#!/usr/bin/env python3
"""Does the test suite actually notice when the code is wrong?

A green suite says the tests pass, not that they would fail if the code broke. The only
way to find out is to break it: change one line to something a person could plausibly
have written, run the suite, and see whether anything fails. A defect nothing fails on is
a rule the suite is not checking, and it is free to come back in the next refactor.

Every mutation below was written as a realistic mistake rather than a random edit — an
inclusive comparison where it should be strict, a distance measured from the wrong end of
a codon, a correction applied to the wrong array. Twenty-eight of these passed the suite
as it stood when each was first tried: the 50-nucleotide rule the tool is built on, the
last-exon exception, whether the reported q-value is Benjamini-Hochberg-adjusted at all,
seven of the eight features the cryptic-junction classifier learns from, and three
defects reachable only from the minus strand, which the consequence layer was hardly
tested on. Three of the twenty-eight are equivalent mutants; the other twenty-five are
now each pinned by a test.

Three mutants survive and are expected to: they are equivalent, not uncaught.

  - ``find_ptc`` stopping at ``len(sequence)`` instead of ``len(sequence) - 2`` only ever
    compares slices shorter than three characters, which cannot equal a stop codon.
  - Leaving the likelihood-ratio statistic unclamped gives the same p-value, because
    ``chi2.sf`` of a negative value is 1.0. Only the reported statistic would look wrong.
  - Testing the far end of an exon rather than its near end in
    ``_downstream_exon_lengths`` agrees for every exon, because the event always lies
    inside an intron of the transcript and so no exon spans it.

    python validation/mutation_survey.py                 # every mutation
    python validation/mutation_survey.py --module enrich  # one module's
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src" / "splicescope"

#: ``(module, label, old, new, expected)``. ``expected`` is "caught" for a real defect
#: and "equivalent" for a mutant that provably cannot change any answer.
MUTATIONS: list[tuple[str, str, str, str, str]] = [
    # ---- consequence: the NMD rule ------------------------------------------------
    ("consequence", "50-nt rule boundary made inclusive",
     "distance, distance > NMD_DISTANCE_RULE",
     "distance, distance >= NMD_DISTANCE_RULE", "caught"),
    ("consequence", "last exon's length counted into the distance",
     "last_junction = sum(lengths[:-1]) if len(lengths) > 1 else None",
     "last_junction = sum(lengths) if len(lengths) > 1 else None", "caught"),
    ("consequence", "distance measured from the start of the stop codon",
     "distance = last_junction - (ptc + 3)",
     "distance = last_junction - ptc", "caught"),
    ("consequence", "last-exon exception removed",
     "last_junction = sum(lengths[:-1]) if len(lengths) > 1 else None",
     "last_junction = sum(lengths[:-1]) if len(lengths) > 0 else None", "caught"),
    ("consequence", "reading frame not inherited across the junction",
     "    begin = (3 - frame_offset) % 3",
     "    begin = frame_offset % 3", "caught"),
    ("consequence", "frame ignored entirely",
     "    begin = (3 - frame_offset) % 3", "    begin = 0", "caught"),
    ("consequence", "native stop offset off by one codon",
     "return remaining - 3 if remaining >= 3 else None",
     "return remaining if remaining >= 3 else None", "caught"),
    ("consequence", "a partial codon at the end may be a stop",
     "for i in range(begin, len(sequence) - 2, 3):",
     "for i in range(begin, len(sequence), 3):", "equivalent"),
    # ---- consequence: transcript structure ----------------------------------------
    ("consequence", "GTF phase added instead of subtracted",
     "return (self.cds_length_before(position) - self.cds_phase) % 3",
     "return (self.cds_length_before(position) + self.cds_phase) % 3", "caught"),
    ("consequence", "GTF phase ignored, so a 5'-incomplete CDS starts mid-codon",
     "return (self.cds_length_before(position) - self.cds_phase) % 3",
     "return self.cds_length_before(position) % 3", "caught"),
    ("consequence", "minus-strand blocks not put in transcription order",
     "        self.exons = sorted(self.exons, reverse=reverse)",
     "        self.exons = sorted(self.exons)", "caught"),
    ("consequence", "intron bounds taken as the exon bounds",
     "out = [(a[1] + 1, b[0] - 1) for a, b in zip(blocks, blocks[1:], strict=False)]",
     "out = [(a[1], b[0]) for a, b in zip(blocks, blocks[1:], strict=False)]", "caught"),
    ("consequence", "introns of a minus-strand transcript left in genomic order",
     'return out[::-1] if self.strand == "-" else out', "return out", "caught"),
    ("consequence", "coding length before a position off by one on the plus strand",
     "                    total += position - cstart",
     "                    total += position - cstart + 1", "caught"),
    ("consequence", "coding length before a position off by one on the minus strand",
     "                    total += cend - position",
     "                    total += cend - position + 1", "caught"),
    # ---- consequence: the sequence a ribosome would read ---------------------------
    ("consequence", "downstream sequence not clipped at the event",
     "            start = max(estart, resume)",
     "            start = estart", "caught"),
    ("consequence", "the exon holding the resume point dropped",
     "            if eend < resume:\n                continue",
     "            if eend <= resume:\n                continue", "caught"),
    ("consequence", "downstream exon lengths off by one",
     "        lengths.append(end - start + 1)",
     "        lengths.append(end - start)", "caught"),
    ("consequence", "exons after the event found by the wrong end",
     "        after = estart > end if tx.strand == \"+\" else eend < start",
     "        after = eend > end if tx.strand == \"+\" else estart < start", "equivalent"),
    ("consequence", "transcription resumes at the furthest exon, not the next one",
     "        return min(later) if later else None",
     "        return max(later) if later else None", "caught"),
    # ---- consequence: what the junction changes ------------------------------------
    ("consequence", "extension and truncation swapped at the 3' end",
     '            if end < iend:\n                return "extension", end + 1, iend\n'
     '            return "truncation", iend + 1, end',
     '            if end < iend:\n                return "truncation", end + 1, iend\n'
     '            return "extension", iend + 1, end', "caught"),
    ("consequence", "the exon-skip guard removed",
     "        if start in starts and end in ends:\n            return None",
     "        if start in starts and end in ends:\n            pass", "caught"),
    ("consequence", "the extended stretch off by one at the 5' end",
     '                return "extension", istart, start - 1',
     '                return "extension", istart, start', "caught"),
    ("consequence", "a host intron need only overlap the event, not contain it",
     "            if istart <= start and end <= iend:",
     "            if istart <= end and start <= iend:", "caught"),
    # ---- betabinom: the statistics ------------------------------------------------
    ("betabinom", "ungrouped samples left in the null",
     "    mask = mask & (group_a | group_b)[None, :]\n", "", "caught"),
    ("betabinom", "likelihood ratio not doubled",
     "statistic = np.maximum(2.0 * (ll_alt - ll_null), 0.0)",
     "statistic = np.maximum(ll_alt - ll_null, 0.0)", "caught"),
    ("betabinom", "chi-square on two degrees of freedom",
     "pvalue = stats.chi2.sf(statistic, df=1)",
     "pvalue = stats.chi2.sf(statistic, df=2)", "caught"),
    ("betabinom", "the group means swapped in the alternative",
     "ll_alt = loglik(k, n, mask_a, mu_a, s) + loglik(k, n, mask_b, mu_b, s)",
     "ll_alt = loglik(k, n, mask_a, mu_b, s) + loglik(k, n, mask_b, mu_a, s)", "caught"),
    ("betabinom", "rank-test floor made one-sided",
     "return 2.0 / comb(n_a + n_b, n_a)",
     "return 1.0 / comb(n_a + n_b, n_a)", "caught"),
    ("betabinom", "degrees-of-freedom correction dropped",
     "    corrected = pearson * n_total / df", "    corrected = pearson", "caught"),
    ("betabinom", "trials-excess forgets the -1",
     "trials_excess = float(np.where(keep, n - 1.0, 0.0).sum())",
     "trials_excess = float(np.where(keep, n, 0.0).sum())", "caught"),
    ("betabinom", "precision bounds removed",
     "    return float(np.clip(s, min_precision, max_precision))",
     "    return float(s)", "caught"),
    ("betabinom", "statistic left unclamped",
     "statistic = np.maximum(2.0 * (ll_alt - ll_null), 0.0)",
     "statistic = 2.0 * (ll_alt - ll_null)", "equivalent"),
    # ---- diff ---------------------------------------------------------------------
    ("diff", "Psi pivoted by sum instead of mean", 'aggfunc="mean"', 'aggfunc="sum"',
     "caught"),
    ("diff", "raw p-values reported as q-values",
     'res["qvalue"] = benjamini_hochberg(res["pvalue"].to_numpy())',
     'res["qvalue"] = res["pvalue"].to_numpy()', "caught"),
    # ---- quantify -----------------------------------------------------------------
    ("quantify", "coverage threshold made inclusive",
     "usage[totals < min_reads] = np.nan",
     "usage[totals <= min_reads] = np.nan", "caught"),
    ("quantify", "an omitted junction filled with one read",
     'missing["count"] = 0', 'missing["count"] = 1', "caught"),
    ("quantify", "zeros filled from the donor site only",
     'for site in ("donor", "acceptor")', 'for site in ("donor",)', "caught"),
    ("quantify", "the zero grid ignores which sample lacks the junction",
     'grid = pd.concat(by_site, ignore_index=True).drop_duplicates(JUNCTION_KEY + ["sample"])',
     "grid = pd.concat(by_site, ignore_index=True).drop_duplicates(JUNCTION_KEY)", "caught"),
    # ---- annotate -----------------------------------------------------------------
    ("annotate", "novel_donor and novel_acceptor swapped",
     '    if d_known and not a_known:\n        return "novel_acceptor"',
     '    if d_known and not a_known:\n        return "novel_donor"', "caught"),
    ("annotate", "is_novel inverted",
     'out["is_novel"] = out["sclass"] != "annotated"',
     'out["is_novel"] = out["sclass"] == "annotated"', "caught"),
    ("annotate", "an unplaceable junction guessed onto the + strand",
     '            or "."', '            or "+"', "caught"),
    ("annotate", "the chromosome dropped from the annotation index",
     "    junctions = set(zip(chroms, starts, ends, strands, strict=True))",
     "    junctions = set(zip(starts, ends, strands, strict=True))", "caught"),
    # ---- io -----------------------------------------------------------------------
    ("io", "STAR strand code 0 read as + rather than unknown",
     '_STAR_STRAND = {0: ".", 1: "+", 2: "-"}',
     '_STAR_STRAND = {0: "+", 1: "+", 2: "-"}', "caught"),
    ("io", "GTF exon boundaries used as intron boundaries",
     "            intron_start = e_end + 1\n            intron_end = n_start - 1",
     "            intron_start = e_end\n            intron_end = n_start", "caught"),
    ("io", "donor and acceptor not swapped on the minus strand",
     '    if strand == "-":\n        return end, start',
     '    if strand == "-":\n        return start, end', "caught"),
    # ---- enrich -------------------------------------------------------------------
    ("enrich", "hypergeometric tail off by one",
     "p = float(stats.hypergeom.sf(k - 1, M, n, N))",
     "p = float(stats.hypergeom.sf(k, M, n, N))", "caught"),
    ("enrich", "Wallenius tail off by one",
     "p = float(stats.nchypergeom_wallenius.sf(k - 1, M, n, N, odds))",
     "p = float(stats.nchypergeom_wallenius.sf(k, M, n, N, odds))", "caught"),
    ("enrich", "fold enrichment inverted",
     "fold = (k / N) / (n / M)", "fold = (n / M) / (k / N)", "caught"),
    ("enrich", "Jeffreys smoothing replaced by clipping",
     "        rate = (observed + _PSEUDO_HITS) / (len(members) + _PSEUDO_TOTAL)",
     "        rate = min(observed / len(members), 1.0 - 1e-6)", "caught"),
    ("enrich", "the remainder bin folded into the previous one",
     "        chunks.append(current)\n\n    propensity",
     "        chunks[-1].extend(current) if chunks else chunks.append(current)"
     "\n\n    propensity", "caught"),
    ("enrich", "probabilities averaged instead of odds",
     "        inside_sum = math.fsum(odds[g] for g in members)",
     "        inside_sum = math.fsum(odds[g] / (1 + odds[g]) for g in members)", "caught"),
    # ---- events -------------------------------------------------------------------
    ("events", "the chromosome dropped from the cassette key",
     "        by_start[(row.chrom, row.strand, row.start)].append(row.end)",
     "        by_start[(row.strand, row.start)].append(row.end)", "caught"),
    ("events", "the chromosome dropped from the MXE key",
     "        by_cs[(r.chrom, r.strand)].append(r)",
     "        by_cs[(r.strand,)].append(r)", "caught"),
    ("events", "the chromosome dropped from the alternative-site key",
     "        groups[(row.chrom, getattr(row, shared), row.strand)].append(row)",
     "        groups[(getattr(row, shared), row.strand)].append(row)", "caught"),
    ("events", "mutually exclusive exons allowed to overlap",
     "if a[1] < b[0]:  # exon A strictly upstream of exon B",
     "if a[1] <= b[0] + 50:  # exon A strictly upstream of exon B", "caught"),
    ("events", "the longer exon chosen as the inclusion isoform",
     "incl = min(members, key=lambda m: abs(getattr(m, variable) - site))",
     "incl = max(members, key=lambda m: abs(getattr(m, variable) - site))", "caught"),
    ("events", "the candidate-exon window widened without bound",
     "hi = bisect_right(starts, j5.end + 1 + max_exon)", "hi = len(starts)", "caught"),
    ("events", "A5SS and A3SS labels swapped",
     '    shared, variable = ("acceptor", "donor") if kind == "A5SS" else ("donor", "acceptor")',
     '    shared, variable = ("donor", "acceptor") if kind == "A5SS" else ("acceptor", "donor")',
     "caught"),
    ("events", "crowded anchors truncated to a single pair",
     "        plist = plist[:max_candidates]", "        plist = plist[:1]", "caught"),
    # ---- cryptic: the features the classifier learns from -------------------------
    ("cryptic", "nearest known site only searched to the right",
     "    for j in (i - 1, i):", "    for j in (i,):", "caught"),
    ("cryptic", "intron length off by one",
     '"intron_length": int(end - start + 1),',
     '"intron_length": int(end - start),', "caught"),
    ("cryptic", "supporting samples counted including the zeros",
     '"n_samples_support": int((counts > 0).sum()),',
     '"n_samples_support": int((counts >= 0).sum()),', "caught"),
    ("cryptic", "canonical-motif flag inverted",
     '"canonical_motif": int(motif in _CANONICAL),',
     '"canonical_motif": int(motif not in _CANONICAL),', "caught"),
    ("cryptic", "the minus-strand canonical motif dropped",
     '_CANONICAL = {"GT/AG", "CT/AC"}', '_CANONICAL = {"GT/AG"}', "caught"),
    ("cryptic", "only donors indexed as known sites",
     "sites.setdefault((chrom, strand), set()).update((d, a))",
     "sites.setdefault((chrom, strand), set()).update((d,))", "caught"),
    ("cryptic", "truth label taken as the minimum over samples",
     'df.groupby(key, observed=True)["is_cryptic_truth"].max().reset_index(drop=True)',
     'df.groupby(key, observed=True)["is_cryptic_truth"].min().reset_index(drop=True)',
     "caught"),
    ("cryptic", "novel_only keeps the annotated junctions instead",
     '        df = df[df["sclass"] != "annotated"]',
     '        df = df[df["sclass"] == "annotated"]', "caught"),
    # ---- ml: the classifier and how it is assessed ---------------------------------
    ("ml", "probability taken for the wrong class in evaluation",
     'method="predict_proba", n_jobs=None\n        )[:, 1]',
     'method="predict_proba", n_jobs=None\n        )[:, 0]', "caught"),
    ("ml", "probability taken for the wrong class in prediction",
     "        return self.pipeline.predict_proba(x)[:, 1]",
     "        return self.pipeline.predict_proba(x)[:, 0]", "caught"),
    ("ml", "ROC-AUC computed against inverted scores",
     '"roc_auc": float(roc_auc_score(y, proba)),',
     '"roc_auc": float(roc_auc_score(y, 1.0 - proba)),', "caught"),
    ("ml", "more folds than the minority class has members",
     "n_splits = min(self.n_splits, int(np.bincount(y).min()))",
     "n_splits = max(self.n_splits, int(np.bincount(y).min()))", "caught"),
    ("ml", "class imbalance left uncorrected",
     '    class_weight: str | None = "balanced"',
     "    class_weight: str | None = None", "caught"),
    ("ml", "the model card repeats a literal instead of reading the model",
     '            "hyperparameters": self._hyperparameters(),',
     '            "hyperparameters": {"class_weight": "balanced"},', "caught"),
    ("ml", "score table sorted least-likely first",
     'return out.sort_values("cryptic_score", ascending=False)[',
     'return out.sort_values("cryptic_score", ascending=True)[', "caught"),
    # ---- cli ------------------------------------------------------------------------
    ("cli", "splice-site shifts restricted to one of the two classes",
     'SHIFT_CLASSES = ("novel_donor", "novel_acceptor")',
     'SHIFT_CLASSES = ("novel_donor",)', "caught"),
    ("cli", "junctions already explained by an event are not excluded",
     "                & ~_junction_index(annotated).isin(_event_junctions(evs))",
     "                & _junction_index(annotated).notna()", "caught"),
    ("cli", "exon-level consequences run on every event type but SE",
     '        table = events[events["event_type"] == "SE"].copy()',
     '        table = events[events["event_type"] != "SE"].copy()', "caught"),
    # ---- plotting: the figures a reader looks at first ------------------------------
    ("plotting", "volcano significance needs only one of the two thresholds",
     '    sig = (d["qvalue"] <= q) & (d["delta_psi"].abs() >= min_delta)',
     '    sig = (d["qvalue"] <= q) | (d["delta_psi"].abs() >= min_delta)', "caught"),
    ("plotting", "volcano q-value threshold inverted",
     '    sig = (d["qvalue"] <= q) & (d["delta_psi"].abs() >= min_delta)',
     '    sig = (d["qvalue"] >= q) & (d["delta_psi"].abs() >= min_delta)', "caught"),
    ("plotting", "the highlighted volcano points are the wrong ones",
     'ax.scatter(d.loc[sig, "delta_psi"], y[sig], s=22, color=_ACCENT,',
     'ax.scatter(d.loc[~sig, "delta_psi"], y[~sig], s=22, color=_ACCENT,', "caught"),
    ("plotting", "volcano y axis not negated, so the best hits sink",
     '    y = -np.log10(d["qvalue"].clip(lower=1e-300))\n    sig =',
     '    y = np.log10(d["qvalue"].clip(lower=1e-300))\n    sig =', "caught"),
    ("plotting", "the event volcano draws one event type in every colour",
     '        m = d["event_type"] == etype',
     '        m = d["event_type"] == next(iter(_EVENT_COLORS))', "caught"),
]


def survey(mutations, quiet: bool) -> int:
    originals = {m[0]: (SRC / f"{m[0]}.py").read_text() for m in mutations}
    surprises = []
    try:
        for module, label, old, new, expected in mutations:
            path = SRC / f"{module}.py"
            source = originals[module]
            if source.count(old) != 1:
                verdict = f"NOT APPLIED ({source.count(old)} matches)"
                surprises.append((module, label, verdict))
            else:
                path.write_text(source.replace(old, new))
                done = subprocess.run(
                    [sys.executable, "-m", "pytest", "-x", "-q", "tests"],
                    cwd=ROOT, capture_output=True, text=True,
                )
                path.write_text(source)
                caught = done.returncode != 0
                got = "caught" if caught else "survived"
                want = "caught" if expected == "caught" else "survived"
                verdict = got if got == want else f"{got.upper()}, expected {want}"
                if got != want:
                    surprises.append((module, label, verdict))
            if not quiet:
                print(f"  {module:<12} {label:<58} {verdict}")
    finally:
        for module, source in originals.items():
            (SRC / f"{module}.py").write_text(source)

    print(f"\n{len(mutations)} mutations, {len(surprises)} surprise(s)")
    for module, label, verdict in surprises:
        print(f"  {module}: {label} -> {verdict}")
    return 1 if surprises else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--module", help="restrict to one module (e.g. enrich)")
    parser.add_argument("--quiet", action="store_true", help="only report surprises")
    args = parser.parse_args()

    chosen = [m for m in MUTATIONS if args.module in (None, m[0])]
    if not chosen:
        modules = sorted({m[0] for m in MUTATIONS})
        parser.error(f"no mutations for {args.module!r}; have {modules}")
    print(f"{len(chosen)} mutation(s); each runs the whole suite, so this takes a while.\n")
    return survey(chosen, args.quiet)


if __name__ == "__main__":
    raise SystemExit(main())
