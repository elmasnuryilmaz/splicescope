#!/usr/bin/env python3
"""Does the test suite actually notice when the code is wrong?

A green suite says the tests pass, not that they would fail if the code broke. The only
way to find out is to break it: change one line to something a person could plausibly
have written, run the suite, and see whether anything fails. A defect nothing fails on is
a rule the suite is not checking, and it is free to come back in the next refactor.

Every mutation below was written as a realistic mistake rather than a random edit — an
inclusive comparison where it should be strict, a distance measured from the wrong end of
a codon, a correction applied to the wrong array. Thirteen of these passed all 165 tests
when the survey was first run, including the 50-nucleotide rule the tool is built on,
the last-exon exception, and whether the reported q-value is Benjamini-Hochberg-adjusted
at all. Those are now covered.

Two mutants survive and are expected to: they are equivalent, not uncaught.

  - ``find_ptc`` stopping at ``len(sequence)`` instead of ``len(sequence) - 2`` only ever
    compares slices shorter than three characters, which cannot equal a stop codon.
  - Leaving the likelihood-ratio statistic unclamped gives the same p-value, because
    ``chi2.sf`` of a negative value is 1.0. Only the reported statistic would look wrong.

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
