"""Synthetic ground-truth data so the whole toolkit runs anywhere — no downloads.

The generator builds a tiny annotated "genome" and then emits per-sample splice
junctions that mimic what STAR reports, including:

* **canonical** introns (the known annotation),
* **true cryptic-exon events** — biologically faithful: a cryptic exon inside an
  intron produces a ``novel_acceptor`` junction (sharing the upstream *known*
  donor) and a ``novel_donor`` junction (sharing the downstream *known*
  acceptor). These recur across replicates, use a canonical GT/AG motif, and are
  **up-regulated in condition B**, and
* **noise** novel junctions — sporadic, low-support, often non-canonical.

Every junction carries an ``is_cryptic_truth`` label, which is what the ML module
is trained to recover. Fully reproducible via ``seed``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

_NONCANONICAL = ["non-canonical", "GC/AG", "AT/AC", "GT/AT"]


@dataclass
class SimulatedDataset:
    known: pd.DataFrame
    observed: pd.DataFrame
    groups: dict[str, str] = field(default_factory=dict)


def simulate_dataset(
    n_genes: int = 8,
    exons_per_gene: int = 5,
    n_per_group: int = 4,
    cryptic_fraction: float = 0.5,
    alt_ss_fraction: float = 0.0,
    label_noise: float = 0.0,
    seed: int = 0,
) -> SimulatedDataset:
    """Generate a reproducible splicing dataset with known ground truth.

    ``label_noise`` flips ``is_cryptic_truth`` for that fraction of junctions,
    mimicking imperfect curation. Left at 0 the labels are exact (deterministic);
    a small value (e.g. 0.1) makes the downstream ML task realistically hard.
    """
    rng = np.random.default_rng(seed)
    chrom, strand = "chr1", "+"
    exon_len, intron_len, gene_gap = 120, 500, 3000

    known_rows: list[tuple] = []
    gene_introns: dict[str, list[tuple[int, int]]] = {}
    pos = 1000
    for g in range(n_genes):
        gene_id = f"g{g:02d}"
        exons = []
        for _ in range(exons_per_gene):
            start, end = pos, pos + exon_len - 1
            exons.append((start, end))
            pos = end + intron_len + 1
        introns = []
        for (_, e1), (s2, _) in zip(exons[:-1], exons[1:], strict=False):
            i_start, i_end = e1 + 1, s2 - 1
            known_rows.append((chrom, i_start, i_end, strand, gene_id))
            introns.append((i_start, i_end))
        gene_introns[gene_id] = introns
        pos += gene_gap

    known = pd.DataFrame(known_rows, columns=["chrom", "start", "end", "strand", "gene_id"])

    samples = [f"A{i}" for i in range(n_per_group)] + [f"B{i}" for i in range(n_per_group)]
    groups = {s: ("A" if s.startswith("A") else "B") for s in samples}

    records: list[dict] = []

    def emit(start, end, motif, truth, per_sample_counts):
        for s, c in per_sample_counts.items():
            if c <= 0:
                continue
            records.append(
                {
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "motif": motif,
                    "annotated_star": 0,
                    "count": int(c),
                    "sample": s,
                    "is_cryptic_truth": int(truth),
                }
            )

    # 1) canonical introns — well expressed everywhere
    canonical_base = {}
    for _, i_start, i_end, _, _ in known.itertuples(index=False):
        counts = {s: int(rng.poisson(200)) for s in samples}
        canonical_base[(i_start, i_end)] = counts
        emit(i_start, i_end, "GT/AG", 0, counts)

    # 2) true cryptic-exon events in a subset of genes (B-upregulated).
    # Effect sizes vary — some events are strong, some are subtle — so the
    # classification task is non-trivial rather than perfectly separable.
    cryptic_genes = [g for g in gene_introns if rng.random() < cryptic_fraction]
    for gene_id in cryptic_genes:
        introns = gene_introns[gene_id]
        if not introns:
            continue
        i_start, i_end = introns[rng.integers(len(introns))]
        if i_end - i_start < 200:
            continue
        c_start = i_start + 100
        c_end = c_start + 60
        strength = rng.uniform(0.35, 1.0)  # subtle -> strong
        b_lambda = 15 + 55 * strength
        incl = {
            s: int(rng.poisson(b_lambda if groups[s] == "B" else 6)) for s in samples
        }
        # a minority of genuine events use a non-canonical motif
        motif = "GT/AG" if rng.random() > 0.15 else _NONCANONICAL[rng.integers(len(_NONCANONICAL))]
        emit(i_start, c_start - 1, motif, 1, incl)  # novel_acceptor (shares known donor)
        emit(c_end + 1, i_end, motif, 1, incl)      # novel_donor (shares known acceptor)

    # 2b) alternative 5′/3′ splice-site events (B-upregulated alternative usage).
    alt_genes = [g for g in gene_introns if rng.random() < alt_ss_fraction]
    for gene_id in alt_genes:
        introns = gene_introns[gene_id]
        if not introns:
            continue
        i_start, i_end = introns[rng.integers(len(introns))]
        if i_end - i_start < 120:
            continue
        delta = int(rng.integers(20, 60))
        alt = {s: int(rng.poisson(50 if groups[s] == "B" else 12)) for s in samples}
        if rng.random() < 0.5:
            emit(i_start + delta, i_end, "GT/AG", 0, alt)  # A5SS: alt donor, shared acceptor
        else:
            emit(i_start, i_end - delta, "GT/AG", 0, alt)  # A3SS: shared donor, alt acceptor

    # 3) noise novel junctions — mostly sporadic/low, but a fraction mimic real
    # events (canonical motif, recurrent support) so classes overlap, truth=0.
    n_noise = max(6, n_genes * 3)
    for _ in range(n_noise):
        gene_id = list(gene_introns)[rng.integers(n_genes)]
        introns = gene_introns[gene_id]
        i_start, i_end = introns[rng.integers(len(introns))]
        if i_end - i_start < 80:
            continue
        ns = i_start + int(rng.integers(20, max(21, i_end - i_start - 40)))
        ne = ns + int(rng.integers(20, 60))
        deceptive = rng.random() < 0.35
        if deceptive:
            # a hard negative: shares the intron's known donor (so it looks like a
            # cryptic novel_acceptor), canonical motif, recurrent moderate support
            motif = "GT/AG"
            hit = rng.choice(samples, size=int(rng.integers(2, len(samples))), replace=False)
            counts = {s: (int(rng.poisson(9)) if s in hit else 0) for s in samples}
            emit(i_start, ne, motif, 0, counts)  # known donor -> novel acceptor
        else:
            motif = _NONCANONICAL[rng.integers(len(_NONCANONICAL))]
            hit = rng.choice(samples, size=int(rng.integers(1, 3)), replace=False)
            counts = {s: (int(rng.poisson(3)) if s in hit else 0) for s in samples}
            emit(ns, ne, motif, 0, counts)  # both-novel sporadic artefact

    observed = pd.DataFrame.from_records(records)

    if label_noise > 0 and not observed.empty:
        jid = (
            observed["chrom"].astype(str)
            + ":"
            + observed["start"].astype(str)
            + "-"
            + observed["end"].astype(str)
            + observed["strand"].astype(str)
        )
        observed = observed.assign(_jid=jid)
        unique_ids = observed["_jid"].drop_duplicates().to_numpy()
        n_flip = int(round(len(unique_ids) * label_noise))
        if n_flip > 0:
            flip_ids = set(rng.choice(unique_ids, size=n_flip, replace=False))
            mask = observed["_jid"].isin(flip_ids)
            observed.loc[mask, "is_cryptic_truth"] = 1 - observed.loc[mask, "is_cryptic_truth"]
        observed = observed.drop(columns="_jid")

    return SimulatedDataset(known=known, observed=observed, groups=groups)


def write_dataset(ds: SimulatedDataset, outdir: str | Path) -> Path:
    """Write a simulated dataset to disk as a GTF, per-sample SJ.out.tab and groups.tsv."""
    outdir = Path(outdir)
    (outdir / "sj").mkdir(parents=True, exist_ok=True)

    # minimal GTF (exons implied by introns: reconstruct exon blocks per gene)
    gtf_lines = []
    for gene_id, sub in ds.known.groupby("gene_id"):
        introns = sorted(zip(sub["start"], sub["end"], strict=False))
        # reconstruct exon blocks that flank the known introns (120 bp flanks)
        exon_coords = []
        prev_end = introns[0][0] - 121
        for i_start, i_end in introns:
            exon_coords.append((prev_end + 1, i_start - 1))
            prev_end = i_end
        exon_coords.append((prev_end + 1, prev_end + 120))
        for es, ee in exon_coords:
            attrs = f'gene_id "{gene_id}"; transcript_id "{gene_id}.t1";'
            gtf_lines.append(f"chr1\tsim\texon\t{es}\t{ee}\t.\t+\t.\t{attrs}")
    (outdir / "annotation.gtf").write_text("\n".join(gtf_lines) + "\n")

    star_cols = [
        "chrom", "start", "end", "strand", "motif", "annotated", "n_unique", "n_multi", "oh",
    ]
    strand_code = {"+": 1, "-": 2, ".": 0}
    motif_code = {"GT/AG": 1, "CT/AC": 2, "GC/AG": 3, "AT/AC": 5, "non-canonical": 0, "GT/AT": 6}
    for sample, sub in ds.observed.groupby("sample"):
        rows = []
        for r in sub.itertuples(index=False):
            rows.append(
                [
                    r.chrom,
                    r.start,
                    r.end,
                    strand_code.get(r.strand, 0),
                    motif_code.get(r.motif, 0),
                    0,
                    r.count,
                    0,
                    30,
                ]
            )
        pd.DataFrame(rows, columns=star_cols).to_csv(
            outdir / "sj" / f"{sample}.SJ.out.tab", sep="\t", header=False, index=False
        )

    pd.DataFrame(
        {"sample": list(ds.groups), "condition": list(ds.groups.values())}
    ).to_csv(outdir / "groups.tsv", sep="\t", index=False)
    return outdir
