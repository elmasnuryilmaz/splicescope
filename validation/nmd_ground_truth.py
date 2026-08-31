#!/usr/bin/env python3
"""Score consequence predictions against experimentally measured NMD sensitivity.

The prediction says a cryptic exon introduces a premature termination codon and
that the transcript is therefore degraded. Whether that is true is not something
our own statistics can answer, so the label is taken from a different kind of
experiment: cryptic exons in TDP-43-depleted i3Neurons were re-measured after
knocking down NMD factors (XRN1, UPF1, SMG6). If inclusion rises when decay is
blocked, that transcript was being degraded.

Ground truth: Table 1 of doi:10.1101/2025.06.28.661837 (supplementary xlsx),
which lists 421 cryptic exons with hg38 exon coordinates and per-condition PSI.

    python validation/nmd_ground_truth.py \\
        --table Table1.xlsx --gtf gencode.gtf.gz --genome genome.fa \\
        --out results/nmd_validation.tsv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

from splicescope.consequence import GenomeFasta, annotate_consequences, load_transcripts

#: Classes in which no premature stop is predicted, so NMD should not apply.
NO_PTC = ["in_frame_insertion", "utr_insertion", "non_coding_host"]
KD_COLUMN = "avgPSI_sgTARDBP"


def load_ground_truth(path: Path, event_type: str = "cassette") -> pd.DataFrame:
    d = pd.read_excel(path, sheet_name="Table1")
    location = next(c for c in d.columns if "Cryptic Exon Location" in c)
    coords = d[location].astype(str).str.extract(r"(chr[\w]+):(\d+)-(\d+)")
    d["chrom"] = coords[0]
    d["exon_start"] = pd.to_numeric(coords[1])
    d["exon_end"] = pd.to_numeric(coords[2])
    d["gene_id"] = d["Ensembl ID"]

    # Coordinates are 1-based inclusive; verified against the stated exon length.
    stated = pd.to_numeric(d["Cryptic Exon (bp)"], errors="coerce")
    length = d.exon_end - d.exon_start + 1
    agree = (stated == length).sum() / stated.notna().sum()
    if agree < 0.95:
        raise ValueError(f"exon coordinates do not look 1-based inclusive ({agree:.0%} agree)")

    nmd_columns = [c for c in d.columns if c.startswith(KD_COLUMN + "_")]
    blocked = d[nmd_columns].apply(pd.to_numeric, errors="coerce").mean(axis=1)
    d["nmd_gain"] = blocked - pd.to_numeric(d[KD_COLUMN], errors="coerce")
    kept = d[(d["type"] == event_type) & d.chrom.notna()].copy()
    return kept


def evaluate(d: pd.DataFrame) -> dict:
    d = d.dropna(subset=["nmd_gain"])
    predicted = d[d.consequence_class == "ptc_nmd"].nmd_gain
    other = d[d.consequence_class.isin(NO_PTC)].nmd_gain
    scored = d[d.consequence_class.isin(["ptc_nmd"] + NO_PTC)]
    labels = scored.consequence_class.eq("ptc_nmd").to_numpy()
    _, pvalue = stats.mannwhitneyu(predicted, other, alternative="greater")

    rng = np.random.default_rng(0)
    boot = []
    for _ in range(2000):
        idx = rng.choice(len(scored), len(scored), replace=True)
        if labels[idx].std() > 0:
            boot.append(roc_auc_score(labels[idx], scored.nmd_gain.to_numpy()[idx]))
    low, high = np.percentile(boot, [2.5, 97.5])

    return {
        "n_ptc_nmd": len(predicted),
        "median_ptc_nmd": predicted.median(),
        "n_no_ptc": len(other),
        "median_no_ptc": other.median(),
        "pvalue": pvalue,
        "auroc": roc_auc_score(labels, scored.nmd_gain),
        "auroc_ci": (low, high),
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--table", required=True, type=Path, help="supplementary Table 1 (.xlsx)")
    p.add_argument("--gtf", required=True, type=Path)
    p.add_argument("--genome", required=True, type=Path, help="indexed FASTA (.fai required)")
    p.add_argument("--out", required=True, type=Path)
    args = p.parse_args()

    truth = load_ground_truth(args.table)
    transcripts = load_transcripts(args.gtf, genes=set(truth.gene_id.dropna().astype(str)))
    strands = {}
    for tx in transcripts.values():
        for key in {tx.gene_id, tx.gene_id.split(".")[0]}:
            strands.setdefault(key, tx.strand)
    truth["strand"] = truth.gene_id.astype(str).map(strands)
    truth = truth[truth.strand.notna()]

    with GenomeFasta(args.genome) as genome:
        scored = annotate_consequences(
            truth, transcripts, genome, start_col="exon_start", end_col="exon_end"
        )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    scored.to_csv(args.out, sep="\t", index=False)

    m = evaluate(scored)
    print(f"wrote {args.out}")
    print(f"  ptc_nmd  n={m['n_ptc_nmd']:>3}  median ΔPSI on NMD block = {m['median_ptc_nmd']:.2f}")
    print(f"  no PTC   n={m['n_no_ptc']:>3}  median ΔPSI on NMD block = {m['median_no_ptc']:.2f}")
    print(f"  Mann-Whitney one-sided p = {m['pvalue']:.3g}")
    print(f"  AUROC = {m['auroc']:.3f}  95% CI [{m['auroc_ci'][0]:.3f}, {m['auroc_ci'][1]:.3f}]")


if __name__ == "__main__":
    main()
