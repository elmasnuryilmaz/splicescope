#!/usr/bin/env python3
"""Reproduce every number reported in the NMD-prediction manuscript.

The three tables in ``manuscript_tables/`` are the per-event outputs of
``splicescope consequence`` run against the two published ground-truth sets, plus
the LeafCutter2 labels from the full-transcriptome run described in Methods 2.4.
This script takes those tables and emits every quantity that appears in the text,
in section order, so that a reader can check the manuscript line by line.

    python validation/manuscript_numbers.py

Every bootstrap uses ``numpy.random.default_rng(0)`` and 2000 resamples, so the
confidence intervals printed here are exactly those in the manuscript.
"""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

SEED = 0
BOOT = 2000

# Classes in which no premature termination codon can arise. The pooled contrast
# reported throughout is ptc_nmd against this group.
NO_STOP = ["in_frame_insertion", "utr_insertion", "non_coding_host"]

KD_COLUMNS = [
    "avgPSI_sgTARDBP_shUPF1/shSMG6",
    "avgPSI_sgTARDBP_shXRN1",
    "avgPSI_sgTARDBP_shXRN1/shSMG6",
    "avgPSI_sgTARDBP_shXRN1/shUPF1",
]
BASELINE = "avgPSI_sgTARDBP"


def auroc_ci(pos: np.ndarray, neg: np.ndarray) -> tuple[float, float, float]:
    """AUROC of measured dPSI by predicted label, with a stratified bootstrap CI.

    The score is the *measured* dPSI and the label is the *predicted* class, so
    this equals U / (n1 * n2) from the corresponding Mann-Whitney test: the
    probability that a randomly drawn predicted-NMD event gained more inclusion
    on NMD blockade than a randomly drawn no-stop event.
    """
    rng = np.random.default_rng(SEED)
    y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
    point = roc_auc_score(y, np.r_[pos, neg])
    draws = [
        roc_auc_score(
            y,
            np.r_[
                pos[rng.integers(0, len(pos), len(pos))],
                neg[rng.integers(0, len(neg), len(neg))],
            ],
        )
        for _ in range(BOOT)
    ]
    lo, hi = np.percentile(draws, [2.5, 97.5])
    return point, lo, hi


def median_ci(values: np.ndarray) -> tuple[float, float, float]:
    rng = np.random.default_rng(SEED)
    draws = [np.median(rng.choice(values, len(values), replace=True)) for _ in range(BOOT)]
    lo, hi = np.percentile(draws, [2.5, 97.5])
    return float(np.median(values)), lo, hi


def head(title: str) -> None:
    print("\n" + "=" * 74)
    print(title)
    print("=" * 74)


def section_3_1(d: pd.DataFrame) -> pd.DataFrame:
    head("Methods 2.3 / Section 3.1 - cassette cryptic exons (Sinha et al.)")
    print(f"cassette exons in Table 1          : {len(d)}")
    d = d[d.nmd_gain.notna()].copy()
    print(f"with a measurable dPSI (analysed)  : {len(d)}")

    print("\nTable 2 - median dPSI on NMD blockade by predicted class:")
    for cls in [
        "ptc_nmd",
        "ptc_escape",
        "in_frame_insertion",
        "utr_insertion",
        "non_coding_host",
        "no_host_transcript",
    ]:
        v = d.loc[d.consequence_class == cls, "nmd_gain"].to_numpy()
        if not len(v):
            continue
        m, lo, hi = median_ci(v)
        print(f"  {cls:<20} n={len(v):>3}  median {m:7.2f}  95% CI {lo:7.2f} - {hi:6.2f}")

    pos = d.loc[d.consequence_class == "ptc_nmd", "nmd_gain"].to_numpy()
    neg = d.loc[d.consequence_class.isin(NO_STOP), "nmd_gain"].to_numpy()
    auc, lo, hi = auroc_ci(pos, neg)
    p = stats.mannwhitneyu(pos, neg, alternative="two-sided").pvalue
    print(f"\npooled contrast: ptc_nmd n={len(pos)} vs no-stop n={len(neg)}")
    print(f"  medians                : {np.median(pos):.2f} against {np.median(neg):.2f}")
    print(f"  Mann-Whitney two-sided : p = {p:.2e}")
    print(f"  AUROC                  : {auc:.3f}  95% CI {lo:.3f} - {hi:.3f}")

    print("\nrobustness - alternative definitions of NMD sensitivity:")
    aucs, worst_p = [], 0.0
    definitions = [(c.replace("avgPSI_sgTARDBP_", ""), d[c] - d[BASELINE]) for c in KD_COLUMNS]
    definitions += [
        ("mean of the four", d[KD_COLUMNS].mean(axis=1) - d[BASELINE]),
        ("max of the four", d[KD_COLUMNS].max(axis=1) - d[BASELINE]),
    ]
    for label, gain in definitions:
        a = gain[d.consequence_class == "ptc_nmd"].dropna().to_numpy()
        b = gain[d.consequence_class.isin(NO_STOP)].dropna().to_numpy()
        val = roc_auc_score(np.r_[np.ones(len(a)), np.zeros(len(b))], np.r_[a, b])
        pv = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
        aucs.append(val)
        worst_p = max(worst_p, pv)
        print(f"  {label:<20} AUROC {val:.3f}   p = {pv:.2e}")
    print(f"  -> range {min(aucs):.3f} - {max(aucs):.3f}; largest p = {worst_p:.1e}")

    one = d.sort_values("nmd_gain", ascending=False).drop_duplicates("gene_name")
    a = one.loc[one.consequence_class == "ptc_nmd", "nmd_gain"].to_numpy()
    b = one.loc[one.consequence_class.isin(NO_STOP), "nmd_gain"].to_numpy()
    val = roc_auc_score(np.r_[np.ones(len(a)), np.zeros(len(b))], np.r_[a, b])
    print(f"  {'one exon per gene':<20} AUROC {val:.3f}   (n = {len(a)} vs {len(b)})")
    return d


def section_3_4(d: pd.DataFrame) -> None:
    head("Section 3.4 - sensitive, not specific")
    sub = d[d.consequence_class.isin(["ptc_nmd"] + NO_STOP)]
    y = sub.consequence_class.eq("ptc_nmd").to_numpy()
    print(f"pooled-contrast set                : n = {len(sub)}")
    called = (d.consequence_class == "ptc_nmd").sum()
    print(f"called ptc_nmd of all analysed      : {called}/{len(d)} = {called / len(d):.1%}")

    sens, spec = [], []
    for t in np.arange(0, 26, 1.0):
        truth = (sub.nmd_gain > t).to_numpy()
        sens.append((y & truth).sum() / truth.sum())
        spec.append(((~y) & (~truth)).sum() / (~truth).sum())
    print("over dPSI thresholds 0-25 step 1 (the grid plotted in Figure 1B):")
    print(f"  sensitivity {min(sens):.2f} - {max(sens):.2f}")
    print(f"  specificity {min(spec):.2f} - {max(spec):.2f}")
    frac = (d.nmd_gain > 10).mean()
    print(f"  measurably NMD-sensitive at dPSI > 10: {frac:.1%}")

    pos = d.loc[d.consequence_class == "ptc_nmd", "nmd_gain"].to_numpy()
    neg = d.loc[d.consequence_class.isin(NO_STOP), "nmd_gain"].to_numpy()
    wide = d.loc[d.consequence_class.isin(["ptc_nmd", "ptc_escape"]), "nmd_gain"].to_numpy()
    narrow_auc = roc_auc_score(np.r_[np.ones(len(pos)), np.zeros(len(neg))], np.r_[pos, neg])
    wide_auc = roc_auc_score(np.r_[np.ones(len(wide)), np.zeros(len(neg))], np.r_[wide, neg])
    print(f"\nbroadening the predicted-NMD class : AUROC {narrow_auc:.3f} -> {wide_auc:.3f}")

    s = d[(d.consequence_class == "ptc_nmd") & d.distance_to_last_junction.notna()]
    rho, p = stats.spearmanr(s.distance_to_last_junction, s.nmd_gain)
    print(f"distance to last junction vs dPSI  : Spearman rho = {rho:.3f}, "
          f"p = {p:.2f} (n = {len(s)})")
    print(f"stop-codon probability 1-(61/64)^42: {1 - (61 / 64) ** 42:.4f}")


def section_3_2(z: pd.DataFrame) -> None:
    head("Section 3.2 - splice-site shifts (Zeng et al.)")
    print(f"junctions after the gene-name filter: {len(z)}")
    interpretable = z.consequence_class != "no_host_transcript"
    print(f"interpretable                       : {interpretable.sum()} "
          f"= {interpretable.mean():.1%}")
    for verdict, g in z.groupby("verdict"):
        rate = (g.consequence_class != "no_host_transcript").mean()
        print(f"    {verdict:<22} n = {len(g):>4}   interpretable {rate:.1%}")

    hit = int(((z.consequence_class == "ptc_nmd") & (z.nmd_sensitive == 1)).sum())
    n1 = int((z.consequence_class == "ptc_nmd").sum())
    print(f"\nptc_nmd NMD-sensitive               : {hit}/{n1} = {hit / n1:.1%}")

    print("odds ratio against every grouping of the no-stop classes:")
    ors = []
    for size in (1, 2, 3):
        for combo in itertools.combinations(NO_STOP, size):
            m = z.consequence_class.isin(combo)
            c, n2 = int((m & (z.nmd_sensitive == 1)).sum()), int(m.sum())
            odds = (hit * (n2 - c)) / ((n1 - hit) * c)
            p = stats.fisher_exact([[hit, n1 - hit], [c, n2 - c]], alternative="two-sided")[1]
            ors.append(odds)
            star = "  <- reported" if size == 3 else ""
            print(
                f"  {'+'.join(combo):<50} n={n2:>4} sens={c / n2:5.1%} "
                f"OR {odds:5.2f} p={p:.2g}{star}"
            )
    print(f"  -> OR range {min(ors):.1f} - {max(ors):.1f}")


def section_3_3(h: pd.DataFrame) -> None:
    head("Section 3.3 - head to head with LeafCutter2")
    print(f"exons carrying both calls           : {len(h)}  (of 357; {357 - len(h)} unlabelled)")
    print(f"LeafCutter2 labels                  : {h.lc2.value_counts().to_dict()}")

    pos = h.loc[h.consequence_class == "ptc_nmd", "nmd_gain"].to_numpy()
    neg = h.loc[h.consequence_class != "ptc_nmd", "nmd_gain"].to_numpy()
    auc, lo, hi = auroc_ci(pos, neg)
    p = stats.mannwhitneyu(pos, neg, alternative="two-sided").pvalue
    print(f"\n  splicescope ptc_nmd n={len(pos)} median {np.median(pos):.2f} | "
          f"other n={len(neg)} median {np.median(neg):.2f}")
    print(f"    AUROC {auc:.3f} ({lo:.3f}-{hi:.3f})  p = {p:.1e}")

    up = h.loc[h.lc2 == "UP", "nmd_gain"].to_numpy()
    pr = h.loc[h.lc2 == "PR", "nmd_gain"].to_numpy()
    ne = h.loc[h.lc2 == "NE", "nmd_gain"].to_numpy()
    auc, lo, hi = auroc_ci(up, pr)
    p = stats.mannwhitneyu(up, pr, alternative="two-sided").pvalue
    print(f"  LeafCutter2 UP n={len(up)} median {np.median(up):.2f} | "
          f"PR n={len(pr)} median {np.median(pr):.2f} | NE n={len(ne)} median {np.median(ne):.2f}")
    print(f"    AUROC {auc:.3f} ({lo:.3f}-{hi:.3f})  p = {p:.1e}")

    resolved = h[h.lc2.isin(["UP", "PR"])]
    agree = ((resolved.consequence_class == "ptc_nmd") == (resolved.lc2 == "UP")).mean()
    print(f"  agreement on the {len(resolved)} exons LeafCutter2 resolves: {agree:.1%}")
    print(f"  called intergenic, in neither comparison: {(h.lc2 == 'IN').sum()}")


def main() -> None:
    here = Path(__file__).resolve().parent / "manuscript_tables"
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tables", type=Path, default=here)
    args = p.parse_args()

    cassettes = pd.read_csv(args.tables / "i3n_prediction_vs_truth.tsv", sep="\t")
    shifts = pd.read_csv(args.tables / "smg1i_junction_consequence.tsv", sep="\t")
    head_to_head = pd.read_csv(args.tables / "headtohead_fulltranscriptome.tsv", sep="\t")

    analysed = section_3_1(cassettes)
    section_3_2(shifts)
    section_3_3(head_to_head)
    section_3_4(analysed)
    print()


if __name__ == "__main__":
    main()
