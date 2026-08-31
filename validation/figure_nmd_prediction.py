#!/usr/bin/env python3
"""Figure: does a sequence-based PTC/NMD call predict measured NMD sensitivity?

Ground truth comes from an orthogonal experiment, not from our own statistics:
cryptic exons in TDP-43-depleted i3Neurons were re-measured after knocking down
NMD factors, so the rise in inclusion when decay is blocked is a direct readout
of how NMD-sensitive each exon was. We ask whether ``consequence`` recovers that
from sequence and annotation alone.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ORDER = [
    "ptc_nmd",
    "frameshift",
    "ptc_escape",
    "in_frame_insertion",
    "utr_insertion",
    "non_coding_host",
]
LABEL = {
    "ptc_nmd": "PTC → NMD",
    "frameshift": "frameshift",
    "ptc_escape": "PTC, escapes",
    "in_frame_insertion": "in frame",
    "utr_insertion": "UTR",
    "non_coding_host": "non-coding host",
}
PREDICTED_NMD = "#1b6ca8"
PREDICTED_NOT = "#b23a48"
NEUTRAL = "#8a8a8a"


def panel_a(ax, d: pd.DataFrame) -> None:
    groups = [g for g in ORDER if (d.consequence_class == g).sum() >= 3]
    data = [d.loc[d.consequence_class == g, "nmd_gain"].to_numpy() for g in groups]
    colours = [
        PREDICTED_NMD
        if g == "ptc_nmd"
        else (NEUTRAL if g in ("frameshift", "ptc_escape") else PREDICTED_NOT)
        for g in groups
    ]
    bp = ax.boxplot(data, vert=False, widths=0.6, patch_artist=True, showfliers=False)
    for patch, colour in zip(bp["boxes"], colours, strict=False):
        patch.set_facecolor(colour)
        patch.set_alpha(0.35)
        patch.set_edgecolor(colour)
    for element in ("medians", "whiskers", "caps"):
        for line in bp[element]:
            line.set_color("0.3")
    rng = np.random.default_rng(0)
    for i, values in enumerate(data, start=1):
        ax.scatter(values, i + rng.uniform(-0.16, 0.16, len(values)), s=5, alpha=0.4,
                   color=colours[i - 1], linewidths=0)
    ax.axvline(0, color="0.6", lw=1, ls="--")
    ax.set_yticks(range(1, len(groups) + 1))
    ax.set_yticklabels([f"{LABEL[g]}\n(n={len(v)})" for g, v in zip(groups, data, strict=False)],
                       fontsize=8)
    ax.set_xlabel("measured ΔPSI when NMD is blocked")
    ax.set_title("A  Predicted consequence vs. measured NMD sensitivity", loc="left", fontsize=11)


def panel_b(ax, d: pd.DataFrame) -> None:
    no_ptc = ["in_frame_insertion", "utr_insertion", "non_coding_host"]
    sub = d[d.consequence_class.isin(["ptc_nmd"] + no_ptc)]
    y = sub.consequence_class.eq("ptc_nmd").to_numpy()
    thresholds = np.arange(0, 26, 1.0)
    sens, spec = [], []
    for t in thresholds:
        truth = (sub.nmd_gain > t).to_numpy()
        sens.append((y & truth).sum() / max(truth.sum(), 1))
        spec.append((~y & ~truth).sum() / max((~truth).sum(), 1))
    ax.plot(thresholds, sens, "-", color=PREDICTED_NMD, lw=2, label="sensitivity")
    ax.plot(thresholds, spec, "-", color=PREDICTED_NOT, lw=2, label="specificity")
    ax.set_ylim(0, 1)
    ax.set_xlabel("ΔPSI threshold defining “NMD-sensitive”")
    ax.set_ylabel("fraction")
    ax.set_title("B  A high-recall screen, not a filter", loc="left", fontsize=11)
    ax.legend(frameon=False, fontsize=8, loc="center right")
    ax.annotate("the call catches most NMD-sensitive exons\nbut also flags many that are not",
                xy=(13, 0.55), fontsize=7.5, color="0.35", ha="center")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--table", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    args = p.parse_args()

    d = pd.read_csv(args.table, sep="\t").dropna(subset=["nmd_gain"])
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
    panel_a(axes[0], d)
    panel_b(axes[1], d)
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(args.out, dpi=200, bbox_inches="tight")

    a = d.loc[d.consequence_class == "ptc_nmd", "nmd_gain"]
    b = d.loc[d.consequence_class.isin(["in_frame_insertion", "utr_insertion", "non_coding_host"]),
              "nmd_gain"]
    _, pval = stats.mannwhitneyu(a, b, alternative="greater")
    print(f"wrote {args.out}")
    print(f"  ptc_nmd n={len(a)} median={a.median():.2f}")
    print(f"  no-PTC  n={len(b)} median={b.median():.2f}")
    print(f"  Mann-Whitney one-sided p = {pval:.3g}")


if __name__ == "__main__":
    main()
