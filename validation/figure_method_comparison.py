#!/usr/bin/env python3
"""Figure: two independent NMD predictions against the same measured ground truth.

Neither method's operating characteristics had been published. Scoring both on
the exons they can each classify shows they separate the ground truth to a
similar degree, and that the difference between them is one of resolution --
LeafCutter2 declines to decide on a fifth of the events -- rather than accuracy.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

SPLICESCOPE = "#1b6ca8"
LEAFCUTTER2 = "#c9772f"
UNDECIDED = "#9a9a9a"


def panel_a(ax, d: pd.DataFrame) -> None:
    """Measured NMD gain, split by each method's call."""
    groups = [
        ("splicescope\nPTC → NMD", d[d.consequence_class.eq("ptc_nmd")].nmd_gain, SPLICESCOPE),
        ("splicescope\nother", d[~d.consequence_class.eq("ptc_nmd")].nmd_gain, SPLICESCOPE),
        ("LeafCutter2\nunproductive", d[d.lc2.eq("UP")].nmd_gain, LEAFCUTTER2),
        ("LeafCutter2\nproductive", d[d.lc2.eq("PR")].nmd_gain, LEAFCUTTER2),
        ("LeafCutter2\nundecided", d[d.lc2.eq("NE")].nmd_gain, UNDECIDED),
    ]
    data = [g[1].to_numpy() for g in groups]
    bp = ax.boxplot(data, vert=False, widths=0.6, patch_artist=True, showfliers=False)
    for patch, (_, _, colour) in zip(bp["boxes"], groups, strict=False):
        patch.set_facecolor(colour)
        patch.set_alpha(0.3)
        patch.set_edgecolor(colour)
    for element in ("medians", "whiskers", "caps"):
        for line in bp[element]:
            line.set_color("0.3")
    rng = np.random.default_rng(0)
    for i, (values, (_, _, colour)) in enumerate(zip(data, groups, strict=False), start=1):
        ax.scatter(values, i + rng.uniform(-0.15, 0.15, len(values)), s=6, alpha=0.4,
                   color=colour, linewidths=0)
    ax.axvline(0, color="0.6", lw=1, ls="--")
    ax.set_yticks(range(1, len(groups) + 1))
    ax.set_yticklabels([f"{g[0]}\n(n={len(g[1])})" for g in groups], fontsize=7.5)
    ax.set_xlabel("measured ΔPSI when NMD is blocked")
    ax.set_title("A  Both methods separate the ground truth", loc="left", fontsize=11)


def panel_b(ax, d: pd.DataFrame) -> None:
    """Discrimination, with bootstrap intervals, on each method's own decided set."""
    no_ptc = ["in_frame_insertion", "utr_insertion", "non_coding_host"]
    ss = d[d.consequence_class.isin(["ptc_nmd"] + no_ptc)]
    lc = d[d.lc2.isin(["UP", "PR"])]
    entries = [
        ("splicescope", ss.consequence_class.eq("ptc_nmd").to_numpy(), ss.nmd_gain.to_numpy(),
         SPLICESCOPE),
        ("LeafCutter2", lc.lc2.eq("UP").to_numpy(), lc.nmd_gain.to_numpy(), LEAFCUTTER2),
    ]
    rng = np.random.default_rng(0)
    for i, (_name, y, score, colour) in enumerate(entries):
        auc = roc_auc_score(y, score)
        boot = []
        for _ in range(2000):
            idx = rng.choice(len(y), len(y), replace=True)
            if y[idx].std() > 0:
                boot.append(roc_auc_score(y[idx], score[idx]))
        low, high = np.percentile(boot, [2.5, 97.5])
        ax.errorbar(auc, i, xerr=[[auc - low], [high - auc]], fmt="o", ms=8,
                    color=colour, capsize=4, lw=2)
        ax.text(auc, i + 0.18, f"{auc:.3f}", ha="center", fontsize=8.5, color=colour)
    ax.axvline(0.5, color="0.5", ls="--", lw=1.2)
    ax.text(0.505, -0.45, "chance", fontsize=7.5, color="0.45", ha="left")
    ax.set_yticks(range(len(entries)))
    ax.set_yticklabels([e[0] for e in entries], fontsize=9)
    ax.set_ylim(-0.6, 1.6)
    ax.set_xlim(0.45, 0.85)
    ax.set_xlabel("AUROC against measured NMD sensitivity (95% CI)")
    ax.set_title("B  Comparable, and both modest", loc="left", fontsize=11)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--table", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    args = p.parse_args()
    d = pd.read_csv(args.table, sep="\t").dropna(subset=["nmd_gain", "lc2"])

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
    panel_a(axes[0], d)
    panel_b(axes[1], d)
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(args.out, dpi=200, bbox_inches="tight")
    print(f"wrote {args.out}  (n={len(d)} exons)")


if __name__ == "__main__":
    main()
