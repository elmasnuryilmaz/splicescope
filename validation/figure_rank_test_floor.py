#!/usr/bin/env python3
"""Figure: why a rank test cannot call differential splicing genome-wide.

Panel A is arithmetic — the smallest p-value a two-sided Mann-Whitney can return
for a given group size, against the p-value the top hit needs to survive
Benjamini-Hochberg across the junctions actually tested. Panel B is the observed
consequence on real data: every p-value the rank test produced, none of them
close to the threshold, beside the p-values the count-based test produced on the
same junctions.
"""

from __future__ import annotations

import argparse
from math import comb
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

RANK = "#b23a48"
COUNT = "#1b6ca8"


def panel_a(ax, n_tests: int) -> None:
    sizes = np.arange(3, 15)
    floors = np.array([2 / comb(2 * n, n) for n in sizes])
    needed = 0.05 / n_tests
    # Where the floor finally drops under the threshold, so the claim can be read
    # off the figure rather than taken on trust.
    crossing = int(sizes[np.argmax(floors < needed)]) if (floors < needed).any() else None

    ax.semilogy(sizes, floors, "o-", color=RANK, lw=2, ms=4,
                label="smallest p a rank test can return")
    ax.axhline(needed, color="0.25", ls="--", lw=1.5,
               label=f"needed for BH q$\\leq$0.05\n({n_tests:,} tests)")
    unreachable = floors > needed
    ax.fill_between(sizes[unreachable], floors[unreachable], needed,
                    color=RANK, alpha=0.08)
    ax.set_xlabel("replicates per group")
    ax.set_ylabel("two-sided p-value")
    ax.set_title("A  No usable design clears the threshold", loc="left", fontsize=11)
    ax.set_xticks(sizes[::2])
    ax.set_ylim(needed / 30, 3)
    ax.legend(fontsize=8, loc="lower left", frameon=True, framealpha=0.92,
              edgecolor="none", facecolor="white")
    ax.annotate("3 vs 3\n$p_{min}=0.1$", xy=(3, 0.1), xytext=(4.6, 0.5),
                fontsize=8, color=RANK, ha="center",
                arrowprops=dict(arrowstyle="->", color=RANK, lw=1))
    if crossing:
        ax.annotate(f"reachable only\nfrom {crossing} vs {crossing}",
                    xy=(crossing, needed), xytext=(crossing - 1.4, needed / 12),
                    fontsize=8, color="0.25", ha="center",
                    arrowprops=dict(arrowstyle="->", color="0.35", lw=1))


def panel_b(ax, rank_p: np.ndarray, count_p: np.ndarray, n_tests: int) -> None:
    floor = 1e-30
    bins = np.logspace(np.log10(floor), 0, 60)
    needed = 0.05 / n_tests
    n_past = int((count_p <= needed).sum())
    ax.hist(np.clip(count_p, floor, 1), bins=bins, color=COUNT, alpha=0.75,
            label=f"beta-binomial ({n_past:,} past threshold)")
    ax.hist(np.clip(rank_p, floor, 1), bins=bins, color=RANK, alpha=0.75,
            label="Mann-Whitney (0 past threshold)")
    ax.axvline(needed, color="0.25", ls="--", lw=1.5)
    ax.annotate("values $<10^{-30}$\ncollected here",
                xy=(floor * 1.5, 40), fontsize=7, color="0.35")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("p-value")
    ax.set_ylabel("junctions")
    ax.set_title("B  The same junctions, both tests", loc="left", fontsize=11)
    ax.legend(frameon=False, fontsize=8, loc="upper left")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--rank", required=True, type=Path, help="differential_splicing.tsv, rank test")
    p.add_argument(
        "--count", required=True, type=Path, help="differential_splicing.tsv, count test"
    )
    p.add_argument("--out", required=True, type=Path)
    args = p.parse_args()

    rank = pd.read_csv(args.rank, sep="\t", low_memory=False)["pvalue"].dropna().to_numpy()
    count = pd.read_csv(args.count, sep="\t", low_memory=False)["pvalue"].dropna().to_numpy()

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8))
    panel_a(axes[0], len(rank))
    panel_b(axes[1], rank, count, len(rank))
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=200, bbox_inches="tight")
    print(f"wrote {args.out}")
    print(f"  rank test:  n={len(rank):,}  min p={rank.min():.4g}")
    print(f"  count test: n={len(count):,}  min p={count.min():.4g}")


if __name__ == "__main__":
    main()
