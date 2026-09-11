#!/usr/bin/env python3
"""How badly gene-level ORA misreads splicing hits, and how much the weighting helps.

A gene contributes as many chances of being a hit as it has tested junctions, so a long,
many-exon gene is over-represented among differentially spliced genes for reasons that
are not biological. The hypergeometric null assumes the opposite — that every gene was one
equally likely draw.

This measures the cost. Each simulated universe gives every gene a number of tested units
and then draws hits *purely* in proportion to that number: there is no biology anywhere in
the data. Gene sets are then built from the largest genes, so every one of them is
biologically null and every rejection is a false positive.

The numbers quoted in docs/METHODS.md §7b come from this script.

    python validation/ora_bias_calibration.py
"""

from __future__ import annotations

import argparse

import numpy as np

from splicescope.enrich import over_representation

N_GENES = 2000
N_SETS = 40
SET_SIZE = 150


def universe(seed: int) -> tuple[list[str], list[str], dict[str, float], list[str]]:
    """Genes, hits drawn only in proportion to opportunity, weights, and a size ordering."""
    rng = np.random.default_rng(seed)
    genes = [f"G{i:04d}" for i in range(N_GENES)]
    units = rng.integers(1, 60, size=N_GENES)
    hit_probability = units / units.max() * 0.5
    hits = [g for g, p in zip(genes, hit_probability, strict=True) if rng.random() < p]
    weights = dict(zip(genes, units.astype(float), strict=True))
    by_size = [g for _, g in sorted(zip(units, genes, strict=True), reverse=True)]
    return genes, hits, weights, by_size


def null_gene_sets(rng: np.random.Generator, by_size: list[str]) -> dict[str, list[str]]:
    """Sets drawn from progressively wider slices of the largest genes."""
    return {
        f"S{i}": list(rng.choice(by_size[: 200 + i * 20], size=SET_SIZE, replace=False))
        for i in range(N_SETS)
    }


def report(name: str, pvalues: np.ndarray) -> None:
    print(
        f"  {name:24s} p<=0.05: {(pvalues <= 0.05).mean():6.1%}   "
        f"p<=1e-6: {(pvalues <= 1e-6).mean():6.1%}   median p: {np.median(pvalues):.3f}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--universes", type=int, default=12)
    args = parser.parse_args()

    plain, weighted = [], []
    for seed in range(args.universes):
        genes, hits, weights, by_size = universe(seed)
        gene_sets = null_gene_sets(np.random.default_rng(seed), by_size)
        plain.append(over_representation(hits, genes, gene_sets)["pvalue"].to_numpy())
        weighted.append(
            over_representation(hits, genes, gene_sets, weights=weights)["pvalue"].to_numpy()
        )
    plain, weighted = np.concatenate(plain), np.concatenate(weighted)

    print(
        f"{len(plain)} biologically null gene sets over {args.universes} universes of "
        f"{N_GENES} genes.\nEvery rejection below is a false positive.\n"
    )
    report("plain hypergeometric", plain)
    report("opportunity-weighted", weighted)
    print(
        "\n  The weighting is not free of assumptions — the propensity is estimated from\n"
        "  the very hits being scored — but the uncorrected test is unusable here, and\n"
        "  the corrected one sits near its nominal rate."
    )


if __name__ == "__main__":
    main()
