#!/usr/bin/env python3
"""What one shared dispersion costs, and what the per-unit floor buys back.

The beta-binomial test needs a precision `s`. Estimating one for every unit from three
replicates is far too noisy, so `differential_splicing` estimates a single shared value
by default. That is well calibrated only when dispersion is actually homogeneous, and
the aggregate false-positive rate hides what happens when it is not: the shared value
lands between the extremes, so the test over-calls everything noisier than it and
under-calls everything tighter.

`dispersion="per_unit_floor"` also estimates each unit's own precision and takes the
smaller of the two. A per-unit estimate is too noisy to use directly, but a minimum can
only *widen* the null, so it can only remove false positives — at a cost in power.

This measures both, on simulated data with no differences (for the false-positive rates)
and with a real one (for the power). Every parameter is stated below rather than tuned;
the numbers in docs/METHODS.md §5.4 come from this script, and a test compares them.

    python validation/dispersion_trade.py
    python validation/dispersion_trade.py --seeds 20   # tighter Monte-Carlo error
"""

from __future__ import annotations

import argparse

import numpy as np

from splicescope.betabinom import estimate_precision, estimate_precision_per_unit, lrt

#: 3 against 3 at 60 reads is the design this whole module exists for — it is where a
#: rank test cannot return a p-value below 0.1 however large the effect.
N_PER_GROUP = 3
DEPTH = 60
#: Units per stratum. The heterogeneous null uses one stratum at each precision.
N_UNITS = 1000
#: The two precisions of the heterogeneous null, and the single one of the homogeneous.
S_TIGHT, S_LOOSE, S_HOMOGENEOUS = 200.0, 5.0, 50.0
#: The difference used for power, at the homogeneous precision.
PSI_A, PSI_B = 0.30, 0.45
ALPHA = 0.05

GROUP_A = np.array([True] * N_PER_GROUP + [False] * N_PER_GROUP)


def _counts(mu_per_sample: np.ndarray, s_values: np.ndarray, seed: int):
    """Beta-binomial counts: one row per unit, its own precision, its own per-sample Ψ."""
    rng = np.random.default_rng(seed)
    n_samples = 2 * N_PER_GROUP
    n = np.full((len(s_values), n_samples), float(DEPTH))
    k = np.empty_like(n)
    for i, s in enumerate(s_values):
        mu = mu_per_sample[i]
        k[i] = rng.binomial(DEPTH, rng.beta(mu * s, (1.0 - mu) * s))
    return k, n, np.ones_like(n, dtype=bool)


def _rates(k, n, mask, s_true=None):
    """Rejection rate under the shared precision and under the per-unit floor."""
    shared = estimate_precision(k, n, mask, groups=[GROUP_A, ~GROUP_A])
    per_unit = estimate_precision_per_unit(k, n, mask, groups=[GROUP_A, ~GROUP_A])
    floored = np.minimum(shared, per_unit)[:, None]
    _, _, _, p_shared = lrt(k, n, mask, GROUP_A, ~GROUP_A, shared)
    _, _, _, p_floor = lrt(k, n, mask, GROUP_A, ~GROUP_A, floored)

    out = {"shared estimate": shared}
    for label, p in (("shared", p_shared), ("per_unit_floor", p_floor)):
        out[f"{label} all"] = float((p <= ALPHA).mean())
        if s_true is not None:
            loose = s_true == S_LOOSE
            out[f"{label} loose"] = float((p[loose] <= ALPHA).mean())
            out[f"{label} tight"] = float((p[~loose] <= ALPHA).mean())
    return out


def measure(seeds: int = 5) -> dict[str, float]:
    """Every figure in METHODS §5.4, averaged over independent simulations."""
    flat = np.ones((N_UNITS * 2, 2 * N_PER_GROUP))
    s_hetero = np.array([S_TIGHT] * N_UNITS + [S_LOOSE] * N_UNITS)
    s_homo = np.full(N_UNITS * 2, S_HOMOGENEOUS)
    difference = np.where(GROUP_A[None, :], PSI_A, PSI_B) * np.ones((N_UNITS * 2, 1))

    runs: list[dict[str, float]] = []
    for seed in range(seeds):
        hetero = _rates(*_counts(flat * 0.3, s_hetero, seed), s_true=s_hetero)
        homo = _rates(*_counts(flat * 0.3, s_homo, 1000 + seed))
        power = _rates(*_counts(difference, s_homo, 2000 + seed))
        runs.append(
            {
                "heterogeneous shared estimate": hetero["shared estimate"],
                "heterogeneous loose, shared": hetero["shared loose"],
                "heterogeneous tight, shared": hetero["shared tight"],
                "heterogeneous all, shared": hetero["shared all"],
                "heterogeneous loose, floor": hetero["per_unit_floor loose"],
                "heterogeneous all, floor": hetero["per_unit_floor all"],
                "homogeneous, shared": homo["shared all"],
                "homogeneous, floor": homo["per_unit_floor all"],
                "power, shared": power["shared all"],
                "power, floor": power["per_unit_floor all"],
            }
        )
    return {key: float(np.mean([r[key] for r in runs])) for key in runs[0]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=int, default=5)
    args = parser.parse_args()

    values = measure(args.seeds)
    print(
        f"{N_UNITS} units per stratum, {N_PER_GROUP} against {N_PER_GROUP} at {DEPTH} reads, "
        f"averaged over {args.seeds} simulations.\n"
    )
    width = max(len(k) for k in values)
    for key, value in values.items():
        digits = 1 if "estimate" in key else 3
        print(f"  {key:<{width}}  {value:.{digits}f}")


if __name__ == "__main__":
    main()
