#!/usr/bin/env python3
"""What a splice site with nothing to choose between does to everyone else's p-value.

Most splice sites in a genome are constitutive and carry a single junction. Their Ψ is
1.0 in every sample — not because splicing is precise there, but because the site has
nothing else to splice to. `compute_psi` measures them like any other site, and
`differential_splicing` then tests them, twice over:

  * each one enters the Benjamini-Hochberg denominator, though its p-value is exactly 1
    and no threshold could ever reject it; and
  * each one enters the *shared dispersion estimate*, contributing a residual of exactly
    zero, though it holds no information about replicate-to-replicate variability at all.

The first costs power. The second is worse: it pulls the estimated precision up, which
narrows the null for every real unit and makes its p-value too small.

This measures both, on data where the truth is known. `filter_invariant=True` sets such
units aside before either step; the numbers in docs/METHODS.md §5.5 come from here and a
test compares them.

    python validation/invariant_units.py
    python validation/invariant_units.py --seeds 20   # tighter Monte-Carlo error
"""

from __future__ import annotations

import argparse

import numpy as np

from splicescope.betabinom import estimate_precision, lrt


def _estimate_counting_degenerate_groups(k, n, mask, groups):
    """The estimate as it was before groups pinned to 0 or 1 were excluded.

    Kept here, and only here, so the "before" row of the table in METHODS §5.5 is a
    measurement rather than a remembered number. It is the same moment estimator with one
    difference: an observation whose fitted group mean sits at a boundary is counted in
    the degrees of freedom, as it used to be, while still contributing a residual of zero.
    """
    mu = np.zeros_like(k)
    covered = np.zeros_like(mask)
    n_params = np.zeros(k.shape[0])
    for block in groups:
        block_mask = mask & block[None, :]
        total = np.where(block_mask, n, 0.0).sum(axis=1)
        included = np.where(block_mask, k, 0.0).sum(axis=1)
        fitted = np.zeros_like(total)
        np.divide(included, total, out=fitted, where=total > 0)
        mu = np.where(block_mask, fitted[:, None], mu)
        covered |= block_mask
        n_params += (total > 0).astype(float)

    usable = (covered.sum(axis=1) - n_params) > 0
    mu = np.clip(mu, 1e-9, 1.0 - 1e-9)
    resid_sq = (k - n * mu) ** 2 / np.maximum(n * mu * (1.0 - mu), 1e-9)
    keep = usable[:, None] & covered
    pearson = float(np.where(keep, resid_sq, 0.0).sum())
    n_total = float(keep.sum())
    df = n_total - float(n_params[usable].sum())
    excess = pearson * n_total / df - n_total
    if excess <= 0:
        return 1e5
    return float(np.clip(float(np.where(keep, n - 1.0, 0.0).sum()) / excess - 1.0, 1.0, 1e5))

#: The design the beta-binomial exists for, matching validation/dispersion_trade.py.
N_PER_GROUP = 3
DEPTH = 60
#: Units that vary. Everything below is measured on these; the invariant ones are only
#: ever present or absent.
N_REAL = 1000
#: Invariant units per real unit. A human annotation is far more lopsided than this —
#: constitutive sites outnumber alternative ones by more than an order of magnitude —
#: so 1:1 is a deliberately mild version of the problem.
INVARIANT_RATIO = 1.0
#: The precision the real units actually have.
S_TRUE = 50.0
#: Null Ψ, and the difference used for power.
PSI_NULL = 0.30
PSI_A, PSI_B = 0.30, 0.45
ALPHA = 0.05

GROUP_A = np.array([True] * N_PER_GROUP + [False] * N_PER_GROUP)


def _real_units(mu_per_sample: np.ndarray, seed: int):
    """Beta-binomial counts for units that genuinely vary."""
    rng = np.random.default_rng(seed)
    n = np.full((N_REAL, 2 * N_PER_GROUP), float(DEPTH))
    k = np.empty_like(n)
    for i in range(N_REAL):
        mu = mu_per_sample[i]
        k[i] = rng.binomial(DEPTH, rng.beta(mu * S_TRUE, (1.0 - mu) * S_TRUE))
    return k, n


def _invariant_units(count: int):
    """A constitutive site: every read on the one junction, so Ψ = 1 in every sample."""
    n = np.full((count, 2 * N_PER_GROUP), float(DEPTH))
    return n.copy(), n


def _measure(k_real, n_real, seed: int) -> dict[str, float]:
    """Rejection rates for the real units, with the invariant ones in and out."""
    n_invariant = int(round(N_REAL * INVARIANT_RATIO))
    k_inv, n_inv = _invariant_units(n_invariant)

    k_all = np.vstack([k_real, k_inv])
    n_all = np.vstack([n_real, n_inv])
    mask_all = np.ones_like(n_all, dtype=bool)
    mask_real = np.ones_like(n_real, dtype=bool)

    blocks = [GROUP_A, ~GROUP_A]
    estimates = {
        # what the estimator does now, with the invariant units present and absent
        "in": estimate_precision(k_all, n_all, mask_all, groups=blocks),
        "out": estimate_precision(k_real, n_real, mask_real, groups=blocks),
        # and what it did before they were excluded from the degrees of freedom
        "counted": _estimate_counting_degenerate_groups(k_all, n_all, mask_all, blocks),
    }
    out = {"n invariant": float(n_invariant)}
    for label, s in estimates.items():
        # every p-value set is for the *real* units; only the precision behind them differs
        _, _, _, p = lrt(k_real, n_real, mask_real, GROUP_A, ~GROUP_A, s)
        out[f"precision {label}"] = float(s)
        out[f"rate {label}"] = float((p <= ALPHA).mean())
    return out


def measure(seeds: int = 5) -> dict[str, float]:
    """Every figure in METHODS §5.5, averaged over independent simulations."""
    flat = np.ones((N_REAL, 2 * N_PER_GROUP))
    difference = np.where(GROUP_A[None, :], PSI_A, PSI_B) * np.ones((N_REAL, 1))

    runs: list[dict[str, float]] = []
    for seed in range(seeds):
        null = _measure(*_real_units(flat * PSI_NULL, seed), seed)
        power = _measure(*_real_units(difference, 2000 + seed), seed)
        n_total = N_REAL + null["n invariant"]
        runs.append(
            {
                "true precision": S_TRUE,
                "estimated precision, invariant absent": null["precision out"],
                "estimated precision, invariant present": null["precision in"],
                "estimated precision, before the fix": null["precision counted"],
                "false-positive rate, invariant absent": null["rate out"],
                "false-positive rate, invariant present": null["rate in"],
                "false-positive rate, before the fix": null["rate counted"],
                "power, invariant present": power["rate in"],
                "power, before the fix": power["rate counted"],
                # what the correction costs: the BH threshold for the best unit
                "BH denominator, invariant in": n_total,
                "BH denominator, invariant out": float(N_REAL),
            }
        )
    return {key: float(np.mean([r[key] for r in runs])) for key in runs[0]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=int, default=5)
    args = parser.parse_args()

    values = measure(args.seeds)
    print(
        f"{N_REAL} units that vary and {int(N_REAL * INVARIANT_RATIO)} that cannot, "
        f"{N_PER_GROUP} against {N_PER_GROUP} at {DEPTH} reads,\n"
        f"averaged over {args.seeds} simulations. Every rate is measured on the units "
        f"that vary.\n"
    )
    width = max(len(k) for k in values)
    for key, value in values.items():
        digits = 3 if "rate" in key or "power" in key else 1
        print(f"  {key:<{width}}  {value:.{digits}f}")


if __name__ == "__main__":
    main()
