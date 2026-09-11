"""Beta-binomial estimation, calibration and the rank-test floor it replaces."""

from __future__ import annotations

import numpy as np
import pytest

from splicescope import betabinom as bb
from splicescope.diff import benjamini_hochberg

N_SAMPLES = 6
GROUP_A = np.array([True] * 3 + [False] * 3)
GROUP_B = ~GROUP_A


def simulate(precision, n_events=2000, delta=0.0, depth=60, seed=0):
    """Beta-binomial counts with genuine per-replicate variability."""
    rng = np.random.default_rng(seed)
    n = (rng.poisson(depth, size=(n_events, N_SAMPLES)) + 5).astype(float)
    mu_a = rng.uniform(0.2, 0.6, n_events)
    mu = np.where(GROUP_A[None, :], mu_a[:, None], np.clip(mu_a + delta, 0.01, 0.99)[:, None])
    p = rng.beta(mu * precision, (1.0 - mu) * precision)
    k = rng.binomial(n.astype(int), p).astype(float)
    return k, n, np.ones_like(n, dtype=bool)


def test_rank_test_floor_is_what_motivates_this_module():
    # No 3-vs-3 rank test can drop below 0.1, so nothing survives BH genome-wide.
    assert bb.min_achievable_rank_pvalue(3, 3) == pytest.approx(0.1)
    assert bb.min_achievable_rank_pvalue(4, 4) == pytest.approx(2 / 70)
    assert bb.min_achievable_rank_pvalue(10, 10) < 1e-4


def test_fit_mu_recovers_the_simulated_psi():
    k, n, mask = simulate(precision=200.0, seed=3)
    mu = bb.fit_mu(k, n, mask, 200.0)
    observed = k.sum(axis=1) / n.sum(axis=1)
    assert np.corrcoef(mu, observed)[0, 1] > 0.99
    assert np.abs(mu - observed).mean() < 0.02


def test_fit_mu_handles_all_or_nothing_events():
    n = np.full((2, 4), 20.0)
    k = np.array([[0.0] * 4, [20.0] * 4])
    mask = np.ones_like(n, dtype=bool)
    mu = bb.fit_mu(k, n, mask, 50.0)
    assert mu[0] == pytest.approx(0.0, abs=1e-3)
    assert mu[1] == pytest.approx(1.0, abs=1e-3)


@pytest.mark.parametrize("precision", [10.0, 40.0, 200.0])
def test_precision_estimator_is_roughly_unbiased(precision):
    k, n, mask = simulate(precision, n_events=4000, seed=5)
    assert bb.estimate_precision(k, n, mask) == pytest.approx(precision, rel=0.35)


def test_pooling_differing_groups_would_inflate_dispersion():
    """Estimating about a pooled mean charges real signal to dispersion."""
    k, n, mask = simulate(200.0, n_events=3000, delta=0.35, seed=17)
    pooled = bb.estimate_precision(k, n, mask)
    by_group = bb.estimate_precision(k, n, mask, groups=[GROUP_A, GROUP_B])
    assert by_group > 5 * pooled


def test_binomial_data_reports_no_overdispersion():
    rng = np.random.default_rng(11)
    n = np.full((3000, N_SAMPLES), 50.0)
    p = rng.uniform(0.3, 0.7, 3000)[:, None]
    k = rng.binomial(50, np.repeat(p, N_SAMPLES, axis=1)).astype(float)
    assert bb.estimate_precision(k, n, np.ones_like(n, dtype=bool)) > 1e4


@pytest.mark.parametrize("precision", [20.0, 200.0])
def test_null_pvalues_are_calibrated(precision):
    k, n, mask = simulate(precision, n_events=4000, delta=0.0, seed=7)
    s = bb.estimate_precision(k, n, mask, groups=[GROUP_A, GROUP_B])
    *_, pvalue = bb.lrt(k, n, mask, GROUP_A, GROUP_B, s)
    # A miscalibrated test is worse than no test: false positives look like hits.
    assert np.mean(pvalue < 0.05) < 0.075
    assert np.mean(pvalue < 0.01) < 0.02


def test_real_difference_clears_bh_where_a_rank_test_cannot():
    k, n, mask = simulate(200.0, n_events=2000, delta=0.35, seed=13)
    s = bb.estimate_precision(k, n, mask, groups=[GROUP_A, GROUP_B])
    mu_a, mu_b, statistic, pvalue = bb.lrt(k, n, mask, GROUP_A, GROUP_B, s)
    assert np.all(statistic >= 0)
    assert (mu_b - mu_a).mean() > 0.2
    # The whole point: p-values far below the 0.1 floor of a 3-vs-3 rank test.
    assert pvalue.min() < bb.min_achievable_rank_pvalue(3, 3) / 100
    assert (benjamini_hochberg(pvalue) <= 0.05).sum() > 0


def _beta_binomial_null(s_values, mu=0.3, depth=60, n_per_group=3, seed=0):
    """Counts with no difference between groups, at a chosen precision per unit."""
    rng = np.random.default_rng(seed)
    n_samples = 2 * n_per_group
    n = np.full((len(s_values), n_samples), float(depth))
    k = np.empty_like(n)
    for i, s in enumerate(s_values):
        p = rng.beta(mu * s, (1 - mu) * s, size=n_samples)
        k[i] = rng.binomial(depth, p)
    group_a = np.array([True] * n_per_group + [False] * n_per_group)
    return k, n, np.ones_like(n, dtype=bool), group_a


def test_the_shared_precision_is_well_calibrated_when_dispersion_is_homogeneous():
    from splicescope.betabinom import estimate_precision, lrt

    k, n, mask, group_a = _beta_binomial_null(np.full(2000, 50.0), seed=0)
    s = estimate_precision(k, n, mask, groups=[group_a, ~group_a])
    assert 40 < s < 65, f"expected to recover s near 50, got {s}"
    _, _, _, p = lrt(k, n, mask, group_a, ~group_a, s)
    assert 0.03 < (p <= 0.05).mean() < 0.075


def test_the_shared_precision_is_anti_conservative_when_dispersion_is_heterogeneous():
    """Recorded, not endorsed: one precision for every unit is not uniformly
    conservative, and the overall rate hides it."""
    from splicescope.betabinom import estimate_precision, lrt

    s_true = np.array([200.0] * 1000 + [5.0] * 1000)
    k, n, mask, group_a = _beta_binomial_null(s_true, seed=0)
    s = estimate_precision(k, n, mask, groups=[group_a, ~group_a])
    _, _, _, p = lrt(k, n, mask, group_a, ~group_a, s)

    loose = s_true == 5.0
    assert (p[loose] <= 0.05).mean() > 0.12, "the loosely dispersed half should over-call"
    assert (p[~loose] <= 0.05).mean() < 0.01, "the tight half should under-call"


def test_the_per_unit_floor_never_raises_the_precision():
    """It can only widen the null, so it can only remove false positives."""
    from splicescope.betabinom import estimate_precision, estimate_precision_per_unit

    s_true = np.array([200.0] * 300 + [5.0] * 300)
    k, n, mask, group_a = _beta_binomial_null(s_true, seed=1)
    shared = estimate_precision(k, n, mask, groups=[group_a, ~group_a])
    per_unit = estimate_precision_per_unit(k, n, mask, groups=[group_a, ~group_a])
    assert (np.minimum(shared, per_unit) <= shared + 1e-9).all()


def test_the_per_unit_floor_reduces_the_heterogeneous_false_positive_rate():
    from splicescope.betabinom import estimate_precision, estimate_precision_per_unit, lrt

    s_true = np.array([200.0] * 1000 + [5.0] * 1000)
    k, n, mask, group_a = _beta_binomial_null(s_true, seed=0)
    shared = estimate_precision(k, n, mask, groups=[group_a, ~group_a])
    per_unit = estimate_precision_per_unit(k, n, mask, groups=[group_a, ~group_a])
    floored = np.minimum(shared, per_unit)[:, None]

    _, _, _, p_shared = lrt(k, n, mask, group_a, ~group_a, shared)
    _, _, _, p_floor = lrt(k, n, mask, group_a, ~group_a, floored)

    loose = s_true == 5.0
    assert (p_floor[loose] <= 0.05).mean() < (p_shared[loose] <= 0.05).mean()
    assert (p_floor <= 0.05).mean() < 0.06


def test_per_unit_precision_matches_the_shared_estimator_applied_row_by_row():
    from splicescope.betabinom import estimate_precision, estimate_precision_per_unit

    k, n, mask, group_a = _beta_binomial_null(np.full(40, 30.0), seed=2)
    blocks = [group_a, ~group_a]
    vectorised = estimate_precision_per_unit(k, n, mask, groups=blocks)
    row_by_row = np.array(
        [
            estimate_precision(k[i : i + 1], n[i : i + 1], mask[i : i + 1], groups=blocks)
            for i in range(k.shape[0])
        ]
    )
    assert np.allclose(vectorised, row_by_row)
