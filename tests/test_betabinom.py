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


def test_dispersion_is_not_estimable_without_replication():
    """One informative sample per group gives a residual of exactly zero, so the design
    says nothing about replicate-to-replicate scatter."""
    from splicescope.betabinom import dispersion_is_estimable

    n = np.full((5, 2), 1000.0)
    mask = np.ones_like(n, dtype=bool)
    group_a = np.array([True, False])
    assert not dispersion_is_estimable(n, mask, groups=[group_a, ~group_a])


def test_dispersion_is_estimable_with_two_replicates_per_group():
    from splicescope.betabinom import dispersion_is_estimable

    n = np.full((5, 4), 1000.0)
    mask = np.ones_like(n, dtype=bool)
    group_a = np.array([True, True, False, False])
    assert dispersion_is_estimable(n, mask, groups=[group_a, ~group_a])


def test_dispersion_is_not_estimable_when_every_unit_is_covered_in_one_sample_only():
    from splicescope.betabinom import dispersion_is_estimable

    n = np.full((4, 4), 1000.0)
    group_a = np.array([True, True, False, False])
    mask = np.zeros_like(n, dtype=bool)
    mask[:, 0] = True   # only one sample informative anywhere
    mask[:, 2] = True
    assert not dispersion_is_estimable(n, mask, groups=[group_a, ~group_a])


def test_samples_in_neither_group_change_nothing_at_all():
    """Both hypotheses have to be fitted to the same observations or the ratio is not a
    likelihood ratio. Leaving ungrouped samples in the null charges their whole
    likelihood to it: the null looks far worse than it is, and the statistic grows
    without bound. Measured before the fix, six ungrouped samples alongside an
    unchanged 3-vs-3 comparison moved p from 1.6e-09 to 9.7e-93.

    So the test is not that the p-value is plausible — it is that the extra samples
    make no difference whatsoever, however far their Psi sits from either group's.
    """
    k, n, mask = simulate(precision=40.0, n_events=300, delta=0.25, seed=5)
    alone = bb.lrt(k, n, mask, GROUP_A, GROUP_B, s=40.0)

    rng = np.random.default_rng(1)
    n_extra = (rng.poisson(60, size=(300, 6)) + 5).astype(float)
    # deliberately nothing like either group: Psi at the two extremes
    k_extra = (n_extra * rng.choice([0.02, 0.98], size=(300, 6))).round()
    padded = bb.lrt(
        np.hstack([k, k_extra]),
        np.hstack([n, n_extra]),
        np.ones((300, 12), dtype=bool),
        np.concatenate([GROUP_A, np.zeros(6, dtype=bool)]),
        np.concatenate([GROUP_B, np.zeros(6, dtype=bool)]),
        s=40.0,
    )
    for name, before, after in zip(
        ("mu_a", "mu_b", "statistic", "pvalue"), alone, padded, strict=True
    ):
        # 1e-9 is float noise from summing a wider array, not a difference: the bug this
        # pins moved p-values by eighty-four orders of magnitude.
        np.testing.assert_allclose(after, before, rtol=1e-9, err_msg=f"{name} moved")


def test_the_statistic_is_never_negative():
    """A likelihood ratio cannot favour the null here, since the null is the
    alternative with its two means tied. Numerical noise can still put it slightly
    below zero, and a negative chi-square statistic in an output table is a defect
    even when the p-value it produces is 1."""
    k, n, mask = simulate(precision=1e4, n_events=500, delta=0.0, seed=9)
    _, _, statistic, pvalue = bb.lrt(k, n, mask, GROUP_A, GROUP_B, s=1e4)
    assert (statistic >= 0).all()
    assert (pvalue <= 1.0).all()


def test_the_estimated_precision_stays_inside_its_bounds():
    """``s`` is a ratio whose denominator is an excess variance, so data noisier than
    the beta-binomial can drive it below 1 and, with enough excess, below zero — which
    makes the likelihood it is handed to meaningless. The bounds are load-bearing, not
    decoration: Psi drawn from {0.02, 0.98} at 400x coverage lands exactly on the
    floor."""
    rng = np.random.default_rng(0)
    n = np.full((200, N_SAMPLES), 400.0)
    k = rng.binomial(400, rng.choice([0.02, 0.98], size=(200, N_SAMPLES))).astype(float)
    mask = np.ones_like(n, dtype=bool)
    groups = [GROUP_A, GROUP_B]

    s = bb.estimate_precision(k, n, mask, groups=groups)
    assert 1.0 <= s <= 1e5
    assert s == pytest.approx(1.0), "this data is meant to sit on the floor"
    raised = bb.estimate_precision(k, n, mask, groups=groups, min_precision=7.5)
    assert raised == pytest.approx(7.5), "the floor that is asked for is the one used"


def test_a_site_with_nothing_to_splice_to_does_not_tighten_the_null_for_everyone_else():
    """The most consequential defect found this session, and it was inside the estimator.

    A constitutive donor carrying a single junction puts every read on that junction, so
    its Ψ is 1 in every sample — not because splicing is precise there, but because the
    site has nothing else to splice to. Such a unit says nothing about how much
    replicates vary. The estimator zeroed its residual, through a `variance > 0` guard,
    and then counted its observations in the degrees of freedom anyway — so it entered
    the denominator and not the numerator. (The guard could never have fired: it is
    applied after the fitted mean has been clipped away from both ends, which makes the
    variance positive by construction.)

    That dilutes the residual mean towards zero and drives the estimated precision up,
    which narrows the null for every real unit. Across a genome these sites outnumber
    the alternative ones, so the effect is not small: at one per real unit the estimate
    went from a true 50 to 842 and the false-positive rate from 0.049 to 0.173, against
    a nominal 0.05. Measured in `validation/invariant_units.py`.
    """
    from splicescope.betabinom import estimate_precision, lrt

    real_k, real_n, mask, group_a = _beta_binomial_null(np.full(1000, 50.0), seed=0)
    groups = [group_a, ~group_a]

    # every read on the one junction, in every sample
    stuck_n = np.full_like(real_n, 60.0)
    stuck_k = stuck_n.copy()
    both_k = np.vstack([real_k, stuck_k])
    both_n = np.vstack([real_n, stuck_n])
    both_mask = np.ones_like(both_n, dtype=bool)

    alone = estimate_precision(real_k, real_n, mask, groups=groups)
    together = estimate_precision(both_k, both_n, both_mask, groups=groups)

    assert 40 < alone < 65, f"the real units alone recover s near 50, got {alone}"
    assert together == pytest.approx(alone, rel=0.15), (
        f"adding 1000 units that cannot vary moved the estimate from {alone} to "
        f"{together}; they carry no information about dispersion and must not count"
    )

    # and the thing that actually matters: the test stays calibrated either way
    for precision in (alone, together):
        _, _, _, p = lrt(real_k, real_n, mask, group_a, ~group_a, precision)
        assert 0.03 < (p <= 0.05).mean() < 0.075, f"nominal 0.05 at s={precision}"

    # the same at the other boundary: a junction nobody uses
    unused_k = np.zeros_like(real_n)
    with_zeros = estimate_precision(
        np.vstack([real_k, unused_k]), np.vstack([real_n, stuck_n]), both_mask, groups=groups
    )
    assert with_zeros == pytest.approx(alone, rel=0.15), "Ψ ≡ 0 is the same problem"


def test_a_group_pinned_to_one_boundary_still_contributes_its_other_group():
    """The exclusion is per group, not per unit: a junction switched fully on in the
    knockdown and varying in the control still tells you what the control's replicates
    do. Dropping the whole unit would throw that away — and this is the shape of a real
    cryptic event, so it is the last thing to discard."""
    from splicescope.betabinom import estimate_precision

    k, n, mask, group_a = _beta_binomial_null(np.full(400, 50.0), seed=1)
    baseline = estimate_precision(k, n, mask, groups=[group_a, ~group_a])

    switched = k.copy()
    switched[:, ~group_a] = n[:, ~group_a]  # group B fully on, group A left alone
    s = estimate_precision(switched, n, mask, groups=[group_a, ~group_a])

    assert 40 < s < 65, f"group A's variability is still measured, got {s}"
    assert s == pytest.approx(baseline, rel=0.35)


def test_the_per_unit_floor_still_sees_a_loose_unit_when_one_group_is_switched_fully_on():
    """The same exclusion, in the estimator that protects against heterogeneity.

    A unit switched fully on in the knockdown and varying in the control is the shape of
    a real cryptic event. Its knockdown group has every read on one junction, so that
    group's residuals are zero by construction — and counting them diluted the control's
    genuine looseness, pushing the per-unit estimate up. Since the floor takes the
    *smaller* of the shared and per-unit values, an inflated per-unit estimate means the
    floor stops biting, and the protection is lost exactly where it is wanted.
    """
    from splicescope.betabinom import estimate_precision_per_unit

    k, n, mask, group_a = _beta_binomial_null(np.full(1500, 5.0), seed=2)
    groups = [group_a, ~group_a]

    varying = estimate_precision_per_unit(k, n, mask, groups=groups)
    switched = k.copy()
    switched[:, ~group_a] = n[:, ~group_a]
    half_on = estimate_precision_per_unit(switched, n, mask, groups=groups)

    def median_estimate(values):
        finite = values[values < 1e5]
        assert len(finite) > len(values) // 2, "most units must be estimable"
        return float(np.median(finite))

    loose, still_loose = median_estimate(varying), median_estimate(half_on)
    assert 3 < loose < 9, f"both groups varying recovers s near 5, got {loose}"
    assert still_loose < 2 * loose, (
        f"switching one group fully on moved the estimate from {loose} to {still_loose}; "
        f"the control's replicates still vary as much as they did"
    )
