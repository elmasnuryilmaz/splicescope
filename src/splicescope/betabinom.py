"""Beta-binomial machinery for count-based differential splicing.

Why not a rank test
-------------------
Comparing per-sample Ψ with Mann-Whitney throws away the read counts that Ψ was
computed from, and a rank test on ``n`` replicates per group cannot return a
p-value below ``2 / C(2n, n)`` — 0.1 for the 3-vs-3 designs that dominate
RNA-seq. Nothing survives genome-wide multiple testing, however large the
effect. Modelling the counts removes that floor: evidence grows with coverage,
not only with replicate count.

The model
---------
Inclusion reads ``k`` out of ``n`` informative reads are beta-binomial with mean
Ψ and precision ``s = α + β``, so that ``α = Ψs`` and ``β = (1 - Ψ)s``. The beta
layer absorbs the replicate-to-replicate variability a plain binomial would
mistake for signal. ``s`` is shared across events and estimated once by profile
likelihood; per event we then fit Ψ under

* ``H0`` — one Ψ for every sample, and
* ``H1`` — a separate Ψ per group,

and compare ``2(ℓ₁ − ℓ₀)`` to a chi-square with one degree of freedom.

Everything is vectorised over events: the fits are bisections on a monotone
score function, so a few hundred thousand events cost a handful of seconds.
"""

from __future__ import annotations

import numpy as np
from scipy import stats
from scipy.special import digamma, gammaln

#: Ψ is fitted inside this open interval to keep the likelihood finite.
_EPS = 1e-6
#: Bisection steps; 60 takes the bracket below double precision.
_BISECT_STEPS = 60


def loglik(k: np.ndarray, n: np.ndarray, mask: np.ndarray, mu: np.ndarray, s: float) -> np.ndarray:
    """Per-event beta-binomial log-likelihood, up to a constant.

    The binomial coefficient is dropped: it does not depend on ``mu`` and
    cancels in the likelihood ratio.
    """
    mu = np.clip(mu, _EPS, 1.0 - _EPS)[:, None]
    a, b = mu * s, (1.0 - mu) * s
    terms = (
        gammaln(k + a)
        + gammaln(n - k + b)
        - gammaln(n + s)
        - gammaln(a)
        - gammaln(b)
        + gammaln(s)
    )
    return np.where(mask, terms, 0.0).sum(axis=1)


def _score(k: np.ndarray, n: np.ndarray, mask: np.ndarray, mu: np.ndarray, s: float) -> np.ndarray:
    """Derivative of the log-likelihood with respect to ``mu`` (drops the factor ``s``)."""
    mu = np.clip(mu, _EPS, 1.0 - _EPS)[:, None]
    a, b = mu * s, (1.0 - mu) * s
    terms = digamma(k + a) - digamma(a) - digamma(n - k + b) + digamma(b)
    return np.where(mask, terms, 0.0).sum(axis=1)


def fit_mu(k: np.ndarray, n: np.ndarray, mask: np.ndarray, s: float) -> np.ndarray:
    """Maximum-likelihood Ψ per event, by vectorised bisection on the score.

    The score is decreasing in ``mu``, so a sign change brackets the root. Events
    whose score never changes sign are at a boundary (all reads included or all
    skipped) and are returned as 0 or 1.
    """
    n_events = k.shape[0]
    lo = np.full(n_events, _EPS)
    hi = np.full(n_events, 1.0 - _EPS)

    score_lo = _score(k, n, mask, lo, s)
    score_hi = _score(k, n, mask, hi, s)
    interior = (score_lo > 0) & (score_hi < 0)

    for _ in range(_BISECT_STEPS):
        mid = 0.5 * (lo + hi)
        positive = _score(k, n, mask, mid, s) > 0
        lo = np.where(positive, mid, lo)
        hi = np.where(positive, hi, mid)

    mu = 0.5 * (lo + hi)
    # Boundary cases: no interior optimum, so take the edge the score points to.
    mu = np.where(interior, mu, np.where(score_lo <= 0, 0.0, 1.0))
    return mu


def dispersion_is_estimable(
    n: np.ndarray, mask: np.ndarray, groups: list[np.ndarray] | None = None
) -> bool:
    """Whether the design carries any information about replicate-to-replicate scatter.

    Dispersion is measured by how far replicates fall from their own group's mean, so a
    group holding a single informative sample contributes a residual of exactly zero and
    says nothing. When no unit anywhere has more informative samples than fitted group
    means there is no residual degree of freedom at all, and :func:`estimate_precision`
    falls back to ``max_precision`` — which asserts *no overdispersion* rather than
    *unknown dispersion*, narrowing the test to a plain binomial one.

    That failure is quiet and severe: a 1-vs-1 comparison of Ψ 0.300 against 0.360 at
    1000 reads returns ``q = 4.5e-03`` on no replication whatsoever. Callers should check
    this before trusting a p-value.
    """
    blocks = groups if groups else [np.ones(n.shape[1], dtype=bool)]
    covered = np.zeros_like(mask)
    n_params = np.zeros(n.shape[0])
    for block in blocks:
        block_mask = mask & block[None, :]
        covered |= block_mask
        n_params += (np.where(block_mask, n, 0.0).sum(axis=1) > 0).astype(float)

    usable = (covered.sum(axis=1) - n_params) > 0
    if not usable.any():
        return False
    keep = usable[:, None] & covered
    return float(keep.sum()) - float(n_params[usable].sum()) > 0


def estimate_precision(
    k: np.ndarray,
    n: np.ndarray,
    mask: np.ndarray,
    groups: list[np.ndarray] | None = None,
    min_precision: float = 1.0,
    max_precision: float = 1e5,
) -> float:
    """Shared precision ``s``, from df-corrected Pearson residuals.

    Profile likelihood is not used here: fitting one Ψ per event on a handful of
    replicates absorbs part of the very variance we are trying to measure, which
    biases ``s`` upward and makes the test anti-conservative. The moment
    estimator below corrects for the fitted means through the residual degrees
    of freedom instead.

    Under the model ``E[r²] = 1 + (n - 1) / (s + 1)`` for Pearson residuals
    ``r``, so a df-corrected residual sum determines ``s`` in one pass.

    A single shared ``s`` is a deliberate simplification, and it is **not uniformly
    conservative**. When dispersion is homogeneous it is well calibrated (simulated at
    s = 50: estimate 52.3, false-positive rate 0.051 against a nominal 0.05). When
    dispersion is heterogeneous the shared value sits between the extremes and the test
    is anti-conservative for everything noisier than it: with half the units at s = 200
    and half at s = 5 the shared estimate is 10.7, and the loosely dispersed half runs
    at a **false-positive rate of 0.173** while the tight half runs at 0.000. The overall
    rate of 0.087 hides both. :func:`estimate_precision_per_unit` and
    ``differential_splicing(dispersion="per_unit_floor")`` trade power for protection
    against this; see ``docs/METHODS.md`` §5.4.

    ``groups`` must be the same sample masks the test will use. Residuals are
    then taken about each group's own mean: pooling groups that genuinely differ
    would charge that difference to dispersion and destroy the power the test is
    supposed to have.

    A group whose fitted mean sits exactly at 0 or 1 is excluded — from the residual sum
    and from the degrees of freedom both. It has every read on one side, so its residuals
    are zero by construction and it carries no information about replicate-to-replicate
    variability. Most such groups are structural rather than biological: a constitutive
    donor with a single junction has Ψ = 1 in every sample because the site has nothing
    else to splice to, and across a genome those sites outnumber the alternative ones.
    Counting them on one side only drove the estimate from a true 50 to 842 and the
    false-positive rate from 0.049 to 0.173 — see ``validation/invariant_units.py`` and
    ``docs/METHODS.md`` §5.5. The exclusion is per group, so a junction switched fully on
    in one condition still contributes what the other condition's replicates do, which is
    the shape of a real cryptic event.

    When the design cannot support an estimate at all this returns ``max_precision``,
    which asserts no overdispersion rather than admitting ignorance — check
    :func:`dispersion_is_estimable` first. ``differential_splicing`` does.
    """
    blocks = groups if groups else [np.ones(k.shape[1], dtype=bool)]

    mu = np.zeros_like(k)
    covered = np.zeros_like(mask)
    n_params = np.zeros(k.shape[0])
    for block in blocks:
        block_mask = mask & block[None, :]
        total = np.where(block_mask, n, 0.0).sum(axis=1)
        included = np.where(block_mask, k, 0.0).sum(axis=1)
        fitted = np.zeros_like(total)
        np.divide(included, total, out=fitted, where=total > 0)
        mu = np.where(block_mask, fitted[:, None], mu)
        covered |= block_mask
        n_params += (total > 0).astype(float)

    # A group whose fitted mean is exactly 0 or 1 put every read on one side, so its
    # residuals are zero by construction and it says nothing about how much replicates
    # vary. Most of them are structural: a constitutive donor carrying a single junction
    # has Psi = 1 in every sample because the site has nothing else to splice to, and
    # such sites outnumber the alternative ones across a genome. Counting them in the
    # degrees of freedom but not in the residual sum -- which is what happened, because
    # the `variance > 0` guard below fires after `mu` has been clipped away from both
    # ends and so never fires at all -- dilutes the residual mean towards zero and drives
    # the estimated precision up. At one such unit per real one the estimate went from a
    # true 50 to 842, and the false-positive rate from 0.049 to 0.173 against a nominal
    # 0.05. They are excluded from both sides here; see validation/invariant_units.py.
    informative = covered & (mu > 0.0) & (mu < 1.0)
    n_params = np.zeros(k.shape[0])
    for block in blocks:
        n_params += (informative & block[None, :]).any(axis=1).astype(float)

    usable = (informative.sum(axis=1) - n_params) > 0
    if not usable.any():
        return max_precision

    mu = np.clip(mu, _EPS, 1.0 - _EPS)
    variance = n * mu * (1.0 - mu)
    resid_sq = np.where(
        informative & (variance > 0), (k - n * mu) ** 2 / np.maximum(variance, _EPS), 0.0
    )

    keep = usable[:, None] & informative
    pearson = float(resid_sq[keep].sum())
    n_total = float(keep.sum())
    df = n_total - float(n_params[usable].sum())
    if df <= 0 or pearson <= 0:
        return max_precision

    corrected = pearson * n_total / df
    excess = corrected - n_total
    if excess <= 0:
        return max_precision  # no overdispersion beyond binomial

    trials_excess = float(np.where(keep, n - 1.0, 0.0).sum())
    s = trials_excess / excess - 1.0
    return float(np.clip(s, min_precision, max_precision))


def estimate_precision_per_unit(
    k: np.ndarray,
    n: np.ndarray,
    mask: np.ndarray,
    groups: list[np.ndarray] | None = None,
    min_precision: float = 1.0,
    max_precision: float = 1e5,
) -> np.ndarray:
    """The same moment estimator applied to each unit on its own.

    One unit's handful of replicates gives a very noisy estimate, so this is not a
    replacement for the shared value — it is the input to a floor. A unit that looks
    *more* dispersed than the shared estimate says is the case the shared value gets
    wrong, and taking the smaller of the two can only widen the null.

    Units that cannot be estimated return ``max_precision``, which leaves the shared
    value in charge of them.
    """
    blocks = groups if groups else [np.ones(k.shape[1], dtype=bool)]

    mu = np.zeros_like(k)
    covered = np.zeros_like(mask)
    n_params = np.zeros(k.shape[0])
    for block in blocks:
        block_mask = mask & block[None, :]
        total = np.where(block_mask, n, 0.0).sum(axis=1)
        included = np.where(block_mask, k, 0.0).sum(axis=1)
        fitted = np.zeros_like(total)
        np.divide(included, total, out=fitted, where=total > 0)
        mu = np.where(block_mask, fitted[:, None], mu)
        covered |= block_mask
        n_params += (total > 0).astype(float)

    mu = np.clip(mu, _EPS, 1.0 - _EPS)
    variance = n * mu * (1.0 - mu)
    resid_sq = np.where(
        covered & (variance > 0), (k - n * mu) ** 2 / np.maximum(variance, _EPS), 0.0
    )

    pearson = resid_sq.sum(axis=1)
    n_total = covered.sum(axis=1).astype(float)
    df = n_total - n_params

    out = np.full(k.shape[0], max_precision)
    ok = (df > 0) & (pearson > 0) & (n_total > 0)
    corrected = np.zeros_like(pearson)
    np.divide(pearson * n_total, df, out=corrected, where=ok)
    excess = corrected - n_total
    ok &= excess > 0

    trials_excess = np.where(covered, n - 1.0, 0.0).sum(axis=1)
    estimate = np.zeros_like(excess)
    np.divide(trials_excess, excess, out=estimate, where=ok)
    out = np.where(ok, np.clip(estimate - 1.0, min_precision, max_precision), out)
    return out


def lrt(
    k: np.ndarray,
    n: np.ndarray,
    mask: np.ndarray,
    group_a: np.ndarray,
    group_b: np.ndarray,
    s: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Likelihood-ratio test of "one Ψ" against "one Ψ per group".

    ``group_a`` and ``group_b`` are boolean masks over the sample axis. Returns
    ``(mu_a, mu_b, statistic, pvalue)``; the statistic is compared to a
    chi-square with one degree of freedom.

    Samples in neither group take no part. Both hypotheses must be fitted to the
    same observations or the ratio is not a likelihood ratio: leaving them in the
    null while the alternative sees only the two groups charges their entire
    likelihood to the null and inflates the statistic without bound — with six
    ungrouped samples alongside a 3-vs-3 comparison, p went from 1.6e-09 to 9.7e-93
    on unchanged group data.
    """
    mask = mask & (group_a | group_b)[None, :]
    mask_a = mask & group_a[None, :]
    mask_b = mask & group_b[None, :]

    mu_null = fit_mu(k, n, mask, s)
    mu_a = fit_mu(k, n, mask_a, s)
    mu_b = fit_mu(k, n, mask_b, s)

    ll_null = loglik(k, n, mask, mu_null, s)
    ll_alt = loglik(k, n, mask_a, mu_a, s) + loglik(k, n, mask_b, mu_b, s)

    statistic = np.maximum(2.0 * (ll_alt - ll_null), 0.0)
    pvalue = stats.chi2.sf(statistic, df=1)
    return mu_a, mu_b, statistic, pvalue


def min_achievable_rank_pvalue(n_a: int, n_b: int) -> float:
    """Smallest two-sided Mann-Whitney p-value possible for these group sizes.

    Provided so the limitation this module exists to fix can be shown rather
    than asserted: with 3 against 3 the floor is 0.1, so no rank test can clear
    genome-wide correction no matter how large the effect.

    This is the *exact* test's floor, and not what ``scipy.stats.mannwhitneyu`` returns
    once there are ties: it then uses the normal approximation, which gives 0.081 for an
    untied 3 against 3, 0.064 with one group tied, and 0.047 when Psi is 0 in every
    control and 1 in every knockdown — the signature this module exists to find, under
    the nominal 0.05. What stops the rank test calling it is Benjamini-Hochberg across a
    genome's worth of junctions, not the floor. See METHODS section 5.1.
    """
    from math import comb

    return 2.0 / comb(n_a + n_b, n_a)
