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
    ``r``, so a df-corrected residual sum determines ``s`` in one pass. A single
    shared ``s`` is a deliberate simplification, and a conservative one: it is
    dominated by the many low-coverage events rather than by the few loud ones.

    ``groups`` must be the same sample masks the test will use. Residuals are
    then taken about each group's own mean: pooling groups that genuinely differ
    would charge that difference to dispersion and destroy the power the test is
    supposed to have.
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

    usable = (covered.sum(axis=1) - n_params) > 0
    if not usable.any():
        return max_precision

    mu = np.clip(mu, _EPS, 1.0 - _EPS)
    variance = n * mu * (1.0 - mu)
    resid_sq = np.where(
        covered & (variance > 0), (k - n * mu) ** 2 / np.maximum(variance, _EPS), 0.0
    )

    keep = usable[:, None] & covered
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
    """
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
    """
    from math import comb

    return 2.0 / comb(n_a + n_b, n_a)
