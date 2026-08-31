"""Differential splicing between two groups of samples.

Two tests are available, and the difference between them decides whether the
tool can find anything at all on a real dataset.

``betabinom`` (default where read counts are available)
    A beta-binomial likelihood-ratio test on the inclusion and total counts Ψ
    was computed from. Evidence accumulates with coverage, so a well-covered
    event can reach very small p-values even with three replicates per group.

``ranksum`` (legacy, used only when counts are absent)
    Two-sided Mann-Whitney U on the per-sample Ψ values. This discards the
    counts, and a rank test on ``n`` replicates per group cannot return a
    p-value below ``2 / C(2n, n)`` — 0.1 for a 3-vs-3 design. After
    Benjamini-Hochberg across a genome's worth of junctions nothing can ever be
    called significant, whatever the effect size. Keep it only for Ψ tables that
    carry no counts.

Both report ``delta_psi`` (group B minus group A), ``pvalue`` and BH ``qvalue``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from .betabinom import estimate_precision, lrt

#: Count column pairs tried in order when ``inc_col``/``total_col`` are not given.
_COUNT_COLUMNS = {
    "psi_donor": ("count", "donor_total"),
    "psi_acceptor": ("count", "acceptor_total"),
    "psi": ("inc_reads", "total_reads"),
    "psi_cassette": ("inc_reads", "total_reads"),
}


def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    """Return BH-adjusted q-values for a 1-D array of p-values (NaNs preserved)."""
    p = np.asarray(pvals, dtype=float)
    q = np.full_like(p, np.nan)
    mask = ~np.isnan(p)
    m = int(mask.sum())
    if m == 0:
        return q
    idx = np.where(mask)[0]
    order = idx[np.argsort(p[idx])]
    ranked = p[order] * m / (np.arange(1, m + 1))
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]  # enforce monotonicity
    q[order] = np.clip(ranked, 0, 1)
    return q


def _resolve_count_columns(
    df: pd.DataFrame, value: str, inc_col: str | None, total_col: str | None
) -> tuple[str, str] | None:
    if inc_col and total_col:
        return (inc_col, total_col) if {inc_col, total_col} <= set(df.columns) else None
    guess = _COUNT_COLUMNS.get(value)
    if guess and set(guess) <= set(df.columns):
        return guess
    return None


def _wide(df: pd.DataFrame, key: list[str], column: str, samples: list[str]) -> np.ndarray:
    table = df.pivot_table(index=key, columns="sample", values=column, aggfunc="sum")
    return table.reindex(columns=samples)


def _finalize(res: pd.DataFrame) -> pd.DataFrame:
    if res.empty:
        return res
    res["qvalue"] = benjamini_hochberg(res["pvalue"].to_numpy())
    res["abs_delta_psi"] = res["delta_psi"].abs()
    return res.sort_values(["qvalue", "abs_delta_psi"], ascending=[True, False]).reset_index(
        drop=True
    )


def _betabinom_test(
    df: pd.DataFrame,
    key: list[str],
    extra: list[str],
    value: str,
    counts: tuple[str, str],
    a_name: str,
    b_name: str,
    groups: dict[str, str],
    min_samples: int,
) -> pd.DataFrame:
    inc_col, total_col = counts
    samples = sorted(df["sample"].unique())
    psi_wide = _wide(df, key, value, samples)
    index = psi_wide.index

    valid = np.array(psi_wide.notna().to_numpy(), dtype=bool)  # writable copy
    k = np.nan_to_num(_wide(df, key, inc_col, samples).to_numpy(), nan=0.0)
    n = np.nan_to_num(_wide(df, key, total_col, samples).to_numpy(), nan=0.0)
    valid &= n > 0

    is_a = np.array([groups.get(s) == a_name for s in samples])
    is_b = np.array([groups.get(s) == b_name for s in samples])
    n_a = (valid & is_a[None, :]).sum(axis=1)
    n_b = (valid & is_b[None, :]).sum(axis=1)

    keep = (n_a >= min_samples) & (n_b >= min_samples)
    if not keep.any():
        return pd.DataFrame()
    k, n, valid = k[keep], n[keep], valid[keep]

    precision = estimate_precision(k, n, valid, groups=[is_a, is_b])
    mu_a, mu_b, statistic, pvalue = lrt(k, n, valid, is_a, is_b, precision)

    frame = index[keep].to_frame(index=False)[key]
    frame[f"mean_{a_name}"] = mu_a
    frame[f"mean_{b_name}"] = mu_b
    frame["delta_psi"] = mu_b - mu_a
    frame["pvalue"] = pvalue
    frame["lrt_statistic"] = statistic
    frame["precision"] = precision
    frame["n_a"] = n_a[keep]
    frame["n_b"] = n_b[keep]

    if extra:
        meta = df.groupby(key, observed=True)[extra].first().reset_index()
        frame = frame.merge(meta, on=key, how="left")
    return frame


def _ranksum_test(
    df: pd.DataFrame,
    key: list[str],
    extra: list[str],
    value: str,
    a_name: str,
    b_name: str,
    min_samples: int,
) -> pd.DataFrame:
    records = []
    for junc, sub in df.groupby(key, observed=True):
        junc = junc if isinstance(junc, tuple) else (junc,)
        a = sub.loc[sub["condition"] == a_name, value].dropna().to_numpy()
        b = sub.loc[sub["condition"] == b_name, value].dropna().to_numpy()
        if len(a) < min_samples or len(b) < min_samples:
            continue
        try:
            _, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        except ValueError:
            p = np.nan
        rec = dict(zip(key, junc, strict=False))
        first = sub.iloc[0]
        for c in extra:
            rec[c] = first[c]
        rec.update(
            {
                f"mean_{a_name}": float(np.mean(a)),
                f"mean_{b_name}": float(np.mean(b)),
                "delta_psi": float(np.mean(b) - np.mean(a)),
                "pvalue": float(p) if p == p else np.nan,
                "n_a": len(a),
                "n_b": len(b),
            }
        )
        records.append(rec)
    return pd.DataFrame.from_records(records)


def differential_splicing(
    psi_df: pd.DataFrame,
    groups: dict[str, str],
    value: str = "psi_donor",
    min_samples: int = 2,
    key: list[str] | None = None,
    test: str = "auto",
    inc_col: str | None = None,
    total_col: str | None = None,
) -> pd.DataFrame:
    """Test each junction (or event) for differential Ψ between two conditions.

    Parameters
    ----------
    psi_df : a long Ψ table (e.g. from :func:`splicescope.quantify.compute_psi`
        or :func:`splicescope.events.event_psi`).
    groups : mapping ``sample -> condition``; exactly two conditions expected.
    value : which Ψ column to test.
    min_samples : minimum informative (non-NaN) samples required *per group*.
    key : grouping columns identifying a testable unit; defaults to the junction
        coordinates ``[chrom, start, end, strand]``. Pass ``["event_id"]`` for
        event-level tests.
    test : ``"auto"`` (beta-binomial when counts are present, else rank sum),
        ``"betabinom"`` or ``"ranksum"``.
    inc_col, total_col : count columns to model; inferred from ``value`` when
        omitted.

    Returns one row per unit, sorted by q-value.
    """
    conditions = sorted(set(groups.values()))
    if len(conditions) != 2:
        raise ValueError(f"expected exactly 2 conditions, got {conditions}")
    a_name, b_name = conditions
    if test not in ("auto", "betabinom", "ranksum"):
        raise ValueError(f"unknown test: {test!r}")

    df = psi_df.copy()
    df["condition"] = df["sample"].map(groups)
    key = key or ["chrom", "start", "end", "strand"]
    extra = [
        c for c in ("gene_id", "sclass", "event_type") if c in df.columns and c not in key
    ]

    counts = _resolve_count_columns(df, value, inc_col, total_col)
    if test == "betabinom" and counts is None:
        raise ValueError(
            "the beta-binomial test needs inclusion and total count columns; "
            "pass inc_col/total_col or use test='ranksum'"
        )
    use_counts = counts is not None and test in ("auto", "betabinom")

    if use_counts:
        res = _betabinom_test(
            df, key, extra, value, counts, a_name, b_name, groups, min_samples
        )
    else:
        res = _ranksum_test(df, key, extra, value, a_name, b_name, min_samples)
    return _finalize(res)


def significant(diff: pd.DataFrame, q: float = 0.05, min_delta: float = 0.1) -> pd.DataFrame:
    """Filter a differential table to confident hits (``qvalue`` and |ΔΨ|)."""
    if diff.empty:
        return diff
    return diff[(diff["qvalue"] <= q) & (diff["abs_delta_psi"] >= min_delta)].reset_index(drop=True)
