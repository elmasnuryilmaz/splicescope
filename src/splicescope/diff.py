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

import warnings

import numpy as np
import pandas as pd
from scipy import stats

from .betabinom import (
    dispersion_is_estimable,
    estimate_precision,
    estimate_precision_per_unit,
    lrt,
)

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


def _wide(
    df: pd.DataFrame, key: list[str], column: str, samples: list[str], aggfunc: str = "sum"
) -> pd.DataFrame:
    """Pivot one column to a unit × sample matrix.

    ``aggfunc`` matters more than it looks: pandas sums an all-NaN group to ``0.0``,
    so summing Ψ turns "uninformative" into "zero usage" and erases the ``min_reads``
    coverage filter entirely. Counts are summed; Ψ is averaged, which leaves NaN as NaN.
    """
    table = df.pivot_table(index=key, columns="sample", values=column, aggfunc=aggfunc)
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
    dispersion: str = "shared",
) -> pd.DataFrame:
    inc_col, total_col = counts
    samples = sorted(df["sample"].unique())
    # The summed count pivot keeps every unit, so use its index as the canonical one
    # and align the others to it; the averaged Ψ pivot can drop an all-NaN unit.
    n_wide = _wide(df, key, total_col, samples)
    index = n_wide.index
    psi_wide = _wide(df, key, value, samples, aggfunc="mean").reindex(index)

    valid = np.array(psi_wide.notna().to_numpy(), dtype=bool)  # writable copy
    k = np.nan_to_num(_wide(df, key, inc_col, samples).reindex(index).to_numpy(), nan=0.0)
    n = np.nan_to_num(n_wide.to_numpy(), nan=0.0)
    valid &= n > 0

    is_a = np.array([groups.get(s) == a_name for s in samples])
    is_b = np.array([groups.get(s) == b_name for s in samples])
    n_a = (valid & is_a[None, :]).sum(axis=1)
    n_b = (valid & is_b[None, :]).sum(axis=1)

    keep = (n_a >= min_samples) & (n_b >= min_samples)
    if not keep.any():
        return pd.DataFrame()
    k, n, valid = k[keep], n[keep], valid[keep]

    if not dispersion_is_estimable(n, valid, groups=[is_a, is_b]):
        warnings.warn(
            "dispersion could not be estimated: no unit has more informative replicates "
            "than fitted group means, so replicate-to-replicate variability is unknown. "
            "The beta-binomial test falls back to assuming none, which makes it as narrow "
            "as a binomial test and its p-values far too small — a 1-vs-1 comparison of "
            "Psi 0.300 against 0.360 at 1000 reads returns q=4.5e-03 on no replication at "
            "all. Add replicates, raise min_samples, or do not treat these p-values as "
            "evidence.",
            stacklevel=3,
        )
    precision = estimate_precision(k, n, valid, groups=[is_a, is_b])
    if dispersion == "per_unit_floor":
        per_unit = estimate_precision_per_unit(k, n, valid, groups=[is_a, is_b])
        precision = np.minimum(precision, per_unit)
        mu_a, mu_b, statistic, pvalue = lrt(k, n, valid, is_a, is_b, precision[:, None])
    else:
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


def _check_inputs(
    psi_df: pd.DataFrame, groups: dict[str, str], value: str, key: list[str]
) -> None:
    """Refuse the inputs that would otherwise produce a clean run and no result.

    The dangerous one is the last: if no sample in the table is named in ``groups``,
    every row is unassigned, every unit fails ``min_samples``, and the function returns
    an empty table. Nothing raises, the caller reports "0 significant junctions" and
    believes it. For a tool whose job is to find differential splicing that is the worst
    outcome available, so it is an error here rather than a silence. The CLI catches the
    same mistake from the filenames; this catches it for everyone using the library.

    The others only replace a ``KeyError`` raised somewhere inside pandas with a
    sentence naming the column and what was available instead.
    """
    if "sample" not in psi_df.columns:
        raise ValueError(
            "the Psi table has no 'sample' column, so no row can be assigned to a "
            "condition. Pass the long table from compute_psi or event_psi, not a matrix."
        )
    if value not in psi_df.columns:
        candidates = [c for c in psi_df.columns if str(c).startswith("psi")]
        hint = f"; did you mean one of {candidates}?" if candidates else ""
        raise ValueError(f"no column {value!r} in the Psi table{hint}")
    absent = [c for c in key if c not in psi_df.columns]
    if absent:
        raise ValueError(
            f"the grouping key names {absent}, which the Psi table does not have. "
            f"Its columns are {list(psi_df.columns)}."
        )

    named = set(groups)
    present = set(psi_df["sample"].dropna().unique())
    if not named & present:
        raise ValueError(
            "no sample in the Psi table is named in `groups`, so every row would be "
            "unassigned and this would return an empty table instead of failing.\n"
            f"  in the table: {sorted(map(str, present))[:8]}\n"
            f"  in `groups`:  {sorted(map(str, named))[:8]}"
        )
    for missing, where in (
        (sorted(map(str, present - named)), "in the Psi table but not in `groups`"),
        (sorted(map(str, named - present)), "in `groups` but not in the Psi table"),
    ):
        if missing:
            shown = ", ".join(missing[:6]) + (" ..." if len(missing) > 6 else "")
            warnings.warn(
                f"{len(missing)} sample(s) {where} and will be ignored: {shown}",
                stacklevel=3,
            )


def differential_splicing(
    psi_df: pd.DataFrame,
    groups: dict[str, str],
    value: str = "psi_donor",
    min_samples: int = 2,
    key: list[str] | None = None,
    test: str = "auto",
    inc_col: str | None = None,
    total_col: str | None = None,
    dispersion: str = "shared",
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
    dispersion : ``"shared"`` (default) uses one precision for every unit;
        ``"per_unit_floor"`` additionally estimates each unit's own precision and
        takes the smaller of the two, which can only widen the null. Use it when
        dispersion is likely heterogeneous — the shared value is anti-conservative
        for anything noisier than itself. It costs power: on simulated homogeneous
        data (Ψ 0.30 vs 0.45, 3 vs 3, 60 reads) power falls from 0.52 to 0.44, while
        on heterogeneous null data the false-positive rate of the loosely dispersed
        units falls from 0.173 to 0.090 and the overall rate from 0.087 to 0.045.

    Returns one row per unit, sorted by q-value.
    """
    conditions = sorted(set(groups.values()))
    if len(conditions) != 2:
        raise ValueError(f"expected exactly 2 conditions, got {conditions}")
    a_name, b_name = conditions
    if test not in ("auto", "betabinom", "ranksum"):
        raise ValueError(f"unknown test: {test!r}")
    if dispersion not in ("shared", "per_unit_floor"):
        raise ValueError(f"unknown dispersion: {dispersion!r}")
    key = key or ["chrom", "start", "end", "strand"]
    _check_inputs(psi_df, groups, value, key)

    df = psi_df.copy()
    df["condition"] = df["sample"].map(groups)
    extra = [
        c
        for c in ("gene_id", "gene_name", "sclass", "event_type")
        if c in df.columns and c not in key
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
            df, key, extra, value, counts, a_name, b_name, groups, min_samples, dispersion
        )
    else:
        res = _ranksum_test(df, key, extra, value, a_name, b_name, min_samples)
    return _finalize(res)


def significant(diff: pd.DataFrame, q: float = 0.05, min_delta: float = 0.1) -> pd.DataFrame:
    """Filter a differential table to confident hits (``qvalue`` and |ΔΨ|)."""
    if diff.empty:
        return diff
    return diff[(diff["qvalue"] <= q) & (diff["abs_delta_psi"] >= min_delta)].reset_index(drop=True)
