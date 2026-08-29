"""Differential splicing between two groups of samples.

For every junction we compare per-sample Ψ (donor usage) between two conditions:

* ``delta_psi`` = mean(Ψ | group B) − mean(Ψ | group A)
* ``pvalue``    = two-sided Mann–Whitney U test on the per-sample Ψ values
* ``qvalue``    = Benjamini–Hochberg FDR across all tested junctions

Junctions without enough informative samples in both groups are skipped. The
Mann–Whitney test is used because per-sample Ψ is a bounded ratio and we do not
want to assume normality on a handful of replicates.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


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


def differential_splicing(
    psi_df: pd.DataFrame,
    groups: dict[str, str],
    value: str = "psi_donor",
    min_samples: int = 2,
) -> pd.DataFrame:
    """Test each junction for differential donor usage between two conditions.

    Parameters
    ----------
    psi_df : output of :func:`splicescope.quantify.compute_psi`.
    groups : mapping ``sample -> condition``; exactly two conditions expected.
    value : which Ψ column to test.
    min_samples : minimum informative (non-NaN) samples required *per group*.

    Returns one row per junction, sorted by q-value.
    """
    conditions = sorted(set(groups.values()))
    if len(conditions) != 2:
        raise ValueError(f"expected exactly 2 conditions, got {conditions}")
    a_name, b_name = conditions

    df = psi_df.copy()
    df["condition"] = df["sample"].map(groups)
    key = ["chrom", "start", "end", "strand"]
    extra = [c for c in ("gene_id", "sclass") if c in df.columns]

    records = []
    for junc, sub in df.groupby(key, observed=True):
        a = sub.loc[sub["condition"] == a_name, value].dropna().to_numpy()
        b = sub.loc[sub["condition"] == b_name, value].dropna().to_numpy()
        if len(a) < min_samples or len(b) < min_samples:
            continue
        delta = float(np.mean(b) - np.mean(a))
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
                "delta_psi": delta,
                "pvalue": float(p) if p == p else np.nan,  # keep NaN as NaN
                "n_a": len(a),
                "n_b": len(b),
            }
        )
        records.append(rec)

    res = pd.DataFrame.from_records(records)
    if res.empty:
        return res
    res["qvalue"] = benjamini_hochberg(res["pvalue"].to_numpy())
    res["abs_delta_psi"] = res["delta_psi"].abs()
    return res.sort_values(["qvalue", "abs_delta_psi"], ascending=[True, False]).reset_index(
        drop=True
    )


def significant(diff: pd.DataFrame, q: float = 0.05, min_delta: float = 0.1) -> pd.DataFrame:
    """Filter a differential table to confident hits (``qvalue`` and |ΔΨ|)."""
    if diff.empty:
        return diff
    return diff[(diff["qvalue"] <= q) & (diff["abs_delta_psi"] >= min_delta)].reset_index(drop=True)
