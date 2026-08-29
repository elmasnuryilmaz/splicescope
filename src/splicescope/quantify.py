"""Quantify splice-site usage (Ψ) from junction read counts.

We use a well-defined, model-free metric: **donor usage**. For a junction *j*
with 5' splice site (donor) *D*, its usage in a sample is

    Ψ(j) = count(j) / Σ count(all junctions sharing donor D)

i.e. the fraction of that donor's spliced reads that flow through *j*. This is
the quantity differential-splicing methods build on and needs no event model.
Acceptor usage is defined symmetrically and returned alongside.

A minimum coverage threshold avoids dividing tiny counts: sites with fewer than
``min_reads`` total reads in a sample yield ``NaN`` (uninformative) rather than a
noisy ratio.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .io import donor_acceptor


def _usage(df: pd.DataFrame, site_cols: list[str], min_reads: int) -> pd.Series:
    totals = df.groupby(site_cols + ["sample"])["count"].transform("sum")
    usage = df["count"] / totals
    usage[totals < min_reads] = np.nan
    return usage


def compute_psi(annotated: pd.DataFrame, min_reads: int = 10) -> pd.DataFrame:
    """Add ``psi_donor`` and ``psi_acceptor`` columns (per sample) to junctions.

    ``annotated`` must have ``[chrom, start, end, strand, sample, count]``.
    """
    df = annotated.copy()
    da = [
        donor_acceptor(s, e, st)
        for s, e, st in zip(df["start"], df["end"], df["strand"], strict=False)
    ]
    df["donor"] = [d for d, _ in da]
    df["acceptor"] = [a for _, a in da]
    df["psi_donor"] = _usage(df, ["chrom", "donor", "strand"], min_reads)
    df["psi_acceptor"] = _usage(df, ["chrom", "acceptor", "strand"], min_reads)
    return df


def psi_matrix(psi_df: pd.DataFrame, value: str = "psi_donor") -> pd.DataFrame:
    """Pivot to a junction × sample matrix of Ψ values (junctions as the index)."""
    keyed = psi_df.assign(
        junction=psi_df["chrom"].astype(str)
        + ":"
        + psi_df["start"].astype(str)
        + "-"
        + psi_df["end"].astype(str)
        + ":"
        + psi_df["strand"].astype(str)
    )
    return keyed.pivot_table(index="junction", columns="sample", values=value, aggfunc="mean")
