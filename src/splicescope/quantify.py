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


def _usage(df: pd.DataFrame, site_cols: list[str], min_reads: int) -> tuple[pd.Series, pd.Series]:
    """Return ``(usage, totals)``; the totals are the denominator Ψ was formed from."""
    totals = df.groupby(site_cols + ["sample"])["count"].transform("sum")
    usage = df["count"] / totals
    usage[totals < min_reads] = np.nan
    return usage, totals


JUNCTION_KEY = ["chrom", "start", "end", "strand"]


def add_unobserved_zeros(df: pd.DataFrame) -> pd.DataFrame:
    """Make the zeros an aligner reports by omission explicit.

    STAR and its kin list only junctions carrying at least one read, so a junction
    missing from a sample arrives looking like missing data when it is in fact a
    measurement: zero reads out of a denominator that is known. For cryptic splicing
    that is the strongest evidence there is — off in every control, on in every
    knockdown — and treating it as missing drops exactly the events the tool exists
    to find.

    A zero row is added only where the junction's donor or acceptor site has reads in
    that sample, so Ψ has a denominator to be a fraction of; a site nobody sequenced
    stays genuinely unobserved. Requires ``donor`` and ``acceptor`` columns.
    """
    if df.empty or "sample" not in df.columns:
        return df
    junctions = df.drop_duplicates(JUNCTION_KEY).drop(columns=["sample", "count"])
    by_site = [
        junctions.merge(
            df[["chrom", site, "strand", "sample"]].drop_duplicates(),
            on=["chrom", site, "strand"],
        )
        for site in ("donor", "acceptor")
    ]
    grid = pd.concat(by_site, ignore_index=True).drop_duplicates(JUNCTION_KEY + ["sample"])

    seen = df[JUNCTION_KEY + ["sample"]].drop_duplicates()
    missing = grid.merge(seen, on=JUNCTION_KEY + ["sample"], how="left", indicator=True)
    missing = missing[missing["_merge"] == "left_only"].drop(columns="_merge")
    if missing.empty:
        return df
    missing["count"] = 0
    return pd.concat([df, missing[df.columns]], ignore_index=True)


def compute_psi(
    annotated: pd.DataFrame, min_reads: int = 10, fill_unobserved: bool = True
) -> pd.DataFrame:
    """Add ``psi_donor`` and ``psi_acceptor`` columns (per sample) to junctions.

    ``annotated`` must have ``[chrom, start, end, strand, sample, count]``. The
    matching ``donor_total`` and ``acceptor_total`` denominators are kept so a
    count-based test can use the reads Ψ was computed from, not just the ratio.

    ``fill_unobserved`` adds the zero-count rows the aligner left out — see
    :func:`add_unobserved_zeros`. Turn it off only to reproduce the pre-0.9.0
    behaviour, which silently discarded junctions absent from a whole group.
    """
    df = annotated.copy()
    da = [
        donor_acceptor(s, e, st)
        for s, e, st in zip(df["start"], df["end"], df["strand"], strict=False)
    ]
    df["donor"] = [d for d, _ in da]
    df["acceptor"] = [a for _, a in da]
    if fill_unobserved:
        df = add_unobserved_zeros(df)
    df["psi_donor"], df["donor_total"] = _usage(df, ["chrom", "donor", "strand"], min_reads)
    df["psi_acceptor"], df["acceptor_total"] = _usage(
        df, ["chrom", "acceptor", "strand"], min_reads
    )
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
