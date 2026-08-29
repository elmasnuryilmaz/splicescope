"""Feature engineering for cryptic-event calling.

Novel junctions are a mix of genuine cryptic splicing and technical noise
(mis-alignment, template switching, low-level transcriptional leakage). We turn
each junction into a feature vector that a classifier can use to tell them apart:

======================  ==================================================
feature                 intuition
======================  ==================================================
``intron_length``       cryptic introns have a plausible length distribution
``log_max_count``       real events are reproducibly supported by reads
``n_samples_support``   real events recur across replicates; noise does not
``mean_psi_donor``      real events take a non-trivial share of their donor
``canonical_motif``     genuine splicing tends to use GT/AG (or CT/AC)
``dist_known_donor``    cryptic sites sit near, but not on, annotated sites
``dist_known_acceptor`` same, for the 3' site
``is_novel_both``       both sites unannotated (the strongest cryptic prior)
======================  ==================================================
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .io import donor_acceptor

_CANONICAL = {"GT/AG", "CT/AC"}


def _safe_nanmean(values: np.ndarray) -> float:
    """np.nanmean without the empty-slice warning (all-NaN -> NaN)."""
    values = np.asarray(values, dtype=float)
    if values.size == 0 or np.all(np.isnan(values)):
        return float("nan")
    return float(np.nanmean(values))


def _known_site_index(known: pd.DataFrame) -> dict[tuple[str, str], np.ndarray]:
    """Sorted arrays of known donor+acceptor positions per (chrom, strand)."""
    sites: dict[tuple[str, str], set[int]] = {}
    for chrom, start, end, strand in known[["chrom", "start", "end", "strand"]].itertuples(
        index=False
    ):
        d, a = donor_acceptor(start, end, strand)
        sites.setdefault((chrom, strand), set()).update((d, a))
    return {k: np.array(sorted(v)) for k, v in sites.items()}


def _nearest_distance(pos: int, sorted_sites: np.ndarray) -> float:
    if sorted_sites.size == 0:
        return np.nan
    i = int(np.searchsorted(sorted_sites, pos))
    best = np.inf
    for j in (i - 1, i):
        if 0 <= j < sorted_sites.size:
            best = min(best, abs(int(sorted_sites[j]) - pos))
    return float(best)


def extract_features(
    psi_annotated: pd.DataFrame, known: pd.DataFrame, novel_only: bool = True
) -> pd.DataFrame:
    """Build a per-junction feature table from an annotated+quantified long table.

    ``psi_annotated`` is the output of :func:`quantify.compute_psi` run on an
    annotated table (so it carries ``sclass`` and ``psi_donor``). One row per
    unique junction is returned, ready for :mod:`splicescope.ml`.
    """
    df = psi_annotated
    if novel_only:
        df = df[df["sclass"] != "annotated"]

    site_index = _known_site_index(known)
    key = ["chrom", "start", "end", "strand"]
    rows = []
    for (chrom, start, end, strand), sub in df.groupby(key, observed=True):
        d, a = donor_acceptor(start, end, strand)
        sites = site_index.get((chrom, strand), np.array([], dtype=int))
        motif = sub["motif"].iloc[0] if "motif" in sub else "non-canonical"
        sclass = sub["sclass"].iloc[0]
        counts = sub["count"].to_numpy()
        rows.append(
            {
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "gene_id": sub["gene_id"].iloc[0] if "gene_id" in sub else None,
                "sclass": str(sclass),
                "intron_length": int(end - start + 1),
                "log_max_count": float(np.log1p(counts.max())),
                "n_samples_support": int((counts > 0).sum()),
                "mean_psi_donor": _safe_nanmean(sub["psi_donor"].to_numpy())
                if "psi_donor" in sub
                else np.nan,
                "canonical_motif": int(motif in _CANONICAL),
                "dist_known_donor": _nearest_distance(d, sites),
                "dist_known_acceptor": _nearest_distance(a, sites),
                "is_novel_both": int(str(sclass) == "cryptic"),
            }
        )
    feats = pd.DataFrame(rows)
    if "is_cryptic_truth" in df.columns and not feats.empty:
        truth = (
            df.groupby(key, observed=True)["is_cryptic_truth"].max().reset_index(drop=True)
        )
        feats["is_cryptic_truth"] = truth.astype(int)
    return feats


FEATURE_COLUMNS = [
    "intron_length",
    "log_max_count",
    "n_samples_support",
    "mean_psi_donor",
    "canonical_motif",
    "dist_known_donor",
    "dist_known_acceptor",
    "is_novel_both",
]
