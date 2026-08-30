"""Gene-set over-representation analysis (ORA) for differentially spliced genes.

After differential splicing, a natural question is *which pathways are affected?*
This module runs the standard hypergeometric (Fisher's exact, one-sided) test used
by tools like clusterProfiler: given a set of "hit" genes, a background, and a
collection of gene sets (GO terms, KEGG pathways, MSigDB, …), it asks whether each
set is over-represented among the hits, and corrects across sets with BH-FDR.

The method is generic — bring any ``{term: [genes]}`` mapping (e.g. via
:func:`splicescope.io.read_gmt`). No gene sets are bundled, because meaningful
enrichment needs real annotations rather than synthetic ones.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence

import pandas as pd
from scipy import stats

from .diff import benjamini_hochberg


def over_representation(
    hits: Iterable[str],
    background: Iterable[str],
    gene_sets: Mapping[str, Sequence[str]],
    min_size: int = 2,
    max_size: int | None = None,
) -> pd.DataFrame:
    """Hypergeometric over-representation of ``gene_sets`` among ``hits``.

    Parameters
    ----------
    hits : the genes of interest (e.g. differentially spliced).
    background : all genes that could have been a hit (the universe).
    gene_sets : mapping of term -> member genes.
    min_size, max_size : restrict to sets of this size *within the background*.

    Returns one row per tested set with the 2×2 counts, fold enrichment, p-value
    and BH q-value, sorted by q-value. Uses the survival function
    ``P(X ≥ k) = hypergeom.sf(k-1, M, n, N)`` with ``M`` background size, ``n`` set
    size in background, ``N`` number of hits in background, ``k`` the overlap.
    """
    bg = set(background)
    hit_set = set(hits) & bg
    M, N = len(bg), len(hit_set)
    if M == 0 or N == 0:
        return pd.DataFrame(
            columns=[
                "term", "set_size", "overlap", "n_hits", "n_background",
                "fold_enrichment", "pvalue", "qvalue", "genes",
            ]
        )

    records = []
    for term, genes in gene_sets.items():
        in_bg = set(genes) & bg
        n = len(in_bg)
        if n < min_size or (max_size is not None and n > max_size):
            continue
        overlap = hit_set & in_bg
        k = len(overlap)
        if k == 0:
            continue
        p = float(stats.hypergeom.sf(k - 1, M, n, N))
        fold = (k / N) / (n / M)
        records.append(
            {
                "term": term,
                "set_size": n,
                "overlap": k,
                "n_hits": N,
                "n_background": M,
                "fold_enrichment": fold,
                "pvalue": p,
                "genes": ",".join(sorted(overlap)),
            }
        )

    res = pd.DataFrame.from_records(records)
    if res.empty:
        return res
    res["qvalue"] = benjamini_hochberg(res["pvalue"].to_numpy())
    return res.sort_values(["qvalue", "pvalue"]).reset_index(drop=True)


def enrich_differential(
    diff_table: pd.DataFrame,
    gene_sets: Mapping[str, Sequence[str]],
    q: float = 0.05,
    min_delta: float = 0.1,
    gene_col: str = "gene_id",
    **kwargs,
) -> pd.DataFrame:
    """Convenience: ORA of significant genes from a differential table.

    Hits are the genes of significant units (``q ≤`` threshold and ``|ΔΨ|`` ≥
    ``min_delta``); the background is every gene that was tested.
    """
    from .diff import significant

    if diff_table.empty or gene_col not in diff_table.columns:
        return over_representation([], [], gene_sets, **kwargs)
    background = diff_table[gene_col].dropna().unique().tolist()
    hits = significant(diff_table, q=q, min_delta=min_delta)[gene_col].dropna().unique().tolist()
    return over_representation(hits, background, gene_sets, **kwargs)
