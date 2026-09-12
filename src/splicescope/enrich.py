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

import re
from collections.abc import Iterable, Mapping, Sequence

import pandas as pd
from scipy import stats

from .diff import benjamini_hochberg

_VERSION_SUFFIX = re.compile(r"\.\d+$")
#: What a de-versioned stem must look like before the suffix is treated as a version:
#: Ensembl (``ENSG``, ``ENSMUST``, …) or RefSeq-style (``NM_``, ``XP_``, …) accessions.
_ACCESSION = re.compile(r"^(?:ENS[A-Z]*[EGTP]\d+|[A-Z]{2}_\d+)$", re.IGNORECASE)
#: Pseudo-counts added to each bin's hit rate (a Jeffreys prior). Without them a bin whose
#: members all are, or all are not, hits gives a rate of exactly 1 or 0, and the odds
#: ``p/(1-p)`` that :func:`_bias_odds` averages become infinite in one direction. Clipping
#: to a small epsilon instead is worse than useless: it replaces infinity with an arbitrary
#: 1e6 that still swamps every other gene in the average.
_PSEUDO_HITS = 0.5
_PSEUDO_TOTAL = 1.0


def normalize_gene_id(gene: str) -> str:
    """Put an identifier in the form both sides of an intersection can agree on.

    A GTF gives versioned accessions (``ENSG00000141510.16``, ``NM_000546.6``); gene-set
    files give either symbols or *unversioned* accessions, so a verbatim comparison
    matches nothing at all. The version suffix is dropped and the rest upper-cased.

    Only from accessions, though. Stripping ``.<digits>`` from everything looks safe for
    human symbols and is not safe in general: *C. elegans* sequence names — the standard
    identifier for the majority of worm genes, which have no CGC name — are exactly that
    shape, so ``C42D8.1``, ``C42D8.2`` and ``C42D8.3`` would collapse to one gene and take
    the background size, the set size, the overlap and the p-value down with them.
    """
    text = str(gene).strip()
    stem = _VERSION_SUFFIX.sub("", text)
    if stem != text and _ACCESSION.match(stem):
        text = stem
    return text.upper()


def _normalized_index(values: Iterable[str]) -> dict[str, str]:
    """Map each normalised identifier to the first original spelling seen."""
    index: dict[str, str] = {}
    for value in values:
        index.setdefault(normalize_gene_id(value), str(value))
    return index


def selection_propensity(
    background: Iterable[str],
    hits: Iterable[str],
    weights: Mapping[str, float],
    bins: int = 10,
) -> dict[str, float]:
    """Estimated P(hit) for each background gene, given how much opportunity it had.

    In splicing data "opportunity" is the number of units tested in that gene: a long,
    many-exon gene offers far more chances to contain a significant junction than a
    two-exon one, for reasons that have nothing to do with the biology. Genes are binned
    by that measure and each bin's observed hit rate becomes the propensity of its
    members — the same device ``goseq`` uses for gene length in RNA-seq.
    """
    bg = sorted({normalize_gene_id(g) for g in background})
    hit_set = {normalize_gene_id(h) for h in hits}
    if not bg:
        return {}
    baseline = len(hit_set & set(bg)) / len(bg)

    # Bin by weight, never splitting a tied group. Most genes contribute exactly one
    # tested unit, so a plain equal-count split would cut that block at an arbitrary
    # point and two genes with identical opportunity would get different propensities
    # according to their names.
    by_weight: dict[float, list[str]] = {}
    for gene in bg:
        by_weight.setdefault(float(weights.get(gene, 0.0)), []).append(gene)
    target = max(1, len(bg) // max(1, bins))
    chunks, current = [], []
    for weight in sorted(by_weight):
        current.extend(by_weight[weight])
        if len(current) >= target:
            chunks.append(current)
            current = []
    if current:
        # Its own bin, never folded into the previous one. A single large tied group can
        # fill a bin on its own — most genes contribute exactly one unit — and folding the
        # remainder into it would put the highest-opportunity genes in the same bin as the
        # lowest, which is the distinction this whole estimate exists to make.
        chunks.append(current)

    propensity: dict[str, float] = {}
    for members in chunks:
        if not members:
            continue
        observed = sum(1 for g in members if g in hit_set)
        # Shrink towards 1/2 by the bin's own size, so a small bin cannot assert certainty.
        # A 7-gene bin whose members are all hits gives 0.94 (odds 15), not 1 (odds 1e6);
        # a 200-gene bin with none gives 0.0025. The smaller the bin, the harder it shrinks.
        rate = (observed + _PSEUDO_HITS) / (len(members) + _PSEUDO_TOTAL)
        for gene in members:
            propensity[gene] = rate
    if not propensity:
        smoothed = (baseline * len(bg) + _PSEUDO_HITS) / (len(bg) + _PSEUDO_TOTAL)
        return dict.fromkeys(bg, smoothed)
    return propensity


def _bias_odds(in_set: set, bg: set, propensity: Mapping[str, float]) -> float | None:
    """Wallenius odds for a gene set: how much likelier its genes were to be drawn.

    Wallenius' ``ω`` is a ratio of sampling *weights*, so the quantity to average is the
    odds ``p / (1 - p)``, not the probability. ``goseq`` averages the probabilities, which
    is the same thing only while ``p`` is small; here propensities reach 0.5. Measured on
    biologically null gene sets biased toward large genes, at identical power (80 % at a
    0.10 enrichment, 100 % above it): averaging probabilities leaves **38–40 %** of them
    called at p ≤ 0.05, averaging odds leaves **0 %**.

    The odds diverge as ``p`` approaches 1, which is why :func:`selection_propensity`
    smooths its rates rather than clipping them. A bin whose members all happened to be
    hits would otherwise contribute an arbitrary 1e6 and decide the mean by itself — and
    since the same genes sit in the denominator for every other set, one of them rewrote
    the whole table.
    """
    outside = bg - in_set
    if not in_set or not outside:
        return None
    inside = sum(propensity[g] / (1.0 - propensity[g]) for g in in_set) / len(in_set)
    beyond = sum(propensity[g] / (1.0 - propensity[g]) for g in outside) / len(outside)
    if beyond <= 0 or inside <= 0:
        return None
    return inside / beyond


def over_representation(
    hits: Iterable[str],
    background: Iterable[str],
    gene_sets: Mapping[str, Sequence[str]],
    min_size: int = 2,
    max_size: int | None = None,
    weights: Mapping[str, float] | None = None,
) -> pd.DataFrame:
    """Hypergeometric over-representation of ``gene_sets`` among ``hits``.

    Parameters
    ----------
    hits : the genes of interest (e.g. differentially spliced).
    background : all genes that could have been a hit (the universe).
    gene_sets : mapping of term -> member genes.
    min_size, max_size : restrict to sets of this size *within the background*.

    weights : per-gene opportunity (for splicing, the number of units tested in that
        gene). Supplied, the null stops assuming every gene was equally likely to be a
        hit: genes are binned by opportunity, each bin's observed hit rate becomes its
        members' propensity, and the set's p-value comes from Wallenius' non-central
        hypergeometric with odds = (mean of ``p/(1-p)`` inside) / (mean outside) — see
        :func:`_bias_odds`. This is ``goseq``'s device for gene length in RNA-seq, and it
        matters here because a long, many-exon gene has many more chances to contain a
        significant junction. Omitted, the plain hypergeometric is used, which is what
        Wallenius becomes at odds 1.

    Returns one row per tested set with the 2×2 counts, fold enrichment, p-value
    and BH q-value, sorted by q-value. Uses the survival function
    ``P(X ≥ k) = hypergeom.sf(k-1, M, n, N)`` with ``M`` background size, ``n`` set
    size in background, ``N`` number of hits in background, ``k`` the overlap, plus
    ``bias_odds`` when weights were given.
    """
    bg_index = _normalized_index(background)
    bg = set(bg_index)
    hit_set = {normalize_gene_id(h) for h in hits} & bg
    M, N = len(bg), len(hit_set)
    propensity = (
        selection_propensity(bg, hit_set, {normalize_gene_id(g): w for g, w in weights.items()})
        if weights
        else None
    )
    if M == 0 or N == 0:
        return pd.DataFrame(
            columns=[
                "term", "set_size", "overlap", "n_hits", "n_background",
                "fold_enrichment", "bias_odds", "pvalue", "qvalue", "genes",
            ]
        )

    records = []
    for term, genes in gene_sets.items():
        in_bg = {normalize_gene_id(g) for g in genes} & bg
        n = len(in_bg)
        if n < min_size or (max_size is not None and n > max_size):
            continue
        overlap = hit_set & in_bg
        k = len(overlap)
        if k == 0:
            continue
        odds = _bias_odds(in_bg, bg, propensity) if propensity else None
        if odds is None:
            p = float(stats.hypergeom.sf(k - 1, M, n, N))
        else:
            p = float(stats.nchypergeom_wallenius.sf(k - 1, M, n, N, odds))
        fold = (k / N) / (n / M)
        records.append(
            {
                "term": term,
                "set_size": n,
                "overlap": k,
                "n_hits": N,
                "n_background": M,
                "fold_enrichment": fold,
                "bias_odds": odds if odds is not None else 1.0,
                "pvalue": p,
                "genes": ",".join(sorted(bg_index[g] for g in overlap)),
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
    name_col: str = "gene_name",
    weight_by_units: bool = True,
    **kwargs,
) -> pd.DataFrame:
    """Convenience: ORA of significant genes from a differential table.

    Hits are the genes of significant units (``q ≤`` threshold and ``|ΔΨ|`` ≥
    ``min_delta``); the background is every gene that was tested.

    Gene sets are keyed by symbols (MSigDB ``*.symbols.gmt``, GO, KEGG) about as often
    as by accessions, and a GTF supplies both, so whichever of ``gene_col`` and
    ``name_col`` overlaps the sets more is the one used. Identifiers are matched through
    :func:`normalize_gene_id` rather than verbatim.

    ``weight_by_units`` corrects the bias that makes gene-level ORA on splicing data
    misleading: a gene contributes as many chances of being a hit as it has tested
    junctions, so large genes are over-represented among the hits for reasons that are
    not biological. The number of rows per gene becomes its opportunity measure — see
    ``weights`` in :func:`over_representation`. Set it to ``False`` for the plain
    hypergeometric, which is what every ORA tool does and what this did before.
    """
    from .diff import significant

    columns = [c for c in (gene_col, name_col) if c and c in diff_table.columns]
    if diff_table.empty or not columns:
        return over_representation([], [], gene_sets, **kwargs)

    members = {normalize_gene_id(g) for genes in gene_sets.values() for g in genes}
    chosen = max(
        columns,
        key=lambda c: len(
            {normalize_gene_id(g) for g in diff_table[c].dropna().unique()} & members
        ),
    )
    background = diff_table[chosen].dropna().unique().tolist()
    hits = significant(diff_table, q=q, min_delta=min_delta)[chosen].dropna().unique().tolist()
    if weight_by_units and "weights" not in kwargs:
        kwargs["weights"] = diff_table[chosen].dropna().value_counts().to_dict()
    return over_representation(hits, background, gene_sets, **kwargs)
