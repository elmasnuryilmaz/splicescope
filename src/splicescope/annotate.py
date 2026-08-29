"""Classify observed junctions against a reference annotation.

Each observed junction is labelled by how it relates to the known splice sites:

======================  ==================================================
class                   meaning
======================  ==================================================
``annotated``           donor+acceptor pair is a known intron
``novel_combination``   both sites known, but not as a pair (e.g. exon skip)
``novel_donor``         acceptor known, 5' site novel
``novel_acceptor``      donor known, 3' site novel
``cryptic``             both sites novel — candidate cryptic splicing
======================  ==================================================
"""

from __future__ import annotations

import pandas as pd

from .io import donor_acceptor

CLASSES = ["annotated", "novel_combination", "novel_donor", "novel_acceptor", "cryptic"]


def _site_sets(known: pd.DataFrame):
    junctions = set(
        zip(known["chrom"], known["start"], known["end"], known["strand"], strict=False)
    )
    donors, acceptors, gene_of_site = set(), set(), {}
    gene_of_junction = {}
    for chrom, start, end, strand, gene in known[
        ["chrom", "start", "end", "strand", "gene_id"]
    ].itertuples(index=False):
        d, a = donor_acceptor(start, end, strand)
        donors.add((chrom, d, strand))
        acceptors.add((chrom, a, strand))
        gene_of_site.setdefault((chrom, d, strand), gene)
        gene_of_site.setdefault((chrom, a, strand), gene)
        gene_of_junction[(chrom, start, end, strand)] = gene
    return junctions, donors, acceptors, gene_of_site, gene_of_junction


def classify_one(chrom, start, end, strand, junctions, donors, acceptors) -> str:
    """Classify a single junction. Pure function for easy testing."""
    if (chrom, start, end, strand) in junctions:
        return "annotated"
    d, a = donor_acceptor(start, end, strand)
    d_known = (chrom, d, strand) in donors
    a_known = (chrom, a, strand) in acceptors
    if d_known and a_known:
        return "novel_combination"
    if d_known and not a_known:
        return "novel_acceptor"
    if a_known and not d_known:
        return "novel_donor"
    return "cryptic"


def annotate_junctions(observed: pd.DataFrame, known: pd.DataFrame) -> pd.DataFrame:
    """Annotate a table of observed junctions.

    ``observed`` needs columns ``[chrom, start, end, strand]`` (extra columns are
    preserved). Returns a copy with two new columns: ``sclass`` (one of
    :data:`CLASSES`) and ``gene_id`` (best-effort assignment by shared splice site).
    """
    junctions, donors, acceptors, gene_of_site, gene_of_junction = _site_sets(known)

    sclass, genes = [], []
    for chrom, start, end, strand in observed[["chrom", "start", "end", "strand"]].itertuples(
        index=False
    ):
        cls = classify_one(chrom, start, end, strand, junctions, donors, acceptors)
        sclass.append(cls)
        d, a = donor_acceptor(start, end, strand)
        gene = (
            gene_of_junction.get((chrom, start, end, strand))
            or gene_of_site.get((chrom, d, strand))
            or gene_of_site.get((chrom, a, strand))
        )
        genes.append(gene)

    out = observed.copy()
    out["sclass"] = pd.Categorical(sclass, categories=CLASSES)
    out["gene_id"] = genes
    out["is_novel"] = out["sclass"] != "annotated"
    return out


def annotation_summary(annotated: pd.DataFrame) -> pd.DataFrame:
    """Counts of unique junctions per class (deduplicated across samples)."""
    uniq = annotated.drop_duplicates(subset=["chrom", "start", "end", "strand"])
    counts = uniq["sclass"].value_counts().reindex(CLASSES, fill_value=0)
    return counts.rename_axis("sclass").reset_index(name="n_junctions")
