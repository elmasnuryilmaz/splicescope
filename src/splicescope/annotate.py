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
    has_name = "gene_name" in known.columns
    columns = ["chrom", "start", "end", "strand", "gene_id"] + (["gene_name"] if has_name else [])
    # Each column is converted once. pandas 3 backs string columns with Arrow, and
    # reading one a value at a time costs several times what is then done with the
    # values — on a GENCODE-sized annotation it was the largest single cost in
    # annotating a whole experiment.
    chroms, starts, ends, strands, genes, *rest = (known[c].to_numpy() for c in columns)
    # keep the symbol alongside the accession: gene-set files use one or the other
    names = rest[0] if has_name else genes

    junctions = set(zip(chroms, starts, ends, strands, strict=True))
    donors, acceptors, gene_of_site = set(), set(), {}
    gene_of_junction = {}
    for chrom, start, end, strand, gene, name in zip(
        chroms, starts, ends, strands, genes, names, strict=True
    ):
        value = (gene, name)
        d, a = donor_acceptor(start, end, strand)
        donors.add((chrom, d, strand))
        acceptors.add((chrom, a, strand))
        gene_of_site.setdefault((chrom, d, strand), value)
        gene_of_site.setdefault((chrom, a, strand), value)
        gene_of_junction[(chrom, start, end, strand)] = value
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


def resolve_unstranded(observed: pd.DataFrame, known: pd.DataFrame) -> pd.DataFrame:
    """Give strand-undefined junctions the strand their annotation implies.

    STAR writes strand code 0 — which this package reads as ``"."`` — whenever it cannot
    infer the strand from the intron motif, which is routine for non-canonical junctions.
    Such a junction can never equal a stranded annotation, so it is classified ``cryptic``
    however ordinary it is, loses its ``gene_id``, and — being the only junction at its own
    ``(chrom, position, ".")`` site — is handed Ψ ≡ 1.0 in every sample. It then occupies a
    row in the differential table that can never show a difference.

    Where the annotation knows the intron, or either of its splice sites, the strand is
    recoverable. Junctions it cannot place are left as ``"."``: the strand genuinely is
    unknown, and guessing one would be worse than saying so.
    """
    if "strand" not in observed.columns:
        return observed
    unstranded = observed["strand"] == "."
    if not unstranded.any() or known.empty:
        return observed

    by_intron = {
        (c, s, e): st
        for c, s, e, st in known[["chrom", "start", "end", "strand"]].itertuples(index=False)
    }
    by_site: dict[tuple, str] = {}
    for c, s, e, st in known[["chrom", "start", "end", "strand"]].itertuples(index=False):
        by_site.setdefault((c, s), st)
        by_site.setdefault((c, e), st)

    out = observed.copy()
    resolved = []
    for chrom, start, end in out.loc[unstranded, ["chrom", "start", "end"]].itertuples(
        index=False
    ):
        resolved.append(
            by_intron.get((chrom, start, end))
            or by_site.get((chrom, start))
            or by_site.get((chrom, end))
            or "."
        )
    out.loc[unstranded, "strand"] = resolved
    return out


def annotate_junctions(observed: pd.DataFrame, known: pd.DataFrame) -> pd.DataFrame:
    """Annotate a table of observed junctions.

    ``observed`` needs columns ``[chrom, start, end, strand]`` (extra columns are
    preserved). Returns a copy with two new columns: ``sclass`` (one of
    :data:`CLASSES`) and ``gene_id`` (best-effort assignment by shared splice site).

    Strand-undefined junctions are first placed against the annotation where they can
    be — see :func:`resolve_unstranded`.
    """
    observed = resolve_unstranded(observed, known)
    junctions, donors, acceptors, gene_of_site, gene_of_junction = _site_sets(known)

    # A junction's class and gene depend only on its coordinates, so each distinct one
    # is resolved once and the answer reused. A six-sample experiment lists every
    # junction six times, and on a human-sized table four fifths of the rows are repeats.
    # The keys come from to_numpy() rather than itertuples(): pandas 3 backs string
    # columns with Arrow, and pulling 1.3 million values out of one an element at a time
    # cost more than the classification did.
    columns = [observed[c].to_numpy() for c in ("chrom", "start", "end", "strand")]
    rows = list(zip(*columns, strict=True))

    resolved: dict[tuple, tuple] = {}
    for key in rows:
        if key in resolved:
            continue
        chrom, start, end, strand = key
        d, a = donor_acceptor(start, end, strand)
        gene = (
            gene_of_junction.get((chrom, start, end, strand))
            or gene_of_site.get((chrom, d, strand))
            or gene_of_site.get((chrom, a, strand))
        )
        resolved[key] = (
            classify_one(chrom, start, end, strand, junctions, donors, acceptors),
            gene[0] if gene else None,
            gene[1] if gene else None,
        )

    answers = [resolved[key] for key in rows]
    out = observed.copy()
    out["sclass"] = pd.Categorical([a[0] for a in answers], categories=CLASSES)
    out["gene_id"] = [a[1] for a in answers]
    out["gene_name"] = [a[2] for a in answers]
    out["is_novel"] = out["sclass"] != "annotated"
    return out


def annotation_summary(annotated: pd.DataFrame) -> pd.DataFrame:
    """Counts of unique junctions per class (deduplicated across samples)."""
    uniq = annotated.drop_duplicates(subset=["chrom", "start", "end", "strand"])
    counts = uniq["sclass"].value_counts().reindex(CLASSES, fill_value=0)
    return counts.rename_axis("sclass").reset_index(name="n_junctions")
