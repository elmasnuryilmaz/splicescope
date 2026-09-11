import numpy as np
import pandas as pd
import pytest

from splicescope.enrich import enrich_differential, over_representation
from splicescope.io import read_gmt

_GTF_LINE = (
    'chr1\tsrc\texon\t{start}\t{end}\t.\t+\t.\t'
    'gene_id "ENSG00000141510.16"; transcript_id "T1"; gene_name "TP53";\n'
)


@pytest.fixture
def gtf_with_symbols(tmp_path):
    """A two-exon GTF spelled the way GENCODE spells things: versioned id + symbol."""
    path = tmp_path / "symbols.gtf"
    path.write_text(_GTF_LINE.format(start=100, end=200) + _GTF_LINE.format(start=401, end=500))
    return path


def test_over_representation_math_and_ranking():
    background = [f"g{i}" for i in range(20)]
    hits = [f"g{i}" for i in range(8)]  # g0..g7
    gene_sets = {
        "enriched": ["g0", "g1", "g2", "g3", "g4"],   # entirely within hits
        "partial": ["g0", "g1", "g10", "g11", "g12"],  # 2/5 in hits
        "random": ["g10", "g11", "g12", "g13", "g14"],  # none in hits (dropped)
    }
    res = over_representation(hits, background, gene_sets)
    terms = set(res["term"])
    assert "enriched" in terms and "partial" in terms
    assert "random" not in terms  # zero overlap is not tested

    enr = res.set_index("term")
    # fold enrichment = (k/N) / (n/M)
    assert np.isclose(enr.loc["enriched", "fold_enrichment"], (5 / 8) / (5 / 20))
    assert np.isclose(enr.loc["partial", "fold_enrichment"], (2 / 8) / (5 / 20))
    # the fully-overlapping set is more significant than the partial one
    assert enr.loc["enriched", "qvalue"] <= enr.loc["partial", "qvalue"]
    assert (res["pvalue"] >= 0).all() and (res["pvalue"] <= 1).all()


def test_over_representation_empty_hits():
    res = over_representation([], ["g0", "g1"], {"s": ["g0"]})
    assert res.empty


def test_read_gmt_roundtrip(tmp_path):
    gmt = tmp_path / "sets.gmt"
    gmt.write_text(
        "PATHWAY_A\ta description\tg0\tg1\tg2\n"
        "PATHWAY_B\t\tg3\tg4\n"
        "junk_line_without_genes\n"
    )
    sets = read_gmt(gmt)
    assert sets == {"PATHWAY_A": ["g0", "g1", "g2"], "PATHWAY_B": ["g3", "g4"]}


def test_enrich_differential_end_to_end():
    # a toy differential table: g0..g4 significant, g10..g14 not
    rows = []
    for i in range(5):
        rows.append({"gene_id": f"g{i}", "qvalue": 0.001, "delta_psi": 0.3, "abs_delta_psi": 0.3})
    for i in range(10, 20):
        rows.append({"gene_id": f"g{i}", "qvalue": 0.9, "delta_psi": 0.0, "abs_delta_psi": 0.0})
    diff = pd.DataFrame(rows)
    gene_sets = {"hit_pathway": ["g0", "g1", "g2", "g3"], "other": ["g15", "g16", "g17"]}
    res = enrich_differential(diff, gene_sets, q=0.05, min_delta=0.1)
    assert res.iloc[0]["term"] == "hit_pathway"
    assert res.iloc[0]["fold_enrichment"] > 1


def _p53_table():
    """A differential table spelled the way a real GENCODE annotation spells things."""
    return pd.DataFrame(
        {
            "gene_id": [
                "ENSG00000141510.16", "ENSG00000012048.23",
                "ENSG00000139618.15", "ENSG00000146648.18",
                "ENSG00000171862.10", "ENSG00000105221.15",
            ],
            "gene_name": ["TP53", "BRCA1", "BRCA2", "EGFR", "PTEN", "AKT2"],
            "qvalue": [0.001, 0.002, 0.9, 0.9, 0.9, 0.9],
            "delta_psi": [0.5, 0.5, 0.01, 0.01, 0.01, 0.01],
            "abs_delta_psi": [0.5, 0.5, 0.01, 0.01, 0.01, 0.01],
        }
    )


def test_versioned_accessions_still_match_unversioned_gene_sets():
    """MSigDB's Ensembl GMTs are unversioned; a GTF's ids are not. Matching them
    verbatim meant every set had a zero-size background and nothing was ever tested."""
    from splicescope.enrich import enrich_differential

    sets = {"P53": ["ENSG00000141510", "ENSG00000012048", "ENSG00000139618"]}
    result = enrich_differential(_p53_table(), sets)
    assert len(result) == 1
    assert result.iloc[0]["set_size"] == 3
    assert result.iloc[0]["overlap"] == 2
    # the original spelling is reported back, not the normalised one
    assert "ENSG00000141510.16" in result.iloc[0]["genes"]


def test_symbol_gene_sets_match_through_the_gene_name_column():
    """Every MSigDB *.symbols.gmt, GO and KEGG set is keyed by symbol."""
    from splicescope.enrich import enrich_differential

    result = enrich_differential(_p53_table(), {"P53": ["TP53", "BRCA1", "BRCA2"]})
    assert len(result) == 1
    assert result.iloc[0]["overlap"] == 2
    assert set(result.iloc[0]["genes"].split(",")) == {"TP53", "BRCA1"}


def test_identifier_matching_is_case_insensitive():
    from splicescope.enrich import enrich_differential

    result = enrich_differential(_p53_table(), {"P53": ["tp53", "brca1", "brca2"]})
    assert len(result) == 1
    assert result.iloc[0]["overlap"] == 2


def test_the_column_that_matches_the_gene_sets_is_the_one_used():
    """With both columns present, picking the wrong one silently tests nothing."""
    from splicescope.enrich import enrich_differential

    table = _p53_table()
    by_symbol = enrich_differential(table, {"P53": ["TP53", "BRCA1", "BRCA2"]})
    by_accession = enrich_differential(
        table, {"P53": ["ENSG00000141510", "ENSG00000012048", "ENSG00000139618"]}
    )
    assert by_symbol.iloc[0]["overlap"] == by_accession.iloc[0]["overlap"] == 2


def test_gene_symbols_reach_the_differential_table_from_the_gtf(gtf_with_symbols):
    """The chain GTF -> annotate -> diff has to carry gene_name, or symbol sets have
    nothing to match against."""
    from splicescope.annotate import annotate_junctions
    from splicescope.io import read_gtf_junctions

    known = read_gtf_junctions(gtf_with_symbols)
    assert "gene_name" in known.columns
    assert set(known["gene_name"]) == {"TP53"}

    observed = pd.DataFrame(
        [dict(chrom="chr1", start=201, end=400, strand="+", sample="s1", count=10)]
    )
    annotated = annotate_junctions(observed, known)
    assert annotated.loc[0, "gene_name"] == "TP53"
    assert annotated.loc[0, "gene_id"] == "ENSG00000141510.16"
