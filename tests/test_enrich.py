import numpy as np
import pandas as pd

from splicescope.enrich import enrich_differential, over_representation
from splicescope.io import read_gmt


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
