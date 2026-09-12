import numpy as np
import pandas as pd
import pytest
from scipy import stats

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


def test_sequence_names_ending_in_digits_are_not_de_versioned():
    """C. elegans sequence names are the standard identifier for most worm genes and are
    exactly `.<digits>` shaped. Stripping that collapsed whole gene families into one,
    shrinking the background, the set size, the overlap and the p-value with them."""
    from splicescope.enrich import normalize_gene_id

    worm = ["C42D8.1", "C42D8.2", "C42D8.3", "Y110A7A.10", "ZK1067.1", "F14H3.4"]
    assert len({normalize_gene_id(g) for g in worm}) == len(worm)


def test_accessions_are_still_de_versioned():
    from splicescope.enrich import normalize_gene_id

    assert normalize_gene_id("ENSG00000141510.16") == "ENSG00000141510"
    assert normalize_gene_id("ENSMUSG00000017146.8") == "ENSMUSG00000017146"
    assert normalize_gene_id("NM_000546.6") == "NM_000546"
    assert normalize_gene_id("XP_011527858.1") == "XP_011527858"


def test_a_worm_gene_family_keeps_every_member_in_the_enrichment():
    """Each collapsed member was silently dropped from the background and from the
    reported leading-edge gene list."""
    from splicescope.enrich import over_representation

    family = ["C42D8.1", "C42D8.2", "C42D8.3"]
    background = family + [f"gene-{i}" for i in range(97)]
    result = over_representation(family, background, {"SET": family + ["gene-0", "gene-1"]})

    assert len(result) == 1
    row = result.iloc[0]
    assert row["n_background"] == 100
    assert row["set_size"] == 5
    assert row["overlap"] == 3
    assert set(row["genes"].split(",")) == set(family)


# --- the gene-opportunity bias ---------------------------------------------------


def _size_driven_universe(seed=0, n_genes=800):
    """Genes whose only distinguishing feature is how many units were tested in them,
    and hits drawn purely in proportion to that. There is no biology here at all."""
    rng = np.random.default_rng(seed)
    genes = [f"G{i:04d}" for i in range(n_genes)]
    units = rng.integers(1, 60, size=n_genes)
    hits = [
        g
        for g, u in zip(genes, units, strict=True)
        if rng.random() < u / units.max() * 0.5
    ]
    weights = dict(zip(genes, units.astype(float), strict=True))
    largest = [g for _, g in sorted(zip(units, genes, strict=True), reverse=True)]
    return genes, hits, weights, largest


def test_unweighted_over_representation_is_still_the_hypergeometric():
    """Wallenius at odds 1 is the hypergeometric; omitting weights must change nothing."""
    genes, hits, _, largest = _size_driven_universe()
    result = over_representation(hits, genes, {"SET": largest[:150]})
    row = result.iloc[0]
    expected = stats.hypergeom.sf(
        row["overlap"] - 1, row["n_background"], row["set_size"], row["n_hits"]
    )
    assert row["pvalue"] == pytest.approx(expected)
    assert row["bias_odds"] == 1.0


def test_weighting_by_opportunity_defuses_a_pure_size_effect():
    """A gene set that is only "the largest genes" is a textbook false positive for ORA
    on splicing data: a long, many-exon gene has many more chances to contain a
    significant junction than a two-exon one."""
    genes, hits, weights, largest = _size_driven_universe()
    gene_sets = {"BIG_GENES": largest[:150]}

    plain = over_representation(hits, genes, gene_sets).iloc[0]
    weighted = over_representation(hits, genes, gene_sets, weights=weights).iloc[0]

    assert plain["pvalue"] < 1e-6, "the uncorrected test should be fooled outright"
    assert weighted["pvalue"] > plain["pvalue"] * 1e6
    assert weighted["bias_odds"] > 1.5, "the set should be flagged as opportunity-rich"


def test_the_correction_removes_the_false_positives():
    """Measured rather than asserted: across biologically null but size-biased sets the
    plain test calls essentially all of them, the weighted one calls almost none."""
    plain_p, weighted_p = [], []
    for seed in range(4):
        genes, hits, weights, largest = _size_driven_universe(seed=seed)
        rng = np.random.default_rng(100 + seed)
        gene_sets = {
            f"S{s}": list(rng.choice(largest[: 150 + s * 20], size=100, replace=False))
            for s in range(10)
        }
        plain_p.append(over_representation(hits, genes, gene_sets)["pvalue"].to_numpy())
        weighted_p.append(
            over_representation(hits, genes, gene_sets, weights=weights)["pvalue"].to_numpy()
        )
    plain_p = np.concatenate(plain_p)
    weighted_p = np.concatenate(weighted_p)

    assert (plain_p <= 0.05).mean() > 0.9, "the uncorrected null should fail almost always"
    assert (weighted_p <= 0.05).mean() < 0.10, "the corrected null should be near nominal"
    assert (weighted_p <= 1e-6).mean() == 0.0, "no confident false call should survive"


def test_propensity_rises_with_opportunity():
    from splicescope.enrich import selection_propensity

    genes, hits, weights, largest = _size_driven_universe()
    propensity = selection_propensity(genes, hits, weights)
    smallest_third = [g for g in largest[-200:]]
    largest_third = [g for g in largest[:200]]
    assert (
        sum(propensity[g] for g in largest_third)
        > sum(propensity[g] for g in smallest_third) * 2
    )


def test_enrich_differential_weights_by_units_per_gene_by_default():
    """The opportunity measure is sitting right there in the differential table: one row
    per tested unit."""
    rows = []
    for gene, n_units in (("BIG", 40), ("SMALL", 2)):
        for i in range(n_units):
            rows.append(
                {
                    "gene_id": gene,
                    "qvalue": 0.001 if i < 2 else 0.9,
                    "delta_psi": 0.5 if i < 2 else 0.01,
                    "abs_delta_psi": 0.5 if i < 2 else 0.01,
                }
            )
    for i in range(60):                      # filler genes, one unit each, not hits
        rows.append(
            {"gene_id": f"F{i}", "qvalue": 0.9, "delta_psi": 0.01, "abs_delta_psi": 0.01}
        )
    table = pd.DataFrame(rows)
    sets = {"SET": ["BIG", "SMALL", "F0", "F1"]}

    weighted = enrich_differential(table, sets)
    unweighted = enrich_differential(table, sets, weight_by_units=False)
    assert not weighted.empty and not unweighted.empty
    assert weighted.iloc[0]["bias_odds"] > 1.0
    assert unweighted.iloc[0]["bias_odds"] == 1.0


def test_the_odds_are_averaged_not_the_probabilities():
    """Wallenius' ω is a ratio of sampling weights, so the quantity to average is
    p/(1-p). goseq averages the probabilities instead; the two agree only while p is
    small, and a selection propensity here reaches 0.5."""
    from splicescope.enrich import _bias_odds

    # two genes inside at p = 0.5, two outside at p = 0.25
    propensity = {"A": 0.5, "B": 0.5, "C": 0.25, "D": 0.25}
    odds = _bias_odds({"A", "B"}, set(propensity), propensity)

    ratio_of_probabilities = 0.5 / 0.25                      # what goseq would use
    ratio_of_odds = (0.5 / 0.5) / (0.25 / 0.75)              # what Wallenius asks for
    assert odds == pytest.approx(ratio_of_odds)
    assert odds == pytest.approx(3.0)
    assert odds != pytest.approx(ratio_of_probabilities)


def test_weighting_keeps_the_power_to_detect_a_real_enrichment():
    """A correction that removed the false positives by refusing to call anything would
    be no use. A genuinely enriched set must survive it."""
    rng = np.random.default_rng(7)
    n_genes = 800
    genes = [f"G{i:04d}" for i in range(n_genes)]
    units = rng.integers(1, 60, size=n_genes)
    true_set = set(rng.choice(genes, size=120, replace=False))
    hits = [
        g
        for g, u in zip(genes, units, strict=True)
        if rng.random() < min(0.95, u / units.max() * 0.5 + (0.25 if g in true_set else 0.0))
    ]
    weights = dict(zip(genes, units.astype(float), strict=True))

    result = over_representation(hits, genes, {"REAL": sorted(true_set)}, weights=weights)
    assert len(result) == 1
    assert result.iloc[0]["pvalue"] < 0.01, "a real enrichment must still be found"


def test_a_bin_can_never_assert_certainty():
    """A bin whose members all are (or all are not) hits gave a rate of exactly 1 or 0.
    `_bias_odds` averages p/(1-p), so a rate of 1 contributed 1e6 after clipping and
    decided the mean by itself — in one direction for the set that held the gene, and in
    the other for every set that did not."""
    from splicescope.enrich import selection_propensity

    rng = np.random.default_rng(21)
    genes = [f"G{i:04d}" for i in range(2000)]
    units = np.maximum(1, rng.lognormal(1.6, 1.0, size=2000).astype(int))
    hits = [
        g
        for g, u in zip(genes, units, strict=True)
        if rng.random() < 1 - (1 - 0.05) ** u
    ]
    weights = dict(zip(genes, units.astype(float), strict=True))
    propensity = selection_propensity(genes, hits, weights)

    assert 0.0 < min(propensity.values())
    assert max(propensity.values()) < 1.0
    # smoothing is by bin size, so even an all-hit bin of 7 cannot exceed 7.5/8
    assert max(propensity.values()) <= 0.95


def test_one_high_opportunity_gene_cannot_destroy_a_real_enrichment():
    """The user-facing consequence of a saturated bin: adding a single gene — one that
    made the set *more* enriched — took its p-value from 4e-15 to 0.22."""
    rng = np.random.default_rng(21)
    genes = [f"G{i:04d}" for i in range(2000)]
    units = np.maximum(1, rng.lognormal(1.6, 1.0, size=2000).astype(int))
    hits = [
        g
        for g, u in zip(genes, units, strict=True)
        if rng.random() < 1 - (1 - 0.05) ** u
    ]
    weights = dict(zip(genes, units.astype(float), strict=True))
    hit_set = set(hits)
    biggest_hit = max(hit_set, key=lambda g: weights[g])

    members = list(rng.choice(sorted(hit_set - {biggest_hit}), size=60, replace=False))
    members += list(
        rng.choice([g for g in genes if g not in hit_set], size=40, replace=False)
    )
    with_it = members[:-1] + [biggest_hit]

    before = over_representation(hits, genes, {"S": members}, weights=weights).iloc[0]
    after = over_representation(hits, genes, {"S": with_it}, weights=weights).iloc[0]

    assert after["overlap"] > before["overlap"], "the swap must add a hit"
    assert before["pvalue"] < 1e-6
    # the weighted null may legitimately soften the call, but not annihilate it
    assert after["pvalue"] < 1e-3, (
        f"one gene moved p from {before.pvalue:.1e} to {after.pvalue:.1e}"
    )
