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


def test_the_odds_table_is_built_once_for_the_whole_collection(monkeypatch):
    """The p/(1-p) of every background gene does not depend on which set is being
    tested, so it must be computed once, not once per set. Rebuilding it per set is
    what made the weighted path cost one division per background gene per gene set:
    on a human background and an MSigDB-sized collection, 35 s against 2.5 s."""
    from splicescope import enrich

    calls = []
    original = enrich._odds_table
    monkeypatch.setattr(
        enrich, "_odds_table", lambda bg, prop: (calls.append(len(bg)), original(bg, prop))[1]
    )

    genes = [f"G{i:04d}" for i in range(300)]
    weights = {g: float(i % 17 + 1) for i, g in enumerate(genes)}
    hits = genes[:60]
    sets = {f"S{j}": genes[j : j + 40] for j in range(120)}

    result = enrich.over_representation(hits, genes, sets, weights=weights)
    assert not result.empty
    assert calls == [300], f"the odds table was rebuilt {len(calls)} times"


def test_taking_the_outside_sum_by_complement_changes_no_answer():
    """The outside sum is the total minus the inside one. That is only safe if the
    subtraction never has to produce a *small* number out of two near-equal ones, which
    is why the smaller side is always the one summed outright."""
    from splicescope.enrich import _bias_odds, _odds_table, selection_propensity

    rng = np.random.default_rng(11)
    genes = [f"G{i:04d}" for i in range(600)]
    units = np.maximum(1, rng.lognormal(1.5, 1.0, size=600).astype(int))
    weights = dict(zip(genes, units.astype(float), strict=True))
    hits = [g for g, u in zip(genes, units, strict=True) if rng.random() < 1 - 0.97**u]
    propensity = selection_propensity(genes, hits, weights)
    bg = set(genes)
    table = _odds_table(bg, propensity)

    def summed_both_ways(in_set):
        inside = sum(propensity[g] / (1 - propensity[g]) for g in in_set) / len(in_set)
        rest = bg - in_set
        beyond = sum(propensity[g] / (1 - propensity[g]) for g in rest) / len(rest)
        return inside / beyond

    # one gene, half the background, and the case a complement is worst at: all but one
    for size in (1, 2, 300, 598, 599):
        in_set = set(rng.choice(genes, size=size, replace=False))
        assert _bias_odds(in_set, bg, None, table) == pytest.approx(
            summed_both_ways(in_set), rel=1e-9
        ), f"the two disagree for a set of {size}"


def test_a_gene_set_member_missing_from_the_background_is_not_a_crash():
    """Only genes that were tested can be drawn, so a set member outside the background
    is simply not in the set for this purpose. Indexing the propensity by it raised
    KeyError before."""
    from splicescope.enrich import _bias_odds, _odds_table

    propensity = {"A": 0.4, "B": 0.2, "C": 0.2}
    table = _odds_table(set(propensity), propensity)
    assert _bias_odds({"A", "ZZZ"}, set(propensity), None, table) == pytest.approx(
        _bias_odds({"A"}, set(propensity), None, table)
    )
    assert _bias_odds({"ZZZ"}, set(propensity), None, table) is None


def test_two_spellings_of_one_gene_have_their_opportunity_added_up(monkeypatch):
    """A merged annotation can carry one accession at two versions. The background
    already counts those as one gene, so its opportunity is the sum of both. Keeping
    whichever the dictionary yielded last under-weights exactly the genes that are
    split, and can drop one into a bin of far less active genes."""
    from splicescope import enrich

    seen = {}
    original = enrich.selection_propensity
    monkeypatch.setattr(
        enrich,
        "selection_propensity",
        lambda bg, hits, w, **kw: (seen.update(w), original(bg, hits, w, **kw))[1],
    )

    genes = [f"ENSG{i:011d}.1" for i in range(30)]
    weights = {g: 1.0 for g in genes}
    weights["ENSG00000000000.1"] = 20.0
    weights["ENSG00000000000.7"] = 30.0  # the same gene, a second version

    enrich.over_representation(genes[:10], genes, {"S": genes[:8]}, weights=weights)
    assert seen["ENSG00000000000"] == 50.0, "the two spellings were not added up"


def test_the_weighted_test_reduces_to_the_plain_one_when_no_gene_is_favoured():
    """Wallenius' distribution at odds 1 *is* the hypergeometric, so giving every gene
    the same opportunity must reproduce the uncorrected p-value exactly. That identity
    is also the only thing that pins the tail convention on the weighted branch: the
    plain branch had a test for using `sf(k-1)` rather than `sf(k)`, and the weighted
    one did not."""
    genes = [f"G{i:04d}" for i in range(120)]
    hits = genes[:30]
    sets = {"S1": genes[:20], "S2": genes[10:40], "S3": genes[60:100]}

    plain = over_representation(hits, genes, sets).set_index("term")
    flat = over_representation(hits, genes, sets, weights=dict.fromkeys(genes, 3.0))
    flat = flat.set_index("term")

    assert set(flat.index) == set(plain.index)
    for term in plain.index:
        assert flat.loc[term, "bias_odds"] == pytest.approx(1.0), term
        assert flat.loc[term, "pvalue"] == pytest.approx(plain.loc[term, "pvalue"], rel=1e-6), term


def test_gene_sets_that_name_no_tested_gene_say_so_rather_than_finding_nothing():
    """An empty enrichment reads as "no pathway is enriched" — a conclusion. When no
    identifier is shared the truth is that the question was never asked, and gene-set
    files are keyed by symbols about as often as by accessions, so the mistake is easy.
    The CLI has warned since 0.8.1; the library returned an empty frame in silence."""
    import warnings

    diff = pd.DataFrame(
        {
            "chrom": ["chr1"] * 6, "start": range(6), "end": range(100, 106),
            "strand": ["+"] * 6,
            "gene_id": [f"ENSG{i:011d}" for i in range(6)],
            "gene_name": [f"GENE{i}" for i in range(6)],
            "qvalue": [1e-8] * 3 + [0.9] * 3,
            "delta_psi": [0.4] * 3 + [0.01] * 3,
            "abs_delta_psi": [0.4] * 3 + [0.01] * 3,
        }
    )

    with pytest.warns(UserWarning, match="no gene set shares an identifier"):
        empty = enrich_differential(diff, {"S": ["YFL001C", "YFL002C", "YFL003C"]})
    assert empty.empty

    # a collection that does name tested genes says nothing, and neither does no
    # collection at all — there is no mistake to report in either
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert len(enrich_differential(diff, {"S": ["GENE0", "GENE1", "GENE2", "GENE4"]})) == 1
        assert enrich_differential(diff, {}).empty


def test_the_warning_names_both_spellings_so_the_mismatch_is_visible():
    import warnings

    diff = pd.DataFrame(
        {
            "chrom": ["chr1"] * 4, "start": range(4), "end": range(100, 104),
            "strand": ["+"] * 4,
            "gene_id": [f"ENSG{i:011d}" for i in range(4)],
            "qvalue": [1e-8] * 2 + [0.9] * 2, "delta_psi": [0.4] * 2 + [0.01] * 2,
            "abs_delta_psi": [0.4] * 2 + [0.01] * 2,
        }
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        enrich_differential(diff, {"S": ["TARDBP", "STMN2", "UNC13A"]})
    message = str(caught[0].message)
    assert "STMN2" in message or "TARDBP" in message, "what the gene sets look like"
    assert "ENSG" in message, "what the annotation looks like"


def test_a_result_with_no_rows_still_has_the_columns_a_result_has():
    """Found by a property test, which could not even index the empty frame it got back.

    `over_representation` skips a gene set that contains no hit, so a collection where
    none of them overlaps left `records` empty and `from_records([])` returned a frame
    with no columns at all. The M == 0 / N == 0 path next to it already returned the
    full schema, so the shape of an empty result depended on *which* way it came out
    empty: selecting a column raised KeyError on one and worked on the other, and
    `to_csv` wrote a headerless file for a run that tested sets and enriched none.
    """
    from splicescope.enrich import RESULT_COLUMNS

    background = [f"G{i:04d}" for i in range(8)]
    hits = background[:3]

    populated = over_representation(hits, background, {"S": background[:4]})
    # every set is real and sized, none of them contains a hit
    no_overlap = over_representation(hits, background, {"S": background[4:6]})
    # nothing is a hit at all — the other empty path
    no_hits = over_representation([], background, {"S": background[:4]})

    assert len(populated) == 1 and no_overlap.empty and no_hits.empty
    for empty in (no_overlap, no_hits):
        assert list(empty.columns) == list(populated.columns) == RESULT_COLUMNS
        # the things a caller does next, which used to raise on this frame
        assert empty["term"].empty and empty[["pvalue", "qvalue", "genes"]].empty
        assert empty.to_csv(index=False).splitlines()[0] == (
            populated.to_csv(index=False).splitlines()[0]
        ), "same header, so a run that found nothing still writes a readable file"
