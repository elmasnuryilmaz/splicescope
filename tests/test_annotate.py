import pandas as pd
import pytest

from splicescope.annotate import annotate_junctions, classify_one
from splicescope.io import donor_acceptor


def _known():
    # one known intron on + strand: chr1:200-300
    return pd.DataFrame(
        [("chr1", 200, 300, "+", "g0")], columns=["chrom", "start", "end", "strand", "gene_id"]
    )


def test_donor_acceptor_strand():
    assert donor_acceptor(200, 300, "+") == (200, 300)
    assert donor_acceptor(200, 300, "-") == (300, 200)


def test_classify_annotated_and_novel():
    junctions = {("chr1", 200, 300, "+")}
    donors = {("chr1", 200, "+")}
    acceptors = {("chr1", 300, "+")}

    assert classify_one("chr1", 200, 300, "+", junctions, donors, acceptors) == "annotated"
    # known donor, novel acceptor
    assert classify_one("chr1", 200, 260, "+", junctions, donors, acceptors) == "novel_acceptor"
    # novel donor, known acceptor
    assert classify_one("chr1", 250, 300, "+", junctions, donors, acceptors) == "novel_donor"
    # both novel -> cryptic
    assert classify_one("chr1", 250, 260, "+", junctions, donors, acceptors) == "cryptic"
    # both known but not a known pair -> novel_combination
    donors2 = donors | {("chr1", 500, "+")}
    acceptors2 = acceptors | {("chr1", 600, "+")}
    assert (
        classify_one("chr1", 500, 300, "+", {("chr1", 200, 300, "+")}, donors2, acceptors2)
        == "novel_combination"
    )


def test_annotate_junctions_adds_columns():
    known = _known()
    obs = pd.DataFrame(
        [
            ("chr1", 200, 300, "+", 100),
            ("chr1", 200, 260, "+", 5),
        ],
        columns=["chrom", "start", "end", "strand", "count"],
    )
    out = annotate_junctions(obs, known)
    assert list(out["sclass"]) == ["annotated", "novel_acceptor"]
    assert out.loc[0, "gene_id"] == "g0"
    assert out.loc[0, "is_novel"] is False or out.loc[0, "is_novel"] == False  # noqa: E712


def test_strand_undefined_junctions_are_placed_against_the_annotation():
    """STAR writes strand code 0 whenever the intron motif does not reveal a strand.
    Such a junction can never equal a stranded annotation, so it was classified cryptic
    however ordinary it was, lost its gene, and — alone at its own (chrom, pos, '.')
    site — was handed Ψ = 1.0 in every sample."""
    from splicescope.quantify import compute_psi

    known = pd.DataFrame(
        [
            dict(chrom="chr1", start=1000, end=2000, strand="+", gene_id="G1"),
            dict(chrom="chr1", start=1000, end=3000, strand="+", gene_id="G1"),
        ]
    )
    observed = pd.DataFrame(
        [
            dict(chrom="chr1", start=1000, end=2000, strand=".", sample="s1", count=40),
            dict(chrom="chr1", start=1000, end=2500, strand=".", sample="s1", count=10),
            dict(chrom="chr1", start=7000, end=7500, strand=".", sample="s1", count=5),
        ]
    )
    out = annotate_junctions(observed, known)

    # the exact intron, and a junction sharing its donor, are both placed on '+'
    assert list(out["strand"]) == ["+", "+", "."]
    assert list(out["sclass"]) == ["annotated", "novel_acceptor", "cryptic"]
    assert list(out["gene_id"])[:2] == ["G1", "G1"]

    # and they now share a donor, so Ψ is a real fraction rather than 1.0 each
    psi = compute_psi(out, min_reads=5)
    assert psi.loc[0, "psi_donor"] == pytest.approx(0.8)
    assert psi.loc[1, "psi_donor"] == pytest.approx(0.2)


def test_an_unplaceable_junction_keeps_its_undefined_strand():
    """Guessing a strand the annotation cannot support would be worse than saying so."""
    known = pd.DataFrame([dict(chrom="chr1", start=1000, end=2000, strand="+", gene_id="G1")])
    observed = pd.DataFrame(
        [dict(chrom="chr2", start=500, end=900, strand=".", sample="s1", count=3)]
    )
    assert annotate_junctions(observed, known).loc[0, "strand"] == "."


def test_an_annotation_with_no_junctions_is_refused():
    """Every class but `cryptic` is defined by agreeing with the annotation, so an
    empty one does not fail — it reports a genome of novel splicing. A GTF carrying only
    gene or transcript records produces exactly that, because introns are derived from
    exons."""
    import pytest

    observed = pd.DataFrame(
        [dict(chrom="chr1", start=200, end=399, strand="+", sample="s1", count=40)]
    )
    empty = pd.DataFrame(columns=["chrom", "start", "end", "strand", "gene_id"])
    with pytest.raises(ValueError, match="contains no junctions"):
        annotate_junctions(observed, empty)


@pytest.mark.parametrize(
    ("junction_chrom", "annotation_chrom"),
    [("chr1", "1"), ("1", "chr1")],
)
def test_chromosome_names_from_different_sources_are_refused(junction_chrom, annotation_chrom):
    """GENCODE writes `chr1` where Ensembl writes `1`. Mixing them matches nothing, and
    a run that calls every junction cryptic is the most exciting wrong answer this tool
    can give."""
    observed = pd.DataFrame(
        [dict(chrom=junction_chrom, start=200, end=399, strand="+", sample="s1", count=40)]
    )
    known = pd.DataFrame(
        [dict(chrom=annotation_chrom, start=200, end=399, strand="+", gene_id="G1")]
    )
    with pytest.raises(ValueError, match="name the same chromosomes differently"):
        annotate_junctions(observed, known)


def test_a_junction_on_a_contig_the_annotation_does_not_cover_is_ordinary():
    """Scaffolds, decoys and a chromosome left out of a small analysis all produce
    junctions the annotation says nothing about. That is not a mistake, so it is not an
    error — only the naming mismatch is, because it is recognisable."""
    observed = pd.DataFrame(
        [
            dict(chrom="chrUn_KI270742v1", start=200, end=399, strand="+",
                 sample="s1", count=40),
            dict(chrom="chr1", start=200, end=399, strand="+", sample="s1", count=40),
        ]
    )
    known = pd.DataFrame([dict(chrom="chr1", start=200, end=399, strand="+", gene_id="G1")])
    out = annotate_junctions(observed, known)
    assert list(out["sclass"]) == ["cryptic", "annotated"]


def test_every_junction_is_classified_by_the_rule_the_methods_table_states():
    """METHODS §3's five rules, written out here and applied to simulated data.

    The existing tests build a junction per class by hand, which checks each rule once on
    an input chosen to exercise it. This applies all five to whatever the simulator emits,
    from the annotation alone — a second reading of the table rather than a second call
    into the code that implements it, so a rule and its documentation cannot drift apart
    without one of them failing.
    """
    from splicescope.annotate import annotate_junctions
    from splicescope.io import donor_acceptor
    from splicescope.simulate import simulate_dataset

    #: (junction known, donor known, acceptor known) -> the class the table names
    rules = {
        (True, True, True): "annotated",
        (False, True, True): "novel_combination",
        (False, True, False): "novel_acceptor",
        (False, False, True): "novel_donor",
        (False, False, False): "cryptic",
    }

    seen = {}
    for seed in (5, 13):
        ds = simulate_dataset(n_genes=30, n_per_group=3, cryptic_fraction=0.7,
                              alt_ss_fraction=0.5, seed=seed)
        junctions, donors, acceptors = set(), set(), set()
        for row in ds.known.itertuples(index=False):
            junctions.add((row.chrom, row.start, row.end, row.strand))
            donor, acceptor = donor_acceptor(row.start, row.end, row.strand)
            donors.add((row.chrom, donor, row.strand))
            acceptors.add((row.chrom, acceptor, row.strand))

        annotated = annotate_junctions(ds.observed, ds.known)
        unique = annotated.drop_duplicates(subset=["chrom", "start", "end", "strand"])
        for row in unique.itertuples(index=False):
            donor, acceptor = donor_acceptor(row.start, row.end, row.strand)
            expected = rules[(
                (row.chrom, row.start, row.end, row.strand) in junctions,
                (row.chrom, donor, row.strand) in donors,
                (row.chrom, acceptor, row.strand) in acceptors,
            )]
            assert expected == row.sclass, (
                f"{row.chrom}:{row.start}-{row.end}{row.strand}: "
                f"the table says {expected}, the code says {row.sclass}"
            )
            seen[expected] = seen.get(expected, 0) + 1

    # the simulator makes no novel_combination — see METHODS §8 — so four of the five
    assert set(seen) == {"annotated", "novel_acceptor", "novel_donor", "cryptic"}
    assert min(seen.values()) >= 20, f"too few of some class to mean much: {seen}"
