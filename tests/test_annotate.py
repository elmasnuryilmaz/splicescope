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
