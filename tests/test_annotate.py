import pandas as pd

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
