"""`validation/bam2sj.awk`, run on hand-written SAM records.

The script turns alignments into the `SJ.out.tab` that the rest of the package reads, so
its output has to look like STAR's: one row per intron. It used to key a row on the
XS:A strand tag as well as the coordinates, and HISAT2 tags some of an intron's reads
and not others, so one intron came out as two rows. `compute_psi` refuses such a
table, and before it did, each row held only part of the intron's reads.

The SAM text is piped straight into awk, so the test needs neither samtools nor a BAM.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from splicescope.annotate import annotate_junctions
from splicescope.io import read_star_sj
from splicescope.quantify import compute_psi

AWK_SCRIPT = Path(__file__).resolve().parent.parent / "validation" / "bam2sj.awk"

pytestmark = pytest.mark.skipif(shutil.which("awk") is None, reason="needs awk on PATH")


def _read(qname, chrom, pos, cigar, *tags, flag=0):
    """One SAM record; sequence and qualities are irrelevant to the script."""
    fields = [qname, str(flag), chrom, str(pos), "60", cigar, "*", "0", "0", "*", "*"]
    return "\t".join(fields + list(tags))


# Each block below is one intron; the comment gives its 1-based inclusive coordinates.
SAM = "\n".join(
    [
        "@HD\tVN:1.6\tSO:coordinate",
        "@SQ\tSN:chr1\tLN:100000",
        "@SQ\tSN:chr2\tLN:100000",
        # chr1:110-159, one read tagged + and two with no XS:A at all. The old script
        # wrote this intron twice: strand 1 with one read, strand 0 with two.
        _read("a1", "chr1", 100, "10M50N20M", "NH:i:1", "XS:A:+"),
        _read("a2", "chr1", 95, "15M50N25M", "NH:i:1"),
        _read("a3", "chr1", 104, "6M50N30M", "NH:i:1"),
        # chr1:320-399, tagged + twice and - once: the majority is +
        _read("b1", "chr1", 300, "20M80N20M", "NH:i:1", "XS:A:+"),
        _read("b2", "chr1", 300, "20M80N20M", "NH:i:1", "XS:A:+"),
        _read("b3", "chr1", 300, "20M80N20M", "NH:i:1", "XS:A:-"),
        # chr1:620-659, one + and one -, and an untagged read that does not vote: a tie
        _read("c1", "chr1", 600, "20M40N20M", "NH:i:1", "XS:A:+"),
        _read("c2", "chr1", 600, "20M40N20M", "NH:i:1", "XS:A:-"),
        _read("c3", "chr1", 600, "20M40N20M", "NH:i:1"),
        # chr1:925-994, - strand, one unique read and two multi-mappers
        _read("d1", "chr1", 900, "25M70N25M", "NH:i:1", "XS:A:-"),
        _read("d2", "chr1", 900, "25M70N25M", "NH:i:3", "XS:A:-"),
        _read("d3", "chr1", 900, "25M70N25M", "NH:i:2"),
        # chr2:1020-1119 and chr2:1150-1349 from one read with two N operations, and a
        # second read across the first intron only
        _read("e1", "chr2", 1000, "20M100N30M200N25M", "NH:i:1", "XS:A:-"),
        _read("e2", "chr2", 1005, "15M100N12M", "NH:i:1", "XS:A:-"),
        # chr2:5020-5079, reached through a soft clip, an insertion and a deletion. The
        # overhangs are 8 (8M | 12M), 3 (8M | 3M) and 9 (30M | 9M), so the column is 9.
        _read("f1", "chr2", 5000, "5S6M1I4M2D8M60N12M3S", "NH:i:1", "XS:A:+"),
        _read("f2", "chr2", 5012, "8M60N3M", "NH:i:1", "XS:A:+"),
        _read("f3", "chr2", 4990, "30M60N9M", "NH:i:1"),
        # neither spliced nor mapped, so neither may appear
        _read("u1", "chr1", 50, "50M", "NH:i:1", "XS:A:+"),
        _read("u2", "*", 0, "*", flag=4),
    ]
) + "\n"

#: (chrom, start, end) -> (strand_code, motif, annotated, n_unique, n_multi, max_overhang)
EXPECTED = {
    ("chr1", 110, 159): (1, 0, 0, 3, 0, 15),
    ("chr1", 320, 399): (1, 0, 0, 3, 0, 20),
    ("chr1", 620, 659): (0, 0, 0, 3, 0, 20),
    ("chr1", 925, 994): (2, 0, 0, 1, 2, 25),
    ("chr2", 1020, 1119): (2, 0, 0, 2, 0, 20),
    ("chr2", 1150, 1349): (2, 0, 0, 1, 0, 25),
    ("chr2", 5020, 5079): (1, 0, 0, 3, 0, 9),
}


def _run_awk(sam: str = SAM) -> list[list[str]]:
    result = subprocess.run(
        ["awk", "-f", str(AWK_SCRIPT)], input=sam, capture_output=True, text=True, check=True
    )
    return [line.split("\t") for line in result.stdout.splitlines()]


def _by_intron(rows: list[list[str]]) -> dict[tuple, list[tuple]]:
    out: dict[tuple, list[tuple]] = {}
    for row in rows:
        chrom, start, end, *rest = row
        out.setdefault((chrom, int(start), int(end)), []).append(tuple(int(v) for v in rest))
    return out


def test_every_intron_is_one_row():
    """The regression itself: an intron crossed by tagged and untagged reads, or by reads
    tagged + and -, is still one intron."""
    rows = _run_awk()
    assert all(len(row) == 9 for row in rows)
    repeated = {k: v for k, v in _by_intron(rows).items() if len(v) > 1}
    assert not repeated, f"introns written more than once: {repeated}"


def test_tagged_and_untagged_reads_are_summed_under_the_tagged_strand():
    rows = _by_intron(_run_awk())
    assert rows[("chr1", 110, 159)] == [EXPECTED[("chr1", 110, 159)]]


def test_conflicting_tags_take_the_majority_and_a_tie_is_undefined():
    rows = _by_intron(_run_awk())
    assert rows[("chr1", 320, 399)] == [EXPECTED[("chr1", 320, 399)]]
    assert rows[("chr1", 620, 659)] == [EXPECTED[("chr1", 620, 659)]]


def test_multimappers_are_counted_apart_from_unique_reads():
    rows = _by_intron(_run_awk())
    assert rows[("chr1", 925, 994)] == [EXPECTED[("chr1", 925, 994)]]


def test_a_read_with_two_introns_counts_towards_both():
    rows = _by_intron(_run_awk())
    assert rows[("chr2", 1020, 1119)] == [EXPECTED[("chr2", 1020, 1119)]]
    assert rows[("chr2", 1150, 1349)] == [EXPECTED[("chr2", 1150, 1349)]]


def test_coordinates_and_overhang_through_clips_insertions_and_deletions():
    rows = _by_intron(_run_awk())
    assert rows[("chr2", 5020, 5079)] == [EXPECTED[("chr2", 5020, 5079)]]


def test_the_whole_table():
    """Nothing beyond the seven introns: headers, an unspliced read and an unmapped one
    contribute no row."""
    assert {k: v[0] for k, v in _by_intron(_run_awk()).items()} == EXPECTED


def test_a_read_without_an_nh_tag_counts_as_unique():
    rows = _by_intron(_run_awk(_read("x", "chr3", 10, "5M20N5M") + "\n"))
    assert rows == {("chr3", 15, 34): [(0, 0, 0, 1, 0, 5)]}


def test_the_table_is_read_and_passes_the_repeat_guard(tmp_path):
    """The output goes through splicescope's own reader, annotation and Ψ, which is where
    the doubled rows were refused.

    The tie at chr1:620-659 is written as strand 0, read as ".", and placed on the strand
    the annotation gives it. Under the old script chr1:110-159 had a strand-0 row beside
    its + row, and the annotation turned it into a second + row of the same junction.
    """
    sj = tmp_path / "S1.SJ.out.tab"
    sj.write_text("".join("\t".join(row) + "\n" for row in _run_awk()))

    observed = read_star_sj(sj)
    known = pd.DataFrame(
        [
            ("chr1", 110, 159, "+", "G1", "GENE1"),
            ("chr1", 620, 659, "+", "G1", "GENE1"),
            ("chr2", 1020, 1119, "-", "G2", "GENE2"),
        ],
        columns=["chrom", "start", "end", "strand", "gene_id", "gene_name"],
    )
    annotated = annotate_junctions(observed, known)
    psi = compute_psi(annotated, min_reads=1)  # raises on a junction listed twice

    assert set(observed["sample"]) == {"S1"}
    as_read = observed.set_index(["chrom", "start", "end"])["strand"]
    assert as_read[("chr1", 110, 159)] == "+"
    assert as_read[("chr1", 620, 659)] == "."
    assert as_read[("chr1", 925, 994)] == "-"

    rows = psi[psi["count"] > 0].set_index(["chrom", "start", "end"])
    assert len(rows) == len(EXPECTED)
    for intron in [("chr1", 110, 159), ("chr1", 620, 659)]:
        strand, count, sclass = rows.loc[intron, ["strand", "count", "sclass"]]
        assert (strand, count, sclass) == ("+", 3, "annotated"), intron
