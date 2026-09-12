"""The consequence layer on the minus strand.

Every coordinate in this module means something different depending on the strand: the
donor and acceptor swap ends, exons are walked in decreasing genomic order, coding
length accumulates downwards, and each piece of sequence is reverse-complemented on the
way out. Half of any genome is on that strand.

It was barely tested. `test_consequence_pipeline.py` never used it at all, and of five
deliberate defects that the whole suite failed to notice, two were reachable only from
the minus strand: exon blocks left in genomic order rather than transcription order, and
the coding length before a position counted one base too many.

The main test here is a mirror. One gene is built on the plus strand, and a second is
its exact reverse complement: the chromosome is reverse-complemented, every coordinate
is reflected, and the strand is flipped. The two therefore describe the same transcript,
the same open reading frame and the same event, so every prediction has to agree. Any
asymmetry left in the arithmetic shows up as a disagreement, whatever the right answer
happens to be.
"""

from __future__ import annotations

import pandas as pd
import pytest

from splicescope.consequence import (
    GenomeFasta,
    Transcript,
    predict_consequence,
    predict_junction_consequence,
    reverse_complement,
)

CHROM_LENGTH = 1_200
LINE = 60
EXONS = [(101, 200), (401, 500), (701, 800)]  # introns at 201-400 and 501-700


def _forward_sequence() -> str:
    """A chromosome whose exons carry an open reading frame and whose introns do not."""
    bases = ["C"] * CHROM_LENGTH
    for start, end in EXONS:
        codons = "AAGGCTACC" * 12  # 108 nt of sense codons, no stop in any frame
        for offset, position in enumerate(range(start, end + 1)):
            bases[position - 1] = codons[offset]
    for start, end, filler in ((201, 400, "TAA"), (501, 700, "TAG")):
        for offset, position in enumerate(range(start, end + 1)):
            bases[position - 1] = filler[offset % 3]
    return "".join(bases)


def _reflect(start: int, end: int) -> tuple[int, int]:
    """The same interval on a reverse-complemented chromosome of the same length."""
    return CHROM_LENGTH - end + 1, CHROM_LENGTH - start + 1


def _write_fasta(tmp_path, sequences: dict[str, str]):
    fasta = tmp_path / "genome.fa"
    rows, offset = [], 0
    with open(fasta, "w") as handle:
        for name, seq in sequences.items():
            header = f">{name}\n"
            handle.write(header)
            offset += len(header)
            rows.append(f"{name}\t{len(seq)}\t{offset}\t{LINE}\t{LINE + 1}")
            for i in range(0, len(seq), LINE):
                handle.write(seq[i : i + LINE] + "\n")
                offset += LINE + 1
    fasta.with_suffix(".fa.fai").write_text("\n".join(rows) + "\n")
    return fasta


def _transcript(chrom: str, strand: str, exons: list[tuple[int, int]]) -> Transcript:
    return Transcript(
        transcript_id=f"T_{strand}",
        gene_id="G1",
        gene_name="GENE1",
        chrom=chrom,
        strand=strand,
        exons=list(exons),
        cds=list(exons),
    )


@pytest.fixture
def mirrored(tmp_path):
    """``(genome, plus_transcript, minus_transcript)`` describing one gene twice."""
    forward = _forward_sequence()
    fasta = _write_fasta(tmp_path, {"fwd": forward, "rev": reverse_complement(forward)})
    plus = _transcript("fwd", "+", EXONS)
    minus = _transcript("rev", "-", [_reflect(s, e) for s, e in EXONS])
    return GenomeFasta(fasta), plus, minus


#: Everything a Consequence claims. Compared field by field, so a disagreement names
#: itself rather than arriving as "the two objects differ".
FIELDS = (
    "consequence_class",
    "insert_length",
    "frame_offset",
    "frameshift",
    "ptc_offset",
    "distance_to_last_junction",
    "nmd_predicted",
)


def _same(forward, reverse, label):
    for name in FIELDS:
        assert getattr(forward, name) == getattr(reverse, name), f"{name} differs: {label}"


def test_the_two_strands_describe_the_same_reading_frame(mirrored):
    """The premise everything below rests on: reflecting the coordinates and reverse
    complementing the chromosome leaves the spliced coding sequence unchanged."""
    genome, plus, minus = mirrored
    with genome as fasta:
        forward_cds = "".join(fasta.fetch("fwd", s, e, "+") for s, e in plus.exons)
        reverse_cds = "".join(fasta.fetch("rev", s, e, "-") for s, e in minus.exons)
    assert forward_cds == reverse_cds
    assert len(forward_cds) == 300


@pytest.mark.parametrize(
    ("start", "end", "label"),
    [
        (201, 380, "acceptor shifted earlier, extending the downstream exon"),
        (201, 420, "acceptor shifted later, truncating the downstream exon"),
        (230, 400, "donor shifted later, extending the upstream exon"),
        (501, 640, "second intron, acceptor shifted earlier"),
        (501, 725, "second intron, acceptor shifted past the exon start"),
    ],
)
def test_a_shifted_splice_site_is_called_the_same_on_either_strand(mirrored, start, end, label):
    genome, plus, minus = mirrored
    with genome as fasta:
        forward = predict_junction_consequence(plus, fasta, "fwd", start, end)
        reverse = predict_junction_consequence(minus, fasta, "rev", *_reflect(start, end))
    assert forward is not None, label
    assert reverse is not None, f"the minus strand could not interpret it: {label}"
    _same(forward, reverse, label)


@pytest.mark.parametrize(("exon_start", "exon_end"), [(250, 300), (250, 299), (301, 390)])
def test_a_cassette_exon_is_called_the_same_on_either_strand(mirrored, exon_start, exon_end):
    genome, plus, minus = mirrored
    with genome as fasta:
        forward = predict_consequence(plus, fasta, "fwd", exon_start, exon_end)
        reverse = predict_consequence(minus, fasta, "rev", *_reflect(exon_start, exon_end))
    _same(forward, reverse, f"cassette {exon_start}-{exon_end}")


def test_the_mirror_is_not_passing_because_nothing_is_predicted(mirrored):
    """A mirror test agrees trivially if both sides give up, so at least one of these
    events has to reach a real consequence with real numbers in it."""
    genome, plus, _ = mirrored
    with genome as fasta:
        extension = predict_junction_consequence(plus, fasta, "fwd", 201, 380)
        truncation = predict_junction_consequence(plus, fasta, "fwd", 201, 420)
        cassette = predict_consequence(plus, fasta, "fwd", 250, 299)

    assert extension.insert_length == 20 and extension.frameshift  # 381..400
    assert truncation.insert_length == -20 and truncation.frameshift  # 401..420
    assert cassette.insert_length == 50 and cassette.frameshift
    for call in (extension, truncation, cassette):
        assert call.consequence_class not in {"no_host_transcript", "non_coding_host"}
    assert any(c.ptc_offset is not None for c in (extension, truncation, cassette))


def test_coding_length_before_a_position_counts_downwards_on_the_minus_strand():
    """On the minus strand the bases transcribed *before* a position are the ones at
    higher coordinates, and the base at the position itself is not among them."""
    tx = _transcript("rev", "-", [(101, 200), (401, 500)])
    assert tx.exons == [(401, 500), (101, 200)], "transcription order, not genomic"
    assert tx.cds_length_before(500) == 0
    assert tx.cds_length_before(450) == 50   # 451..500
    assert tx.cds_length_before(401) == 99   # 402..500
    assert tx.cds_length_before(150) == 150  # all of 401..500, then 151..200


def test_the_annotated_phase_is_subtracted_from_the_frame():
    """A 5'-incomplete CDS (GENCODE ``cds_start_NF``) begins part-way through a codon,
    and the GTF phase says how many bases to drop. Adding it instead of subtracting it
    gives a different frame for every phase but 0, and so a different stop codon."""
    for phase, expected in ((0, 1), (1, 0), (2, 2)):
        tx = Transcript(
            transcript_id="T1", gene_id="G1", gene_name="GENE1", chrom="fwd", strand="+",
            exons=[(101, 200)], cds=[(101, 200)], cds_phase=phase,
        )
        assert tx.cds_length_before(150) == 49
        assert tx.frame_at(150) == expected, f"phase {phase}"


def test_consequences_run_end_to_end_on_a_minus_strand_table(mirrored):
    """The pipeline helper, not just the per-transcript function. Nothing in
    `test_consequence_pipeline.py` used the minus strand at all."""
    from splicescope.consequence import annotate_junction_consequences

    genome, _, minus = mirrored
    rows = []
    for start, end, sclass in ((201, 380, "novel_acceptor"), (230, 400, "novel_donor")):
        reflected = _reflect(start, end)
        rows.append(
            {
                "chrom": "rev", "strand": "-", "gene_id": "G1", "sclass": sclass,
                "intron_start": reflected[0], "intron_end": reflected[1],
            }
        )
    with genome as fasta:
        out = annotate_junction_consequences(
            pd.DataFrame(rows), {minus.transcript_id: minus}, fasta
        )
    assert len(out) == 2
    assert out["consequence_class"].notna().all()
    assert (out["insert_length"] > 0).all(), "both of these extend an exon"


@pytest.mark.parametrize(
    ("strand", "resume"),
    [("+", 200), ("-", 401)],
)
def test_the_exon_that_holds_the_resume_point_is_kept(mirrored, strand, resume):
    """Transcription resumes *at* a base, not after it, so the exon containing that
    base contributes it. On the plus strand the boundary case is an exon ending there
    and on the minus strand one starting there; a `<=` in place of the `<` drops the
    exon entirely and shortens the sequence a ribosome is modelled as reading."""
    from splicescope.consequence import downstream_sequence

    genome, _, _ = mirrored
    tx = _transcript("fwd", strand, [(101, 200), (401, 500)])
    with genome as fasta:
        sequence, lengths = downstream_sequence(tx, fasta, "fwd", resume)
    assert lengths == [1, 100], "one base from the boundary exon, then the whole next one"
    assert len(sequence) == 101
