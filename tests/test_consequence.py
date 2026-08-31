"""Consequence prediction: frame inheritance, PTC detection and the NMD rule."""

from __future__ import annotations

import pandas as pd
import pytest

from splicescope.consequence import (
    FRAMESHIFT,
    IN_FRAME,
    NO_HOST,
    NON_CODING_HOST,
    PTC_ESCAPE,
    PTC_NMD,
    UTR_INSERTION,
    GenomeFasta,
    Transcript,
    annotate_consequences,
    find_ptc,
    host_transcripts,
    predict_consequence,
    reverse_complement,
)

LINE = 60


def write_fasta(tmp_path, sequences: dict[str, str]):
    """Write a FASTA plus the ``.fai`` index GenomeFasta needs."""
    fasta = tmp_path / "genome.fa"
    fai_rows, offset = [], 0
    with open(fasta, "w") as handle:
        for name, seq in sequences.items():
            header = f">{name}\n"
            handle.write(header)
            offset += len(header)
            fai_rows.append(f"{name}\t{len(seq)}\t{offset}\t{LINE}\t{LINE + 1}")
            for i in range(0, len(seq), LINE):
                chunk = seq[i : i + LINE] + "\n"
                handle.write(chunk)
                offset += len(chunk)
    fasta.with_suffix(".fa.fai").write_text("\n".join(fai_rows) + "\n")
    return fasta


def make_transcript(strand="+", cds=True, exons=None):
    exons = exons or [(101, 200), (401, 500), (701, 800)]
    return Transcript(
        transcript_id="T1",
        gene_id="G1",
        gene_name="GENE1",
        chrom="chr1",
        strand=strand,
        exons=list(exons),
        cds=list(exons) if cds else [],
    )


def test_reverse_complement_round_trips():
    assert reverse_complement("ATGCN") == "NGCAT"
    assert reverse_complement(reverse_complement("ACGTACGT")) == "ACGTACGT"


def test_find_ptc_respects_inherited_frame():
    # TAA sits at offset 0; with one base of the codon already upstream the
    # reading frame starts at offset 2 and the stop is no longer in frame.
    assert find_ptc("TAAGGGCCC", 0) == 0
    assert find_ptc("TAAGGGCCC", 1) is None


def test_find_ptc_returns_none_without_stop():
    assert find_ptc("GGGCCCGGGCCC", 0) is None


def test_introns_are_in_transcription_order():
    assert make_transcript("+").introns == [(201, 400), (501, 700)]
    assert make_transcript("-").introns == [(501, 700), (201, 400)]


def test_host_transcripts_requires_containment_and_strand():
    tx = {"T1": make_transcript("+")}
    assert len(host_transcripts(tx, "chr1", 250, 300, "+")) == 1
    assert host_transcripts(tx, "chr1", 250, 300, "-") == []  # wrong strand
    assert host_transcripts(tx, "chr1", 150, 300, "+") == []  # overlaps an exon
    assert host_transcripts(tx, "chr1", 250, 300, "chr2") == []


def test_ptc_far_from_last_junction_is_called_nmd(tmp_path):
    # A 99 nt first exon leaves the frame at 0, so the stop codon three bases
    # into the insert is read in frame, with a whole exon downstream of it.
    exons = [(101, 199), (401, 500), (701, 800)]
    insert = "TTT" + "TAA" + "GGG" * 40
    seq = "A" * 200 + insert + "A" * (1000 - 200 - len(insert))
    with GenomeFasta(write_fasta(tmp_path, {"chr1": seq})) as fa:
        result = predict_consequence(
            make_transcript("+", exons=exons), fa, "chr1", 201, 200 + len(insert)
        )
    assert result.frame_offset == 0
    assert result.consequence_class == PTC_NMD
    assert result.ptc_offset == 3
    assert result.nmd_predicted


def test_ptc_close_to_last_junction_escapes_nmd(tmp_path):
    # Single downstream exon and a stop at the very end of the insert, so the
    # 50-nt rule is not met and the transcript should escape decay.
    exons = [(101, 199), (401, 500)]
    insert = "GGG" * 20 + "TGA"
    seq = "A" * 200 + insert + "A" * (1000 - 200 - len(insert))
    with GenomeFasta(write_fasta(tmp_path, {"chr1": seq})) as fa:
        result = predict_consequence(
            make_transcript("+", exons=exons), fa, "chr1", 201, 200 + len(insert)
        )
    assert result.consequence_class == PTC_ESCAPE
    assert not result.nmd_predicted


def test_stopless_insert_is_frameshift_or_in_frame(tmp_path):
    seq = "A" * 200 + "GGC" * 30 + "A" * 500
    with GenomeFasta(write_fasta(tmp_path, {"chr1": seq})) as fa:
        tx = make_transcript("+")
        in_frame = predict_consequence(tx, fa, "chr1", 201, 290)  # 90 nt
        shifted = predict_consequence(tx, fa, "chr1", 201, 291)  # 91 nt
    assert in_frame.consequence_class == IN_FRAME
    assert shifted.consequence_class == FRAMESHIFT
    assert shifted.frameshift


def test_non_coding_host_and_utr_are_separated(tmp_path):
    seq = "A" * 1000
    with GenomeFasta(write_fasta(tmp_path, {"chr1": seq})) as fa:
        no_cds = predict_consequence(make_transcript(cds=False), fa, "chr1", 250, 300)
        # CDS restricted to the first exon, so a later intron is outside it.
        utr = make_transcript()
        utr.cds = [(101, 200)]
        outside = predict_consequence(utr, fa, "chr1", 550, 600)
    assert no_cds.consequence_class == NON_CODING_HOST
    assert outside.consequence_class == UTR_INSERTION


def test_annotate_consequences_reports_no_host(tmp_path):
    seq = "A" * 1000
    events = pd.DataFrame(
        [("chr1", 250, 300, "+"), ("chr9", 250, 300, "+")],
        columns=["chrom", "start", "end", "strand"],
    )
    with GenomeFasta(write_fasta(tmp_path, {"chr1": seq})) as fa:
        out = annotate_consequences(events, {"T1": make_transcript()}, fa)
    assert list(out["consequence_class"])[1] == NO_HOST
    assert out.loc[0, "n_host_transcripts"] == 1
    assert len(out) == 2


def test_missing_fasta_index_is_reported(tmp_path):
    fasta = tmp_path / "nope.fa"
    fasta.write_text(">chr1\nACGT\n")
    with pytest.raises(FileNotFoundError):
        GenomeFasta(fasta)


def test_index_by_gene_restricts_the_transcript_search(tmp_path):
    from splicescope.consequence import index_by_gene

    tx = make_transcript()
    other = Transcript(
        transcript_id="T2",
        gene_id="G2.3",
        gene_name="GENE2",
        chrom="chr1",
        strand="+",
        exons=[(101, 200), (401, 500)],
        cds=[(101, 200), (401, 500)],
    )
    index = index_by_gene({"T1": tx, "T2": other})
    # Versioned id, bare id and gene name all reach the same transcript.
    assert index["G2.3"] == [other]
    assert index["G2"] == [other]
    assert index["GENE2"] == [other]
    # Searching only one gene's transcripts gives the same hosts as searching all.
    seq = "A" * 1000
    events = pd.DataFrame(
        [("chr1", 250, 300, "+", "G1")], columns=["chrom", "start", "end", "strand", "gene_id"]
    )
    with GenomeFasta(write_fasta(tmp_path, {"chr1": seq})) as fa:
        scoped = annotate_consequences(events, {"T1": tx, "T2": other}, fa)
    assert scoped.loc[0, "n_host_transcripts"] == 1


# --- splice-site shifts (novel donor / novel acceptor) -----------------------


def test_junction_change_detects_extension_at_either_end():
    from splicescope.consequence import junction_change

    tx = make_transcript("+")  # exons 101-200, 401-500, 701-800; intron 201-400
    # 3' end of the intron moved earlier: the downstream exon gains 351-400.
    assert junction_change(tx, 201, 350) == ("extension", 351, 400)
    # 5' end moved later: the upstream exon gains 201-250.
    assert junction_change(tx, 251, 400) == ("extension", 201, 250)


def test_junction_change_detects_truncation():
    from splicescope.consequence import junction_change

    tx = make_transcript("+")
    # Intron reaches further than annotated, so the exon loses sequence.
    assert junction_change(tx, 201, 450) == ("truncation", 401, 450)
    assert junction_change(tx, 151, 400) == ("truncation", 151, 200)


def test_junction_change_needs_one_annotated_site():
    from splicescope.consequence import junction_change

    tx = make_transcript("+")
    assert junction_change(tx, 250, 350) is None  # neither end is annotated
    assert junction_change(tx, 201, 400) is None  # the annotated intron itself


def test_junction_change_works_on_the_minus_strand():
    from splicescope.consequence import junction_change

    tx = make_transcript("-")
    assert junction_change(tx, 201, 350) == ("extension", 351, 400)
    assert junction_change(tx, 201, 450) == ("truncation", 401, 450)


def test_extension_is_scored_like_a_cassette_exon(tmp_path):
    from splicescope.consequence import predict_junction_consequence

    # A 99 nt first exon leaves the frame at 0, so the stop placed at the first
    # base of the extension is read in frame.
    exons = [(101, 199), (401, 500), (701, 800)]  # intron 200-400
    seq = "A" * 199 + "TAA" + "GGG" * 19 + "A" * (1000 - 259)
    with GenomeFasta(write_fasta(tmp_path, {"chr1": seq})) as fa:
        # The intron's 5' end moved to 261, so the exon gains 200-260.
        result = predict_junction_consequence(
            make_transcript("+", exons=exons), fa, "chr1", 261, 400
        )
    assert result.insert_length == 61
    assert result.ptc_offset == 0
    assert result.consequence_class == PTC_NMD


def test_truncation_reports_frame_only(tmp_path):
    from splicescope.consequence import EXON_TRUNCATION, predict_junction_consequence

    seq = "A" * 1000
    with GenomeFasta(write_fasta(tmp_path, {"chr1": seq})) as fa:
        tx = make_transcript("+")
        shifted = predict_junction_consequence(tx, fa, "chr1", 201, 450)  # loses 50 nt
        in_frame = predict_junction_consequence(tx, fa, "chr1", 201, 451)  # loses 51 nt
    assert shifted.consequence_class == EXON_TRUNCATION
    assert shifted.insert_length == -50
    assert shifted.frameshift
    assert not in_frame.frameshift


def test_annotate_junction_consequences_end_to_end(tmp_path):
    from splicescope.consequence import annotate_junction_consequences

    seq = "A" * 200 + "GGC" * 30 + "A" * 500
    events = pd.DataFrame(
        [("chr1", 261, 400, "+", "G1"), ("chr1", 250, 350, "+", "G1")],
        columns=["chrom", "intron_start", "intron_end", "strand", "gene_id"],
    )
    with GenomeFasta(write_fasta(tmp_path, {"chr1": seq})) as fa:
        out = annotate_junction_consequences(events, {"T1": make_transcript()}, fa)
    assert out.loc[0, "n_host_transcripts"] == 1
    assert out.loc[1, "consequence_class"] == NO_HOST  # unanchored junction


# --- premature stops that lie downstream of the event ------------------------


def test_frameshift_insert_finds_the_stop_in_a_downstream_exon(tmp_path):
    # The insert shifts the frame without carrying a stop itself; the first stop
    # is in the next exon, and must still be found.
    exons = [(101, 199), (401, 500), (701, 800)]
    # A 40 nt insert leaves the downstream frame at 1, so the first codon of the
    # next exon starts two bases in; the stop is placed there.
    seq = "A" * 200 + "G" * 40 + "A" * 160 + "GG" + "TAA" + "G" * 595
    with GenomeFasta(write_fasta(tmp_path, {"chr1": seq})) as fa:
        result = predict_consequence(
            make_transcript("+", exons=exons), fa, "chr1", 201, 240
        )
    assert result.frameshift
    assert result.consequence_class in (PTC_NMD, PTC_ESCAPE)
    assert result.ptc_offset is not None


def test_truncation_finds_a_downstream_stop(tmp_path):
    from splicescope.consequence import predict_junction_consequence

    exons = [(101, 199), (401, 500), (701, 800)]  # intron 200-400
    # Losing 10 bases keeps the frame at 0, and transcription resumes at 411,
    # where the stop sits — after the removed stretch, not inside it.
    seq = "A" * 410 + "TAA" + "G" * 587
    with GenomeFasta(write_fasta(tmp_path, {"chr1": seq})) as fa:
        # The intron reaches to 410, so the second exon loses its first 10 bases.
        result = predict_junction_consequence(
            make_transcript("+", exons=exons), fa, "chr1", 200, 410
        )
    assert result.insert_length == -10
    assert result.frameshift
    assert result.ptc_offset is not None


def test_in_frame_truncation_without_a_stop_stays_a_truncation(tmp_path):
    from splicescope.consequence import EXON_TRUNCATION, predict_junction_consequence

    exons = [(101, 199), (401, 500), (701, 800)]
    seq = "A" * 1000  # poly-A downstream: no stop codon anywhere
    with GenomeFasta(write_fasta(tmp_path, {"chr1": seq})) as fa:
        result = predict_junction_consequence(
            make_transcript("+", exons=exons), fa, "chr1", 200, 409
        )
    assert result.consequence_class == EXON_TRUNCATION
    assert result.ptc_offset is None


def test_downstream_sequence_is_clipped_and_strand_aware(tmp_path):
    from splicescope.consequence import downstream_sequence

    seq = "A" * 100 + "C" * 100 + "G" * 100 + "T" * 700
    with GenomeFasta(write_fasta(tmp_path, {"chr1": seq})) as fa:
        tx = make_transcript("+")  # exons 101-200, 401-500, 701-800
        tail, lengths = downstream_sequence(tx, fa, "chr1", 401)
        assert lengths == [100, 100]
        assert len(tail) == 200
        minus = make_transcript("-")
        tail_m, lengths_m = downstream_sequence(minus, fa, "chr1", 500)
        assert lengths_m == [100, 100]  # exons 401-500 and 101-200, in that order
