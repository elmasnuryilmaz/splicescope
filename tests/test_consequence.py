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


# --------------------------------------------------------------------------------
# a stop where the protein natively ends is not a premature stop
# --------------------------------------------------------------------------------


def test_an_in_frame_exon_truncation_is_not_a_premature_stop(tmp_path):
    """Removing whole codons leaves the downstream frame untouched, so the first stop
    found is the transcript's own. Reporting it as a PTC turned an ordinary in-frame
    deletion into a predicted termination event."""
    from splicescope.consequence import EXON_TRUNCATION, predict_junction_consequence

    # exons 101-200, 401-500, 701-800; the native stop is the last codon of the CDS
    tx = make_transcript()
    chrom = ["A"] * 900
    for i in range(100, 200):
        chrom[i] = "AAC"[(i - 100) % 3]
    for i in range(400, 500):
        chrom[i] = "AAC"[(i - 400) % 3]
    for i in range(700, 800):
        chrom[i] = "AAC"[(i - 700) % 3]
    # place the native stop as the final codon of the last CDS block
    mature_len = 300
    assert mature_len % 3 == 0
    chrom[797], chrom[798], chrom[799] = "T", "A", "A"
    fasta = write_fasta(tmp_path, {"chr1": "".join(chrom)})

    with GenomeFasta(fasta) as fa:
        # intron 1 is 201-400; a donor shifted 3 nt earlier removes one codon
        in_frame = predict_junction_consequence(tx, fa, "chr1", 198, 400)
        assert in_frame is not None
        assert in_frame.consequence_class == EXON_TRUNCATION
        assert in_frame.ptc_offset is None
        assert in_frame.frameshift is False


def test_a_frame_shifting_truncation_still_reports_its_premature_stop(tmp_path):
    from splicescope.consequence import predict_junction_consequence

    tx = make_transcript()
    chrom = ["A"] * 900
    for start in (100, 400, 700):
        for i in range(start, start + 100):
            chrom[i] = "AAC"[(i - start) % 3]
    # Removing 2 nt leaves frame 2, so the scan starts at offset 1 of the retained
    # sequence and steps by 3: offset 19 is read as a codon, offset 20 is not.
    # Exon 2 starts at genomic 401, so offset 19 is chrom index 419.
    chrom[419], chrom[420], chrom[421] = "T", "G", "A"
    fasta = write_fasta(tmp_path, {"chr1": "".join(chrom)})

    with GenomeFasta(fasta) as fa:
        shifted = predict_junction_consequence(tx, fa, "chr1", 199, 400)  # removes 2 nt
        assert shifted.frameshift is True
        assert shifted.consequence_class in (PTC_NMD, PTC_ESCAPE)
        assert shifted.ptc_offset == 19


# --------------------------------------------------------------------------------
# GTF phase: 5'-incomplete coding sequences
# --------------------------------------------------------------------------------


def test_cds_phase_is_read_from_the_gtf_and_shifts_the_frame(tmp_path):
    """GENCODE marks 5'-incomplete transcripts (cds_start_NF) with a non-zero phase on
    the first CDS record. Ignoring it translates the whole transcript out of frame."""
    from splicescope.consequence import load_transcripts

    def gtf(phase):
        attrs = 'gene_id "G1"; transcript_id "T1"; gene_name "GENE1";'
        return (
            f"chr1\tsrc\texon\t101\t200\t.\t+\t.\t{attrs}\n"
            f"chr1\tsrc\tCDS\t101\t200\t.\t+\t{phase}\t{attrs}\n"
            f"chr1\tsrc\texon\t401\t500\t.\t+\t.\t{attrs}\n"
            f"chr1\tsrc\tCDS\t401\t500\t.\t+\t2\t{attrs}\n"
        )

    frames = {}
    for phase in (0, 1, 2):
        path = tmp_path / f"p{phase}.gtf"
        path.write_text(gtf(phase))
        tx = load_transcripts(path)["T1"]
        assert tx.cds_phase == phase
        frames[phase] = tx.frame_at(401)

    # the three phases must give three different reading frames at the same position
    assert len(set(frames.values())) == 3


def test_cds_phase_defaults_to_zero_when_the_column_is_a_dot(tmp_path):
    from splicescope.consequence import load_transcripts

    attrs = 'gene_id "G1"; transcript_id "T1";'
    path = tmp_path / "dot.gtf"
    path.write_text(
        f"chr1\tsrc\texon\t101\t200\t.\t+\t.\t{attrs}\n"
        f"chr1\tsrc\tCDS\t101\t200\t.\t+\t.\t{attrs}\n"
    )
    assert load_transcripts(path)["T1"].cds_phase == 0


def test_cds_phase_is_taken_from_the_first_block_in_transcription_order(tmp_path):
    """On the minus strand the first coding block is the genomically last one."""
    from splicescope.consequence import load_transcripts

    attrs = 'gene_id "G1"; transcript_id "T1";'
    path = tmp_path / "minus.gtf"
    path.write_text(
        f"chr1\tsrc\texon\t101\t200\t.\t-\t.\t{attrs}\n"
        f"chr1\tsrc\tCDS\t101\t200\t.\t-\t1\t{attrs}\n"
        f"chr1\tsrc\texon\t401\t500\t.\t-\t.\t{attrs}\n"
        f"chr1\tsrc\tCDS\t401\t500\t.\t-\t2\t{attrs}\n"
    )
    assert load_transcripts(path)["T1"].cds_phase == 2


def test_a_ptc_in_an_extension_of_the_final_exon_escapes_nmd(tmp_path):
    """An extension is contiguous with the exon it joins — there is no junction between
    them. When that exon is the last one, the PTC lies past the final exon-exon junction
    and NMD cannot be triggered, however far it sits from the transcript's end."""
    from splicescope.consequence import predict_junction_consequence

    tx = make_transcript()  # exons 101-200, 401-500, 701-800; last intron 501-700
    seq = list("A" * 900)
    # An acceptor shifted 150 nt earlier extends the final exon to 551-700. The frame
    # there is 2, so codons are read from offset 1: put the stop at offset 4 (genomic 555).
    for i in range(550, 700):
        seq[i] = "GGC"[(i - 550) % 3]
    seq[554], seq[555], seq[556] = "T", "A", "A"
    fasta = write_fasta(tmp_path, {"chr1": "".join(seq)})

    with GenomeFasta(fasta) as fa:
        final = predict_junction_consequence(tx, fa, "chr1", 501, 550)
        assert final is not None
        assert final.consequence_class == PTC_ESCAPE
        assert final.distance_to_last_junction is None
        assert final.nmd_predicted is False


def test_a_ptc_in_an_extension_of_an_internal_exon_can_still_trigger_nmd(tmp_path):
    """The other half: an internal exon does have a junction after it."""
    from splicescope.consequence import predict_junction_consequence

    tx = make_transcript()  # first intron 201-400 precedes exon 401-500
    seq = list("A" * 900)
    # extension 251-400; frame there is 1, so codons start at offset 2 — stop at offset 5
    for i in range(250, 400):
        seq[i] = "GGC"[(i - 250) % 3]
    seq[255], seq[256], seq[257] = "T", "A", "A"
    fasta = write_fasta(tmp_path, {"chr1": "".join(seq)})

    with GenomeFasta(fasta) as fa:
        internal = predict_junction_consequence(tx, fa, "chr1", 201, 250)
        assert internal.consequence_class == PTC_NMD
        assert internal.distance_to_last_junction > 50


def test_an_exon_skipping_junction_is_not_read_as_a_splice_site_shift():
    """Both of its sites are annotated, just not as a pair. Interpreting it as a shift
    reports the skipped exon *and the intron beyond it* as one contiguous deletion."""
    from splicescope.consequence import junction_change

    tx = Transcript(
        transcript_id="T1", gene_id="G1", gene_name="G", chrom="chr1", strand="+",
        exons=[(101, 200), (401, 500), (701, 800), (1001, 1100)],
        cds=[(101, 200), (401, 500), (701, 800), (1001, 1100)],
    )
    assert tx.introns == [(201, 400), (501, 700), (801, 1000)]

    # skips exon 401-500: donor of intron 1, acceptor of intron 2
    assert junction_change(tx, 201, 700) is None
    # the annotated intron itself is not a change either
    assert junction_change(tx, 201, 400) is None
    # genuine shifts are still interpreted
    assert junction_change(tx, 201, 350) == ("extension", 351, 400)
    assert junction_change(tx, 251, 400) == ("extension", 201, 250)


def test_the_fifty_nucleotide_boundary_is_exclusive_and_measured_from_the_stop_codon():
    """The rule this tool is built on: a premature stop triggers decay when it sits
    more than 50 nucleotides upstream of the last exon-exon junction. Three things
    about that sentence are separately wrong-able — which junction ("last", so every
    exon length but the final one), which end of the stop codon the distance is
    measured from (its end, so ptc + 3), and whether 50 itself counts (it does not).

    One sequence is used twice here, with the downstream exons split one base
    differently. Nothing about the stop changes; only the distance to the junction
    does, from exactly 50 to exactly 51, and the call has to flip there and nowhere
    else.
    """
    from splicescope.consequence import NMD_DISTANCE_RULE, _nmd_from_downstream

    assert NMD_DISTANCE_RULE == 50
    sequence = "AAA" * 20 + "TAA" + "AAA" * 47  # first in-frame stop at offset 60
    assert len(sequence) == 204

    offset, distance, nmd = _nmd_from_downstream(sequence, [113, 91], frame=0, offset_before=0)
    assert (offset, distance) == (60, 50)
    assert not nmd, "a stop exactly 50 nt from the last junction escapes decay"

    offset, distance, nmd = _nmd_from_downstream(sequence, [114, 90], frame=0, offset_before=0)
    assert (offset, distance) == (60, 51)
    assert nmd, "one nucleotide further and the same stop triggers decay"


def test_a_stop_in_the_last_exon_has_no_junction_to_be_upstream_of():
    """The last-exon exception, which is not a special case bolted on but the rule
    read literally: decay needs a junction downstream of the stop, and a stop in the
    final exon has none. The distance is then not a large number or a negative one, it
    does not exist — and a caller that prints it must show that."""
    from splicescope.consequence import _nmd_from_downstream

    sequence = "AAA" * 20 + "TAA" + "AAA" * 47
    offset, distance, nmd = _nmd_from_downstream(sequence, [204], frame=0, offset_before=0)
    assert offset == 60
    assert distance is None, "there is no last junction, so there is no distance"
    assert not nmd


def test_the_reported_stop_offset_is_on_the_scale_of_the_event():
    """``offset_before`` is how many transcript bases precede the searched sequence, so
    the offset a caller sees stays comparable with the event's own coordinates."""
    from splicescope.consequence import _nmd_from_downstream

    sequence = "AAA" * 20 + "TAA" + "AAA" * 47
    plain = _nmd_from_downstream(sequence, [114, 90], frame=0, offset_before=0)
    shifted = _nmd_from_downstream(sequence, [114, 90], frame=0, offset_before=300)
    assert shifted[0] == plain[0] + 300
    assert shifted[1:] == plain[1:], "only the offset moves; the distance and call do not"


def test_every_consequence_class_is_explained_in_a_sentence():
    """The table says `ptc_nmd`, `insert_length=61`, `distance_to_last_junction=395`.
    That is an answer only to a reader who already knows the rule, and the dashboard
    shows it to people who do not."""
    from splicescope.consequence import CONSEQUENCE_CLASSES, describe

    base = dict(
        transcript_id="T1", gene_name="G", insert_length=61, frame_offset=0,
        frameshift=True, ptc_offset=23, distance_to_last_junction=395,
        nmd_predicted=True,
    )
    for name in CONSEQUENCE_CLASSES:
        sentence = describe({**base, "consequence_class": name})
        assert sentence and sentence[0].isupper() and sentence.endswith("."), name
        assert "consequence_class" not in sentence, f"{name}: it printed the column name"
        assert len(sentence.split()) > 8, f"{name}: too terse to explain anything"


def test_the_fifty_nucleotide_rule_is_stated_in_words_both_ways():
    from splicescope.consequence import describe

    base = dict(insert_length=61, frameshift=True, ptc_offset=23, nmd_predicted=True)
    decayed = describe(
        {**base, "consequence_class": "ptc_nmd", "distance_to_last_junction": 395}
    )
    assert "395" in decayed and "50" in decayed
    assert "more than" in decayed and "decay" in decayed

    escaped = describe(
        {**base, "consequence_class": "ptc_escape", "distance_to_last_junction": 20}
    )
    assert "20" in escaped and "within the 50" in escaped
    assert "truncated protein" in escaped


def test_a_stop_in_the_final_exon_is_explained_by_the_missing_junction():
    """Not "the distance was too small" — there is no junction downstream at all, which
    is the last-exon exception read literally."""
    from splicescope.consequence import describe

    sentence = describe(
        {
            "consequence_class": "ptc_escape", "insert_length": 61, "frameshift": True,
            "ptc_offset": 104, "distance_to_last_junction": None, "nmd_predicted": False,
        }
    )
    assert "final exon" in sentence
    assert "no exon-exon junction downstream" in sentence
    assert "truncated protein" in sentence


def test_a_removed_stretch_is_described_as_removed():
    """A splice-site shift that truncates an exon carries a negative insert length."""
    from splicescope.consequence import describe

    sentence = describe(
        {
            "consequence_class": "exon_truncation", "insert_length": -20,
            "frameshift": True, "ptc_offset": None,
            "distance_to_last_junction": None, "nmd_predicted": False,
        }
    )
    assert sentence.startswith("Removing these 20 nucleotides")


def test_an_in_frame_change_says_whether_the_protein_gains_or_loses_residues():
    """A cassette exon adds codons; a shifted splice site that trims one removes them.
    The sentence has to say which, and swapping the two words was invisible to the
    suite until this test."""
    from splicescope.consequence import describe

    added = describe(
        {
            "consequence_class": "in_frame_insertion", "insert_length": 42,
            "frameshift": False, "ptc_offset": None,
            "distance_to_last_junction": None, "nmd_predicted": False,
        }
    )
    assert "Including these 42 nucleotides" in added
    assert "gains 14 amino acids" in added and "loses" not in added

    removed = describe(
        {
            "consequence_class": "in_frame_insertion", "insert_length": -42,
            "frameshift": False, "ptc_offset": None,
            "distance_to_last_junction": None, "nmd_predicted": False,
        }
    )
    assert "Removing these 42 nucleotides" in removed
    assert "loses 14 amino acids" in removed and "gains" not in removed


def _tiny_genome(tmp_path, name: str, length: int):
    """A FASTA whose index says exactly `length`, so the check has something to read."""
    fasta = tmp_path / f"{name}.fa"
    fasta.write_text(f">{name}\n" + "ACGT" * (length // 4) + "\n")
    index = f"{name}\t{length}\t{len(name) + 2}\t{length}\t{length + 1}\n"
    fasta.with_suffix(".fa.fai").write_text(index)
    return fasta


def test_a_genome_from_another_assembly_is_refused_rather_than_read_off_the_end(tmp_path):
    """Reading past the end of a contig returns a truncated string, not an error, so a
    mismatched genome does not fail — it finds no stop codon anywhere and calls
    everything a frameshift. On the simulator that turned six `ptc_nmd` calls into zero
    while reporting the same nine events as confidently as before."""
    events = pd.DataFrame([("chr1", 250, 300, "+")], columns=["chrom", "start", "end", "strand"])
    with GenomeFasta(_tiny_genome(tmp_path, "chr1", 200)) as fa:
        with pytest.raises(ValueError, match="only 200 bases long"):
            annotate_consequences(events, {"T1": make_transcript()}, fa)


def test_a_genome_named_the_other_way_round_says_which_convention(tmp_path):
    """GENCODE writes `chr1`, Ensembl writes `1`, and the genome has to match the GTF."""
    events = pd.DataFrame([("chr1", 250, 300, "+")], columns=["chrom", "start", "end", "strand"])
    with GenomeFasta(_tiny_genome(tmp_path, "1", 4000)) as fa:
        with pytest.raises(ValueError, match="has no contig"):
            annotate_consequences(events, {"T1": make_transcript()}, fa)
    with GenomeFasta(_tiny_genome(tmp_path, "1", 4000)) as fa:
        try:
            annotate_consequences(events, {"T1": make_transcript()}, fa)
        except ValueError as exc:
            assert "GENCODE writes 'chr1' where Ensembl writes '1'" in str(exc)


def test_a_contig_the_annotation_does_not_use_is_not_the_genome_s_problem(tmp_path):
    """The check covers only the chromosomes where a prediction will actually be made —
    those with both an event and a transcript to host it. An event on a contig the
    annotation says nothing about already has an honest answer, and a genome that covers
    more than the analysis needs is nobody's mistake."""
    events = pd.DataFrame(
        [("chr1", 250, 300, "+"), ("chr9", 250, 300, "+")],
        columns=["chrom", "start", "end", "strand"],
    )
    with GenomeFasta(write_fasta(tmp_path, {"chr1": "A" * 1000})) as fa:
        out = annotate_consequences(events, {"T1": make_transcript()}, fa)
    assert len(out) == 2
    assert out.loc[1, "consequence_class"] == NO_HOST
