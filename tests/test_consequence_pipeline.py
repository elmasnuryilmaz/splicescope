"""The protein-consequence layer, end to end on the built-in synthetic data.

Until the simulator emitted a genome these paths could only be tested against
hand-built transcripts. Here the same code runs on a generated chromosome that
carries real open reading frames, and is driven through the command line the way
a user drives it.
"""

import pandas as pd
import pytest

from splicescope.cli import main
from splicescope.consequence import STOP_CODONS, GenomeFasta, load_transcripts
from splicescope.simulate import gene_exons, simulate_dataset, simulate_genome, write_dataset

CONSEQUENCE_CLASSES = {
    "ptc_nmd",
    "ptc_escape",
    "frameshift",
    "exon_truncation",
    "in_frame_insertion",
    "utr_insertion",
    "non_coding_host",
    "no_host_transcript",
}


@pytest.fixture(scope="module")
def demo(tmp_path_factory):
    """A simulated dataset written to disk, with its genome, run once per module."""
    data = tmp_path_factory.mktemp("demo")
    ds = simulate_dataset(n_genes=9, n_per_group=4, cryptic_fraction=0.7, seed=5)
    write_dataset(ds, data, seed=5)
    return data


def _codons(seq):
    return [seq[i : i + 3] for i in range(0, len(seq), 3)]


def test_simulated_genome_encodes_real_open_reading_frames(demo):
    """Every simulated gene must translate cleanly, or consequence calls mean nothing."""
    transcripts = load_transcripts(demo / "annotation.gtf")
    assert transcripts
    with GenomeFasta(demo / "genome.fa") as fasta:
        for tx in transcripts.values():
            cds = "".join(fasta.fetch(tx.chrom, s, e, tx.strand) for s, e in tx.cds)
            assert len(cds) % 3 == 0, f"{tx.transcript_id} CDS is not a whole number of codons"
            codons = _codons(cds)
            assert not [c for c in codons[:-1] if c in STOP_CODONS], (
                f"{tx.transcript_id} has an internal stop codon"
            )
            assert codons[-1] in STOP_CODONS, f"{tx.transcript_id} does not end in a stop"


def test_simulated_introns_carry_canonical_splice_motifs(demo):
    transcripts = load_transcripts(demo / "annotation.gtf")
    with GenomeFasta(demo / "genome.fa") as fasta:
        motifs = {
            fasta.fetch(tx.chrom, i, i + 1) + "/" + fasta.fetch(tx.chrom, j - 1, j)
            for tx in transcripts.values()
            for i, j in tx.introns
        }
    assert motifs == {"GT/AG"}


def test_simulated_genes_span_all_three_reading_frames(demo):
    """A fixed UTR would put every intron on a codon boundary and never test frame
    inheritance, which is the subtle half of the consequence prediction."""
    transcripts = load_transcripts(demo / "annotation.gtf")
    frames = {tx.cds_length_before(tx.introns[0][0]) % 3 for tx in transcripts.values()}
    assert frames == {0, 1, 2}


def test_write_dataset_writes_a_genome_its_index_can_address(demo):
    """The .fai offsets must agree with the file, or every fetched base is shifted."""
    raw = (demo / "genome.fa").read_text().splitlines()
    assert raw[0].startswith(">")
    sequence = "".join(raw[1:])

    name, length, offset, linebases, linewidth = (
        (demo / "genome.fa.fai").read_text().split()
    )
    assert name == "chr1"
    assert int(length) == len(sequence)
    assert int(offset) == len(raw[0]) + 1
    assert int(linewidth) == int(linebases) + 1

    with GenomeFasta(demo / "genome.fa") as fasta:
        # spot-check positions that straddle FASTA line boundaries
        for start in (1, 59, 60, 61, 119, len(sequence) - 5):
            assert fasta.fetch("chr1", start, start + 4) == sequence[start - 1 : start + 4]


def test_simulate_genome_leaves_the_junction_data_untouched():
    """The genome uses its own random stream; published junction numbers must not move."""
    ds = simulate_dataset(n_genes=6, seed=3)
    before = ds.observed.copy()
    simulate_genome(ds, seed=3)
    pd.testing.assert_frame_equal(before, ds.observed)


def test_gene_exons_reconstructs_blocks_that_flank_every_known_intron():
    ds = simulate_dataset(n_genes=4, seed=1)
    exons = gene_exons(ds)
    for gene_id, sub in ds.known.groupby("gene_id"):
        blocks = exons[str(gene_id)]
        assert len(blocks) == len(sub) + 1
        gaps = {(a[1] + 1, b[0] - 1) for a, b in zip(blocks, blocks[1:], strict=False)}
        assert gaps == set(zip(sub["start"], sub["end"], strict=False))


def _run(demo, outdir, *extra):
    return main(
        [
            "run",
            "--sj-dir", str(demo / "sj"),
            "--gtf", str(demo / "annotation.gtf"),
            "--groups", str(demo / "groups.tsv"),
            "--outdir", str(outdir),
            *extra,
        ]
    )


@pytest.fixture(scope="module")
def results(demo, tmp_path_factory):
    out = tmp_path_factory.mktemp("results")
    assert _run(demo, out, "--genome", str(demo / "genome.fa")) == 0
    return out


def test_run_with_a_genome_predicts_consequences_for_cassettes_and_shifts(results):
    cassettes = pd.read_csv(results / "consequence.tsv", sep="\t")
    shifts = pd.read_csv(results / "junction_consequence.tsv", sep="\t")
    assert not cassettes.empty and not shifts.empty
    for table in (cassettes, shifts):
        assert set(table["consequence_class"]) <= CONSEQUENCE_CLASSES
    # a real prediction was made, not a table of "could not interpret"
    assert (shifts["consequence_class"] != "no_host_transcript").any()
    assert (results / "figures" / "consequence.png").exists()
    assert (results / "figures" / "junction_consequence.png").exists()
    assert (results / "figures" / "event_summary.png").exists()


def test_run_only_reports_shifts_for_anchored_novel_junctions(results):
    shifts = pd.read_csv(results / "junction_consequence.tsv", sep="\t")
    assert set(shifts["sclass"]) <= {"novel_donor", "novel_acceptor"}


def test_cassette_junctions_are_not_also_reported_as_splice_site_shifts(results):
    """Read alone, a cassette inclusion junction looks like an exon extension running
    to the end of the intron. It is already explained as an event, so it must not be
    re-interpreted — otherwise the same signal is reported twice, wrongly the second time."""
    events = pd.read_csv(results / "events.tsv", sep="\t")
    se = events[events["event_type"] == "SE"]
    used = set()
    for row in se.itertuples(index=False):
        for s, e in (
            (row.inc1_start, row.inc1_end),
            (row.inc2_start, row.inc2_end),
            (row.skip_start, row.skip_end),
        ):
            used.add((row.chrom, int(s), int(e), row.strand))
    assert used, "the demo must contain cassette events for this test to mean anything"

    shifts = pd.read_csv(results / "junction_consequence.tsv", sep="\t")
    reported = {
        (r.chrom, int(r.start), int(r.end), r.strand) for r in shifts.itertuples(index=False)
    }
    assert reported & used == set()


def test_run_without_a_genome_skips_consequence(demo, tmp_path):
    out = tmp_path / "nogenome"
    assert _run(demo, out) == 0
    assert (out / "events.tsv").exists()
    assert not (out / "consequence.tsv").exists()
    assert not (out / "junction_consequence.tsv").exists()


def test_run_reports_a_missing_genome_index_instead_of_crashing(demo, tmp_path):
    bare = tmp_path / "unindexed.fa"
    bare.write_text(">chr1\nACGT\n")  # no .fai beside it
    assert _run(demo, tmp_path / "out", "--genome", str(bare)) == 2


def test_consequence_command_auto_detects_junction_input(demo, tmp_path):
    """A plain junction table has start/end, not exon_start/exon_end."""
    transcripts = load_transcripts(demo / "annotation.gtf")
    tx = next(iter(transcripts.values()))
    intron_start, intron_end = tx.introns[0]
    events = pd.DataFrame(
        [
            {
                "chrom": tx.chrom,
                "start": intron_start,
                "end": intron_end - 30,  # annotated donor, novel acceptor 30 nt early
                "strand": tx.strand,
                "gene_id": tx.gene_id,
            }
        ]
    )
    src = tmp_path / "junctions.tsv"
    events.to_csv(src, sep="\t", index=False)
    out = tmp_path / "jc.tsv"

    rc = main(
        [
            "consequence",
            "--events", str(src),
            "--gtf", str(demo / "annotation.gtf"),
            "--genome", str(demo / "genome.fa"),
            "--out", str(out),
        ]
    )
    assert rc == 0
    table = pd.read_csv(out, sep="\t")
    assert len(table) == 1
    # the junction is interpretable against the annotation, so it is not "no host"
    assert table.loc[0, "consequence_class"] in CONSEQUENCE_CLASSES - {"no_host_transcript"}


def test_consequence_command_auto_detects_cassette_input(demo, tmp_path):
    transcripts = load_transcripts(demo / "annotation.gtf")
    tx = next(iter(transcripts.values()))
    intron_start, _ = tx.introns[0]
    events = pd.DataFrame(
        [
            {
                "chrom": tx.chrom,
                "strand": tx.strand,
                "exon_start": intron_start + 100,
                "exon_end": intron_start + 160,
                "gene_id": tx.gene_id,
            }
        ]
    )
    src = tmp_path / "exons.tsv"
    events.to_csv(src, sep="\t", index=False)
    out = tmp_path / "ec.tsv"

    rc = main(
        [
            "consequence",
            "--events", str(src),
            "--gtf", str(demo / "annotation.gtf"),
            "--genome", str(demo / "genome.fa"),
            "--out", str(out),
        ]
    )
    assert rc == 0
    table = pd.read_csv(out, sep="\t")
    assert table.loc[0, "insert_length"] == 61
    assert table.loc[0, "consequence_class"] in CONSEQUENCE_CLASSES


def test_consequence_command_rejects_a_table_it_cannot_interpret(demo, tmp_path):
    src = tmp_path / "nonsense.tsv"
    pd.DataFrame([{"chrom": "chr1", "strand": "+", "foo": 1}]).to_csv(
        src, sep="\t", index=False
    )
    rc = main(
        [
            "consequence",
            "--events", str(src),
            "--gtf", str(demo / "annotation.gtf"),
            "--genome", str(demo / "genome.fa"),
            "--out", str(tmp_path / "out.tsv"),
        ]
    )
    assert rc == 2
