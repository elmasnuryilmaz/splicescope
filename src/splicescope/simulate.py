"""Synthetic ground-truth data so the whole toolkit runs anywhere — no downloads.

The generator builds a tiny annotated "genome" and then emits per-sample splice
junctions that mimic what STAR reports, including:

* **canonical** introns (the known annotation),
* **true cryptic-exon events** — biologically faithful: a cryptic exon inside an
  intron produces a ``novel_acceptor`` junction (sharing the upstream *known*
  donor) and a ``novel_donor`` junction (sharing the downstream *known*
  acceptor). These recur across replicates, use a canonical GT/AG motif, and are
  **up-regulated in condition B**, and
* **noise** novel junctions — sporadic, low-support, often non-canonical.

Every junction carries an ``is_cryptic_truth`` label, which is what the ML module
is trained to recover. Fully reproducible via ``seed``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

_NONCANONICAL = ["non-canonical", "GC/AG", "AT/AC", "GT/AT"]


#: Columns of :attr:`SimulatedDataset.truth`. Type-specific fields are NaN elsewhere,
#: the same shape :func:`splicescope.events.detect_events` returns.
TRUTH_COLUMNS = [
    "event_type", "gene_id", "chrom", "strand",
    "intron_start", "intron_end",
    "exonA_start", "exonA_end", "exonB_start", "exonB_end",
    "site_pos", "alt_pos",
]


@dataclass
class SimulatedDataset:
    known: pd.DataFrame
    observed: pd.DataFrame
    groups: dict[str, str] = field(default_factory=dict)
    #: One row per event actually injected — see :data:`TRUTH_COLUMNS`. Recovering this
    #: from the coordinates is not equivalent: a gene drawn for an event is skipped when
    #: its intron is too short to hold one, so the arithmetic gives every geometry the
    #: simulator *could* have used rather than the ones it did. Recall measured against
    #: the former is a lower bound, which is what the MXE test used to assert.
    truth: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=TRUTH_COLUMNS))


def simulate_dataset(
    n_genes: int = 8,
    exons_per_gene: int = 5,
    n_per_group: int = 4,
    cryptic_fraction: float = 0.5,
    alt_ss_fraction: float = 0.0,
    mxe_fraction: float = 0.0,
    label_noise: float = 0.0,
    seed: int = 0,
) -> SimulatedDataset:
    """Generate a reproducible splicing dataset with known ground truth.

    ``label_noise`` flips ``is_cryptic_truth`` for that fraction of junctions,
    mimicking imperfect curation. Left at 0 the labels are exact (deterministic);
    a small value (e.g. 0.1) makes the downstream ML task realistically hard.
    """
    rng = np.random.default_rng(seed)
    chrom, strand = "chr1", "+"
    exon_len, intron_len, gene_gap = 120, 500, 3000

    known_rows: list[tuple] = []
    gene_introns: dict[str, list[tuple[int, int]]] = {}
    pos = 1000
    for g in range(n_genes):
        gene_id = f"g{g:02d}"
        exons = []
        for _ in range(exons_per_gene):
            start, end = pos, pos + exon_len - 1
            exons.append((start, end))
            pos = end + intron_len + 1
        introns = []
        for (_, e1), (s2, _) in zip(exons[:-1], exons[1:], strict=False):
            i_start, i_end = e1 + 1, s2 - 1
            known_rows.append((chrom, i_start, i_end, strand, gene_id))
            introns.append((i_start, i_end))
        gene_introns[gene_id] = introns
        pos += gene_gap

    known = pd.DataFrame(known_rows, columns=["chrom", "start", "end", "strand", "gene_id"])

    samples = [f"A{i}" for i in range(n_per_group)] + [f"B{i}" for i in range(n_per_group)]
    groups = {s: ("A" if s.startswith("A") else "B") for s in samples}

    records: list[dict] = []
    truth_rows: list[dict] = []

    def note(event_type, gene_id, intron, **fields):
        truth_rows.append(
            {
                "event_type": event_type,
                "gene_id": gene_id,
                "chrom": chrom,
                "strand": strand,
                "intron_start": intron[0],
                "intron_end": intron[1],
                **fields,
            }
        )

    def emit(start, end, motif, truth, per_sample_counts):
        for s, c in per_sample_counts.items():
            if c <= 0:
                continue
            records.append(
                {
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "motif": motif,
                    "annotated_star": 0,
                    "count": int(c),
                    "sample": s,
                    "is_cryptic_truth": int(truth),
                }
            )

    # choose MXE introns up front so their canonical (skipping) junction can be
    # suppressed — mutually-exclusive exons have no exon-skipping isoform. Guarded
    # so that mxe_fraction == 0 draws no randomness and leaves the stream untouched.
    mxe_intron_map: dict[str, tuple[int, int]] = {}
    if mxe_fraction > 0:
        for gene_id, introns in gene_introns.items():
            if rng.random() >= mxe_fraction:
                continue
            candidates = [iv for iv in introns if iv[1] - iv[0] >= 300]
            if candidates:
                mxe_intron_map[gene_id] = candidates[rng.integers(len(candidates))]
    mxe_introns = set(mxe_intron_map.values())

    # One intron carries at most one injected event, so the truth table is unambiguous.
    # Two events in one intron do not merely crowd each other, they change what the reads
    # mean: an MXE intron has its skipping junction suppressed, so a cryptic exon placed
    # in it has no cassette to be detected as, and an alternative donor there is also a
    # leg of an MXE pair and is reported as that instead. Both are the right reading of
    # the data and both make the truth table claim events that are not there.
    taken: set[tuple[int, int]] = set(mxe_introns)

    # 1) canonical introns — well expressed everywhere (except MXE introns)
    canonical_base = {}
    for _, i_start, i_end, _, _ in known.itertuples(index=False):
        if (i_start, i_end) in mxe_introns:
            continue
        counts = {s: int(rng.poisson(200)) for s in samples}
        canonical_base[(i_start, i_end)] = counts
        emit(i_start, i_end, "GT/AG", 0, counts)

    # 2) true cryptic-exon events in a subset of genes (B-upregulated).
    # Effect sizes vary — some events are strong, some are subtle — so the
    # classification task is non-trivial rather than perfectly separable.
    cryptic_genes = [g for g in gene_introns if rng.random() < cryptic_fraction]
    for gene_id in cryptic_genes:
        free = [iv for iv in gene_introns[gene_id] if iv not in taken and iv[1] - iv[0] >= 200]
        if not free:
            continue
        i_start, i_end = free[rng.integers(len(free))]
        taken.add((i_start, i_end))
        c_start = i_start + 100
        c_end = c_start + 60
        strength = rng.uniform(0.35, 1.0)  # subtle -> strong
        b_lambda = 15 + 55 * strength
        incl = {
            s: int(rng.poisson(b_lambda if groups[s] == "B" else 6)) for s in samples
        }
        # a minority of genuine events use a non-canonical motif
        motif = "GT/AG" if rng.random() > 0.15 else _NONCANONICAL[rng.integers(len(_NONCANONICAL))]
        emit(i_start, c_start - 1, motif, 1, incl)  # novel_acceptor (shares known donor)
        emit(c_end + 1, i_end, motif, 1, incl)      # novel_donor (shares known acceptor)
        note("cryptic_exon", gene_id, (i_start, i_end), exonA_start=c_start, exonA_end=c_end)

    # 2b) alternative 5′/3′ splice-site events (B-upregulated alternative usage).
    alt_genes = [g for g in gene_introns if rng.random() < alt_ss_fraction]
    for gene_id in alt_genes:
        free = [iv for iv in gene_introns[gene_id] if iv not in taken and iv[1] - iv[0] >= 120]
        if not free:
            continue
        i_start, i_end = free[rng.integers(len(free))]
        taken.add((i_start, i_end))
        delta = int(rng.integers(20, 60))
        alt = {s: int(rng.poisson(50 if groups[s] == "B" else 12)) for s in samples}
        if rng.random() < 0.5:
            emit(i_start + delta, i_end, "GT/AG", 0, alt)  # A5SS: alt donor, shared acceptor
            note("A5SS", gene_id, (i_start, i_end), site_pos=i_end, alt_pos=i_start + delta)
        else:
            emit(i_start, i_end - delta, "GT/AG", 0, alt)  # A3SS: shared donor, alt acceptor
            note("A3SS", gene_id, (i_start, i_end), site_pos=i_start, alt_pos=i_end - delta)

    # 2c) mutually-exclusive-exon (MXE) events: exon A favoured in condition A,
    # exon B in condition B, between shared flanking exons (no skipping isoform).
    for gene_id, (i_start, i_end) in mxe_intron_map.items():
        a_start, a_end = i_start + 50, i_start + 90
        b_start, b_end = i_start + 150, i_start + 190
        if b_end + 1 >= i_end:
            continue
        incl_a = {s: int(rng.poisson(55 if groups[s] == "A" else 15)) for s in samples}
        incl_b = {s: int(rng.poisson(15 if groups[s] == "A" else 55)) for s in samples}
        emit(i_start, a_start - 1, "GT/AG", 0, incl_a)  # up -> exon A
        emit(a_end + 1, i_end, "GT/AG", 0, incl_a)      # exon A -> down
        emit(i_start, b_start - 1, "GT/AG", 0, incl_b)  # up -> exon B
        emit(b_end + 1, i_end, "GT/AG", 0, incl_b)      # exon B -> down
        note(
            "MXE", gene_id, (i_start, i_end),
            exonA_start=a_start, exonA_end=a_end, exonB_start=b_start, exonB_end=b_end,
        )

    # 3) noise novel junctions — mostly sporadic/low, but a fraction mimic real
    # events (canonical motif, recurrent support) so classes overlap, truth=0.
    n_noise = max(6, n_genes * 3)
    for _ in range(n_noise):
        gene_id = list(gene_introns)[rng.integers(n_genes)]
        introns = gene_introns[gene_id]
        i_start, i_end = introns[rng.integers(len(introns))]
        if i_end - i_start < 80:
            continue
        ns = i_start + int(rng.integers(20, max(21, i_end - i_start - 40)))
        ne = ns + int(rng.integers(20, 60))
        deceptive = rng.random() < 0.35
        if deceptive:
            # a hard negative: shares the intron's known donor (so it looks like a
            # cryptic novel_acceptor), canonical motif, recurrent moderate support
            motif = "GT/AG"
            hit = rng.choice(samples, size=int(rng.integers(2, len(samples))), replace=False)
            counts = {s: (int(rng.poisson(9)) if s in hit else 0) for s in samples}
            emit(i_start, ne, motif, 0, counts)  # known donor -> novel acceptor
        else:
            motif = _NONCANONICAL[rng.integers(len(_NONCANONICAL))]
            hit = rng.choice(samples, size=int(rng.integers(1, 3)), replace=False)
            counts = {s: (int(rng.poisson(3)) if s in hit else 0) for s in samples}
            emit(ns, ne, motif, 0, counts)  # both-novel sporadic artefact

    observed = pd.DataFrame.from_records(records)

    if label_noise > 0 and not observed.empty:
        jid = (
            observed["chrom"].astype(str)
            + ":"
            + observed["start"].astype(str)
            + "-"
            + observed["end"].astype(str)
            + observed["strand"].astype(str)
        )
        observed = observed.assign(_jid=jid)
        unique_ids = observed["_jid"].drop_duplicates().to_numpy()
        n_flip = int(round(len(unique_ids) * label_noise))
        if n_flip > 0:
            flip_ids = set(rng.choice(unique_ids, size=n_flip, replace=False))
            mask = observed["_jid"].isin(flip_ids)
            observed.loc[mask, "is_cryptic_truth"] = 1 - observed.loc[mask, "is_cryptic_truth"]
        observed = observed.drop(columns="_jid")

    truth = pd.DataFrame(truth_rows, columns=TRUTH_COLUMNS)
    # nullable integers, so a coordinate that does not apply to this event type is blank
    # rather than NaN — which would otherwise make every coordinate a float on the way out
    for column in TRUTH_COLUMNS[4:]:
        truth[column] = truth[column].astype("Int64")
    return SimulatedDataset(known=known, observed=observed, groups=groups, truth=truth)


_BASES = ("A", "C", "G", "T")
_STOP_CODONS = ("TAA", "TAG", "TGA")
#: The 61 sense codons. Coding exons are built from these so that every simulated
#: gene carries a genuine open reading frame and a stop only where one is intended.
_SENSE_CODONS = [
    a + b + c for a in _BASES for b in _BASES for c in _BASES if a + b + c not in _STOP_CODONS
]


def gene_exons(ds: SimulatedDataset) -> dict[str, list[tuple[int, int]]]:
    """Exon blocks per gene, reconstructed from the known introns (120 bp flanks)."""
    out: dict[str, list[tuple[int, int]]] = {}
    for gene_id, sub in ds.known.groupby("gene_id"):
        introns = sorted(zip(sub["start"], sub["end"], strict=False))
        exon_coords = []
        prev_end = introns[0][0] - 121
        for i_start, i_end in introns:
            exon_coords.append((prev_end + 1, i_start - 1))
            prev_end = i_end
        exon_coords.append((prev_end + 1, prev_end + 120))
        out[str(gene_id)] = exon_coords
    return out


def _mature_to_genomic(
    blocks: list[tuple[int, int]], offset: int, length: int
) -> list[tuple[int, int]]:
    """Genomic intervals covering ``[offset, offset+length)`` of the mature transcript."""
    out, pos = [], 0
    for s, e in blocks:
        n = e - s + 1
        lo, hi = max(offset, pos), min(offset + length, pos + n)
        if hi > lo:
            out.append((s + lo - pos, s + hi - 1 - pos))
        pos += n
    return out


def simulate_genome(
    ds: SimulatedDataset, utr5: int = 30, utr3: int = 60, seed: int = 0
) -> tuple[str, dict[str, list[tuple[int, int]]]]:
    """Build a synthetic chromosome consistent with the simulated annotation.

    Exonic sequence is drawn from sense codons, so every gene carries a real open
    reading frame that ends in a single stop; intron ends carry the canonical
    ``GT``/``AG`` dinucleotides; everything else is random. That is enough for the
    protein-consequence layer to be exercised end to end — a cryptic exon spliced
    into one of these transcripts genuinely shifts the frame and genuinely has (or
    has not) an in-frame stop, rather than being asserted to.

    Uses its own random stream, so the junction data in ``ds`` is untouched.
    Returns the chromosome sequence and the CDS blocks of each gene.
    """
    rng = np.random.default_rng(seed)
    exons_by_gene = gene_exons(ds)
    span = max(e for blocks in exons_by_gene.values() for _, e in blocks)
    seq = list(rng.choice(_BASES, size=span + 120))

    cds_blocks: dict[str, list[tuple[int, int]]] = {}
    for gene_index, (gene_id, blocks) in enumerate(exons_by_gene.items()):
        # Stagger the 5'UTR by one base per gene so the coding frame at an intron
        # boundary cycles through 0, 1 and 2. Exons here are all 120 bp, so a fixed
        # UTR would put every intron on a codon boundary and frame inheritance —
        # the subtle part of the consequence layer — would never be exercised.
        gene_utr5 = utr5 - 1 + gene_index % 3
        mature_len = sum(e - s + 1 for s, e in blocks)
        cds_len = max(0, (mature_len - gene_utr5 - utr3) // 3 * 3)
        if cds_len < 6:
            continue
        utr5_len = gene_utr5
        codons = [*rng.choice(_SENSE_CODONS, size=cds_len // 3 - 1), "TAA"]
        mature = (
            "".join(rng.choice(_BASES, size=utr5_len))
            + "".join(codons)
            + "".join(rng.choice(_BASES, size=mature_len - utr5_len - cds_len))
        )
        offset = 0
        for s, e in blocks:
            n = e - s + 1
            seq[s - 1 : e] = list(mature[offset : offset + n])
            offset += n
        cds_blocks[gene_id] = _mature_to_genomic(blocks, utr5_len, cds_len)

    for _, i_start, i_end, _, _ in ds.known.itertuples(index=False):
        seq[i_start - 1], seq[i_start] = "G", "T"  # donor GT
        seq[i_end - 2], seq[i_end - 1] = "A", "G"  # acceptor AG

    return "".join(seq), cds_blocks


def write_fasta(sequence: str, path: str | Path, name: str = "chr1", width: int = 60) -> Path:
    """Write ``sequence`` as a FASTA plus the ``.fai`` index :class:`GenomeFasta` needs."""
    path = Path(path)
    header = f">{name}\n"
    lines = [sequence[i : i + width] for i in range(0, len(sequence), width)]
    path.write_text(header + "\n".join(lines) + "\n")
    fai = path.with_suffix(path.suffix + ".fai")
    fai.write_text(f"{name}\t{len(sequence)}\t{len(header)}\t{width}\t{width + 1}\n")
    return path


def write_dataset(ds: SimulatedDataset, outdir: str | Path, seed: int = 0) -> Path:
    """Write a simulated dataset: GTF, per-sample SJ.out.tab, groups.tsv and a genome.

    The GTF carries CDS records and the genome is written with its ``.fai``, so the
    protein-consequence layer runs on the built-in data with no downloads.
    """
    outdir = Path(outdir)
    (outdir / "sj").mkdir(parents=True, exist_ok=True)

    genome, cds_blocks = simulate_genome(ds, seed=seed)
    write_fasta(genome, outdir / "genome.fa")

    gtf_lines = []
    for gene_id, exon_coords in gene_exons(ds).items():
        attrs = (
            f'gene_id "{gene_id}"; transcript_id "{gene_id}.t1"; gene_name "{gene_id}";'
        )
        for es, ee in exon_coords:
            gtf_lines.append(f"chr1\tsim\texon\t{es}\t{ee}\t.\t+\t.\t{attrs}")
        written = 0
        for cs, ce in cds_blocks.get(gene_id, []):
            phase = (3 - written % 3) % 3  # GTF frame: bases to remove to reach a codon start
            gtf_lines.append(f"chr1\tsim\tCDS\t{cs}\t{ce}\t.\t+\t{phase}\t{attrs}")
            written += ce - cs + 1
    (outdir / "annotation.gtf").write_text("\n".join(gtf_lines) + "\n")

    star_cols = [
        "chrom", "start", "end", "strand", "motif", "annotated", "n_unique", "n_multi", "oh",
    ]
    strand_code = {"+": 1, "-": 2, ".": 0}
    motif_code = {"GT/AG": 1, "CT/AC": 2, "GC/AG": 3, "AT/AC": 5, "non-canonical": 0, "GT/AT": 6}
    for sample, sub in ds.observed.groupby("sample"):
        rows = []
        for r in sub.itertuples(index=False):
            rows.append(
                [
                    r.chrom,
                    r.start,
                    r.end,
                    strand_code.get(r.strand, 0),
                    motif_code.get(r.motif, 0),
                    0,
                    r.count,
                    0,
                    30,
                ]
            )
        pd.DataFrame(rows, columns=star_cols).to_csv(
            outdir / "sj" / f"{sample}.SJ.out.tab", sep="\t", header=False, index=False
        )

    pd.DataFrame(
        {"sample": list(ds.groups), "condition": list(ds.groups.values())}
    ).to_csv(outdir / "groups.tsv", sep="\t", index=False)
    # the events that were injected, so recall can be measured against them rather
    # than guessed at from the coordinates
    ds.truth.to_csv(outdir / "truth.tsv", sep="\t", index=False)
    return outdir
