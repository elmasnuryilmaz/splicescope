"""Predict what a cryptic exon does to the protein its host transcript encodes.

Detecting a cryptic exon says nothing about whether it matters. The same event
can be silently tolerated, shift the reading frame, or introduce a premature
termination codon (PTC) that sends the transcript to nonsense-mediated decay --
which is how TDP-43 cryptic exons deplete proteins such as STMN2 and UNC13A.
This module resolves that, one event at a time:

1. locate the host transcripts whose intron contains the event,
2. inherit the reading frame from the coding sequence upstream of the intron,
3. translate the inserted sequence in that frame and look for a stop codon,
4. apply the 50-nucleotide rule to decide whether a PTC triggers NMD.

Sequence comes from an indexed genome FASTA read directly through its ``.fai``,
so nothing here needs a compiled sequence library.

Consequence classes
-------------------
======================  ==================================================
class                   meaning
======================  ==================================================
``ptc_nmd``             PTC introduced, far enough upstream to trigger NMD
``ptc_escape``          PTC introduced but predicted to escape NMD
``frameshift``          insert length not a multiple of 3, no PTC inside it
``in_frame_insertion``  multiple of 3 and no stop codon: protein extended
``utr_insertion``       event falls outside the coding sequence
``non_coding_host``     host transcript has no annotated CDS
``exon_truncation``     a novel splice site shortens an annotated exon
``no_host_transcript``  no annotated intron contains the event
======================  ==================================================

Two shapes of event are handled. A **cassette** cryptic exon is given directly as
an exon interval. A **splice-site shift** — a novel donor or acceptor paired with
an annotated site on the other side — is given as a junction, and the sequence it
adds to or removes from the neighbouring exon is derived from the annotation.
Splice-site shifts are the majority of cryptic events reported by junction-level
callers, so a layer that only understood cassettes would miss most of them.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})

#: Minimum distance (nt) from a PTC to the last exon-exon junction for NMD.
NMD_DISTANCE_RULE = 50

PTC_NMD = "ptc_nmd"
PTC_ESCAPE = "ptc_escape"
FRAMESHIFT = "frameshift"
IN_FRAME = "in_frame_insertion"
UTR_INSERTION = "utr_insertion"
EXON_TRUNCATION = "exon_truncation"
NON_CODING_HOST = "non_coding_host"
NO_HOST = "no_host_transcript"

#: Classes in which the event is expected to reduce functional protein.
PROTEIN_DISRUPTING = (PTC_NMD, PTC_ESCAPE, FRAMESHIFT, EXON_TRUNCATION)

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    """Reverse complement of a DNA string."""
    return seq.translate(_COMPLEMENT)[::-1]


class GenomeFasta:
    """Random access to an indexed FASTA using only its ``.fai`` sidecar.

    The index gives, per contig, the byte offset of its first base plus how many
    bases and how many bytes each line holds, which is enough to seek straight to
    any coordinate without reading the file.
    """

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        fai = self.path.with_suffix(self.path.suffix + ".fai")
        if not fai.exists():
            raise FileNotFoundError(f"missing FASTA index: {fai}")
        self.index: dict[str, tuple[int, int, int, int]] = {}
        for line in fai.read_text().splitlines():
            name, length, offset, linebases, linewidth = line.split("\t")[:5]
            self.index[name] = (int(length), int(offset), int(linebases), int(linewidth))
        self._handle = None

    def __enter__(self) -> GenomeFasta:
        self._handle = open(self.path, "rb")
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def close(self) -> None:
        if self._handle is not None:
            self._handle.close()
            self._handle = None

    def fetch(self, chrom: str, start: int, end: int, strand: str = "+") -> str:
        """Sequence for ``chrom:start-end``, 1-based inclusive, strand-aware."""
        if chrom not in self.index:
            return ""
        length, offset, linebases, linewidth = self.index[chrom]
        start = max(1, start)
        end = min(length, end)
        if end < start:
            return ""
        if self._handle is None:
            self._handle = open(self.path, "rb")
        begin = offset + (start - 1) // linebases * linewidth + (start - 1) % linebases
        stop = offset + (end - 1) // linebases * linewidth + (end - 1) % linebases
        self._handle.seek(begin)
        raw = self._handle.read(stop - begin + 1).decode("ascii")
        seq = "".join(raw.split()).upper()
        return reverse_complement(seq) if strand == "-" else seq


@dataclass
class Transcript:
    """Exon and CDS structure of one transcript, in transcription order."""

    transcript_id: str
    gene_id: str
    gene_name: str
    chrom: str
    strand: str
    exons: list[tuple[int, int]]
    cds: list[tuple[int, int]] = field(default_factory=list)

    def __post_init__(self) -> None:
        reverse = self.strand == "-"
        self.exons = sorted(self.exons, reverse=reverse)
        self.cds = sorted(self.cds, reverse=reverse)

    @property
    def introns(self) -> list[tuple[int, int]]:
        """Introns in transcription order, 1-based inclusive."""
        blocks = sorted(self.exons)
        out = [(a[1] + 1, b[0] - 1) for a, b in zip(blocks, blocks[1:], strict=False)]
        return out[::-1] if self.strand == "-" else out

    def cds_length_before(self, position: int) -> int:
        """Length of coding sequence transcribed before ``position``."""
        total = 0
        for cstart, cend in self.cds:
            if self.strand == "+":
                if cend < position:
                    total += cend - cstart + 1
                elif cstart <= position:
                    total += position - cstart
            else:
                if cstart > position:
                    total += cend - cstart + 1
                elif cend >= position:
                    total += cend - position
        return total


def _strip_version(gene_id: str) -> str:
    """``ENSG00000130477.16`` -> ``ENSG00000130477``."""
    return str(gene_id).split(".", 1)[0]


def _attr(attributes: str, key: str) -> str:
    marker = key + ' "'
    i = attributes.find(marker)
    if i < 0:
        return ""
    i += len(marker)
    j = attributes.find('"', i)
    return attributes[i:j]


def load_transcripts(
    gtf_path: str | Path, genes: set[str] | None = None
) -> dict[str, Transcript]:
    """Read exon and CDS records from a GTF into transcript models.

    The file may be plain or gzipped. ``genes`` restricts parsing to a set of
    gene names or IDs, which matters:
    a full GENCODE GTF holds millions of records and a cryptic exon report
    normally concerns a few hundred genes.
    """
    # Callers pass gene names or Ensembl IDs, versioned or not; accept all three.
    wanted: set[str] | None = None
    if genes:
        wanted = set()
        for gene in genes:
            upper = str(gene).upper()
            wanted.add(upper)
            wanted.add(upper.split(".")[0])
    records: dict[str, dict] = {}

    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in ("exon", "CDS"):
                continue
            attributes = fields[8]
            gene_name = _attr(attributes, "gene_name")
            gene_id = _attr(attributes, "gene_id")
            if wanted is not None and not (
                gene_name.upper() in wanted or gene_id.split(".")[0].upper() in wanted
            ):
                continue
            tid = _attr(attributes, "transcript_id")
            if not tid:
                continue
            entry = records.setdefault(
                tid,
                {
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "chrom": fields[0],
                    "strand": fields[6],
                    "exon": [],
                    "CDS": [],
                },
            )
            entry[fields[2]].append((int(fields[3]), int(fields[4])))

    return {
        tid: Transcript(
            transcript_id=tid,
            gene_id=entry["gene_id"],
            gene_name=entry["gene_name"],
            chrom=entry["chrom"],
            strand=entry["strand"],
            exons=entry["exon"],
            cds=entry["CDS"],
        )
        for tid, entry in records.items()
        if entry["exon"]
    }


def index_by_gene(transcripts: dict[str, Transcript]) -> dict[str, list[Transcript]]:
    """Group transcripts by gene, under both the versioned and bare gene id."""
    index: dict[str, list[Transcript]] = {}
    for tx in transcripts.values():
        for key in {tx.gene_id, _strip_version(tx.gene_id), tx.gene_name.upper()}:
            if key:
                index.setdefault(key, []).append(tx)
    return index


def host_transcripts(
    transcripts: dict[str, Transcript] | list[Transcript],
    chrom: str,
    start: int,
    end: int,
    strand: str,
) -> list[Transcript]:
    """Transcripts with an intron fully containing ``chrom:start-end``.

    Pass only the candidate transcripts — the host gene's, via
    :func:`index_by_gene` — when the caller knows which gene the event belongs
    to; scanning every transcript for every event does not scale.
    """
    candidates = transcripts.values() if isinstance(transcripts, dict) else transcripts
    hosts = []
    for tx in candidates:
        if tx.chrom != chrom or tx.strand != strand:
            continue
        for istart, iend in tx.introns:
            if istart <= start and end <= iend:
                hosts.append(tx)
                break
    return hosts


def downstream_sequence(
    tx: Transcript, genome: GenomeFasta, chrom: str, resume: int
) -> tuple[str, list[int]]:
    """Mature sequence retained from ``resume`` onward, and its exon lengths.

    ``resume`` is the genomic position at which transcription continues after
    the event. Exons are walked in transcription order and clipped to that
    point, so the result is the sequence a ribosome would read downstream, with
    the length of each remaining exon kept so the last exon-exon junction can be
    located.
    """
    pieces: list[str] = []
    lengths: list[int] = []
    for estart, eend in tx.exons:  # already in transcription order
        if tx.strand == "+":
            if eend < resume:
                continue
            start = max(estart, resume)
            end = eend
        else:
            if estart > resume:
                continue
            start = estart
            end = min(eend, resume)
        if end < start:
            continue
        pieces.append(genome.fetch(chrom, start, end, tx.strand))
        lengths.append(end - start + 1)
    return "".join(pieces), lengths


def _nmd_from_downstream(
    sequence: str, lengths: list[int], frame: int, offset_before: int
) -> tuple[int | None, int | None, bool]:
    """Find the first in-frame stop downstream and apply the 50-nucleotide rule.

    ``offset_before`` is how many transcript bases precede ``sequence``, so the
    reported PTC offset stays on the same scale as the event itself.
    """
    ptc = find_ptc(sequence, frame)
    if ptc is None:
        return None, None, False
    last_junction = sum(lengths[:-1]) if len(lengths) > 1 else None
    if last_junction is None:
        return offset_before + ptc, None, False
    distance = last_junction - (ptc + 3)
    return offset_before + ptc, distance, distance > NMD_DISTANCE_RULE


def junction_change(tx: Transcript, start: int, end: int) -> tuple[str, int, int] | None:
    """What a novel junction adds to, or removes from, a neighbouring exon.

    ``start`` and ``end`` are the first and last intronic base, 1-based
    inclusive. The junction has to share one splice site with an annotated
    intron of ``tx``; the other site is the novel one, and the sequence between
    the annotated and the novel position is what changes.

    Returns ``("extension" | "truncation", first_base, last_base)``, or ``None``
    when neither site matches the annotation (an unanchored junction, which the
    annotation cannot interpret).
    """
    for istart, iend in tx.introns:
        if istart == start and iend != end:
            # The 3' end of the intron moved.
            if end < iend:
                return "extension", end + 1, iend
            return "truncation", iend + 1, end
        if iend == end and istart != start:
            # The 5' end of the intron moved.
            if start > istart:
                return "extension", istart, start - 1
            return "truncation", start, istart - 1
    return None


def junction_hosts(
    transcripts: dict[str, Transcript] | list[Transcript],
    chrom: str,
    start: int,
    end: int,
    strand: str,
) -> list[tuple[Transcript, str, int, int]]:
    """Transcripts whose annotation can interpret this junction, with the change."""
    candidates = transcripts.values() if isinstance(transcripts, dict) else transcripts
    hosts = []
    for tx in candidates:
        if tx.chrom != chrom or tx.strand != strand:
            continue
        change = junction_change(tx, start, end)
        if change is not None:
            hosts.append((tx, *change))
    return hosts


def predict_junction_consequence(
    tx: Transcript,
    genome: GenomeFasta,
    chrom: str,
    start: int,
    end: int,
) -> Consequence | None:
    """Predict the effect of a novel junction, via the exon change it implies.

    An extension adds sequence to an exon and is evaluated exactly like a
    cassette exon. A truncation removes sequence, which cannot introduce a stop
    inside the removed interval, so only the frame consequence is reported.
    """
    change = junction_change(tx, start, end)
    if change is None:
        return None
    kind, cstart, cend = change
    if kind == "extension":
        return predict_consequence(tx, genome, chrom, cstart, cend)

    length = cend - cstart + 1
    bounds = _cds_bounds(tx)
    base = dict(
        transcript_id=tx.transcript_id,
        gene_name=tx.gene_name,
        insert_length=-length,
        frame_offset=0,
        frameshift=length % 3 != 0,
        ptc_offset=None,
        distance_to_last_junction=None,
        nmd_predicted=False,
    )
    if bounds is None:
        return Consequence(**base, consequence_class=NON_CODING_HOST)
    if cend < bounds[0] or cstart > bounds[1]:
        return Consequence(**base, consequence_class=UTR_INSERTION)

    boundary = cstart if tx.strand == "+" else cend
    frame = tx.cds_length_before(boundary) % 3
    base["frame_offset"] = frame
    # Removing bases shifts everything downstream; the stop, if any, is there.
    resume = cend + 1 if tx.strand == "+" else cstart - 1
    tail, lengths = downstream_sequence(tx, genome, chrom, resume)
    offset, distance, nmd = _nmd_from_downstream(tail, lengths, frame, 0)
    if offset is None:
        return Consequence(**base, consequence_class=EXON_TRUNCATION)
    base["ptc_offset"] = offset
    base["distance_to_last_junction"] = distance
    base["nmd_predicted"] = nmd
    return Consequence(**base, consequence_class=PTC_NMD if nmd else PTC_ESCAPE)


@dataclass
class Consequence:
    """Predicted protein-level effect of including one cryptic exon."""

    transcript_id: str
    gene_name: str
    insert_length: int
    frame_offset: int
    frameshift: bool
    ptc_offset: int | None
    distance_to_last_junction: int | None
    nmd_predicted: bool
    consequence_class: str

    def as_dict(self) -> dict:
        return {
            "transcript_id": self.transcript_id,
            "gene_name": self.gene_name,
            "insert_length": self.insert_length,
            "frame_offset": self.frame_offset,
            "frameshift": self.frameshift,
            "ptc_offset": self.ptc_offset,
            "distance_to_last_junction": self.distance_to_last_junction,
            "nmd_predicted": self.nmd_predicted,
            "consequence_class": self.consequence_class,
        }


def _downstream_exon_lengths(tx: Transcript, start: int, end: int) -> list[int]:
    """Lengths of the exons transcribed after the intron holding the event."""
    blocks = tx.exons
    out = []
    for estart, eend in blocks:
        after = estart > end if tx.strand == "+" else eend < start
        if after:
            out.append(eend - estart + 1)
    return out


def _resume_after(tx: Transcript, start: int, end: int) -> int | None:
    """Genomic position where transcription continues after an event."""
    if tx.strand == "+":
        later = [s for s, _ in tx.exons if s > end]
        return min(later) if later else None
    earlier = [e for _, e in tx.exons if e < start]
    return max(earlier) if earlier else None


def _cds_bounds(tx: Transcript) -> tuple[int, int] | None:
    if not tx.cds:
        return None
    starts = [c[0] for c in tx.cds]
    ends = [c[1] for c in tx.cds]
    return min(starts), max(ends)


def find_ptc(sequence: str, frame_offset: int) -> int | None:
    """Offset of the first in-frame stop codon, or ``None``.

    ``frame_offset`` is how many bases of the codon spanning the exon boundary
    already lie upstream, so the first complete codon inside ``sequence`` starts
    at ``(3 - frame_offset) % 3``.
    """
    begin = (3 - frame_offset) % 3
    for i in range(begin, len(sequence) - 2, 3):
        if sequence[i : i + 3] in STOP_CODONS:
            return i
    return None


def predict_consequence(
    tx: Transcript,
    genome: GenomeFasta,
    chrom: str,
    start: int,
    end: int,
) -> Consequence:
    """Predict the effect of inserting ``chrom:start-end`` into ``tx``."""
    insert_length = end - start + 1
    bounds = _cds_bounds(tx)
    base = dict(
        transcript_id=tx.transcript_id,
        gene_name=tx.gene_name,
        insert_length=insert_length,
        frame_offset=0,
        frameshift=insert_length % 3 != 0,
        ptc_offset=None,
        distance_to_last_junction=None,
        nmd_predicted=False,
    )

    if bounds is None:
        return Consequence(**base, consequence_class=NON_CODING_HOST)
    cds_start, cds_end = bounds
    if end < cds_start or start > cds_end:
        return Consequence(**base, consequence_class=UTR_INSERTION)

    boundary = start if tx.strand == "+" else end
    frame_offset = tx.cds_length_before(boundary) % 3
    base["frame_offset"] = frame_offset

    sequence = genome.fetch(chrom, start, end, tx.strand)
    ptc_offset = find_ptc(sequence, frame_offset)

    if ptc_offset is None:
        if not base["frameshift"]:
            return Consequence(**base, consequence_class=IN_FRAME)
        # The frame is shifted but the insert itself carries no stop, so the
        # first premature stop lies downstream in the retained exons.
        resume = _resume_after(tx, start, end)
        if resume is None:
            return Consequence(**base, consequence_class=FRAMESHIFT)
        tail, lengths = downstream_sequence(tx, genome, chrom, resume)
        offset, distance, nmd = _nmd_from_downstream(
            tail, lengths, (frame_offset + insert_length) % 3, insert_length
        )
        if offset is None:
            return Consequence(**base, consequence_class=FRAMESHIFT)
        base["ptc_offset"] = offset
        base["distance_to_last_junction"] = distance
        base["nmd_predicted"] = nmd
        return Consequence(**base, consequence_class=PTC_NMD if nmd else PTC_ESCAPE)

    base["ptc_offset"] = ptc_offset
    downstream = _downstream_exon_lengths(tx, start, end)
    # Distance from the PTC to the final exon-exon junction of the transcript:
    # what is left of the cryptic exon, plus every downstream exon but the last.
    distance = (insert_length - ptc_offset - 3) + sum(downstream[:-1])
    base["distance_to_last_junction"] = distance
    nmd = bool(downstream) and distance > NMD_DISTANCE_RULE
    base["nmd_predicted"] = nmd
    return Consequence(**base, consequence_class=PTC_NMD if nmd else PTC_ESCAPE)


#: Most disruptive first; used to pick one prediction when several transcripts host an event.
_PRIORITY = {
    PTC_NMD: 0,
    FRAMESHIFT: 1,
    EXON_TRUNCATION: 2,
    PTC_ESCAPE: 3,
    IN_FRAME: 4,
    UTR_INSERTION: 5,
    NON_CODING_HOST: 6,
    NO_HOST: 7,
}


def annotate_junction_consequences(
    events: pd.DataFrame,
    transcripts: dict[str, Transcript],
    genome: GenomeFasta,
    chrom_col: str = "chrom",
    start_col: str = "intron_start",
    end_col: str = "intron_end",
    strand_col: str = "strand",
    gene_col: str | None = "gene_id",
) -> pd.DataFrame:
    """Add consequence columns to a table of novel junctions.

    Coordinates are the first and last intronic base. Each junction is
    interpreted against its host gene's annotation, and where several
    transcripts interpret it differently the most disruptive reading is kept.
    """
    by_gene = index_by_gene(transcripts) if gene_col in events.columns else None
    columns = [chrom_col, start_col, end_col, strand_col]
    if by_gene is not None:
        columns.append(gene_col)

    rows = []
    for record in events[columns].itertuples(index=False):
        chrom, start, end, strand = record[0], int(record[1]), int(record[2]), record[3]
        if by_gene is not None:
            gene = str(record[4])
            candidates = by_gene.get(gene) or by_gene.get(_strip_version(gene)) or []
        else:
            candidates = transcripts
        hosts = junction_hosts(candidates, chrom, start, end, strand)
        predictions = [
            p
            for tx, _, _, _ in hosts
            if (p := predict_junction_consequence(tx, genome, chrom, start, end)) is not None
        ]
        if not predictions:
            rows.append(
                {
                    "transcript_id": "",
                    "gene_name": "",
                    "insert_length": 0,
                    "frame_offset": 0,
                    "frameshift": False,
                    "ptc_offset": None,
                    "distance_to_last_junction": None,
                    "nmd_predicted": False,
                    "consequence_class": NO_HOST,
                    "n_host_transcripts": 0,
                }
            )
            continue
        predictions.sort(key=lambda c: _PRIORITY[c.consequence_class])
        best = predictions[0].as_dict()
        best["n_host_transcripts"] = len(predictions)
        rows.append(best)

    return pd.concat([events.reset_index(drop=True), pd.DataFrame(rows)], axis=1)


def annotate_consequences(
    events: pd.DataFrame,
    transcripts: dict[str, Transcript],
    genome: GenomeFasta,
    chrom_col: str = "chrom",
    start_col: str = "start",
    end_col: str = "end",
    strand_col: str = "strand",
    gene_col: str | None = "gene_id",
) -> pd.DataFrame:
    """Add consequence columns to a table of cryptic exon candidates.

    Where several transcripts host the event the most disruptive prediction is
    reported, since that is the one that determines whether protein is lost.
    ``n_host_transcripts`` records how many were considered.

    ``gene_col`` names a column holding the event's gene, which restricts the
    transcript search to that gene. Without it every transcript is scanned for
    every event, which is quadratic and impractical at genome scale.
    """
    priority = _PRIORITY
    by_gene = index_by_gene(transcripts) if gene_col in events.columns else None
    columns = [chrom_col, start_col, end_col, strand_col]
    if by_gene is not None:
        columns.append(gene_col)

    rows = []
    cache: dict[tuple, list[Consequence]] = {}
    for record in events[columns].itertuples(index=False):
        chrom, start, end, strand = record[0], record[1], record[2], record[3]
        if by_gene is not None:
            gene = str(record[4])
            candidates = by_gene.get(gene) or by_gene.get(_strip_version(gene)) or []
        else:
            candidates = transcripts
        hosts = host_transcripts(candidates, chrom, int(start), int(end), strand)
        if not hosts:
            rows.append(
                {
                    "transcript_id": "",
                    "gene_name": "",
                    "insert_length": int(end) - int(start) + 1,
                    "frame_offset": 0,
                    "frameshift": False,
                    "ptc_offset": None,
                    "distance_to_last_junction": None,
                    "nmd_predicted": False,
                    "consequence_class": NO_HOST,
                    "n_host_transcripts": 0,
                }
            )
            continue
        key = (chrom, int(start), int(end), strand)
        predictions = cache.get(key)
        if predictions is None:
            predictions = [
                predict_consequence(tx, genome, chrom, int(start), int(end)) for tx in hosts
            ]
            cache[key] = predictions
        predictions.sort(key=lambda c: priority[c.consequence_class])
        best = predictions[0].as_dict()
        best["n_host_transcripts"] = len(hosts)
        rows.append(best)

    return pd.concat(
        [events.reset_index(drop=True), pd.DataFrame(rows)], axis=1
    )


__all__ = [
    "Consequence",
    "GenomeFasta",
    "Transcript",
    "annotate_consequences",
    "annotate_junction_consequences",
    "find_ptc",
    "index_by_gene",
    "junction_change",
    "junction_hosts",
    "predict_junction_consequence",
    "host_transcripts",
    "load_transcripts",
    "predict_consequence",
    "reverse_complement",
]
