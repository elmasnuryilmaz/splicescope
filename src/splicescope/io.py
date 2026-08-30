"""Reading splice junctions from STAR and reference annotations from GTF.

Coordinate convention (matches STAR ``SJ.out.tab``)
---------------------------------------------------
A junction is the *intron* it spans: ``start`` is the first intronic base and
``end`` the last intronic base, both 1-based inclusive. The 5' splice site
(donor) and 3' splice site (acceptor) therefore depend on strand:

* ``+`` strand: donor = ``start``, acceptor = ``end``
* ``-`` strand: donor = ``end``,   acceptor = ``start``

Every junction is keyed by ``(chrom, start, end, strand)``.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import pandas as pd

# STAR SJ.out.tab strand codes and intron-motif codes.
_STAR_STRAND = {0: ".", 1: "+", 2: "-"}
# STAR motif code 0 = non-canonical; odd codes are + strand, even are - strand.
_STAR_MOTIF = {
    0: "non-canonical",
    1: "GT/AG",
    2: "CT/AC",
    3: "GC/AG",
    4: "CT/GC",
    5: "AT/AC",
    6: "GT/AT",
}

JUNCTION_KEY = ["chrom", "start", "end", "strand"]


def donor_acceptor(start: int, end: int, strand: str) -> tuple[int, int]:
    """Return ``(donor_pos, acceptor_pos)`` for an intron given its strand."""
    if strand == "-":
        return end, start
    return start, end


def read_star_sj(path: str | Path, sample: str | None = None) -> pd.DataFrame:
    """Read a STAR ``SJ.out.tab`` file into a tidy junction table.

    Parameters
    ----------
    path : path to ``SJ.out.tab``.
    sample : sample label; defaults to the file's stem.

    Returns a DataFrame with columns
    ``[chrom, start, end, strand, motif, annotated_star, count, sample]`` where
    ``count`` is the number of uniquely-mapping reads crossing the junction.
    """
    path = Path(path)
    sample = sample or path.stem.replace(".SJ.out", "").replace("SJ.out", "") or path.stem
    cols = [
        "chrom",
        "start",
        "end",
        "strand_code",
        "motif_code",
        "annotated_star",
        "n_unique",
        "n_multi",
        "max_overhang",
    ]
    df = pd.read_csv(path, sep="\t", header=None, names=cols, dtype={"chrom": str})
    df["strand"] = df["strand_code"].map(_STAR_STRAND).fillna(".")
    df["motif"] = df["motif_code"].map(_STAR_MOTIF).fillna("non-canonical")
    df["count"] = df["n_unique"].astype("int64")
    df["sample"] = sample
    return df[["chrom", "start", "end", "strand", "motif", "annotated_star", "count", "sample"]]


def read_many_star_sj(paths: dict[str, str | Path]) -> pd.DataFrame:
    """Read several ``SJ.out.tab`` files into one long table.

    ``paths`` maps sample name -> file path.
    """
    frames = [read_star_sj(p, sample=name) for name, p in paths.items()]
    if not frames:
        cols = ["chrom", "start", "end", "strand", "motif", "annotated_star", "count", "sample"]
        return pd.DataFrame(columns=cols)
    return pd.concat(frames, ignore_index=True)


def _parse_attributes(field: str) -> dict[str, str]:
    """Parse the 9th GTF column into a dict (handles GTF ``key "value";`` style)."""
    out: dict[str, str] = {}
    for chunk in field.strip().split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        if " " in chunk:
            key, _, value = chunk.partition(" ")
            out[key] = value.strip().strip('"')
    return out


def read_gtf_junctions(path: str | Path) -> pd.DataFrame:
    """Derive the set of *known* introns (junctions) from a GTF annotation.

    Exons are grouped per transcript, sorted, and the gaps between consecutive
    exons become known introns. The result is a DataFrame with columns
    ``[chrom, start, end, strand, gene_id]`` (one row per unique junction).
    """
    path = Path(path)
    exons_by_tx: dict[str, list[tuple[int, int]]] = defaultdict(list)
    tx_meta: dict[str, tuple[str, str, str]] = {}  # tx -> (chrom, strand, gene_id)

    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "exon":
                continue
            chrom, _, _, start, end, _, strand, _, attrs = parts
            a = _parse_attributes(attrs)
            tx = a.get("transcript_id")
            if tx is None:
                continue
            exons_by_tx[tx].append((int(start), int(end)))
            tx_meta[tx] = (chrom, strand, a.get("gene_id", tx))

    rows: dict[tuple[str, int, int, str], str] = {}
    for tx, exons in exons_by_tx.items():
        chrom, strand, gene_id = tx_meta[tx]
        exons.sort()
        for (_, e_end), (n_start, _) in zip(exons[:-1], exons[1:], strict=False):
            intron_start = e_end + 1
            intron_end = n_start - 1
            if intron_end < intron_start:
                continue  # overlapping/degenerate exons
            rows[(chrom, intron_start, intron_end, strand)] = gene_id

    if not rows:
        return pd.DataFrame(columns=["chrom", "start", "end", "strand", "gene_id"])
    data = [(c, s, e, st, g) for (c, s, e, st), g in rows.items()]
    return pd.DataFrame(data, columns=["chrom", "start", "end", "strand", "gene_id"])


def read_gmt(path: str | Path) -> dict[str, list[str]]:
    """Read a GMT gene-set file into ``{term: [genes]}``.

    GMT lines are tab-separated: ``term<TAB>description<TAB>gene1<TAB>gene2...``.
    The description column is ignored. Compatible with MSigDB / GO / KEGG exports.
    """
    sets: dict[str, list[str]] = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            term = parts[0].strip()
            genes = [g.strip() for g in parts[2:] if g.strip()]
            if term and genes:
                sets[term] = genes
    return sets
