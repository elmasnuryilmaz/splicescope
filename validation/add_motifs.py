#!/usr/bin/env python3
"""Fill in the intron-motif column of SJ.out.tab files from an indexed genome.

`bam2sj.awk` recovers coordinates, strand and counts from a BAM, but the intron
motif is a property of the reference rather than the alignment, so it is left as
0 (non-canonical). Downstream features such as `canonical_motif` are useless
until it is filled in, which is what this script does — reading the two bases at
each end of every intron and assigning STAR's motif code.

    python validation/add_motifs.py --sj-dir sj/ --genome genome.fa
"""

from __future__ import annotations

import argparse
from pathlib import Path

from splicescope.consequence import GenomeFasta

#: STAR motif codes, keyed by the (donor, acceptor) dinucleotide on the + strand.
MOTIF_CODES = {
    ("GT", "AG"): 1,
    ("CT", "AC"): 2,
    ("GC", "AG"): 3,
    ("CT", "GC"): 4,
    ("AT", "AC"): 5,
    ("GT", "AT"): 6,
}
#: Odd codes are + strand introns, even codes are - strand.
STRAND_OF_MOTIF = {1: 1, 3: 1, 5: 1, 2: 2, 4: 2, 6: 2}


def annotate_file(path: Path, genome: GenomeFasta, fix_strand: bool = True) -> dict:
    stats = {"rows": 0, "canonical": 0, "strand_fixed": 0}
    out_lines = []
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) < 9:
            continue
        chrom, start, end = fields[0], int(fields[1]), int(fields[2])
        donor = genome.fetch(chrom, start, start + 1)
        acceptor = genome.fetch(chrom, end - 1, end)
        code = MOTIF_CODES.get((donor, acceptor), 0)
        fields[4] = str(code)
        # The motif also determines strand, which is more reliable than a missing
        # or ambiguous aligner tag.
        if fix_strand and code and fields[3] != str(STRAND_OF_MOTIF[code]):
            fields[3] = str(STRAND_OF_MOTIF[code])
            stats["strand_fixed"] += 1
        stats["rows"] += 1
        stats["canonical"] += code > 0
        out_lines.append("\t".join(fields))
    path.write_text("\n".join(out_lines) + "\n")
    return stats


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sj-dir", required=True, type=Path)
    parser.add_argument("--genome", required=True, type=Path)
    parser.add_argument("--keep-strand", action="store_true", help="do not correct strand")
    args = parser.parse_args()

    with GenomeFasta(args.genome) as genome:
        for path in sorted(args.sj_dir.glob("*.SJ.out.tab")):
            s = annotate_file(path, genome, fix_strand=not args.keep_strand)
            share = s["canonical"] / s["rows"] if s["rows"] else 0
            print(
                f"{path.name}: {s['rows']:,} junctions, "
                f"{share:.1%} canonical, {s['strand_fixed']:,} strands corrected"
            )


if __name__ == "__main__":
    main()
