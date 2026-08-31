#!/usr/bin/env bash
# Derive STAR-style SJ.out.tab files from coordinate-sorted BAMs.
# Usage: extract_sj.sh <bam_dir> <out_dir>
set -euo pipefail
BAMDIR="$1"; OUTDIR="$2"
AWK="$(dirname "$0")/bam2sj.awk"
mkdir -p "$OUTDIR"
for bam in "$BAMDIR"/*.bam; do
    s=$(basename "$bam" .sorted.bam); s=${s%.bam}
    samtools view -F 0x104 "$bam" \
      | awk -f "$AWK" \
      | sort -k1,1 -k2,2n > "$OUTDIR/${s}.SJ.out.tab"
    echo "done $s $(wc -l < "$OUTDIR/${s}.SJ.out.tab") junctions"
done
