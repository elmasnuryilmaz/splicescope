# Extract splice junctions from SAM records into STAR SJ.out.tab format.
#
# Input is SAM records, as `samtools view -F 0x104` writes them; header lines are
# skipped if present. Output is one row per intron, which is what STAR writes and what
# splicescope's reader expects:
#   chrom  intron_start  intron_end  strand_code  motif_code  annotated
#   n_unique  n_multi  max_overhang
#
# Junctions are the N operations of the CIGAR, in 1-based inclusive intron coordinates.
# One read with two N operations contributes to two introns. A row is keyed by
# (chrom, intron_start, intron_end) alone, and every column is aggregated over all the
# reads that cross that intron:
#   n_unique      reads with NH:i:1 (a read without an NH tag counts as unique)
#   n_multi       reads with NH:i greater than 1
#   max_overhang  the largest, over those reads, of the shorter of the two M blocks
#                 immediately flanking the N
#   strand_code   1 (+) or 2 (-) from HISAT2's XS:A tag, decided by majority over the
#                 reads that carry the tag, unique and multi alike; 0 (undefined) when
#                 no read carries it or when + and - are tied
#   motif_code    always 0: the motif is a property of the reference, and this script
#                 reads no genome. `add_motifs.py` fills it in afterwards, and resets
#                 strand_code from the motif wherever the motif is canonical.
#   annotated     always 0: this script reads no annotation. splicescope does not use
#                 the column; it compares junctions with the GTF itself.
#
# Why the strand is not part of the key. HISAT2 puts XS:A on some of the reads that
# cross an intron and leaves it off others, so keying on the tag wrote one intron as two
# rows, one with strand 0 and one with 1 or 2, each holding part of its reads.
# splicescope reads strand 0 as "." and then gives a "." junction the strand its
# annotation implies, at which point the two rows are the same junction twice in one
# sample and `compute_psi` refuses the table. Before that refusal existed, each of the
# two rows carried only part of the intron's reads. STAR names each intron once, and so
# does this script.
#
# Why a majority, and 0 on a tie. A row is one intron and carries one strand. Where the
# reads crossing it disagree, the larger share is the better estimate, and the
# minority's reads still count towards n_unique and n_multi, because they did cross the
# intron. When the two sides are tied the reads give no answer, and 0 says so:
# splicescope then takes the strand from the annotation where it knows the intron or
# one of its splice sites, and leaves the junction unstranded where it does not, rather
# than guessing.
/^@/ { next }
{
    xs = 0
    nh = 1
    for (i = 12; i <= NF; i++) {
        if ($i ~ /^XS:A:/) { s = substr($i, 6); xs = (s == "+") ? 1 : ((s == "-") ? 2 : 0) }
        else if ($i ~ /^NH:i:/) { nh = substr($i, 6) + 0 }
    }

    cigar = $6
    if (cigar !~ /N/) next

    pos = $4                 # 1-based leftmost reference position
    n = 0                    # number of ops parsed
    delete oplen; delete opch
    while (match(cigar, /^[0-9]+[MIDNSHP=X]/)) {
        tok = substr(cigar, 1, RLENGTH)
        n++
        oplen[n] = substr(tok, 1, RLENGTH - 1) + 0
        opch[n] = substr(tok, RLENGTH, 1)
        cigar = substr(cigar, RLENGTH + 1)
    }

    for (i = 1; i <= n; i++) {
        op = opch[i]; len = oplen[i]
        if (op == "N") {
            istart = pos
            iend = pos + len - 1
            # Anchor length: the aligned block immediately flanking the intron.
            left = (i > 1 && opch[i-1] == "M") ? oplen[i-1] : 0
            right = (i < n && opch[i+1] == "M") ? oplen[i+1] : 0
            overhang = (left < right) ? left : right
            key = $3 "\t" istart "\t" iend
            if (nh == 1) uniq[key]++ ; else multi[key]++
            if (xs == 1) plus[key]++ ; else if (xs == 2) minus[key]++
            if (!(key in oh) || overhang > oh[key]) oh[key] = overhang
        }
        if (op == "M" || op == "D" || op == "N" || op == "=" || op == "X") pos += len
    }
}
END {
    for (k in oh) {
        p = plus[k] + 0; m = minus[k] + 0
        strand = (p > m) ? 1 : ((m > p) ? 2 : 0)
        printf "%s\t%d\t0\t0\t%d\t%d\t%d\n", k, strand, uniq[k] + 0, multi[k] + 0, oh[k]
    }
}
