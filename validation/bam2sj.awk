# Extract splice junctions from SAM records into STAR SJ.out.tab format.
# Junctions are the N operations of the CIGAR; strand comes from HISAT2's XS:A
# tag and uniqueness from NH:i. Emits one row per junction:
#   chrom  intron_start  intron_end  strand_code  motif_code  annotated
#   n_unique  n_multi  max_overhang
{
    strand = 0
    nh = 1
    for (i = 12; i <= NF; i++) {
        if ($i ~ /^XS:A:/) { s = substr($i, 6); strand = (s == "+") ? 1 : ((s == "-") ? 2 : 0) }
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
            key = $3 "\t" istart "\t" iend "\t" strand
            if (nh == 1) uniq[key]++ ; else multi[key]++
            if (overhang > oh[key]) oh[key] = overhang
        }
        if (op == "M" || op == "D" || op == "N" || op == "=" || op == "X") pos += len
    }
}
END {
    for (k in oh) printf "%s\t0\t0\t%d\t%d\t%d\n", k, uniq[k] + 0, multi[k] + 0, oh[k]
}
