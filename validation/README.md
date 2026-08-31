# Real-data validation

The bundled demo dataset is simulated, which proves the code runs but not that it
recovers biology. This directory records a validation on a real experiment.

## Dataset

**GSE245332** (PRJNA1027826) — TDP-43 knockdown, human, GRCh38, three replicates per group.

| group | runs |
|---|---|
| TDP-43 KD (shT2) | SRR26383188, SRR26383189, SRR26383190 |
| control (shNC3) | SRR26383185, SRR26383186, SRR26383187 |

Reads were aligned with HISAT2; the BAMs total ~24 GB.

## Getting splice junctions from BAM

`splicescope` reads STAR `SJ.out.tab`. When only alignments are available,
`extract_sj.sh` derives the same table from any coordinate-sorted BAM: junctions are the
`N` operations of the CIGAR, strand comes from the `XS:A` tag and uniqueness from `NH:i`.

```bash
./validation/extract_sj.sh /path/to/bam_dir /path/to/sj_dir
```

Six BAMs yielded 1,977,980 junction records (848,749 distinct junctions).

## Running

```bash
splicescope run \
    --sj-dir  sj_dir \
    --gtf     gencode.v47.basic.annotation.gtf \
    --groups  sj_dir/groups.tsv \
    --genome  GRCh38.primary_assembly.genome.fa \
    --outdir  results
```

Runtime for the whole pipeline on these six samples — junctions in, consequences out —
is **~1.5 minutes** on a laptop (2 GB peak). It was ~45 minutes before the quadratic
steps in MXE detection and consequence lookup were replaced; both changes are exact, and
the outputs are unchanged apart from row order.

## What it showed

### The rank test could not work

| | junction tests | smallest p | smallest q | significant |
|---|---:|---:|---:|---:|
| Mann–Whitney (≤ 0.6.0) | 138,643 | 0.0594 | 1.00 | **0** |
| beta-binomial (0.7.0) | 189,700 | ~1e-132 | ~1e-127 | **1,925** |

Not a threshold problem: with three replicates per group no rank test can return a
p-value below 0.1, so nothing survives BH across ~10⁵ junctions. Meanwhile 46 junctions
had |ΔΨ| ≥ 0.5 — the signal was present the whole time.

At event level the count-based test calls **4,699** differentially spliced events, against
**3,724** SE events called by rMATS on the same BAMs: the same order of magnitude from an
independent implementation.

### Known TDP-43 cryptic exons are recovered

Significant at q ≤ 0.05 and |ΔΨ| ≥ 0.1:

| gene | smallest q | largest \|ΔΨ\| |
|---|---:|---:|
| STMN2 | 9.7e-36 | 0.500 |
| PFKP | 9.4e-28 | 0.510 |
| HDGFL2 | 3.0e-14 | 0.700 |
| ARHGAP32 | 1.2e-12 | 0.665 |
| ATG4B | 3.8e-08 | 0.258 |
| AGRN | 5.2e-07 | 0.474 |
| KALRN | 1.6e-06 | 0.510 |
| RAP1GAP | 1.5e-02 | 0.469 |

Not called, with reasons: **UNC13A** (26 events tested, best q 0.31 — expressed too
weakly here for its cryptic exon to be assessed), **CAMK2B** (|ΔΨ| = 0.502 but only 5
events and 11 junctions, so the count model correctly declines to call it), **SYNJ2**
(q 0.12), **ELAVL3** (|ΔΨ| = 0.09, below the effect threshold).

CAMK2B is worth keeping in view: a large ratio on a handful of reads is exactly what a
count-based test should refuse to call, and what a Ψ-only test would have no way to judge.

### What the consequence layer showed — and did not

Run over all 296,034 cassette exons with the full GENCODE annotation:

| | events with a coding host | protein-disrupting | of those, `ptc_nmd` |
|---|---:|---:|---:|
| significant (q ≤ 0.05, \|ΔΨ\| ≥ 0.1) | 1,322 | 61.6% | 53.6% |
| all tested events | 116,387 | 66.9% | 58.4% |

**The consequence class does not distinguish real events from background** (OR 0.79,
Fisher p 5.6e-05 — if anything slightly depleted). The reason is arithmetic: a random
128 nt interval read in a fixed frame contains a stop codon with probability
`1 − (61/64)^42 ≈ 0.87`. Almost any intronic sequence looks disruptive, so a PTC call
carries little information about whether an event is genuine.

Use the consequence layer to **interpret events that are already significant**, not to
rank candidates. That is how it is applied above, and it is what makes STMN2 and HDGFL2
interpretable rather than merely detected.

One genuinely informative contrast did appear: significant events are far **less** likely
to have an annotated host intron at all (39.7% vs 62.5%), which is what one expects if
they are truly cryptic — sitting in sequence the annotation does not model as a
transcribed intron.

### Does the classifier add anything?

A second experiment, **GSE122069** (human, 3 vs 3, same direction convention), gives a
label that is not circular: an event or junction is "replicating" if it is also called in
that independent experiment, with the same sign. 140,401 junctions were tested in both.

Ranking by each candidate, scored against that label (classifier trained with
gene-held-out 5-fold CV, so no gene appears in both train and test):

| ranker | strict label (96 pos.) | loose label (1,578 pos.) |
|---|---:|---:|
| beta-binomial −log10(q) | **0.062** | 0.040 |
| \|ΔΨ\| | 0.020 | 0.049 |
| cryptic classifier | 0.006 | 0.051 |
| classifier + statistic | 0.048 | **0.057** |
| random | 0.0007 | 0.011 |

(average precision; the strict label is `q ≤ 0.05, |ΔΨ| ≥ 0.1` in the replication set,
the loose one `p < 0.05, |ΔΨ| ≥ 0.05`.)

**The answer depends on the label, and that is the finding.** For confidently replicating
events — the ones that matter — the differential statistic dominates and the classifier
is 11× worse; adding it to the statistic makes the ranking worse, not better. Only under
a permissive definition of replication does the classifier become competitive.

This is what its features predict. `intron_length`, `log_max_count`, `canonical_motif`
and distance to known splice sites describe whether a junction is *plausible*, not
whether it *changes between conditions* — none of them sees the group structure at all.
That is a useful thing to know about a junction, and a different question from the one
the differential test answers. It should not be presented as the way to find cryptic
events.

Caveats: 96 positives under the strict label is modest, and GSE122069 carries a weaker
effect than GSE245332 (261 vs 4,699 differential events), so replication is a
conservative proxy that will miss real events rather than invent them.

### Predicting NMD sensitivity

Whether a cryptic exon is actually degraded is not something our own statistics can
answer, so the label comes from a different kind of experiment: cryptic exons in
TDP-43-depleted neurons re-measured after NMD is blocked. Two independent studies,
different cell models, different ways of blocking decay.

**Cassette exons** — i3Neurons, NMD factors knocked down (doi:10.1101/2025.06.28.661837,
357 exons with hg38 coordinates and a quantitative ΔPSI label):

| predicted class | n | median ΔPSI when NMD is blocked |
|---|---:|---:|
| `ptc_nmd` | 275 | **12.70** |
| `ptc_escape` | 9 | 8.60 |
| `in_frame_insertion` | 31 | 6.25 |
| `utr_insertion` | 31 | 4.17 |
| `non_coding_host` | 6 | **−1.57** |

Mann-Whitney one-sided p = 3.8e-06, AUROC 0.675 (95% CI 0.596–0.746, 2000 bootstraps).
The last row is the internal control: a transcript with no coding sequence cannot carry a
premature stop, and shows no gain.

**Splice-site shifts** — iNeurons and iMNs, SMG1 inhibitor
(doi:10.1101/2025.07.09.664014, 3,419 interpretable junctions, binary label):
`ptc_nmd` is 48.5% NMD-sensitive against 23.5% where no stop is possible,
odds ratio 3.07, Fisher p = 2.3e-12.

Robustness on the first dataset: the result holds under all six ways of defining the
label (each NMD knockdown alone, their mean, their maximum; AUROC 0.65–0.68, every one
significant), and is unchanged when restricted to one exon per gene (AUROC 0.687).

**What this does and does not support.** Sensitivity is high (0.85–0.88) and specificity
low (0.26–0.46): most events are called `ptc_nmd`, because a random intronic interval
read in a fixed frame contains a stop with probability ≈0.87. Use the call to interpret
and prioritise events, not to filter them. The negative set in the second dataset is also
imperfect — it is defined as "not in the NMD-sensitive table", which was itself filtered
at ΔPSI ≥ 0.1, so mildly sensitive events sit among the negatives and bias toward the
null.

**A coordinate trap worth recording.** The junction table reports the flanking exon
boundaries, not the intron: the true intron is `[start + 1, end − 1]`. Taken literally,
interpretable coverage is 0%; corrected, it is 96%. The offset was found by measuring the
distance from each junction to the nearest annotated intron, which is systematically ±1.

## Reproducing

Raw data are public (SRA). Everything after alignment is in this repository; the
reference GTF and genome are the standard GENCODE and GRCh38 primary assembly files.
