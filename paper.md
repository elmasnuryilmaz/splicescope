---
title: 'splicescope: detecting cryptic splicing from splice junctions and predicting its effect on the encoded protein'
tags:
  - Python
  - bioinformatics
  - RNA-seq
  - alternative splicing
  - cryptic exons
  - nonsense-mediated decay
authors:
  - name: Elmasnur Yılmaz
    orcid: 0000-0001-9711-0203
    affiliation: 1
    corresponding: true
affiliations:
  - name: Department of Biotechnology, Ege University, İzmir, Türkiye
    index: 1
date: 31 October 2026
bibliography: paper.bib
---

# Summary

Most human multi-exon genes are cut and re-joined in more than one way before they are translated
into protein, a process called splicing. When a regulator of splicing is lost, cells can
start using splice sites that are normally silent and include stretches of intron in the
mature message. These "cryptic exons" are a hallmark of neurodegenerative diseases such as
amyotrophic lateral sclerosis (ALS) and frontotemporal dementia, where loss of the RNA
binding protein TDP-43 switches them on in genes needed by neurons
[@Ling2015; @Klim2019; @Melamed2019; @Ma2022; @Brown2022]. Two questions follow for every candidate:
is the event real rather than sequencing noise, and what does it do to the protein?

`splicescope` is a Python package that answers both from data most RNA sequencing
experiments already produce. It reads the splice-junction tables written by aligners such
as STAR [@Dobin2013], together with a reference annotation, and returns a classification
of every junction against the annotation, splice-site usage and event-level
percent-spliced-in values, a count-based test for differential splicing between two
conditions, and, given the genome sequence, a prediction of whether each event shifts the
reading frame, introduces a premature stop codon and exposes the transcript to
nonsense-mediated decay (NMD). It is distributed with a command-line interface, a Nextflow
pipeline [@DiTommaso2017], an interactive Streamlit dashboard, a simulator that writes
ground-truth datasets, and a documented validation on public data.

# Statement of need

`splicescope` is for RNA biologists and neurogeneticists who have junction tables from a
standard alignment and need an interpretable, statistically sound shortlist of cryptic
events. Three gaps motivated it.

First, novel junctions are dominated by technical noise, and separating real events often
relies on ad hoc thresholds. `splicescope` makes the categories explicit: each junction is
labelled `annotated`, `novel_donor`, `novel_acceptor`, `novel_combination` or `cryptic`
relative to the annotation, so that the classes a user filters on are defined rather than
implied.

Second, small experiments are common in this field, and some convenient tests cannot
succeed on them. A two-sided rank test on `n` replicates per group cannot return a p-value
below $2/\binom{2n}{n}$, which is 0.1 for a three-versus-three design, so nothing survives
Benjamini-Hochberg correction [@Benjamini1995] across tens of thousands of junctions,
however large the effect. `splicescope` therefore tests the read counts that splice-site
usage is computed from, with a beta-binomial likelihood-ratio test, so that evidence grows
with coverage and not only with the number of replicates.

Third, detecting an event says nothing about its consequence. Whether a cryptic exon is
tolerated, shifts the frame or triggers NMD decides whether it depletes the protein, as for
the TDP-43 target *UNC13A*, whose cryptic exon introduces a premature stop codon
[@Ma2022; @Brown2022]. Other routes exist: *STMN2* is lost through premature
polyadenylation [@Klim2019; @Melamed2019], which this prediction does not model.
`splicescope` predicts this for cassette exons and also for splice-site shifts (novel
donors and acceptors that extend or truncate an exon), using the position of the first
premature stop relative to the last exon-exon junction [@NagyMaquat1998]. The
consequence step also runs as a standalone command on candidates produced by other tools,
so it can be added to an existing analysis rather than replacing it.

# State of the field

Several established tools detect and test alternative splicing. rMATS
[@Shen2014; @Wang2024] tests differential percent-spliced-in for annotated and novel event
types from alignments using a hierarchical model across replicates. LeafCutter [@Li2018]
quantifies intron excision within clusters of overlapping introns without relying on an
annotation. MAJIQ [@VaqueroGarcia2016] models local splicing variations that go beyond
binary event types. SUPPA2 [@Trincado2018] derives percent-spliced-in from transcript
abundance estimates and tests differences quickly across conditions. These tools are
mature, and `splicescope` does not aim to replace them for discovery.

Prediction of functional consequence is less well served. IsoformSwitchAnalyzeR
[@VittingSeerup2019] predicts NMD sensitivity and protein features for switches between
full-length isoforms, which requires transcript-level quantification and isoforms that
contain the event. LeafCutter2 [@BuenAbadNajar2025] labels junctions as productive or
unproductive from annotated start and stop codons. `splicescope` addresses a narrower
question with a different input: a single junction or exon interval, observed but possibly
absent from every annotated transcript, resolved against the host transcripts it could
splice into, with the frame, the stop-codon offset and the distance to the last exon-exon
junction reported for each event alongside one of eight consequence classes.

The build-versus-contribute choice followed from this. The consequence logic needs the
genome sequence, annotated coding sequences and per-event frame arithmetic, none of which
is part of the statistical engines above, and those engines are written in different
languages around different data models (event files, intron clusters, splice graphs). A
small Python package that reads their outputs adds this layer to all of them, which is why
`splicescope consequence` accepts event tables produced elsewhere.

# Software design

The package is a set of small modules, one per step (input, junction taxonomy,
quantification, event detection, differential testing, consequence, enrichment,
classification, plotting and simulation), that exchange plain `pandas` [@McKinney2010]
data frames. Numerical work uses NumPy, SciPy and scikit-learn, and figures use Matplotlib
[@Harris2020; @Virtanen2020; @Hunter2007]. The command-line interface, the Nextflow
pipeline and the dashboard are thin layers over the same functions, so every interface
produces the same numbers.

Junction tables, rather than alignments, are the input: they are small, and a bundled
script derives the same table from any coordinate-sorted BAM file, for example from HISAT2
[@Kim2019]. The documented cost is that intron retention, which is measured from coverage
inside the intron, is out of reach.
Splice-site usage is a deliberately model-free proxy, reported beside event-level
percent-spliced-in rather than instead of it. The differential test estimates one
dispersion parameter shared across all units with a moment estimator, which is nearly
unbiased in simulation where fitting each unit on a handful of replicates is not; the
methods document also measures the price, an inflated false-positive rate for units whose
true dispersion is much larger than the shared estimate. Unobserved
junctions are filled in as measured zeros rather than treated as missing, because an
aligner omits a junction when it has no reads.

The classifier is included but held to the same standard as the statistics. It is
evaluated by cross-validation with feature scaling fitted inside each fold and is shipped
with a model card. The validation, which held whole genes out of each training fold,
reports that it ranks confidently replicating events far worse than the differential
statistic does (average precision 0.006 against 0.062) and becomes competitive only under
a permissive definition of replication, because its features describe whether a junction
is plausible rather than whether it changes between conditions. It is documented as a
noise filter, not as the way to find cryptic events.

Correctness is checked at several levels. The simulator writes a genome whose coding exons
are built from sense codons with canonical intron ends, so the consequence layer is tested
against reading frames that genuinely shift and stop codons that genuinely occur. The test
suite combines unit tests, end-to-end command-line runs and property-based tests written
with Hypothesis [@MacIver2019]. Continuous integration runs it on Python 3.10 to 3.13 and
on the oldest declared dependency versions, executes the tutorial notebook, runs the
Nextflow pipeline on its demonstration and real-data input paths and renders the dashboard
headlessly. A mutation survey deliberately breaks the code in plausible ways and reports
any change the suite fails to detect, and numbers quoted in the documentation are checked
by tests against the code that produces them.

# Research impact statement

`splicescope` is early-stage software. Its evidence of significance is a reproducible
validation on public data rather than adoption by other groups. The figures below were
measured with versions up to 0.8.1; `validation/README.md` lists the later corrections
that will move them when re-run.

On a TDP-43 knockdown in human cells, GEO GSE245332 [@Hou2024], with three replicates per
group aligned with HISAT2 and analysed with the GENCODE annotation [@Frankish2019], a
rank test found no significant junctions among 138,643 tested, as the bound above
predicts, whereas the count-based test recovered known TDP-43 cryptic splicing targets
including *STMN2*, *HDGFL2*, *PFKP*, *ARHGAP32*, *KALRN*, *AGRN*, *ATG4B* and *RAP1GAP*.
The whole pipeline ran on the six samples in about 1.5 minutes on a laptop. Replication in
a second experiment, GSE122069 [@Melamed2019], provided the non-circular label used to
evaluate the classifier.

The NMD prediction was compared with external measurements from two independent
experiments in which NMD was blocked in TDP-43-depleted neurons [@Sinha2025; @Zeng2025].
Among 356 cassette exons with a measured change, the 275 predicted to carry an
NMD-triggering stop codon gained more on NMD blockade than the 68 where no stop was
possible (area under the ROC curve, AUROC, 0.675; 95% bootstrap interval 0.601 to 0.752).
On the 348 of these exons that LeafCutter2 [@BuenAbadNajar2025] also labels, the
`splicescope` NMD call separated NMD sensitivity about as well as LeafCutter2's
unproductive-versus-productive label (AUROC 0.684 and 0.697; the latter leaves out the 87
exons LeafCutter2 places in neither class), and the two agreed on 85% of the 261 exons it
does place. For 3,419 interpretable splice-site shifts, predicted NMD targets were
NMD-sensitive more often than events where no stop was possible (48.5% against 23.5%, odds
ratio 3.07). The same analyses show the limit of the call: sensitivity is high and
specificity low, because a random intronic interval read in a fixed frame usually contains
a stop codon. The documentation therefore recommends using the prediction to interpret
events that are already significant, not to rank candidates. The per-event tables and the
script that regenerates every one of these numbers are in the repository.

<!-- Before submission: add concrete evidence of research use (a thesis chapter, preprint
or analysis that used splicescope, and any users outside the author's group), and state
whether a related manuscript on the NMD validation is published, in review or in
preparation. JOSS screens for this; see the eligibility notes in pull request #1. -->

# AI usage disclosure

Generative AI (Claude Code, Anthropic) was used to assist with software development,
testing, documentation and the preparation of this paper. The author reviewed and edited
all AI-assisted output. The correctness of that output was checked by the automated test
suite and continuous integration described above, by regenerating the quoted
NMD-validation numbers from the committed tables, and by checking every reference against
a bibliographic record. The author takes full responsibility for the software and the
paper.

# Acknowledgements

<!-- Before submission: acknowledge colleagues and funding sources, say whether any funder
was involved, and declare any conflicts of interest or state that there are none. -->

# References
