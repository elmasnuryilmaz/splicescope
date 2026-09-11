# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/) and this project adheres to
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Fixed
- **`normalize_gene_id` collapsed whole gene families in a model organism.** 0.9.0 began
  stripping a trailing `.<digits>` from every identifier, on the docstring's own premise
  that "gene symbols do not end in `.<digits>`". *C. elegans* sequence names — the standard
  identifier for the majority of worm genes, which have no CGC name — are exactly that
  shape, so `C42D8.1`, `C42D8.2` and `C42D8.3` normalised to one key. The background size,
  the set size, the overlap, the fold enrichment, the p-value and the reported
  leading-edge gene list all shrank with them, silently and with exit 0. The suffix is now
  removed only when the stem looks like an Ensembl or RefSeq accession.
- **`run` crashed on digit-only sample names.** pandas reads a `sample` column of `101`,
  `102`, … as int64, which can never equal the strings derived from filenames, so the new
  mismatch check took its error branch on a perfectly consistent cohort — and then died
  with `TypeError: sequence item 0: expected str instance, int found` and exit 1 rather
  than printing its diagnostic. The column is now read as text, which also keeps `007`
  from becoming `7`.
- **MXE detection could grow as the fourth power of a locus's exons.** Removing the
  one-pair-per-anchor truncation in 0.9.0 fixed a correctness bug and introduced a
  scalability one: every `(upstream, downstream)` junction combination at an anchor is a
  candidate, so `n` real exons arrive with about `n**2` of them — the real exons plus
  every span from one exon's start to a later exon's end — and pairing all of those is
  quadratic again. **80 exons under one anchor produced 1,677,140 rows in 5.9 s.**
  Candidates are now ranked by the read support of their weaker flanking junction and the
  best `max_candidates` (default 50) are kept, with a warning naming the busiest anchor.
  Reachable as `detect_events(max_candidates=...)` and `run --max-mxe-candidates`.
  Two geometric rules were tried first and both **deleted 15–22 % of the injected events**
  along with the spans, because a noise junction sharing the anchor's donor is
  indistinguishable in shape from a real exon; read support distinguishes them. All 200
  injected events are recovered at every seed tested.
- **An unestimable dispersion was silently treated as no dispersion.** When no unit has
  more informative replicates than fitted group means there is no residual degree of
  freedom anywhere, and `estimate_precision` returned `max_precision` — asserting *no
  overdispersion* rather than admitting the design says nothing. The beta-binomial test
  then narrows to a plain binomial one: a **1-vs-1 comparison of Ψ 0.300 against 0.360 at
  1000 reads returns q = 4.5e-03**, on no replication whatsoever, with nothing in the
  output to say so. `dispersion_is_estimable` now reports whether the design can support
  an estimate and `differential_splicing` warns when it cannot. The default
  `min_samples=2` already avoided this; passing `min_samples=1` walked straight into it.

## [0.9.0] — 2026-09-11

### Added
- **The protein-consequence layer now runs on the built-in data, with no downloads.**
  `simulate` writes the chromosome its annotation describes: coding exons drawn from the
  61 sense codons so every gene carries a real open reading frame ending in a single stop,
  canonical `GT`/`AG` intron ends, CDS records in the GTF, and a `.fai` index. The 5′UTR is
  staggered by one base per gene so the coding frame at an intron boundary cycles through
  0, 1 and 2 — with the uniform 120 bp exons a fixed UTR put every intron on a codon
  boundary, and frame inheritance, the subtle half of the prediction, was never exercised.
  `simulate_genome` uses its own random stream, so no junction-level number moves.
- **Splice-site shifts are reachable from the command line.** `run --genome` now writes
  `junction_consequence.tsv` beside `consequence.tsv`, and `consequence` gained
  `--mode {auto,exon,junction}`, defaulting to reading the mode off the columns present.
  `annotate_junction_consequences` — the headline of 0.8.0, and already described in the
  README — existed only as a Python API and could not be invoked by a user.
- A junction already explained by a detected cassette or MXE event is excluded from the
  splice-site-shift pass. Read on its own, a cassette inclusion junction looks like an exon
  extension running to the end of the intron, which is the wrong reading of it; on the demo
  this removes 12 of 31 candidate junctions.
- `plot_consequence_summary`, and the `event_summary` / `event_volcano` figures that
  `plotting` had provided since 0.4.0 but `run` never called.
- 14 tests for the above, including that the simulated ORFs translate cleanly, that the
  `.fai` offsets address the FASTA correctly across line boundaries, that all three reading
  frames occur, and that cassette junctions are not double-reported as shifts.
- CI smoke-tests the consequence path end to end.

### Fixed
- **`detect_mxe_events` emitted one event per anchor pair and abandoned the rest.** Inside
  each shared (donor, acceptor) context the first non-overlapping exon pair set a flag and
  broke out of both loops, so a locus with three or more mutually exclusive exons — the
  interesting case — reported one pair and lost the others. The `seen` set already handled
  deduplication, so the breaks bought nothing. The loss was total rather than partial: the
  discarded exons' junctions share the emitted pair's donor and acceptor, so once
  `detect_events` marked those as used, each site had a single remaining alternative and
  A5SS/A3SS detection rejected them too. The exons disappeared from `events.tsv` with no
  warning. Worse than under-reporting, the survivor was the genomically *leftmost* pair
  rather than the best-supported one, so noise junctions displaced real events: on
  `simulate_dataset(n_genes=200, mxe_fraction=1.0)`, **31–44 of 200 emitted MXE rows
  carried wrong exon coordinates** while the injected pair sat unused in the same group
  (seeds 4, 6, 11). Every pair is now enumerated; all 200 injected events are recovered at
  every seed, and detection still takes 0.02 s.
- `max_exon` is reachable through `detect_events` and as `run --max-exon`; the 1 kb
  candidate-exon window was hard-coded for anyone using the documented entry points.
- **The `min_reads` coverage filter was inert for the beta-binomial test** — the default
  test, on every run since 0.7.0. `_wide` pivoted Ψ with `aggfunc="sum"`, and pandas sums
  an all-NaN group to `0.0`, so `psi_wide.notna()` was true for every cell that had a row
  at all and the NaN masking `compute_psi` performs was discarded. `--min-reads` changed
  nothing: on `simulate_dataset(n_genes=20, n_per_group=4, seed=3)` the test returned the
  same 105 units and 18 significant hits at `min_reads=10` and at `min_reads=100000`,
  where every Ψ is NaN. Consequences were not cosmetic — the BH denominator was inflated
  by every sub-threshold unit, so **every q-value was wrong**, and `estimate_precision`
  was fitting dispersion to 1–2 read residuals. Ψ is now pivoted with `mean`; counts are
  still summed. The rank-sum path was never affected because it uses `dropna`, which is
  why `validation/README.md` records 138,643 rank tests against 189,700 count tests on
  identical input: the 51,057 difference *was* the bug.
- **A junction absent from a sample is a measured zero, not missing data.** Aligners list
  only junctions with at least one read, so `compute_psi` never produced a row for a
  junction in a sample where it had none, the pivot left NaN, and the unit was dropped for
  want of informative replicates. The events layer already defaulted a missing count to
  zero, so the two halves of the tool disagreed about the same data. The casualties were
  exactly the cleanest cryptic events: a junction with 0 reads in all three controls and
  100 in all three knockdowns was **absent from `differential_splicing.tsv` entirely**,
  while leakier events were reported. `compute_psi` now completes the junction × sample
  grid (`add_unobserved_zeros`), adding a zero only where that junction's donor or
  acceptor site has coverage in the sample, so a site nobody sequenced stays unobserved.
  `fill_unobserved=False` restores the old behaviour.
- 8 regression tests; 5 of them fail against the previous code.

### Changed
- Every junction- and event-level p-value and q-value moves as a result of the two fixes
  above, in both directions: sub-threshold units leave the BH family, and junctions absent
  from a group enter it. The real-data numbers in `validation/README.md` were measured
  with 0.8.1 and are marked as awaiting a re-run; the qualitative conclusions there do not
  depend on the exact values. On simulated data the cryptic classifier moves from
  ROC-AUC 0.785 to 0.795 (AP 0.702 → 0.700), because `mean_psi_donor` now averages over
  the samples where a junction was measured absent as well.

### Also added
- Test coverage for the plotting module: threshold behaviour of both volcano
  plots, fixed event-type ordering with zero-fill, top-N truncation in the
  enrichment plot, and that `savefig` creates missing parent directories.
- `CODE_OF_CONDUCT.md`, issue templates (bug report, feature request) and a pull
  request template that requires stating the effect on published numbers.
- Python 3.13 to the CI matrix, and the corresponding trove classifiers.
- `Issues`, `Changelog`, `Documentation` and `Archive` (Zenodo DOI) entries under
  `[project.urls]`.
- Zenodo DOI badge, and the DOI in `CITATION.cff` (concept and version).
- Gzipped input: `read_gtf_junctions` and `read_gmt` now open `.gz` directly. GENCODE and
  Ensembl ship `.gtf.gz`, and the README's own example passes one, which until now failed
  with a `UnicodeDecodeError`.
- `io.sample_name_from_path`, and a test that `__version__`, `CITATION.cff` and the
  CHANGELOG agree — they had drifted to 0.5.0, 0.8.1 and 0.8.1 respectively, so
  `splicescope --version` reported a release three versions old.

### Changed (breaking)
- `events.cassette_psi` is removed; use `events.event_psi`, which keys by `event_id`.
  `cassette_psi` keyed its rows by the *skipping* junction, so several cassette exons
  between the same pair of flanking exons collapsed into one test — and since they can
  move in opposite directions, their ΔΨ cancelled. A worked example with three events at
  one locus produced a single junction key. Nothing in the pipeline used it; `run` has
  been on `event_psi` since 0.3.0.

### Added (statistics)
- `betabinom.estimate_precision_per_unit` and
  `differential_splicing(dispersion="per_unit_floor")`, which takes the smaller of each
  unit's own precision and the shared one. One shared dispersion is **not uniformly
  conservative**, contrary to what the docstring claimed: on a null with half the units at
  `s = 200` and half at `s = 5`, the shared estimate is 10.7 and the loosely dispersed half
  runs at a false-positive rate of **0.173** against a nominal 0.05 while the tight half
  runs at 0.000 — the pooled 0.087 hides both. The floor brings the loose half to 0.090 and
  the overall rate to 0.045, at a cost of roughly a sixth of the power (0.52 → 0.44 on a
  homogeneous Ψ 0.30-vs-0.45 contrast). **The default is unchanged**, so no published
  number moves; `docs/METHODS.md` §5.4 records the trade-off and the tests record the
  behaviour.
- More broadly: `predict_consequence` gained `contiguous_downstream`; `Transcript` gained
  `cds_phase` and `frame_at`; `enrich.normalize_gene_id` is public.

### Fixed (real data)
- **Strand-undefined junctions became phantom cryptic calls.** STAR writes strand code 0
  whenever the intron motif does not reveal a strand — routine for non-canonical junctions
  — and such a junction can never equal a stranded annotation. It was therefore classified
  `cryptic` however ordinary it was, lost its `gene_id`, and, being the only junction at
  its own `(chrom, position, ".")` site, was handed **Ψ ≡ 1.0 in every sample**: a row in
  the differential table that can never show a difference, inflating the BH denominator and
  the headline `cryptic` count alike. `resolve_unstranded` now places such a junction on
  the strand its annotation implies, by exact intron or by either splice site. Junctions
  the annotation cannot place keep `"."` — the strand really is unknown.
- **Alternative-splice-site events were deleted for sharing an anchor with a cassette.**
  Cassette and MXE junctions are excluded from A5SS/A3SS detection so the same signal is
  not told twice, but they were dropped before the alternatives at a site were counted, so
  a site with a genuine third alternative fell to one and its event vanished. A site is now
  skipped only when *every* junction at it is already explained.
- `splicescope consequence` raised `KeyError: 'consequence_class'` when filtering an event
  table to cassette exons left nothing; it now writes an empty result and exits 0.

### Fixed (statistics and biology)
- **Samples in neither group inflated the likelihood-ratio test without bound.** `lrt`
  fitted the null over every valid sample but the alternative over the two groups only,
  so the two hypotheses were not fitted to the same observations and the whole likelihood
  of the ungrouped samples was charged to the null. Adding six samples that belong to
  neither group moved an unchanged 3-vs-3 comparison from **p = 1.6e-09 to p = 9.7e-93**,
  with `n_a` and `n_b` still reported as 3 and 3. Samples outside both groups now take no
  part. This is not exotic: an SJ directory routinely holds more samples than a given
  contrast lists, which `run` now also warns about.
- **An in-frame exon truncation was reported as a premature stop.** Removing a whole
  number of codons leaves the downstream reading frame untouched, so the first in-frame
  stop found is the transcript's *own* — reported as `ptc_escape` with a PTC offset. On the
  demo genome, deletions of 3, 6, 9 and 30 nt all returned the identical offset 296, which
  is the annotated stop codon. A stop at or beyond where the protein natively ends is no
  longer called premature; these are `exon_truncation`, and frameshifting truncations still
  report their genuinely premature stops.
- **GTF phase was discarded, so every 5'-incomplete CDS was translated in the wrong
  frame.** GENCODE marks these `cds_start_NF` and gives the first CDS record a non-zero
  phase; `load_transcripts` ignored column 8 entirely. `Transcript.cds_phase` now carries
  it — taken from the first block in *transcription* order, so the minus strand is right —
  and `Transcript.frame_at` applies it wherever frame is inherited.
- **A PTC in an extension of the final exon was called `ptc_nmd`.** An extension is
  contiguous with the exon it joins, so there is no junction between them; when that exon
  is the last one, the stop lies past the final exon-exon junction and NMD cannot be
  triggered. The old code measured the distance to a junction that does not exist — on the
  demo genome it reported `ptc_nmd` at a distance of 120 for an event that escapes.
- **An exon-skipping junction was read as a splice-site shift.** Both its sites are
  annotated, just not as a pair, so `junction_change` matched one of them and reported the
  skipped exon *and the intron beyond it* as a single contiguous deletion. Such junctions
  are now left uninterpreted, which is what they are: a shift has exactly one novel site.

### Fixed (Nextflow)
- **The pipeline's real-data mode could not have worked.** One channel was passed for all
  three inputs, and the process script then looked for `<sj_dir>/sj`,
  `<sj_dir>/annotation.gtf` and `<sj_dir>/groups.tsv` — a layout no real dataset has. The
  three inputs are now staged independently, `--genome` is exposed (with its `.fai` staged
  alongside), `--min-reads` is passed through, and demo mode runs the consequence layer
  too now that `simulate` writes a genome. *Not executed locally: Nextflow is not
  installed on the development machine. The CLI invocation the process builds was
  verified against a replica of Nextflow's staging layout.*

### Fixed (enrichment)
- **Pathway ORA could never test a single gene set against a real annotation.**
  `over_representation` intersected identifiers verbatim, but a GTF gives versioned
  accessions (`ENSG00000141510.16`) while every real GMT gives symbols or *unversioned*
  accessions. No set ever met `min_size`, so the result was always empty — and the CLI
  reported `0 tested` rather than an error, so the run looked successful. A four-gene set
  overlapping the background perfectly returned nothing, lost entirely to the `.16`.
  Identifiers now go through `normalize_gene_id` (version suffix dropped, case-folded).
- Gene symbols reach the differential table: `read_gtf_junctions` returns `gene_name`,
  `annotate_junctions` carries it, and `differential_splicing` keeps it. Symbol-keyed
  gene sets — MSigDB `*.symbols.gmt`, GO, KEGG, i.e. most of them — had nothing to match
  against before. `enrich_differential` uses whichever of `gene_id`/`gene_name` overlaps
  the sets more.
- `run --gene-sets` now warns, with an example identifier from each side, when no gene set
  shares an identifier with the tested genes.

### Fixed (packaging and CLI)
- **`run` reported "0 significant junctions" and exited 0 when sample names did not
  match.** Names were derived with `stem.replace(".SJ.out", "")`, but STAR writes
  `{prefix}SJ.out.tab`, so the near-universal real filename `SampleA_SJ.out.tab` yielded
  the sample name `SampleA_SJ.out` while `groups.tsv` said `SampleA`. Nothing checked the
  overlap: conditions mapped to NaN, no sample joined either group, the test returned an
  empty frame, and the user got a clean run, an empty table and the conclusion "no
  differential splicing". Names now tolerate any separator before `SJ.out`, a total
  mismatch is an error naming both sets, and a partial one warns.
- The single source of truth for the version is `splicescope.__version__`; `pyproject.toml`
  reads it dynamically.
- `savefig` now closes the figure it wrote. pyplot keeps every figure alive until closed,
  and `run` writes up to eight per invocation, so a run leaked all of them and matplotlib
  warned past twenty. `close=False` keeps the old behaviour.

## [0.8.1] — 2026-09-03

### Fixed
- Streamlit dashboard deployment: the app now imports from `src`, declares its runtime
  dependencies explicitly, guards against load errors, and drops the removed
  `use_container_width` argument.

## [0.8.0] — 2026-08-30

### Added
- **Splice-site shifts are now interpreted**, not just cassette exons. A novel donor
  or acceptor paired with an annotated site on the other side is resolved against the
  annotation into the sequence it adds to, or removes from, the neighbouring exon
  (`junction_change`, `junction_hosts`, `predict_junction_consequence`,
  `annotate_junction_consequences`). These are the majority of cryptic events reported
  by junction-level callers: on a published set of 3,548 such junctions, interpretable
  coverage went from **0% to 96%** (98% of novel donors, 99% of novel acceptors).
- **Premature stops downstream of the event are found.** Previously an insert that
  shifted the frame without carrying a stop was reported as `frameshift` with no PTC,
  and a truncation got no PTC search at all. Both now assemble the retained downstream
  exons and continue the in-frame scan (`downstream_sequence`), so the 50-nucleotide
  rule can be applied where the stop actually is.
- New class `exon_truncation` for a splice-site shift that shortens an exon without
  producing a stop.

### Changed
- Against experimentally measured NMD sensitivity the class ordering is now monotonic
  and the effect is stronger. On cassette exons (i3Neuron, quantitative label) the
  ambiguous `frameshift` bucket disappears and classes rank
  `ptc_nmd` 12.7 > `ptc_escape` 8.6 > in-frame 6.3 > UTR 4.2 > non-coding −1.6 ΔPSI.
  On splice-site shifts (SMG1i, binary label) `ptc_nmd` is 48.5% NMD-sensitive against
  23.5% for events where no stop is possible — **odds ratio 2.31 → 3.07,
  two-sided Fisher p 1e-04 → 4.5e-12**.

## [0.7.0] — 2026-08-30

### Changed
- **The differential test now models read counts** (`betabinom`), replacing the
  Mann–Whitney U test on per-sample Ψ as the default wherever counts are available.
  A rank test on `n` replicates per group cannot return a p-value below `2/C(2n,n)`
  — 0.1 at 3-vs-3 — so nothing could ever survive genome-wide BH correction. On a real
  human TDP-43 knockdown (GSE245332, 3 vs 3) the old test returned **0 significant
  junctions with min(q) = 1.0**, while rMATS called 3 724 events on the same BAMs, and
  46 junctions had |ΔΨ| ≥ 0.5. The signal was there; the statistic could not reach it.
- `differential_splicing` gained `test` (`"auto"`, `"betabinom"`, `"ranksum"`) plus
  `inc_col`/`total_col`. `"auto"` picks the count-based test when counts are present and
  falls back to the rank test otherwise, so existing calls keep working.
- `compute_psi` now also returns `donor_total` and `acceptor_total`, the denominators Ψ
  was formed from, so junction-level tests can use counts too.

### Added
- `betabinom`: vectorised beta-binomial fitting (`fit_mu`, `loglik`), a df-corrected
  moment estimator for the shared dispersion (`estimate_precision`) and the
  likelihood-ratio test (`lrt`). Fits are bisections on a monotone score, so hundreds of
  thousands of events cost seconds.
- `min_achievable_rank_pvalue` so the limitation above can be shown rather than asserted.
- 11 tests: Ψ recovery, boundary events, dispersion bias, null calibration
  (nominal 0.05 → 0.045–0.053) and power past the rank-test floor.

### Fixed
- `detect_mxe_events` compared every junction against every other one on a chromosome to
  find candidate exons, which is quadratic and dominated the whole pipeline (81% of
  runtime on a real dataset). A candidate exon is at most `max_exon` long, so only
  downstream junctions starting inside that window can pair with a given upstream
  junction; sorting by start and binary-searching the window is exact and much cheaper.
  **89.5 s → 0.5 s** on two chromosomes, with identical output. Together with the
  consequence fix below this takes a full six-sample run from ~45 minutes to **~1.5
  minutes**; only the row order of `events.tsv` changes.
- `annotate_consequences` scanned every transcript for every event, which is quadratic:
  296,034 cassette exons took ~34 minutes. Transcripts are now indexed by gene
  (`index_by_gene`) and only the event's own gene is searched — the same run takes
  **19 seconds**, and events can no longer be attributed to an overlapping neighbour.

### Documentation
- `validation/` records an end-to-end run on real data (GSE245332) plus a replication
  experiment against a second, independent dataset (GSE122069): `bam2sj.awk` and
  `extract_sj.sh` derive STAR-format junctions from any sorted BAM, and `add_motifs.py`
  fills in intron motifs from an indexed genome (85–90% canonical, and the motif agrees
  with the aligner's `XS:A` tag for 99.94% of junctions).
- The classifier's positioning in the README and METHODS is corrected. Scored against
  replication in the independent dataset it reaches an average precision of 0.006 for
  confidently replicating junctions where the beta-binomial statistic reaches 0.062, and
  the combination is worse than the statistic alone. Its features describe whether a
  junction is plausible, not whether it changes between conditions; it is a noise filter,
  not the way to find cryptic events.

### Notes
- Dispersion is estimated about each group's own mean. Pooling groups that genuinely
  differ charges the difference to dispersion — at ΔΨ = 0.35 the pooled estimate collapses
  from `s ≈ 200` to `s ≈ 5` and power goes to zero.

## [0.6.0] — 2026-08-30

### Added
- **Protein-consequence prediction** (`consequence`): for each cassette exon, inherit the
  reading frame from the coding sequence upstream of the host intron, translate the
  insert, locate the first in-frame premature termination codon (PTC) and apply the
  50-nucleotide rule to decide whether it triggers nonsense-mediated decay. Events are
  classified as `ptc_nmd`, `ptc_escape`, `frameshift`, `in_frame_insertion`,
  `utr_insertion`, `non_coding_host` or `no_host_transcript`. Detecting a cryptic exon
  says nothing about whether it matters; this is what closes that gap.
- `GenomeFasta`: random access to an indexed genome through its `.fai` alone, so
  sequence resolution adds no compiled dependency.
- `load_transcripts`: exon/CDS models from a plain or gzipped GTF, optionally restricted
  to a set of gene names or Ensembl IDs (versioned or not).
- CLI `consequence --events --gtf --genome --out` for candidates produced elsewhere
  (e.g. an existing rMATS run), and `run --genome` to fold the same prediction into the
  end-to-end pipeline as `consequence.tsv`.
- 11 tests covering frame inheritance, PTC detection, the NMD distance rule, strand
  handling and the FASTA index reader.

### Notes
- Use the **full** GENCODE annotation rather than `basic`: on 300 real cryptic exons the
  reduced set left 43% of events without a host transcript, against 20% with the full
  annotation.

## [0.5.0] — 2026-08-30

### Added
- **Pathway over-representation analysis** (`enrich`): a hypergeometric/BH-FDR test for
  gene sets enriched among differentially spliced genes — the same statistic as
  clusterProfiler, on any `{term: [genes]}` mapping. No gene sets are bundled (real
  enrichment needs real annotations).
- GMT gene-set reader (`io.read_gmt`, MSigDB/GO/KEGG-compatible) and an enrichment plot.
- CLI `run --gene-sets sets.gmt` writes `enrichment.tsv` and a figure.
- An **executed tutorial notebook** (`examples/tutorial.ipynb`, rebuilt via
  `examples/_build_tutorial.py`) with a Colab badge.
- Tests for the ORA math, GMT parsing and end-to-end enrichment.

## [0.4.0] — 2026-08-29

### Added
- **Mutually-exclusive-exon (MXE) events** — two non-overlapping exons between shared
  flanking exons, with PSI = exon-A share. splicescope now covers 4 of the 5 canonical
  rMATS event classes (SE, MXE, A5SS, A3SS; RI remains out of scope — needs coverage).
- MXE junctions are excluded from downstream A5SS/A3SS detection; the simulator can inject
  MXE events (`mxe_fraction`) and suppresses the skipping isoform for those introns.
- CLI `simulate --mxe`; the showcase gains MXE in the event-type and event-volcano panels.
- Tests for MXE detection and event-level differential inclusion.

## [0.3.0] — 2026-08-29

### Added
- **Alternative 5′/3′ splice-site events** (A5SS, A3SS) alongside cassette exons, with a
  unified `detect_events` / `event_psi` API and rMATS-style PSI per type.
- Cassette-exon junctions are excluded from alt-splice-site detection so shared signal is
  not double-reported.
- `differential_splicing` gained a `key` argument for event-level tests (`key=["event_id"]`).
- Simulator can inject A5SS/A3SS events (`alt_ss_fraction`); CLI `run` writes `events.tsv`,
  `event_differential.tsv` and an event volcano; the showcase figure is now a 2×3 panel.
- Tests for A5SS/A3SS detection and event-level differential inclusion.

### Notes
- Intron retention (RI) is intentionally out of scope: it cannot be quantified from
  junctions alone (it needs intronic read coverage). Documented in METHODS.

## [0.2.0] — 2026-08-29

### Added
- **Event-level splicing** (`events`): reconstruct cassette (skipped) exons directly from
  junctions and compute rMATS-style percent-spliced-in (Ψ). Strand-independent detection.
- `differential_splicing` now works on any Ψ column (e.g. `psi_cassette`), so cassette
  events get ΔΨ + Mann–Whitney + FDR for free.
- CLI `run` writes `cassette_events.tsv`, `cassette_differential.tsv` and a
  `cassette_volcano.png`.
- Tests for cassette detection, the PSI formula, and end-to-end differential inclusion.

## [0.1.0] — 2026-08-27

### Added
- Junction I/O for STAR `SJ.out.tab` and known-intron extraction from GTF (`io`).
- Junction taxonomy: annotated / novel_donor / novel_acceptor / novel_combination /
  cryptic (`annotate`).
- Model-free splice-site usage Ψ (`quantify`).
- Differential splicing: ΔΨ, Mann–Whitney U, Benjamini–Hochberg FDR (`diff`).
- Feature engineering for cryptic-event calling (`cryptic`).
- Cross-validated `CrypticClassifier` with permutation importance and a model card (`ml`).
- Biologically faithful synthetic data generator with ground-truth labels and optional
  label noise (`simulate`).
- Publication-quality, headless-safe plots (`plotting`).
- `splicescope` CLI (`simulate`, `run`).
- Nextflow (DSL2) pipeline and a Streamlit dashboard.
- Test suite (unit + property + end-to-end) and GitHub Actions CI across Python 3.10–3.12.
