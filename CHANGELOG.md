# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/) and this project adheres to
[Semantic Versioning](https://semver.org/).

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
  p 1e-04 → 2.3e-12**.

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
