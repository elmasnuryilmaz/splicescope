# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/) and this project adheres to
[Semantic Versioning](https://semver.org/).

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
