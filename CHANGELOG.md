# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/) and this project adheres to
[Semantic Versioning](https://semver.org/).

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
