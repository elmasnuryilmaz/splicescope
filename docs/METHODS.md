# Methods

This document specifies exactly what `splicescope` computes, so results are
interpretable and reproducible. Notation: a junction *j* is an intron with 1-based
inclusive coordinates `[s, e]` on chromosome `c` and strand `σ`.

## 1. Junction ingestion and coordinate convention

Junctions are read from STAR `SJ.out.tab`. Following STAR, `s` is the first intronic
base and `e` the last intronic base. The 5′ splice site (donor) `D(j)` and 3′ splice
site (acceptor) `A(j)` are strand-dependent:

$$D(j),\,A(j) = \begin{cases}(s,\,e) & \sigma = +\\ (e,\,s) & \sigma = -\end{cases}$$

Only uniquely-mapping reads (`n_unique`) are used as the junction count, to avoid
multi-mapping inflation.

## 2. Known junctions from annotation

Exons are grouped per transcript from a GTF and sorted by coordinate. Each gap between
consecutive exons `(…​,e_k)` and `(s_{k+1},…)` defines a known intron
`[e_k + 1,\; s_{k+1} - 1]`. The union over transcripts gives the known-junction set
`𝒦`, plus the known-donor set `𝒟` and known-acceptor set `𝒜`.

## 3. Junction classification

Each observed junction is assigned exactly one class:

| class | rule |
|-------|------|
| `annotated` | `j ∈ 𝒦` |
| `novel_combination` | `D(j) ∈ 𝒟` **and** `A(j) ∈ 𝒜`, but `j ∉ 𝒦` (e.g. exon skipping) |
| `novel_acceptor` | `D(j) ∈ 𝒟`, `A(j) ∉ 𝒜` |
| `novel_donor` | `A(j) ∈ 𝒜`, `D(j) ∉ 𝒟` |
| `cryptic` | `D(j) ∉ 𝒟` **and** `A(j) ∉ 𝒜` |

Genes are assigned by shared splice site. Note that a **cryptic exon** manifests as a
`novel_acceptor` junction sharing the upstream known donor *plus* a `novel_donor`
junction sharing the downstream known acceptor — both are detectable here.

## 4. Splice-site usage (Ψ)

We use a deliberately model-free quantity. For a junction *j* with donor `D`,

$$\Psi_D(j) = \frac{\mathrm{count}(j)}{\sum_{j' :\, D(j') = D} \mathrm{count}(j')}$$

i.e. the fraction of that donor's spliced reads that flow through *j*. Acceptor usage
`Ψ_A` is symmetric. To avoid noisy ratios from thinly-covered sites, if the denominator
is `< min_reads` (default 10) in a sample, Ψ is set to `NaN` (uninformative). Ψ needs no
event model, which is what makes it robust to incomplete annotation.

## 5. Differential splicing

For two conditions A and B and each junction, we compare the per-sample Ψ vectors:

- **Effect size:** `ΔΨ = mean(Ψ | B) − mean(Ψ | A)`.
- **Test:** two-sided **Mann–Whitney U** on the per-sample Ψ values. Ψ is a bounded
  ratio and replicate counts are small, so a rank test is preferred over a *t*-test.
- **Multiple testing:** **Benjamini–Hochberg** FDR across all tested junctions, with the
  standard monotonicity enforcement.

Junctions with fewer than `min_samples` informative replicates per group are skipped.
"Significant" defaults to `q ≤ 0.05` and `|ΔΨ| ≥ 0.1`.

## 6. Feature engineering for cryptic calling

Each novel junction becomes a feature vector:

| feature | rationale |
|---------|-----------|
| `intron_length` | genuine introns have a plausible length distribution |
| `log_max_count` | real events are reproducibly supported |
| `n_samples_support` | real events recur across replicates; noise does not |
| `mean_psi_donor` | real events take a non-trivial share of their donor |
| `canonical_motif` | genuine splicing favours GT/AG (or CT/AC) |
| `dist_known_donor` | cryptic sites sit near, but not on, annotated sites |
| `dist_known_acceptor` | as above, 3′ side |
| `is_novel_both` | both sites unannotated (strongest cryptic prior) |

Distances to the nearest annotated site are computed by binary search over the sorted
per-`(chrom,strand)` known-site positions.

## 7. Classifier and evaluation

A `StandardScaler` → `RandomForestClassifier` pipeline (`class_weight="balanced"` to
handle the noise-dominated negative class). Performance is estimated with **stratified
k-fold cross-validation** using out-of-fold probabilities, reporting **ROC-AUC** and
**average precision** (PR-AUC). Feature attributions use **permutation importance** on
the fitted model. Every fitted model exports a **model card** (JSON + Markdown) recording
hyper-parameters, the evaluation protocol, metrics, importances, intended use and
limitations. Scaling is fit inside each CV fold via the pipeline, so there is no
train/test leakage.

## 8. Simulation model

To keep the toolkit self-contained and testable, `simulate` generates a small annotated
genome and per-sample junctions:

1. **Canonical introns** — Poisson-distributed counts (λ ≈ 200) in every sample.
2. **True cryptic-exon events** in a fraction of genes: a cryptic exon inside an intron
   yields a `novel_acceptor` junction (shares the known donor) and a `novel_donor`
   junction (shares the known acceptor), with inclusion up-regulated in condition B
   (variable effect size). A minority use a non-canonical motif.
3. **Noise** — sporadic, low-support novel junctions; a fraction are *hard negatives*
   that share a known donor, use a canonical motif and recur (so the classes overlap).
4. **Label noise** (optional) — a fraction of ground-truth labels are flipped to mimic
   imperfect curation, so the ML task is realistically hard rather than trivially
   separable.

Generation is fully seeded and therefore reproducible.

## 9. Limitations

- Ψ here is splice-site *usage*, not event-level PSI (cassette exon, A5SS/A3SS, IR); the
  latter is a natural extension built on the same junction graph.
- The bundled data is simulated. On real data the classifier must be retrained on curated
  labels and validated on **held-out genes** (not just held-out junctions) to avoid
  optimistic estimates from shared-gene leakage.
- Classification uses splice-site identity, not sequence models of the splice site; a
  learned donor/acceptor motif score would refine the `canonical_motif` feature.
