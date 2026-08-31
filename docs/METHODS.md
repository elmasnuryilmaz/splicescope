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

For two conditions A and B and each junction or event we test the **read counts** Ψ was
formed from, not Ψ itself.

### 5.1 Why not a rank test

A two-sided Mann–Whitney U on `n` replicates per group cannot return a p-value below

```
p_min = 2 / C(2n, n)
```

which is `0.1` for the 3-vs-3 designs that dominate RNA-seq (`0.029` at 4-vs-4, `1.1e-5`
even at 10-vs-10). Surviving Benjamini–Hochberg across ~1.4 × 10⁵ junctions needs the top
hit at `p ≤ 3.6e-7`, so a rank test reports **nothing significant at any realistic
replicate count**, however large the effect. On a real TDP-43 knockdown (3 vs 3) this
produced 0 hits with `min(q) = 1.0`, while rMATS called 3 724 events on the same BAMs.
Modelling counts removes the floor: evidence grows with coverage, not only with replicates.

### 5.2 Beta-binomial likelihood-ratio test

Inclusion reads `k` out of `n` informative reads are beta-binomial with mean Ψ and
precision `s = α + β`, so `α = Ψs` and `β = (1 − Ψ)s`. The beta layer absorbs the
replicate-to-replicate variability a plain binomial would mistake for signal. Per unit we
fit Ψ by maximum likelihood under

- **H₀** — one Ψ shared by every sample, and
- **H₁** — a separate Ψ per group,

and refer `G = 2(ℓ₁ − ℓ₀)` to a χ² with one degree of freedom. Ψ is fitted by bisection on
the score function, which is monotone in Ψ, so the whole step vectorises over events.

### 5.3 Dispersion

`s` is shared across events and estimated once from **df-corrected Pearson residuals**:
under the model `E[r²] = 1 + (n − 1)/(s + 1)`, so one pass over the residual sum
determines `s`. Two choices matter:

- **Not profile likelihood.** Fitting one Ψ per event on a handful of replicates absorbs
  part of the variance being measured, biasing `s` upward (200 → 631 in simulation) and
  making the test anti-conservative. The moment estimator corrects for the fitted means
  through the residual degrees of freedom and is essentially unbiased (200 → 202).
- **Residuals about each group's own mean.** Pooling groups that genuinely differ charges
  that difference to dispersion: at ΔΨ = 0.35 the pooled estimate collapses from `s ≈ 200`
  to `s ≈ 5` and power drops to zero. `estimate_precision` therefore takes the same group
  masks the test will use.

Simulation gives well-calibrated null p-values (nominal 0.05 → observed 0.045–0.053 across
dispersion levels) and, at 3-vs-3 with ~60× coverage, 94 % of ΔΨ = 0.2 events clearing
BH `q ≤ 0.05`.

- **Effect size:** `ΔΨ = Ψ̂(B) − Ψ̂(A)`, the fitted group means.
- **Multiple testing:** **Benjamini–Hochberg** FDR across all tested units, with the
  standard monotonicity enforcement.

Units with fewer than `min_samples` informative replicates per group are skipped.
"Significant" defaults to `q ≤ 0.05` and `|ΔΨ| ≥ 0.1`. The rank test remains available as
`test="ranksum"` for Ψ tables that carry no counts.

## 5b. Event-level PSI (cassette exons)

Splice-site usage (§4) is model-free but not an *event*. `events` additionally
reconstructs the canonical **cassette-exon** event from three junctions that share
endpoints — two inclusion junctions flanking the exon and one skipping junction across it
— and computes rMATS-style percent-spliced-in:

$$\text{inclusion} = \tfrac{1}{2}\big(\mathrm{count}(\text{inc}_1) + \mathrm{count}(\text{inc}_2)\big),\qquad
\Psi = \frac{\text{inclusion}}{\text{inclusion} + \mathrm{count}(\text{skip})}$$

Detection uses the coordinate pattern skip `(s₃, e₃)`, inc₁ `(s₃, e₁)` with `e₁ < e₃`,
inc₂ `(s₂, e₃)` with `s₂ > s₃` and `e₁ < s₂`; the cassette exon is `[e₁+1, s₂-1]`. This
pattern is strand-independent (only the strand must match). Event-level ΔΨ between
conditions reuses the differential test of §5 (keyed by `event_id`).

**Alternative 5′/3′ splice sites (A5SS / A3SS).** An A5SS event is a splice *acceptor*
served by two or more *donors*; A3SS is the mirror (one donor, several acceptors). The
inclusion isoform is the one whose variable site is closest to the shared site (the longer
exon), per the rMATS convention, and Ψ is that isoform's share of all reads at the shared
site. Junctions already used by a cassette event are excluded from A5SS/A3SS detection so
the same signal is not reported twice.

**Mutually exclusive exons (MXE).** An MXE is two non-overlapping exons, A (upstream) and
B (downstream), that lie between the *same* upstream donor and downstream acceptor and are
normally never included together (so, unlike a cassette, there is no exon-skipping
junction). Each exon is reached by its own inclusion-junction pair; Ψ reports exon A's
share, `Ψ = incl(A) / (incl(A) + incl(B))`. Detection groups every exon that connects a
given donor→acceptor pair and emits an MXE when two of them are non-overlapping; junctions
already claimed by a cassette event are excluded first.

**Intron retention (RI) is out of scope by design.** Retention is defined by *reads inside
the intron*, which junction counts do not carry — quantifying it requires intronic
coverage from the BAM/bigWig. Rather than approximate it unreliably, `splicescope` reports
the four junction-quantifiable classes (SE, MXE, A5SS, A3SS); RI is a natural extension
once coverage is available.

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

### 7c. What the classifier is and is not for

Its features (`intron_length`, `log_max_count`, `n_samples_support`, `mean_psi_donor`,
`canonical_motif`, distance to the nearest known donor/acceptor, `is_novel_both`) describe
a junction in isolation; none of them sees the experimental groups. It therefore scores
*plausibility*, not differential usage.

Measured against replication in an independent experiment (GSE122069), it reaches an
average precision of 0.006 for confidently replicating junctions where the beta-binomial
statistic reaches 0.062 — an order of magnitude worse — and combining the two is worse
than the statistic alone. Under a permissive definition of replication the ordering
reverses. Use it to filter noise before testing, not to rank cryptic candidates.

## 7b. Pathway over-representation (ORA)

To ask *which pathways are affected*, `enrich` runs the standard one-sided
hypergeometric test (equivalent to Fisher's exact) used by clusterProfiler. With `M`
background genes, a gene set of size `n` (within the background), `N` differentially
spliced genes, and `k` of them in the set, the enrichment p-value is
`P(X ≥ k) = hypergeom.sf(k−1, M, n, N)`, with fold enrichment `(k/N)/(n/M)` and
Benjamini–Hochberg correction across sets. The background is every gene that was tested;
hits are the genes of significant units. Gene sets are supplied by the user (GMT via
`io.read_gmt`) — none are bundled, because meaningful enrichment requires real
annotations, not synthetic ones.

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

## 8b. Protein consequence

For a cassette exon the reading frame is inherited from the coding sequence upstream of
the host intron, the insert is translated in that frame, the first in-frame stop is
located and the **50-nucleotide rule** decides whether it triggers NMD: a PTC more than
50 nt upstream of the last exon–exon junction is predicted to be degraded. Where several
transcripts host the event, the most disruptive prediction is reported.

Two shapes of event are handled. A **cassette** cryptic exon is supplied as an exon
interval. A **splice-site shift** — a novel donor or acceptor paired with an annotated
site on the other side — is supplied as a junction, and the annotation determines what
it adds to (extension) or removes from (truncation) the neighbouring exon. Shifts are the
majority of cryptic events reported by junction-level callers, so a layer restricted to
cassettes would miss most of them.

When the event shifts the frame without carrying a stop itself, the retained downstream
exons are assembled and the in-frame scan continues there, because that is where the
first premature stop then lies. The same applies to truncations, which cannot contain a
stop inside the removed interval by construction.

**This is interpretation, not evidence of reality.** A random 128 nt interval read in a
fixed frame contains a stop codon with probability `1 − (61/64)^42 ≈ 0.87`, so most
intronic sequence looks disruptive whether or not the event is real. On GSE245332 the
protein-disrupting fraction among significant events (61.6%) was no higher than
background (66.9%). Apply the consequence layer to events that already passed the
differential test.

## 9. Limitations

- Ψ here is splice-site *usage*, not event-level PSI (cassette exon, A5SS/A3SS, IR); the
  latter is a natural extension built on the same junction graph.
- The bundled data is simulated. On real data the classifier must be retrained on curated
  labels and validated on **held-out genes** (not just held-out junctions) to avoid
  optimistic estimates from shared-gene leakage.
- Classification uses splice-site identity, not sequence models of the splice site; a
  learned donor/acceptor motif score would refine the `canonical_motif` feature.
