# Quick start

This walks through the whole analysis on the built-in synthetic dataset. Nothing is
downloaded, and it takes a few seconds.

## 1. Make a dataset

```bash
splicescope simulate --outdir demo_data
```

`demo_data/` now holds what a real experiment would provide, plus the answer key:

| file | contents |
|---|---|
| `sj/*.SJ.out.tab` | one STAR junction table per sample, four replicates in each of two groups |
| `annotation.gtf` | the reference annotation, with coding sequences |
| `groups.tsv` | `sample` and `condition` for every sample |
| `genome.fa`, `genome.fa.fai` | an indexed genome whose genes carry real open reading frames |
| `truth.tsv` | the events that were injected, so that recall can be measured |

## 2. Run the pipeline

```bash
splicescope run \
    --sj-dir demo_data/sj \
    --gtf    demo_data/annotation.gtf \
    --groups demo_data/groups.tsv \
    --genome demo_data/genome.fa \
    --outdir results
```

Leave out `--genome` and everything except the consequence step still runs.

## 3. Read the results

| file | one row per |
|---|---|
| `annotation_summary.tsv` | junction class: how many junctions are `annotated`, `novel_donor`, `novel_acceptor`, `novel_combination` or `cryptic` |
| `differential_splicing.tsv` | junction: mean Ψ in each group, ΔΨ, beta-binomial p-value and Benjamini-Hochberg q-value |
| `events.tsv` | splicing event (SE, MXE, A5SS, A3SS) with its coordinates |
| `event_differential.tsv` | event: event-level percent-spliced-in per group, ΔPSI, p and q |
| `consequence.tsv` | cassette exon: frame, premature stop, NMD prediction and consequence class |
| `junction_consequence.tsv` | splice-site shift: the same, for novel donors and acceptors |
| `enrichment.tsv` | gene set, when `--gene-sets` is given |
| `cryptic_scores.tsv`, `model_card.json` | junction classifier output, when curated `--labels` are given |
| `figures/` | volcano plots, class summaries and consequence panels |

The columns and the statistics behind them are defined in [Methods](METHODS.md).

## 4. Or do it from Python

Every step is a function that takes and returns a `pandas` data frame:

```python
from splicescope import annotate, diff, quantify
from splicescope.simulate import simulate_dataset

ds = simulate_dataset(n_genes=20, n_per_group=6, seed=11)
ann = annotate.annotate_junctions(ds.observed, ds.known)   # junction taxonomy
psi = quantify.compute_psi(ann)                             # splice-site usage
dpsi = diff.differential_splicing(psi, ds.groups)           # ΔΨ, LRT, q-values
print(dpsi.sort_values("qvalue").head())
```

`splicescope.demo.run_demo()` runs all of it, protein consequences included, in one call.
The [tutorial notebook](https://github.com/elmasnuryilmaz/splicescope/blob/main/examples/tutorial.ipynb)
goes through the same analysis step by step, with its outputs, and the
[API reference](api.md) lists every public function.

## Next steps

- Real data: point `--sj-dir`, `--gtf` and `--groups` at your own files. Use the full
  GENCODE annotation rather than the `basic` subset; see the note in the README.
- Many samples or a cluster: the [Nextflow pipeline](nextflow.md).
- Exploring one event interactively: the [dashboard](dashboard.md).
