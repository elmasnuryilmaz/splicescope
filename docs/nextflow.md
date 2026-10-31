# Nextflow pipeline

`nextflow/main.nf` is a DSL2 pipeline that runs the same analysis as
`splicescope run`, as processes that Nextflow can schedule on a laptop, a cluster or the
cloud. It needs Nextflow 23.04 or newer and `splicescope` installed in the environment
the processes run in.

## Demo mode

With no inputs, the pipeline simulates a dataset (genome included) and analyses it:

```bash
nextflow run nextflow/main.nf -profile test
```

Results are published to `results/`: the simulated inputs in `results/data/` and the
analysis in `results/out/`, with the same files `splicescope run` writes (see
[Quick start](quickstart.md#3-read-the-results)). An execution report and a timeline
are written next to them as HTML.

## Real data

Give the three inputs, and optionally an indexed genome for the consequence layer. Each
is staged on its own, so they can live in different places:

```bash
nextflow run nextflow/main.nf \
    --sj_dir  /path/to/SJ_out_tabs \
    --gtf     gencode.v47.annotation.gtf.gz \
    --groups  groups.tsv \
    --genome  GRCh38.primary_assembly.genome.fa \
    --outdir  results
```

The `.fai` index must sit next to the FASTA (`GRCh38.primary_assembly.genome.fa.fai`).
Without `--genome` the consequence step is skipped and everything else runs.

## Parameters

| parameter | default | meaning |
|---|---|---|
| `--sj_dir` | none | directory of STAR `*.SJ.out.tab` files; with `--gtf` and `--groups`, switches to real-data mode |
| `--gtf` | none | reference annotation, plain or gzipped |
| `--groups` | none | TSV with columns `sample` and `condition` |
| `--genome` | none | indexed genome FASTA; enables protein-consequence prediction |
| `--min_reads` | 10 | minimum reads at a splice site in a sample for its Ψ to be computed |
| `--outdir` | `results` | where results are published |
| `--genes`, `--replicates` | 12, 6 | size of the simulated dataset in demo mode |

## Profiles

| profile | effect |
|---|---|
| `test` | a small simulated dataset on one CPU; what CI runs |
| `docker` | runs each process in the `python:3.11-bookworm` container, after installing the splicescope release that matches the pipeline's version from PyPI inside it |
| `conda` | runs each process in the conda environment `nextflow/environment.yml` (Python 3.11 and the matching splicescope release from PyPI) |
| none | runs each process in your current environment, where `splicescope` must already be installed |

Both the `docker` and the `conda` profile need network access to PyPI. With `docker`, the
install is part of each task's script, so it happens inside the container and never
touches your own Python environment. The image is the full Debian one rather than
`-slim` because Nextflow needs `ps` inside the container to record task metrics.

## What CI checks

Every push runs the pipeline three ways: demo mode, real-data mode with the inputs in
separate directories, and real-data mode without a genome (which must skip the
consequence step). See the `nextflow` job in
[`.github/workflows/ci.yml`](https://github.com/elmasnuryilmaz/splicescope/blob/main/.github/workflows/ci.yml).
