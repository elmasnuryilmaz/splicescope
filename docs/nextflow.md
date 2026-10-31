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
| `docker` | runs each process in the `python:3.11-slim` container |
| `conda` | runs each process in a conda environment with Python 3.11 and pip |

!!! warning "The `docker` and `conda` profiles provide Python, not splicescope"
    The container and the conda environment named in `nextflow/nextflow.config` do not
    contain `splicescope`, so with those profiles the processes cannot find the
    `splicescope` command. Until a container image with the package is published, run
    the pipeline without them, in an environment where `splicescope` is installed.

## What CI checks

Every push runs the pipeline three ways: demo mode, real-data mode with the inputs in
separate directories, and real-data mode without a genome (which must skip the
consequence step). See the `nextflow` job in
[`.github/workflows/ci.yml`](https://github.com/elmasnuryilmaz/splicescope/blob/main/.github/workflows/ci.yml).
