# Command-line interface

Installing the package puts one command on your `PATH`, `splicescope`, with three
subcommands. The help text below is generated from the parser when the site is built,
so it always matches the installed version.

| subcommand | what it does |
|---|---|
| `simulate` | write a synthetic dataset with known ground truth |
| `run` | the whole pipeline on a directory of STAR junction tables |
| `consequence` | protein-consequence prediction for events produced by any tool |

## `splicescope`

<!-- cli-help: -->

## `splicescope simulate`

Writes STAR `SJ.out.tab` files, a GTF annotation with coding sequences, `groups.tsv`, an
indexed genome FASTA and `truth.tsv`, which lists the injected events so that recall can
be measured. Nothing is downloaded.

<!-- cli-help: simulate -->

## `splicescope run`

Annotates, quantifies, detects events, tests for differential splicing and, when
`--genome` is given, predicts protein consequences. See [Quick start](quickstart.md) for
the files it writes.

<!-- cli-help: run -->

The groups file is a tab-separated table with the columns `sample` and `condition`, and
exactly two conditions. Sample names come from the junction file names with a trailing
`SJ.out` removed: `A1_SJ.out.tab`, `A1.SJ.out.tab`, `A1SJ.out.tab` and `A1.tab` all give
`A1`. The run stops with an error if no sample name is shared between the two.

## `splicescope consequence`

Runs only the consequence layer, on candidates produced elsewhere (for example an rMATS
or LeafCutter run). Cassette exons are given as exon intervals (`exon_start`,
`exon_end`); novel donors and acceptors as junctions (`start`, `end`). The genome needs
its `samtools faidx` index next to it.

<!-- cli-help: consequence -->

The consequence classes, the frame arithmetic and the NMD rule are defined in
[Methods, section 8b](METHODS.md#8b-protein-consequence).
