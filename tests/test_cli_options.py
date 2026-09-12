"""The `run` options that had no test, and the outputs only they produce.

Coverage put `cli.py` at 80 %, and two of the three uncovered blocks were whole
features: pathway enrichment behind `--gene-sets`, and the supervised classifier. The
second could not run at all — it was guarded on a truth column that nothing an aligner
writes ever carries, so `model_card.json` and `cryptic_scores.tsv` were unreachable from
the command line while the model card told every reader to retrain on curated labels.
`--labels` is how they are handed in.
"""

from __future__ import annotations

import json

import pandas as pd
import pytest

from splicescope.cli import main
from splicescope.simulate import simulate_dataset, write_dataset


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    """A written-out simulated dataset, plus the labels it knows the truth of."""
    root = tmp_path_factory.mktemp("cli")
    data = root / "data"
    ds = simulate_dataset(n_genes=12, n_per_group=4, cryptic_fraction=0.7, seed=5)
    write_dataset(ds, data)
    labels = (
        ds.observed.groupby(["chrom", "start", "end", "strand"], observed=True)[
            "is_cryptic_truth"
        ]
        .max()
        .reset_index()
    )
    labels.to_csv(root / "labels.tsv", sep="\t", index=False)
    genes = sorted(ds.known["gene_id"].unique())
    return root, data, genes


def _run(data, outdir, *extra):
    return main(
        [
            "run",
            "--sj-dir", str(data / "sj"),
            "--gtf", str(data / "annotation.gtf"),
            "--groups", str(data / "groups.tsv"),
            "--outdir", str(outdir),
            *extra,
        ]
    )


def test_curated_labels_train_the_classifier_and_write_its_model_card(dataset, tmp_path):
    root, data, _ = dataset
    out = tmp_path / "labelled"
    assert _run(data, out, "--labels", str(root / "labels.tsv")) == 0

    assert (out / "cryptic_scores.tsv").exists()
    scores = pd.read_csv(out / "cryptic_scores.tsv", sep="\t")
    assert "cryptic_score" in scores.columns
    assert scores["cryptic_score"].between(0, 1).all()
    assert scores["cryptic_score"].is_monotonic_decreasing, "most likely first"

    card = json.loads((out / "model_card.json").read_text())
    assert card["cv_metrics"]["n_positive"] > 0
    assert card["hyperparameters"]["class_weight"] == "balanced"

    # the human-readable card must carry the same hyper-parameters as the JSON one,
    # which is the thing a reviewer reads
    markdown = (out / "model_card.md").read_text()
    assert "## Hyper-parameters" in markdown
    for key, value in card["hyperparameters"].items():
        assert f"- **{key}**: {value}" in markdown


def test_a_label_file_missing_its_columns_says_which(dataset, tmp_path, capsys):
    root, data, _ = dataset
    bad = tmp_path / "bad.tsv"
    pd.DataFrame({"chrom": ["chr1"], "start": [100]}).to_csv(bad, sep="\t", index=False)
    assert _run(data, tmp_path / "out", "--labels", str(bad)) == 2
    error = capsys.readouterr().err
    assert "is_cryptic_truth" in error and "missing" in error


def test_labels_that_match_nothing_are_refused_rather_than_silently_ignored(
    dataset, tmp_path, capsys
):
    """The same failure the sample/group check exists for: a coordinate convention off
    by one would leave nothing labelled, and a classifier trained on nothing is worse
    than no classifier."""
    root, data, _ = dataset
    wrong = tmp_path / "wrong.tsv"
    pd.DataFrame(
        {
            "chrom": ["chrZ"] * 3, "start": [1, 2, 3], "end": [9, 10, 11],
            "strand": ["+"] * 3, "is_cryptic_truth": [1, 0, 1],
        }
    ).to_csv(wrong, sep="\t", index=False)
    assert _run(data, tmp_path / "out", "--labels", str(wrong)) == 2
    error = capsys.readouterr().err
    assert "no junction" in error and "1-based" in error


def test_gene_sets_produce_an_enrichment_table_and_figure(dataset, tmp_path):
    root, data, genes = dataset
    gmt = tmp_path / "sets.gmt"
    gmt.write_text(
        "\n".join(
            [
                "\t".join(["SET_ALL", "every simulated gene", *genes]),
                "\t".join(["SET_HALF", "half of them", *genes[: len(genes) // 2]]),
            ]
        )
        + "\n"
    )
    out = tmp_path / "enriched"
    assert _run(data, out, "--gene-sets", str(gmt)) == 0

    table = pd.read_csv(out / "enrichment.tsv", sep="\t")
    assert not table.empty
    assert {"term", "overlap", "fold_enrichment", "pvalue", "qvalue"} <= set(table.columns)
    assert table["qvalue"].between(0, 1).all()
    assert (out / "figures" / "enrichment.png").exists()


def test_gene_sets_naming_other_identifiers_say_so_instead_of_reporting_nothing(
    dataset, tmp_path, capsys
):
    """The enrichment equivalent of the sample/group mismatch: a GMT keyed by symbols
    against a GTF keyed by accessions tests nothing, and "0 sets enriched" is a
    conclusion rather than an error."""
    root, data, _ = dataset
    gmt = tmp_path / "foreign.gmt"
    gmt.write_text("SET_X\tnot this organism\tYFL001C\tYFL002C\tYFL003C\n")
    assert _run(data, tmp_path / "foreign", "--gene-sets", str(gmt)) == 0

    captured = capsys.readouterr()
    assert "no gene set shares an identifier" in captured.err
    assert "gene sets name genes like" in captured.err
    assert "0 tested" in captured.out


def test_the_simulate_subcommand_writes_a_dataset_the_run_command_accepts(tmp_path, capsys):
    """The README's very first command, and it had no test — the suite built its
    datasets by calling the library directly. CI's smoke test runs it, so a break would
    have shown up there, but nothing pinned what it writes."""
    data = tmp_path / "demo"
    assert main(["simulate", "--outdir", str(data), "--genes", "6", "--seed", "1"]) == 0

    for name in ("annotation.gtf", "groups.tsv", "genome.fa", "genome.fa.fai", "truth.tsv"):
        assert (data / name).exists(), name
    assert len(list((data / "sj").glob("*.SJ.out.tab"))) == 8, "4 replicates per group"
    printed = capsys.readouterr().out
    assert "unique junctions" in printed and "truth.tsv" in printed

    # and what it wrote is what `run` consumes, genome and all
    out = tmp_path / "results"
    assert _run(data, out, "--genome", str(data / "genome.fa")) == 0
    assert (out / "differential_splicing.tsv").exists()
    assert (out / "consequence.tsv").exists()


def test_an_empty_sj_directory_is_an_error_not_an_empty_result(dataset, tmp_path, capsys):
    """Pointing --sj-dir at the wrong place is an easy mistake, and a run that reports
    nothing is the one outcome this tool must never produce quietly."""
    _, data, _ = dataset
    empty = tmp_path / "nothing"
    empty.mkdir()
    rc = main(
        [
            "run",
            "--sj-dir", str(empty),
            "--gtf", str(data / "annotation.gtf"),
            "--groups", str(data / "groups.tsv"),
            "--outdir", str(tmp_path / "out"),
        ]
    )
    assert rc == 2
    assert "no *.tab files" in capsys.readouterr().err


def test_a_sample_listed_in_groups_with_no_file_warns_and_the_run_continues(
    dataset, tmp_path, capsys
):
    """A dropped sample is a real experiment, not a typo, so this must not be fatal —
    but it must not be silent either. Only the total-mismatch case was tested."""
    _, data, _ = dataset
    groups = pd.read_csv(data / "groups.tsv", sep="\t")
    extra = pd.concat(
        [groups, pd.DataFrame([{"sample": "B99", "condition": groups["condition"].iloc[-1]}])],
        ignore_index=True,
    )
    path = tmp_path / "groups_with_a_ghost.tsv"
    extra.to_csv(path, sep="\t", index=False)

    rc = main(
        [
            "run",
            "--sj-dir", str(data / "sj"),
            "--gtf", str(data / "annotation.gtf"),
            "--groups", str(path),
            "--outdir", str(tmp_path / "ghost"),
        ]
    )
    assert rc == 0
    error = capsys.readouterr().err
    assert "B99" in error and "no *.tab file" in error


def test_the_consequence_subcommand_names_the_columns_it_needs(dataset, tmp_path, capsys):
    """Run on its own, it takes a table of events from anywhere — including rMATS or a
    hand-made one — so the column check is the first thing a new user meets."""
    _, data, _ = dataset
    events = tmp_path / "events.tsv"
    pd.DataFrame({"chrom": ["chr1"], "gene_id": ["g00"]}).to_csv(events, sep="\t", index=False)
    def run_consequence(*extra):
        return main(
            [
                "consequence",
                "--events", str(events),
                "--gtf", str(data / "annotation.gtf"),
                "--genome", str(data / "genome.fa"),
                "--out", str(tmp_path / "cons.tsv"),
                *extra,
            ]
        )

    # with no coordinate columns at all it cannot even tell which kind of table this is
    assert run_consequence() == 2
    assert "pass --mode explicitly" in capsys.readouterr().err

    # told which kind, it names the columns that kind needs
    assert run_consequence("--mode", "junction") == 2
    error = capsys.readouterr().err
    assert "missing columns" in error and "strand" in error


def test_a_gtf_from_the_wrong_source_is_an_error_with_a_readable_message(dataset, tmp_path, capsys):
    """The mistake is easy and the wrong answer is exciting: mix an Ensembl GTF with
    chr-prefixed alignments and every junction comes out cryptic. The CLI must say so
    the way it says everything else, not as a traceback."""
    _, data, _ = dataset
    ensembl = tmp_path / "ensembl.gtf"
    ensembl.write_text(
        "\n".join(
            line.replace("chr1\t", "1\t", 1)
            for line in (data / "annotation.gtf").read_text().splitlines()
        )
        + "\n"
    )
    rc = main(
        [
            "run",
            "--sj-dir", str(data / "sj"),
            "--gtf", str(ensembl),
            "--groups", str(data / "groups.tsv"),
            "--outdir", str(tmp_path / "out"),
        ]
    )
    assert rc == 2
    error = capsys.readouterr().err
    assert error.startswith("error: ")
    assert "name the same chromosomes differently" in error
    assert "GENCODE writes 'chr1' where Ensembl writes '1'" in error


def test_the_cli_says_which_step_could_not_be_run(dataset, tmp_path, capsys):
    """Library warnings reach a terminal as `warning: …` rather than as a Python warning
    with a file and a line number, and labelled, because the same sentence comes from the
    junction-level test and the event-level one."""
    _, data, _ = dataset
    rc = main(
        [
            "run",
            "--sj-dir", str(data / "sj"),
            "--gtf", str(data / "annotation.gtf"),
            "--groups", str(data / "groups.tsv"),
            "--outdir", str(tmp_path / "starved"),
            "--min-reads", "1000000",
        ]
    )
    assert rc == 0, "an impossible threshold is a choice, not a failure"
    error = capsys.readouterr().err
    assert "warning: junctions: no unit had" in error
    assert "warning: events: no unit had" in error
    assert "UserWarning" not in error, "not raw Python warnings"
