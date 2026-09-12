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
