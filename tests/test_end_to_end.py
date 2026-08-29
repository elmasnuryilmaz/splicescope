from splicescope.cli import main
from splicescope.simulate import simulate_dataset, write_dataset


def test_cli_simulate_then_run(tmp_path):
    data = tmp_path / "data"
    ds = simulate_dataset(n_genes=10, n_per_group=4, cryptic_fraction=0.7, seed=5)
    write_dataset(ds, data)
    assert (data / "annotation.gtf").exists()
    assert (data / "groups.tsv").exists()
    assert list((data / "sj").glob("*.tab"))

    out = tmp_path / "results"
    rc = main(
        [
            "run",
            "--sj-dir",
            str(data / "sj"),
            "--gtf",
            str(data / "annotation.gtf"),
            "--groups",
            str(data / "groups.tsv"),
            "--outdir",
            str(out),
        ]
    )
    assert rc == 0
    assert (out / "annotation_summary.tsv").exists()
    assert (out / "differential_splicing.tsv").exists()
    assert (out / "figures" / "annotation_summary.png").exists()
    # STAR SJ.out.tab files carry no truth labels, so the file-based run
    # correctly skips the (supervised) cryptic ML step — see test_simulate_ml
    # for the labelled, in-memory ML validation.
    assert not (out / "model_card.json").exists()
