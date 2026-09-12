"""Version consistency, real-world filenames and compressed inputs.

None of this is science, and all of it is how a tool loses a user's trust before the
science is reached: a `--version` three releases out of date, a GENCODE `.gtf.gz` that
raises a decode error, or a run that reports nothing because STAR's filenames were not
the ones the CLI expected.
"""

import gzip
import re
from pathlib import Path

import pandas as pd
import pytest

import splicescope
from splicescope.cli import main
from splicescope.io import read_gmt, read_gtf_junctions, sample_name_from_path
from splicescope.simulate import simulate_dataset, write_dataset

ROOT = Path(__file__).resolve().parent.parent


def test_version_is_consistent_across_the_repository():
    """__version__ was 0.5.0 while pyproject and CITATION.cff said 0.8.1, so the CLI
    reported a release three versions old."""
    version = splicescope.__version__

    citation = (ROOT / "CITATION.cff").read_text()
    assert re.search(rf'^version: "{re.escape(version)}"$', citation, re.M), (
        f"CITATION.cff does not declare {version}"
    )

    changelog = (ROOT / "CHANGELOG.md").read_text()
    released = re.findall(r"^## \[(\d+\.\d+\.\d+)\]", changelog, re.M)
    assert released, "no released version headings in CHANGELOG.md"
    assert released[0] == version, (
        f"CHANGELOG's newest release is {released[0]}, but __version__ is {version}"
    )


def test_pyproject_reads_the_version_from_the_package():
    """A literal version in pyproject.toml is a second source of truth, and it drifted."""
    pyproject = (ROOT / "pyproject.toml").read_text()
    assert 'dynamic = ["version"]' in pyproject
    assert 'attr = "splicescope.__version__"' in pyproject
    assert not re.search(r"^version = \"\d", pyproject, re.M)


@pytest.mark.parametrize(
    ("filename", "expected"),
    [
        ("A1_SJ.out.tab", "A1"),      # STAR with an underscore-terminated prefix
        ("A1.SJ.out.tab", "A1"),      # STAR with a dot-terminated prefix
        ("A1SJ.out.tab", "A1"),       # STAR with a bare prefix
        ("A1-SJ.out.tab", "A1"),
        ("A1.tab", "A1"),             # already-renamed files
        ("sample_02_SJ.out.tab", "sample_02"),
    ],
)
def test_sample_names_survive_stars_filename_conventions(filename, expected):
    assert sample_name_from_path(filename) == expected


def test_run_fails_loudly_when_no_sample_name_matches(tmp_path):
    """Previously this printed "0 significant junctions" and exited 0 — a clean run and
    a wrong conclusion."""
    data = tmp_path / "data"
    write_dataset(simulate_dataset(n_genes=4, n_per_group=2, seed=1), data)
    groups = pd.read_csv(data / "groups.tsv", sep="\t")
    groups["sample"] = "prefix_" + groups["sample"]
    groups.to_csv(data / "renamed.tsv", sep="\t", index=False)

    rc = main(
        [
            "run",
            "--sj-dir", str(data / "sj"),
            "--gtf", str(data / "annotation.gtf"),
            "--groups", str(data / "renamed.tsv"),
            "--outdir", str(tmp_path / "out"),
        ]
    )
    assert rc == 2


def test_run_works_with_stars_actual_filenames(tmp_path):
    data = tmp_path / "data"
    write_dataset(simulate_dataset(n_genes=8, n_per_group=3, seed=1), data)
    for path in (data / "sj").glob("*.SJ.out.tab"):
        path.rename(path.with_name(path.name.replace(".SJ.out.tab", "_SJ.out.tab")))

    out = tmp_path / "out"
    rc = main(
        [
            "run",
            "--sj-dir", str(data / "sj"),
            "--gtf", str(data / "annotation.gtf"),
            "--groups", str(data / "groups.tsv"),
            "--outdir", str(out),
        ]
    )
    assert rc == 0
    result = pd.read_csv(out / "differential_splicing.tsv", sep="\t")
    assert not result.empty, "renaming the files must not empty the differential table"


def test_gtf_can_be_gzipped(tmp_path):
    """GENCODE and Ensembl ship .gtf.gz, and the README's own example passes one."""
    data = tmp_path / "data"
    write_dataset(simulate_dataset(n_genes=4, seed=1), data)
    plain = data / "annotation.gtf"
    packed = data / "annotation.gtf.gz"
    with gzip.open(packed, "wt") as fh:
        fh.write(plain.read_text())

    pd.testing.assert_frame_equal(read_gtf_junctions(plain), read_gtf_junctions(packed))


def test_gmt_can_be_gzipped(tmp_path):
    content = "TERM_A\tdescription\tG1\tG2\nTERM_B\tdescription\tG2\tG3\n"
    plain = tmp_path / "sets.gmt"
    plain.write_text(content)
    packed = tmp_path / "sets.gmt.gz"
    with gzip.open(packed, "wt") as fh:
        fh.write(content)

    assert read_gmt(plain) == read_gmt(packed) == {"TERM_A": ["G1", "G2"], "TERM_B": ["G2", "G3"]}


def _numeric_cohort(tmp_path, names):
    """A dataset whose samples are named with digits, as sequencing run IDs often are."""
    data = tmp_path / "data"
    write_dataset(simulate_dataset(n_genes=8, n_per_group=3, seed=1), data)
    renamed = data / "numeric"
    renamed.mkdir()
    for name, path in zip(names, sorted((data / "sj").glob("*.tab")), strict=True):
        (renamed / f"{name}_SJ.out.tab").write_bytes(path.read_bytes())
    groups = data / "numeric_groups.tsv"
    groups.write_text(
        "sample\tcondition\n"
        + "".join(f"{n}\t{'A' if i < 3 else 'B'}\n" for i, n in enumerate(names))
    )
    return data, renamed, groups


@pytest.mark.parametrize(
    "names",
    [
        ["101", "102", "103", "104", "105", "106"],   # plain digits
        ["001", "002", "003", "004", "005", "006"],   # zero-padded: "007" must not become 7
    ],
)
def test_digit_only_sample_names_are_labels_not_numbers(tmp_path, names):
    """pandas reads a digit-only sample column as int64, which can never equal the strings
    derived from filenames. The mismatch check then crashed on `str.join` over ints."""
    data, sj_dir, groups = _numeric_cohort(tmp_path, names)
    out = tmp_path / f"out_{names[0]}"

    rc = main(
        [
            "run",
            "--sj-dir", str(sj_dir),
            "--gtf", str(data / "annotation.gtf"),
            "--groups", str(groups),
            "--outdir", str(out),
        ]
    )
    assert rc == 0
    result = pd.read_csv(out / "differential_splicing.tsv", sep="\t")
    assert not result.empty, "a consistently named numeric cohort must produce results"


def test_stars_strand_and_motif_codes_are_read_as_star_writes_them(tmp_path):
    """Column 4 of SJ.out.tab is 0, 1 or 2, and 0 means STAR could not tell — which is
    routine for non-canonical junctions, and the reason `resolve_unstranded` exists.
    Reading 0 as "+" would silently invent a strand for every one of them, place them
    against the wrong half of the annotation, and leave nothing to resolve. Nothing in
    the suite read a real SJ.out.tab and checked the mapping.
    """
    import pandas as pd

    from splicescope.io import read_star_sj

    # chrom start end strand motif annotated n_unique n_multi overhang
    path = Path(tmp_path) / "S1.SJ.out.tab"
    path.write_text(
        "chr1\t100\t200\t1\t1\t1\t30\t2\t40\n"   # + strand, GT/AG
        "chr1\t300\t400\t2\t2\t1\t25\t0\t38\n"   # - strand, CT/AC
        "chr1\t500\t600\t0\t0\t0\t7\t1\t20\n"    # strand unknown, non-canonical
        "chr1\t700\t800\t1\t3\t0\t9\t0\t22\n"    # + strand, GC/AG
    )
    df = read_star_sj(path)

    assert list(df["strand"]) == ["+", "-", ".", "+"]
    assert list(df["motif"]) == ["GT/AG", "CT/AC", "non-canonical", "GC/AG"]
    # `count` is the uniquely-mapping reads, not the total
    assert list(df["count"]) == [30, 25, 7, 9]
    assert list(df["sample"]) == ["S1"] * 4
    assert pd.api.types.is_integer_dtype(df["count"])
