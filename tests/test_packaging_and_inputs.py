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

    # The README carries the citation a reader will copy, and it was left at v0.8.1
    # through two releases because nothing checked it.
    readme = (ROOT / "README.md").read_text()
    cited = re.findall(r"splice junctions\* \(v(\d+\.\d+\.\d+)\)", readme)
    assert cited, "no version in the README's suggested citation"
    assert cited == [version], (
        f"the README cites {cited}, but __version__ is {version}"
    )


def test_the_readme_does_not_offer_an_older_release_as_this_version():
    """A reader who wants to cite this exact version was sent to v0.8.1's DOI.

    `CITATION.cff` lists two identifiers: a concept DOI that resolves to the latest
    release, and a version DOI whose own description says it belongs to v0.8.1. The
    README offered the second as "this exact version" while the first line of the same
    citation says v0.9.1, so the two disagreed about what the reader would be citing.
    Only a release mints a new version DOI, so the honest fix is to say which version the
    existing one belongs to.

    This checks the two files agree: a DOI that `CITATION.cff` attributes to some other
    version must not be offered by the README as the current one.
    """
    import re

    from splicescope import __version__ as version

    citation = (ROOT / "CITATION.cff").read_text()
    readme = (ROOT / "README.md").read_text()

    #: DOI -> the version its description names, for the entries that name one.
    attributed = {
        doi: found.group(1)
        for doi, description in re.findall(
            r'value: "(10\.\d+/[^"]+)"\s*\n\s*description: "([^"]*)"', citation
        )
        if (found := re.search(r"v(\d+\.\d+\.\d+)", description))
    }
    assert attributed, "CITATION.cff names no version for any DOI"

    for doi, belongs_to in attributed.items():
        if belongs_to == version:
            continue
        for sentence in re.split(r"(?<=[.])\s", readme):
            if doi in sentence and "exact version" in sentence:
                assert belongs_to in sentence, (
                    f"the README offers {doi} as this exact version, but CITATION.cff "
                    f"says it belongs to v{belongs_to} and this is v{version}"
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


def test_the_oldest_constraints_match_the_declared_floors():
    """`constraints-oldest.txt` is what CI installs to prove the declared support range
    is real. It only proves anything while it stays in step with `pyproject.toml`: raise
    a floor without repinning and the job silently tests a version the package no longer
    claims to support, or a dependency added without a pin is never exercised at its
    floor at all.

    Every *runtime* dependency has to be pinned. A pin may also name something from an
    optional extra — `hypothesis` is, because the property tests import it at module
    level and a floor nobody installs is a promise rather than a fact — but those are not
    required, since pinning a linter at its floor would test which rules existed rather
    than anything about this package."""
    from packaging.requirements import Requirement
    from packaging.version import Version

    try:  # tomllib is stdlib from 3.11; pytest brings tomli on 3.10, which we support
        import tomllib
    except ModuleNotFoundError:  # pragma: no cover - only on Python 3.10
        import tomli as tomllib

    project = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]
    runtime = {
        (r := Requirement(spec)).name: r.specifier for spec in project["dependencies"]
    }
    optional = {
        (r := Requirement(spec)).name: r.specifier
        for specs in project.get("optional-dependencies", {}).values()
        for spec in specs
    }
    declared = {**optional, **runtime}
    pinned = {}
    for line in (ROOT / "constraints-oldest.txt").read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            name, _, version = line.partition("==")
            assert version, f"constraints must pin exactly, got {line!r}"
            pinned[name] = Version(version)

    assert not set(pinned) - set(declared), (
        f"pinned but declared nowhere: {sorted(set(pinned) - set(declared))}"
    )
    assert not set(runtime) - set(pinned), (
        f"a runtime dependency with no floor anybody installs: "
        f"{sorted(set(runtime) - set(pinned))}"
    )
    for name, version in pinned.items():
        assert declared[name].contains(version), (
            f"{name}=={version} does not satisfy the declared {name}{declared[name]}"
        )


def test_the_package_tells_type_checkers_its_annotations_are_real():
    """PEP 561: without a `py.typed` marker shipped in the wheel, a type checker in
    someone else's project ignores every annotation in this package. Sixty of the
    sixty-five public functions carried them and none of it was visible."""
    try:
        import tomllib
    except ModuleNotFoundError:  # pragma: no cover - only on Python 3.10
        import tomli as tomllib

    marker = ROOT / "src" / "splicescope" / "py.typed"
    assert marker.exists(), "the marker itself"

    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text())
    package_data = pyproject["tool"]["setuptools"].get("package-data", {})
    assert "py.typed" in package_data.get("splicescope", []), (
        "the marker exists but is not declared as package data, so it never reaches "
        "the wheel and nothing changes for anyone installing this"
    )


def test_every_public_function_is_annotated():
    """What the marker above promises. A partially annotated package is worse than an
    unannotated one: the checker trusts what is there and infers `Any` for the rest."""
    import ast

    unannotated = []
    for path in sorted((ROOT / "src" / "splicescope").glob("*.py")):
        for node in ast.parse(path.read_text()).body:
            if not isinstance(node, ast.FunctionDef) or node.name.startswith("_"):
                continue
            arguments = [a for a in node.args.args if a.arg not in {"self", "cls"}]
            if node.returns is None and not any(a.annotation for a in arguments):
                unannotated.append(f"{path.name}:{node.name}")
    assert not unannotated, f"public functions with no annotations at all: {unannotated}"


def _gtf_line(feature, start, end, tx="T1", gene="G1", chrom="chr1", extra=""):
    attrs = f'gene_id "{gene}"; transcript_id "{tx}"; gene_name "{gene}";{extra}'
    return f"{chrom}\tsrc\t{feature}\t{start}\t{end}\t.\t+\t.\t{attrs}"


def test_a_real_gtf_s_awkward_lines_are_skipped_rather_than_parsed(tmp_path):
    """Everything here occurs in annotations people actually download. GENCODE files
    open with several `##` lines; a transcript with one exon has no intron to
    contribute; some annotations carry exons that overlap within one transcript, which
    would otherwise produce an intron running backwards."""
    path = tmp_path / "awkward.gtf"
    path.write_text(
        "\n".join(
            [
                "##description: evidence-based annotation",
                "#!genome-build GRCh38",
                "",
                _gtf_line("gene", 100, 900),          # not an exon
                _gtf_line("transcript", 100, 900),    # not an exon either
                _gtf_line("exon", 100, 200),          # T1 exon 1
                _gtf_line("exon", 401, 500),          # T1 exon 2  -> intron 201-400
                # a transcript whose exons overlap: the "intron" would run backwards
                _gtf_line("exon", 100, 300, tx="T2", gene="G2"),
                _gtf_line("exon", 250, 400, tx="T2", gene="G2"),
                # a single-exon transcript contributes no intron at all
                _gtf_line("exon", 600, 700, tx="T3", gene="G3"),
                "chr1\tsrc\texon\t800\t900\t.\t+\t.\tgene_id \"G4\";",  # no transcript_id
                "chr1\tsrc\texon\t950",                                 # truncated line
            ]
        )
        + "\n"
    )
    known = read_gtf_junctions(path)

    assert len(known) == 1, f"only T1 has an intron, got {known.to_dict('records')}"
    row = known.iloc[0]
    assert (row["chrom"], row["start"], row["end"], row["strand"]) == ("chr1", 201, 400, "+")
    assert row["gene_id"] == "G1"


def test_a_gtf_of_single_exon_transcripts_yields_no_junctions(tmp_path):
    """Which is not an error here, but it is one downstream: with nothing to match
    against, every observed junction would be called cryptic. `annotate_junctions` is
    what refuses it."""
    import pandas as pd
    import pytest

    from splicescope.annotate import annotate_junctions

    path = tmp_path / "single.gtf"
    path.write_text(
        _gtf_line("exon", 100, 200, tx="T1") + "\n" + _gtf_line("exon", 600, 700, tx="T2") + "\n"
    )
    known = read_gtf_junctions(path)
    assert known.empty
    assert list(known.columns) == ["chrom", "start", "end", "strand", "gene_id", "gene_name"]

    observed = pd.DataFrame(
        [dict(chrom="chr1", start=201, end=400, strand="+", sample="s1", count=30)]
    )
    with pytest.raises(ValueError, match="contains no junctions"):
        annotate_junctions(observed, known)


def test_no_sj_files_gives_an_empty_table_with_the_right_columns(tmp_path):
    """`read_many_star_sj({})` is reachable from the library even though the CLI
    refuses an empty directory first, and an empty frame with no columns is the shape
    that breaks every caller downstream."""
    from splicescope.io import read_many_star_sj

    out = read_many_star_sj({})
    assert out.empty
    assert list(out.columns) == [
        "chrom", "start", "end", "strand", "motif", "annotated_star", "count", "sample",
    ]


def test_everything_the_suite_imports_is_a_declared_dependency():
    """Six pushes ran red before anyone looked at CI.

    `tests/test_properties.py` imports `hypothesis` at module level and nothing declared
    it, so the module failed to collect on every Python in the matrix and on the oldest
    job. Locally it was installed and the suite was green, which is the whole problem: a
    dependency you already have is invisible until someone else installs the package.

    This walks every import in the repository and checks it against `pyproject.toml`.
    Modules with a documented fallback are listed below with the reason, so adding one is
    a deliberate act rather than an omission.
    """
    import ast
    import sys

    try:
        import tomllib
    except ModuleNotFoundError:  # pragma: no cover - only on Python 3.10
        import tomli as tomllib

    from packaging.requirements import Requirement

    #: Third-party names that need no declaration, and why.
    ALLOWED = {
        "splicescope": "this package",
        "tomllib": "stdlib from 3.11; on 3.10 the import is guarded and falls back to tomli",
        "tomli": "stdlib tomllib from 3.11; pytest brings tomli on 3.10, guarded by try",
        "packaging": "a pip and setuptools dependency, present wherever this installs",
        "pytest_cov": "a pytest plugin, never imported",
    }
    #: import name -> distribution name, where they differ.
    IMPORT_NAMES = {"sklearn": "scikit-learn", "PIL": "pillow", "yaml": "pyyaml"}

    project = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]
    declared = {
        Requirement(spec).name.lower().replace("-", "_")
        for spec in project["dependencies"]
        + [s for specs in project.get("optional-dependencies", {}).values() for s in specs]
    }

    stdlib = set(sys.stdlib_module_names)
    missing: dict[str, set[str]] = {}
    for directory in ("src", "tests", "app", "examples", "validation", "docs"):
        for path in sorted((ROOT / directory).rglob("*.py")):
            for node in ast.walk(ast.parse(path.read_text())):
                if isinstance(node, ast.Import):
                    names = [alias.name.split(".")[0] for alias in node.names]
                elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                    names = [node.module.split(".")[0]]
                else:
                    continue
                for name in names:
                    if name in stdlib or name in ALLOWED:
                        continue
                    dist = IMPORT_NAMES.get(name, name).lower().replace("-", "_")
                    if dist not in declared:
                        missing.setdefault(name, set()).add(
                            str(path.relative_to(ROOT))
                        )

    assert not missing, "imported but declared nowhere in pyproject.toml:\n  " + "\n  ".join(
        f"{name} — {sorted(where)[:3]}" for name, where in sorted(missing.items())
    )


def test_the_job_that_runs_the_suite_installs_what_the_suite_imports():
    """The other half of the test above, and the half that was actually red.

    Declaring a dependency is not the same as installing it where the tests run. CI
    installed `.[dev]` while the suite imports `nbformat` from `.[docs]` — to check that
    the committed tutorial is the one its builder produces — so that test failed on every
    Python in the matrix and passed on any machine that had ever rebuilt the notebook.

    This maps each import to the extra that declares it and compares that against what
    the workflow installs, so an extra added to `pyproject.toml` without being added to
    CI is caught here rather than on the next push.
    """
    import ast
    import re
    import sys

    try:
        import tomllib
    except ModuleNotFoundError:  # pragma: no cover - only on Python 3.10
        import tomli as tomllib

    from packaging.requirements import Requirement

    IMPORT_NAMES = {"sklearn": "scikit-learn", "PIL": "pillow", "yaml": "pyyaml"}
    IGNORE = {"splicescope", "tomli", "tomllib", "packaging", "pytest_cov"}
    #: pillow arrives with matplotlib, so `docs` is not what puts it there.
    TRANSITIVE = {"pillow": "matplotlib"}
    #: `app/` is never imported by the suite — the page is executed through Streamlit's
    #: own harness, in a module that skips itself when that extra is absent.
    SEARCHED = ("src", "tests", "examples", "validation", "docs")

    project = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]
    extra_of = {}
    for name in (Requirement(s).name for s in project["dependencies"]):
        extra_of[name.lower()] = None  # runtime: always installed
    for extra, specs in project.get("optional-dependencies", {}).items():
        for spec in specs:
            extra_of.setdefault(Requirement(spec).name.lower(), extra)

    stdlib = set(sys.stdlib_module_names)
    needed: dict[str, str] = {}
    for directory in SEARCHED:
        for path in sorted((ROOT / directory).rglob("*.py")):
            tree = ast.parse(path.read_text())
            # a module that calls `pytest.importorskip` for something is asking to be
            # skipped where it is missing, which is a declaration that it is optional
            skipped = {
                ast.literal_eval(call.args[0])
                for call in ast.walk(tree)
                if isinstance(call, ast.Call)
                and isinstance(call.func, ast.Attribute)
                and call.func.attr == "importorskip"
                and call.args
                and isinstance(call.args[0], ast.Constant)
            }
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    names = [a.name.split(".")[0] for a in node.names]
                elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                    names = [node.module.split(".")[0]]
                else:
                    continue
                for name in names:
                    if name in stdlib or name in IGNORE or name in skipped:
                        continue
                    dist = IMPORT_NAMES.get(name, name).lower()
                    dist = TRANSITIVE.get(dist, dist)
                    extra = extra_of.get(dist)
                    if extra is not None:
                        needed[extra] = f"{name} in {path.relative_to(ROOT)}"

    workflow = (ROOT / ".github" / "workflows" / "ci.yml").read_text()
    step = workflow.index("- name: Test (pytest)")
    installs = re.findall(r'pip install -e "\.\[([^\]]+)\]"', workflow[:step])
    assert installs, "no editable install found before the step that runs the suite"
    provided = {e.strip() for e in installs[-1].split(",")}

    uncovered = {e: why for e, why in needed.items() if e not in provided}
    assert not uncovered, (
        f"CI installs {sorted(provided)} before running the suite, but it also imports "
        + "; ".join(f"{why} (extra {extra!r})" for extra, why in sorted(uncovered.items()))
    )

    # and what a person is told to install has to be at least what CI installs, or the
    # first thing they see after following the README is a failure CI never shows them
    for doc in ("README.md", "CONTRIBUTING.md"):
        documented = re.findall(r'pip install -e "\.\[([^\]]+)\]"', (ROOT / doc).read_text())
        assert documented, f"{doc} documents no editable install"
        offered = {e.strip() for line in documented for e in line.split(",")}
        assert needed.keys() <= offered, (
            f"{doc} tells a reader to install {sorted(offered)}, which leaves out "
            + ", ".join(sorted(needed.keys() - offered))
        )
