"""Numbers in the documentation, checked against the code that produces them.

A measurement quoted in `docs/METHODS.md` and the script advertised as reproducing it
are two sources of truth, and they drift: the weighted ORA's false-positive rate moved
from 1.9 % to 2.1 % when the propensity estimate gained its Jeffreys smoothing, and the
table went on saying 1.9 % because nothing compared them.
"""

from __future__ import annotations

import importlib.util
import re
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent

#: One row of the METHODS §7b table: name, rate at p<=0.05, rate at p<=1e-6, median p.
_TABLE_ROW = re.compile(
    r"^\| (plain hypergeometric|opportunity-weighted) \| \*\*([\d.]+) %\*\* "
    r"\| ([\d.]+) % \| ([\d.]+) \|$",
    re.M,
)


def _load_docs(script: str):
    """Import a `docs/` generator by path, the same way as a validation script."""
    return _load_from(ROOT / "docs" / f"{script}.py", script)


def _load(script: str):
    """Import a `validation/` script by path — they are programs, not a package."""
    return _load_from(ROOT / "validation" / f"{script}.py", script)


def _load_from(path, script: str):
    spec = importlib.util.spec_from_file_location(script, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _at_documented_precision(value: float, documented: str) -> float:
    """Round a computed value to however many decimals the document quotes."""
    return round(value, len(documented.partition(".")[2]))


def test_the_ora_calibration_table_matches_the_script_that_produces_it():
    """METHODS §7b quotes a false-positive rate for the uncorrected and the corrected
    test over 480 biologically null gene sets, and points the reader at
    `validation/ora_bias_calibration.py`. This runs it and compares."""
    methods = (ROOT / "docs" / "METHODS.md").read_text()
    rows = {name: rest for name, *rest in _TABLE_ROW.findall(methods)}
    assert set(rows) == {"plain hypergeometric", "opportunity-weighted"}, (
        "could not parse the calibration table in METHODS §7b"
    )

    results = _load("ora_bias_calibration").measure()
    assert len(results["plain hypergeometric"]) == 480, "the table says 480 gene sets"

    for name, (at_05, at_1e6, median) in rows.items():
        pvalues = results[name]
        for computed, documented, where in (
            ((pvalues <= 0.05).mean() * 100, at_05, "p <= 0.05"),
            ((pvalues <= 1e-6).mean() * 100, at_1e6, "p <= 1e-6"),
            (float(np.median(pvalues)), median, "the median p"),
        ):
            assert _at_documented_precision(computed, documented) == float(documented), (
                f"{name} at {where}: METHODS says {documented}, the script gives {computed:.4f}"
            )


#: The two tables of METHODS §5.4, keyed by the name `dispersion_trade.measure` uses.
_DISPERSION_ROWS = (
    (r"\| tightly dispersed \| 200 \| ([\d.]+) \| ([\d.]+) \|",
     ("heterogeneous shared estimate", "heterogeneous tight, shared")),
    (r"\| loosely dispersed \| 5 \| [\d.]+ \| \*\*([\d.]+)\*\* \|",
     ("heterogeneous loose, shared",)),
    (r"\| all \| — \| [\d.]+ \| ([\d.]+) \|",
     ("heterogeneous all, shared",)),
    (r"\| `shared` \(default\) \| ([\d.]+) \| ([\d.]+) \| ([\d.]+) \| ([\d.]+) \|",
     ("heterogeneous loose, shared", "heterogeneous all, shared",
      "homogeneous, shared", "power, shared")),
    (r"\| `per_unit_floor` \| ([\d.]+) \| ([\d.]+) \| ([\d.]+) \| ([\d.]+) \|",
     ("heterogeneous loose, floor", "heterogeneous all, floor",
      "homogeneous, floor", "power, floor")),
)


def test_the_dispersion_trade_tables_match_the_script_that_produces_them():
    """METHODS §5.4 is the argument for leaving `dispersion="shared"` as the default:
    it costs roughly a sixth of the power to halve the false-positive rate on loosely
    dispersed units. Those are the numbers a reader would weigh, and the suite's own
    tests only check the direction — that the loose half over-calls and the floor
    lowers it — with bounds loose enough for any of the figures to drift.

    The tolerance is absolute 0.005. The simulations are seeded, so this reproduces
    exactly on one machine; the slack is for NumPy's right to change a generator's
    stream between versions, and is still far below the drift that would matter (the
    figures this replaced were off by 0.02).
    """
    methods = (ROOT / "docs" / "METHODS.md").read_text()
    documented: dict[str, str] = {}
    for pattern, names in _DISPERSION_ROWS:
        found = re.search(pattern, methods)
        assert found, f"could not parse §5.4 with {pattern}"
        for name, value in zip(names, found.groups(), strict=True):
            if name in documented:
                assert documented[name] == value, (
                    f"§5.4 states {name} twice and disagrees: {documented[name]} vs {value}"
                )
            documented[name] = value

    measured = _load("dispersion_trade").measure()
    assert set(documented) <= set(measured), sorted(set(documented) - set(measured))
    for name, value in documented.items():
        tolerance = 0.05 if "estimate" in name else 0.005
        assert abs(measured[name] - float(value)) <= tolerance, (
            f"{name}: METHODS says {value}, the script gives {measured[name]:.4f}"
        )


def test_the_readme_s_python_examples_run():
    """A snippet that stopped working is a visible defect: it is the first code a reader
    tries, and nothing else in the suite touches it. Both of the README's examples had
    something wrong when this was written — one assigned `evaluate`'s metrics dict to a
    variable called `clf`, so the obvious next line, `clf.fit(...)`, would have raised."""
    import re

    blocks = re.findall(r"```python\n(.*?)```", (ROOT / "README.md").read_text(), re.S)
    assert len(blocks) >= 2, "the quickstart and the one-call example"

    for index, source in enumerate(blocks):
        namespace: dict = {}
        try:
            exec(compile(source, f"README.md[python block {index}]", "exec"), namespace)
        except Exception as exc:  # noqa: BLE001 - the point is to report any failure
            raise AssertionError(
                f"README python block {index} failed: {type(exc).__name__}: {exc}\n"
                f"{source}"
            ) from exc
        assert namespace, f"block {index} defined nothing — is it really runnable?"


def test_the_nmd_animation_flips_where_the_rule_says_it_does():
    """The animation on the README is the tool's central claim in one picture, and a
    picture can say what the code does not. Every frame's verdict comes from
    `_nmd_from_downstream`, and this checks the sequence of verdicts: the call has to
    change exactly once, between a distance above 50 and one below it."""
    from splicescope.consequence import NMD_DISTANCE_RULE

    generator = _load_docs("make_nmd_rule_gif")
    walk = range(300, generator.FIRST_EXON - 9, 12)
    verdicts = [generator.verdict(stop) for stop in walk]

    distances = [distance for _, distance, _ in verdicts]
    calls = [nmd for _, _, nmd in verdicts]
    assert distances == sorted(distances, reverse=True), "the stop walks towards the junction"
    assert calls[0] is True and calls[-1] is False, "and the fate changes along the way"

    flips = [i for i in range(1, len(calls)) if calls[i] != calls[i - 1]]
    assert len(flips) == 1, f"the call must change exactly once, changed at {flips}"
    before, after = distances[flips[0] - 1], distances[flips[0]]
    assert before > NMD_DISTANCE_RULE >= after, (
        f"the flip straddles {NMD_DISTANCE_RULE}: {before} nt then {after} nt"
    )


def test_the_committed_tutorial_is_the_one_the_builder_produces():
    """`examples/tutorial.ipynb` ships its executed outputs, and those outputs are a
    claim: this is what the current code does. The claim is only checkable if rebuilding
    the notebook and diffing it means something, and for a long time it did not — every
    rebuild rewrote 24 random cell ids and four execution timestamps per cell, so a real
    change arrived as four lines inside two hundred lines of noise. That is how the
    tutorial came to report 43 significant junctions where the code found 46.

    Executing the notebook takes minutes, so this does not do it. It checks the two
    things that make the rebuild-and-diff check work: the file carries nothing that
    records *when* it ran, and its cells are still the builder's cells — a notebook
    edited in Jupyter and committed would pass neither, and would be silently discarded
    by the next rebuild.
    """
    import json

    builder = _load_from(ROOT / "examples" / "_build_tutorial.py", "_build_tutorial")
    committed = json.loads((ROOT / "examples" / "tutorial.ipynb").read_text())
    expected = builder.build()

    assert len(committed["cells"]) == len(expected.cells), "a cell was added or removed"
    for i, (got, want) in enumerate(zip(committed["cells"], expected.cells, strict=True)):
        assert got["source"] == want["source"].splitlines(keepends=True), (
            f"cell {i} is not the builder's — edit examples/_build_tutorial.py, not the "
            f".ipynb, then rerun it"
        )
        assert got["id"] == want["id"] == f"cell-{i:02d}", "cell ids must be stable"
        assert "execution" not in got.get("metadata", {}), (
            f"cell {i} carries wall-clock timestamps, so every rebuild is a diff"
        )

    assert committed["metadata"]["kernelspec"]["name"] == "python3", (
        "a venv-specific kernel name would not open on anyone else's machine"
    )


def test_the_invariant_unit_tables_match_the_script_that_produces_them():
    """METHODS §5.5 quotes what a splice site with nothing to choose between costs, and
    what excluding it recovers. Both rows are measurements, so both must still be
    reproducible — including the "before the fix" one, which is why the script keeps the
    old arithmetic rather than a remembered number."""
    text = (ROOT / "docs" / "METHODS.md").read_text()
    section = text[text.index("### 5.5 Sites with nothing to choose between") :]
    section = section[: section.index("## 5b.")]

    values = _load("invariant_units").measure(seeds=5)

    rows = re.findall(
        r"^\| (none present|one per unit that can, counted as before|"
        r"one per unit that can, excluded) \| ([\d.]+) \| \*{0,2}([\d.]+)\*{0,2} \|$",
        section,
        re.M,
    )
    assert len(rows) == 3, f"expected three rows in the §5.5 table, parsed {rows}"
    keys = {
        "none present": "absent",
        "one per unit that can, counted as before": "before the fix",
        "one per unit that can, excluded": "invariant present",
    }
    for label, precision, rate in rows:
        key = keys[label]
        assert float(precision) == pytest.approx(
            values[f"estimated precision, invariant {key}" if key == "absent"
                   else f"estimated precision, {key}"], rel=0.05
        ), f"{label}: precision"
        assert float(rate) == pytest.approx(
            values[f"false-positive rate, invariant {key}" if key == "absent"
                   else f"false-positive rate, {key}"], abs=0.01
        ), f"{label}: false-positive rate"

    # the prose numbers, which drift just as readily as a table
    assert "48 %" in section and "1.91" in section
    assert str(int(values["true precision"])) in section


def test_the_readme_counts_what_the_repository_actually_contains():
    """Every number in the README that counts something here, checked against it.

    These rot faster than anything else in the documentation and nothing noticed: the
    README claimed the code was broken 85 ways when the survey held 111 mutations, and
    97 % coverage when the suite reached 98 %. A reader has no way to tell a stale count
    from a current one, and a count that is quietly wrong makes every other number in the
    file worth less.

    The collected-test total comes from `tests/conftest.py`, which records what pytest
    collected, so this costs nothing. It is skipped on a partial run, where the number
    would be meaningless.
    """
    import sys

    # pytest has already imported it; the name depends on the import mode
    conftest = sys.modules.get("conftest") or sys.modules.get("tests.conftest")

    readme = (ROOT / "README.md").read_text()
    survey = _load("mutation_survey")

    mutations = len(survey.MUTATIONS)
    assert f"broken {mutations} ways" in readme, f"the survey holds {mutations} mutations"

    workflow = (ROOT / ".github" / "workflows" / "mutation.yml").read_text()
    assert f"Break the code {mutations} ways" in workflow, "the CI step name too"

    equivalent = sum(1 for m in survey.MUTATIONS if m[4] == "equivalent")
    spelled = {2: "two", 3: "three", 4: "four", 5: "five"}[equivalent]
    assert f"{spelled} documented equivalent mutants" in workflow
    assert survey.__doc__.count("equivalent") >= 1

    properties = sum(
        1 for line in (ROOT / "tests" / "test_properties.py").read_text().splitlines()
        if line.startswith("def test_")
    )
    spelled = {14: "Fourteen", 20: "Twenty", 21: "Twenty-one", 22: "Twenty-two"}.get(properties)
    assert spelled and f"{spelled} of those are **property-based**" in readme, (
        f"test_properties.py holds {properties} tests"
    )

    if conftest is not None and conftest.COLLECTED is not None:
        assert f"# {conftest.COLLECTED} tests:" in readme, (
            f"pytest collected {conftest.COLLECTED}"
        )


def test_every_mutation_still_has_something_to_break():
    """A mutation whose target text has moved is a rule that has quietly stopped being
    checked, and nothing said so for a week.

    `describe` was fixed to handle a missing distance on the `ptc_nmd` branch, which
    rewrote the line one mutation replaces. The survey does report `NOT APPLIED`, and
    that is how it was found — but the survey runs weekly, so the rule went unchecked
    until then. This is the same check, in the suite that runs on every push, and it
    costs a few string searches.

    It also rules out the quieter failure: a target that now matches more than once,
    where the survey would rewrite two places and test something other than what its
    label says.
    """
    survey = _load("mutation_survey")
    src = ROOT / "src" / "splicescope"

    sources = {name: (src / f"{name}.py").read_text() for name in
               sorted({m[0] for m in survey.MUTATIONS})}
    problems = []
    for module, label, old, new, _expected in survey.MUTATIONS:
        found = sources[module].count(old)
        if found != 1:
            problems.append(f"{module}: {label!r} matches {found} places, needs exactly 1")
        if old == new:
            problems.append(f"{module}: {label!r} does not change anything")
    assert not problems, "\n  " + "\n  ".join(problems)

    # and the mutated source still has to be importable, or "caught" means a SyntaxError
    # rather than a test noticing anything
    import ast

    for module, label, old, new, _expected in survey.MUTATIONS:
        mutated = sources[module].replace(old, new)
        try:
            ast.parse(mutated)
        except SyntaxError as exc:  # pragma: no cover - only on a broken mutation
            raise AssertionError(f"{module}: {label!r} does not parse: {exc}") from exc
