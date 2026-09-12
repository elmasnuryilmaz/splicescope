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

ROOT = Path(__file__).resolve().parent.parent

#: One row of the METHODS §7b table: name, rate at p<=0.05, rate at p<=1e-6, median p.
_TABLE_ROW = re.compile(
    r"^\| (plain hypergeometric|opportunity-weighted) \| \*\*([\d.]+) %\*\* "
    r"\| ([\d.]+) % \| ([\d.]+) \|$",
    re.M,
)


def _load(script: str):
    """Import a `validation/` script by path — they are programs, not a package."""
    path = ROOT / "validation" / f"{script}.py"
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
