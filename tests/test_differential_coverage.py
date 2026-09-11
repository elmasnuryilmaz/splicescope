"""Which measurements reach the differential test, and which are correctly withheld.

Two bugs lived here and neither was visible to the suite: the ``min_reads`` coverage
filter was inert on the default code path, and a junction an aligner omitted because it
had no reads was treated as missing data rather than as a measured zero. Both are about
the boundary between "we did not measure this" and "we measured this to be zero", so
they are tested together.
"""

import numpy as np
import pandas as pd

from splicescope.diff import _wide, differential_splicing
from splicescope.quantify import add_unobserved_zeros, compute_psi
from splicescope.simulate import simulate_dataset

SAMPLES = ["C1", "C2", "C3", "K1", "K2", "K3"]
GROUPS = {s: ("ctrl" if s.startswith("C") else "kd") for s in SAMPLES}


def _junction(start, end, sample, count, chrom="chr1", strand="+"):
    return dict(
        chrom=chrom, start=start, end=end, strand=strand, sample=sample, count=count
    )


# --------------------------------------------------------------------------------
# the min_reads coverage filter
# --------------------------------------------------------------------------------


def test_all_nan_psi_survives_the_pivot_as_nan():
    """pandas sums an all-NaN group to 0.0, which silently turned "uninformative"
    into "zero usage" and erased the coverage filter."""
    df = pd.DataFrame(
        {"k": ["j1", "j1"], "sample": ["s1", "s2"], "psi": [np.nan, np.nan]}
    )
    wide = _wide(df, ["k"], "psi", ["s1", "s2"], aggfunc="mean")
    assert wide.isna().all().all()


def test_raising_min_reads_removes_units_from_the_differential_test():
    """The regression: min_reads had no effect whatsoever on the default test."""
    ds = simulate_dataset(n_genes=20, n_per_group=4, seed=3)
    counts = []
    for min_reads in (10, 200, 100_000):
        psi = compute_psi(ds.observed, min_reads=min_reads)
        counts.append(len(differential_splicing(psi, ds.groups)))
    assert counts[0] > counts[1] > counts[2]
    assert counts[2] == 0, "no unit can be informative at a 100k-read threshold"


def test_a_unit_below_the_coverage_threshold_is_never_called_significant():
    """Three reads per sample must not produce a genome-wide-significant call."""
    rows = []
    for s in SAMPLES:
        kd = s.startswith("K")
        rows.append(_junction(200, 300, s, 3 if kd else 0))
        rows.append(_junction(200, 260, s, 0 if kd else 3))
    df = pd.DataFrame([r for r in rows if r["count"] > 0])
    psi = compute_psi(df, min_reads=10)
    assert psi["psi_donor"].isna().all()
    assert differential_splicing(psi, GROUPS).empty


# --------------------------------------------------------------------------------
# junctions an aligner omits because they carry no reads
# --------------------------------------------------------------------------------


def test_a_junction_absent_from_one_group_is_tested_as_zero_usage():
    """The cleanest cryptic event there is — absent in every control, present in every
    knockdown — was dropped from the differential table entirely."""
    rows = []
    for s in SAMPLES:
        rows.append(_junction(1000, 2000, s, 100))  # annotated partner, always present
        if s.startswith("K"):
            rows.append(_junction(1000, 1500, s, 100))  # cryptic, KD only
    res = differential_splicing(compute_psi(pd.DataFrame(rows)), GROUPS)

    cryptic = res[res["end"] == 1500]
    assert len(cryptic) == 1, "the cryptic junction must be tested, not dropped"
    row = cryptic.iloc[0]
    assert row["mean_ctrl"] == 0.0
    assert row["delta_psi"] > 0.4
    assert row["qvalue"] < 0.01


def test_zeros_are_filled_only_where_the_splice_site_was_actually_covered():
    """A site nobody sequenced stays unobserved; inventing a zero there would claim a
    measurement that was never made."""
    rows = [
        _junction(1000, 2000, "s1", 50),  # donor 1000 covered in s1 only
        _junction(1000, 1500, "s1", 50),
        _junction(5000, 6000, "s2", 50),  # s2 covers an unrelated locus
    ]
    psi = compute_psi(pd.DataFrame(rows))
    filled = set(zip(psi["start"], psi["end"], psi["sample"], strict=False))
    # s1 saw both junctions at donor 1000; nothing to fill there
    assert (5000, 6000, "s1") not in filled, "donor 5000 has no coverage in s1"
    assert (1000, 2000, "s2") not in filled, "donor 1000 has no coverage in s2"
    assert (1000, 1500, "s2") not in filled


def test_filled_zero_rows_carry_the_junctions_own_annotation():
    """A filled row that lost sclass or gene_id would corrupt every downstream filter."""
    rows = []
    for s in ("s1", "s2"):
        rows.append({**_junction(1000, 2000, s, 40), "motif": "GT/AG", "sclass": "annotated",
                     "gene_id": "G1"})
    rows.append({**_junction(1000, 1500, "s1", 10), "motif": "GT/AG",
                 "sclass": "novel_acceptor", "gene_id": "G1"})
    out = add_unobserved_zeros(compute_psi(pd.DataFrame(rows), fill_unobserved=False))

    added = out[(out["end"] == 1500) & (out["sample"] == "s2")]
    assert len(added) == 1
    assert added.iloc[0]["count"] == 0
    assert added.iloc[0]["sclass"] == "novel_acceptor"
    assert added.iloc[0]["gene_id"] == "G1"
    assert added.iloc[0]["motif"] == "GT/AG"


def test_filling_can_be_turned_off_to_reproduce_the_old_behaviour():
    rows = [_junction(1000, 2000, s, 100) for s in SAMPLES]
    rows += [_junction(1000, 1500, s, 100) for s in SAMPLES if s.startswith("K")]
    df = pd.DataFrame(rows)
    assert len(compute_psi(df, fill_unobserved=False)) == len(df)
    assert len(compute_psi(df, fill_unobserved=True)) == len(SAMPLES) * 2


def test_junction_and_event_levels_agree_about_an_absent_junction():
    """events.py already defaulted a missing count to zero; the junction level did not,
    so the two halves of the same tool disagreed about the same data."""
    rows = []
    for s in SAMPLES:
        rows.append(_junction(1000, 2000, s, 100))
        if s.startswith("K"):
            rows.append(_junction(1000, 1500, s, 100))
    psi = compute_psi(pd.DataFrame(rows))
    tested = psi[psi["sample"] == "C1"]
    assert set(tested["end"]) == {1500, 2000}
    assert tested.loc[tested["end"] == 1500, "psi_donor"].iloc[0] == 0.0
