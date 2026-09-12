"""Which measurements reach the differential test, and which are correctly withheld.

Two bugs lived here and neither was visible to the suite: the ``min_reads`` coverage
filter was inert on the default code path, and a junction an aligner omitted because it
had no reads was treated as missing data rather than as a measured zero. Both are about
the boundary between "we did not measure this" and "we measured this to be zero", so
they are tested together.
"""

import warnings

import numpy as np
import pandas as pd
import pytest

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
        with warnings.catch_warnings():
            # the 100k threshold leaves nothing testable, which the function now says;
            # this test is about the counts falling, not about that message
            warnings.simplefilter("ignore", UserWarning)
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
    # and the caller is told that nothing was tested, rather than left to read an empty
    # table as "no difference"
    with pytest.warns(UserWarning, match="nothing was tested"):
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


def _small_psi():
    rows = [
        _junction(100, 199, s, c)
        for s, c in zip(SAMPLES, [50, 48, 52, 20, 22, 18], strict=True)
    ] + [
        _junction(100, 299, s, c)
        for s, c in zip(SAMPLES, [10, 12, 8, 40, 38, 42], strict=True)
    ]
    return compute_psi(pd.DataFrame(rows))


def test_groups_that_name_no_sample_in_the_table_is_an_error_not_an_empty_result():
    """The worst outcome this tool can produce is a clean run and a wrong conclusion.
    With none of the sample names matching, every row is unassigned, every unit fails
    min_samples, and an empty table comes back — so the caller reports "0 significant
    junctions" and believes it. The CLI already refused this from the filenames; the
    library returned the empty table."""
    import pytest

    psi = _small_psi()
    renamed = {f"other_{s}": g for s, g in GROUPS.items()}
    with pytest.raises(ValueError, match="no sample in the Psi table is named"):
        differential_splicing(psi, renamed)


def test_a_partial_sample_mismatch_warns_and_keeps_going():
    """Half a mismatch is usually a real experiment with one sample dropped, so it must
    not be fatal — but it must not be silent either."""
    import warnings

    psi = _small_psi()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        out = differential_splicing(psi, {**GROUPS, "K4": "kd"})
    assert not out.empty
    assert any("K4" in str(w.message) for w in caught)


def test_a_misspelled_value_column_says_which_ones_exist():
    import pytest

    psi = _small_psi()
    with pytest.raises(ValueError, match=r"psi_donor"):
        differential_splicing(psi, GROUPS, value="psi")


def test_a_grouping_key_the_table_does_not_have_is_named_in_the_error():
    """Event-level tests are run by passing key=["event_id"]; handing that to a
    junction-level table used to fail as a KeyError from inside pandas."""
    import pytest

    psi = _small_psi()
    with pytest.raises(ValueError, match=r"event_id"):
        differential_splicing(psi, GROUPS, key=["event_id"])


def test_a_site_with_exactly_min_reads_is_kept_not_dropped():
    """`min_reads` is the coverage a splice site needs before Psi means anything, so a
    site carrying exactly that many reads is measured, not withheld. The boundary was
    not pinned, and this filter has been wrong before — it was inert on the default code
    path until 0.9.0."""
    rows = [_junction(100, 199, s, c) for s, c in zip(SAMPLES, [6, 6, 6, 6, 6, 6], strict=True)]
    rows += [_junction(100, 299, s, c) for s, c in zip(SAMPLES, [4, 4, 4, 4, 4, 4], strict=True)]
    obs = pd.DataFrame(rows)  # every donor site totals exactly 10

    at_threshold = compute_psi(obs, min_reads=10)
    assert at_threshold["psi_donor"].notna().all(), "10 reads is not fewer than 10"
    assert at_threshold["donor_total"].eq(10).all()

    just_above = compute_psi(obs, min_reads=11)
    assert just_above["psi_donor"].isna().all(), "11 is more than the site has"


def test_a_zero_is_filled_from_whichever_site_was_covered():
    """The rule is that a junction gets its measured zero where *either* of its splice
    sites has reads in that sample, because either one gives Psi a denominator. Filling
    only from the donor loses every junction whose acceptor was the covered end — half
    the cases, and on the minus strand the other half."""
    shared_acceptor = [
        _junction(100, 300, "C1", 30),   # only C1 uses this junction
        _junction(200, 300, "C1", 10),   # shares its acceptor, 300
        _junction(200, 300, "K1", 40),   # K1 covers the acceptor but nothing at donor 100
    ]
    filled = add_unobserved_zeros(compute_psi(pd.DataFrame(shared_acceptor), fill_unobserved=False))
    added = filled[(filled["start"] == 100) & (filled["sample"] == "K1")]
    assert len(added) == 1, "the junction's acceptor had reads in K1, so its zero is measurable"
    assert added["count"].iloc[0] == 0

    psi = compute_psi(pd.DataFrame(shared_acceptor))
    row = psi[(psi["start"] == 100) & (psi["sample"] == "K1")]
    assert row["psi_acceptor"].iloc[0] == 0.0
    assert row["acceptor_total"].iloc[0] == 40.0


def test_a_psi_matrix_instead_of_a_long_table_says_so():
    """`psi_matrix` returns junctions by samples, which is the shape a person reaches
    for — and with no `sample` column nothing can be assigned to a condition."""
    import pytest

    from splicescope.quantify import psi_matrix

    psi = _small_psi()
    wide = psi_matrix(psi)
    with pytest.raises(ValueError, match="no 'sample' column"):
        differential_splicing(wide.reset_index(), GROUPS)


def test_a_threshold_nothing_can_meet_says_so_instead_of_reporting_no_difference():
    """`--min-reads 100000` printed "0 significant junctions" and stopped. That reads as
    "no differential splicing", which is a conclusion; the truth is that not one unit
    was testable. Ψ is NaN wherever a splice site is under the threshold, so raising it
    far enough empties the table with no other sign."""
    import pytest

    psi = _small_psi()
    with pytest.warns(UserWarning, match="nothing was tested"):
        out = differential_splicing(compute_psi(psi, min_reads=10**6), GROUPS)
    assert out.empty
    assert list(out.columns), "and it still comes back with the schema of a result"


def test_a_threshold_that_leaves_units_testable_says_nothing():
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        out = differential_splicing(_small_psi(), GROUPS)
    assert not out.empty
