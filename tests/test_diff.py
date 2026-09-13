import numpy as np
import pandas as pd
import pytest

from splicescope.diff import benjamini_hochberg, differential_splicing, significant
from splicescope.quantify import compute_psi
from splicescope.simulate import simulate_dataset


def test_bh_monotone_and_bounds():
    p = np.array([0.001, 0.01, 0.5, 0.7, np.nan])
    q = benjamini_hochberg(p)
    finite = q[~np.isnan(q)]
    assert np.all(finite >= 0) and np.all(finite <= 1)
    # q-values are >= their p-values
    assert np.all(finite >= p[~np.isnan(p)] - 1e-9)
    assert np.isnan(q[-1])


def test_differential_detects_injected_cryptic():
    ds = simulate_dataset(n_genes=10, n_per_group=5, cryptic_fraction=1.0, seed=1)
    psi = compute_psi(
        ds.observed.assign(sclass="x"), min_reads=5
    )  # sclass placeholder; not needed here
    diff = differential_splicing(psi, ds.groups)
    assert not diff.empty
    hits = significant(diff, q=0.1, min_delta=0.05)
    # at least one true cryptic junction should surface as differential
    truth_juncs = set(
        map(
            tuple,
            ds.observed.loc[ds.observed["is_cryptic_truth"] == 1, ["start", "end"]]
            .drop_duplicates()
            .to_numpy(),
        )
    )
    hit_juncs = set(map(tuple, hits[["start", "end"]].drop_duplicates().to_numpy()))
    assert truth_juncs & hit_juncs, "no injected cryptic event recovered as differential"


def _one_junction(samples_and_counts, depth=1000):
    rows = []
    for sample, inc in samples_and_counts.items():
        rows.append(dict(chrom="chr1", start=1000, end=1500, strand="+",
                         sample=sample, count=inc))
        rows.append(dict(chrom="chr1", start=1000, end=2000, strand="+",
                         sample=sample, count=depth - inc))
    return pd.DataFrame(rows)


def test_unestimable_dispersion_is_warned_about_not_assumed_away():
    """Falling back to max_precision asserts *no* overdispersion rather than *unknown*
    dispersion, which narrows the test to a binomial one. A 1-vs-1 comparison then calls
    Psi 0.300 against 0.360 significant on no replication at all."""
    import warnings

    psi = compute_psi(_one_junction({"ctrl1": 300, "kd1": 360}))
    groups = {"ctrl1": "ctrl", "kd1": "kd"}

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = differential_splicing(psi, groups, min_samples=1)
    messages = [str(w.message) for w in caught]
    assert any("dispersion could not be estimated" in m for m in messages), messages
    # the number is still produced — the point is that the user is told not to trust it
    assert not result.empty


def test_a_replicated_design_does_not_warn():
    import warnings

    counts = {"C1": 300, "C2": 310, "C3": 295, "K1": 360, "K2": 355, "K3": 370}
    psi = compute_psi(_one_junction(counts))
    groups = {s: ("ctrl" if s.startswith("C") else "kd") for s in counts}

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        differential_splicing(psi, groups)
    assert [str(w.message) for w in caught] == []


def test_the_reported_qvalue_is_benjamini_hochberg_of_the_pvalue():
    """The FDR claim in every output table. `benjamini_hochberg` is tested on its own,
    but nothing checked that `differential_splicing` actually applies it — handing the
    raw p-values through as q-values passed the whole suite."""
    import numpy as np

    from splicescope.diff import benjamini_hochberg, differential_splicing

    rng = np.random.default_rng(4)
    n_units, samples = 200, ["C1", "C2", "C3", "K1", "K2", "K3"]
    rows = []
    for u in range(n_units):
        shift = 0.3 if u < 20 else 0.0  # a handful of real differences
        for s in samples:
            base = 0.35 + (shift if s.startswith("K") else 0.0)
            total = int(rng.integers(60, 140))
            rows.append(
                dict(
                    chrom="chr1", start=1000 + 10 * u, end=1200 + 10 * u, strand="+",
                    sample=s, psi_donor=base, donor_total=float(total),
                    count=round(base * total),
                )
            )
    psi = pd.DataFrame(rows)
    out = differential_splicing(psi, {s: ("ctrl" if s[0] == "C" else "kd") for s in samples})

    assert len(out) == n_units
    expected = benjamini_hochberg(out["pvalue"].to_numpy())
    np.testing.assert_allclose(out["qvalue"].to_numpy(), expected, rtol=1e-12)
    assert (out["qvalue"] >= out["pvalue"] - 1e-12).all(), "a q-value cannot beat its p-value"
    assert not np.allclose(out["qvalue"], out["pvalue"]), (
        "with 200 units the correction has to move something"
    )


def _psi_without_counts(n_units=40, delta=0.35, seed=0, tied=False):
    """A long Ψ table carrying no read counts, which is what forces the rank test.

    ``tied`` saturates the knockdown group at 1.0, which is what a cryptic junction
    actually looks like — off in every control, on in every knockdown — and which sends
    `scipy.stats.mannwhitneyu` down its normal-approximation path.
    """
    rng = np.random.default_rng(seed)
    samples = ["C1", "C2", "C3", "K1", "K2", "K3"]
    rows = []
    for unit in range(n_units):
        base = 0.20 + 0.005 * unit
        for sample in samples:
            knocked = sample.startswith("K")
            if tied:
                psi = 1.0 if knocked else 0.0
            else:
                psi = base + (delta if knocked else 0.0) + rng.normal(0, 0.005)
            rows.append(
                {
                    "chrom": "chr1", "start": 1000 + 10 * unit, "end": 1200 + 10 * unit,
                    "strand": "+", "sample": sample, "psi_donor": psi,
                }
            )
    return pd.DataFrame(rows)


def _groups(psi):
    return {s: ("ctrl" if s[0] == "C" else "kd") for s in psi["sample"].unique()}


def test_the_rank_test_cannot_beat_the_exact_floor():
    """The legacy path, and the whole argument of METHODS §5.1: a two-sided
    Mann-Whitney on n replicates per group cannot return a p-value below 2/C(2n,n),
    which is 0.1 for 3 against 3. Nothing tested the implementation — only the formula
    for the floor, in `test_betabinom.py`."""
    from splicescope.betabinom import min_achievable_rank_pvalue
    from splicescope.diff import differential_splicing

    psi = _psi_without_counts()
    out = differential_splicing(psi, _groups(psi), test="ranksum")

    assert len(out) == 40
    assert out["n_a"].eq(3).all() and out["n_b"].eq(3).all()
    assert out["delta_psi"].min() > 0.3, "the difference is real and one-directional"

    floor = min_achievable_rank_pvalue(3, 3)
    assert floor == pytest.approx(0.1)
    assert out["pvalue"].eq(floor).all(), (
        "separated, untied groups of three give exactly the floor and never less"
    )
    assert out["qvalue"].min() > 0.05, "so nothing survives correction, at any effect size"


def test_ties_take_the_rank_test_below_its_exact_floor_and_correction_holds_anyway():
    """The floor is the *exact* test's. scipy computes it exactly only for a small,
    tie-free sample; with ties it uses the normal approximation, which returns less. A
    junction that is off in every control and on in every knockdown is fully tied on
    both sides and gives 0.047 — under the nominal 0.05, and the signature this tool
    exists to find. So the rank test is not floored at 0.1 in practice.

    What keeps it from calling is Benjamini-Hochberg across everything else tested. That
    is a weaker guarantee than the floor, and it is worth having in a test: 40 switching
    junctions among 2 000 need the 40th at p <= 0.001, and 0.047 is nowhere near it.
    """
    from splicescope.diff import differential_splicing

    rng = np.random.default_rng(3)
    switching, tested = 40, 2000
    samples = ["C1", "C2", "C3", "K1", "K2", "K3"]
    rows = []
    for unit in range(tested):
        for sample in samples:
            knocked = sample.startswith("K")
            psi = (1.0 if knocked else 0.0) if unit < switching else rng.uniform(0.2, 0.8)
            rows.append(
                {
                    "chrom": "chr1", "start": 1000 + 10 * unit, "end": 1200 + 10 * unit,
                    "strand": "+", "sample": sample, "psi_donor": psi,
                }
            )
    psi_table = pd.DataFrame(rows)
    out = differential_splicing(psi_table, _groups(psi_table), test="ranksum")
    switched = out[out["delta_psi"] == 1.0]

    assert len(switched) == switching
    assert np.allclose(switched["pvalue"], 0.046854, atol=1e-5), (
        "a fully tied 3-vs-3 goes below both the exact floor and the nominal 0.05"
    )
    assert switched["qvalue"].min() > 0.05, (
        "correction is what refuses it, and this is the guarantee that actually holds"
    )


def test_auto_falls_back_to_the_rank_test_only_when_counts_are_absent():
    from splicescope.diff import differential_splicing

    psi = _psi_without_counts()
    without = differential_splicing(psi, _groups(psi), test="auto")
    assert "lrt_statistic" not in without.columns, "no counts, so no likelihood ratio"

    with_counts = psi.assign(donor_total=80.0, count=(psi["psi_donor"] * 80).round())
    chosen = differential_splicing(with_counts, _groups(psi), test="auto")
    assert "lrt_statistic" in chosen.columns, "counts present, so the beta-binomial is used"
    assert chosen["pvalue"].min() < 1e-6, "which clears correction where the rank test cannot"
    assert chosen["qvalue"].min() < 0.05


def test_the_beta_binomial_cannot_be_asked_for_without_counts():
    from splicescope.diff import differential_splicing

    psi = _psi_without_counts()
    with pytest.raises(ValueError, match="count"):
        differential_splicing(psi, _groups(psi), test="betabinom")


def _counts_table(n_units=60, delta=0.15, seed=0):
    """A Ψ table with the counts it was computed from, so the beta-binomial can run."""
    rng = np.random.default_rng(seed)
    samples = ["C1", "C2", "C3", "K1", "K2", "K3"]
    rows = []
    for unit in range(n_units):
        for sample in samples:
            psi = 0.35 + (delta if sample.startswith("K") else 0.0) + rng.normal(0, 0.04)
            total = int(rng.integers(60, 120))
            rows.append(
                {
                    "chrom": "chr1", "start": 1000 + 10 * unit, "end": 1200 + 10 * unit,
                    "strand": "+", "gene_id": f"G{unit // 10}", "sample": sample,
                    "psi_donor": psi, "donor_total": float(total),
                    "count": round(psi * total),
                }
            )
    return pd.DataFrame(rows)


def test_the_per_unit_floor_option_reaches_the_test_and_only_widens_the_null():
    """METHODS §5.4 is an argument for a public option, and the option's own branch in
    `differential_splicing` had no test — only the estimator underneath it. Taking the
    smaller of the shared and per-unit precisions can only widen the null, so every
    p-value must move up or stay, and the reported precision must not exceed the shared
    one for any unit."""
    from splicescope.diff import differential_splicing

    psi = _counts_table()
    groups = _groups(psi)
    shared = differential_splicing(psi, groups, dispersion="shared").set_index("start")
    floored = differential_splicing(psi, groups, dispersion="per_unit_floor").set_index("start")

    assert len(floored) == len(shared) == 60
    assert shared["precision"].nunique() == 1, "one precision for every unit"
    assert floored["precision"].nunique() > 1, "per unit, so they differ"
    assert (floored["precision"] <= shared["precision"].iloc[0] + 1e-9).all()

    common = shared.index
    assert (floored.loc[common, "pvalue"] >= shared.loc[common, "pvalue"] - 1e-9).all(), (
        "a wider null cannot make any unit more significant"
    )
    # ΔΨ is a difference of *fitted* group means, and fit_mu depends on the precision it
    # is given, so shrinking the precision moves the effect size a little. It must not
    # move much, and it must never change direction.
    moved = (floored.loc[common, "delta_psi"] - shared.loc[common, "delta_psi"]).abs()
    assert moved.max() < 0.01, f"ΔΨ moved by {moved.max():.4f}"
    assert (
        np.sign(floored.loc[common, "delta_psi"]) == np.sign(shared.loc[common, "delta_psi"])
    ).all()


def test_an_unknown_dispersion_is_refused():
    from splicescope.diff import differential_splicing

    psi = _counts_table(n_units=4)
    with pytest.raises(ValueError, match="unknown dispersion"):
        differential_splicing(psi, _groups(psi), dispersion="empirical_bayes")


def test_a_unit_without_enough_replicates_is_skipped_by_the_rank_test():
    """`min_samples` is per group, and a unit that cannot meet it in both is dropped
    rather than tested on one side."""
    from splicescope.diff import differential_splicing

    psi = _psi_without_counts(n_units=3)
    # leave the second unit with a single knockdown sample
    thin = psi[~((psi["start"] == 1010) & psi["sample"].isin(["K2", "K3"]))]
    out = differential_splicing(thin, _groups(psi), test="ranksum", min_samples=2)

    assert set(out["start"]) == {1000, 1020}, "the thinned unit is absent, not half-tested"
    assert out["n_b"].eq(3).all()


def test_the_rank_test_carries_the_annotation_columns_through():
    """`gene_id` and `sclass` are what make a result table usable, and the rank test
    copies them from the unit's first row."""
    from splicescope.diff import differential_splicing

    psi = _counts_table(n_units=12).drop(columns=["count", "donor_total"])
    psi = psi.assign(sclass="novel_donor")
    out = differential_splicing(psi, _groups(psi), test="ranksum")
    assert {"gene_id", "sclass"} <= set(out.columns)
    assert out["gene_id"].nunique() == 2 and out["sclass"].eq("novel_donor").all()


@pytest.mark.parametrize("test", ["betabinom", "ranksum"])
def test_a_result_with_no_rows_still_has_the_columns_a_result_has(test):
    """Found by a property test, which could not even index the empty frame it got back.

    `differential_splicing` returned a bare `pd.DataFrame()` when no unit met
    `min_samples`, so the empty case had a different shape from every other: selecting
    columns raised KeyError where the same code worked on a result with rows, and
    `to_csv` wrote a file with no header line at all — which is what the CLI produced for
    a run that found nothing.
    """
    from splicescope.diff import differential_splicing, significant

    psi = _counts_table(n_units=6)
    if test == "ranksum":
        psi = psi.drop(columns=["count", "donor_total"])
    groups = _groups(psi)

    populated = differential_splicing(psi, groups, test=test)
    with pytest.warns(UserWarning, match="nothing was tested"):
        empty = differential_splicing(psi, groups, test=test, min_samples=99)

    assert len(populated) == 6 and empty.empty
    assert list(empty.columns) == list(populated.columns), "same schema, same order"
    # the things a caller does next, which used to raise on this frame
    assert empty[["chrom", "start", "delta_psi", "qvalue"]].empty
    def header(frame):
        return frame.to_csv(index=False).splitlines()[0]

    assert header(empty) == header(populated)
    assert significant(empty).empty


def test_a_unit_that_cannot_vary_is_reported_as_untestable_and_can_be_set_aside():
    """A splice site with nothing to choose between gets p = 1 from either test.

    A constitutive donor carrying a single junction has Ψ = 1 in every sample, because
    the site has nothing else to splice to. The rank test has no ordering to work with
    and the beta-binomial's two groups have the same pooled proportion, so both return
    exactly 1.0 whichever way the samples are labelled. It cannot be rejected at any
    threshold — yet it still enters the Benjamini-Hochberg denominator, where it makes
    every real unit's q-value worse. `filter_invariant` sets such units aside; it is off
    by default because it changes every q-value in a run.
    """
    from splicescope.diff import differential_splicing, significant

    psi = _counts_table(n_units=4).drop(columns=["count", "donor_total"])
    flat = psi["start"] == 1000
    psi.loc[flat, "psi_donor"] = 0.5   # the same Ψ in all six samples
    groups = _groups(psi)

    kept = differential_splicing(psi, groups, test="ranksum")
    stuck = kept[kept["start"] == 1000].iloc[0]
    assert len(kept) == 4 and stuck["delta_psi"] == 0.0
    # Exactly 1.0 on most SciPy versions, missing on those that refuse a rank test with
    # no ordering to work with. Both say the same thing and the difference is not this
    # package's to fix, so the invariant is asserted rather than the value — which is
    # what CI said when 3.12 and 3.13 disagreed with everything else.
    assert pd.isna(stuck["pvalue"]) or stuck["pvalue"] == 1.0, (
        f"no evidence either way, got {stuck['pvalue']}"
    )
    # q = 0.5 is already absurdly lax, and it is still not called at it
    assert stuck["start"] not in set(significant(kept, q=0.5, min_delta=0.0)["start"])

    # a unit taking exactly two values — one per group — is the clearest real difference
    # there is, and the filter must never reach it
    two_valued = psi["start"] == 1010
    psi.loc[two_valued, "psi_donor"] = np.where(
        psi.loc[two_valued, "sample"].str.startswith("C"), 0.4, 0.7
    )
    assert psi.loc[two_valued, "psi_donor"].nunique() == 2

    with pytest.warns(UserWarning, match="set aside"):
        dropped = differential_splicing(psi, groups, test="ranksum", filter_invariant=True)

    kept = differential_splicing(psi, groups, test="ranksum")
    assert len(dropped) == 3 and 1000 not in set(dropped["start"])
    assert 1010 in set(dropped["start"]), "two values is variation, not invariance"
    # the units that could be tested keep their p-values; only the correction moves
    merged = kept.merge(dropped, on=["chrom", "start", "end", "strand"], suffixes=("", "_f"))
    assert len(merged) == 3
    assert np.allclose(merged["pvalue"], merged["pvalue_f"], equal_nan=True)
    assert (merged["qvalue_f"] <= merged["qvalue"] + 1e-12).all(), "never worse"
    if stuck["pvalue"] == 1.0:
        # There is only something to gain where the unit reached the correction at all.
        # On a SciPy that reports nothing for it, Benjamini-Hochberg already leaves it
        # out of the denominator and the filter has nothing left to remove — which is
        # not true of the beta-binomial, whose likelihood ratio is a real 1.0.
        assert merged["qvalue_f"].min() < merged["qvalue"].min(), "better where it counts"


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"test": "wilcoxon"}, "unknown test"),
        ({"dispersion": "empirical"}, "unknown dispersion"),
    ],
)
def test_an_option_that_does_not_exist_is_refused_by_name(kwargs, message):
    """Silently falling back to the default would run a different analysis from the one
    that was asked for, and say nothing about it."""
    from splicescope.diff import differential_splicing

    psi = _counts_table(n_units=2)
    with pytest.raises(ValueError, match=message):
        differential_splicing(psi, _groups(psi), **kwargs)


def test_more_than_two_conditions_is_refused_and_the_message_names_them():
    """A three-arm design — control and two doses, or a timecourse — is a real experiment
    and a natural thing to point this at. It is not what a two-group test computes, so it
    is refused, and the message lists what was found so the fix is obvious."""
    from splicescope.diff import differential_splicing

    psi = _counts_table(n_units=2)
    groups = {s: {"C1": "ctrl", "C2": "ctrl", "C3": "low"}.get(s, "high") for s in
              psi["sample"].unique()}

    with pytest.raises(ValueError, match="expected exactly 2 conditions") as excinfo:
        differential_splicing(psi, groups)
    for condition in ("ctrl", "low", "high"):
        assert condition in str(excinfo.value), "say which conditions were found"


def test_the_invariant_filter_leaves_a_table_with_nothing_to_drop_alone():
    """The ordinary case for a table of alternative splice sites, where every unit varies.
    The filter must then be a no-op in both senses: the same rows, and no warning — one
    that fired on a run it changed nothing about would be noise, and noise is how people
    learn to ignore warnings."""
    import warnings

    from splicescope.diff import _drop_invariant_units, differential_splicing

    psi = _counts_table(n_units=5).drop(columns=["count", "donor_total"])
    key = ["chrom", "start", "end", "strand"]
    kept, set_aside = _drop_invariant_units(psi, key, "psi_donor")
    assert set_aside == 0 and kept is psi, "not even a copy is made"

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        filtered = differential_splicing(psi, _groups(psi), test="ranksum",
                                         filter_invariant=True)
    plain = differential_splicing(psi, _groups(psi), test="ranksum")
    assert len(filtered) == len(plain) == 5
    assert np.allclose(filtered["qvalue"], plain["qvalue"])
