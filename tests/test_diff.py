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
