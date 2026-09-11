import numpy as np
import pandas as pd

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
