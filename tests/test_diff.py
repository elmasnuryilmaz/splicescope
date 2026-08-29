import numpy as np

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
