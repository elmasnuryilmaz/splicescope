import numpy as np
import pandas as pd

from splicescope.quantify import compute_psi, psi_matrix


def test_psi_donor_sums_to_one_per_site():
    # two junctions sharing donor 200 in one sample
    df = pd.DataFrame(
        [
            ("chr1", 200, 300, "+", 30, "s1"),
            ("chr1", 200, 260, "+", 10, "s1"),
        ],
        columns=["chrom", "start", "end", "strand", "count", "sample"],
    )
    out = compute_psi(df, min_reads=1)
    # 30/40 and 10/40
    assert np.isclose(out.loc[0, "psi_donor"], 0.75)
    assert np.isclose(out.loc[1, "psi_donor"], 0.25)
    assert np.isclose(out["psi_donor"].sum(), 1.0)


def test_psi_nan_below_min_reads():
    df = pd.DataFrame(
        [("chr1", 200, 300, "+", 3, "s1")],
        columns=["chrom", "start", "end", "strand", "count", "sample"],
    )
    out = compute_psi(df, min_reads=10)
    assert np.isnan(out.loc[0, "psi_donor"])


def test_psi_matrix_shape():
    df = pd.DataFrame(
        [
            ("chr1", 200, 300, "+", 30, "s1"),
            ("chr1", 200, 300, "+", 20, "s2"),
        ],
        columns=["chrom", "start", "end", "strand", "count", "sample"],
    )
    out = compute_psi(df, min_reads=1)
    m = psi_matrix(out)
    assert m.shape == (1, 2)
    assert set(m.columns) == {"s1", "s2"}


def test_a_junction_listed_twice_for_one_sample_is_refused():
    """Found by re-deriving the beta-binomial's fitted group means by hand, on a unit
    where they disagreed.

    The denominator below sums every row at the splice site, so a junction listed twice
    counts the site twice and Ψ becomes its share of a total that is not the site's. The
    figures still plot, the test still runs, and the number is wrong. A STAR
    `SJ.out.tab` names each junction once per sample, so a repeat is two files read under
    one name or a table concatenated with itself — a mistake worth naming rather than
    averaging away.
    """
    import pytest

    from splicescope.quantify import compute_psi

    shared_site = [
        {"chrom": "chr1", "start": 100, "end": 200, "strand": "+", "sample": "S1",
         "count": 5},
        {"chrom": "chr1", "start": 100, "end": 300, "strand": "+", "sample": "S1",
         "count": 15},
    ]
    clean = pd.DataFrame(shared_site)
    psi = compute_psi(clean, min_reads=1)
    assert psi.loc[psi["end"] == 200, "psi_donor"].iloc[0] == pytest.approx(0.25)

    doubled = pd.DataFrame(shared_site + [shared_site[0]])
    with pytest.raises(ValueError, match="more than once for the same sample") as excinfo:
        compute_psi(doubled, min_reads=1)
    message = str(excinfo.value)
    assert "chr1:100-200+" in message and "'S1'" in message, "name the one it found"
    assert "counted twice" in message, "and what it would have done"


def test_the_simulator_never_writes_a_junction_twice_for_one_sample():
    """It used to, on a quarter of seeds. The deceptive-noise branch runs from an
    intron's known donor to a novel acceptor, which is the shape of a cryptic exon's
    upstream junction, and the two coincided exactly — leaving one junction labelled both
    cryptic and noise, and a splice site whose total was counted twice. Events are written
    before noise, so the real one wins.
    """
    from splicescope.simulate import simulate_dataset

    key = ["chrom", "start", "end", "strand", "sample"]
    collisions = 0
    for seed in range(1, 13):
        observed = simulate_dataset(
            n_genes=40, n_per_group=3, cryptic_fraction=0.7, seed=seed
        ).observed
        sizes = observed.groupby(key, observed=True).size()
        collisions += int((sizes > 1).sum())
    assert collisions == 0, f"{collisions} junctions written twice for one sample"
