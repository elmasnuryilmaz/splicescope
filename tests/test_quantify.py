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
