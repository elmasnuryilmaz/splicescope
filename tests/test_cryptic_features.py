"""The features the cryptic-junction classifier learns from.

These were almost entirely unchecked: of eight deliberate defects introduced here,
seven passed the whole suite. A feature that is quietly wrong does not make the model
fail, it makes it learn the wrong thing and report a good score for doing so, which is
the failure mode a classifier is worst at revealing.

One hand-built input, with every expected value worked out from the coordinates rather
than copied from the code's output.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from splicescope.cryptic import FEATURE_COLUMNS, extract_features

# one known intron: its donor is at 200 and its acceptor at 399 (plus strand), so the
# known-site index for (chr1, +) is exactly {200, 399}
KNOWN = pd.DataFrame(
    [("chr1", 200, 399, "+", "G1", "GENE1")],
    columns=["chrom", "start", "end", "strand", "gene_id", "gene_name"],
)


def _rows(start, end, sclass, motif, counts, truths, psi=None, chrom="chr1"):
    psi = psi if psi is not None else [np.nan] * len(counts)
    return [
        {
            "chrom": chrom, "start": start, "end": end, "strand": "+",
            "sample": f"S{i}", "count": c, "motif": motif, "sclass": sclass,
            "gene_id": "G1", "psi_donor": p, "is_cryptic_truth": t,
        }
        for i, (c, t, p) in enumerate(zip(counts, truths, psi, strict=True))
    ]


@pytest.fixture
def features():
    rows = []
    # a novel donor 50 nt inside the intron, sharing the known acceptor
    rows += _rows(250, 399, "novel_donor", "GT/AG", [5, 0, 12], [1, 1, 1], [0.2, np.nan, 0.4])
    # a novel acceptor 49 nt before the known one, sharing the known donor,
    # with the minus-strand canonical motif
    rows += _rows(200, 350, "novel_acceptor", "CT/AC", [3, 3, 3], [0, 0, 0])
    # both ends novel, well past the annotation, non-canonical, labelled in one sample
    rows += _rows(500, 600, "cryptic", "non-canonical", [1, 0, 0], [0, 1, 0])
    # a chromosome the annotation says nothing about
    rows += _rows(900, 999, "cryptic", "GT/AG", [4, 4, 4], [0, 0, 0], chrom="chr9")
    feats = extract_features(pd.DataFrame(rows), KNOWN)
    return feats.set_index("start")


def test_the_distance_to_a_known_site_looks_both_ways(features):
    """`searchsorted` gives the insertion point, so the nearest site may be the one
    before it. Searching only forwards turns a 50 nt shift into a 149 nt one, which is
    the difference between a plausible cryptic site and an implausible one."""
    assert features.loc[250, "dist_known_donor"] == 50   # 200 is nearer than 399
    assert features.loc[200, "dist_known_acceptor"] == 49  # 399 is nearer than 200
    assert features.loc[500, "dist_known_donor"] == 101
    assert features.loc[500, "dist_known_acceptor"] == 201


def test_both_splice_sites_of_a_known_intron_count_as_known(features):
    """A known intron contributes two known positions, not one. Indexing only its donor
    leaves every real acceptor looking novel."""
    assert features.loc[250, "dist_known_acceptor"] == 0  # 399 is the known acceptor
    assert features.loc[200, "dist_known_donor"] == 0     # 200 is the known donor


def test_a_chromosome_with_no_annotation_gives_no_distance(features):
    """Not zero, and not the distance to something on another chromosome."""
    assert np.isnan(features.loc[900, "dist_known_donor"])
    assert np.isnan(features.loc[900, "dist_known_acceptor"])


def test_intron_length_is_inclusive_of_both_ends(features):
    assert features.loc[250, "intron_length"] == 150   # 399 - 250 + 1
    assert features.loc[200, "intron_length"] == 151   # 350 - 200 + 1


def test_only_samples_with_reads_count_as_support(features):
    """The measured zeros an aligner reports by omission are filled in upstream, so a
    count of 0 is a real row. Counting it as support makes every junction look
    recurrent, which is the feature that separates events from artefacts."""
    assert features.loc[250, "n_samples_support"] == 2  # 5, 0, 12
    assert features.loc[200, "n_samples_support"] == 3
    assert features.loc[500, "n_samples_support"] == 1


def test_both_canonical_motifs_are_canonical(features):
    """GT/AG on the plus strand reads as CT/AC on the minus. Dropping the second makes
    every minus-strand junction look non-canonical, and the motif is one of the
    strongest features the model has."""
    assert features.loc[250, "canonical_motif"] == 1   # GT/AG
    assert features.loc[200, "canonical_motif"] == 1   # CT/AC
    assert features.loc[500, "canonical_motif"] == 0   # non-canonical


def test_a_junction_is_labelled_cryptic_if_any_of_its_rows_says_so(features):
    """Labels belong to junctions, not to rows, so the aggregation over a junction's
    samples has to be total rather than partial."""
    assert features.loc[500, "is_cryptic_truth"] == 1  # 0, 1, 0
    assert features.loc[250, "is_cryptic_truth"] == 1
    assert features.loc[200, "is_cryptic_truth"] == 0


def test_each_label_lands_on_its_own_junction(features):
    """The labels are aggregated in a second pass over the same grouping and attached by
    position. That is only correct while both passes agree on the order, and a silent
    misalignment would train the model on shuffled labels while every metric still
    looked healthy."""
    assert dict(features["is_cryptic_truth"]) == {200: 0, 250: 1, 500: 1, 900: 0}


def test_mean_psi_skips_the_samples_that_have_none(features):
    assert features.loc[250, "mean_psi_donor"] == pytest.approx(0.3)  # 0.2, -, 0.4
    assert np.isnan(features.loc[200, "mean_psi_donor"])


def test_both_ends_novel_is_its_own_feature(features):
    assert features.loc[500, "is_novel_both"] == 1
    assert features.loc[250, "is_novel_both"] == 0
    assert set(FEATURE_COLUMNS) <= set(features.reset_index().columns)
