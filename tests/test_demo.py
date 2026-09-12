"""What the dashboard runs.

`app/streamlit_app.py` is the first thing a reader of this repository tries, and its
analysis had no test: a change to any signature it touched would have appeared as a
traceback on a public page rather than as a red build. The analysis now lives in
`splicescope.demo` and this covers it, including the corners of the sliders the page
exposes — one of which crashed.
"""

from __future__ import annotations

import itertools

import matplotlib
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from splicescope import plotting  # noqa: E402
from splicescope.demo import run_demo  # noqa: E402

#: The slider ranges the dashboard offers, at their extremes.
CORNERS = list(itertools.product((5, 40), (3, 10), (0.1, 1.0), (0.0, 0.3)))


@pytest.fixture(scope="module")
def default_view():
    """The dashboard's own settings, at the classifier size it really uses. Shared,
    because training the full forest is the most expensive thing in this file."""
    return run_demo()


def test_the_default_dashboard_view_has_everything_in_it(default_view):
    result = default_view
    assert not result.annotation_summary.empty
    assert not result.differential.empty
    assert result.classifier_note is None
    assert result.metrics["roc_auc"] > 0.75
    assert result.scores is not None and not result.scores.empty
    assert len(result.importances) == 8


@pytest.mark.parametrize(("n_genes", "n_rep", "fraction", "noise"), CORNERS)
def test_no_slider_position_raises(n_genes, n_rep, fraction, noise):
    """Four of these sixteen combinations used to raise `IndexError: index 1 is out of
    bounds for axis 1 with size 1`, from `predict_proba` on a model fitted to a single
    class. All four had the cryptic fraction at its minimum and the label noise at zero,
    which is a reasonable thing for a reader to ask for."""
    # a small forest: what is being tested is that nothing raises, which does not
    # depend on how well the model is trained
    result = run_demo(n_genes=n_genes, n_per_group=n_rep, cryptic_fraction=fraction,
                      label_noise=noise, n_estimators=20)
    assert not result.annotation_summary.empty
    assert not result.differential.empty
    # either the classifier ran, or it said why it did not — never both, never neither
    assert (result.metrics is None) == (result.classifier_note is not None)
    if result.classifier_note is not None:
        assert result.scores is None and result.importances is None
        assert "class" in result.classifier_note


def test_a_one_class_dataset_explains_itself_instead_of_crashing():
    result = run_demo(n_genes=5, n_per_group=3, cryptic_fraction=0.1, label_noise=0.0,
                      n_estimators=20)
    assert result.metrics is None
    assert "no two classes to separate" in result.classifier_note
    assert "Raise the fraction" in result.classifier_note
    # the rest of the page still has content to show
    assert len(result.differential) > 0
    assert result.annotation_summary["n_junctions"].sum() > 0


def test_every_figure_the_dashboard_draws_renders(default_view):
    """The page draws three. A plotting signature change would otherwise surface as a
    traceback in the deployed app."""
    result = default_view
    for draw in (
        lambda ax: plotting.plot_annotation_summary(result.annotation_summary, ax=ax),
        lambda ax: plotting.plot_volcano(result.differential, q=0.05, min_delta=0.1, ax=ax),
        lambda ax: plotting.plot_importance(result.importances, ax=ax),
    ):
        fig, ax = plt.subplots()
        draw(ax)
        plt.close(fig)


def test_the_classifier_refuses_a_single_class_rather_than_indexing_past_it():
    """Directly, at the level the bug lived: `np.bincount` cannot see a class that is
    absent, so the fold-count guard passed and the forest fitted a one-class model."""
    import numpy as np
    import pandas as pd

    from splicescope.cryptic import FEATURE_COLUMNS
    from splicescope.ml import CrypticClassifier

    feats = pd.DataFrame({c: np.arange(20, dtype=float) for c in FEATURE_COLUMNS})
    feats["is_cryptic_truth"] = 0
    with pytest.raises(ValueError, match="no two classes"):
        CrypticClassifier(n_estimators=10).fit(feats)
    with pytest.raises(ValueError, match="no two classes"):
        CrypticClassifier(n_estimators=10).evaluate(feats)

    feats["is_cryptic_truth"] = 1  # the all-positive case was caught only by accident
    with pytest.raises(ValueError, match="no two classes"):
        CrypticClassifier(n_estimators=10).fit(feats)
