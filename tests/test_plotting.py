import matplotlib
import numpy as np
import pandas as pd
import pytest

from splicescope import plotting


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    matplotlib.pyplot.close("all")


def test_backend_is_headless():
    # importing the module must not require a display (CI, servers)
    assert matplotlib.get_backend().lower() == "agg"


def test_annotation_summary_bar_heights_match_input():
    summary = pd.DataFrame(
        {"sclass": ["annotated", "novel_donor", "novel_both"], "n_junctions": [120, 30, 7]}
    )
    ax = plotting.plot_annotation_summary(summary)
    heights = [p.get_height() for p in ax.patches]
    assert heights == [120, 30, 7]
    assert [t.get_text() for t in ax.get_xticklabels()] == list(summary["sclass"])


def test_volcano_splits_points_on_both_thresholds():
    diff = pd.DataFrame(
        {
            "qvalue": [0.001, 0.001, 0.20, 0.30, np.nan],
            "delta_psi": [0.50, 0.02, 0.60, 0.01, 0.9],
        }
    )
    # significant needs qvalue <= 0.05 AND |delta_psi| >= 0.1 -> only the first row.
    # the NaN row is dropped entirely.
    ax = plotting.plot_volcano(diff, q=0.05, min_delta=0.1)
    ns, sig = ax.collections[0], ax.collections[1]
    assert len(sig.get_offsets()) == 1
    assert len(ns.get_offsets()) == 3
    assert sig.get_offsets()[0][0] == pytest.approx(0.50)


def test_volcano_thresholds_are_configurable():
    diff = pd.DataFrame({"qvalue": [0.04, 0.04], "delta_psi": [0.30, 0.05]})
    lenient = plotting.plot_volcano(diff, q=0.05, min_delta=0.01)
    assert len(lenient.collections[1].get_offsets()) == 2
    strict = plotting.plot_volcano(diff, q=0.05, min_delta=0.5)
    assert len(strict.collections[1].get_offsets()) == 0


def test_volcano_clips_zero_qvalues_instead_of_infinity():
    diff = pd.DataFrame({"qvalue": [0.0], "delta_psi": [0.4]})
    ax = plotting.plot_volcano(diff)
    y = ax.collections[1].get_offsets()[0][1]
    assert np.isfinite(y) and y == pytest.approx(300.0)


def test_event_summary_uses_fixed_order_and_zero_fills():
    events = pd.DataFrame({"event_type": ["SE", "SE", "A3SS"]})
    ax = plotting.plot_event_summary(events)
    assert [t.get_text() for t in ax.get_xticklabels()] == ["SE", "MXE", "A5SS", "A3SS"]
    # MXE and A5SS are absent from the data but must still be drawn as zero
    assert [p.get_height() for p in ax.patches] == [2, 0, 0, 1]


def test_event_volcano_draws_one_series_per_event_type():
    ediff = pd.DataFrame(
        {
            "qvalue": [0.01, 0.01, 0.01],
            "delta_psi": [0.3, -0.4, 0.2],
            "event_type": ["SE", "SE", "A5SS"],
        }
    )
    ax = plotting.plot_event_volcano(ediff)
    assert len(ax.collections) == 4  # SE, MXE, A5SS, A3SS
    counts = [len(c.get_offsets()) for c in ax.collections]
    assert counts == [2, 0, 1, 0]
    assert [t.get_text() for t in ax.get_legend().get_texts()] == ["SE", "MXE", "A5SS", "A3SS"]


def test_enrichment_keeps_only_top_n_and_drops_missing_qvalues():
    enrich = pd.DataFrame(
        {
            "term": [f"set{i}" for i in range(6)],
            "qvalue": [0.001, 0.002, 0.003, 0.004, 0.005, np.nan],
        }
    )
    ax = plotting.plot_enrichment(enrich, top=3)
    assert len(ax.patches) == 3
    # bars are reversed so the most significant term sits at the top of the axis
    assert [t.get_text() for t in ax.get_yticklabels()] == ["set2", "set1", "set0"]


def test_importance_bars_carry_error_bars():
    imp = pd.DataFrame(
        {"feature": ["a", "b"], "importance": [0.4, 0.1], "std": [0.05, 0.02]}
    )
    ax = plotting.plot_importance(imp)
    assert [p.get_width() for p in ax.patches] == [0.1, 0.4]  # reversed for display
    assert ax.containers[0].has_xerr


def test_roc_reports_the_auc_it_computed():
    y_true = np.array([0, 0, 1, 1])
    scores = np.array([0.1, 0.4, 0.35, 0.8])  # AUC = 0.75
    ax = plotting.plot_roc(y_true, scores)
    assert "AUC = 0.75" in ax.get_legend().get_texts()[0].get_text()


def test_savefig_creates_missing_parent_directories(tmp_path):
    fig, ax = matplotlib.pyplot.subplots()
    ax.plot([0, 1], [0, 1])
    out = tmp_path / "figures" / "nested" / "roc.png"
    written = plotting.savefig(fig, out)
    assert written == out
    assert out.exists() and out.stat().st_size > 0
