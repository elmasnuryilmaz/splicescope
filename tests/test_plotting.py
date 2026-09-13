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


def test_savefig_releases_the_figure_so_pipelines_do_not_leak(tmp_path):
    """pyplot keeps every figure alive until closed; `run` writes up to eight per
    invocation and kept them all, which matplotlib warns about past twenty."""
    before = len(plotting.plt.get_fignums())
    for i in range(25):
        fig, ax = plotting.plt.subplots()
        ax.plot([0, 1], [0, 1])
        plotting.savefig(fig, tmp_path / f"f{i}.png")
    assert len(plotting.plt.get_fignums()) == before

    fig, _ = plotting.plt.subplots()
    plotting.savefig(fig, tmp_path / "kept.png", close=False)
    assert fig.number in plotting.plt.get_fignums()
    plotting.plt.close(fig)


def test_event_psi_shows_every_replicate_and_both_group_means():
    """The number a differential test reports is a summary; this is the measurement
    behind it. A reader who cannot see the replicates has to take ΔΨ on trust."""
    import matplotlib.pyplot as plt
    import pandas as pd

    from splicescope.plotting import plot_event_psi

    psi = pd.DataFrame(
        {
            "event_id": ["E1"] * 6,
            "sample": ["C1", "C2", "C3", "K1", "K2", "K3"],
            "psi": [0.04, 0.05, 0.03, 0.22, 0.26, 0.24],
        }
    )
    groups = {"C1": "ctrl", "C2": "ctrl", "C3": "ctrl", "K1": "kd", "K2": "kd", "K3": "kd"}

    fig, ax = plt.subplots()
    plot_event_psi(psi, groups, ax=ax)
    points = [c for c in ax.collections if len(getattr(c, "get_offsets", lambda: [])())]
    plotted = sum(len(c.get_offsets()) for c in points)
    assert plotted >= 6, "every replicate is drawn, not just a mean"
    assert [t.get_text() for t in ax.get_xticklabels()] == ["ctrl", "kd"]
    assert ax.get_ylabel().startswith("Ψ")
    plt.close(fig)


def test_event_psi_ignores_a_sample_the_design_does_not_name():
    """A Ψ table can carry a sample the groups file left out; it belongs to neither
    condition and must not be drawn under one of them."""
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    from splicescope.plotting import plot_event_psi

    psi = pd.DataFrame(
        {"event_id": ["E1"] * 3, "sample": ["C1", "K1", "STRAY"], "psi": [0.1, 0.5, 0.9]}
    )
    fig, ax = plt.subplots()
    plot_event_psi(psi, {"C1": "ctrl", "K1": "kd"}, ax=ax)
    drawn = np.concatenate(
        [c.get_offsets()[:, 1] for c in ax.collections if len(c.get_offsets())]
    )
    assert 0.9 not in set(np.round(drawn, 6)), "the unnamed sample was plotted anyway"
    plt.close(fig)


def test_every_plot_makes_its_own_figure_when_given_no_axis():
    """`ax=None` is the default and so the first thing anyone calls, and for two of these
    it had never run. The figure it creates is reachable through `ax.figure`, which is how
    a caller saves or closes it — `savefig` takes a figure, not an axis — so that is what
    this checks rather than merely that nothing raised.
    """
    import matplotlib.pyplot as plt

    from splicescope.consequence import CONSEQUENCE_CLASSES

    diff = pd.DataFrame(
        {
            "chrom": ["chr1"] * 4, "start": [10, 20, 30, 40], "end": [90, 80, 70, 60],
            "strand": ["+"] * 4, "sclass": ["cryptic"] * 4,
            "delta_psi": [0.4, -0.3, 0.02, 0.5], "qvalue": [1e-8, 1e-3, 0.9, 1e-12],
            "event_type": ["SE", "MXE", "A5SS", "A3SS"],
        }
    )
    calls = {
        "plot_annotation_summary": (
            pd.DataFrame({"sclass": ["annotated", "cryptic"], "n_junctions": [7, 3]}),
        ),
        "plot_volcano": (diff,),
        "plot_event_summary": (pd.DataFrame({"event_type": ["SE", "SE", "MXE"]}),),
        "plot_event_volcano": (diff.assign(event_id=["a", "b", "c", "d"]),),
        "plot_importance": (
            pd.DataFrame({"feature": ["a", "b"], "importance": [0.6, 0.4],
                          "std": [0.05, 0.04]}),
        ),
        "plot_consequence_summary": (
            pd.DataFrame({"consequence_class": list(CONSEQUENCE_CLASSES)[:3]}),
        ),
        "plot_event_psi": (
            pd.DataFrame({"event_id": ["e"] * 4, "sample": ["c1", "c2", "k1", "k2"],
                          "psi": [0.2, 0.25, 0.7, 0.75]}),
            {"c1": "ctrl", "c2": "ctrl", "k1": "kd", "k2": "kd"},
        ),
        "plot_enrichment": (
            pd.DataFrame({"term": ["t1", "t2"], "qvalue": [1e-4, 1e-2],
                          "overlap": [3, 2], "set_size": [10, 8]}),
        ),
        "plot_roc": (np.array([0, 0, 1, 1]), np.array([0.1, 0.4, 0.6, 0.9])),
    }
    assert set(calls) == {n for n in dir(plotting) if n.startswith("plot_")}, (
        "a plot function was added or renamed and is not exercised here"
    )

    for name, args in calls.items():
        before = set(plt.get_fignums())
        ax = getattr(plotting, name)(*args)
        made = set(plt.get_fignums()) - before
        assert len(made) == 1, f"{name} should create exactly one figure, made {len(made)}"
        assert ax.figure.number in made, f"{name}'s axis must belong to the figure it made"
        assert ax.get_children(), f"{name} drew nothing"
        plt.close(ax.figure)
