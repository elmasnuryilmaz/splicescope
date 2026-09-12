"""The dashboard itself, rendered.

`tests/test_demo.py` covers the analysis the page runs; this runs the page. Streamlit
ships a headless harness for exactly this, so the Streamlit calls — which no amount of
testing the analysis can check — are exercised too. Skipped when the optional `app`
extra is not installed.
"""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("streamlit", reason="the dashboard needs the optional 'app' extra")

from streamlit.testing.v1 import AppTest  # noqa: E402

APP = Path(__file__).resolve().parent.parent / "app" / "streamlit_app.py"


def _run(sliders: dict[int, float] | None = None):
    app = AppTest.from_file(str(APP), default_timeout=300)
    app.run()
    if sliders:
        for index, value in sliders.items():
            app.slider[index].set_value(value)
        app.run()
    return app


def test_the_dashboard_renders_with_its_default_settings():
    app = _run()
    assert not app.exception, app.exception
    # the page's own catch-all would have fired if any step had raised
    assert not app.error, [e.value for e in app.error]
    assert [m.label for m in app.metric] == [
        "junctions tested",
        "significant (ΔΨ)",
        "cryptic classifier ROC-AUC",
    ]
    assert app.metric[2].value != "—", "the classifier should have run on the defaults"
    assert len(app.dataframe) == 1, "the ranked cryptic candidates"
    assert not app.info


def test_the_sliders_that_used_to_crash_the_page_now_explain_themselves():
    """Sliders 0-3 are genes, replicates, cryptic fraction and label noise. The minimum
    cryptic fraction with no label noise leaves one class, which raised `IndexError`
    from inside `predict_proba` and put a traceback on a public page."""
    app = _run({0: 5, 1: 3, 2: 0.1, 3: 0.0})
    assert not app.exception, app.exception
    assert not app.error, [e.value for e in app.error]
    assert app.info, "the page must say why the classifier is missing"
    assert "no two classes to separate" in app.info[0].value
    # and the rest of the page is still there
    assert app.metric[2].value == "—"
    assert len(app.metric) == 3
    assert not app.dataframe


def test_the_thresholds_change_the_significant_count_without_reanalysing():
    """The q-value and ΔΨ sliders are applied after the cached analysis, so moving them
    must change what is reported without touching the pipeline."""
    app = AppTest.from_file(str(APP), default_timeout=300)
    app.run()
    strict = int(app.metric[1].value)
    app.slider[4].set_value(0.2)   # q-value threshold
    app.slider[5].set_value(0.0)   # min |ΔΨ|
    app.run()
    assert not app.exception, app.exception
    assert int(app.metric[1].value) >= strict


def test_every_figure_the_page_draws_is_released():
    """This is why the deployed app went over its memory limit.

    `st.pyplot` does *not* close the figure it renders — its `clear_figure` default is
    False — and `pyplot` keeps every figure in a global registry until something does.
    A Streamlit script reruns on each widget change, so four figures accumulate per
    slider move, without bound, against an analysis that costs 47 MB and imports that
    cost 194.

    The check is structural rather than a count, and deliberately so: `AppTest` runs the
    script with its own module state, so the registry this process can see is not the one
    the page uses. A test that counted `plt.get_fignums()` from here passed happily with
    the leak reinstated, which is how this test came to be written this way.
    """
    import ast

    tree = ast.parse(APP.read_text())
    functions = [n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)]

    def calls(node, dotted):
        owner, attr = dotted.split(".")
        return [
            c for c in ast.walk(node)
            if isinstance(c, ast.Call) and isinstance(c.func, ast.Attribute)
            and c.func.attr == attr
            and isinstance(c.func.value, ast.Name) and c.func.value.id == owner
        ]

    renderers = [f for f in functions if calls(f, "st.pyplot")]
    assert renderers, "no st.pyplot call found — has the page stopped drawing?"
    for function in renderers:
        assert calls(function, "plt.close"), (
            f"{function.name}() renders a figure and never closes it"
        )

    # and nothing may render outside those functions, where no close would apply
    inside = {id(c) for f in renderers for c in calls(f, "st.pyplot")}
    assert all(id(c) in inside for c in calls(tree, "st.pyplot")), (
        "a figure is rendered at module level, outside any function that closes it"
    )
