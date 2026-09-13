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
        "cassette exons found",
        "predicted to trigger decay",
    ]
    assert int(app.metric[2].value) > 0, "the page is built around a cassette exon"
    assert int(app.metric[3].value) > 0, "at least one should be predicted to trigger decay"
    assert [s.value for s in app.subheader] == [
        "One event, end to end",
        "The whole dataset",
        "Telling real cryptic junctions from noise",
    ]
    assert len(app.selectbox) == 1, "one event is chosen and shown in full"
    assert app.table, "its numbers"
    assert app.dataframe, "the ranked cryptic candidates"
    assert not app.info


def test_the_chosen_event_is_explained_in_words_not_only_in_columns():
    """The point of the page. A table saying `ptc_nmd`, `insert_length=61`,
    `distance_to_last_junction=395` is an answer only to someone who already knows the
    rule; the page has to say what it means."""
    app = _run()
    prose = " ".join(m.value for m in app.markdown) + " ".join(
        getattr(w, "value", "") for w in app.get("write")
    )
    assert "premature stop" in prose, "the event's consequence must be spelled out"
    assert "nonsense-mediated decay" in prose or "truncated protein" in prose
    # and the argument for the test it used
    assert "Mann-Whitney" in prose and "count-based" in prose


def test_the_sliders_that_leave_nothing_to_show_explain_themselves():
    """Sliders 0-3 are genes, replicates, cryptic fraction and label noise. The minimum
    cryptic fraction with no label noise leaves no cassette exon to walk through and one
    class for the classifier. Both used to be failures; both are now sentences."""
    app = _run({0: 5, 1: 3, 2: 0.1, 3: 0.0})
    assert not app.exception, app.exception
    assert not app.error, [e.value for e in app.error]

    messages = " ".join(i.value for i in app.info)
    assert "No cassette exon was detected" in messages
    assert "no two classes to separate" in messages
    assert not app.selectbox, "nothing to pick from"
    assert int(app.metric[2].value) == 0


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


def test_the_page_tests_for_a_missing_value_in_a_way_that_cannot_raise():
    """`x == x` is a compact way to ask "is this not NaN", and it works right up until
    the column is nullable — `pd.NA == pd.NA` is `pd.NA`, whose truth value raises rather
    than being False. The page used it on two values it takes from a merge, and a
    `TypeError` there does not produce a missing line: it replaces the whole page with an
    error box, because the analysis runs inside one `try`.

    That is not hypothetical. Event coordinates became a nullable integer column on the
    day this was written, and the same idiom in a property test raised immediately. The
    check is structural, because the values happen to be floating point today and a test
    that ran the page would pass either way.
    """
    import ast

    tree = ast.parse(APP.read_text())
    offenders = [
        ast.unparse(node)
        for node in ast.walk(tree)
        if isinstance(node, ast.Compare)
        and len(node.ops) == 1
        and isinstance(node.ops[0], ast.Eq)
        and ast.unparse(node.left) == ast.unparse(node.comparators[0])
    ]
    assert not offenders, (
        "compare with pd.isna instead of a value against itself: " + "; ".join(offenders)
    )
    assert "pd.isna(" in APP.read_text(), "and it is asked somewhere"
