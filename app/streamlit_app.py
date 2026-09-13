"""Interactive splicescope dashboard.

    streamlit run app/streamlit_app.py

The page answers one question end to end: a knockdown switches on splice junctions the
annotation does not contain — which of them are real, and what do they do to the protein?

The analysis lives in :func:`splicescope.demo.run_demo`, so the test suite exercises
exactly what this page shows. This file is presentation only.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Make the app importable whether or not the package was pip-installed. On Streamlit
# Community Cloud the repo is cloned and the source lives in ../src; adding it to the
# path means the app never depends on a local-package build step succeeding.
_SRC = Path(__file__).resolve().parent.parent / "src"
if _SRC.is_dir() and str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

st.set_page_config(page_title="splicescope", page_icon="🧬", layout="wide")
st.title("🧬 splicescope — cryptic splicing explorer")
st.caption(
    "A knockdown switches on splice junctions the annotation does not contain. "
    "Which of them are real, and what do they do to the protein they land in? "
    "Everything below is computed live on reproducible synthetic data."
)

try:
    from splicescope import diff, plotting
    from splicescope.betabinom import min_achievable_rank_pvalue
    from splicescope.consequence import describe
    from splicescope.demo import run_demo
except Exception as exc:  # pragma: no cover - defensive import guard for deployment
    st.error("Could not import splicescope. Please try again shortly.")
    st.exception(exc)
    st.stop()


with st.sidebar:
    st.header("Dataset")
    st.caption("The data is simulated, so the truth is known and every call can be scored.")
    n_genes = st.slider("genes", 5, 40, 20)
    n_rep = st.slider("replicates / group", 3, 10, 6)
    cryptic_fraction = st.slider("fraction of genes with a cryptic event", 0.1, 1.0, 0.6)
    label_noise = st.slider("label noise (curation error)", 0.0, 0.3, 0.12)
    seed = st.number_input("seed", value=11, step=1)
    st.header("Calling")
    q_thr = st.slider("q-value threshold", 0.01, 0.2, 0.05)
    delta_thr = st.slider("min |ΔΨ|", 0.0, 0.5, 0.1)


def draw(figsize, plot):
    """Render one figure and let go of it.

    `st.pyplot` does not close the figure it renders, and `pyplot` keeps every figure in
    a global registry until something does. A Streamlit script reruns on every widget
    change, so left open they accumulate without bound — which is how this app went over
    its memory limit on Community Cloud.
    """
    fig, ax = plt.subplots(figsize=figsize)
    try:
        plot(ax)
        st.pyplot(fig)
    finally:
        plt.close(fig)


# Bounded, so a session spent dragging sliders cannot grow the cache without limit.
@st.cache_data(show_spinner="Simulating, annotating, testing, predicting…", max_entries=8, ttl=3600)
def analyse(n_genes, n_rep, cryptic_fraction, label_noise, seed):
    """Cached wrapper. Returns plain, serialisable pieces rather than the fitted model."""
    result = run_demo(
        n_genes=n_genes,
        n_per_group=n_rep,
        cryptic_fraction=cryptic_fraction,
        label_noise=label_noise,
        seed=int(seed),
    )
    return (
        result.annotation_summary,
        result.differential,
        result.metrics,
        result.importances,
        result.scores,
        result.classifier_note,
        result.events,
        result.event_psi,
        result.event_differential,
        result.consequences,
        result.groups,
    )


def event_label(row, consequence):
    """`g07 · chr1:20760-20820 · predicted to trigger decay`"""
    reading = {
        "ptc_nmd": "predicted to trigger decay",
        "ptc_escape": "premature stop, decay escaped",
        "frameshift": "frame shifted, no premature stop",
        "in_frame_insertion": "in frame",
        "exon_truncation": "shorter protein",
        "utr_insertion": "outside the coding sequence",
        "non_coding_host": "non-coding host",
        "no_host_transcript": "no host transcript",
    }.get(consequence, consequence)
    return f"{row.gene_id} · {row.chrom}:{int(row.exon_start)}-{int(row.exon_end)} · {reading}"


try:
    (summary, dsplice, metrics, importances, scores, note,
     events, epsi, ediff, consequences, groups) = analyse(
        n_genes, n_rep, cryptic_fraction, label_noise, seed
    )

    hits = diff.significant(dsplice, q=q_thr, min_delta=delta_thr)
    decay = 0 if consequences is None else int(consequences["nmd_predicted"].sum())
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("junctions tested", len(dsplice))
    c2.metric("significant (ΔΨ)", len(hits))
    c3.metric("cassette exons found", 0 if consequences is None else len(consequences))
    c4.metric("predicted to trigger decay", decay)

    # ---------------------------------------------------------------- one event
    st.markdown("---")
    st.subheader("One event, end to end")
    if consequences is None or consequences.empty:
        st.info(
            "No cassette exon was detected in this dataset. Raise the fraction of genes "
            "with a cryptic event, or the number of genes."
        )
    else:
        table = consequences.copy()
        if ediff is not None and "event_id" in table.columns:
            table = table.merge(
                ediff[["event_id", "delta_psi", "pvalue", "qvalue"]], on="event_id", how="left"
            )
        table = table.sort_values(
            ["nmd_predicted", "qvalue"], ascending=[False, True], na_position="last"
        ).reset_index(drop=True)

        labels = [
            event_label(row, row.consequence_class) for row in table.itertuples(index=False)
        ]
        chosen = st.selectbox(
            "Pick a cassette exon — the strongest call is first",
            range(len(labels)),
            format_func=lambda i: labels[i],
        )
        row = table.iloc[chosen]

        left, right = st.columns([1, 1.35])
        with left:
            if epsi is not None and groups:
                rows = epsi[epsi["event_id"] == row.event_id]
                if not rows.empty:
                    draw(
                        (4.6, 3.4),
                        lambda ax: plotting.plot_event_psi(rows, groups, ax=ax),
                    )
        with right:
            st.markdown("**What it does to the protein**")
            st.write(describe(row))
            numbers = {
                "gene": row.get("gene_name") or row.get("gene_id"),
                "exon": f"{row.chrom}:{int(row.exon_start)}-{int(row.exon_end)}",
                "length added": f"{int(row.insert_length)} nt",
                "reading frame": "shifted" if row.frameshift else "preserved",
                "call": row.consequence_class,
            }
            # pd.isna, not `x == x`: the latter raises rather than returning False on a
            # nullable column, and this page is where a user would see it happen
            if not pd.isna(row.get("delta_psi")):
                numbers["ΔΨ"] = f"{row.delta_psi:+.3f}"
                numbers["q-value"] = f"{row.qvalue:.2e}"
            st.table({"": list(numbers), " ": list(numbers.values())})

        # ------------------------------------------------ why count-based testing
        if not pd.isna(row.get("pvalue")):
            floor = min_achievable_rank_pvalue(n_rep, n_rep)
            st.markdown("**Why the read counts, and not the ranks**")
            st.write(
                f"A two-sided Mann-Whitney on {n_rep} against {n_rep} replicates cannot "
                f"return a p-value below **{floor:.3g}** however large the difference — it "
                "has only the ordering of the samples to work with, and there are just "
                f"{2 * n_rep} of them. This event's count-based p-value is "
                f"**{row.pvalue:.1e}**, because evidence accumulates with coverage as well "
                "as with replicates. After correction across a genome's worth of junctions, "
                "that is the difference between reporting this exon and reporting nothing."
            )

    # ---------------------------------------------------------------- everything else
    st.markdown("---")
    st.subheader("The whole dataset")
    left, right = st.columns(2)
    with left:
        draw((5, 3.4), lambda ax: plotting.plot_annotation_summary(summary, ax=ax))
    with right:
        draw(
            (5, 3.8),
            lambda ax: plotting.plot_volcano(dsplice, q=q_thr, min_delta=delta_thr, ax=ax),
        )

    if events is not None and not events.empty:
        left, right = st.columns(2)
        with left:
            draw((5, 3.2), lambda ax: plotting.plot_event_summary(events, ax=ax))
        with right:
            if consequences is not None and not consequences.empty:
                draw(
                    (5.5, 3.2),
                    lambda ax: plotting.plot_consequence_summary(consequences, ax=ax),
                )

    st.markdown("---")
    st.subheader("Telling real cryptic junctions from noise")
    if note:
        # a legitimate slider position, not an error: say what to change
        st.info(note)
    else:
        st.caption(
            f"A classifier trained on the junctions' own features — length, support, "
            f"motif, distance to the nearest annotated site — and scored by "
            f"cross-validation, never on the data it was fitted to. ROC-AUC "
            f"{metrics['roc_auc']:.3f}."
        )
        left, right = st.columns([1.2, 1])
        with left:
            st.dataframe(scores.head(20), height=340)
        with right:
            draw((5, 3.4), lambda ax: plotting.plot_importance(importances, ax=ax))
except Exception as exc:  # pragma: no cover - keep the demo from showing a blank error
    st.error("Something went wrong while running the analysis.")
    st.exception(exc)
