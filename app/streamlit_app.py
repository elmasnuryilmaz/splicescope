"""Interactive splicescope dashboard.

    streamlit run app/streamlit_app.py

Explore junction classes, differential splicing and ranked cryptic candidates on a
reproducible synthetic dataset — adjust the sidebar and everything recomputes live.
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
import streamlit as st

st.set_page_config(page_title="splicescope", page_icon="🧬", layout="wide")
st.title("🧬 splicescope — cryptic splicing explorer")
st.caption("From splice junctions to cryptic-event calls. Data below is reproducible & synthetic.")

try:
    from splicescope import annotate, cryptic, diff, plotting, quantify
    from splicescope.ml import CrypticClassifier
    from splicescope.simulate import simulate_dataset
except Exception as exc:  # pragma: no cover - defensive import guard for deployment
    st.error("Could not import splicescope. Please try again shortly.")
    st.exception(exc)
    st.stop()


with st.sidebar:
    st.header("Dataset")
    n_genes = st.slider("genes", 5, 40, 20)
    n_rep = st.slider("replicates / group", 3, 10, 6)
    cryptic_fraction = st.slider("fraction of genes with a cryptic event", 0.1, 1.0, 0.6)
    label_noise = st.slider("label noise (curation error)", 0.0, 0.3, 0.12)
    seed = st.number_input("seed", value=11, step=1)
    q_thr = st.slider("q-value threshold", 0.01, 0.2, 0.05)
    delta_thr = st.slider("min |ΔΨ|", 0.0, 0.5, 0.1)


@st.cache_data(show_spinner=True)
def run_pipeline(n_genes, n_rep, cryptic_fraction, label_noise, seed):
    ds = simulate_dataset(
        n_genes=n_genes,
        n_per_group=n_rep,
        cryptic_fraction=cryptic_fraction,
        label_noise=label_noise,
        seed=int(seed),
    )
    ann = annotate.annotate_junctions(ds.observed, ds.known)
    summary = annotate.annotation_summary(ann)
    psi = quantify.compute_psi(ann, min_reads=5)
    dsplice = diff.differential_splicing(psi, ds.groups)
    feats = cryptic.extract_features(psi, ds.known)
    clf = CrypticClassifier(random_state=0).fit(feats)
    metrics = clf.evaluate(feats)
    scores = clf.score_table(feats)
    # return plain, serialisable objects (not the fitted model) for the cache
    return summary, dsplice, clf.importances, metrics, scores


try:
    summary, dsplice, importances, metrics, scores = run_pipeline(
        n_genes, n_rep, cryptic_fraction, label_noise, seed
    )

    c1, c2, c3 = st.columns(3)
    c1.metric("junctions tested", len(dsplice))
    c2.metric("significant (ΔΨ)", len(diff.significant(dsplice, q=q_thr, min_delta=delta_thr)))
    c3.metric("cryptic classifier ROC-AUC", f"{metrics['roc_auc']:.3f}")

    left, right = st.columns(2)
    with left:
        st.subheader("Junction classes")
        fig, ax = plt.subplots(figsize=(5, 3.4))
        plotting.plot_annotation_summary(summary, ax=ax)
        st.pyplot(fig)
    with right:
        st.subheader("Differential splicing")
        fig, ax = plt.subplots(figsize=(5, 3.8))
        plotting.plot_volcano(dsplice, q=q_thr, min_delta=delta_thr, ax=ax)
        st.pyplot(fig)

    st.subheader("Top cryptic candidates")
    st.dataframe(scores.head(25))

    st.subheader("What the classifier keys on")
    fig, ax = plt.subplots(figsize=(6, 3.2))
    plotting.plot_importance(importances, ax=ax)
    st.pyplot(fig)
except Exception as exc:  # pragma: no cover - keep the demo from showing a blank error
    st.error("Something went wrong while running the analysis.")
    st.exception(exc)
