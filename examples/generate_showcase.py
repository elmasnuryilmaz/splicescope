"""End-to-end demo that produces the showcase figure used in the README.

Runs the whole splicescope pipeline on a reproducible synthetic dataset and
writes a 2x2 panel: junction classes, a ΔΨ volcano, permutation importances and
the cryptic-classifier ROC curve.

    python examples/generate_showcase.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold, cross_val_predict

from splicescope import annotate, cryptic, diff, plotting, quantify
from splicescope.ml import CrypticClassifier
from splicescope.simulate import simulate_dataset

OUT = Path(__file__).resolve().parent.parent / "docs" / "showcase.png"


def main() -> None:
    # label_noise mimics imperfect curation, so the ML task is realistically hard
    ds = simulate_dataset(
        n_genes=20, n_per_group=6, cryptic_fraction=0.6, label_noise=0.12, seed=11
    )

    annotated = annotate.annotate_junctions(ds.observed, ds.known)
    summary = annotate.annotation_summary(annotated)
    psi = quantify.compute_psi(annotated, min_reads=5)
    dsplice = diff.differential_splicing(psi, ds.groups)

    feats = cryptic.extract_features(psi, ds.known)
    clf = CrypticClassifier(random_state=0)
    metrics = clf.evaluate(feats)
    clf.fit(feats)

    # out-of-fold probabilities for an honest ROC
    x, y = clf._xy(feats)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
    oof = cross_val_predict(clf._make_pipeline(), x, y, cv=cv, method="predict_proba")[:, 1]

    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
    fig.suptitle("splicescope — from splice junctions to cryptic-event calls", fontsize=14, y=0.98)
    plotting.plot_annotation_summary(summary, ax=axes[0, 0])
    plotting.plot_volcano(dsplice, ax=axes[0, 1])
    plotting.plot_importance(clf.importances, ax=axes[1, 0])
    plotting.plot_roc(y, oof, ax=axes[1, 1])
    plotting.savefig(fig, OUT)

    print(
        f"cryptic classifier: ROC-AUC={metrics['roc_auc']:.3f} "
        f"AP={metrics['average_precision']:.3f}"
    )
    print(f"significant junctions: {len(diff.significant(dsplice))}")
    print(f"showcase figure -> {OUT}")


if __name__ == "__main__":
    main()
