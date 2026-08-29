"""Render docs/demo.gif — how the cryptic classifier holds up as label noise grows.

For a sweep of ``label_noise`` we recompute the whole pipeline and draw the
cross-validated ROC. The curve degrades gracefully from near-perfect to clearly
imperfect, which is the honest robustness story. Dark theme matches the website.

    python docs/make_demo_gif.py
"""

from __future__ import annotations

import io
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import StratifiedKFold, cross_val_predict

from splicescope import annotate, cryptic, quantify
from splicescope.ml import CrypticClassifier
from splicescope.simulate import simulate_dataset

OUT = Path(__file__).resolve().parent / "demo.gif"

BG = "#07060f"
FG = "#f4f1ff"
MUTED = "#a39fc4"
ACCENT = "#ff3d81"


def roc_for_noise(label_noise: float):
    ds = simulate_dataset(
        n_genes=18, n_per_group=6, cryptic_fraction=0.6, label_noise=label_noise, seed=11
    )
    ann = annotate.annotate_junctions(ds.observed, ds.known)
    psi = quantify.compute_psi(ann, min_reads=5)
    feats = cryptic.extract_features(psi, ds.known)
    clf = CrypticClassifier(n_estimators=200, random_state=0)
    x, y = clf._xy(feats)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
    proba = cross_val_predict(clf._make_pipeline(), x, y, cv=cv, method="predict_proba")[:, 1]
    fpr, tpr, _ = roc_curve(y, proba)
    return fpr, tpr, roc_auc_score(y, proba)


def frame(label_noise: float) -> Image.Image:
    fpr, tpr, auc = roc_for_noise(label_noise)
    fig, ax = plt.subplots(figsize=(6, 5), dpi=110)
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)
    ax.plot([0, 1], [0, 1], ls="--", lw=1, color=MUTED, alpha=0.5)
    ax.plot(fpr, tpr, lw=3, color=ACCENT, solid_capstyle="round")
    ax.fill_between(fpr, tpr, color=ACCENT, alpha=0.12)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("false positive rate", color=MUTED)
    ax.set_ylabel("true positive rate", color=MUTED)
    ax.set_title("splicescope — cryptic classifier vs. label noise", color=FG, fontsize=13, pad=14)
    ax.text(0.97, 0.10, f"AUC = {auc:.2f}", color=FG, ha="right", fontsize=15, weight="bold")
    ax.text(
        0.97, 0.03, f"label noise = {label_noise:.0%}", color=ACCENT, ha="right", fontsize=11
    )
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(MUTED)
    ax.tick_params(colors=MUTED)
    fig.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="png", facecolor=BG)
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).convert("RGB")


def main() -> None:
    noises = [round(x, 3) for x in np.linspace(0.0, 0.30, 11)]
    print("rendering frames:", noises)
    frames = [frame(n) for n in noises]
    frames = frames + frames[-2:0:-1]  # ping-pong loop
    frames[0].save(
        OUT,
        save_all=True,
        append_images=frames[1:],
        duration=260,
        loop=0,
        optimize=True,
    )
    print(f"wrote {OUT} ({OUT.stat().st_size // 1024} KB, {len(frames)} frames)")


if __name__ == "__main__":
    main()
