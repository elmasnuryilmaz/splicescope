"""Publication-quality plots (matplotlib only, theme-consistent)."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless-safe for CI and servers
import matplotlib.pyplot as plt
import numpy as np

_ACCENT = "#ff3d81"
_ACCENT2 = "#8b5cf6"
_CYAN = "#22d3ee"
_MUTED = "#9aa0b4"


def _style(ax):
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(colors="#444")
    ax.grid(axis="y", alpha=0.15)


def plot_annotation_summary(summary, ax=None):
    """Bar chart of unique junctions per class."""
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 3.2))
    ax.bar(summary["sclass"], summary["n_junctions"], color=_ACCENT2)
    ax.set_ylabel("unique junctions")
    ax.set_title("Junction classes")
    ax.set_xticks(range(len(summary)))
    ax.set_xticklabels(summary["sclass"], rotation=30, ha="right", fontsize=8)
    _style(ax)
    return ax


def plot_volcano(diff, q: float = 0.05, min_delta: float = 0.1, ax=None):
    """ΔΨ volcano: effect size vs −log10 q-value."""
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 4))
    d = diff.dropna(subset=["qvalue", "delta_psi"]).copy()
    y = -np.log10(d["qvalue"].clip(lower=1e-300))
    sig = (d["qvalue"] <= q) & (d["delta_psi"].abs() >= min_delta)
    ax.scatter(d.loc[~sig, "delta_psi"], y[~sig], s=14, color=_MUTED, alpha=0.6, label="ns")
    ax.scatter(d.loc[sig, "delta_psi"], y[sig], s=22, color=_ACCENT, label="significant")
    ax.axhline(-np.log10(q), ls="--", lw=0.8, color="#888")
    ax.axvline(min_delta, ls="--", lw=0.8, color="#888")
    ax.axvline(-min_delta, ls="--", lw=0.8, color="#888")
    ax.set_xlabel("ΔΨ (condition B − A)")
    ax.set_ylabel("−log10 q-value")
    ax.set_title("Differential splicing")
    ax.legend(frameon=False, fontsize=8)
    _style(ax)
    return ax


def plot_importance(importances, ax=None):
    """Horizontal bar of permutation importances."""
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 3.4))
    imp = importances.iloc[::-1]
    ax.barh(imp["feature"], imp["importance"], xerr=imp["std"], color=_CYAN)
    ax.set_xlabel("permutation importance")
    ax.set_title("What drives cryptic calls")
    _style(ax)
    ax.grid(axis="x", alpha=0.15)
    return ax


def plot_roc(y_true, scores, ax=None):
    """ROC curve from labels and scores."""
    from sklearn.metrics import roc_auc_score, roc_curve

    if ax is None:
        _, ax = plt.subplots(figsize=(4, 4))
    fpr, tpr, _ = roc_curve(y_true, scores)
    auc = roc_auc_score(y_true, scores)
    ax.plot(fpr, tpr, color=_ACCENT, lw=2, label=f"AUC = {auc:.2f}")
    ax.plot([0, 1], [0, 1], ls="--", color="#888", lw=0.8)
    ax.set_xlabel("false positive rate")
    ax.set_ylabel("true positive rate")
    ax.set_title("Cryptic classifier ROC")
    ax.legend(frameon=False, fontsize=9)
    _style(ax)
    return ax


def savefig(fig, path: str | Path, dpi: int = 150):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    return path
