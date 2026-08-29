"""Learn to separate true cryptic-splicing junctions from noise.

Novel junctions are abundant and mostly artefactual. Rather than hand-tuning
thresholds, we train a classifier on the engineered features from
:mod:`splicescope.cryptic` and evaluate it honestly with stratified
cross-validation. The model ships with a *model card* (metrics + permutation
importances + provenance) so a reviewer can see exactly how it was assessed.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .cryptic import FEATURE_COLUMNS

TRUTH_COLUMN = "is_cryptic_truth"


@dataclass
class CrypticClassifier:
    """A cross-validated classifier for cryptic-junction calling."""

    n_estimators: int = 300
    random_state: int = 0
    n_splits: int = 5
    features: list[str] = field(default_factory=lambda: list(FEATURE_COLUMNS))
    pipeline: Pipeline | None = None
    cv_metrics: dict = field(default_factory=dict)
    importances: pd.DataFrame | None = None

    def _make_pipeline(self) -> Pipeline:
        return Pipeline(
            [
                ("scale", StandardScaler()),
                (
                    "rf",
                    RandomForestClassifier(
                        n_estimators=self.n_estimators,
                        random_state=self.random_state,
                        class_weight="balanced",
                        n_jobs=-1,
                    ),
                ),
            ]
        )

    def _xy(self, feats: pd.DataFrame):
        x = feats[self.features].to_numpy(dtype=float)
        x = np.nan_to_num(x, nan=0.0)
        y = feats[TRUTH_COLUMN].to_numpy(dtype=int)
        return x, y

    def evaluate(self, feats: pd.DataFrame) -> dict:
        """Stratified-CV evaluation. Returns ROC-AUC and average precision."""
        x, y = self._xy(feats)
        n_splits = min(self.n_splits, int(np.bincount(y).min()))
        if n_splits < 2:
            raise ValueError("need at least 2 examples of each class for cross-validation")
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=self.random_state)
        proba = cross_val_predict(
            self._make_pipeline(), x, y, cv=cv, method="predict_proba", n_jobs=None
        )[:, 1]
        self.cv_metrics = {
            "n_junctions": int(len(y)),
            "n_positive": int(y.sum()),
            "n_splits": int(n_splits),
            "roc_auc": float(roc_auc_score(y, proba)),
            "average_precision": float(average_precision_score(y, proba)),
        }
        return self.cv_metrics

    def fit(self, feats: pd.DataFrame) -> CrypticClassifier:
        """Fit on all data and compute permutation importances."""
        x, y = self._xy(feats)
        self.pipeline = self._make_pipeline().fit(x, y)
        imp = permutation_importance(
            self.pipeline, x, y, n_repeats=20, random_state=self.random_state, n_jobs=-1
        )
        self.importances = (
            pd.DataFrame(
                {
                    "feature": self.features,
                    "importance": imp.importances_mean,
                    "std": imp.importances_std,
                }
            )
            .sort_values("importance", ascending=False)
            .reset_index(drop=True)
        )
        return self

    def predict_proba(self, feats: pd.DataFrame) -> np.ndarray:
        if self.pipeline is None:
            raise RuntimeError("call fit() first")
        x = np.nan_to_num(feats[self.features].to_numpy(dtype=float), nan=0.0)
        return self.pipeline.predict_proba(x)[:, 1]

    def score_table(self, feats: pd.DataFrame) -> pd.DataFrame:
        """Return junctions with a ``cryptic_score`` column, most-likely first."""
        out = feats.copy()
        out["cryptic_score"] = self.predict_proba(feats)
        cols = ["chrom", "start", "end", "strand", "gene_id", "sclass", "cryptic_score"]
        cols = [c for c in cols if c in out.columns]
        return out.sort_values("cryptic_score", ascending=False)[
            cols + [c for c in out.columns if c not in cols]
        ].reset_index(drop=True)

    def model_card(self) -> dict:
        return {
            "model": "RandomForestClassifier (StandardScaler pipeline)",
            "task": "binary classification — true cryptic junction vs noise",
            "features": self.features,
            "hyperparameters": {
                "n_estimators": self.n_estimators,
                "class_weight": "balanced",
                "random_state": self.random_state,
            },
            "evaluation": "stratified k-fold cross-validation (no leakage)",
            "cv_metrics": self.cv_metrics,
            "permutation_importance": (
                self.importances.to_dict(orient="records") if self.importances is not None else []
            ),
            "intended_use": "prioritising candidate cryptic-splicing junctions for review",
            "limitations": (
                "trained/validated on simulated data in this demo; on real data, retrain with "
                "curated labels and validate on held-out genes to avoid optimistic estimates"
            ),
        }

    def write_model_card(self, path: str | Path) -> Path:
        path = Path(path)
        card = self.model_card()
        path.write_text(json.dumps(card, indent=2))
        md = path.with_suffix(".md")
        lines = [
            "# Model card — cryptic junction classifier",
            "",
            f"**Model:** {card['model']}  ",
            f"**Task:** {card['task']}  ",
            f"**Evaluation:** {card['evaluation']}",
            "",
            "## Cross-validated metrics",
        ]
        for k, v in card["cv_metrics"].items():
            lines.append(f"- **{k}**: {v}")
        lines += ["", "## Permutation importance"]
        for row in card["permutation_importance"]:
            lines.append(f"- `{row['feature']}`: {row['importance']:.3f} ± {row['std']:.3f}")
        lines += ["", "## Limitations", card["limitations"], ""]
        md.write_text("\n".join(lines))
        return path
