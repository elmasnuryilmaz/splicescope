import pandas as pd

from splicescope.annotate import annotate_junctions
from splicescope.cryptic import extract_features
from splicescope.ml import CrypticClassifier
from splicescope.quantify import compute_psi
from splicescope.simulate import simulate_dataset


def test_simulation_is_reproducible():
    a = simulate_dataset(seed=42)
    b = simulate_dataset(seed=42)
    pd.testing.assert_frame_equal(a.observed, b.observed)


def test_simulation_has_both_truth_classes():
    ds = simulate_dataset(n_genes=10, cryptic_fraction=1.0, seed=3)
    assert ds.observed["is_cryptic_truth"].nunique() == 2


def test_ml_recovers_signal():
    ds = simulate_dataset(n_genes=14, n_per_group=5, cryptic_fraction=0.7, seed=7)
    annotated = annotate_junctions(ds.observed, ds.known)
    psi = compute_psi(annotated, min_reads=5)
    feats = extract_features(psi, ds.known)
    assert feats["is_cryptic_truth"].nunique() == 2

    clf = CrypticClassifier(random_state=0)
    metrics = clf.evaluate(feats)
    # signal is strong by construction; the model must beat chance comfortably
    assert metrics["roc_auc"] > 0.75

    clf.fit(feats)
    scored = clf.score_table(feats)
    assert "cryptic_score" in scored.columns
    assert scored["cryptic_score"].is_monotonic_decreasing
    card = clf.model_card()
    assert card["cv_metrics"]["roc_auc"] == metrics["roc_auc"]
