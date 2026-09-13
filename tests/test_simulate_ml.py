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


def test_the_injected_events_are_written_out_next_to_the_data(tmp_path):
    """A ground-truth dataset whose ground truth stays in memory is only useful from
    inside the test suite. `truth.tsv` lands beside the SJ files so anyone can measure
    recall on the simulated data without reimplementing the generator's arithmetic."""
    import pandas as pd

    from splicescope.simulate import TRUTH_COLUMNS, simulate_dataset, write_dataset

    ds = simulate_dataset(n_genes=14, mxe_fraction=0.5, alt_ss_fraction=0.5, seed=2)
    out = write_dataset(ds, tmp_path / "sim", seed=2)

    written = pd.read_csv(out / "truth.tsv", sep="\t")
    assert list(written.columns) == TRUTH_COLUMNS
    assert len(written) == len(ds.truth)
    assert set(written["event_type"]) <= {"cryptic_exon", "A5SS", "A3SS", "MXE"}
    # coordinates survive as integers, not as 1220.0
    assert (written["intron_start"] % 1 == 0).all()
    assert written["exonA_start"].dropna().astype(int).gt(0).all()
    # every row names an intron the annotation actually contains
    known = {(r.start, r.end) for r in ds.known.itertuples(index=False)}
    for row in written.itertuples(index=False):
        assert (row.intron_start, row.intron_end) in known, row


def _separable_features(n_positive=20, n_negative=20, seed=0):
    """A feature table the classifier cannot get wrong, so what is tested is the
    plumbing rather than the learning."""
    import numpy as np

    from splicescope.cryptic import FEATURE_COLUMNS

    rng = np.random.default_rng(seed)
    frames = []
    for label, n, count, support, motif in (
        (1, n_positive, 5.0, 6, 1),
        (0, n_negative, 1.0, 1, 0),
    ):
        block = pd.DataFrame(
            {c: rng.normal(0.0, 0.1, n) for c in FEATURE_COLUMNS}
        )
        block["log_max_count"] = count + rng.normal(0, 0.05, n)
        block["n_samples_support"] = support
        block["canonical_motif"] = motif
        block["is_cryptic_truth"] = label
        frames.append(block)
    return pd.concat(frames, ignore_index=True)


def test_the_score_is_the_probability_of_being_cryptic_not_of_being_noise():
    """`predict_proba` returns a column per class, and taking the wrong one inverts
    every score. Nothing noticed: the score table's own ordering stays self-consistent
    either way, and `evaluate` reaches the probabilities by a different route."""
    feats = _separable_features()
    clf = CrypticClassifier(n_estimators=50, random_state=0).fit(feats)

    scores = clf.predict_proba(feats)
    positive = scores[feats["is_cryptic_truth"] == 1]
    negative = scores[feats["is_cryptic_truth"] == 0]
    assert positive.min() > 0.5, "a known cryptic junction must score above a half"
    assert negative.max() < 0.5, "a known artefact must score below it"
    assert positive.mean() > negative.mean()

    # and the same orientation must survive into the table a user reads
    table = clf.score_table(feats)
    assert table["cryptic_score"].iloc[0] > table["cryptic_score"].iloc[-1]


def test_the_fold_count_cannot_exceed_the_rarer_class():
    """True cryptic junctions are the rare class, and asking for more folds than there
    are positives either throws from inside scikit-learn or silently evaluates on folds
    with no positive in them."""
    feats = _separable_features(n_positive=3, n_negative=30)
    clf = CrypticClassifier(n_estimators=20, n_splits=5, random_state=0)
    metrics = clf.evaluate(feats)
    assert metrics["n_splits"] == 3
    assert metrics["n_positive"] == 3


def test_the_model_card_reports_the_model_that_was_actually_fitted():
    """The card exists so a reviewer does not have to take the method on trust, which
    it cannot do while it repeats literals from the source. `class_weight` was written
    into it as the string "balanced" regardless of what the pipeline had."""
    feats = _separable_features()
    clf = CrypticClassifier(n_estimators=50, class_weight=None, random_state=3).fit(feats)

    card = clf.model_card()["hyperparameters"]
    fitted = clf.pipeline.named_steps["rf"]
    assert card["class_weight"] == fitted.class_weight is None
    assert card["n_estimators"] == fitted.n_estimators == 50
    assert card["random_state"] == fitted.random_state == 3

    balanced = CrypticClassifier(n_estimators=50, random_state=0).fit(feats)
    assert balanced.model_card()["hyperparameters"]["class_weight"] == "balanced"


def test_the_classifier_refuses_a_dataset_too_small_to_cross_validate():
    """Coverage found this message had never been produced. It is the one a user meets
    first on a small pilot dataset, so it has to say what is wrong and what would fix
    it — a traceback out of scikit-learn would not."""
    import pandas as pd
    import pytest

    from splicescope.cryptic import FEATURE_COLUMNS
    from splicescope.ml import CrypticClassifier

    feats = pd.DataFrame(
        {c: [float(i) for i in range(12)] for c in FEATURE_COLUMNS}
        | {"is_cryptic_truth": [1] + [0] * 11}
    )
    with pytest.raises(ValueError, match="at least 2 examples of each class") as excinfo:
        CrypticClassifier().evaluate(feats)
    assert "has 1" in str(excinfo.value), "say how many the rarer class has"

    # two of the rarer class is enough to fold, and it says nothing
    feats.loc[1, "is_cryptic_truth"] = 1
    metrics = CrypticClassifier().evaluate(feats)
    assert metrics["n_splits"] == 2 and metrics["n_positive"] == 2


def test_scoring_before_fitting_says_so_rather_than_failing_inside_sklearn():
    """`score_table` on an unfitted classifier used to reach `None.predict_proba`."""
    import pandas as pd
    import pytest

    from splicescope.cryptic import FEATURE_COLUMNS
    from splicescope.ml import CrypticClassifier

    feats = pd.DataFrame({c: [1.0, 2.0] for c in FEATURE_COLUMNS})
    clf = CrypticClassifier()
    for call in (clf.predict_proba, clf.score_table):
        with pytest.raises(RuntimeError, match="call fit"):
            call(feats)


def test_the_score_table_orders_its_ties_the_same_way_every_time():
    """A forest that is certain gives many junctions exactly the same score, and sorting
    on the score alone left those rows in whatever order they arrived in. The top of this
    table is what a reader looks at and what the tutorial prints, so an order that moves
    between runs — or between one machine and another — is one they cannot cite.

    The fixture gives ten junctions one feature vector and ten another, so each block of
    ten scores identically and only the tie-break decides their order.
    """
    import numpy as np
    import pandas as pd

    from splicescope.cryptic import FEATURE_COLUMNS
    from splicescope.ml import CrypticClassifier

    rng = np.random.default_rng(0)
    cryptic = {c: 2.0 for c in FEATURE_COLUMNS}
    noise = {c: -2.0 for c in FEATURE_COLUMNS}
    feats = pd.DataFrame([cryptic] * 10 + [noise] * 10)
    feats["chrom"] = ["chr1"] * 20
    feats["start"] = rng.permutation(np.arange(1000, 1000 + 20 * 10, 10))
    feats["end"] = feats["start"] + 100
    feats["strand"] = ["+"] * 20
    feats["is_cryptic_truth"] = [1] * 10 + [0] * 10

    clf = CrypticClassifier(random_state=0).fit(feats)
    ordered = clf.score_table(feats)
    assert ordered["cryptic_score"].value_counts().max() >= 10, "the case this is about"

    shuffled = clf.score_table(feats.sample(frac=1.0, random_state=7).reset_index(drop=True))
    assert list(ordered["start"]) == list(shuffled["start"]), (
        "the same junctions in the same order, whichever order they were given in"
    )
    # and within one score, the coordinates rise
    for _, block in ordered.groupby("cryptic_score", sort=False):
        assert list(block["start"]) == sorted(block["start"]), "ties break on coordinates"
