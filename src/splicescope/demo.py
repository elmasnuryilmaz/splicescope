"""The whole pipeline on synthetic data, in one call.

This is what the Streamlit dashboard runs, and it lives here rather than in the app so
that it is exercised by the test suite. An interactive demo is the first thing a reader
tries, and its analysis had no test at all: a change to any signature it touches would
have shown up as a traceback on a public page rather than as a red build.

    from splicescope.demo import run_demo
    result = run_demo(n_genes=20, cryptic_fraction=0.6)
    result.differential.head()
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass
class DemoResult:
    """Everything the dashboard displays, and why any of it is missing."""

    annotation_summary: pd.DataFrame
    differential: pd.DataFrame
    features: pd.DataFrame
    #: ``None`` when the classifier could not be trained — see :attr:`classifier_note`.
    metrics: dict | None = None
    importances: pd.DataFrame | None = None
    scores: pd.DataFrame | None = None
    #: Why the classifier was skipped, in a sentence fit to show a reader.
    classifier_note: str | None = None
    #: SE / MXE / A5SS / A3SS events, and their differential inclusion between groups.
    events: pd.DataFrame | None = None
    #: Per-sample PSI for every event, which is what a reader has to see to believe a call.
    event_psi: pd.DataFrame | None = None
    event_differential: pd.DataFrame | None = None
    #: ``{sample: condition}``, so a caller can split the PSI by group.
    groups: dict[str, str] | None = None
    #: What each cassette exon does to the protein it lands in: frame, PTC, NMD.
    consequences: pd.DataFrame | None = None


def run_demo(
    n_genes: int = 20,
    n_per_group: int = 6,
    cryptic_fraction: float = 0.6,
    label_noise: float = 0.12,
    seed: int = 11,
    min_reads: int = 5,
    n_estimators: int = 300,
    with_consequences: bool = True,
) -> DemoResult:
    """Simulate a dataset and run annotation, Ψ, the differential test and the classifier.

    The classifier is skipped, rather than raising, when the simulated data leaves every
    junction in one class: a small number of genes with a low cryptic fraction and no
    label noise is a perfectly reasonable thing to ask for, and the junction classes and
    differential results are still worth showing. :attr:`DemoResult.classifier_note` then
    says why the rest is absent.

    ``with_consequences`` also detects events and predicts what each cassette exon does
    to the protein — the part that answers *so what?* about a cryptic junction. It needs
    a genome, so the simulated dataset is written to a temporary directory and read back
    through the same reader a real run uses.

    ``n_estimators`` is exposed only so that a test sweeping many parameter combinations
    does not have to train a full forest for each; the dashboard leaves it alone.
    """
    from .annotate import annotate_junctions, annotation_summary
    from .cryptic import extract_features
    from .diff import differential_splicing
    from .ml import TRUTH_COLUMN, CrypticClassifier
    from .quantify import compute_psi
    from .simulate import simulate_dataset

    ds = simulate_dataset(
        n_genes=n_genes,
        n_per_group=n_per_group,
        cryptic_fraction=cryptic_fraction,
        label_noise=label_noise,
        seed=int(seed),
    )
    annotated = annotate_junctions(ds.observed, ds.known)
    psi = compute_psi(annotated, min_reads=min_reads)
    result = DemoResult(
        annotation_summary=annotation_summary(annotated),
        differential=differential_splicing(psi, ds.groups),
        features=extract_features(psi, ds.known),
    )

    labels = result.features.get(TRUTH_COLUMN)
    if labels is None or labels.nunique() < 2:
        only = "cryptic" if labels is not None and labels.eq(1).all() else "not cryptic"
        result.classifier_note = (
            f"Every one of the {len(result.features)} novel junctions here is "
            f"'{only}', so there are no two classes to separate. Raise the fraction of "
            "genes with a cryptic event, or add a little label noise."
        )
        return result

    clf = CrypticClassifier(random_state=0, n_estimators=n_estimators).fit(result.features)
    result.metrics = clf.evaluate(result.features)
    result.importances = clf.importances
    result.scores = clf.score_table(result.features)
    if with_consequences:
        _add_events_and_consequences(result, ds, annotated, min_reads, seed)
    return result


def _add_events_and_consequences(result, ds, annotated, min_reads: int, seed: int) -> None:
    """Detect events, test them, and predict what the cassette exons do to the protein."""
    import tempfile

    from .consequence import GenomeFasta, annotate_consequences, load_transcripts
    from .diff import differential_splicing
    from .events import detect_events, event_psi
    from .simulate import write_dataset

    events = detect_events(annotated)
    if events.empty:
        return
    result.events = events
    result.groups = dict(ds.groups)
    psi = event_psi(annotated, events, min_reads=min_reads)
    result.event_psi = psi
    result.event_differential = differential_splicing(
        psi, ds.groups, value="psi", key=["event_id"]
    )

    cassettes = events[events["event_type"] == "SE"]
    if cassettes.empty:
        return
    with tempfile.TemporaryDirectory() as tmp:
        outdir = write_dataset(ds, tmp, seed=seed)
        transcripts = load_transcripts(outdir / "annotation.gtf")
        with GenomeFasta(outdir / "genome.fa") as fasta:
            result.consequences = annotate_consequences(
                cassettes.copy(), transcripts, fasta,
                start_col="exon_start", end_col="exon_end",
            )
