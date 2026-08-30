"""splicescope — detect, quantify and characterize splicing events from splice junctions.

A small, reproducible toolkit built around a single idea: most of what we want to know
about alternative and *cryptic* splicing can be recovered from the splice junctions that
aligners such as STAR already report, plus a reference annotation.

Public API
----------
- :func:`splicescope.io.read_star_sj`      – read STAR ``SJ.out.tab`` files
- :func:`splicescope.io.read_gtf_junctions` – derive known junctions from a GTF
- :func:`splicescope.annotate.annotate_junctions` – classify observed junctions
- :func:`splicescope.quantify.compute_psi` – splice-site usage (Ψ)
- :func:`splicescope.diff.differential_splicing` – ΔΨ between conditions
- :func:`splicescope.events.detect_events` – SE / MXE / A5SS / A3SS events + event-level PSI
- :func:`splicescope.cryptic.extract_features` – features for cryptic-event calling
- :func:`splicescope.ml.CrypticClassifier` – learn to separate cryptic events from noise
- :func:`splicescope.enrich.enrich_differential` – pathway over-representation (ORA)
- :func:`splicescope.simulate.simulate_dataset` – synthetic ground-truth data
"""

from __future__ import annotations

__version__ = "0.5.0"

__all__ = ["__version__"]
