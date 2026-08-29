"""Build and execute examples/tutorial.ipynb (reproducible, outputs embedded).

    python examples/_build_tutorial.py

Constructs the tutorial notebook cell-by-cell with nbformat, runs it in the
project kernel so every figure and table is captured, then writes the .ipynb.
"""

from __future__ import annotations

from pathlib import Path

import nbformat
from nbclient import NotebookClient
from nbformat.v4 import new_code_cell, new_markdown_cell, new_notebook

OUT = Path(__file__).resolve().parent / "tutorial.ipynb"

CELLS: list[tuple[str, str]] = [
    ("md", """\
# splicescope — a guided tour

This notebook walks through the whole splicescope pipeline end to end, on a small
**synthetic** dataset so it runs anywhere with no downloads. We'll go from raw splice
junctions to:

1. a **junction taxonomy** (annotated / novel / cryptic),
2. **splice-site usage** (Ψ),
3. **event-level PSI** for the four rMATS-style event classes (SE / MXE / A5SS / A3SS),
4. **differential splicing** between two conditions, and
5. a cross-validated **classifier** that prioritises genuine cryptic events over noise.

The last section shows how to point the exact same code at your own STAR `SJ.out.tab`
files instead of the simulator.

> **Running on Google Colab?** Run this first, then restart the runtime:
> `!pip install git+https://github.com/elmasnuryilmaz/splicescope`
"""),
    ("code", """\
import pandas as pd

from splicescope import annotate, cryptic, diff, events, plotting, quantify
from splicescope.ml import CrypticClassifier
from splicescope.simulate import simulate_dataset

# activate the inline backend *after* importing plotting (which sets Agg), so
# figures are captured in the notebook
%matplotlib inline
import matplotlib.pyplot as plt

pd.set_option("display.max_columns", 12)
"""),
    ("md", """\
## 1. A dataset with known ground truth

`simulate_dataset` builds a tiny annotated genome and emits per-sample junctions that
mimic STAR: canonical introns, true cryptic-exon events (up-regulated in condition **B**),
alternative 5′/3′ splice sites, mutually exclusive exons, and noise. `label_noise` flips a
few ground-truth labels to make the ML task realistically hard.
"""),
    ("code", """\
ds = simulate_dataset(
    n_genes=20, n_per_group=6,
    cryptic_fraction=0.6, alt_ss_fraction=0.4, mxe_fraction=0.3,
    label_noise=0.08, seed=1,
)
n_junc = ds.observed[["chrom", "start", "end", "strand"]].drop_duplicates().shape[0]
print(f"{ds.observed['sample'].nunique()} samples, {n_junc} unique junctions")
ds.observed.head()
"""),
    ("md", """\
## 2. Classify every junction

Each observed junction is labelled by how it relates to the reference: `annotated`,
`novel_donor`, `novel_acceptor`, `novel_combination`, or `cryptic` (both splice sites
unannotated — the strongest cryptic prior).
"""),
    ("code", """\
annotated = annotate.annotate_junctions(ds.observed, ds.known)
summary = annotate.annotation_summary(annotated)
summary
"""),
    ("code", """\
fig, ax = plt.subplots(figsize=(5.2, 3.2))
plotting.plot_annotation_summary(summary, ax=ax)
plt.show()
"""),
    ("md", """\
## 3. Splice-site usage (Ψ)

Ψ is model-free: for a junction it is the fraction of its splice site's reads that flow
through it. Thinly covered sites become `NaN` rather than a noisy ratio.
"""),
    ("code", """\
psi = quantify.compute_psi(annotated, min_reads=5)
psi.loc[psi["sclass"] != "annotated",
        ["chrom", "start", "end", "sclass", "sample", "count", "psi_donor"]].head()
"""),
    ("md", """\
## 4. Event-level PSI — SE / MXE / A5SS / A3SS

Beyond usage, splicescope reconstructs canonical splicing *events* from junctions and
computes rMATS-style percent-spliced-in for each. It covers four of the five rMATS event
classes (intron retention needs read coverage and is out of scope).
"""),
    ("code", """\
evs = events.detect_events(annotated)
evs["event_type"].value_counts()
"""),
    ("code", """\
fig, ax = plt.subplots(figsize=(5.2, 3.2))
plotting.plot_event_summary(evs, ax=ax)
plt.show()
"""),
    ("md", """\
## 5. Differential splicing between conditions

For each junction (and each event) we compare per-sample Ψ between conditions with a
Mann–Whitney U test and Benjamini–Hochberg FDR. Positive ΔΨ means more inclusion in B.
"""),
    ("code", """\
dsplice = diff.differential_splicing(psi, ds.groups)
print(len(diff.significant(dsplice)), "significant junctions (q<=0.05, |ΔΨ|>=0.1)")

fig, ax = plt.subplots(figsize=(5.2, 4))
plotting.plot_volcano(dsplice, ax=ax)
plt.show()
"""),
    ("code", """\
# the same differential test, at the event level
epsi = events.event_psi(annotated, evs, min_reads=5)
ediff = diff.differential_splicing(epsi, ds.groups, value="psi", key=["event_id"])
ediff.head(6)[["event_id", "event_type", "delta_psi", "qvalue"]]
"""),
    ("code", """\
fig, ax = plt.subplots(figsize=(5.2, 4))
plotting.plot_event_volcano(ediff, ax=ax)
plt.show()
"""),
    ("md", """\
## 6. Prioritising cryptic events with a classifier

Novel junctions are dominated by noise. Rather than hand-tuning thresholds, we train a
RandomForest on engineered features and evaluate it with **stratified cross-validation**
(no leakage). Every model ships with a model card.
"""),
    ("code", """\
feats = cryptic.extract_features(psi, ds.known)
clf = CrypticClassifier(random_state=0)
metrics = clf.evaluate(feats)
clf.fit(feats)
metrics
"""),
    ("code", """\
fig, ax = plt.subplots(figsize=(6, 3.2))
plotting.plot_importance(clf.importances, ax=ax)
plt.show()
"""),
    ("code", """\
from sklearn.model_selection import StratifiedKFold, cross_val_predict

x, y = clf._xy(feats)
oof = cross_val_predict(
    clf._make_pipeline(), x, y,
    cv=StratifiedKFold(5, shuffle=True, random_state=0), method="predict_proba",
)[:, 1]

fig, ax = plt.subplots(figsize=(4.2, 4.2))
plotting.plot_roc(y, oof, ax=ax)
plt.show()
"""),
    ("code", """\
# highest-scoring candidate cryptic junctions
clf.score_table(feats).head(8)[
    ["chrom", "start", "end", "sclass", "cryptic_score"]
]
"""),
    ("md", """\
## 7. Using your own data

Everything above works on real STAR output — just replace the simulator with the readers.
Given a folder of `*.SJ.out.tab` files, a reference GTF, and a `sample,condition` table:

```python
from splicescope import io

observed = io.read_many_star_sj({
    "ctrl1": "sj/ctrl1.SJ.out.tab",
    "ctrl2": "sj/ctrl2.SJ.out.tab",
    "kd1":   "sj/kd1.SJ.out.tab",
    "kd2":   "sj/kd2.SJ.out.tab",
})
known  = io.read_gtf_junctions("annotation.gtf")
groups = {"ctrl1": "ctrl", "ctrl2": "ctrl", "kd1": "kd", "kd2": "kd"}

annotated = annotate.annotate_junctions(observed, known)
psi       = quantify.compute_psi(annotated)
result    = diff.differential_splicing(psi, groups)
```

Or from the command line:

```bash
splicescope run --sj-dir sj/ --gtf annotation.gtf --groups groups.tsv --outdir results
```

That's the whole tour — from junctions to events to cryptic calls, reproducibly.
"""),
]


def build() -> nbformat.NotebookNode:
    nb = new_notebook()
    for kind, src in CELLS:
        nb.cells.append(new_markdown_cell(src) if kind == "md" else new_code_cell(src))
    return nb


def main() -> None:
    nb = build()
    print(f"executing {len(nb.cells)} cells ...")
    NotebookClient(nb, timeout=600, kernel_name="splicescope-venv").execute()
    # keep the file portable: a generic kernelspec, no venv-specific name
    nb.metadata["kernelspec"] = {
        "display_name": "Python 3",
        "language": "python",
        "name": "python3",
    }
    nbformat.write(nb, OUT)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
