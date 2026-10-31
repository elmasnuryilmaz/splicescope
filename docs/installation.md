# Installation

`splicescope` needs Python 3.10 or newer and runs on Linux, macOS and Windows. Its
runtime dependencies are NumPy, pandas, SciPy, Matplotlib and scikit-learn, and pip
installs them automatically.

## From PyPI

Once the package has been released on PyPI:

```bash
pip install splicescope
```

Optional extras add the pieces that not every user needs:

```bash
pip install "splicescope[app]"    # the Streamlit dashboard
```

## From source

Until the first PyPI release, or to work on the code, install from GitHub:

```bash
pip install "git+https://github.com/elmasnuryilmaz/splicescope.git"
```

or clone it and install in editable mode, with the test and notebook tools:

```bash
git clone https://github.com/elmasnuryilmaz/splicescope.git
cd splicescope
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,docs]"     # add ",app" for the dashboard, ",site" for this website
```

| extra | what it adds | needed for |
|---|---|---|
| `app` | Streamlit | `streamlit run app/streamlit_app.py` |
| `dev` | pytest, pytest-cov, Hypothesis, ruff | running the test suite and the linter |
| `docs` | nbformat, nbclient, ipykernel, Pillow | rebuilding the tutorial notebook and the figures in `docs/` |
| `site` | MkDocs, Material for MkDocs, mkdocstrings | building this website with `mkdocs serve` |

## Check the installation

```bash
splicescope --version
splicescope simulate --outdir demo_data
```

From a source checkout, the full test suite is the most thorough check:

```bash
pytest
```

## Other tools you may need

- **An aligner that writes splice junctions.** `splicescope` reads STAR's `SJ.out.tab`.
  For alignments from another aligner, `validation/extract_sj.sh` derives the same table
  from any coordinate-sorted BAM file.
- **A FASTA index for the genome** (`.fai`, for example from `samtools faidx`), only for
  protein-consequence prediction.
- **Nextflow 23.04 or newer**, only for the [pipeline](nextflow.md).
