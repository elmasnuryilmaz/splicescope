# Streamlit dashboard

`app/streamlit_app.py` is an interactive page that answers one question end to end: a
knockdown switches on junctions the annotation does not contain, so which of them are
real, and what do they do to the protein?

It picks a cryptic cassette exon, shows its Ψ in every replicate, and states in words what
including it does to the reading frame, whether a premature stop codon appears and
whether nonsense-mediated decay is predicted, next to the reason a rank test could not
have called it.

## Run it locally

```bash
pip install "splicescope[app]"          # or: pip install -e ".[app]" in a checkout
streamlit run app/streamlit_app.py
```

The app needs the `app/` directory from the repository, so run it from a clone. A hosted
copy is linked from the README.

## How it relates to the package

The page is presentation only. The analysis behind it is
[`splicescope.demo.run_demo`](api.md#splicescope.demo.run_demo), so the test suite
exercises exactly what the page shows, and a separate CI job renders the page headlessly
with Streamlit's test harness on every push (`tests/test_dashboard.py`).
