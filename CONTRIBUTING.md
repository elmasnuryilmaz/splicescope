# Contributing

Thanks for your interest! This is a research-grade project; contributions and issues are
welcome.

## Reporting problems and asking questions

Bugs, feature requests and questions all go to the
[issue tracker](https://github.com/elmasnuryilmaz/splicescope/issues); the bug and
feature templates ask for what is needed to act on them. [`SUPPORT.md`](SUPPORT.md) says
where to look first, what response to expect, and how decisions about the project are
made.

## Development setup

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,docs]"
```

## Before opening a PR

```bash
ruff check .      # lint (and `ruff format .` to auto-format)
pytest -q         # the suite must stay green
```

If you changed anything the analysis computes, rebuild what the repository ships with it:

```bash
python examples/_build_tutorial.py     # the executed notebook
python examples/generate_showcase.py   # docs/showcase.png
python docs/make_demo_gif.py           # docs/demo.gif
python docs/make_nmd_rule_gif.py       # docs/nmd_rule.gif
```

To preview the documentation website after changing a page, a docstring or the CLI:

```bash
pip install -e ".[site]"
mkdocs serve                   # http://127.0.0.1:8000, rebuilt on every save
mkdocs build --strict          # what CI runs: any broken link or docstring fails it
```

CI runs the notebook and fails if the committed one no longer prints what the code
prints. The figures it cannot check that way — font rendering differs between machines,
so their bytes would differ for no reason worth failing over — which is why they are
listed here.

Behaviour worth having a test is usually worth a line in
[`validation/mutation_survey.py`](validation/mutation_survey.py): the survey breaks the
code one line at a time and reports anything the suite fails to notice. It runs weekly in
CI and takes about an hour locally, so `--module <name>` is the usual way to run it.

[`docs/REPRODUCING.md`](docs/REPRODUCING.md) is the map of what is checked and by what —
read it before adding a number to the documentation, because most of them are pinned to
the thing they describe and a new one probably should be.

## Guidelines

- Keep functions pure and testable where possible; every new behaviour needs a test.
- New biological metrics should come with a one-paragraph docstring stating exactly what
  is computed (and its assumptions).
- The project must stay runnable with **no external downloads** — extend
  `splicescope.simulate` rather than depending on real datasets in tests.
- Public API changes go in the `CHANGELOG.md` under *Unreleased*.
