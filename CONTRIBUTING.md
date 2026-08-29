# Contributing

Thanks for your interest! This is a research-grade project; contributions and issues are
welcome.

## Development setup

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

## Before opening a PR

```bash
ruff check .      # lint (and `ruff format .` to auto-format)
pytest -q         # the suite must stay green
```

## Guidelines

- Keep functions pure and testable where possible; every new behaviour needs a test.
- New biological metrics should come with a one-paragraph docstring stating exactly what
  is computed (and its assumptions).
- The project must stay runnable with **no external downloads** — extend
  `splicescope.simulate` rather than depending on real datasets in tests.
- Public API changes go in the `CHANGELOG.md` under *Unreleased*.
