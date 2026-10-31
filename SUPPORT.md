# Support and governance

## Getting help

- **Questions about using splicescope**, including how to run it on your own data:
  open an issue at <https://github.com/elmasnuryilmaz/splicescope/issues> and include the
  command you ran, the version (`splicescope --version`) and, if you can, a small input
  that shows the problem. The built-in simulator (`splicescope simulate`) is often enough
  to make a shareable example.
- **What a number means**: the definitions of Ψ, the count-based test, the consequence
  classes and the classifier are in [`docs/METHODS.md`](docs/METHODS.md); how every
  documented result is reproduced is in [`docs/REPRODUCING.md`](docs/REPRODUCING.md).
- **Bugs and feature requests**: use the issue templates on GitHub. A bug report is most
  useful with the full error message and the operating system and Python version.

Please do not report security-sensitive problems in a public issue; contact the
maintainer by email instead (the address is in `pyproject.toml` and `CITATION.cff`).

## What to expect

splicescope is maintained by one person alongside research work, so issues and pull
requests are answered on a best-effort basis. Confirmed bugs that change results are
fixed with priority over new features. A fix that changes any number in the
documentation is released with a changelog entry that says which numbers moved and why.

## Governance

- **Maintainer.** The project is led by its author, Elmasnur Yılmaz, who reviews and
  merges changes and makes releases.
- **How decisions are made.** Changes to methods, defaults or thresholds are proposed in
  an issue or pull request with the evidence that motivates them (see the pull request
  template), discussed in public, and recorded in [`CHANGELOG.md`](CHANGELOG.md).
- **Contributors.** Anyone may open issues and pull requests; the process is described in
  [`CONTRIBUTING.md`](CONTRIBUTING.md), and everyone taking part is expected to follow
  the [Code of Conduct](CODE_OF_CONDUCT.md). Regular contributors may be invited to
  become co-maintainers with merge rights.
- **Releases and archiving.** Releases follow semantic versioning, are tagged on GitHub,
  published to PyPI and archived on Zenodo, so that every released version stays citable
  and installable even if development stops.
