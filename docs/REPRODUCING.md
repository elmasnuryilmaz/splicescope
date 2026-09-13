# Reproducing this project

Everything in this repository that can be checked, and the command that checks it.

The point of this document is that you do not have to take any number here on trust. Each
section names a claim, the command that reproduces it, roughly how long that takes, and
what you should see. Where something *cannot* be reproduced without downloading data, it
says so plainly rather than pretending otherwise.

Read it in order the first time. After that, the tables are the index.

---

## 1. Setting up

```bash
git clone https://github.com/elmasnuryilmaz/splicescope.git
cd splicescope
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,docs]"        # add ",app" for the Streamlit dashboard
```

`dev` brings the test suite, the linter and `hypothesis`; `docs` brings the notebook
tools, which the suite needs because it checks that the committed tutorial is the one its
builder produces. Nothing downloads biological data — the bundled dataset is simulated,
and that is deliberate: the whole pipeline has to run on a laptop with no network.

**The environment these numbers were measured on.** Results should not depend on it, and
the `oldest` CI job exists to prove the declared floors are real, but a difference in the
fourth decimal of a cross-validated metric usually traces back to here.

| | measured on | declared floor |
|---|---|---|
| Python | 3.11.5 | 3.10 |
| numpy | 2.4.6 | 1.23.5 |
| pandas | 3.0.5 | 1.5.0 |
| scipy | 1.17.1 | 1.9.3 |
| scikit-learn | 1.9.0 | 1.2.0 |
| matplotlib | 3.11.1 | 3.6.0 |
| hypothesis | 6.168.0 | 6.80.0 |

The floors are pinned in [`constraints-oldest.txt`](../constraints-oldest.txt) and
installed by a CI job, because a floor nobody installs is a promise rather than a fact.

---

## 2. The evidence, layer by layer

Seven layers, weakest claim first. Each answers a different question, and none of them
substitutes for another.

| layer | question it answers | command | time |
|---|---|---|---|
| suite | does the code do what the tests say? | `pytest -q` | ~1 min |
| coverage | is any of it never run? | `pytest --cov=splicescope` | ~1 min |
| properties | does it hold on inputs nobody wrote by hand? | part of the suite | — |
| mutations | would the suite *notice* if the code broke? | `python validation/mutation_survey.py` | ~50 min |
| outputs | do the committed figures still come from this code? | see §2.5 | ~35 s |
| measurements | do the numbers in METHODS hold? | see §2.6 | ~15 s |
| real data | does it recover biology? | see §2.8 | needs downloads |

### 2.1 The test suite

```bash
pytest -q
```

**360 tests, about 61 seconds.** Twenty test files. The suite is the ordinary kind — one
behaviour per test, named for the behaviour rather than the function — with two
exceptions worth knowing about, described in §2.3 and §3.

### 2.2 Coverage

```bash
pytest --cov=splicescope --cov-report=term-missing
```

**98 % of statements**, 39 uncovered out of 2,128. CI fails below 90 %, which is a floor
rather than a target: it fires when a feature lands untested, not when coverage drifts by
a line. The uncovered statements are defensive early returns; `pytest --cov-report=term-missing`
prints them by line number if you want to judge that for yourself.

Coverage is what found five whole features that could not run at all — pathway enrichment
behind `--gene-sets`, the supervised classifier behind `--labels`, and three others. It
does not tell you whether the code is *right*, which is what the next two layers are for.

### 2.3 Property-based tests

[`tests/test_properties.py`](../tests/test_properties.py), **20 tests**, run as part of the
suite. Generators produce junction tables, annotations and DNA; the invariants have to hold
whatever comes out:

- Ψ exhausts its splice site — every site's values sum to 1
- reverse-complement is an involution
- `find_ptc` points only at a real in-frame stop, and never past an earlier one
- Benjamini–Hochberg commutes with permutation and never undercuts its p-values
- annotation does not depend on row order
- no event names a coordinate that was never observed

These see what the other layers cannot: neither coverage nor the mutation survey can
produce an input a person would not have written. One found a real defect on its first
run — an empty result whose columns it could not index — and sweeping the package for that
same shape turned up three more of it.

To search harder than the suite does:

```bash
pytest tests/test_properties.py -q -p hypothesis --hypothesis-seed=0
```

### 2.4 The mutation survey

```bash
python validation/mutation_survey.py                  # every mutation, ~50 min
python validation/mutation_survey.py --module enrich   # one module, ~4 min
```

A green suite says the tests pass, not that they would fail if the code broke. The survey
breaks it on purpose: **122 mutations**, each a mistake someone could plausibly make — an
inclusive comparison where it should be strict, a distance measured from the wrong end of
a codon, a correction applied to the wrong array. Each one changes a line, runs the whole
suite, and reports anything the suite fails to notice.

Expected output ends with:

```
122 mutations, 0 surprise(s)
```

Three mutations are *expected* to survive: they are equivalent, not uncaught, and the
script's docstring says which and why. A fourth outcome, `NOT APPLIED`, means a mutation's
target text has moved — a rule that has quietly stopped being checked. That is also a
surprise, and it has happened.

Two things about the survey are worth copying if you build one:

- **It refuses to start on a suite that is already failing.** A mutation is judged caught
  by the suite failing, so a red suite reports every mutant as caught and says everything
  is fine. That happened here once; the only sign was two mutants documented as equivalent
  turning up as caught.
- **It puts the sources back on SIGTERM, not only on Ctrl-C.** Otherwise a job runner
  killing it leaves a planted bug in the working tree with nothing to say where it came
  from.

### 2.5 The committed outputs

The repository ships an executed notebook and four figures. They are a claim about what
the current code does, so rebuilding them has to mean something:

```bash
python examples/_build_tutorial.py        # examples/tutorial.ipynb   ~11 s
python examples/generate_showcase.py      # docs/showcase.png          ~9 s
python docs/make_demo_gif.py              # docs/demo.gif             ~12 s
python docs/make_nmd_rule_gif.py          # docs/nmd_rule.gif          ~2 s
git status --short                        # expect: nothing
```

**Every rebuild is byte-identical.** The notebook has numbered cell ids and no execution
timestamps, so a rebuild that changed nothing produces the same file and `git diff` shows
the analysis moving or nothing at all. Before that, every rebuild rewrote 24 random cell
ids and four wall-clock stamps per cell, and a real change arrived as four lines inside two
hundred lines of churn — which is how the tutorial spent time reporting 43 significant
junctions where the code found 46.

Two checks keep that true, and they are different:

```bash
pytest tests/test_documented_numbers.py -k committed_tutorial   # the cells are the builder's
python examples/_build_tutorial.py --check                      # the output is the code's
```

The first compares the notebook's *cells* against its builder, which catches a notebook
hand-edited in Jupyter. The second executes it and compares what it *printed*, which
catches a change to the analysis. CI runs the second on every push. Only text is compared:
the embedded figures are PNGs whose bytes depend on the platform's font rendering, so a
laptop and a CI runner differ for no reason worth failing over. Counts are compared
exactly and every other number to two significant figures — enough to let a different
scikit-learn build through, and still enough to report a p-value moving from 1.2e-53 to
2.4e-53.

### 2.6 Measurements on simulated data

Four scripts measure something stated in [`docs/METHODS.md`](METHODS.md), on data they
generate themselves. Each prints a table; a test compares every figure in METHODS against
the script that produced it, so the two cannot drift apart.

```bash
python validation/dispersion_trade.py        # METHODS §5.4    ~3 s
python validation/invariant_units.py         # METHODS §5.5    ~2 s
python validation/ora_bias_calibration.py    # METHODS §7b     ~1 s
python validation/consequence_corrections.py # what 0.9.0 changed
```

Each takes `--seeds N` to tighten the Monte-Carlo error. What they show:

- **`dispersion_trade`** — one shared dispersion is well calibrated when dispersion is
  homogeneous (0.048 against a nominal 0.05) and anti-conservative when it is not: with
  half the units at `s = 200` and half at `s = 5`, the loose half runs at 0.155 while the
  aggregate 0.077 hides it. `dispersion="per_unit_floor"` trades about a sixth of the
  power for protection against that.
- **`invariant_units`** — a splice site with nothing to choose between (a constitutive
  donor carrying one junction, Ψ ≡ 1) used to enter the shared dispersion estimate and
  drive it from a true 50 to 842, taking the false-positive rate from 0.049 to **0.173**
  against a nominal 0.05. The script measures it both ways, including the old arithmetic,
  so the "before" row is a measurement and not a remembered number.
- **`ora_bias_calibration`** — gene-level over-representation on splicing hits is unusable
  uncorrected: 99.4 % of biologically null gene sets reject at p ≤ 0.05, because a long
  gene has more chances to contain a significant junction. The opportunity weighting takes
  that to 2.1 %.

### 2.7 The manuscript numbers

```bash
python validation/manuscript_numbers.py      # ~6 s
```

Reproduces every quantity in the NMD-prediction manuscript, in section order, from the
three per-event tables committed in [`validation/manuscript_tables/`](../validation/manuscript_tables/):
the consequence calls for the cassette exons, the calls for the splice-site shifts, and
the LeafCutter2 labels from the full-transcriptome run. Bootstraps use
`numpy.random.default_rng(0)` with 2,000 resamples, so the confidence intervals it prints
are exactly the reported ones.

Figures: `figure_nmd_prediction.py` (Figure 1), `figure_method_comparison.py` (Figure 2),
`figure_rank_test_floor.py` (Supplementary Figure S1). The first two take `--table` and
`--out`.

### 2.8 The real-data run

This is the layer you cannot run from a clean checkout. It needs the SRA runs, a GENCODE
GTF and the GRCh38 primary assembly — all public, none small.

[`validation/README.md`](../validation/README.md) has the run accessions, the
BAM-to-junction step (`extract_sj.sh`, which derives a STAR-shaped `SJ.out.tab` from any
coordinate-sorted BAM), the exact `splicescope run` invocation, and what it showed.

**Read its caveat before quoting any number from it.** Those figures were measured with
0.8.1, and three defects in the differential path have been fixed since — each changing
*which* units are tested or how wide the null is. They will move when it is re-run. What
does not move is the argument they were collected for: a rank test cannot clear
Benjamini–Hochberg at 3 against 3, and the count model clears it by more than a hundred
orders of magnitude.

---

## 3. What continuous integration runs

Every push, in [`.github/workflows/ci.yml`](../.github/workflows/ci.yml):

| job | what it does |
|---|---|
| `test` | lint, suite and coverage floor on Python 3.10, 3.11, 3.12 and 3.13, then a CLI smoke run |
| `oldest` | installs `constraints-oldest.txt` and runs the suite, so the declared floors are real |
| `dashboard` | renders the Streamlit page headlessly |
| `tutorial` | executes the notebook and compares what it printed |
| `nextflow` | the pipeline on all three of its input paths |

Weekly, in [`.github/workflows/mutation.yml`](../.github/workflows/mutation.yml): the
whole mutation survey.

---

## 4. How a written claim is kept true

Numbers in documentation rot faster than anything else, and nothing noticed: the README
claimed the code was broken 85 ways when the survey held 111 mutations, 97 % coverage
against 98 %, and fourteen property-based tests against twenty.

Twenty-eight tests — under a tenth of the suite — now compare a written claim against the
thing it describes. They live in
[`tests/test_documented_numbers.py`](../tests/test_documented_numbers.py) and
[`tests/test_packaging_and_inputs.py`](../tests/test_packaging_and_inputs.py):

| what is pinned | to what |
|---|---|
| the README's test, mutation and property counts | the suite, the survey, the file |
| the README's coverage floor | `--cov-fail-under` in the workflow |
| the README's count of excluded shift candidates | the dataset `splicescope simulate` writes |
| the README's Python examples | they are executed |
| the README's citation version and DOI | `CITATION.cff` |
| the sentence the README quotes from the dashboard | the sentence the page opens with |
| METHODS §3's taxonomy rules | re-derived from the annotation |
| METHODS §5.4, §5.5 and §7b tables | the scripts that print them |
| METHODS §5b's PSI formulas | re-derived from the junction counts |
| the thresholds the CLI prints | `significant`'s defaults |
| the committed tutorial | its builder |
| the NMD animation's verdicts | `_nmd_from_downstream` |
| the package's advertised public API | the modules |
| every mutation's target | the source it mutates |
| the Python floor | `pyproject.toml`, the CI matrix, the oldest job, the README badge |
| the pipeline's `min_reads` and version | the CLI's default and the package's version |
| what CI installs | what the suite imports |
| every import | `pyproject.toml` |

The last two exist because CI was red for seven pushes and nobody looked: `hypothesis` was
imported at module level and declared nowhere, and `nbformat` was declared but not
installed in the job that needed it. Both were invisible locally, which is the whole
problem — a dependency you already have is invisible until somebody else installs the
package.

---

## 5. What auditing found, and how to find it again

Three defects came out of one technique: **re-derive a computation from its own
definition, in a second implementation, and compare.** Not a second call into the same
code — a second reading of the rule, written out where the test can be read.

Eight computations were re-derived. Six agreed; two did not.

| computation | outcome |
|---|---|
| protein-consequence calls | **found a defect** |
| the beta-binomial's fitted group means | **found a defect** |
| event PSI, all four types | agreed on 1,368 values |
| donor-usage Ψ | agreed on 2,491 values |
| the rank test's ΔΨ and sample counts | agreed on 515 units |
| the junction taxonomy | agreed on 1,063 junctions |
| Benjamini–Hochberg | agreed on 200 random vectors |
| the ORA hypergeometric tail | agreed on 218 rows |

Three times the *re-derivation* was the thing that was wrong, and the comparison found
that instead. That is the method working, not failing.

### The three defects

**A shared dispersion estimate that counted what it excluded.** Dispersion is measured
from how far replicates fall from their group's mean. A group with every read on one
junction has residuals of exactly zero and says nothing about replicate-to-replicate
variability — and most such groups are structural rather than biological: a constitutive
donor carrying a single junction has Ψ = 1 in every sample because the site has nothing
else to splice to, and those sites outnumber the alternative ones across a genome. The
estimator zeroed their residual and counted their observations in the degrees of freedom
anyway. Estimated precision went from a true 50 to 842, and the false-positive rate from
0.049 to 0.173 against a nominal 0.05. Reproduce with `validation/invariant_units.py`;
described in METHODS §5.5.

**A stop codon spanning the splice junction, never examined.** A ribosome reads the mature
mRNA straight through, so a codon can begin in the last one or two bases of a cryptic exon
and finish in the next one. The search looked at the exon and then at the downstream
sequence separately, and the downstream search starts at the frame the exon leaves behind
— which is exactly past that codon. Neither half examined it. Across eight simulator
configurations it is the first premature stop for 11 of 719 cassette predictions, and on
one of those the verdict flips from `in_frame_insertion` — the tool saying the protein
simply gains three amino acids — to `ptc_nmd`. The correction can only move events
*towards* decay, because the codon it finds is always earlier than the one reported
before. Pinned by three tests in `tests/test_consequence.py`.

**A junction listed twice for one sample.** The Ψ denominator sums every row at the splice
site, so a repeated junction counts the site twice and Ψ becomes its share of a total that
is not the site's. The figures still plot, the tests still pass, and the number is wrong.
The simulator made them: its deceptive-noise branch runs from an intron's known donor to a
novel acceptor, which is the shape of a cryptic exon's upstream junction, and on 8 of 30
seeds the two coincided exactly — leaving one junction labelled both cryptic and noise.
`compute_psi` now refuses such a table by name, and the simulator no longer writes one.

### Doing it yourself

Pick a computation whose rule is written down. Implement the rule again, in the test, from
the documentation rather than from the code. Run both over whatever the simulator emits.
Print the disagreements and look at each one by hand — the first few will be your own
arithmetic, and that is the point at which you learn the rule properly.
`tests/test_events.py::test_every_event_psi_is_recomputed_from_the_junction_counts` and
`tests/test_consequence_pipeline.py::test_every_cassette_call_is_recomputed_from_the_annotation`
are worked examples.

---

## 6. What is left

Four things, in the order they matter. None of them is code.

**Reboot the Streamlit app.** <https://splicescope.streamlit.app> is still serving an
older build — Streamlit Community Cloud has not redeployed on its own across many pushes.
The README's *Live demo* badge points there, so a reader currently meets a weaker demo
than the repository contains: the current page walks one cassette exon end to end and says
in words what including it does to the protein. One click from the Streamlit dashboard.

**Cut a release.** `CHANGELOG.md` records 0.9.0 and 0.9.1 as dated releases, but neither
is tagged and neither has a GitHub release; the newest real release is v0.8.1. That is
also why the citation is awkward: a version DOI is minted by archiving a release, so
`CITATION.cff` lists one for v0.8.1 and the README has to say there is none for this
version. Tagging closes both. The `Unreleased` section is large and behaviour-changing —
p-values on count-based runs move — so a minor bump rather than a patch.

**Decide two defaults.** Both are off, both are measured, and both are judgment calls
rather than bugs:

- `dispersion="per_unit_floor"` — costs about a sixth of the power, removes roughly half
  the false positives when dispersion is heterogeneous (METHODS §5.4).
- `filter_invariant=True` — sets aside units that cannot be rejected at any threshold
  before correcting. On a 200-gene simulation that is 48 % of the tested units, worth a
  factor of 1.91 in q (METHODS §5.5).

**Re-measure the real-data numbers.** `validation/README.md` carries counts from 0.8.1 and
three fixes have landed since. Needs the BAMs. Two specific things are outstanding there:
the third consequence correction (exon-skipping junctions) needs the GENCODE annotation to
identify, and the junction-level counts have not been re-run at all.

### Known and deliberate

Not to-dos — recorded so nobody has to rediscover them:

- **The truncation path does not catch the junction-spanning codon.** The same reasoning
  as the cassette fix applies and the same fix would, but the simulator emits one
  truncation per six hundred junction changes — and that one at a codon boundary, where no
  codon can straddle — so there is nothing to validate a fix against. METHODS §8b.
- **The simulator produces no `novel_combination` junction.** That is the model being
  faithful: a cryptic exon is absent from the annotation, so the isoform skipping it *is*
  the annotated intron. `novel_combination` comes from skipping an exon the annotation
  already has, which is alternative splicing of a known cassette rather than a cryptic
  event. `annotate` classifies both; only the simulator is narrower. METHODS §8.
- **The classifier fills a missing feature with 0.0**, and zero is a meaningful value for
  several of them — `dist_known_donor = 0` reads as *exactly on an annotated splice site*.
  It was measured rather than argued about: imputing the median moves cross-validated
  ROC-AUC from 0.900 to 0.898, and filling the distances with the largest distance seen
  gives 0.850 against 0.853 with 63 % of junctions on an unannotated contig. Unchanged
  deliberately, with the numbers at the line and a test pinning the limitation.
- **Intron retention is out of reach.** It is measured from coverage *inside* the intron,
  which splice junctions do not carry.
