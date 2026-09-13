# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/) and this project adheres to
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- **The executed tutorial and both figures in `docs/` are regenerated.** They are built
  from the simulator, whose draw order changed when an intron stopped being allowed to
  hold two events, so their committed outputs no longer matched what the code produces —
  the notebook reported 43 significant junctions where the current code finds 46. The
  junction *count* was identical either way (198 both times, coincidentally), which is
  why nothing looked wrong. `examples/tutorial.ipynb`, `docs/showcase.png` and
  `docs/demo.gif` are rebuilt; the analysis behind them is unchanged.
- **The weighted ORA is 14× faster, and gives the same answers.** `p/(1-p)` does not
  depend on which gene set is being tested, but it was being recomputed for every
  background gene for every set: one division per (gene × set), 95.7 million of them on a
  human-sized run. It is now computed once per run, and each set's *outside* sum comes
  from the total by complement, so a set costs its own size rather than the whole
  background. Measured on a 20 000-gene background with the hits drawn in proportion to
  opportunity:

  | gene sets | before | after |
  |---|---|---|
  | 5 000 (≈ MSigDB H + C2) | 10.8 s | 0.9 s |
  | 16 000 (≈ MSigDB C5) | 35.1 s | 2.5 s |

  The plain hypergeometric takes 1.5 s on the same input, so the correction now costs
  1.8× the uncorrected test instead of 24×. This matters because `weight_by_units=True`
  is the default. Across 480 gene sets spanning 1 gene to all-but-one of the background,
  the largest disagreement with the old arithmetic is 3.5e-14 relative — the two differ
  only in floating-point noise. The smaller side is always the one summed outright, so
  the subtraction never has to recover a small number from two near-equal ones.

- **Annotation resolves each junction once, not once per sample.** A junction's class
  and gene depend only on its coordinates, but they were being recomputed for every row
  of every sample — on a six-sample human-sized table four fifths of the work was
  repeats. Reading the values was the other half of it: pandas 3 backs string columns
  with Arrow, and pulling a million values out of one an element at a time cost more
  than the classification did, so each coordinate column is now converted in one call.
  The same two-line change applies in `compute_psi` and in the event detectors' site
  helper, so the fast spelling is the only one left in the codebase.

  | 998 000 rows over 6 samples | before | after |
  |---|---|---|
  | `annotate_junctions` | 3.6 s | 2.1 s |
  | whole pipeline through events | 19.8 s | 17.1 s |

  At 20 samples the annotation step is 2.3× faster, since that is where the repeats are.
  Checked by output rather than by tests alone: the annotated frame is identical — same
  columns, same index, same values — on a human-scale table, a 20-sample one and the
  simulator's, and every one of the six TSVs a full CLI run writes is byte-identical to
  what the previous code produced.

### Added
- **A test that the event-level comparison uses the read counts.**
  `differential_splicing(value="psi")` resolves its count columns by name, and the names
  it looks for are the ones `event_psi` happens to produce. Nothing checked that they
  still match, and the failure is silent: with the counts unreachable the test falls back
  to Mann-Whitney. At six replicates a side that is nearly invisible — 21 events called
  instead of 22 — but the q-values collapse by twenty-three orders of magnitude, and at
  the three-a-side design METHODS §5.1 is about, it is the difference between **32 of 44
  events called and none at all**.

- **Property tests for the day's new code.** The guards and the sentence generator are
  the least battle-tested things here, so they get generators rather than examples:
  `describe` must be a sentence for any combination of fields a row can carry and must
  never reach its catch-all; the annotation check must object to a naming mismatch and to
  nothing else, checked against every pair of chromosome sets; and the genome check must
  accept any contig long enough, however much longer. The first of them found the bug
  below on its first run.

- **`docs/nmd_rule.gif`: the 50-nucleotide rule, animated.** A premature stop codon
  walks towards the last exon-exon junction and the prediction flips exactly once, where
  the rule says it does. That threshold is the whole of what this toolkit predicts about
  a cryptic exon, and it is one number — worth a picture more than anything else on the
  page. Every verdict in it comes from `_nmd_from_downstream`, the function the pipeline
  uses, rather than from a drawing of the rule, and a test walks the same positions and
  requires the call to change once and to straddle 50. The animation cannot drift from
  the code.

- **The README's Python examples are run by the suite, and the diagram shows the
  consequence step.** A snippet that has stopped working is the first thing a reader
  meets, and nothing touched them. Both had something wrong: one assigned `evaluate`'s
  metrics dict to a variable called `clf`, so the obvious next line — `clf.fit(...)` —
  would have raised on a dict; and the flowchart of how the pipeline fits together left
  out the protein-consequence step entirely, which is the layer the demo is now built
  around. There is also a second example showing the whole pipeline, consequences
  included, in one call.

- **A `py.typed` marker, so the annotations are usable.** Sixty of the sixty-five
  public functions carried type hints and none of them was visible to anyone installing
  this: under PEP 561 a type checker ignores a package's inline annotations unless the
  package ships that marker. It is now in the wheel — verified by building one — and
  declared as package data so it stays there. The remaining five functions, all in
  `plotting`, are annotated too: a partly annotated package is worse than an unannotated
  one, because the checker trusts what is there and infers `Any` for the rest.

- **The showcase figure shows what the events do to the protein.** Six panels, and two
  of them were the classifier — which METHODS §7.1 is explicit is the weakest part of
  this toolkit, worth an order of magnitude less than the statistic when measured against
  replication in an independent experiment. The protein consequence, which is the part
  nothing else in the figure says, was absent. It now replaces the permutation-importance
  panel, so the front page shows a premature stop and whether it is far enough from the
  last junction to trigger decay. The generator writes the simulated data out and reads
  it back through the same reader a real run uses, so the figure comes from the code a
  user would run rather than from a shortcut.

- **The dashboard answers a question instead of redrawing charts.** It used to be a
  simulator with sliders: move one, the data is regenerated, the charts change. That
  demonstrates that the code runs, not what it is for — and the tool's most distinctive
  layer, the protein consequence, was not on the page at all.

  It now walks one cryptic cassette exon from junctions to protein. Its Ψ in every
  replicate, control against knockdown, so the effect is visible rather than summarised.
  Then what including it does, in a sentence: *"Including these 61 nucleotides shifts the
  reading frame, and the first premature stop codon appears 63 nucleotides in. It sits
  355 nucleotides before the last exon-exon junction, more than the 50 the rule allows,
  so the ribosome stops while the junction complex is still downstream and
  nonsense-mediated decay is predicted to degrade the transcript."* Beside it, why a rank
  test could not have called it at all.

  Two pieces of that are in the package rather than the page, so they are tested:
  `consequence.describe` turns a prediction into that sentence for every one of the eight
  classes, and `plotting.plot_event_psi` draws the replicates. `run_demo` now also
  detects events, tests them, and predicts consequences, so the page is presentation only.

- **Fourteen property-based tests, and `hypothesis` in the `dev` extra.** The README
  described the suite as "unit + property + end-to-end". Nineteen tests had
  property-shaped names, but every one asserted an invariant on a hand-picked input —
  there was no generator anywhere, so the claim described a style rather than a
  technique. That technique sees what the others cannot: neither the mutation survey nor
  coverage can produce an input a person would not have written. The invariants now
  checked against generated junction tables, annotations and DNA include Ψ exhausting
  its splice site, reverse-complement being an involution, `find_ptc` pointing only at a
  real in-frame stop and never past an earlier one, BH commuting with permutation,
  annotation not depending on row order, filling the omitted zeros not moving a junction
  that was observed, swapping the group labels flipping ΔΨ while leaving the p-value
  alone, a locus on another chromosome changing nothing, and no event naming a coordinate
  that was never observed. Searched at 600 examples per property as a one-off; the
  committed settings run 40, in about four seconds.

- **`--labels`, which makes the classifier reachable from the command line.** The `run`
  command already had a supervised-classifier branch, guarded on a truth column — and
  nothing an aligner writes carries one, so `model_card.json` and `cryptic_scores.tsv`
  could never be produced from the command line. Fifteen lines of unreachable code in the
  main entry point, while the model card told every reader to retrain on curated labels
  before real-data use and the tool offered no way to do it. `--labels` takes a TSV of
  `chrom, start, end, strand, is_cryptic_truth`, merges it onto the Ψ table, and reports
  how many rows matched. A file missing columns names them; a file that matches nothing is
  an error rather than a classifier trained on nothing.
- **The Markdown model card carries the hyper-parameters too.** The JSON one gained them
  from the fitted estimator this release; the human-readable card, which is the one a
  reviewer reads, still omitted them.
- **Coverage is measured, with a floor in CI.** 97 % of statements, and the `test` job
  fails below 90 %. Five of the blocks that were uncovered turned out to be whole
  features with no test: pathway enrichment behind `--gene-sets`, the classifier branch
  above, the entire rank-sum test, the `simulate` subcommand — the README's first
  command, exercised by CI's smoke test but pinned by nothing — and
  `dispersion="per_unit_floor"`, which is the public option METHODS §5.4 exists to
  justify. `cli.py` went from 80 % to 98 % and `diff.py` from 78 % to 95 %. What is left
  uncovered is scattered defensive lines, not features.

- **`validation/dispersion_trade.py`, and the §5.4 tables checked against it.** METHODS
  §5.4 is the argument for leaving `dispersion="shared"` as the default: it costs roughly
  a sixth of the power to halve the false-positive rate on loosely dispersed units. Those
  figures had no reproduction path at all — no script, and no pointer to one, unlike §7b
  — and the suite's own tests only check the direction, with bounds loose enough for any
  of them to drift. There is now a script that states every parameter and averages over
  five independent simulations, and a test that compares every cell of both tables
  against it. The test also caught that the two tables state some of the same figures,
  and now requires them to agree.

- **The ORA calibration table in METHODS is checked against the script that produces
  it.** `validation/ora_bias_calibration.py` gained a `measure()` function, and a test
  runs it and compares every figure in the §7b table at the precision the table quotes.
  It found the table stale on its first run, which is how it earned its place.

- **`constraints-oldest.txt` and a CI job that installs it.** Every dependency has a
  declared floor in `pyproject.toml`, and the CI matrix resolves to the newest of
  everything on every Python it tests, so nothing ever checked that the floor was real —
  a promise to anyone installing into an existing environment. It is real: the whole
  suite passes on numpy 1.23.5, pandas 1.5.0, scipy 1.9.3, scikit-learn 1.2.0 and
  matplotlib 3.6.0, and equally on pandas 2.3.3 and pandas 3.0.5. (scipy is pinned one
  patch above its floor because 1.9.0 ships no wheel for several platforms.) A test keeps
  the pins in step with the declared specifiers in both directions: raise a floor above
  its pin, or add a dependency without one, and it says which.

- **A `docs` extra, so the executed tutorial can actually be re-executed.** Rebuilding
  `examples/tutorial.ipynb` and the figures in `docs/` needs `nbformat`, `nbclient`,
  `ipykernel` and `pillow`, and none of them was declared anywhere. A notebook advertised
  as reproducible needs the tools that reproduce it to be installable: `pip install -e
  ".[docs]"`.
- **The version check covers the README's citation.** It compared `__version__` against
  `CITATION.cff` and the CHANGELOG but not against the citation a reader copies, which is
  why that one sat at v0.8.1 through two releases.

- **`splicescope.demo.run_demo`, and the dashboard rendered in CI.** The README's first
  badge is a live Streamlit demo, which makes a broken dashboard the most visible failure
  this repository has — and its analysis had no test at all, so a change to any signature
  it touched would have appeared as a traceback on a public page rather than as a red
  build. The analysis now lives in the package as one function, the page is presentation
  only, and both are covered: `tests/test_demo.py` sweeps the sixteen corners of the
  slider ranges the page exposes, and `tests/test_dashboard.py` renders the real page
  through Streamlit's own headless harness. A `dashboard` job in CI installs the `app`
  extra and runs it. `run_demo` is also the shortest way to run the whole pipeline on
  synthetic data from Python.

- **Sixteen tests for the consequence layer on the minus strand.** Every coordinate in
  that module means something different depending on the strand: donor and acceptor swap
  ends, exons are walked in decreasing genomic order, coding length accumulates
  downwards, and each piece of sequence is reverse-complemented on the way out. Half of
  any genome is on that strand, and it was hardly tested — `test_consequence_pipeline.py`
  never used it at all. Of sixteen deliberate defects in the module's structural code,
  five passed the whole suite and three were reachable only from the minus strand: exon
  blocks left in genomic rather than transcription order, the coding length before a
  position counted one base too many, and the GTF phase added to the reading frame
  instead of subtracted from it (which changes the frame, and so the stop codon, for
  every 5'-incomplete CDS).

  The main new test is a mirror. One gene is built on the plus strand and a second is its
  exact reverse complement: the chromosome is reverse-complemented, every coordinate
  reflected, the strand flipped. The two therefore describe the same transcript, the same
  reading frame and the same event, so every field of every prediction has to agree —
  whatever the right answer happens to be. Five splice-site shifts and three cassette
  exons are compared this way, plus a check that the mirror is not agreeing merely
  because both sides gave up.

- **Thirteen tests for the features the classifier learns from.** Extending the
  mutation survey to `cryptic.py` found **seven of eight** deliberate defects passing the
  whole suite: the distance to the nearest known splice site searched only forwards from
  the insertion point instead of both ways (turning a 50 nt shift into a 149 nt one),
  only donors indexed as known sites, the intron length off by one, measured zeros
  counted as supporting samples, the canonical-motif flag inverted, the minus-strand
  canonical motif `CT/AC` dropped, and a junction's label taken as the minimum over its
  rows rather than the maximum. A wrong feature does not make a model fail, it makes it
  learn the wrong thing and report a good score for doing so. Each is now pinned in
  `tests/test_cryptic_features.py` against one hand-built input whose expected values are
  worked out from the coordinates rather than copied from the output.

  `ml.py` gave up three more: `predict_proba` taking the column for *noise* instead of
  the one for *cryptic*, which inverts every score a user reads while the score table's
  own ordering stays self-consistent and `evaluate` reaches the probabilities by another
  route; asking for more cross-validation folds than the rare class has members; and the
  class weighting being dropped.

- **The simulator reports the events it injected.** `SimulatedDataset.truth` is a table
  of every injected cryptic exon, A5SS, A3SS and MXE with its host intron and exon
  coordinates, and `write_dataset` writes it as `truth.tsv` beside the SJ files. Recall
  was previously measured against a set reconstructed from the coordinates, which gives
  every geometry the simulator *could* have used rather than the ones it did — a gene
  drawn for an event is skipped when its intron cannot hold one. That is why the MXE test
  asserted `>= 55 of 60` and said so in a comment. It now asserts all of them.

- **`validation/mutation_survey.py`, and a weekly CI job that runs it.** The survey
  that found the gaps below is now a script, so the claim is reproducible rather than
  asserted — the same shape as the other scripts in `validation/`. It breaks the code 44
  ways, runs the suite against each, and exits non-zero on any surprise: a defect nothing
  failed on, or one of the two documented equivalent mutants suddenly being caught. 42
  are caught, 2 survive by construction, 0 surprises. Too slow for every push at one full
  suite run per mutation, so CI runs it on Mondays and on request.

- **Eleven tests for the rules the suite was not checking.** Found by mutation: change
  a line, run the suite, see whether anything fails. Thirty-seven deliberate defects were
  introduced across the eight modules surveyed at that point, and **thirteen of them
  passed all 165 tests**. Two are
  harmless by construction (a slice shorter than three characters can never equal a stop
  codon; a negative chi-square statistic yields the same p-value of 1) and the other
  eleven are now each pinned by a test. Four of them were outside the consequence and
  statistics modules:

  - **The coverage boundary.** A splice site carrying exactly `min_reads` reads is
    measured, not withheld, and making the comparison inclusive broke nothing. This
    filter was already inert on the default code path until 0.9.0.
  - **Which site a measured zero comes from.** A junction gets its zero where *either*
    of its splice sites has reads in that sample; filling only from the donor loses every
    junction whose acceptor was the covered end, and on the minus strand the other half.
  - **STAR's strand codes.** Column 4 of `SJ.out.tab` is 0, 1 or 2, and 0 means STAR
    could not tell. Reading 0 as `+` invents a strand for every non-canonical junction
    and leaves `resolve_unstranded` nothing to resolve. No test had read a real
    `SJ.out.tab` and checked the mapping.
  - **The weighted tail.** The plain branch had a test for summing from `k-1` rather
    than `k`; the weighted branch did not. Wallenius' distribution at odds 1 *is* the
    hypergeometric, so giving every gene the same opportunity now has to reproduce the
    uncorrected p-value exactly.

  The seven in the consequence and statistics modules: Two mutants turned out to be harmless (a shorter slice can never equal a
  three-character stop codon; a negative chi-square statistic gives the same p-value of
  1), and the rest were real gaps, now closed:

  - The **50-nucleotide rule**, which is the tool's central scientific claim, had nothing
    pinning it. Making the threshold inclusive, measuring the distance from the start of
    the stop codon instead of its end, or counting the last exon's length into the
    distance all left 165 tests green. One sequence is now tested twice with the
    downstream exons split one base apart, so the distance is exactly 50 and then exactly
    51 and the call has to flip there and nowhere else.
  - The **last-exon exception** could be deleted without a single failure. A stop in the
    final exon has no junction downstream to be upstream of, so the distance does not
    exist, and the test now checks that it is reported as absent rather than as a
    negative number.
  - **Samples in neither group.** 0.9.0 fixed a likelihood ratio that charged ungrouped
    samples entirely to the null and moved a p-value from 1.6e-09 to 9.7e-93; the fix
    shipped without a regression test. The test now adds six samples whose Psi sits at
    0.02 and 0.98 to an unchanged 3-vs-3 comparison and requires every returned value to
    be unmoved.
  - **The q-value.** Nothing verified that `differential_splicing` applies
    Benjamini-Hochberg at all — passing the raw p-values straight through as q-values
    was invisible to the suite, on which the FDR claim in every output table rests.
  - **The precision bounds.** `s` is a ratio whose denominator is an excess variance, so
    data noisier than the beta-binomial drives it below 1 and can take it negative, which
    makes the likelihood meaningless. Psi drawn from {0.02, 0.98} at 400x coverage lands
    exactly on the floor, so the clip is load-bearing.

- **Two tests for things the event suite could not see.** Every hand-built event test
  uses one chromosome, so a grouping key that forgot the chromosome would pool junctions
  from unrelated genes and assemble events out of the pieces — and nothing in the suite
  failed when that key was removed, from either the cassette or the MXE detector.
  Coordinates repeat across chromosomes, so this is not a hypothetical. There is now a
  test whose two halves are each deliberately incomplete and can only produce an event
  if they are merged. A second test injects events into 400 genes on both strands that
  share donors and acceptors the way real genes do, and compares the detected event
  *identities* against the injected ones — 200 skipped exons, 60 alternative 5' sites and
  210 alternative 3' sites, with nothing missed and nothing invented. Mutation-checked:
  dropping the chromosome from any of the three grouping keys now fails a test, where
  before it failed none.

### Fixed
- **Rebuilding the tutorial is now a meaningful check.** `examples/tutorial.ipynb`
  ships its executed outputs, which are a claim about what the current code does — but
  the claim was not checkable, because every rebuild rewrote 24 random cell ids and four
  wall-clock timestamps per cell. A real change arrived as a handful of lines inside two
  hundred lines of churn, which is how the tutorial spent time reporting 43 significant
  junctions where the code found 46 with nothing looking wrong. Cell ids are now
  numbered and the execution metadata is stripped, so two rebuilds of unchanged code
  produce a byte-identical file and `git diff` shows the analysis moving or nothing at
  all. A test checks that the committed notebook is still the builder's — a notebook
  hand-edited in Jupyter and committed would otherwise be discarded silently by the next
  rebuild.
- **A splice site with nothing to choose between was tightening the null for every other
  site.** The most consequential defect found by this audit, and it was inside the
  dispersion estimator. Most splice sites in a genome are constitutive and carry a single
  junction, so Ψ there is 1.0 in every sample — not because splicing is precise, but
  because the site has nothing else to splice to. Such a group's residuals are zero by
  construction and it says nothing about replicate-to-replicate variability.
  `estimate_precision` zeroed its residual through a `variance > 0` guard and then
  counted its observations in the degrees of freedom anyway, so it entered the
  denominator and not the numerator. (The guard could never have fired: it is applied
  after the fitted mean has been clipped away from both ends, which makes the variance
  positive by construction.) That dilutes the residual mean towards zero and drives the
  estimated precision up, narrowing the null for every real unit:

  | units that cannot vary | estimated `s` (true 50) | false positives at nominal 0.05 |
  |---|---:|---:|
  | none present | 50.6 | 0.049 |
  | one per unit that can, counted as before | 842.1 | **0.173** |
  | one per unit that can, excluded | 50.6 | 0.049 |

  A human annotation is far more lopsided than one-to-one, so this understates it. Such
  groups are now excluded from both sides of the estimator, which removes the effect
  entirely. The exclusion is *per group*, not per unit, so a junction switched fully on in
  the knockdown and varying in the control still contributes what the control's
  replicates do — and that is the shape of a real cryptic event. **p-values on
  count-based runs move, by roughly a factor of two on the tutorial's data; the committed
  outputs are rebuilt.** Reproduce with `validation/invariant_units.py`; METHODS §5.5 and
  a test compare every figure.
- **A design can pass the replication check while the data hold nothing to estimate
  from.** `dispersion_is_estimable` is given no counts, so it can only weigh informative
  samples against fitted group means. A table of constitutive splice sites passes it with
  replicates to spare — every group is fully covered — and still leaves nothing to measure
  dispersion from, because Ψ is 1 in every sample. The estimate then returned its ceiling,
  asserting *no* overdispersion and narrowing the test to a binomial one, in silence: the
  exact failure the design check exists to catch, reached past it. `dispersion_has_
  information` asks the question of the counts, and `differential_splicing` warns on it.
  It deliberately does not fire on data genuinely tighter than binomial, which reach the
  same ceiling — that is a measurement, not the absence of one.
- **The three dispersion routines now share one block of bookkeeping.** They each
  computed the fitted group means, the coverage and the parameter count for themselves,
  which is how the defect above survived: the residual sum excluded groups pinned at Ψ 0
  or 1 and the degrees of freedom counted them, in two copies of what should have been
  one piece of code.
- **The per-unit estimator had the same inconsistency, where it costs protection.** A
  unit switched fully on in the knockdown and varying in the control is the shape of a
  real cryptic event. Counting the knockdown group's zero residuals diluted the control's
  genuine looseness, so a unit whose replicates vary at `s = 5` was estimated at 13.6
  instead of 5.7 — and since `dispersion="per_unit_floor"` takes the *smaller* of the
  shared and per-unit values, the floor stopped biting on precisely the units it exists
  for. The `dispersion_trade.py` figures are unchanged, because that simulation has no
  group at a boundary.
- **`differential_splicing(filter_invariant=True)`, off by default.** The same units also
  enter the Benjamini–Hochberg denominator, where they make every real unit's q-value
  worse though their own p-value is exactly 1 and no threshold could reject them: 48 % of
  the tested units on a 200-gene simulation, costing a factor of 1.91 in q. Setting them
  aside is the degenerate case of independent filtering (Bourgon, Gentleman & Huber, PNAS
  2010) — the criterion is the spread across all samples, which never looks at the group
  labels, and zero spread forces p = 1 whatever those labels are. Off by default because
  it changes every q-value in a run; a warning reports how many units it set aside.
- **An event table's columns were an accident of the data.** Sweeping for the defect
  above found two more of it, in `events`. `detect_events` returned whatever `pd.concat`
  made of the types it happened to find: a dataset with cassette exons but no mutually
  exclusive ones came back without the MXE columns, and one with no events at all came
  back with five columns where a full result has thirty. So `events["exon_start"]` raised
  a `KeyError` or did not depending on the input — and the showcase script's own first
  move, select the SE rows and hand `exon_start` to the consequence layer, was one
  eventless dataset away from raising. The schema is now pinned to the event types that
  were *asked for*, which is something a caller can reason about, and `event_columns()`
  reports it. `event_psi` had the same split: the full schema when there were no events,
  a column-less frame when there were events but no sample carrying them.
- **`events.tsv` wrote an exon boundary as `1840.0`.** A table holding more than one
  event type has a gap wherever a column belongs to a different type, and float64 is the
  only NumPy dtype that can hold one, so every coordinate in a mixed table picked up a
  decimal point. A genomic position with a `.0` reads as a rounded measurement, and a
  tool downstream parsing the column as an integer fails on it. Coordinates are now a
  nullable integer dtype: same values, same gaps, and `1840` in the file. Read counts
  stay floating point, because an inclusion count is the mean of two junctions and 7.5
  is a real value.
- **A consequence run on an alternative-splice-site table no longer answers a question
  it was not asked.** Pinning the event schema above had a consequence of its own, caught
  before it shipped: `splicescope consequence` picked cassette-exon mode by testing
  whether the `exon_start` column was *present*, and it is now always present. On a file
  holding only A5SS and A3SS events it would have found no cassette exon, written an
  empty `consequence.tsv`, and exited 0 — which reads as *these cassette exons change no
  protein*. Detection now tests whether the columns hold anything, and the refusal says
  how many rows the file had and of which types.
- **An enrichment result that found nothing had no columns.** Found by a property test,
  which could not index the frame it got back. `over_representation` skips a gene set
  that contains no hit, so a collection where none of them overlaps left the record list
  empty and `pd.DataFrame.from_records([])` returned a frame with no columns at all.
  Selecting `term` raised a `KeyError` on that frame and worked on every other one, and
  `to_csv` wrote a headerless file for a run that did test sets and enriched none — the
  ordinary outcome of asking an honest question. The schema is now a single constant,
  `enrich.RESULT_COLUMNS`, returned by both empty paths *and* imposed on the populated
  one: `qvalue` is computed after the records are built, so it used to land last there
  and second-to-last in the early return. An empty result and a full one are now the same
  table with a different number of rows. Same defect class as the `differential_splicing`
  fix below, in the other module that had it.
- **`describe` raised on a premature stop with no recorded distance.** Found by a
  property test over every combination of the fields a row can carry. The `ptc_escape`
  branch handled a missing distance — the last-exon case, where there is no junction
  downstream for the stop to be upstream of — and the `ptc_nmd` branch did not, so
  `int(None)` raised a `TypeError` from the middle of a sentence. A row reaches this
  function from a `consequence.tsv` read back with pandas as readily as from the
  pipeline, and an empty cell there is a NaN.
- **A coverage threshold nothing can meet says so.** `--min-reads 100000` printed "0
  significant junctions" and stopped, which reads as *no differential splicing* — a
  conclusion. Not one unit was testable: Ψ is NaN wherever a splice site carries fewer
  reads than the threshold, so raising it far enough empties the table with no other
  sign. `differential_splicing` now says when no unit had enough informative samples in
  both groups, which covers that cause and every other.
- **Library warnings reach a terminal as the CLI's own.** Three checks now live in the
  library, where everyone importing it is covered, and a terminal wants
  `warning: junctions: …` rather than a Python warning carrying a file and a line
  number. The CLI relays each one, labelled by the step that raised it, since the same
  sentence comes from the junction-level test and the event-level one.
- **Gene sets that name nothing tested say so.** With no identifier in common nothing
  can be tested, and the empty result reads as "no pathway is enriched" — a conclusion,
  where the truth is that the question was never asked. Gene-set files are keyed by
  symbols about as often as by accessions, and a GMT for the wrong organism looks exactly
  like one for the right one. The CLI has warned about this since 0.8.1;
  `enrich_differential` returned the empty frame in silence, so everyone using the
  library got it. The warning is now in the library, with both spellings in it, and the
  CLI relays that one message instead of composing a second.
- **A genome from another assembly is refused rather than read off the end of a contig.**
  `GenomeFasta.fetch` clamps to the contig's length and returns a truncated string, not
  an error, so the wrong genome does not fail — it finds no stop codon anywhere and calls
  everything a frameshift. Measured: swapping the simulator's genome for one 200 bases
  long turned **six `ptc_nmd` calls into zero** and reported the same nine events as
  confidently as before. The `.fai` already carries every contig's name and length, so
  both halves of the mistake are cheap to see, and `check_genome` now refuses a
  chromosome the FASTA does not have and a coordinate past the end of one it does. It
  covers only the chromosomes where a prediction will actually be made — those with both
  an event and a transcript to host it — so an event on an unannotated contig keeps its
  honest `no_host_transcript` answer and a genome larger than the analysis is nobody's
  mistake.
- **An annotation nothing can match is refused instead of reporting a genome of novel
  splicing.** Every junction class but `cryptic` is defined by agreeing with the
  annotation, so an annotation that matches nothing does not fail — it calls everything
  cryptic. That is the most exciting result this tool can produce and the easiest one to
  produce by accident. Two ways in, both silent until now:

  - a GTF carrying only gene or transcript records, which yields no introns at all,
    because introns are the gaps between exons;
  - a GTF and an aligner's output from different sources, since GENCODE writes `chr1`
    where Ensembl writes `1`.

  Both are now errors that name what is wrong, and the CLI prints them the way it prints
  its other refusals rather than as a traceback. Only those two: a junction on a contig
  the annotation does not cover — a scaffold, a decoy, a chromosome left out of a small
  analysis — is ordinary and still runs. The naming mismatch is distinguishable from it
  because stripping the `chr` prefix makes the two agree exactly.
- **`CONSEQUENCE_CLASSES` lived in a test file.** The set of things the consequence
  layer can say was defined in `tests/test_consequence_pipeline.py` and nowhere in the
  package, so anything exhaustive over the eight classes had no list to be exhaustive
  against. It is now in `consequence.py`, and the test imports it.
- **The live demo was over its memory limit, and the README pointed at a dead URL.**
  Two separate things. The badge and the personal site both linked to the
  hash-generated address the app had before it moved to a custom subdomain, so the link
  had been dead since the move; both now point at `splicescope.streamlit.app`. And the
  app itself was answering *"This app has gone over its resource limits — it's using too
  much memory"*.

  The page renders four figures and `st.pyplot` does **not** close them: its
  `clear_figure` default is `False`, and `pyplot` keeps every figure in a global
  registry until something else does. A Streamlit script reruns on each widget change,
  so the figures accumulate without bound. Measured against what does not grow:

  | | RSS |
  |---|---:|
  | bare Python | 8 MB |
  | + streamlit | 49 MB |
  | + matplotlib, pandas, numpy, scikit-learn, scipy | 214 MB |
  | + one analysis at the page's defaults | 255 MB |
  | + a second at the slider maximum | 255 MB |
  | leaked figures | +35 MB per 30 interactions, unbounded |

  Every figure is now drawn through a helper that closes it in a `finally`, and the
  cache is bounded (`max_entries`, `ttl`) rather than keeping one entry per slider
  combination for the life of the container.

  The test for this is structural — it parses the page and requires every `st.pyplot`
  call to sit in a function that also closes — because a test that counted open figures
  was written first and was useless: `AppTest` runs the script with its own module state,
  so the registry the test process can see is not the one the page uses, and it passed
  happily with the leak reinstated.
- **A differential result with no rows had no columns.** Found by the first property test
  to run, which could not index the frame it got back. `differential_splicing` returned a
  bare `pd.DataFrame()` when no unit met `min_samples`, so the empty case had a different
  shape from every other: selecting columns raised `KeyError` where the same code worked
  on a result with rows, and `to_csv` wrote a file with no header line at all — which is
  what the CLI produced for `differential_splicing.tsv` on a run that found nothing, and
  what anything downstream would then try to read. Both tests now return their own full
  schema, in their own column order, which differ from each other.
- **A sample listed in `--groups` with no file warned three times.** The CLI names it,
  and then hands the mapping to `differential_splicing`, which warns about the same
  sample again on each of its two calls. The CLI now narrows the mapping to the samples
  it actually read, after saying which it dropped, so one specific message replaces three.
- **METHODS §5.1 overstated the rank test's floor.** It says a two-sided Mann-Whitney on
  3 against 3 cannot return a p-value below `2/C(6,3)` = 0.1. That is the *exact* test's
  floor, and `scipy.stats.mannwhitneyu` computes the exact p-value only for a small,
  tie-free sample — with ties it uses the normal approximation, which returns less:

  | 3 vs 3, two-sided | p |
  |---|---:|
  | separated, untied — exact | 0.100 |
  | the same data, approximation forced | 0.081 |
  | one group holding tied values | 0.064 |
  | Ψ = 0 in every control, 1 in every knockdown | **0.047** |

  The last row is the signature the tool exists to find, and it lands under the nominal
  0.05. The section's conclusion survives — 40 switching junctions among 2 000 need the
  40th at p ≤ 0.001 — but the reason is Benjamini-Hochberg, not the floor, and a small
  enough experiment breaks it. Both the section and
  `min_achievable_rank_pvalue`'s docstring now say so, and a test pins each number.
- **The §5.4 dispersion figures did not reproduce.** Reconstructing the simulation from
  the parameters the section states gives 0.155 where the table said 0.173, 0.077 where
  it said 0.087, and a shared precision estimate of 10.5 rather than 10.7 — the power
  figures, 0.520 and 0.442, were right. The conclusion is unchanged, so the tables now
  carry what the script produces, and `differential_splicing`'s own docstring is updated
  with them.
- **METHODS §7b quoted a false-positive rate the code no longer produces.** The weighted
  ORA's rate on biologically null gene sets moved from 1.9 % to 2.1 % when the propensity
  estimate gained its Jeffreys smoothing in 0.9.1 — a change the release notes describe —
  but the table in §7b was not updated with it. The 0.9.1 entry that quotes the older
  figure now says which one the release ships.
- **The limitations section denied a feature the README advertises.** METHODS §9 said Ψ
  here is splice-site usage "not event-level PSI (cassette exon, A5SS/A3SS, IR)", while
  §5b specifies event-level PSI, `detect_events` returns all four classes and the
  README's second paragraph leads with them. The one class genuinely out of reach is
  intron retention, for the reason §5b already gives: it is defined by reads *inside* the
  intron, which junction counts do not carry. The same stale sentence was in the README.
- **The README's suggested citation named v0.8.1** while the package was 0.9.1.
- METHODS had a subsection numbered 7c sitting above section 7b; it is a subsection of 7
  and is now numbered as one, matching 5.1-5.4.
- **Four of the dashboard's sixteen slider corners crashed it.** Every one had the
  cryptic-event fraction at its minimum and the label noise at zero, which is a
  reasonable thing for a reader to ask for and which leaves every junction in one class.
  `np.bincount` cannot see a class that is entirely absent — for all-zero labels it
  returns a single bin — so the fold-count guard in `evaluate` passed, the forest fitted
  a one-class model, and `predict_proba` raised `IndexError: index 1 is out of bounds for
  axis 1 with size 1` from far away from the cause. The asymmetry hid it: all-positive
  labels happened to be caught, all-negative ones were not. Both `fit` and `evaluate`
  now refuse a single-class label vector in a sentence that names the problem, and the
  dashboard skips the classifier and says which slider to move instead of showing a
  traceback. The junction classes and the differential results are still drawn.
- **The model card could describe a model that was never fitted.** Its `class_weight` was
  the string `"balanced"` written into the card's own source rather than read from the
  estimator, so removing the weighting from the pipeline left the card still claiming it.
  A model card exists precisely so a reviewer need not take the method on trust, which
  makes this the one thing it must not do. `class_weight` is now a field on the
  classifier, the pipeline is built from it, and the card reports the fitted estimator's
  own parameters, falling back to the field before a fit. It also reports the fold count
  actually used, which `evaluate` lowers when the rare class has fewer members than
  `n_splits`.
- **The simulator no longer put two events in one intron.** It does not just crowd
  them, it changes what the reads mean. An MXE intron has its skipping junction
  suppressed, because mutually exclusive exons have no skipping isoform, so a cryptic
  exon placed in the same intron is not a detectable cassette; and an alternative donor
  there is also a leg of an MXE pair, so it is reported as that and excluded from
  alternative-site detection. Both readings are right, and both made the ground truth
  claim events that were not in the data — 5 of 37 cryptic exons and 1 of 24 A5SS sites
  on one seed. An intron now carries at most one event, and a gene whose chosen intron is
  taken picks another instead of being dropped. Recall against the reported truth is then
  **100 % for all four event types on every seed tested**.
- **`differential_splicing` refuses a sample/group mismatch instead of returning
  nothing.** If no sample in the Psi table is named in `groups`, every row is
  unassigned, every unit fails `min_samples`, and an empty table comes back: the caller
  reports "0 significant junctions" and believes it. The CLI has refused this since
  0.8.1, but only by comparing filenames, so everyone using the library got the silence.
  It is now a `ValueError` naming both lists. A partial mismatch warns and continues,
  since that is usually a real experiment with a sample dropped. A misspelled `value`
  column, a `key` column the table lacks, and a Psi table with no `sample` column now
  say so too, instead of surfacing as a `KeyError` from inside pandas.
- **A gene set with a member outside the background no longer raises `KeyError`.** Only
  genes that were tested can be drawn, so such a member is simply not in the set for the
  purpose of the odds; it was being looked up in the propensity table regardless.
- **Two version spellings of one accession have their opportunity added up.** A merged
  annotation can carry `ENSG….16` and `ENSG….17`; the background already counts them as
  one gene, but the weight table kept whichever came last, under-weighting exactly the
  genes that are split.
- Removed an unreachable fallback in `selection_propensity` — with a non-empty
  background every gene lands in a bin, so the code it guarded could never run.

## [0.9.1] — 2026-09-12

### Added
- **The Nextflow pipeline is executed in CI, both branches.** Its real-data mode was broken
  from the day it was written until 0.9.0 and nobody noticed, because nothing ever ran it;
  the 0.9.0 fix itself shipped with a caveat that it had not been executed either. It has
  now been run: demo mode, real-data mode with the inputs staged in three unrelated places
  as a real dataset has them, and real-data mode without `--genome` (which must skip the
  consequence layer and does). The pre-fix pipeline was run on the same inputs for
  comparison and fails exactly as predicted — `--sj-dir star_junctions/sj --gtf
  star_junctions/annotation.gtf`, then `error: no *.tab files`. A `nextflow` job in CI now
  covers all three.
- **Pathway ORA corrects the gene-opportunity bias.** The hypergeometric null treats every
  gene as one equally likely draw, but a gene contributes as many chances of being a hit as
  it has tested junctions, so long many-exon genes are over-represented among the hits for
  reasons that are not biological. Uncorrected this is catastrophic rather than mild: on
  480 biologically null gene sets whose hits were drawn *purely* in proportion to unit
  count, the plain test called **99.4 %** of them at p ≤ 0.05 and 54 % at p ≤ 1e-6
  (`validation/ora_bias_calibration.py` reproduces it). Genes
  are now binned by unit count, each bin's observed hit rate becomes its members'
  propensity, and the p-value comes from Wallenius' non-central hypergeometric with
  odds = (mean of `p/(1-p)` inside the set) / (mean outside) — `goseq`'s device applied to
  junction count instead of transcript length. The same simulation then gives **1.9 %** at
  p ≤ 0.05 and 0 % at p ≤ 1e-6 (2.1 % once the smoothing below is applied, which is what
  this release actually ships). Tied weights share a bin, so a gene's propensity cannot
  depend on its name. `bias_odds` reports the odds per set and
  `weight_by_units=False` restores the plain test, which Wallenius equals at odds 1.
  One deliberate departure from `goseq`: it averages probabilities where this averages
  odds. They agree while `p` is small, but a propensity here reaches 0.5, and at identical
  power (80 % detection of a set enriched by 0.10, 100 % above) averaging probabilities
  leaves **38–40 %** of the null sets called at p ≤ 0.05 against **0 %** for odds.
  Wallenius' `ω` is a ratio of sampling weights, which is what the odds are.

### Fixed
- **A saturated propensity bin let one gene rewrite the whole enrichment table.** Bins are
  filled in ascending weight order and the remainder becomes its own bin, which at the top
  of a heavy-tailed unit-count distribution is a handful of the highest-opportunity genes —
  exactly the genes that are almost always hits. Its observed rate was then 1.0, clipped to
  `1 − 1e-6`, and `_bias_odds` averages `p/(1-p)`, so each such gene contributed **1e6** and
  decided the mean by itself, in the numerator for the set that held it and the denominator
  for every set that did not. Adding one of those genes to a 100-gene set — which made the
  set *more* enriched, k 60 → 61 — moved its p-value from **3.9e-15 to 0.22**: a real
  pathway containing the longest gene in the background was reported as not enriched at
  all. It occurred in 6 of 40 simulated universes. Rates are now shrunk by bin size
  (`(hits + ½)/(n + 1)`) instead of clipped, so no bin can assert certainty: an all-hit bin
  of 7 gives 0.94, not 1. Calibration is unchanged (1.9 % → 2.1 % of null sets at p ≤ 0.05)
  and the same p-value now moves 6.8e-10 → 1.6e-07 rather than flipping.
- **`normalize_gene_id` collapsed whole gene families in a model organism.** 0.9.0 began
  stripping a trailing `.<digits>` from every identifier, on the docstring's own premise
  that "gene symbols do not end in `.<digits>`". *C. elegans* sequence names — the standard
  identifier for the majority of worm genes, which have no CGC name — are exactly that
  shape, so `C42D8.1`, `C42D8.2` and `C42D8.3` normalised to one key. The background size,
  the set size, the overlap, the fold enrichment, the p-value and the reported
  leading-edge gene list all shrank with them, silently and with exit 0. The suffix is now
  removed only when the stem looks like an Ensembl or RefSeq accession.
- **`run` crashed on digit-only sample names.** pandas reads a `sample` column of `101`,
  `102`, … as int64, which can never equal the strings derived from filenames, so the new
  mismatch check took its error branch on a perfectly consistent cohort — and then died
  with `TypeError: sequence item 0: expected str instance, int found` and exit 1 rather
  than printing its diagnostic. The column is now read as text, which also keeps `007`
  from becoming `7`.
- **MXE detection could grow as the fourth power of a locus's exons.** Removing the
  one-pair-per-anchor truncation in 0.9.0 fixed a correctness bug and introduced a
  scalability one: every `(upstream, downstream)` junction combination at an anchor is a
  candidate, so `n` real exons arrive with about `n**2` of them — the real exons plus
  every span from one exon's start to a later exon's end — and pairing all of those is
  quadratic again. **80 exons under one anchor produced 1,677,140 rows in 5.9 s.**
  Candidates are now ranked by the read support of their weaker flanking junction and the
  best `max_candidates` (default 50) are kept, with a warning naming the busiest anchor.
  Reachable as `detect_events(max_candidates=...)` and `run --max-mxe-candidates`.
  Two purely geometric rules were tried first and both **deleted 15–22 % of the injected
  events** along with the spans, because a noise junction sharing the anchor's donor is
  indistinguishable in shape from a real exon; read support distinguishes them. All 200
  injected events are recovered at every seed tested.
  Support alone is not enough either: in a *tandem array* every candidate is flanked by
  equally deep junctions, so the order under the cap was arbitrary and the cut took real
  exons with it — 11 densely packed exons kept only 6, and 15 kept 5. The shorter exon
  breaks the tie, since a span runs from one exon's start to a later exon's end and is
  always longer than the real exon sharing its start. Arrays of 8, 11, 15 and 25 exons now
  keep every one, the simulator still recovers 200/200, and the noisy anchor stays capped.
- **An unestimable dispersion was silently treated as no dispersion.** When no unit has
  more informative replicates than fitted group means there is no residual degree of
  freedom anywhere, and `estimate_precision` returned `max_precision` — asserting *no
  overdispersion* rather than admitting the design says nothing. The beta-binomial test
  then narrows to a plain binomial one: a **1-vs-1 comparison of Ψ 0.300 against 0.360 at
  1000 reads returns q = 4.5e-03**, on no replication whatsoever, with nothing in the
  output to say so. `dispersion_is_estimable` now reports whether the design can support
  an estimate and `differential_splicing` warns when it cannot. The default
  `min_samples=2` already avoided this; passing `min_samples=1` walked straight into it.

## [0.9.0] — 2026-09-11

### Added
- **The protein-consequence layer now runs on the built-in data, with no downloads.**
  `simulate` writes the chromosome its annotation describes: coding exons drawn from the
  61 sense codons so every gene carries a real open reading frame ending in a single stop,
  canonical `GT`/`AG` intron ends, CDS records in the GTF, and a `.fai` index. The 5′UTR is
  staggered by one base per gene so the coding frame at an intron boundary cycles through
  0, 1 and 2 — with the uniform 120 bp exons a fixed UTR put every intron on a codon
  boundary, and frame inheritance, the subtle half of the prediction, was never exercised.
  `simulate_genome` uses its own random stream, so no junction-level number moves.
- **Splice-site shifts are reachable from the command line.** `run --genome` now writes
  `junction_consequence.tsv` beside `consequence.tsv`, and `consequence` gained
  `--mode {auto,exon,junction}`, defaulting to reading the mode off the columns present.
  `annotate_junction_consequences` — the headline of 0.8.0, and already described in the
  README — existed only as a Python API and could not be invoked by a user.
- A junction already explained by a detected cassette or MXE event is excluded from the
  splice-site-shift pass. Read on its own, a cassette inclusion junction looks like an exon
  extension running to the end of the intron, which is the wrong reading of it; on the demo
  this removes 12 of 31 candidate junctions.
- `plot_consequence_summary`, and the `event_summary` / `event_volcano` figures that
  `plotting` had provided since 0.4.0 but `run` never called.
- 14 tests for the above, including that the simulated ORFs translate cleanly, that the
  `.fai` offsets address the FASTA correctly across line boundaries, that all three reading
  frames occur, and that cassette junctions are not double-reported as shifts.
- CI smoke-tests the consequence path end to end.

### Fixed
- **`detect_mxe_events` emitted one event per anchor pair and abandoned the rest.** Inside
  each shared (donor, acceptor) context the first non-overlapping exon pair set a flag and
  broke out of both loops, so a locus with three or more mutually exclusive exons — the
  interesting case — reported one pair and lost the others. The `seen` set already handled
  deduplication, so the breaks bought nothing. The loss was total rather than partial: the
  discarded exons' junctions share the emitted pair's donor and acceptor, so once
  `detect_events` marked those as used, each site had a single remaining alternative and
  A5SS/A3SS detection rejected them too. The exons disappeared from `events.tsv` with no
  warning. Worse than under-reporting, the survivor was the genomically *leftmost* pair
  rather than the best-supported one, so noise junctions displaced real events: on
  `simulate_dataset(n_genes=200, mxe_fraction=1.0)`, **31–44 of 200 emitted MXE rows
  carried wrong exon coordinates** while the injected pair sat unused in the same group
  (seeds 4, 6, 11). Every pair is now enumerated; all 200 injected events are recovered at
  every seed, and detection still takes 0.02 s.
- `max_exon` is reachable through `detect_events` and as `run --max-exon`; the 1 kb
  candidate-exon window was hard-coded for anyone using the documented entry points.
- **The `min_reads` coverage filter was inert for the beta-binomial test** — the default
  test, on every run since 0.7.0. `_wide` pivoted Ψ with `aggfunc="sum"`, and pandas sums
  an all-NaN group to `0.0`, so `psi_wide.notna()` was true for every cell that had a row
  at all and the NaN masking `compute_psi` performs was discarded. `--min-reads` changed
  nothing: on `simulate_dataset(n_genes=20, n_per_group=4, seed=3)` the test returned the
  same 105 units and 18 significant hits at `min_reads=10` and at `min_reads=100000`,
  where every Ψ is NaN. Consequences were not cosmetic — the BH denominator was inflated
  by every sub-threshold unit, so **every q-value was wrong**, and `estimate_precision`
  was fitting dispersion to 1–2 read residuals. Ψ is now pivoted with `mean`; counts are
  still summed. The rank-sum path was never affected because it uses `dropna`, which is
  why `validation/README.md` records 138,643 rank tests against 189,700 count tests on
  identical input: the 51,057 difference *was* the bug.
- **A junction absent from a sample is a measured zero, not missing data.** Aligners list
  only junctions with at least one read, so `compute_psi` never produced a row for a
  junction in a sample where it had none, the pivot left NaN, and the unit was dropped for
  want of informative replicates. The events layer already defaulted a missing count to
  zero, so the two halves of the tool disagreed about the same data. The casualties were
  exactly the cleanest cryptic events: a junction with 0 reads in all three controls and
  100 in all three knockdowns was **absent from `differential_splicing.tsv` entirely**,
  while leakier events were reported. `compute_psi` now completes the junction × sample
  grid (`add_unobserved_zeros`), adding a zero only where that junction's donor or
  acceptor site has coverage in the sample, so a site nobody sequenced stays unobserved.
  `fill_unobserved=False` restores the old behaviour.
- 8 regression tests; 5 of them fail against the previous code.

### Changed
- Every junction- and event-level p-value and q-value moves as a result of the two fixes
  above, in both directions: sub-threshold units leave the BH family, and junctions absent
  from a group enter it. The real-data numbers in `validation/README.md` were measured
  with 0.8.1 and are marked as awaiting a re-run; the qualitative conclusions there do not
  depend on the exact values. On simulated data the cryptic classifier moves from
  ROC-AUC 0.785 to 0.795 (AP 0.702 → 0.700), because `mean_psi_donor` now averages over
  the samples where a junction was measured absent as well.

### Also added
- Test coverage for the plotting module: threshold behaviour of both volcano
  plots, fixed event-type ordering with zero-fill, top-N truncation in the
  enrichment plot, and that `savefig` creates missing parent directories.
- `CODE_OF_CONDUCT.md`, issue templates (bug report, feature request) and a pull
  request template that requires stating the effect on published numbers.
- Python 3.13 to the CI matrix, and the corresponding trove classifiers.
- `Issues`, `Changelog`, `Documentation` and `Archive` (Zenodo DOI) entries under
  `[project.urls]`.
- Zenodo DOI badge, and the DOI in `CITATION.cff` (concept and version).
- Gzipped input: `read_gtf_junctions` and `read_gmt` now open `.gz` directly. GENCODE and
  Ensembl ship `.gtf.gz`, and the README's own example passes one, which until now failed
  with a `UnicodeDecodeError`.
- `io.sample_name_from_path`, and a test that `__version__`, `CITATION.cff` and the
  CHANGELOG agree — they had drifted to 0.5.0, 0.8.1 and 0.8.1 respectively, so
  `splicescope --version` reported a release three versions old.

### Changed (breaking)
- `events.cassette_psi` is removed; use `events.event_psi`, which keys by `event_id`.
  `cassette_psi` keyed its rows by the *skipping* junction, so several cassette exons
  between the same pair of flanking exons collapsed into one test — and since they can
  move in opposite directions, their ΔΨ cancelled. A worked example with three events at
  one locus produced a single junction key. Nothing in the pipeline used it; `run` has
  been on `event_psi` since 0.3.0.

### Added (statistics)
- `betabinom.estimate_precision_per_unit` and
  `differential_splicing(dispersion="per_unit_floor")`, which takes the smaller of each
  unit's own precision and the shared one. One shared dispersion is **not uniformly
  conservative**, contrary to what the docstring claimed: on a null with half the units at
  `s = 200` and half at `s = 5`, the shared estimate is 10.7 and the loosely dispersed half
  runs at a false-positive rate of **0.173** against a nominal 0.05 while the tight half
  runs at 0.000 — the pooled 0.087 hides both. The floor brings the loose half to 0.090 and
  the overall rate to 0.045, at a cost of roughly a sixth of the power (0.52 → 0.44 on a
  homogeneous Ψ 0.30-vs-0.45 contrast). **The default is unchanged**, so no published
  number moves; `docs/METHODS.md` §5.4 records the trade-off and the tests record the
  behaviour.
- More broadly: `predict_consequence` gained `contiguous_downstream`; `Transcript` gained
  `cds_phase` and `frame_at`; `enrich.normalize_gene_id` is public.

### Fixed (real data)
- **Strand-undefined junctions became phantom cryptic calls.** STAR writes strand code 0
  whenever the intron motif does not reveal a strand — routine for non-canonical junctions
  — and such a junction can never equal a stranded annotation. It was therefore classified
  `cryptic` however ordinary it was, lost its `gene_id`, and, being the only junction at
  its own `(chrom, position, ".")` site, was handed **Ψ ≡ 1.0 in every sample**: a row in
  the differential table that can never show a difference, inflating the BH denominator and
  the headline `cryptic` count alike. `resolve_unstranded` now places such a junction on
  the strand its annotation implies, by exact intron or by either splice site. Junctions
  the annotation cannot place keep `"."` — the strand really is unknown.
- **Alternative-splice-site events were deleted for sharing an anchor with a cassette.**
  Cassette and MXE junctions are excluded from A5SS/A3SS detection so the same signal is
  not told twice, but they were dropped before the alternatives at a site were counted, so
  a site with a genuine third alternative fell to one and its event vanished. A site is now
  skipped only when *every* junction at it is already explained.
- `splicescope consequence` raised `KeyError: 'consequence_class'` when filtering an event
  table to cassette exons left nothing; it now writes an empty result and exits 0.

### Fixed (statistics and biology)
- **Samples in neither group inflated the likelihood-ratio test without bound.** `lrt`
  fitted the null over every valid sample but the alternative over the two groups only,
  so the two hypotheses were not fitted to the same observations and the whole likelihood
  of the ungrouped samples was charged to the null. Adding six samples that belong to
  neither group moved an unchanged 3-vs-3 comparison from **p = 1.6e-09 to p = 9.7e-93**,
  with `n_a` and `n_b` still reported as 3 and 3. Samples outside both groups now take no
  part. This is not exotic: an SJ directory routinely holds more samples than a given
  contrast lists, which `run` now also warns about.
- **An in-frame exon truncation was reported as a premature stop.** Removing a whole
  number of codons leaves the downstream reading frame untouched, so the first in-frame
  stop found is the transcript's *own* — reported as `ptc_escape` with a PTC offset. On the
  demo genome, deletions of 3, 6, 9 and 30 nt all returned the identical offset 296, which
  is the annotated stop codon. A stop at or beyond where the protein natively ends is no
  longer called premature; these are `exon_truncation`, and frameshifting truncations still
  report their genuinely premature stops.
- **GTF phase was discarded, so every 5'-incomplete CDS was translated in the wrong
  frame.** GENCODE marks these `cds_start_NF` and gives the first CDS record a non-zero
  phase; `load_transcripts` ignored column 8 entirely. `Transcript.cds_phase` now carries
  it — taken from the first block in *transcription* order, so the minus strand is right —
  and `Transcript.frame_at` applies it wherever frame is inherited.
- **A PTC in an extension of the final exon was called `ptc_nmd`.** An extension is
  contiguous with the exon it joins, so there is no junction between them; when that exon
  is the last one, the stop lies past the final exon-exon junction and NMD cannot be
  triggered. The old code measured the distance to a junction that does not exist — on the
  demo genome it reported `ptc_nmd` at a distance of 120 for an event that escapes.
- **An exon-skipping junction was read as a splice-site shift.** Both its sites are
  annotated, just not as a pair, so `junction_change` matched one of them and reported the
  skipped exon *and the intron beyond it* as a single contiguous deletion. Such junctions
  are now left uninterpreted, which is what they are: a shift has exactly one novel site.

### Fixed (Nextflow)
- **The pipeline's real-data mode could not have worked.** One channel was passed for all
  three inputs, and the process script then looked for `<sj_dir>/sj`,
  `<sj_dir>/annotation.gtf` and `<sj_dir>/groups.tsv` — a layout no real dataset has. The
  three inputs are now staged independently, `--genome` is exposed (with its `.fai` staged
  alongside), `--min-reads` is passed through, and demo mode runs the consequence layer
  too now that `simulate` writes a genome. *Not executed locally: Nextflow is not
  installed on the development machine. The CLI invocation the process builds was
  verified against a replica of Nextflow's staging layout.*

### Fixed (enrichment)
- **Pathway ORA could never test a single gene set against a real annotation.**
  `over_representation` intersected identifiers verbatim, but a GTF gives versioned
  accessions (`ENSG00000141510.16`) while every real GMT gives symbols or *unversioned*
  accessions. No set ever met `min_size`, so the result was always empty — and the CLI
  reported `0 tested` rather than an error, so the run looked successful. A four-gene set
  overlapping the background perfectly returned nothing, lost entirely to the `.16`.
  Identifiers now go through `normalize_gene_id` (version suffix dropped, case-folded).
- Gene symbols reach the differential table: `read_gtf_junctions` returns `gene_name`,
  `annotate_junctions` carries it, and `differential_splicing` keeps it. Symbol-keyed
  gene sets — MSigDB `*.symbols.gmt`, GO, KEGG, i.e. most of them — had nothing to match
  against before. `enrich_differential` uses whichever of `gene_id`/`gene_name` overlaps
  the sets more.
- `run --gene-sets` now warns, with an example identifier from each side, when no gene set
  shares an identifier with the tested genes.

### Fixed (packaging and CLI)
- **`run` reported "0 significant junctions" and exited 0 when sample names did not
  match.** Names were derived with `stem.replace(".SJ.out", "")`, but STAR writes
  `{prefix}SJ.out.tab`, so the near-universal real filename `SampleA_SJ.out.tab` yielded
  the sample name `SampleA_SJ.out` while `groups.tsv` said `SampleA`. Nothing checked the
  overlap: conditions mapped to NaN, no sample joined either group, the test returned an
  empty frame, and the user got a clean run, an empty table and the conclusion "no
  differential splicing". Names now tolerate any separator before `SJ.out`, a total
  mismatch is an error naming both sets, and a partial one warns.
- The single source of truth for the version is `splicescope.__version__`; `pyproject.toml`
  reads it dynamically.
- `savefig` now closes the figure it wrote. pyplot keeps every figure alive until closed,
  and `run` writes up to eight per invocation, so a run leaked all of them and matplotlib
  warned past twenty. `close=False` keeps the old behaviour.

## [0.8.1] — 2026-09-03

### Fixed
- Streamlit dashboard deployment: the app now imports from `src`, declares its runtime
  dependencies explicitly, guards against load errors, and drops the removed
  `use_container_width` argument.

## [0.8.0] — 2026-08-30

### Added
- **Splice-site shifts are now interpreted**, not just cassette exons. A novel donor
  or acceptor paired with an annotated site on the other side is resolved against the
  annotation into the sequence it adds to, or removes from, the neighbouring exon
  (`junction_change`, `junction_hosts`, `predict_junction_consequence`,
  `annotate_junction_consequences`). These are the majority of cryptic events reported
  by junction-level callers: on a published set of 3,548 such junctions, interpretable
  coverage went from **0% to 96%** (98% of novel donors, 99% of novel acceptors).
- **Premature stops downstream of the event are found.** Previously an insert that
  shifted the frame without carrying a stop was reported as `frameshift` with no PTC,
  and a truncation got no PTC search at all. Both now assemble the retained downstream
  exons and continue the in-frame scan (`downstream_sequence`), so the 50-nucleotide
  rule can be applied where the stop actually is.
- New class `exon_truncation` for a splice-site shift that shortens an exon without
  producing a stop.

### Changed
- Against experimentally measured NMD sensitivity the class ordering is now monotonic
  and the effect is stronger. On cassette exons (i3Neuron, quantitative label) the
  ambiguous `frameshift` bucket disappears and classes rank
  `ptc_nmd` 12.7 > `ptc_escape` 8.6 > in-frame 6.3 > UTR 4.2 > non-coding −1.6 ΔPSI.
  On splice-site shifts (SMG1i, binary label) `ptc_nmd` is 48.5% NMD-sensitive against
  23.5% for events where no stop is possible — **odds ratio 2.31 → 3.07,
  two-sided Fisher p 1e-04 → 4.5e-12**.

## [0.7.0] — 2026-08-30

### Changed
- **The differential test now models read counts** (`betabinom`), replacing the
  Mann–Whitney U test on per-sample Ψ as the default wherever counts are available.
  A rank test on `n` replicates per group cannot return a p-value below `2/C(2n,n)`
  — 0.1 at 3-vs-3 — so nothing could ever survive genome-wide BH correction. On a real
  human TDP-43 knockdown (GSE245332, 3 vs 3) the old test returned **0 significant
  junctions with min(q) = 1.0**, while rMATS called 3 724 events on the same BAMs, and
  46 junctions had |ΔΨ| ≥ 0.5. The signal was there; the statistic could not reach it.
- `differential_splicing` gained `test` (`"auto"`, `"betabinom"`, `"ranksum"`) plus
  `inc_col`/`total_col`. `"auto"` picks the count-based test when counts are present and
  falls back to the rank test otherwise, so existing calls keep working.
- `compute_psi` now also returns `donor_total` and `acceptor_total`, the denominators Ψ
  was formed from, so junction-level tests can use counts too.

### Added
- `betabinom`: vectorised beta-binomial fitting (`fit_mu`, `loglik`), a df-corrected
  moment estimator for the shared dispersion (`estimate_precision`) and the
  likelihood-ratio test (`lrt`). Fits are bisections on a monotone score, so hundreds of
  thousands of events cost seconds.
- `min_achievable_rank_pvalue` so the limitation above can be shown rather than asserted.
- 11 tests: Ψ recovery, boundary events, dispersion bias, null calibration
  (nominal 0.05 → 0.045–0.053) and power past the rank-test floor.

### Fixed
- `detect_mxe_events` compared every junction against every other one on a chromosome to
  find candidate exons, which is quadratic and dominated the whole pipeline (81% of
  runtime on a real dataset). A candidate exon is at most `max_exon` long, so only
  downstream junctions starting inside that window can pair with a given upstream
  junction; sorting by start and binary-searching the window is exact and much cheaper.
  **89.5 s → 0.5 s** on two chromosomes, with identical output. Together with the
  consequence fix below this takes a full six-sample run from ~45 minutes to **~1.5
  minutes**; only the row order of `events.tsv` changes.
- `annotate_consequences` scanned every transcript for every event, which is quadratic:
  296,034 cassette exons took ~34 minutes. Transcripts are now indexed by gene
  (`index_by_gene`) and only the event's own gene is searched — the same run takes
  **19 seconds**, and events can no longer be attributed to an overlapping neighbour.

### Documentation
- `validation/` records an end-to-end run on real data (GSE245332) plus a replication
  experiment against a second, independent dataset (GSE122069): `bam2sj.awk` and
  `extract_sj.sh` derive STAR-format junctions from any sorted BAM, and `add_motifs.py`
  fills in intron motifs from an indexed genome (85–90% canonical, and the motif agrees
  with the aligner's `XS:A` tag for 99.94% of junctions).
- The classifier's positioning in the README and METHODS is corrected. Scored against
  replication in the independent dataset it reaches an average precision of 0.006 for
  confidently replicating junctions where the beta-binomial statistic reaches 0.062, and
  the combination is worse than the statistic alone. Its features describe whether a
  junction is plausible, not whether it changes between conditions; it is a noise filter,
  not the way to find cryptic events.

### Notes
- Dispersion is estimated about each group's own mean. Pooling groups that genuinely
  differ charges the difference to dispersion — at ΔΨ = 0.35 the pooled estimate collapses
  from `s ≈ 200` to `s ≈ 5` and power goes to zero.

## [0.6.0] — 2026-08-30

### Added
- **Protein-consequence prediction** (`consequence`): for each cassette exon, inherit the
  reading frame from the coding sequence upstream of the host intron, translate the
  insert, locate the first in-frame premature termination codon (PTC) and apply the
  50-nucleotide rule to decide whether it triggers nonsense-mediated decay. Events are
  classified as `ptc_nmd`, `ptc_escape`, `frameshift`, `in_frame_insertion`,
  `utr_insertion`, `non_coding_host` or `no_host_transcript`. Detecting a cryptic exon
  says nothing about whether it matters; this is what closes that gap.
- `GenomeFasta`: random access to an indexed genome through its `.fai` alone, so
  sequence resolution adds no compiled dependency.
- `load_transcripts`: exon/CDS models from a plain or gzipped GTF, optionally restricted
  to a set of gene names or Ensembl IDs (versioned or not).
- CLI `consequence --events --gtf --genome --out` for candidates produced elsewhere
  (e.g. an existing rMATS run), and `run --genome` to fold the same prediction into the
  end-to-end pipeline as `consequence.tsv`.
- 11 tests covering frame inheritance, PTC detection, the NMD distance rule, strand
  handling and the FASTA index reader.

### Notes
- Use the **full** GENCODE annotation rather than `basic`: on 300 real cryptic exons the
  reduced set left 43% of events without a host transcript, against 20% with the full
  annotation.

## [0.5.0] — 2026-08-30

### Added
- **Pathway over-representation analysis** (`enrich`): a hypergeometric/BH-FDR test for
  gene sets enriched among differentially spliced genes — the same statistic as
  clusterProfiler, on any `{term: [genes]}` mapping. No gene sets are bundled (real
  enrichment needs real annotations).
- GMT gene-set reader (`io.read_gmt`, MSigDB/GO/KEGG-compatible) and an enrichment plot.
- CLI `run --gene-sets sets.gmt` writes `enrichment.tsv` and a figure.
- An **executed tutorial notebook** (`examples/tutorial.ipynb`, rebuilt via
  `examples/_build_tutorial.py`) with a Colab badge.
- Tests for the ORA math, GMT parsing and end-to-end enrichment.

## [0.4.0] — 2026-08-29

### Added
- **Mutually-exclusive-exon (MXE) events** — two non-overlapping exons between shared
  flanking exons, with PSI = exon-A share. splicescope now covers 4 of the 5 canonical
  rMATS event classes (SE, MXE, A5SS, A3SS; RI remains out of scope — needs coverage).
- MXE junctions are excluded from downstream A5SS/A3SS detection; the simulator can inject
  MXE events (`mxe_fraction`) and suppresses the skipping isoform for those introns.
- CLI `simulate --mxe`; the showcase gains MXE in the event-type and event-volcano panels.
- Tests for MXE detection and event-level differential inclusion.

## [0.3.0] — 2026-08-29

### Added
- **Alternative 5′/3′ splice-site events** (A5SS, A3SS) alongside cassette exons, with a
  unified `detect_events` / `event_psi` API and rMATS-style PSI per type.
- Cassette-exon junctions are excluded from alt-splice-site detection so shared signal is
  not double-reported.
- `differential_splicing` gained a `key` argument for event-level tests (`key=["event_id"]`).
- Simulator can inject A5SS/A3SS events (`alt_ss_fraction`); CLI `run` writes `events.tsv`,
  `event_differential.tsv` and an event volcano; the showcase figure is now a 2×3 panel.
- Tests for A5SS/A3SS detection and event-level differential inclusion.

### Notes
- Intron retention (RI) is intentionally out of scope: it cannot be quantified from
  junctions alone (it needs intronic read coverage). Documented in METHODS.

## [0.2.0] — 2026-08-29

### Added
- **Event-level splicing** (`events`): reconstruct cassette (skipped) exons directly from
  junctions and compute rMATS-style percent-spliced-in (Ψ). Strand-independent detection.
- `differential_splicing` now works on any Ψ column (e.g. `psi_cassette`), so cassette
  events get ΔΨ + Mann–Whitney + FDR for free.
- CLI `run` writes `cassette_events.tsv`, `cassette_differential.tsv` and a
  `cassette_volcano.png`.
- Tests for cassette detection, the PSI formula, and end-to-end differential inclusion.

## [0.1.0] — 2026-08-27

### Added
- Junction I/O for STAR `SJ.out.tab` and known-intron extraction from GTF (`io`).
- Junction taxonomy: annotated / novel_donor / novel_acceptor / novel_combination /
  cryptic (`annotate`).
- Model-free splice-site usage Ψ (`quantify`).
- Differential splicing: ΔΨ, Mann–Whitney U, Benjamini–Hochberg FDR (`diff`).
- Feature engineering for cryptic-event calling (`cryptic`).
- Cross-validated `CrypticClassifier` with permutation importance and a model card (`ml`).
- Biologically faithful synthetic data generator with ground-truth labels and optional
  label noise (`simulate`).
- Publication-quality, headless-safe plots (`plotting`).
- `splicescope` CLI (`simulate`, `run`).
- Nextflow (DSL2) pipeline and a Streamlit dashboard.
- Test suite (unit + property + end-to-end) and GitHub Actions CI across Python 3.10–3.12.
