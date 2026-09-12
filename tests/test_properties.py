"""Invariants, on inputs nobody would think to write.

The README said the suite was "unit + property + end-to-end". Nineteen tests had
property-shaped names, but they all asserted invariants on hand-picked inputs; there was
no generator anywhere. That is a different technique, and this session has been a
demonstration that different techniques see different things — the mutation survey found
rules the suite was not checking, and coverage found five features it was not running at
all. Neither can produce an input a person would not have written.

So: real generators, and the invariants that have to hold whatever they produce.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st

from splicescope.annotate import CLASSES, _check_annotation, annotate_junctions
from splicescope.consequence import (
    CONSEQUENCE_CLASSES,
    Transcript,
    describe,
    find_ptc,
    reverse_complement,
)
from splicescope.diff import benjamini_hochberg
from splicescope.enrich import normalize_gene_id
from splicescope.events import detect_events
from splicescope.io import donor_acceptor
from splicescope.quantify import compute_psi

SLOW = settings(
    max_examples=40,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.data_too_large],
)

CHROMS = st.sampled_from(["chr1", "chr2", "chrX"])
STRANDS = st.sampled_from(["+", "-"])
DNA = st.text(alphabet="ACGT", min_size=0, max_size=90)


@st.composite
def junction_table(draw, min_rows=1, max_rows=25, samples=("S1", "S2")):
    """Junctions with plausible coordinates, repeated across samples."""
    n = draw(st.integers(min_value=min_rows, max_value=max_rows))
    rows = []
    for _ in range(n):
        chrom = draw(CHROMS)
        strand = draw(STRANDS)
        start = draw(st.integers(min_value=101, max_value=9_000))
        length = draw(st.integers(min_value=30, max_value=3_000))
        for sample in samples:
            rows.append(
                {
                    "chrom": chrom, "start": start, "end": start + length, "strand": strand,
                    "sample": sample,
                    "count": draw(st.integers(min_value=0, max_value=500)),
                    "motif": draw(st.sampled_from(["GT/AG", "CT/AC", "non-canonical"])),
                }
            )
    frame = pd.DataFrame(rows)
    return frame.drop_duplicates(subset=["chrom", "start", "end", "strand", "sample"])


@st.composite
def known_table(draw, max_rows=12):
    n = draw(st.integers(min_value=1, max_value=max_rows))
    rows = []
    for i in range(n):
        start = draw(st.integers(min_value=101, max_value=9_000))
        rows.append(
            {
                "chrom": draw(CHROMS),
                "start": start,
                "end": start + draw(st.integers(min_value=30, max_value=3_000)),
                "strand": draw(STRANDS),
                "gene_id": f"G{i}",
                "gene_name": f"GENE{i}",
            }
        )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- pure functions


@given(DNA)
def test_reverse_complement_is_an_involution(seq):
    assert reverse_complement(reverse_complement(seq)) == seq


@given(
    st.integers(min_value=1, max_value=10**6),
    st.integers(min_value=1, max_value=10**6),
    STRANDS,
)
def test_donor_and_acceptor_are_the_intron_ends_in_transcription_order(start, end, strand):
    donor, acceptor = donor_acceptor(start, end, strand)
    assert {donor, acceptor} == {start, end}
    assert (donor, acceptor) == ((start, end) if strand == "+" else (end, start))


@given(DNA, st.integers(min_value=0, max_value=2))
def test_find_ptc_points_at_a_real_in_frame_stop_or_at_nothing(sequence, frame):
    offset = find_ptc(sequence, frame)
    if offset is None:
        return
    assert sequence[offset : offset + 3] in {"TAA", "TAG", "TGA"}
    assert (offset - (3 - frame) % 3) % 3 == 0, "the hit must sit on a codon boundary"
    # and it must be the *first* such codon
    begin = (3 - frame) % 3
    earlier = range(begin, offset, 3)
    assert all(sequence[i : i + 3] not in {"TAA", "TAG", "TGA"} for i in earlier)


@given(st.lists(st.floats(min_value=0.0, max_value=1.0, allow_nan=False), min_size=1, max_size=60))
def test_bh_qvalues_never_undercut_their_pvalues_and_stay_monotone(pvalues):
    p = np.array(pvalues)
    q = benjamini_hochberg(p)
    assert q.shape == p.shape
    assert ((q >= p - 1e-12) | np.isnan(q)).all(), "a q-value cannot beat its p-value"
    assert ((q <= 1.0 + 1e-12) | np.isnan(q)).all()
    order = np.argsort(p, kind="stable")
    ranked = q[order]
    assert (np.diff(ranked) >= -1e-12).all(), "q must not fall as p rises"


@given(st.lists(st.floats(min_value=0.0, max_value=1.0, allow_nan=False), min_size=2, max_size=40))
def test_bh_commutes_with_permutation(pvalues):
    p = np.array(pvalues)
    order = np.random.default_rng(len(pvalues)).permutation(len(p))
    assert np.allclose(benjamini_hochberg(p)[order], benjamini_hochberg(p[order]), atol=1e-12)


@given(st.text(min_size=1, max_size=24))
def test_normalizing_a_gene_id_twice_changes_nothing(gene):
    once = normalize_gene_id(gene)
    assert normalize_gene_id(once) == once


# ---------------------------------------------------------------- the pipeline


@given(junction_table(), known_table())
@SLOW
def test_annotation_labels_every_row_exactly_once_and_keeps_the_table_intact(observed, known):
    out = annotate_junctions(observed, known)
    assert len(out) == len(observed)
    assert out.index.equals(observed.index)
    assert set(observed.columns) <= set(out.columns)
    assert out["sclass"].isin(CLASSES).all()
    assert out["sclass"].notna().all(), "classification is a total function"
    assert (out["is_novel"] == (out["sclass"] != "annotated")).all()


@given(junction_table(), known_table())
@SLOW
def test_annotation_does_not_depend_on_the_order_of_the_rows(observed, known):
    """A junction's class comes from its coordinates and the annotation, nothing else."""
    shuffled = observed.sample(frac=1.0, random_state=0)
    plain = annotate_junctions(observed, known).set_index(
        ["chrom", "start", "end", "strand", "sample"]
    )
    other = annotate_junctions(shuffled, known).set_index(
        ["chrom", "start", "end", "strand", "sample"]
    )
    assert plain["sclass"].sort_index().equals(other["sclass"].sort_index())


@given(junction_table(), known_table())
@SLOW
def test_psi_is_a_fraction_or_unmeasured_and_sums_to_one_per_site(observed, known):
    psi = compute_psi(annotate_junctions(observed, known), min_reads=1)
    for value, site in (("psi_donor", "donor"), ("psi_acceptor", "acceptor")):
        measured = psi[psi[value].notna()]
        assert measured[value].between(-1e-12, 1 + 1e-12).all(), f"{value} left [0, 1]"
        totals = measured.groupby(["chrom", site, "strand", "sample"], observed=True)[value].sum()
        assert np.allclose(totals, 1.0, atol=1e-9), (
            f"{value} must exhaust its site: {totals[~np.isclose(totals, 1.0)].head()}"
        )


@st.composite
def clustered_junctions(draw, samples=("S1", "S2")):
    """Junctions that share anchors, so events actually form.

    Purely random coordinates almost never do: an event needs two or three junctions
    meeting at the same splice site. Each locus here gets one donor anchor with several
    acceptors and one acceptor anchor with several donors, which is what an alternative
    splice site, a cassette exon and a mutually exclusive pair are all made of.
    """
    rows = []
    for locus in range(draw(st.integers(min_value=1, max_value=3))):
        chrom, strand = draw(CHROMS), draw(STRANDS)
        base = 1_000 + locus * 40_000
        anchor_start = base + draw(st.integers(min_value=0, max_value=200))
        anchor_end = anchor_start + draw(st.integers(min_value=600, max_value=3_000))
        ends = draw(
            st.lists(
                st.integers(min_value=200, max_value=590), min_size=2, max_size=4, unique=True
            )
        )
        starts = draw(
            st.lists(
                st.integers(min_value=200, max_value=590), min_size=2, max_size=4, unique=True
            )
        )
        pairs = [(anchor_start, anchor_start + off) for off in ends]
        pairs += [(anchor_end - off, anchor_end) for off in starts]
        pairs.append((anchor_start, anchor_end))
        for start, end in pairs:
            if end <= start:
                continue
            for sample in samples:
                rows.append(
                    {
                        "chrom": chrom, "start": start, "end": end, "strand": strand,
                        "sample": sample,
                        "count": draw(st.integers(min_value=1, max_value=300)),
                        "motif": "GT/AG",
                    }
                )
    frame = pd.DataFrame(rows)
    return frame.drop_duplicates(subset=["chrom", "start", "end", "strand", "sample"])


@given(clustered_junctions(), known_table())
@SLOW
def test_every_event_is_built_from_junctions_that_were_observed(observed, known):
    """A detector may find nothing, and may find more than a person expected, but it may
    not invent a coordinate."""
    annotated = annotate_junctions(observed, known)
    events = detect_events(annotated)
    assume(not events.empty)

    seen = set(
        zip(
            annotated["chrom"], annotated["start"], annotated["end"], annotated["strand"],
            strict=True,
        )
    )
    starts = {s for _, s, _, _ in seen}
    ends = {e for _, _, e, _ in seen}
    for row in events.itertuples(index=False):
        assert row.chrom in {c for c, _, _, _ in seen}
        assert row.strand in {"+", "-"}
        for field in ("skip_start", "inc1_start", "inc2_start", "a_j1_start", "b_j1_start"):
            value = getattr(row, field, None)
            if value is not None and value == value:
                assert int(value) in starts, f"{field}={value} was never observed"
        for field in ("skip_end", "inc1_end", "inc2_end", "a_j1_end", "b_j1_end"):
            value = getattr(row, field, None)
            if value is not None and value == value:
                assert int(value) in ends, f"{field}={value} was never observed"


@given(
    st.lists(
        st.tuples(st.integers(min_value=1, max_value=500), st.integers(min_value=1, max_value=300)),
        min_size=1,
        max_size=8,
    ),
    STRANDS,
)
def test_coding_length_before_a_position_is_monotone_along_transcription(blocks, strand):
    """It counts bases already transcribed, so walking forward can only add to it."""
    exons, cursor = [], 100
    for gap, length in blocks:
        cursor += gap
        exons.append((cursor, cursor + length))
        cursor += length
    tx = Transcript(
        transcript_id="T", gene_id="G", gene_name="G", chrom="chr1", strand=strand,
        exons=list(exons), cds=list(exons),
    )
    total = sum(e - s + 1 for s, e in exons)
    positions = [s for s, _ in exons] + [e for _, e in exons]
    walked = sorted(positions, reverse=strand == "-")
    lengths = [tx.cds_length_before(p) for p in walked]
    assert all(b >= a for a, b in zip(lengths, lengths[1:], strict=False)), lengths
    assert 0 <= min(lengths) and max(lengths) <= total


@given(junction_table(), known_table())
@SLOW
def test_filling_the_omitted_zeros_does_not_move_an_observed_junction(observed, known):
    """The whole justification for `add_unobserved_zeros` is that a junction an aligner
    omitted is a measured zero. A zero adds nothing to its site's total, so the Ψ of
    every junction that *was* observed has to come out the same either way. If filling
    moved them, it would be changing measurements rather than completing them."""
    annotated = annotate_junctions(observed, known)
    key = ["chrom", "start", "end", "strand", "sample"]
    plain = compute_psi(annotated, min_reads=1, fill_unobserved=False).set_index(key)
    filled = compute_psi(annotated, min_reads=1, fill_unobserved=True).set_index(key)

    assert set(plain.index) <= set(filled.index), "filling only adds rows"
    for value in ("psi_donor", "psi_acceptor", "donor_total", "acceptor_total"):
        before = plain[value]
        after = filled.loc[before.index, value]
        both = before.notna() & after.notna()
        assert np.allclose(before[both], after[both], atol=1e-12), value
        assert before.isna().equals(after.isna()), f"{value} changed which rows are measured"


@given(junction_table(min_rows=1, max_rows=10, samples=("A1", "A2", "B1", "B2")), known_table())
@SLOW
def test_swapping_the_group_labels_flips_delta_psi_and_leaves_the_pvalue_alone(observed, known):
    """A symmetry of both tests: which group is called A is a naming choice. ΔΨ is
    B minus A so it must change sign, and a two-sided p-value must not move at all."""
    from splicescope.diff import differential_splicing

    psi = compute_psi(annotate_junctions(observed, known), min_reads=1)
    key = ["chrom", "start", "end", "strand"]
    first = {"A1": "aa", "A2": "aa", "B1": "bb", "B2": "bb"}
    second = {"A1": "bb", "A2": "bb", "B1": "aa", "B2": "aa"}

    # generated junctions are often too sparse to test, and differential_splicing says
    # so; `assume` below is what handles it, so the warning is not this test's subject
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        forward = differential_splicing(psi, first, min_samples=2).set_index(key)
        backward = differential_splicing(psi, second, min_samples=2).set_index(key)
    assume(not forward.empty)

    assert set(forward.index) == set(backward.index)
    reversed_delta = backward.loc[forward.index, "delta_psi"]
    assert np.allclose(forward["delta_psi"], -reversed_delta, atol=1e-9, equal_nan=True)
    assert np.allclose(
        forward["pvalue"], backward.loc[forward.index, "pvalue"], atol=1e-9, equal_nan=True
    )


@given(clustered_junctions(), known_table(), CHROMS)
@SLOW
def test_a_locus_on_another_chromosome_changes_nothing(observed, known, spare):
    """Events are local. Adding an unrelated cluster elsewhere may add events of its own,
    but it must not alter, remove or renumber any event that was already there."""
    assume(spare not in set(observed["chrom"]))
    before = detect_events(annotate_junctions(observed, known))
    assume(not before.empty)

    elsewhere = observed.assign(chrom=spare, start=observed["start"] + 500_000,
                                end=observed["end"] + 500_000)
    after = detect_events(
        annotate_junctions(pd.concat([observed, elsewhere], ignore_index=True), known)
    )

    kept = after[after["chrom"] != spare].reset_index(drop=True)
    original = before.reset_index(drop=True)
    assert list(kept["event_id"]) == list(original["event_id"])
    for column in ("event_type", "skip_start", "skip_end", "site_pos", "n_alternatives"):
        if column in original:
            assert kept[column].equals(original[column]), column


@given(
    st.sampled_from(sorted(CONSEQUENCE_CLASSES)),
    st.integers(min_value=-4000, max_value=4000),
    st.booleans(),
    st.one_of(st.none(), st.integers(min_value=0, max_value=6000)),
    st.one_of(st.none(), st.integers(min_value=-3000, max_value=6000)),
    st.booleans(),
)
def test_describing_a_consequence_never_raises_and_never_says_nothing(
    kind, length, frameshift, ptc_offset, distance, nmd
):
    """`describe` is read by people, from a row that may carry anything the pipeline can
    produce: a negative insert length for a truncation, a missing PTC offset, a missing
    distance where the stop is in the final exon. It has to be a sentence whatever the
    combination — and never the catch-all, which exists only for a class added later."""
    row = {
        "consequence_class": kind, "insert_length": length, "frameshift": frameshift,
        "ptc_offset": ptc_offset, "distance_to_last_junction": distance,
        "nmd_predicted": nmd,
    }
    sentence = describe(row)
    assert sentence.endswith(".") and sentence[0].isupper()
    assert "None" not in sentence and "nan" not in sentence
    assert "does not name" not in sentence, f"{kind} fell through to the catch-all"
    assert str(abs(length)) in sentence or kind in {"no_host_transcript"}


@given(
    st.lists(CHROMS, min_size=1, max_size=4, unique=True),
    st.lists(CHROMS, min_size=1, max_size=4, unique=True),
)
@SLOW
def test_the_annotation_check_only_objects_to_a_naming_mismatch(seen, annotated):
    """It must not fire for a junction on a contig the annotation does not cover, which
    is ordinary, and it must fire whenever the two name the same chromosomes differently,
    which is a mistake. Those are the only two cases."""
    observed = pd.DataFrame(
        [dict(chrom=c, start=200, end=399, strand="+", sample="s1", count=10) for c in seen]
    )
    known = pd.DataFrame(
        [dict(chrom=c, start=200, end=399, strand="+", gene_id="G") for c in annotated]
    )
    try:
        _check_annotation(observed, known)
        objected = False
    except ValueError:
        objected = True

    def bare(names):
        return {n[3:] if n.lower().startswith("chr") else n for n in names}

    shares_a_name = bool(set(seen) & set(annotated))
    only_the_prefix_differs = not shares_a_name and bool(bare(seen) & bare(annotated))
    assert objected == only_the_prefix_differs, (
        f"seen={seen} annotated={annotated}: objected={objected}"
    )


@given(st.integers(min_value=1, max_value=5_000), st.integers(min_value=1, max_value=5_000))
def test_the_genome_check_accepts_any_contig_long_enough(reach, slack):
    """A genome larger than the analysis needs is nobody's mistake, and the check must
    say nothing about it however much larger."""
    import tempfile

    from splicescope.consequence import GenomeFasta, Transcript, check_genome

    length = reach + slack
    with tempfile.TemporaryDirectory() as tmp:
        fasta = Path(tmp) / "g.fa"
        fasta.write_text(">chr1\n" + "A" * length + "\n")
        fasta.with_suffix(".fa.fai").write_text(f"chr1\t{length}\t6\t{length}\t{length + 1}\n")
        tx = {
            "T": Transcript(
                transcript_id="T", gene_id="G", gene_name="G", chrom="chr1", strand="+",
                exons=[(1, 10)], cds=[(1, 10)],
            )
        }
        events = pd.DataFrame([dict(chrom="chr1", start=1, end=reach, strand="+")])
        with GenomeFasta(fasta) as genome:
            check_genome(events, tx, genome)  # must not raise


@st.composite
def gene_universe(draw, n_genes=None):
    """A background, a subset of hits, per-gene opportunity, and some gene sets."""
    n = n_genes or draw(st.integers(min_value=6, max_value=60))
    genes = [f"G{i:04d}" for i in range(n)]
    hits = draw(st.lists(st.sampled_from(genes), min_size=1, max_size=n, unique=True))
    weights = {g: float(draw(st.integers(min_value=1, max_value=40))) for g in genes}
    sets = {
        f"S{j}": draw(st.lists(st.sampled_from(genes), min_size=2, max_size=n, unique=True))
        for j in range(draw(st.integers(min_value=1, max_value=4)))
    }
    return genes, hits, weights, sets


@given(gene_universe())
@SLOW
def test_over_representation_reports_probabilities_and_the_same_sets_either_way(universe):
    """Whatever the overlap, a p-value is a probability and a fold enrichment is
    positive. And weighting changes how surprising a set is, never which sets are
    testable — so the two runs must agree on the rows even where they disagree on the
    numbers."""
    from splicescope.enrich import RESULT_COLUMNS, over_representation

    genes, hits, weights, sets = universe
    plain = over_representation(hits, genes, sets)
    weighted = over_representation(hits, genes, sets, weights=weights)

    assert set(plain["term"]) == set(weighted["term"]), "weighting is not a filter"
    for frame in (plain, weighted):
        # the schema does not depend on whether anything came out
        assert list(frame.columns) == RESULT_COLUMNS
        if frame.empty:
            continue
        assert frame["pvalue"].between(0.0, 1.0).all()
        assert frame["qvalue"].between(0.0, 1.0).all()
        assert (frame["fold_enrichment"] > 0).all()
        assert (frame["overlap"] > 0).all(), "a set with no hit is not reported"
        assert (frame["overlap"] <= frame["set_size"]).all()
        assert (frame["bias_odds"] > 0).all()
        # the docstring promises "sorted by q-value", and a reader takes the top row
        assert frame["qvalue"].is_monotonic_increasing


@given(gene_universe())
@SLOW
def test_every_background_gene_gets_a_propensity_strictly_inside_zero_and_one(universe):
    """The odds `p/(1-p)` that the weighted null averages are infinite at either end, so
    the smoothing has to keep every estimate off both. This is what stopped one gene
    rewriting a whole enrichment table."""
    from splicescope.enrich import selection_propensity

    genes, hits, weights, _ = universe
    propensity = selection_propensity(genes, hits, weights)
    assert set(propensity) == {g.upper() for g in genes}
    assert all(0.0 < p < 1.0 for p in propensity.values()), min(propensity.values())
    assert len(set(propensity.values())) <= len(genes)


@given(
    st.integers(min_value=4, max_value=30),
    st.floats(min_value=0.0, max_value=1.0),
    st.floats(min_value=0.0, max_value=1.0),
    st.integers(min_value=0, max_value=10),
)
@SLOW
def test_the_simulator_only_ever_claims_events_it_could_place(n_genes, cryptic, alt, seed):
    """The truth table is what recall is measured against, so a row in it that the data
    does not support would quietly inflate every recall number in the repository."""
    from splicescope.simulate import simulate_dataset

    ds = simulate_dataset(
        n_genes=n_genes, n_per_group=3, cryptic_fraction=cryptic, alt_ss_fraction=alt,
        mxe_fraction=0.3, seed=seed,
    )
    introns = {(r.start, r.end) for r in ds.known.itertuples(index=False)}
    used = list(zip(ds.truth["intron_start"], ds.truth["intron_end"], strict=True))
    assert len(used) == len(set(used)), "one intron carries at most one event"
    for start, end in used:
        assert (start, end) in introns, "the truth names an intron not in the annotation"

    for row in ds.truth.itertuples(index=False):
        for field in ("exonA_start", "exonA_end", "exonB_start", "exonB_end", "alt_pos"):
            value = getattr(row, field)
            if value is not pd.NA and value == value:
                assert row.intron_start <= value <= row.intron_end, (
                    f"{field}={value} lies outside its host intron"
                )
