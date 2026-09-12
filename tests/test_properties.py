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

import numpy as np
import pandas as pd
from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st

from splicescope.annotate import CLASSES, annotate_junctions
from splicescope.consequence import Transcript, find_ptc, reverse_complement
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
