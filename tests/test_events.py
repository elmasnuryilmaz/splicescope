import numpy as np
import pandas as pd

from splicescope.annotate import annotate_junctions
from splicescope.diff import differential_splicing, significant
from splicescope.events import (
    cassette_psi,
    detect_alt_ss_events,
    detect_cassette_events,
    detect_events,
    detect_mxe_events,
    event_psi,
)
from splicescope.simulate import simulate_dataset


def test_detect_simple_cassette():
    # three junctions forming one cassette event on + strand
    obs = pd.DataFrame(
        [
            ("chr1", 100, 300, "+", 40, "s1"),  # skip: donor 100, acceptor 300
            ("chr1", 100, 149, "+", 30, "s1"),  # inc1: shares donor 100
            ("chr1", 201, 300, "+", 30, "s1"),  # inc2: shares acceptor 300
        ],
        columns=["chrom", "start", "end", "strand", "count", "sample"],
    )
    events = detect_cassette_events(obs)
    assert len(events) == 1
    ev = events.iloc[0]
    assert (ev["exon_start"], ev["exon_end"]) == (150, 200)
    assert (ev["skip_start"], ev["skip_end"]) == (100, 300)


def test_cassette_psi_formula():
    obs = pd.DataFrame(
        [
            ("chr1", 100, 300, "+", 20, "s1"),  # skip = 20
            ("chr1", 100, 149, "+", 60, "s1"),  # inc1 = 60
            ("chr1", 201, 300, "+", 60, "s1"),  # inc2 = 60
        ],
        columns=["chrom", "start", "end", "strand", "count", "sample"],
    )
    events = detect_cassette_events(obs)
    psi = cassette_psi(obs, events, min_reads=1)
    # inclusion = (60+60)/2 = 60; PSI = 60 / (60+20) = 0.75
    assert np.isclose(psi.loc[0, "psi_cassette"], 0.75)


def test_cassette_events_on_simulation():
    ds = simulate_dataset(n_genes=12, n_per_group=6, cryptic_fraction=1.0, seed=2)
    annotated = annotate_junctions(ds.observed, ds.known)
    events = detect_cassette_events(annotated)
    assert not events.empty

    psi = cassette_psi(annotated, events, min_reads=5)
    vals = psi["psi_cassette"].dropna()
    assert ((vals >= 0) & (vals <= 1)).all()

    diff = differential_splicing(psi, ds.groups, value="psi_cassette")
    # cryptic exons are up-regulated in B, so some events are differentially included
    assert not significant(diff, q=0.1, min_delta=0.05).empty


def test_detect_a5ss_and_a3ss():
    # A5SS: two donors share acceptor 300; A3SS: two acceptors share donor 100
    obs = pd.DataFrame(
        [
            ("chr1", 100, 300, "+", 40, "s1"),  # canonical
            ("chr1", 140, 300, "+", 20, "s1"),  # alt donor  -> A5SS at acceptor 300
            ("chr1", 100, 260, "+", 20, "s1"),  # alt acceptor -> A3SS at donor 100
        ],
        columns=["chrom", "start", "end", "strand", "count", "sample"],
    )
    a5 = detect_alt_ss_events(obs, "A5SS")
    a3 = detect_alt_ss_events(obs, "A3SS")
    assert (a5["event_type"] == "A5SS").all() and len(a5) == 1
    assert (a3["event_type"] == "A3SS").all() and len(a3) == 1
    # A5SS inclusion is the proximal donor (140, closer to acceptor 300)
    assert a5.iloc[0]["incl_start"] == 140


def test_detect_mxe_simple():
    # two non-overlapping exons (A: 150-190, B: 250-290) between donor 100 & acceptor 400
    obs = pd.DataFrame(
        [
            ("chr1", 100, 149, "+", 40, "s1"),  # up -> exon A
            ("chr1", 191, 400, "+", 40, "s1"),  # exon A -> down
            ("chr1", 100, 249, "+", 30, "s1"),  # up -> exon B
            ("chr1", 291, 400, "+", 30, "s1"),  # exon B -> down
        ],
        columns=["chrom", "start", "end", "strand", "count", "sample"],
    )
    mxe = detect_mxe_events(obs)
    assert len(mxe) == 1
    ev = mxe.iloc[0]
    assert (ev["exonA_start"], ev["exonA_end"]) == (150, 190)
    assert (ev["exonB_start"], ev["exonB_end"]) == (250, 290)


def test_mxe_is_differential():
    ds = simulate_dataset(n_genes=16, n_per_group=6, mxe_fraction=1.0, seed=6)
    annotated = annotate_junctions(ds.observed, ds.known)
    events = detect_events(annotated)
    assert "MXE" in set(events["event_type"])
    psi = event_psi(annotated, events, min_reads=5)
    diff = differential_splicing(psi, ds.groups, value="psi", key=["event_id"])
    mxe_hits = significant(diff, q=0.1, min_delta=0.1)
    assert (mxe_hits["event_type"] == "MXE").any()


def test_unified_events_and_event_level_differential():
    ds = simulate_dataset(
        n_genes=16, n_per_group=6, cryptic_fraction=0.5, alt_ss_fraction=0.8,
        mxe_fraction=0.5, seed=4,
    )
    annotated = annotate_junctions(ds.observed, ds.known)
    events = detect_events(annotated)
    kinds = set(events["event_type"])
    assert {"SE", "A5SS", "A3SS"} & kinds  # at least the injected types appear
    assert kinds <= {"SE", "MXE", "A5SS", "A3SS"}

    psi = event_psi(annotated, events, min_reads=5)
    vals = psi["psi"].dropna()
    assert ((vals >= 0) & (vals <= 1)).all()

    diff = differential_splicing(psi, ds.groups, value="psi", key=["event_id"])
    assert "event_type" in diff.columns
    assert not significant(diff, q=0.1, min_delta=0.05).empty


def _junction(chrom, start, end, strand="+", sample="s1", count=50, gene="G1"):
    return {
        "chrom": chrom, "start": start, "end": end, "strand": strand,
        "sample": sample, "count": count, "sclass": "cryptic", "gene_id": gene,
    }


def test_mxe_window_keeps_valid_pairs_and_respects_max_exon():
    """The windowed search must find exactly the pairs the length filter allows."""
    import pandas as pd

    from splicescope.events import detect_mxe_events

    # Shared donor at 1000 and shared acceptor at 5000, with two candidate exons
    # (1101-1200 and 2101-2200) reachable through their own junction pairs.
    rows = [
        _junction("chr1", 1000, 1100),   # donor -> exon A
        _junction("chr1", 1201, 5000),   # exon A -> acceptor
        _junction("chr1", 1000, 2100),   # donor -> exon B
        _junction("chr1", 2201, 5000),   # exon B -> acceptor
    ]
    events = detect_mxe_events(pd.DataFrame(rows))
    assert len(events) == 1
    ev = events.iloc[0]
    assert (ev.exonA_start, ev.exonA_end) == (1101, 1200)
    assert (ev.exonB_start, ev.exonB_end) == (2101, 2200)

    # An exon longer than max_exon must fall outside the window and find nothing.
    assert detect_mxe_events(pd.DataFrame(rows), max_exon=50).empty


def test_every_mutually_exclusive_pair_is_emitted_not_just_the_first():
    """Stopping at the first pair lost tandem MXE clusters entirely, and — worse —
    kept the genomically leftmost pair rather than the best-supported one, so a noise
    junction to the left could displace the real event."""
    from splicescope.events import detect_mxe_events

    def j(start, end):
        return dict(
            chrom="chr1", start=start, end=end, strand="+",
            sample="s1", count=50, gene_id="g1",
        )

    # three candidate exons between the same donor (1000) and acceptor (5000)
    df = pd.DataFrame(
        [
            j(1000, 1100), j(1201, 5000),   # exon A 1101-1200
            j(1000, 2100), j(2201, 5000),   # exon B 2101-2200
            j(1000, 3100), j(3201, 5000),   # exon C 3101-3200
        ]
    )
    events = detect_mxe_events(df)
    pairs = {
        (r.exonA_start, r.exonA_end, r.exonB_start, r.exonB_end)
        for r in events.itertuples(index=False)
    }
    assert pairs == {
        (1101, 1200, 2101, 2200),
        (1101, 1200, 3101, 3200),
        (2101, 2200, 3101, 3200),
    }


def test_an_mxe_exon_is_not_lost_from_detect_events_altogether():
    """A dropped MXE exon could not fall through to A5SS/A3SS either: its junctions
    share the donor and acceptor of the emitted pair, so after exclusion each site was
    left with a single alternative and the exon vanished from every event type."""
    from splicescope.events import detect_events

    def j(start, end):
        return dict(
            chrom="chr1", start=start, end=end, strand="+",
            sample="s1", count=50, gene_id="g1",
        )

    df = pd.DataFrame(
        [
            j(1000, 1100), j(1201, 5000),
            j(1000, 2100), j(2201, 5000),
            j(1000, 3100), j(3201, 5000),
        ]
    )
    table = detect_events(df)
    text = table.to_string()
    for boundary in ("1101", "2101", "3101"):
        assert boundary in text, f"exon starting at {boundary} is missing from events"


def test_max_exon_is_reachable_through_detect_events():
    """The 1 kb candidate window was hard-coded for everyone using the public entry point."""
    from splicescope.events import detect_events

    def j(start, end):
        return dict(
            chrom="chr1", start=start, end=end, strand="+",
            sample="s1", count=50, gene_id="g1",
        )

    # both candidate exons are 1,500 bp long — outside the default window
    df = pd.DataFrame(
        [
            j(1000, 1100), j(2601, 9000),   # exon 1101-2600
            j(1000, 3100), j(4601, 9000),   # exon 3101-4600
        ]
    )
    assert detect_events(df, types=("MXE",)).empty
    widened = detect_events(df, types=("MXE",), max_exon=2000)
    assert len(widened) == 1
