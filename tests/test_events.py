import numpy as np
import pandas as pd

from splicescope.annotate import annotate_junctions
from splicescope.diff import differential_splicing, significant
from splicescope.events import (
    cassette_psi,
    detect_alt_ss_events,
    detect_cassette_events,
    detect_events,
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


def test_unified_events_and_event_level_differential():
    ds = simulate_dataset(
        n_genes=16, n_per_group=6, cryptic_fraction=0.5, alt_ss_fraction=0.8, seed=4
    )
    annotated = annotate_junctions(ds.observed, ds.known)
    events = detect_events(annotated)
    kinds = set(events["event_type"])
    assert {"SE", "A5SS", "A3SS"} & kinds  # at least the injected types appear
    assert kinds <= {"SE", "A5SS", "A3SS"}

    psi = event_psi(annotated, events, min_reads=5)
    vals = psi["psi"].dropna()
    assert ((vals >= 0) & (vals <= 1)).all()

    diff = differential_splicing(psi, ds.groups, value="psi", key=["event_id"])
    assert "event_type" in diff.columns
    assert not significant(diff, q=0.1, min_delta=0.05).empty
