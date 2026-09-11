import numpy as np
import pandas as pd

from splicescope.annotate import annotate_junctions
from splicescope.diff import differential_splicing, significant
from splicescope.events import (
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
    psi = event_psi(obs, events, min_reads=1)
    # inclusion = (60+60)/2 = 60; PSI = 60 / (60+20) = 0.75
    assert np.isclose(psi.loc[0, "psi"], 0.75)


def test_cassette_events_on_simulation():
    ds = simulate_dataset(n_genes=12, n_per_group=6, cryptic_fraction=1.0, seed=2)
    annotated = annotate_junctions(ds.observed, ds.known)
    events = detect_cassette_events(annotated)
    assert not events.empty

    psi = event_psi(annotated, events, min_reads=5)
    vals = psi["psi"].dropna()
    assert ((vals >= 0) & (vals <= 1)).all()

    diff = differential_splicing(psi, ds.groups, value="psi", key=["event_id"])
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


def test_an_alt_splice_site_event_survives_sharing_an_anchor_with_a_cassette():
    """Cassette junctions are excluded from alt-splice-site detection so the same signal
    is not reported twice — but dropping them outright took the site's other alternatives
    with them, deleting genuine A5SS/A3SS events."""
    from splicescope.events import detect_events

    def j(start, end):
        return dict(
            chrom="chr1", start=start, end=end, strand="+",
            sample="s1", count=50, gene_id="g1",
        )

    df = pd.DataFrame(
        [
            j(1000, 1100),   # cassette inclusion 1
            j(1201, 2000),   # cassette inclusion 2
            j(1000, 2000),   # cassette skip
            j(1000, 1700),   # a third acceptor at the same donor: a real A3SS alternative
        ]
    )
    types = set(detect_events(df)["event_type"])
    assert "SE" in types
    assert "A3SS" in types, "the extra acceptor at donor 1000 is a genuine event"


def test_a_cassette_alone_is_not_also_reported_as_an_alt_splice_site_event():
    """The other half of the rule: a site explained entirely by a cassette exon must
    not be told twice."""
    from splicescope.events import detect_events

    def j(start, end):
        return dict(
            chrom="chr1", start=start, end=end, strand="+",
            sample="s1", count=50, gene_id="g1",
        )

    df = pd.DataFrame([j(1000, 1100), j(1201, 2000), j(1000, 2000)])
    assert set(detect_events(df)["event_type"]) == {"SE"}


def test_cassette_events_sharing_a_skip_junction_stay_separate():
    """Several cassette exons can sit between the same pair of flanking exons. Keying
    their PSI by the shared skipping junction merged them into one test, where opposite
    changes cancel; event_psi keys by event_id, so each is tested on its own."""

    def j(start, end, sample, count):
        return dict(chrom="chr1", start=start, end=end, strand="+", sample=sample, count=count)

    rows = []
    for sample in ("A1", "A2", "B1", "B2"):
        up = sample.startswith("B")
        rows += [
            j(1000, 2000, sample, 100),                  # the shared skipping junction
            j(1000, 1100, sample, 90 if up else 10),     # exon X, up in B
            j(1201, 2000, sample, 90 if up else 10),
            j(1000, 1500, sample, 10 if up else 90),     # exon Y, down in B
            j(1601, 2000, sample, 10 if up else 90),
        ]
    observed = pd.DataFrame(rows)
    events = detect_cassette_events(observed)
    assert len(events) > 1, "this locus must yield several cassette events"

    psi = event_psi(observed, events, min_reads=1)
    assert psi["event_id"].nunique() == len(events)

    groups = {"A1": "A", "A2": "A", "B1": "B", "B2": "B"}
    diff = differential_splicing(psi, groups, value="psi", key=["event_id"])
    assert len(diff) == len(events)
    # the exons move in opposite directions; merged, they would cancel to nothing
    assert diff["delta_psi"].max() > 0.2
    assert diff["delta_psi"].min() < -0.2


def _mxe_cluster(n_exons, spacing=20, count=50):
    """n candidate exons between one donor and one acceptor, all equally supported."""
    rows = []
    for i in range(n_exons):
        boundary = 1000 + i * spacing
        rows.append(dict(chrom="chr1", start=1000, end=boundary, strand="+",
                         sample="s1", count=count, gene_id="g1"))
        rows.append(dict(chrom="chr1", start=boundary + 10, end=9000, strand="+",
                         sample="s1", count=count, gene_id="g1"))
    return pd.DataFrame(rows).drop_duplicates(["chrom", "start", "end", "strand"])


def test_every_injected_mxe_event_is_recovered():
    """The strongest check available: the simulator knows where it put the exons.
    Both the original one-pair-per-anchor bug and two geometric attempts at bounding
    the search lost 15-22% of these while still emitting a plausible-looking table."""
    from splicescope.annotate import annotate_junctions
    from splicescope.events import detect_mxe_events

    for seed in (4, 6, 11):
        ds = simulate_dataset(n_genes=60, n_per_group=6, mxe_fraction=1.0,
                              alt_ss_fraction=0.8, cryptic_fraction=0.5, seed=seed)
        events = detect_mxe_events(annotate_junctions(ds.observed, ds.known))
        injected = {
            (i + 50, i + 90, i + 150, i + 190)
            for i in (r.start for r in ds.known.itertuples(index=False))
        }
        found = {
            (r.exonA_start, r.exonA_end, r.exonB_start, r.exonB_end)
            for r in events.itertuples(index=False)
        }
        missing = injected - found
        # not every known intron gets an MXE injected, so compare against what did
        assert len(missing) == len(injected) - len(injected & found)
        assert injected & found, f"seed {seed} recovered no injected pair"
        recovered = len(injected & found)
        assert recovered >= 55, f"seed {seed} recovered only {recovered} of the injected pairs"


def test_a_crowded_anchor_is_bounded_and_says_so():
    """Every junction combination at an anchor is a candidate, so n exons arrive with
    about n**2 of them and pairing all of those grew as n**4 — 80 exons produced 1.68
    million rows in 5.9 s. The bound must hold, and must not be silent."""
    import warnings

    from splicescope.events import detect_mxe_events

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        events = detect_mxe_events(_mxe_cluster(80), max_candidates=50)
    assert len(events) <= 50 * 49 // 2
    assert any("candidate exons" in str(w.message) for w in caught)


def test_the_best_supported_candidates_survive_the_bound():
    """Geometry cannot tell a real exon from a span or from a noise junction sharing the
    anchor's donor — every structural rule tried here deleted real events. Read support
    can: an exon is only as good as its weaker flanking junction."""
    from splicescope.events import detect_mxe_events

    rows = []
    real = [(1100, 1200), (1400, 1500), (1700, 1800)]
    for exon_start, exon_end in real:                      # deep support
        rows.append(dict(chrom="chr1", start=1000, end=exon_start - 1, strand="+",
                         sample="s1", count=500, gene_id="g1"))
        rows.append(dict(chrom="chr1", start=exon_end + 1, end=5000, strand="+",
                         sample="s1", count=500, gene_id="g1"))
    for i in range(60):                                    # shallow noise at the same anchor
        rows.append(dict(chrom="chr1", start=1000, end=1210 + i * 3, strand="+",
                         sample="s1", count=1, gene_id="g1"))
        rows.append(dict(chrom="chr1", start=1230 + i * 3, end=5000, strand="+",
                         sample="s1", count=1, gene_id="g1"))
    observed = pd.DataFrame(rows).drop_duplicates(["chrom", "start", "end", "strand"])

    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        events = detect_mxe_events(observed, max_candidates=20)
    pairs = {
        (r.exonA_start, r.exonA_end, r.exonB_start, r.exonB_end)
        for r in events.itertuples(index=False)
    }
    for a in range(len(real)):
        for b in range(a + 1, len(real)):
            expected = (*real[a], *real[b])
            assert expected in pairs, f"lost the well-supported pair {expected}"


def test_max_candidates_reaches_through_detect_events():
    import warnings

    from splicescope.events import detect_events

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        wide = detect_events(_mxe_cluster(60), types=("MXE",), max_candidates=200)
        narrow = detect_events(_mxe_cluster(60), types=("MXE",), max_candidates=10)
    assert len(wide) > len(narrow)


def test_a_dense_tandem_array_keeps_every_real_exon():
    """Support ranking alone is not enough. In a tandem array every candidate is flanked
    by equally deep junctions, so the order under the cap is arbitrary and the cut takes
    real exons with it — 11 densely packed exons kept only 6. A span runs from one exon's
    start to a later exon's end, so it is always longer than the real exon sharing its
    start: shortest-first breaks the tie the right way."""
    from splicescope.events import detect_mxe_events

    def array(n_exons, exon=60, gap=40, count=500):
        rows, pos = [], 1000
        for _ in range(n_exons):
            rows.append(dict(chrom="chr1", start=1000, end=pos, strand="+",
                             sample="s1", count=count, gene_id="g1"))
            rows.append(dict(chrom="chr1", start=pos + exon + 1, end=90000, strand="+",
                             sample="s1", count=count, gene_id="g1"))
            pos += exon + gap
        return pd.DataFrame(rows).drop_duplicates(["chrom", "start", "end", "strand"])

    import warnings

    for n in (8, 11, 15, 25):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            events = detect_mxe_events(array(n))
        expected = {(1000 + i * 100 + 1, 1000 + i * 100 + 60) for i in range(n)}
        seen = {(r.exonA_start, r.exonA_end) for r in events.itertuples(index=False)}
        seen |= {(r.exonB_start, r.exonB_end) for r in events.itertuples(index=False)}
        assert expected <= seen, f"{n}-exon array lost {len(expected - seen)} real exons"
