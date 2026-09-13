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
    """The strongest check available: the simulator reports where it put the exons.
    Both the original one-pair-per-anchor bug and two geometric attempts at bounding
    the search lost 15-22% of these while still emitting a plausible-looking table.

    This asserted `>= 55 of 60` until the simulator gained a truth table, because the
    injected set had to be reconstructed from the coordinates and that arithmetic gives
    every geometry the simulator *could* have used. Recall is exact now.
    """
    from splicescope.annotate import annotate_junctions
    from splicescope.events import detect_mxe_events

    for seed in (4, 6, 11):
        ds = simulate_dataset(n_genes=60, n_per_group=6, mxe_fraction=1.0,
                              alt_ss_fraction=0.8, cryptic_fraction=0.5, seed=seed)
        events = detect_mxe_events(annotate_junctions(ds.observed, ds.known))
        injected = {
            (r.exonA_start, r.exonA_end, r.exonB_start, r.exonB_end)
            for r in ds.truth[ds.truth.event_type == "MXE"].itertuples(index=False)
        }
        found = {
            (r.exonA_start, r.exonA_end, r.exonB_start, r.exonB_end)
            for r in events.itertuples(index=False)
        }
        assert injected, "the simulator was asked for MXE events and reported none"
        assert injected <= found, (
            f"seed {seed} missed {len(injected - found)} of {len(injected)} injected pairs"
        )


def test_every_injected_event_of_every_type_is_recovered():
    """Recall for all four types at once, against what the simulator says it injected
    rather than against what the coordinates allow. Noise junctions can form further
    valid geometries, so the detector may report more than was injected; what it may
    not do is miss one.
    """
    from splicescope.annotate import annotate_junctions

    for seed in (4, 6, 11, 17):
        ds = simulate_dataset(n_genes=60, n_per_group=6, mxe_fraction=0.4,
                              alt_ss_fraction=0.8, cryptic_fraction=0.6, seed=seed)
        events = detect_events(annotate_junctions(ds.observed, ds.known))
        truth = ds.truth

        # a cryptic exon is a cassette: its host intron is the skipping junction
        se = events[events.event_type == "SE"]
        found = set(zip(se.skip_start, se.skip_end, se.exon_start, se.exon_end, strict=True))
        want = {
            (r.intron_start, r.intron_end, r.exonA_start, r.exonA_end)
            for r in truth[truth.event_type == "cryptic_exon"].itertuples(index=False)
        }
        assert want <= found, f"seed {seed}: {len(want - found)} of {len(want)} cryptic exons"

        for kind in ("A5SS", "A3SS"):
            want = set(truth.loc[truth.event_type == kind, "site_pos"])
            found = set(events.loc[events.event_type == kind, "site_pos"])
            assert want and want <= found, (
                f"seed {seed}: {len(want - found)} of {len(want)} {kind} sites"
            )

        mxe = events[events.event_type == "MXE"]
        found = set(
            zip(mxe.exonA_start, mxe.exonA_end, mxe.exonB_start, mxe.exonB_end, strict=True)
        )
        want = {
            (r.exonA_start, r.exonA_end, r.exonB_start, r.exonB_end)
            for r in truth[truth.event_type == "MXE"].itertuples(index=False)
        }
        assert want <= found, f"seed {seed}: {len(want - found)} of {len(want)} MXE pairs"


def test_one_intron_carries_at_most_one_injected_event():
    """What makes the recall above meaningful. Two events in one intron do not just
    crowd each other, they change what the reads mean: an MXE intron has no skipping
    junction, so a cryptic exon placed in it is not a detectable cassette, and an
    alternative donor there is also a leg of an MXE pair and is reported as that. Both
    readings are right and both make the truth table claim events that are not there —
    which is why MXE recall used to be quoted as a floor.
    """
    for seed in (4, 6, 11, 17):
        ds = simulate_dataset(n_genes=60, n_per_group=6, mxe_fraction=0.4,
                              alt_ss_fraction=0.8, cryptic_fraction=0.6, seed=seed)
        introns = list(zip(ds.truth.intron_start, ds.truth.intron_end, strict=True))
        assert len(introns) == len(set(introns)), f"seed {seed} reused an intron"


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


def _annotation_scale_dataset(n_genes=400, n_exons=10, n_samples=3):
    """Many genes sharing splice sites, the way a real annotation does.

    Every junction is present in every sample, so nothing is lost to sampling and the
    expected event set is exact. Each gene gets ten 150 nt exons 2 300 nt apart, and
    depending on its index an alternative donor, a skipped exon, a cryptic acceptor, or
    several of those. The gene index also sets the strand, so each event class appears
    on both — which matters, because donor and acceptor swap roles between them.
    """
    chroms = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]
    rows, expected = [], {"SE": set(), "A5SS": set(), "A3SS": set()}
    for g in range(n_genes):
        chrom, strand = chroms[g % len(chroms)], "+" if g % 2 == 0 else "-"
        base = 1_000 + (g // len(chroms)) * 60_000
        starts = [base + k * 2_300 for k in range(n_exons)]
        ends = [s + 149 for s in starts]
        gid, gname = f"ENSG{g:011d}", f"GENE{g}"
        for k in range(n_exons - 1):
            rows.append((chrom, ends[k] + 1, starts[k + 1] - 1, strand, gid, gname, True))
        if g % 5 == 0:  # an alternative donor 27 nt in, sharing intron 2's other end
            rows.append((chrom, ends[2] - 26, starts[3] - 1, strand, gid, gname, True))
            expected["A5SS" if strand == "+" else "A3SS"].add((chrom, strand, starts[3] - 1))
        if g % 3 == 0:  # exon 4 skipped
            rows.append((chrom, ends[3] + 1, starts[5] - 1, strand, gid, gname, True))
            expected["SE"].add((chrom, strand, ends[3] + 1, starts[5] - 1))
        if g % 4 == 0:  # a cryptic acceptor inside intron 6, on the annotated donor
            rows.append((chrom, ends[6] + 1, ends[6] + 900, strand, gid, gname, False))
            expected["A3SS" if strand == "+" else "A5SS"].add((chrom, strand, ends[6] + 1))

    cols = ["chrom", "start", "end", "strand", "gene_id", "gene_name", "in_annotation"]
    pool = pd.DataFrame(rows, columns=cols)
    known = pool[pool["in_annotation"]].drop(columns="in_annotation").reset_index(drop=True)
    rng = np.random.default_rng(0)
    observed = pd.concat(
        [
            pool[["chrom", "start", "end", "strand"]].assign(
                sample=f"S{s}", count=rng.integers(20, 80, len(pool))
            )
            for s in range(n_samples)
        ],
        ignore_index=True,
    )
    return observed, known, expected


def test_every_event_in_an_annotation_sized_dataset_is_recovered_and_nothing_else_is():
    """Counts can hide a detector that drops one event and invents another, so this
    compares the event *identities* against the ones injected — across 400 genes on
    both strands, sharing donors and acceptors the way real genes do.

    A skipping junction is genuinely an alternative 3' site relative to its own
    inclusion junction, so it is excluded from alt-site detection once its cassette is
    called. That exclusion is what keeps these three sets disjoint, and it is the part
    an exact comparison pins down.
    """
    observed, known, expected = _annotation_scale_dataset()
    events = detect_events(annotate_junctions(observed, known))

    se = events[events["event_type"] == "SE"]
    got_se = set(zip(se["chrom"], se["strand"], se["skip_start"], se["skip_end"], strict=True))
    assert got_se == expected["SE"], (
        f"{len(expected['SE'] - got_se)} skipped exons missed, "
        f"{len(got_se - expected['SE'])} invented"
    )

    for kind in ("A5SS", "A3SS"):
        sub = events[events["event_type"] == kind]
        got = set(zip(sub["chrom"], sub["strand"], sub["site_pos"], strict=True))
        assert got == expected[kind], (
            f"{kind}: {len(expected[kind] - got)} missed, {len(got - expected[kind])} invented"
        )

    # the cassette must name the exon that was actually skipped, not merely some exon
    exon = se.set_index(["chrom", "strand", "skip_start"])[["exon_start", "exon_end"]]
    chrom, strand, skip_start, _ = sorted(expected["SE"])[0]
    row = exon.loc[(chrom, strand, skip_start)]
    assert row["exon_end"] - row["exon_start"] == 149


def test_junctions_on_different_chromosomes_are_never_combined_into_one_event():
    """Coordinates repeat across chromosomes — every chromosome has a position 8 050 —
    so a grouping key that forgets the chromosome pools junctions from genes that have
    nothing to do with each other and assembles events out of the pieces.

    Each half below is deliberately incomplete: chr1 has a skipping junction and its
    upstream inclusion junction, chr2 only the downstream one, and neither pair of
    anchors holds both ends of a mutually exclusive exon. No event exists on either
    chromosome, and the only way to report one is to merge them.
    """
    obs = pd.DataFrame(
        [
            # a cassette missing its downstream inclusion junction, which sits on chr2
            ("chr1", 100, 300, "+", 40),
            ("chr1", 100, 199, "+", 40),
            ("chr2", 250, 300, "+", 40),
            # two MXE donors on chr1, the matching acceptors on chr2
            ("chr1", 2000, 2100, "+", 40),
            ("chr1", 2000, 2200, "+", 40),
            ("chr2", 2150, 9000, "+", 40),
            ("chr2", 2250, 9000, "+", 40),
        ],
        columns=["chrom", "start", "end", "strand", "count"],
    )
    obs["sample"] = "s1"

    cassettes = detect_cassette_events(obs)
    assert cassettes.empty, f"invented {len(cassettes)} cassette(s) across chromosomes"
    mxe = detect_mxe_events(obs)
    assert mxe.empty, f"invented {len(mxe)} MXE event(s) across chromosomes"


def test_the_event_level_test_uses_the_counts_and_not_the_ranks():
    """`differential_splicing(value="psi")` finds its count columns by name, and the
    names it looks for are the ones `event_psi` happens to use. Nothing checked that
    they still match, and the failure is silent: with the counts missing the test falls
    back to Mann-Whitney, which at three replicates per group cannot clear correction at
    all. Measured on the simulator: 32 of 44 events called with the counts, 0 without.
    """
    from splicescope.diff import differential_splicing
    from splicescope.simulate import simulate_dataset

    ds = simulate_dataset(n_genes=24, n_per_group=3, cryptic_fraction=0.8, seed=11)
    annotated = annotate_junctions(ds.observed, ds.known)
    psi = event_psi(annotated, detect_events(annotated), min_reads=5)
    assert {"inc_reads", "total_reads"} <= set(psi.columns), (
        "the column names differential_splicing resolves for value='psi'"
    )

    counted = differential_splicing(psi, ds.groups, value="psi", key=["event_id"])
    assert "lrt_statistic" in counted.columns, "the beta-binomial, not the rank test"
    assert (counted["qvalue"] <= 0.05).sum() > 10, "and it calls a good number of them"

    # the same data with the counts hidden, which is what a rename would do
    ranked = differential_splicing(
        psi.drop(columns=["inc_reads", "total_reads"]),
        ds.groups, value="psi", key=["event_id"],
    )
    assert "lrt_statistic" not in ranked.columns
    assert (ranked["qvalue"] <= 0.05).sum() == 0, (
        "3 vs 3 on ranks alone clears nothing — which is why the fallback must not be "
        "reachable by accident"
    )


def _annotated_with_a_cassette_and_an_alt_site():
    """One cassette exon and one alternative acceptor, in two samples.

    Deliberately no mutually exclusive exon: the point of the tests below is that a
    dataset which contains some event types but not others still gets every column.
    """
    rows = []
    for sample in ("s1", "s2"):
        rows += [
            # cassette: 1000-2000 skipping, 1000-1500 + 1600-2000 including
            _junction("chr1", 1000, 2000, sample=sample, count=20),
            _junction("chr1", 1000, 1500, sample=sample, count=60),
            _junction("chr1", 1601, 2000, sample=sample, count=60),
            # two acceptors reached from one donor at 5000
            _junction("chr1", 5000, 6000, sample=sample, count=40),
            _junction("chr1", 5000, 6030, sample=sample, count=15),
        ]
    return pd.DataFrame(rows)


def test_the_event_schema_comes_from_the_request_not_from_the_data():
    """The third and fourth instances of one defect, found by sweeping for it.

    `detect_events` returned whatever `pd.concat` made of the types it happened to find.
    A dataset with cassettes but no mutually exclusive exons came back without the MXE
    columns, and one with no events at all came back with five columns where a full
    result has thirty — so `events["exon_start"]` raised a KeyError or did not depending
    on the data, and the showcase script's own first move (select the SE rows, hand
    `exon_start` to the consequence layer) was one eventless dataset away from raising.

    What a caller can reason about is which types they asked for, so that is what the
    schema is now pinned to.
    """
    from splicescope.events import detect_events, event_columns

    frame = _annotated_with_a_cassette_and_an_alt_site()
    everything = detect_events(frame)
    nothing = detect_events(frame.iloc[0:0])

    assert not everything.empty and nothing.empty
    assert list(everything.columns) == list(nothing.columns) == event_columns()
    assert "exonA_start" in nothing.columns, "MXE columns, though no MXE was found"
    # the showcase's first move, on the result that used to raise
    assert nothing[nothing["event_type"] == "SE"][["exon_start", "exon_end"]].empty

    # and asking for less gets less, found or not
    for types in [("SE",), ("MXE",), ("A5SS", "A3SS")]:
        found = detect_events(frame, types=types)
        empty = detect_events(frame.iloc[0:0], types=types)
        assert list(found.columns) == list(empty.columns) == event_columns(types)
        assert set(found.columns) < set(event_columns()), f"{types} asks for fewer"

    # A5SS and A3SS describe a site the same way: asking for both adds one set, not two
    assert event_columns(("A5SS", "A3SS")) == event_columns(("A5SS",))


def test_event_psi_has_its_columns_even_with_nothing_to_measure():
    """`event_psi` already returned the full schema when there were no events, and a
    column-less frame when there were events but no sample carrying them — the same
    function, disagreeing with itself about what an empty result looks like."""
    from splicescope.events import PSI_COLUMNS, detect_events, event_psi

    frame = _annotated_with_a_cassette_and_an_alt_site()
    evs = detect_events(frame)

    measured = event_psi(frame, evs, min_reads=1)
    no_events = event_psi(frame, evs.iloc[0:0], min_reads=1)
    no_samples = event_psi(frame.iloc[0:0], evs, min_reads=1)

    assert not measured.empty and no_events.empty and no_samples.empty
    for frame_ in (measured, no_events, no_samples):
        assert list(frame_.columns) == PSI_COLUMNS
    # what a caller does next with a result that found nothing
    assert no_samples.groupby("event_id")["psi"].mean().empty
    assert no_samples.to_csv(index=False).splitlines()[0] == (
        measured.to_csv(index=False).splitlines()[0]
    )


def test_event_coordinates_are_written_as_positions_not_as_measurements():
    """`events.tsv` said an exon began at `1840.0`.

    A table holding more than one event type has a missing value wherever a column
    belongs to a different type, and float64 is the only NumPy dtype that can hold one —
    so every coordinate in a mixed table picked up a decimal point. A genomic position
    with a `.0` reads as a rounded measurement, and a downstream tool parsing the column
    as an integer fails on it. A nullable integer column holds the same values and the
    same gaps, and writes neither.
    """
    from splicescope.events import detect_events

    frame = _annotated_with_a_cassette_and_an_alt_site()
    evs = detect_events(frame)
    assert set(evs["event_type"]) == {"SE", "A3SS"}, "a mixed table, which is the hard case"

    header, *rows = evs.to_csv(sep="\t", index=False).splitlines()
    fields = header.split("\t")
    for row in rows:
        for name, value in zip(fields, row.split("\t"), strict=True):
            if name in ("event_id", "event_type", "chrom", "strand", "gene_id", "site_kind"):
                continue
            assert "." not in value, f"{name}={value} is a position, not a measurement"
            assert value == "" or int(value) >= 0

    # and the values survive the dtype: still usable as numbers, still missing where absent
    se = evs[evs["event_type"] == "SE"].iloc[0]
    assert int(se.exon_end) - int(se.exon_start) > 0
    assert pd.isna(se.site_pos), "an SE has no alternative splice site"
