"""Event-level splicing: cassette (skipped) exons and their PSI.

Junction *usage* (see :mod:`splicescope.quantify`) is model-free but not an event.
Here we reconstruct the canonical **cassette-exon** event directly from junctions
and compute the standard *percent-spliced-in* (Ψ).

A cassette exon is defined by three junctions that share endpoints:

```
      inclusion 1            inclusion 2
   ┌───────────────┐      ┌───────────────┐
[ up-exon ]      [ cassette exon ]      [ down-exon ]
   └──────────────────────────────────────┘
                skipping junction
```

* **inclusion 1** ``(s3, e1)`` shares the *start* of the skip and lands on the
  cassette exon's 5′ boundary,
* **inclusion 2** ``(s2, e3)`` shares the *end* of the skip and leaves the exon's
  3′ boundary, with ``e1 < s2`` (the exon is ``[e1+1, s2-1]``),
* **skipping** ``(s3, e3)`` splices straight from up-exon to down-exon.

This coordinate pattern is strand-independent (the strand only has to match), so
detection works for both strands. PSI follows the junction-count definition used
by rMATS-style methods:

    inclusion = (count(inc1) + count(inc2)) / 2
    Ψ = inclusion / (inclusion + count(skip))
"""

from __future__ import annotations

from collections import defaultdict

import pandas as pd

from .io import donor_acceptor

EVENT_KEY = "event_id"
EVENT_TYPES = ("SE", "MXE", "A5SS", "A3SS")


def _unique_junctions(annotated: pd.DataFrame) -> pd.DataFrame:
    cols = ["chrom", "start", "end", "strand"]
    extra = [c for c in ("gene_id",) if c in annotated.columns]
    return annotated.drop_duplicates(subset=cols)[cols + extra].reset_index(drop=True)


def detect_cassette_events(annotated: pd.DataFrame) -> pd.DataFrame:
    """Find cassette-exon events from a table of (annotated) junctions.

    Returns one row per event with the three junction coordinate pairs, the
    inferred cassette-exon interval, and a best-effort ``gene_id``.
    """
    uniq = _unique_junctions(annotated)
    by_start: dict[tuple, list[int]] = defaultdict(list)
    by_end: dict[tuple, list[int]] = defaultdict(list)
    gene_of = {}
    for row in uniq.itertuples(index=False):
        by_start[(row.chrom, row.strand, row.start)].append(row.end)
        by_end[(row.chrom, row.strand, row.end)].append(row.start)
        gene_of[(row.chrom, row.start, row.end, row.strand)] = getattr(row, "gene_id", None)

    events = []
    seen = set()
    for row in uniq.itertuples(index=False):
        c, st, s3, e3 = row.chrom, row.strand, row.start, row.end
        inc1_ends = [e1 for e1 in by_start[(c, st, s3)] if e1 < e3]
        inc2_starts = [s2 for s2 in by_end[(c, st, e3)] if s2 > s3]
        if not inc1_ends or not inc2_starts:
            continue
        for e1 in inc1_ends:
            for s2 in inc2_starts:
                if e1 >= s2:  # inclusion junctions must not overlap
                    continue
                exon_start, exon_end = e1 + 1, s2 - 1
                if exon_end < exon_start:
                    continue
                key = (c, st, s3, e3, e1, s2)
                if key in seen:
                    continue
                seen.add(key)
                events.append(
                    {
                        "event_id": f"SE:{c}:{s3}-{e3}:{st}@{exon_start}-{exon_end}",
                        "event_type": "SE",
                        "chrom": c,
                        "strand": st,
                        "gene_id": gene_of.get((c, s3, e3, st)),
                        "exon_start": exon_start,
                        "exon_end": exon_end,
                        "inc1_start": s3,
                        "inc1_end": e1,
                        "inc2_start": s2,
                        "inc2_end": e3,
                        "skip_start": s3,
                        "skip_end": e3,
                    }
                )
    return pd.DataFrame(events)


def cassette_psi(
    annotated: pd.DataFrame, events: pd.DataFrame, min_reads: int = 10
) -> pd.DataFrame:
    """Compute per-sample cassette-exon PSI for each detected event.

    Returns a long table keyed by the *skipping* junction coordinates (so it can
    be fed straight into :func:`splicescope.diff.differential_splicing` with
    ``value="psi_cassette"``), plus ``event_id``, inclusion/skip read counts and
    ``gene_id``.
    """
    if events.empty:
        return pd.DataFrame(
            columns=[
                "chrom", "start", "end", "strand", "sample",
                "psi_cassette", "inc_reads", "skip_reads", "event_id", "gene_id",
            ]
        )

    counts = (
        annotated.groupby(["chrom", "start", "end", "strand", "sample"], observed=True)["count"]
        .sum()
        .to_dict()
    )
    samples = sorted(annotated["sample"].unique())

    rows = []
    for ev in events.itertuples(index=False):
        for sample in samples:
            i1 = counts.get((ev.chrom, ev.inc1_start, ev.inc1_end, ev.strand, sample), 0)
            i2 = counts.get((ev.chrom, ev.inc2_start, ev.inc2_end, ev.strand, sample), 0)
            skip = counts.get((ev.chrom, ev.skip_start, ev.skip_end, ev.strand, sample), 0)
            inclusion = (i1 + i2) / 2.0
            denom = inclusion + skip
            psi = inclusion / denom if denom >= min_reads else float("nan")
            rows.append(
                {
                    "chrom": ev.chrom,
                    "start": ev.skip_start,
                    "end": ev.skip_end,
                    "strand": ev.strand,
                    "sample": sample,
                    "psi_cassette": psi,
                    "inc_reads": inclusion,
                    "skip_reads": float(skip),
                    "event_id": ev.event_id,
                    "gene_id": ev.gene_id,
                }
            )
    return pd.DataFrame(rows)


def _with_sites(annotated: pd.DataFrame) -> pd.DataFrame:
    df = annotated.copy()
    da = [
        donor_acceptor(s, e, st)
        for s, e, st in zip(df["start"], df["end"], df["strand"], strict=False)
    ]
    df["donor"] = [d for d, _ in da]
    df["acceptor"] = [a for _, a in da]
    return df


def detect_mxe_events(
    annotated: pd.DataFrame, max_exon: int = 1000, exclude: set | None = None
) -> pd.DataFrame:
    """Detect mutually-exclusive-exon (MXE) events.

    An MXE is two non-overlapping exons, A (upstream) and B (downstream), that sit
    between the *same* upstream donor ``s`` and downstream acceptor ``e`` — each
    reached by its own pair of inclusion junctions, and normally never included
    together. PSI is exon A's share:

        Ψ = incl(A) / (incl(A) + incl(B))

    where ``incl(X) = (count(X_5' junction) + count(X_3' junction)) / 2``.
    """
    exclude = exclude or set()
    uniq = _unique_junctions(annotated)
    gene_of = {
        (r.chrom, r.start, r.end, r.strand): getattr(r, "gene_id", None)
        for r in uniq.itertuples(index=False)
    }
    by_cs: dict[tuple, list] = defaultdict(list)
    for r in uniq.itertuples(index=False):
        if (r.chrom, r.start, r.end, r.strand) in exclude:
            continue
        by_cs[(r.chrom, r.strand)].append(r)

    events, seen = [], set()
    for (chrom, strand), js in by_cs.items():
        # every (exon between donor s and acceptor e), grouped by that (s, e) context
        paths: dict[tuple, list] = defaultdict(list)
        for j5 in js:
            for j3 in js:
                if j3.start <= j5.end:
                    continue
                exon_len = j3.start - j5.end - 1
                if not 1 <= exon_len <= max_exon:
                    continue
                paths[(j5.start, j3.end)].append((j5.end + 1, j3.start - 1, j5, j3))

        for (s, e), plist in paths.items():
            if len(plist) < 2:
                continue
            plist.sort(key=lambda p: (p[0], p[1]))
            emitted = False
            for i in range(len(plist)):
                for k in range(i + 1, len(plist)):
                    a, b = plist[i], plist[k]
                    if a[1] < b[0]:  # exon A strictly upstream of exon B
                        key = (chrom, strand, s, e, a[0], a[1], b[0], b[1])
                        if key in seen:
                            continue
                        seen.add(key)
                        events.append(
                            {
                                "event_id": (
                                    f"MXE:{chrom}:{s}-{e}:{strand}"
                                    f"@{a[0]}-{a[1]}|{b[0]}-{b[1]}"
                                ),
                                "event_type": "MXE",
                                "chrom": chrom,
                                "strand": strand,
                                "gene_id": gene_of.get((chrom, a[2].start, a[2].end, strand)),
                                "exonA_start": a[0],
                                "exonA_end": a[1],
                                "exonB_start": b[0],
                                "exonB_end": b[1],
                                "a_j1_start": a[2].start,
                                "a_j1_end": a[2].end,
                                "a_j2_start": a[3].start,
                                "a_j2_end": a[3].end,
                                "b_j1_start": b[2].start,
                                "b_j1_end": b[2].end,
                                "b_j2_start": b[3].start,
                                "b_j2_end": b[3].end,
                            }
                        )
                        emitted = True
                        break
                if emitted:
                    break
    return pd.DataFrame(events)


def detect_alt_ss_events(
    annotated: pd.DataFrame, kind: str, exclude: set | None = None
) -> pd.DataFrame:
    """Detect alternative 5′ (``kind="A5SS"``) or 3′ (``kind="A3SS"``) splice-site events.

    An A5SS event is a splice acceptor served by two or more donors; A3SS is the
    mirror (one donor, several acceptors). The *inclusion* isoform is the one whose
    variable site sits closest to the shared site (the longer exon), following the
    rMATS convention. One event is emitted per shared site with ≥2 alternatives.

    ``exclude`` is a set of ``(chrom, start, end, strand)`` junctions to ignore —
    used to keep cassette-exon inclusion junctions from being re-reported as
    alt-splice-site events.
    """
    if kind not in ("A5SS", "A3SS"):
        raise ValueError("kind must be 'A5SS' or 'A3SS'")
    shared, variable = ("acceptor", "donor") if kind == "A5SS" else ("donor", "acceptor")
    exclude = exclude or set()

    uniq = _with_sites(_unique_junctions(annotated))
    groups: dict[tuple, list] = defaultdict(list)
    for row in uniq.itertuples(index=False):
        if (row.chrom, row.start, row.end, row.strand) in exclude:
            continue
        groups[(row.chrom, getattr(row, shared), row.strand)].append(row)

    events = []
    for (chrom, site, strand), members in groups.items():
        if len({getattr(m, variable) for m in members}) < 2:
            continue
        incl = min(members, key=lambda m: abs(getattr(m, variable) - site))
        events.append(
            {
                "event_id": f"{kind}:{chrom}:{site}:{strand}",
                "event_type": kind,
                "chrom": chrom,
                "strand": strand,
                "gene_id": getattr(incl, "gene_id", None),
                "site_pos": site,
                "site_kind": shared,
                "incl_start": incl.start,
                "incl_end": incl.end,
                "n_alternatives": len(members),
            }
        )
    return pd.DataFrame(events)


def detect_events(annotated: pd.DataFrame, types: tuple[str, ...] = EVENT_TYPES) -> pd.DataFrame:
    """Detect all requested event types and return them in one typed table.

    Cassette (SE) events are detected first; the junctions they use are then
    excluded from A5SS/A3SS detection so the same signal is not double-reported.
    """
    parts = []
    used: set = set()
    if "SE" in types:
        se = detect_cassette_events(annotated)
        parts.append(se)
        for e in se.itertuples(index=False):
            used.add((e.chrom, e.inc1_start, e.inc1_end, e.strand))
            used.add((e.chrom, e.inc2_start, e.inc2_end, e.strand))
            used.add((e.chrom, e.skip_start, e.skip_end, e.strand))
    if "MXE" in types:
        mxe = detect_mxe_events(annotated, exclude=used)
        parts.append(mxe)
        for e in mxe.itertuples(index=False):
            for js, je in (
                (e.a_j1_start, e.a_j1_end),
                (e.a_j2_start, e.a_j2_end),
                (e.b_j1_start, e.b_j1_end),
                (e.b_j2_start, e.b_j2_end),
            ):
                used.add((e.chrom, js, je, e.strand))
    if "A5SS" in types:
        parts.append(detect_alt_ss_events(annotated, "A5SS", exclude=used))
    if "A3SS" in types:
        parts.append(detect_alt_ss_events(annotated, "A3SS", exclude=used))
    parts = [p for p in parts if not p.empty]
    if not parts:
        return pd.DataFrame(columns=["event_id", "event_type", "chrom", "strand", "gene_id"])
    return pd.concat(parts, ignore_index=True)


def event_psi(annotated: pd.DataFrame, events: pd.DataFrame, min_reads: int = 10) -> pd.DataFrame:
    """Per-sample PSI for every event, dispatched by ``event_type``.

    Returns a long table keyed by ``event_id`` (feed to
    :func:`splicescope.diff.differential_splicing` with ``key=["event_id"]``).
    """
    long_cols = ["event_id", "event_type", "gene_id", "sample", "psi", "inc_reads", "total_reads"]
    if events.empty:
        return pd.DataFrame(columns=long_cols)

    df = _with_sites(annotated)
    counts = (
        df.groupby(["chrom", "start", "end", "strand", "sample"], observed=True)["count"]
        .sum()
        .to_dict()
    )
    acc_tot = df.groupby(["chrom", "acceptor", "strand", "sample"], observed=True)["count"].sum()
    don_tot = df.groupby(["chrom", "donor", "strand", "sample"], observed=True)["count"].sum()
    acc_tot, don_tot = acc_tot.to_dict(), don_tot.to_dict()
    samples = sorted(df["sample"].unique())

    rows = []
    for ev in events.itertuples(index=False):
        for sample in samples:
            if ev.event_type == "SE":
                i1 = counts.get(
                    (ev.chrom, int(ev.inc1_start), int(ev.inc1_end), ev.strand, sample), 0
                )
                i2 = counts.get(
                    (ev.chrom, int(ev.inc2_start), int(ev.inc2_end), ev.strand, sample), 0
                )
                skip = counts.get(
                    (ev.chrom, int(ev.skip_start), int(ev.skip_end), ev.strand, sample), 0
                )
                inclusion = (i1 + i2) / 2.0
                total = inclusion + skip
            elif ev.event_type == "MXE":
                a1 = counts.get(
                    (ev.chrom, int(ev.a_j1_start), int(ev.a_j1_end), ev.strand, sample), 0
                )
                a2 = counts.get(
                    (ev.chrom, int(ev.a_j2_start), int(ev.a_j2_end), ev.strand, sample), 0
                )
                b1 = counts.get(
                    (ev.chrom, int(ev.b_j1_start), int(ev.b_j1_end), ev.strand, sample), 0
                )
                b2 = counts.get(
                    (ev.chrom, int(ev.b_j2_start), int(ev.b_j2_end), ev.strand, sample), 0
                )
                inclusion = (a1 + a2) / 2.0  # exon A
                total = inclusion + (b1 + b2) / 2.0  # + exon B
            else:  # A5SS / A3SS
                inclusion = counts.get(
                    (ev.chrom, int(ev.incl_start), int(ev.incl_end), ev.strand, sample), 0
                )
                site_tot = acc_tot if ev.event_type == "A5SS" else don_tot
                total = site_tot.get((ev.chrom, int(ev.site_pos), ev.strand, sample), 0)
            psi = inclusion / total if total >= min_reads else float("nan")
            rows.append(
                {
                    "event_id": ev.event_id,
                    "event_type": ev.event_type,
                    "gene_id": ev.gene_id,
                    "sample": sample,
                    "psi": psi,
                    "inc_reads": float(inclusion),
                    "total_reads": float(total),
                }
            )
    return pd.DataFrame(rows)
