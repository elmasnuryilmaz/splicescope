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

EVENT_KEY = "event_id"


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
                        "event_id": f"{c}:{s3}-{e3}:{st}@{exon_start}-{exon_end}",
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
