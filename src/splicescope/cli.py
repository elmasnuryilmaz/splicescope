"""Command-line interface for splicescope.

Examples
--------
    splicescope simulate --outdir demo_data          # write a synthetic dataset
    splicescope run --sj-dir demo_data/sj \\
        --gtf demo_data/annotation.gtf \\
        --groups demo_data/groups.tsv --outdir results
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__


def _cmd_simulate(args: argparse.Namespace) -> int:
    from .simulate import simulate_dataset, write_dataset

    ds = simulate_dataset(
        n_genes=args.genes,
        n_per_group=args.replicates,
        alt_ss_fraction=args.alt_ss,
        mxe_fraction=args.mxe,
        seed=args.seed,
    )
    out = write_dataset(ds, args.outdir, seed=args.seed)
    n_junc = ds.observed.drop_duplicates(["chrom", "start", "end", "strand"]).shape[0]
    print(f"[simulate] wrote {n_junc} unique junctions for {len(ds.groups)} samples -> {out}")
    print(f"[simulate] annotation.gtf (with CDS) and genome.fa (+.fai) written to {out}")
    return 0


def _cmd_run(args: argparse.Namespace) -> int:
    import pandas as pd

    from . import annotate as _annot
    from . import cryptic as _cryptic
    from . import diff as _diff
    from . import io as _io
    from . import plotting as _plot
    from . import quantify as _quant
    from .ml import CrypticClassifier

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    sj_paths = {p.stem.replace(".SJ.out", ""): p for p in sorted(Path(args.sj_dir).glob("*.tab"))}
    if not sj_paths:
        print(f"error: no *.tab files in {args.sj_dir}", file=sys.stderr)
        return 2
    observed = _io.read_many_star_sj(sj_paths)
    known = _io.read_gtf_junctions(args.gtf)
    groups_df = pd.read_csv(args.groups, sep="\t")
    groups = dict(zip(groups_df["sample"], groups_df["condition"], strict=False))

    annotated = _annot.annotate_junctions(observed, known)
    summary = _annot.annotation_summary(annotated)
    summary.to_csv(outdir / "annotation_summary.tsv", sep="\t", index=False)

    psi = _quant.compute_psi(annotated, min_reads=args.min_reads)
    diff = _diff.differential_splicing(psi, groups)
    diff.to_csv(outdir / "differential_splicing.tsv", sep="\t", index=False)
    hits = _diff.significant(diff)
    print(f"[run] {len(hits)} significant junctions (q<=0.05, |ΔΨ|>=0.1)")

    figs = outdir / "figures"
    fig, ax = _plot.plt.subplots(figsize=(5, 3.2))
    _plot.plot_annotation_summary(summary, ax=ax)
    _plot.savefig(fig, figs / "annotation_summary.png")
    if not diff.empty:
        fig, ax = _plot.plt.subplots(figsize=(5, 4))
        _plot.plot_volcano(diff, ax=ax)
        _plot.savefig(fig, figs / "volcano.png")

    # event-level: cassette (SE), A5SS and A3SS events with rMATS-style PSI
    from . import events as _events

    evs = _events.detect_events(annotated, max_exon=args.max_exon)
    if not evs.empty:
        evs.to_csv(outdir / "events.tsv", sep="\t", index=False)
        epsi = _events.event_psi(annotated, evs, min_reads=args.min_reads)
        ediff = _diff.differential_splicing(epsi, groups, value="psi", key=["event_id"])
        ediff.to_csv(outdir / "event_differential.tsv", sep="\t", index=False)
        by_type = evs["event_type"].value_counts().to_dict()
        print(
            f"[run] {len(evs)} events {by_type}; "
            f"{len(_diff.significant(ediff))} differentially spliced (ΔΨ)"
        )
        fig, ax = _plot.plt.subplots(figsize=(5, 3.2))
        _plot.plot_event_summary(evs, ax=ax)
        _plot.savefig(fig, figs / "event_summary.png")
        if not ediff.empty and "event_type" in ediff:
            fig, ax = _plot.plt.subplots(figsize=(5, 4))
            _plot.plot_event_volcano(ediff, ax=ax)
            _plot.savefig(fig, figs / "event_volcano.png")

    # protein-level consequence (needs an indexed genome)
    if args.genome:
        if not evs.empty:
            cons = _consequence_table(evs, args.gtf, args.genome, mode="exon")
            if cons is None:
                return 2
            if not cons.empty:
                _write_consequence(cons, outdir / "consequence.tsv", "cassette exons")
                fig, ax = _plot.plt.subplots(figsize=(5.5, 3.2))
                _plot.plot_consequence_summary(cons, ax=ax)
                _plot.savefig(fig, figs / "consequence.png")

        # Splice-site shifts: a novel donor or acceptor anchored on an annotated site.
        # These are the majority of cryptic events junction-level callers report and
        # they are not cassette exons, so they need the junction-level interpretation.
        # Junctions already explained by a detected event are left out: read on its
        # own, a cassette-exon inclusion junction looks like an exon extension that
        # runs to the end of the intron, which is the wrong reading of it.
        shifts = (
            annotated[
                annotated["sclass"].isin(SHIFT_CLASSES)
                & ~_junction_index(annotated).isin(_event_junctions(evs))
            ]
            .drop_duplicates(subset=["chrom", "start", "end", "strand"])
            .loc[:, ["chrom", "start", "end", "strand", "sclass", "gene_id"]]
            .reset_index(drop=True)
        )
        if not shifts.empty:
            jcons = _consequence_table(shifts, args.gtf, args.genome, mode="junction")
            if jcons is None:
                return 2
            _write_consequence(jcons, outdir / "junction_consequence.tsv", "splice-site shifts")
            fig, ax = _plot.plt.subplots(figsize=(5.5, 3.2))
            _plot.plot_consequence_summary(jcons, ax=ax)
            _plot.savefig(fig, figs / "junction_consequence.png")

    # pathway over-representation (only if the user supplies real gene sets)
    if args.gene_sets:
        from . import enrich as _enrich

        gene_sets = _io.read_gmt(args.gene_sets)
        enr = _enrich.enrich_differential(diff, gene_sets)
        enr.to_csv(outdir / "enrichment.tsv", sep="\t", index=False)
        n_sig = int((enr["qvalue"] <= 0.05).sum()) if not enr.empty else 0
        print(f"[run] enrichment: {len(gene_sets)} sets, {len(enr)} tested, {n_sig} sig (q<=0.05)")
        if not enr.empty:
            fig, ax = _plot.plt.subplots(figsize=(6, 3.6))
            _plot.plot_enrichment(enr, ax=ax)
            _plot.savefig(fig, figs / "enrichment.png")

    # cryptic ML (only if truth labels are available, e.g. simulated data)
    if "is_cryptic_truth" in psi.columns:
        feats = _cryptic.extract_features(psi, known)
        if feats["is_cryptic_truth"].nunique() > 1:
            clf = CrypticClassifier(random_state=args.seed)
            metrics = clf.evaluate(feats)
            clf.fit(feats)
            clf.write_model_card(outdir / "model_card.json")
            clf.score_table(feats).to_csv(outdir / "cryptic_scores.tsv", sep="\t", index=False)
            print(
                f"[run] cryptic classifier ROC-AUC={metrics['roc_auc']:.3f} "
                f"AP={metrics['average_precision']:.3f}"
            )
            fig, ax = _plot.plt.subplots(figsize=(5, 3.4))
            _plot.plot_importance(clf.importances, ax=ax)
            _plot.savefig(fig, figs / "importance.png")

    print(f"[run] results written to {outdir}")
    return 0


#: Junction coordinate columns each event type is built from.
_EVENT_JUNCTION_COLUMNS = {
    "SE": (("inc1_start", "inc1_end"), ("inc2_start", "inc2_end"), ("skip_start", "skip_end")),
    "MXE": (
        ("a_j1_start", "a_j1_end"),
        ("a_j2_start", "a_j2_end"),
        ("b_j1_start", "b_j1_end"),
        ("b_j2_start", "b_j2_end"),
    ),
}


def _junction_index(df):
    """A Series of ``(chrom, start, end, strand)`` tuples aligned to ``df``."""
    import pandas as pd

    return pd.Series(
        list(zip(df["chrom"], df["start"], df["end"], df["strand"], strict=False)),
        index=df.index,
    )


def _event_junctions(events) -> set:
    """Every junction coordinate consumed by a detected SE or MXE event."""
    out: set = set()
    if events is None or events.empty or "event_type" not in events:
        return out
    for etype, pairs in _EVENT_JUNCTION_COLUMNS.items():
        sub = events[events["event_type"] == etype]
        if sub.empty:
            continue
        for start_col, end_col in pairs:
            if start_col not in sub.columns:
                continue
            block = sub[["chrom", start_col, end_col, "strand"]].dropna()
            out.update(
                (c, int(s), int(e), st) for c, s, e, st in block.itertuples(index=False)
            )
    return out


#: Columns each consequence mode reads, and the label used in messages.
_CONSEQUENCE_MODES = {
    "exon": (("exon_start", "exon_end"), "cassette exons"),
    "junction": (("start", "end"), "splice-site shifts"),
}

#: Junction classes that a splice-site shift can be resolved against the annotation.
#: ``novel_combination`` is deliberately excluded: both of its sites are annotated, so
#: it is an exon-skipping junction rather than a shifted splice site, and reading it as
#: a shift would report the skipped exon as a "truncation".
SHIFT_CLASSES = ("novel_donor", "novel_acceptor")


def _consequence_table(events, gtf: str, genome: str, mode: str = "exon"):
    """Predict protein consequences for a table of events; ``None`` on a bad genome.

    ``mode="exon"`` treats each row as a cassette exon interval; ``mode="junction"``
    treats it as an intron whose novel splice site is resolved against the annotation
    into the sequence it adds to, or removes from, the neighbouring exon.
    """
    from .consequence import (
        GenomeFasta,
        annotate_consequences,
        annotate_junction_consequences,
        load_transcripts,
    )

    (start_col, end_col), _ = _CONSEQUENCE_MODES[mode]
    if mode == "exon" and "event_type" in events:
        table = events[events["event_type"] == "SE"].copy()
    else:
        table = events.copy()
    if table.empty:
        return table
    genes = set(table["gene_id"].dropna().astype(str)) if "gene_id" in table else None
    transcripts = load_transcripts(gtf, genes=genes or None)
    annotate = annotate_consequences if mode == "exon" else annotate_junction_consequences
    try:
        with GenomeFasta(genome) as fasta:
            return annotate(table, transcripts, fasta, start_col=start_col, end_col=end_col)
    except FileNotFoundError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return None


def _write_consequence(table, out_path, label: str) -> None:
    table.to_csv(out_path, sep="\t", index=False)
    counts = table["consequence_class"].value_counts().to_dict()
    print(f"[consequence] {len(table)} {label} -> {out_path}")
    print(f"[consequence] {counts}")


def _cmd_consequence(args: argparse.Namespace) -> int:
    import pandas as pd

    events = pd.read_csv(args.events, sep="\t")
    mode = args.mode
    if mode == "auto":
        if {"exon_start", "exon_end"} <= set(events.columns):
            mode = "exon"
        elif {"start", "end"} <= set(events.columns):
            mode = "junction"
        else:
            print(
                f"error: {args.events} has neither exon_start/exon_end (cassette exons) "
                "nor start/end (junctions); pass --mode explicitly",
                file=sys.stderr,
            )
            return 2
    (start_col, end_col), label = _CONSEQUENCE_MODES[mode]
    missing = {"chrom", "strand", start_col, end_col} - set(events.columns)
    if missing:
        print(f"error: {args.events} is missing columns {sorted(missing)}", file=sys.stderr)
        return 2
    table = _consequence_table(events, args.gtf, args.genome, mode=mode)
    if table is None:
        return 2
    _write_consequence(table, args.out, label)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="splicescope", description=__doc__.splitlines()[0])
    p.add_argument("--version", action="version", version=f"splicescope {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    s = sub.add_parser("simulate", help="write a synthetic ground-truth dataset")
    s.add_argument("--outdir", required=True)
    s.add_argument("--genes", type=int, default=8)
    s.add_argument("--replicates", type=int, default=4)
    s.add_argument(
        "--alt-ss", type=float, default=0.4, help="fraction of genes with an A5SS/A3SS event"
    )
    s.add_argument("--mxe", type=float, default=0.25, help="fraction of genes with an MXE event")
    s.add_argument("--seed", type=int, default=0)
    s.set_defaults(func=_cmd_simulate)

    r = sub.add_parser("run", help="end-to-end: annotate -> quantify -> diff -> cryptic ML")
    r.add_argument("--sj-dir", required=True, help="directory of STAR *.SJ.out.tab files")
    r.add_argument("--gtf", required=True)
    r.add_argument("--groups", required=True, help="TSV with columns sample,condition")
    r.add_argument("--outdir", required=True)
    r.add_argument("--gene-sets", default=None, help="optional GMT file for pathway enrichment")
    r.add_argument(
        "--genome",
        default=None,
        help="indexed genome FASTA (.fai required); enables protein-consequence prediction",
    )
    r.add_argument("--min-reads", type=int, default=10)
    r.add_argument(
        "--max-exon",
        type=int,
        default=1000,
        help="longest candidate exon considered for mutually-exclusive-exon events",
    )
    r.add_argument("--seed", type=int, default=0)
    r.set_defaults(func=_cmd_run)

    c = sub.add_parser(
        "consequence",
        help="predict frame / PTC / NMD effects for cassette exons or splice-site shifts",
    )
    c.add_argument(
        "--events",
        required=True,
        help="TSV of cassette exons (chrom, strand, exon_start, exon_end) "
        "or of junctions (chrom, strand, start, end)",
    )
    c.add_argument("--gtf", required=True)
    c.add_argument("--genome", required=True, help="indexed genome FASTA (.fai required)")
    c.add_argument("--out", required=True)
    c.add_argument(
        "--mode",
        choices=["auto", "exon", "junction"],
        default="auto",
        help="'exon' reads exon_start/exon_end as a cassette exon; 'junction' reads "
        "start/end as an intron and resolves its novel splice site against the "
        "annotation; 'auto' (default) picks by which columns are present",
    )
    c.set_defaults(func=_cmd_consequence)

    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
