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
        n_genes=args.genes, n_per_group=args.replicates, seed=args.seed
    )
    out = write_dataset(ds, args.outdir)
    n_junc = ds.observed.drop_duplicates(["chrom", "start", "end", "strand"]).shape[0]
    print(f"[simulate] wrote {n_junc} unique junctions for {len(ds.groups)} samples -> {out}")
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

    # event-level: cassette-exon PSI and its differential inclusion
    from . import events as _events

    evs = _events.detect_cassette_events(annotated)
    if not evs.empty:
        evs.to_csv(outdir / "cassette_events.tsv", sep="\t", index=False)
        cpsi = _events.cassette_psi(annotated, evs, min_reads=args.min_reads)
        cdiff = _diff.differential_splicing(cpsi, groups, value="psi_cassette")
        cdiff.to_csv(outdir / "cassette_differential.tsv", sep="\t", index=False)
        print(
            f"[run] {len(evs)} cassette-exon events; "
            f"{len(_diff.significant(cdiff))} differentially spliced (ΔΨ)"
        )
        if not cdiff.empty:
            fig, ax = _plot.plt.subplots(figsize=(5, 4))
            _plot.plot_volcano(cdiff, ax=ax)
            ax.set_title("Differential cassette-exon inclusion")
            _plot.savefig(fig, figs / "cassette_volcano.png")

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


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="splicescope", description=__doc__.splitlines()[0])
    p.add_argument("--version", action="version", version=f"splicescope {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    s = sub.add_parser("simulate", help="write a synthetic ground-truth dataset")
    s.add_argument("--outdir", required=True)
    s.add_argument("--genes", type=int, default=8)
    s.add_argument("--replicates", type=int, default=4)
    s.add_argument("--seed", type=int, default=0)
    s.set_defaults(func=_cmd_simulate)

    r = sub.add_parser("run", help="end-to-end: annotate -> quantify -> diff -> cryptic ML")
    r.add_argument("--sj-dir", required=True, help="directory of STAR *.SJ.out.tab files")
    r.add_argument("--gtf", required=True)
    r.add_argument("--groups", required=True, help="TSV with columns sample,condition")
    r.add_argument("--outdir", required=True)
    r.add_argument("--min-reads", type=int, default=10)
    r.add_argument("--seed", type=int, default=0)
    r.set_defaults(func=_cmd_run)

    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
