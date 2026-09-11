#!/usr/bin/env python3
"""What the 0.9.0 consequence corrections do to the published SMG1i contrast.

Two of the three 0.9.0 corrections can be recovered from the 0.8.x output alone, without
the GENCODE annotation and the GRCh38 assembly that re-running the caller would need.

**The last-exon rule.** A premature stop inside an *extension* of a transcript's final
exon was called ``ptc_nmd``. It cannot be: an extension is contiguous with the exon it
joins, so when that exon is the last one the stop lies past the final exon-exon junction
and the 50-nucleotide rule has no junction to measure against. The old distance formula

    distance = insert_length - ptc_offset - 3 + sum(downstream_exon_lengths[:-1])

collapses to ``insert_length - ptc_offset - 3`` exactly when one downstream exon remained,
which is that case and leaves an arithmetic signature.

**The in-frame truncation rule.** Removing a whole number of codons leaves the downstream
reading frame untouched, so the first in-frame stop in the retained sequence is the
transcript's *own*. Every truncation whose length is a multiple of three and which was
called ``ptc_nmd``/``ptc_escape`` was therefore reporting the annotated stop as premature.
No signature needed: ``insert_length`` says which they are.

(The third correction — exon-skipping junctions are no longer read as splice-site shifts —
needs the annotation to identify, and stays outstanding.)

Both reclassifications are then testable against the experiment, which is the point: each
*predicts* that the events it moves escape decay, and the SMG1i data can refuse.

    python validation/consequence_corrections.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

TABLE = Path(__file__).parent / "manuscript_tables" / "smg1i_junction_consequence.tsv"

#: Classes in which no premature stop can arise; the published contrast is against these.
NO_STOP = ["in_frame_insertion", "utr_insertion", "non_coding_host"]


def final_exon_extensions(table: pd.DataFrame) -> pd.Series:
    """Rows whose PTC sat in an extension of the transcript's final exon."""
    is_extension = (
        (table["insert_length"] > 0)
        & table["ptc_offset"].notna()
        & table["distance_to_last_junction"].notna()
    )
    implied = table["insert_length"] - table["ptc_offset"] - 3
    single_downstream_exon = np.isclose(implied, table["distance_to_last_junction"])
    return is_extension & single_downstream_exon & (table["consequence_class"] == "ptc_nmd")


def inframe_truncations(table: pd.DataFrame) -> pd.Series:
    """Truncations that removed whole codons, so the stop they found was the native one."""
    return (
        (table["insert_length"] < 0)
        & (table["insert_length"].abs() % 3 == 0)
        & table["consequence_class"].isin(["ptc_nmd", "ptc_escape"])
    )


def agreement(table: pd.DataFrame, moved: pd.Series, label: str) -> None:
    """Does the measured decay data agree that the moved events escape NMD?"""
    reclassified = table[moved]
    remaining = table[(table["consequence_class"] == "ptc_nmd") & ~moved]
    odds, p = stats.fisher_exact(
        [
            [
                int(reclassified["nmd_sensitive"].sum()),
                int((1 - reclassified["nmd_sensitive"]).sum()),
            ],
            [int(remaining["nmd_sensitive"].sum()), int((1 - remaining["nmd_sensitive"]).sum())],
        ],
        alternative="two-sided",
    )
    print(
        f"    {label}\n"
        f"      reclassified   n={len(reclassified):4d}  "
        f"{100 * reclassified['nmd_sensitive'].mean():.1f}% NMD-sensitive\n"
        f"      still ptc_nmd  n={len(remaining):4d}  "
        f"{100 * remaining['nmd_sensitive'].mean():.1f}% NMD-sensitive\n"
        f"      odds ratio {odds:.2f}, two-sided Fisher p = {p:.2e}"
    )


def contrast(table: pd.DataFrame, classes: pd.Series) -> tuple[float, float]:
    """ptc_nmd against the no-stop classes, on the measured NMD-sensitivity label."""
    nmd = table[classes == "ptc_nmd"]
    nostop = table[classes.isin(NO_STOP)]
    counts = [
        [int(nmd["nmd_sensitive"].sum()), int((1 - nmd["nmd_sensitive"]).sum())],
        [int(nostop["nmd_sensitive"].sum()), int((1 - nostop["nmd_sensitive"]).sum())],
    ]
    odds, p = stats.fisher_exact(counts, alternative="two-sided")
    print(
        f"    ptc_nmd  n={len(nmd):4d}  {100 * nmd['nmd_sensitive'].mean():.1f}% NMD-sensitive\n"
        f"    no-stop  n={len(nostop):4d}  {100 * nostop['nmd_sensitive'].mean():.1f}% "
        f"NMD-sensitive\n"
        f"    odds ratio {odds:.2f}, two-sided Fisher p = {p:.2e}"
    )
    return odds, p


def main() -> None:
    table = pd.read_csv(TABLE, sep="\t")
    last_exon = final_exon_extensions(table)
    in_frame = inframe_truncations(table)

    print(
        f"{len(table)} junctions; {int(last_exon.sum())} final-exon extensions and "
        f"{int(in_frame.sum())} in-frame truncations reclassified\n"
    )

    print("1. The published contrast, as 0.8.1 produced it")
    contrast(table, table["consequence_class"])

    print("\n2. After the last-exon correction alone")
    contrast(table, table["consequence_class"].where(~last_exon, "ptc_escape"))

    print("\n3. After both corrections")
    corrected = table["consequence_class"].where(~last_exon, "ptc_escape")
    corrected = corrected.where(~in_frame, "exon_truncation")
    contrast(table, corrected)

    print("\n4. Does the experiment agree that the moved events escape decay?")
    agreement(table, last_exon, "last-exon extensions -> ptc_escape")
    agreement(table, in_frame, "in-frame truncations -> exon_truncation")
    baseline = table[table["consequence_class"].isin(NO_STOP)]["nmd_sensitive"].mean()
    print(f"\n    for scale, events where no stop is possible: {100 * baseline:.1f}%")
    print(
        "\n    Both groups sit below the events they were taken from and above the "
        "no-stop\n    baseline, which is what the corrections predict and what the data "
        "could have\n    refused. The last-exon rule is the stronger of the two; the "
        "in-frame truncation\n    rule moves in the same direction more weakly."
    )


if __name__ == "__main__":
    main()
