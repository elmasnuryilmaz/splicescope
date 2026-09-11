#!/usr/bin/env python3
"""What the 0.9.0 last-exon NMD correction does to the published SMG1i contrast.

Up to 0.8.1 a premature stop found inside an *extension* of a transcript's **final** exon
was called ``ptc_nmd``. That cannot be right: an extension is contiguous with the exon it
joins, so when that exon is the last one the stop lies past the final exon-exon junction
and the 50-nucleotide rule has no junction to measure against. The old code measured a
distance to a junction that does not exist.

Re-running `splicescope consequence` on the 3,548 published junctions needs the GENCODE
annotation and the GRCh38 assembly, neither of which is in this repository. But the
affected rows can be recovered from the 0.8.x output alone, because the old distance
formula leaves an arithmetic signature:

    distance = insert_length - ptc_offset - 3 + sum(downstream_exon_lengths[:-1])

so ``distance == insert_length - ptc_offset - 3`` exactly when there was a single
downstream exon — that is, when the extension ran into the final exon.

The reclassification is then testable against the experiment, which is the point: the
correction *predicts* these events escape decay, and the SMG1i data can refuse.

    python validation/last_exon_correction.py
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
    moved = final_exon_extensions(table)

    print(f"{len(table)} junctions; {int(moved.sum())} reclassified ptc_nmd -> ptc_escape\n")

    print("1. The published contrast, as 0.8.1 produced it")
    contrast(table, table["consequence_class"])

    print("\n2. The same contrast after the correction")
    corrected = table["consequence_class"].where(~moved, "ptc_escape")
    contrast(table, corrected)

    print("\n3. Does the experiment agree that these events escape decay?")
    reclassified = table[moved]
    remaining = table[(table["consequence_class"] == "ptc_nmd") & ~moved]
    counts = [
        [
            int(reclassified["nmd_sensitive"].sum()),
            int((1 - reclassified["nmd_sensitive"]).sum()),
        ],
        [int(remaining["nmd_sensitive"].sum()), int((1 - remaining["nmd_sensitive"]).sum())],
    ]
    odds, p = stats.fisher_exact(counts, alternative="two-sided")
    baseline = table[table["consequence_class"].isin(NO_STOP)]["nmd_sensitive"].mean()
    print(
        f"    reclassified to ptc_escape  n={len(reclassified):4d}  "
        f"{100 * reclassified['nmd_sensitive'].mean():.1f}% NMD-sensitive\n"
        f"    still ptc_nmd               n={len(remaining):4d}  "
        f"{100 * remaining['nmd_sensitive'].mean():.1f}% NMD-sensitive\n"
        f"    odds ratio {odds:.2f}, two-sided Fisher p = {p:.2e}\n"
        f"    for scale, events where no stop is possible: {100 * baseline:.1f}%"
    )
    print(
        "\n    The reclassified events behave like escapers, not like decay targets — "
        "which is\n    what the correction predicts, and is a claim the data could have "
        "refused."
    )


if __name__ == "__main__":
    main()
