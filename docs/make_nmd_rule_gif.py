"""Render docs/nmd_rule.gif — the 50-nucleotide rule deciding a transcript's fate.

A premature stop codon walks towards the last exon-exon junction. Past 50 nucleotides
the ribosome stops while the junction complex is still downstream and the transcript is
degraded; inside 50 it does not, and a truncated protein is made instead. That threshold
is the whole of what this toolkit predicts about a cryptic exon, and it is one number.

Every verdict in the animation comes from `_nmd_from_downstream` — the function the
pipeline uses — rather than from a drawing of the rule, so the picture cannot say
something the code does not.

    python docs/make_nmd_rule_gif.py
"""

from __future__ import annotations

import io
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
from PIL import Image

from splicescope.consequence import NMD_DISTANCE_RULE, _nmd_from_downstream

OUT = Path(__file__).resolve().parent / "nmd_rule.gif"

BG = "#07060f"
FG = "#f4f1ff"
MUTED = "#a39fc4"
DECAY = "#ff3d81"
ESCAPE = "#22d3ee"
EXON = "#8b5cf6"

#: The retained transcript downstream of the event: two exons, the junction between them.
FIRST_EXON, LAST_EXON = 460, 150
TOTAL = FIRST_EXON + LAST_EXON


def verdict(stop_at: int):
    """Ask the pipeline's own function where the stop is and what it means."""
    filler = "AAG" * (TOTAL // 3 + 2)
    sequence = filler[:stop_at] + "TAA" + filler[stop_at + 3 : TOTAL]
    offset, distance, nmd = _nmd_from_downstream(
        sequence, [FIRST_EXON, LAST_EXON], frame=0, offset_before=0
    )
    return offset, distance, nmd


def frame(stop_at: int) -> Image.Image:
    offset, distance, nmd = verdict(stop_at)
    colour = DECAY if nmd else ESCAPE
    fig, ax = plt.subplots(figsize=(7.2, 3.6), dpi=110)
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)

    # the two retained exons, and the junction between them
    ax.add_patch(Rectangle((0, 0.42), FIRST_EXON, 0.16, color=EXON, alpha=0.85))
    ax.add_patch(Rectangle((FIRST_EXON + 8, 0.42), LAST_EXON - 8, 0.16, color=EXON, alpha=0.85))
    ax.axvline(FIRST_EXON + 4, color=FG, lw=2, ymin=0.30, ymax=0.72)
    ax.text(FIRST_EXON + 4, 0.74, "last exon–exon junction", color=FG, fontsize=9,
            ha="center", va="bottom")

    # where the rule sits: 50 nt before the junction, which is the line the stop crosses
    threshold = FIRST_EXON + 4 - NMD_DISTANCE_RULE
    ax.axvline(threshold, color=MUTED, lw=1.2, ls=(0, (4, 3)), ymin=0.22, ymax=0.66)
    ax.text(threshold, 0.13, f"{NMD_DISTANCE_RULE} nt", color=MUTED, fontsize=9,
            ha="center", va="top")

    # the premature stop, labelled above so nothing collides with the junction
    ax.plot([offset], [0.50], marker="v", markersize=13, color=colour, zorder=5)
    ax.text(offset, 0.62, "premature stop", color=colour, fontsize=9, ha="center", va="bottom")

    # the distance between them, which is the whole of the rule
    ax.add_patch(
        FancyArrowPatch(
            (offset + 3, 0.28), (FIRST_EXON + 1, 0.28),
            arrowstyle="<->", mutation_scale=11, color=colour, lw=1.4, alpha=0.9,
        )
    )
    ax.text(offset - 12, 0.28, f"{distance} nt", color=colour, fontsize=11, weight="bold",
            ha="right", va="center")

    ax.text(
        0.5, 0.92,
        "degraded by nonsense-mediated decay" if nmd else "escapes decay — truncated protein",
        transform=ax.transAxes, color=colour, fontsize=14, weight="bold", ha="center",
    )
    ax.text(
        0.5, 0.05,
        f"more than {NMD_DISTANCE_RULE} nt from the junction"
        if nmd
        else f"within {NMD_DISTANCE_RULE} nt of the junction",
        transform=ax.transAxes, color=MUTED, fontsize=10, ha="center",
    )

    ax.set_xlim(-25, TOTAL + 25)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.tight_layout()

    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", facecolor=BG)
    plt.close(fig)
    buffer.seek(0)
    return Image.open(buffer).convert("RGB")


def main() -> None:
    # the stop walks towards the junction; the call flips exactly once
    positions = list(range(300, FIRST_EXON - 9, 12))
    frames = [frame(p) for p in positions]
    distances = [verdict(p)[1] for p in positions]
    flips = [
        (a, b) for a, b in zip(distances, distances[1:], strict=False)
        if (a > NMD_DISTANCE_RULE) != (b > NMD_DISTANCE_RULE)
    ]
    print(f"distances {distances[0]} down to {distances[-1]} nt; the call flips at {flips}")

    frames = frames + frames[-2:0:-1]  # ping-pong
    frames[0].save(OUT, save_all=True, append_images=frames[1:], duration=220, loop=0,
                   optimize=True)
    print(f"wrote {OUT} ({OUT.stat().st_size // 1024} KB, {len(frames)} frames)")


if __name__ == "__main__":
    main()
