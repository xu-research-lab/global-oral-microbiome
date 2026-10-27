#!/usr/bin/env python3
"""Plot disease age deviations from HC for each oral site."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from age_common import DISEASE_ABBR, FIGURE_DIR, load_metadata
from age_stats import disease_age_gap_stats

OUT_DIR = FIGURE_DIR


def main():
    stats = disease_age_gap_stats(load_metadata())
    if stats.empty:
        raise ValueError("No disease/HC comparisons are available")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    stats.to_csv(OUT_DIR / "age_disease_microbiome_age_acceleration_1row_stats.csv", index=False)

    from matplotlib.gridspec import GridSpec
    fig = plt.figure(figsize=(12, 4.3))
    gs = GridSpec(3, 2, figure=fig, hspace=0.45, wspace=0.75,
                  left=0.16, right=0.82, top=0.96, bottom=0.08,
                  width_ratios=[1.0, 1.0])

    def draw_panel(ax, j):
        s = stats[stats.site == j].sort_values("delta")
        if len(s) == 0:
            ax.set_axis_off()
            ax.set_title(f"{j} — no data")
            return
        ax.axvline(0, color="k", ls="--", lw=0.9, zorder=1)

        y = np.arange(len(s))
        col = np.where(s.fdr < .05, np.where(s.delta > 0, "#c0392b", "#2471a3"), "0.55")
        ax.hlines(y, s.ci_lo, s.ci_hi, color=col, lw=1.6, zorder=2)
        ax.scatter(s.delta, y, s=40 + 3 * np.sqrt(s.n), c=col,
                   edgecolor="k", lw=.5, zorder=3)

        for yi, row in zip(y, s.itertuples()):
            # Both labels sit OUTSIDE the axes (x>1.0 in axes coords) so they
            # never overlap the data points. n= is left-aligned at a fixed x,
            # and the sig stars are left-aligned at a fixed x just to the right
            # — so within a figure column every n= starts at the same place
            # and every star column starts at the same place, with a small gap.
            ax.text(1.02, yi, f"n={row.n}",
                    transform=ax.get_yaxis_transform(),
                    va="center", ha="left", fontsize=11, color="0.3")
            ax.text(1.25, yi, f"{row.sig}",
                    transform=ax.get_yaxis_transform(),
                    va="center", ha="left", fontsize=11, color="0.3")

        ax.set_yticks(y)
        ax.set_yticklabels([DISEASE_ABBR.get(v, v) for v in s.status], fontsize=12)
        ax.set_ylim(-0.7, len(s) - 0.3)
        # compress y axis vertically: less whitespace above/below points
        ax.margins(y=0.02)
        ax.set_xlabel("Age deviation relative to HC (years)", fontsize=12)
        ax.set_title(j,
                     fontsize=13, fontweight="normal")
        ax.spines[["top", "right"]].set_visible(False)

    # saliva spans the full left column (rows 0:3, col 0)
    ax_saliva = fig.add_subplot(gs[:, 0])
    draw_panel(ax_saliva, "saliva")

    # right column: tooth / oral mucosa / tongue stacked top-to-bottom
    # only the last (bottom) panel shows the x label to avoid overlap
    right_sites = ["tooth", "oral mucosa", "tongue"]
    for i, j in enumerate(right_sites):
        ax = fig.add_subplot(gs[i, 1])
        draw_panel(ax, j)
        if i < len(right_sites) - 1:
            ax.set_xlabel("")

    for ext in ("png", "pdf"):
        fig.savefig(OUT_DIR / f"age_disease_microbiome_age_acceleration_2col.{ext}",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUT_DIR / 'age_disease_microbiome_age_acceleration_2col.png'}")
    print(f"Stats: {OUT_DIR / 'age_disease_microbiome_age_acceleration_1row_stats.csv'}")


if __name__ == "__main__":
    main()
