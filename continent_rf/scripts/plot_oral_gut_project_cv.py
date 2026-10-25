#!/usr/bin/env python3
"""Plot the oral/gut project-level OVO two-panel AUC figure.

Companion to run_oral_gut_project_cv_repeat5.py.
Reads continent_rf/result/oral/.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


OUTDIR = (Path(__file__).resolve().parents[1] / "result" / "oral")
REGIONS = ["North America", "Asia", "Europe", "Oceania", "South America", "Africa"]


def _ovo_pair_frame(summary, site):
    frame = summary[summary["site"] == site]
    classes = [
        r for r in REGIONS if r in set(frame["region_1"]) | set(frame["region_2"])
    ]
    matrix = pd.DataFrame(np.nan, index=classes, columns=classes)
    for row in frame.itertuples(index=False):
        matrix.loc[row.region_1, row.region_2] = row.mean_auc
        matrix.loc[row.region_2, row.region_1] = row.mean_auc
    return matrix


def plot_ovo_saliva_gut_panels_upper_triangle():
    """Oral (saliva) and Gut OVO heatmaps side by side in a single figure.

    Same pair data as the combined triangle figure, one panel per dataset.
    Lower-triangle layout: each panel keeps only pairs (col=earlier, row=later
    continent); North America never appears as the later member and Africa
    never as the earlier one, so the empty top row and right column are
    cropped — leaving a compact 5x5 staircase per panel. Axes are positioned
    manually (no tight_layout) so the two panels sit as close as the row
    labels of the right panel allow. Rows/cols are pair members, no train/test
    roles.
    """
    summary = pd.read_csv(OUTDIR / "region_ovo_summary.tsv", sep="\t")
    sites = [s for s in ("Saliva", "Gut") if s in set(summary["site"])]
    if not sites:
        return
    # colour ramp starts at the lowest displayed AUC (floor to 0.05), not 0
    vmin = np.floor(
        min(summary.loc[summary.site == s, "mean_auc"].min() for s in sites) * 20
    ) / 20
    # cell size just large enough for the printed AUC values (~0.46in)
    fig = plt.figure(figsize=(7.0, 3.3))
    panel = (2.3 / 7.0, 2.3 / 3.3)             # 2.3in square panels, 5 cells each
    x0 = 0.10                                  # left panel's row labels
    y0 = 0.22                                  # room for rotated column labels
    gap = 0.121                                # right panel's row labels
    cax = fig.add_axes((0.93, y0, 0.012, panel[1]))
    for k, site in enumerate(sites):
        ax = fig.add_axes((x0 + k * (panel[0] + gap), y0, panel[0], panel[1]))
        matrix = _ovo_pair_frame(summary, site)
        classes = [r for r in REGIONS if r in matrix.index]
        matrix = matrix.reindex(index=classes, columns=classes)
        matrix = matrix.iloc[1:, :-1]          # drop empty top row / right col
        labels = matrix.map(lambda v: f"{v * 100:.1f}" if pd.notna(v) else "")
        # after cropping the top row / right column, the real pair cells sit
        # on and below the diagonal of the cropped 5x5 grid: hide strictly-above
        mask = np.triu(np.ones(matrix.shape, dtype=bool), k=1)
        sns.heatmap(
            matrix,
            mask=mask,
            annot=labels,
            fmt="",
            annot_kws={"size": 8.5, "color": "white"},
            cmap="Blues",
            vmin=vmin,
            vmax=1,
            linewidths=0.8,
            linecolor="white",
            cbar=False,
            square=True,
            xticklabels=classes[:-1],
            yticklabels=classes[1:],
            ax=ax,
        )
        displayed = "Oral (saliva)" if site == "Saliva" else site
        ax.set_title(displayed, fontsize=11, pad=8)
        ax.tick_params(labelsize=7.5)
        ax.tick_params(axis="x", rotation=45)
        ax.tick_params(axis="y", pad=1)
        # anchor the x labels' right end on the tick (left of centre alignment)
        for lbl in ax.get_xticklabels():
            lbl.set_rotation_mode("anchor")
            lbl.set_ha("right")
        # dark digits on the pale low-AUC cells, white on the dark ones
        norm_vals = (matrix.to_numpy() - vmin) / (1 - vmin)
        shown = norm_vals[~mask]
        for t, frac in zip(ax.texts, shown):
            t.set_color("#333333" if frac < 0.30 else "white")
    sm = plt.cm.ScalarMappable(cmap="Blues", norm=plt.Normalize(vmin, 1))
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax, label="Mean AUC (%)")
    cbar.ax.yaxis.set_major_formatter(
        matplotlib.ticker.PercentFormatter(xmax=1.0, decimals=0, symbol="")
    )
    fig.savefig(OUTDIR / "region_ovo_oral_gut_two_panels_upper_triangle.png",
                dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / "region_ovo_oral_gut_two_panels_upper_triangle.pdf",
                bbox_inches="tight")
    plt.close(fig)


def main():
    if not (OUTDIR / "region_ovo_summary.tsv").exists():
        raise FileNotFoundError(
            f"Missing {OUTDIR / 'region_ovo_summary.tsv'}; "
            "run run_oral_gut_project_cv_repeat5.py first"
        )
    OUTDIR.mkdir(parents=True, exist_ok=True)
    plot_ovo_saliva_gut_panels_upper_triangle()
    print(f"Plots written to {OUTDIR}")


if __name__ == "__main__":
    main()
