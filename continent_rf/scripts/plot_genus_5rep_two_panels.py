#!/usr/bin/env python3
"""Reproduce the region_ovo_oral_gut_two_panels_upper_triangle figure from the
genus-level 5-rep results (result_genus_5rep/), same style as
plot_oral_gut_project_cv.plot_ovo_saliva_gut_panels_upper_triangle().

Reads result_genus_5rep/{oral,gut}/region_ovo_summary.tsv and writes
  result_genus_5rep/region_ovo_oral_gut_two_panels_upper_triangle.{png,pdf}
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

OUTDIR = Path(__file__).resolve().parents[1] / "result_genus_5rep"
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


def main():
    summary = pd.concat([
        pd.read_csv(OUTDIR / d / "region_ovo_summary.tsv", sep="\t")
        for d in ("oral", "gut")
    ])
    sites = [s for s in ("Saliva", "Gut") if s in set(summary["site"])]
    vmin = np.floor(
        min(summary.loc[summary.site == s, "mean_auc"].min() for s in sites) * 20
    ) / 20
    fig = plt.figure(figsize=(7.0, 3.3))
    panel = (2.3 / 7.0, 2.3 / 3.3)
    x0 = 0.10
    y0 = 0.22
    gap = 0.121
    cax = fig.add_axes((0.93, y0, 0.012, panel[1]))
    for k, site in enumerate(sites):
        ax = fig.add_axes((x0 + k * (panel[0] + gap), y0, panel[0], panel[1]))
        matrix = _ovo_pair_frame(summary, site)
        classes = [r for r in REGIONS if r in matrix.index]
        matrix = matrix.reindex(index=classes, columns=classes)
        matrix = matrix.iloc[1:, :-1]
        labels = matrix.map(lambda v: f"{v * 100:.1f}" if pd.notna(v) else "")
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
        for lbl in ax.get_xticklabels():
            lbl.set_rotation_mode("anchor")
            lbl.set_ha("right")
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
    print(f"[done] -> {OUTDIR}/region_ovo_oral_gut_two_panels_upper_triangle.{{png,pdf}}")


if __name__ == "__main__":
    main()