#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_pooled_pcoa_density.py

Re-render of
results/single_pair/03_beta/pcoa_permanova/plots/pooled/pcoa_pooled_paired_adonis2_weighted_unifrac.png
with the scatter replaced by the DENSITY-CONTOUR style of
/home/hexincheng/projects/hxc/oral_project/beta_diversity/beta_pcoa_20260715/04-plot-density-contour-pcoa.R

The original generator script is not in this repo (the PCoA/PERMANOVA run
was executed under /home/zhangcaiyi/project/oral_gut_pair), but all inputs
are present here, so the figure is recomputed from data:

  * pooled paired PCoA coordinates (dataset_label=filter,
    metric=weighted_unifrac, analysis_scope=pooled_all_bioprojects)
  * PERMANOVA (adonis2 D ~ BioProject + site, pair-blocked permutation)
    stats from pcoa_permanova_pooled.tsv
  * axis % variance from pcoa_variance_explained.tsv

Density style (ported from the R script):
  * one filled 2D-KDE contour layer per site, 6 bins, alpha=0.15 per band,
    same-colour 0.3-width contour lines; bands are nested (region ABOVE each
    level is filled), so opacity accumulates toward the mode exactly like
    ggplot2::stat_density_2d(geom="polygon", bins=6, alpha=0.15)
  * Oral drawn first (background), Gut second (foreground)
  * white plot background, thin black frame; zero lines + dotted grid kept
    from the original figure; PERMANOVA annotation box bottom-left

Inputs : results/single_pair/03_beta/pcoa_permanova/pcoa_coordinates_pooled.tsv
         results/single_pair/03_beta/pcoa_permanova/pcoa_permanova_pooled.tsv
         results/single_pair/03_beta/pcoa_permanova/pcoa_variance_explained.tsv
Outputs: results/new/figures/11_pooled_pcoa_density/
         pcoa_pooled_paired_adonis2_weighted_unifrac_density.{png,pdf}
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

ROOT = Path("/home/hexincheng/projects/hxc/oral_project/oral_gut_pair")
IN_DIR = ROOT / "results/single_pair/03_beta/pcoa_permanova"
OUT_DIR = ROOT / "results/new/figures/11_pooled_pcoa_density"

DATASET_LABEL = "filter"
METRIC = "weighted_unifrac"
METRIC_TITLE = "Weighted UniFrac"
SITE_COLORS = {"Oral": "#D65F4B", "Gut": "#63B6C7"}
N_BINS = 6          # R: stat_density_2d(bins = 6)
BAND_ALPHA = 0.15   # R: alpha = 0.15

plt.rcParams.update(
    {
        "pdf.fonttype": 42,
        "font.size": 9,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
    }
)


def draw_density_layer(ax, x: np.ndarray, y: np.ndarray, color: str,
                       n_bins: int = N_BINS, alpha: float = BAND_ALPHA) -> None:
    """ggplot2 stat_density_2d equivalent: nested filled KDE contours.

    Each of the (n_bins - 1) interior levels fills the region where the
    density is ABOVE that level with `alpha`; stacked bands therefore
    accumulate opacity toward the mode, matching the R output.
    """
    kde = gaussian_kde(np.vstack([x, y]))
    xg = np.linspace(x.min(), x.max(), 200)
    yg = np.linspace(y.min(), y.max(), 200)
    xx, yy = np.meshgrid(xg, yg)
    zz = kde(np.vstack([xx.ravel(), yy.ravel()])).reshape(xx.shape)

    zmax = float(zz.max())
    levels = np.linspace(0.0, zmax, n_bins + 1)[1:-1]  # interior breakpoints
    for lv in levels:
        ax.contourf(
            xx, yy, zz,
            levels=[lv, zmax * 1.0001],
            colors=[color],
            alpha=alpha,
            zorder=2,
        )
    ax.contour(
        xx, yy, zz,
        levels=np.concatenate([levels, [zmax * 0.9999]]),
        colors=[color],
        linewidths=0.3,
        zorder=2.5,
    )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    coords = pd.read_csv(
        IN_DIR / "pcoa_coordinates_pooled.tsv", sep="\t",
        usecols=["PC1", "PC2", "site_for_analysis", "dataset_label",
                 "metric", "analysis_scope"],
    )
    d = coords.loc[
        coords["dataset_label"].astype(str).eq(DATASET_LABEL)
        & coords["metric"].astype(str).eq(METRIC)
        & coords["analysis_scope"].astype(str).eq("pooled_all_bioprojects")
    ].copy()
    assert not d.empty, "no pooled rows matched"

    perm = pd.read_csv(IN_DIR / "pcoa_permanova_pooled.tsv", sep="\t")
    row = perm.loc[
        perm["dataset_label"].astype(str).eq(DATASET_LABEL)
        & perm["metric"].astype(str).eq(METRIC)
        & perm["analysis_scope"].astype(str).eq("pooled_all_bioprojects")
    ].iloc[0]

    ve = pd.read_csv(IN_DIR / "pcoa_variance_explained.tsv", sep="\t")
    ve_row = ve.loc[
        ve["dataset_label"].astype(str).eq(DATASET_LABEL)
        & ve["metric"].astype(str).eq(METRIC)
        & ve["analysis_scope"].astype(str).eq("pooled_all_bioprojects")
    ].iloc[0]
    pc1_pct = float(ve_row["pc1_percent"])
    pc2_pct = float(ve_row["pc2_percent"])

    # canvas shrunk to smallest legible size (user request: as small as possible);
    # fonts dropped to the fig 05/06 convention (title 10, axes 8).
    # Unified plot-area geometry (user request 2026-10-20): canvas 3.40 in
    # tall, axes 2.26 in with bottom edge 0.70 in / top edge 2.96 in from
    # the canvas bottom — identical to figs 03/05/06/07 — saved WITHOUT
    # bbox_inches="tight"; title uses pad=4 so its two lines stay inside
    # the 0.44-in top band
    fig = plt.figure(figsize=(3.30, 3.40))
    ax = fig.add_axes([0.60 / 3.30, 0.70 / 3.40,
                       (3.30 - 0.60 - 0.08) / 3.30, 2.26 / 3.40])

    for site in ("Oral", "Gut"):   # Oral = background layer, Gut = foreground
        s = d.loc[d["site_for_analysis"].astype(str).eq(site)]
        draw_density_layer(ax, s["PC1"].to_numpy(float), s["PC2"].to_numpy(float),
                           SITE_COLORS[site])

    ax.axhline(0, color="#c8c8c8", linewidth=0.9, zorder=1)
    ax.axvline(0, color="#c8c8c8", linewidth=0.9, zorder=1)
    ax.grid(linestyle=":", linewidth=0.5, alpha=0.4)
    ax.set_axisbelow(True)

    p_txt = f"P={float(row['p_value']):.4f}"
    ax.text(
        0.03, 0.06,
        (f"adonis2: D ~ BioProject + site\n"
         f"pair-blocked permutation\n"
         f"R²(site)={float(row['r2_term']):.4f}, "
         f"F={float(row['pseudo_F']):.3f}, {p_txt}"),
        transform=ax.transAxes, ha="left", va="bottom", fontsize=6.5,
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
                  edgecolor="#999", alpha=0.92),
    )

    ax.tick_params(labelsize=8)
    ax.set_xlabel(f"PC1 ({pc1_pct:.2f}%)", fontsize=8)
    ax.set_ylabel(f"PC2 ({pc2_pct:.2f}%)", fontsize=8)
    ax.set_title(f"Pooled paired samples\n{METRIC_TITLE} PCoA", fontsize=10,
                 pad=4)

    from matplotlib.patches import Patch

    handles = [
        Patch(facecolor=SITE_COLORS["Oral"], alpha=0.55, edgecolor="none", label="Oral"),
        Patch(facecolor=SITE_COLORS["Gut"], alpha=0.55, edgecolor="none", label="Gut"),
    ]
    ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=7.5,
              title=None)

    for ext in ("png", "pdf"):
        fig.savefig(
            OUT_DIR / f"pcoa_pooled_paired_adonis2_weighted_unifrac_density.{ext}",
            dpi=300,   # no bbox_inches="tight" — unified geometry, see above
        )
    plt.close(fig)
    print(
        f"n={len(d)} samples (Oral {int((d['site_for_analysis'] == 'Oral').sum())}, "
        f"Gut {int((d['site_for_analysis'] == 'Gut').sum())}); "
        f"R2={float(row['r2_term']):.4f}, F={float(row['pseudo_F']):.3f}, "
        f"{p_txt}; PC1 {pc1_pct:.2f}%, PC2 {pc2_pct:.2f}%"
    )
    print(f"outputs -> {OUT_DIR}")


if __name__ == "__main__":
    main()