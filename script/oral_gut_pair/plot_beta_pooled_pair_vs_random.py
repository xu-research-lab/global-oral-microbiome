#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_beta_pooled_pair_vs_random.py

Optimized re-render of the POOLED paired-vs-random beta diversity boxplots:

    results/single_pair/03_beta/pair_vs_random/plots/
        beta_pooled_pair_vs_random_filter_unweighted_unifrac.png
        beta_pooled_pair_vs_random_filter_weighted_unifrac.png

Original generator: s04_2_beta_diversity_pair_vs_random.py (plot_pooled_boxplot).
This script re-plots the SAME precomputed tables — no data recalculation:
every point is one focal oral unit; "Random-unpaired" is the median of 100
random gut draws for that unit (already computed upstream).

Improvements over the original figure (data unchanged):
  * vector PDF output (pdf.fonttype 42) alongside 300-dpi PNG
  * violin plots (KDE clamped to the observed range) with slim inner boxes
    showing median + IQR + whiskers, replacing the boxplot + jitter cloud
  * stats (Wilcoxon P, n, median difference) as a single compact line at the
    top, below the title
  * significance bracket between the two groups (end ticks pointing down)
  * cleaner frame: subtle y-grid, metric-specific y-axis label

Inputs : results/single_pair/03_beta/pair_vs_random/
             beta_pair_vs_random_long_plot_table.tsv
             beta_pair_vs_random_wilcoxon_stats.tsv
Outputs: results/new/figures/05_beta_pair_vs_random/
             beta_pooled_pair_vs_random_filter_unweighted_unifrac.{png,pdf}
             beta_pooled_pair_vs_random_filter_weighted_unifrac.{png,pdf}
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

ROOT = Path("/home/hexincheng/projects/hxc/oral_project/oral_gut_pair")
IN_DIR = ROOT / "results/single_pair/03_beta/pair_vs_random"
IN_LONG = IN_DIR / "beta_pair_vs_random_long_plot_table.tsv"
IN_STATS = IN_DIR / "beta_pair_vs_random_wilcoxon_stats.tsv"
OUT_DIR = ROOT / "results/new/figures/05_beta_pair_vs_random"

DATASET_LABEL = "filter"
METRICS = {
    "unweighted_unifrac": {
        "title": "Unweighted UniFrac distance",
        "ylabel": "Unweighted UniFrac distance",
    },
    "weighted_unifrac": {
        "title": "Weighted UniFrac distance",
        "ylabel": "Weighted UniFrac distance",
    },
}
PAIRED_COLOR = "#D65F4B"
RANDOM_COLOR = "#63B6C7"
X_LABELS = ["Paired\n(same pair_id)", "Random-unpaired\n(median of 100 gut)"]

plt.rcParams.update(
    {
        "pdf.fonttype": 42,
        "font.size": 9,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
    }
)


def format_p(p: float) -> str:
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4f}"


def draw_violin(ax, x: float, vals: np.ndarray, color: str, width: float = 0.6) -> None:
    """Half-free violin (KDE clamped to the observed data range) + slim inner box."""
    kde = gaussian_kde(vals)
    y = np.linspace(vals.min(), vals.max(), 200)
    d = kde(y)
    d = d / d.max() * (width / 2.0)
    ax.fill_betweenx(y, x - d, x + d, facecolor=color, edgecolor="black",
                     linewidth=0.8, alpha=0.85, zorder=2)
    ax.boxplot(
        [vals], positions=[x], widths=0.12, patch_artist=True, showfliers=False,
        boxprops={"facecolor": "white", "edgecolor": "black", "linewidth": 0.8},
        medianprops={"color": "black", "linewidth": 1.4},
        whiskerprops={"linewidth": 0.8, "color": "black"},
        capprops={"linewidth": 0.8, "color": "black"},
        zorder=3,
    )


def draw_panel(long_df: pd.DataFrame, stat_df: pd.DataFrame, metric: str, cfg: dict) -> None:
    d = long_df.loc[
        (long_df["dataset_label"] == DATASET_LABEL) & (long_df["metric"] == metric)
    ].copy()
    if d.empty:
        print(f"WARNING: no pooled data for {DATASET_LABEL} / {metric}")
        return

    paired = d.loc[d["distance_type"] == "Paired", "distance"].astype(float).values
    random = d.loc[d["distance_type"] == "Random-unpaired", "distance"].astype(float).values

    row = stat_df.loc[
        (stat_df["scope"] == "pooled_all_bioprojects")
        & (stat_df["dataset_label"] == DATASET_LABEL)
        & (stat_df["metric"] == metric)
    ].iloc[0]
    p_text = format_p(float(row["p_value"]))
    stars = str(row["p_significance_label"])
    n_units = int(row["n_focal_units"])
    median_diff = float(row["median_difference_paired_minus_random"])

    # canvas shrunk to near-minimum (user request: as small as possible):
    # 3.0x3.0 in with fonts dropped ~1 pt; xlim tightened to (-0.42, 1.08)
    # (violins keep 0.12 data-unit side margins) so the two long x labels
    # still clear each other — with the old (-0.55, 1.3) they collided
    # below 3.6 in width
    fig, ax = plt.subplots(figsize=(3.0, 3.4))

    # --- violins (KDE clamped to the data range) with slim inner boxes ---
    draw_violin(ax, 0, paired, PAIRED_COLOR)
    draw_violin(ax, 0.75, random, RANDOM_COLOR)

    ax.set_xlim(-0.42, 1.17)

    # --- significance bracket ---
    y_min = min(paired.min(), random.min())
    y_max = max(paired.max(), random.max())
    span = y_max - y_min
    # bracket hugging the violins (user request 2026-10-07: the old
    # +0.045·span bracket under a +0.11·span head margin floated too high)
    ax.set_ylim(y_min - 0.03 * span, y_max + 0.12 * span)
    y_br = y_max + 0.028 * span
    tick = 0.016 * span
    # bracket ends back ON the violin centre lines (user request 2026-10-07:
    # 括号要在中线 — the ±0.135 widening put the ticks off the whisker lines)
    x_l = 0
    x_r = 0.75
    ax.plot([x_l, x_l, x_r, x_r], [y_br - tick, y_br, y_br, y_br - tick],
            lw=1.0, color="black", zorder=4)
    ax.text((x_l + x_r) / 2, y_br + 0.008 * span, stars, ha="center",
            va="bottom", fontsize=9.5, color="black")

    # --- stats line removed (user request); bracket + stars remain ---
    ax.set_title(cfg["title"], fontsize=10, pad=6)

    ax.set_xticks([0, 0.75])
    ax.set_xticklabels(X_LABELS, fontsize=8)
    ax.set_ylabel(cfg["ylabel"], fontsize=9)
    ax.grid(axis="y", linestyle=":", linewidth=0.5, alpha=0.4)
    ax.set_axisbelow(True)

    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT_DIR / f"beta_pooled_pair_vs_random_{DATASET_LABEL}_{metric}.{ext}",
                    dpi=300)
    plt.close(fig)
    print(f"done: {metric}  P={p_text} ({stars})  n={n_units}  dmedian={median_diff:.4g}")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    long_df = pd.read_csv(IN_LONG, sep="\t")
    stat_df = pd.read_csv(IN_STATS, sep="\t")
    for metric, cfg in METRICS.items():
        draw_panel(long_df, stat_df, metric, cfg)
    print(f"outputs -> {OUT_DIR}")


if __name__ == "__main__":
    main()