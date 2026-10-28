#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_site_r2_hc_vs_disease.py

Standalone re-render of ONE panel from
results/single_pair/03_beta/pcoa_permanova/plots/pcoa_hc_vs_disease_weighted_unifrac.png
(original: scripts/python/oral_gut_pair/report_figures.py ->
make_pcoa_status_figures, bottom-left subplot): the Weighted UniFrac
PERMANOVA site-R² boxplot comparing HC-only vs Disease-only BioProjects.

Everything else in the original composite (PCoA scatter panels, unweighted
panel, legend) is intentionally NOT drawn.

Data source is the table the original function itself wrote:
  results/single_pair/03_beta/pcoa_permanova/pcoa_site_r2_hc_vs_disease.tsv
(per-BP paired adonis2 D ~ site, R² of the "site" term), so the medians,
MWU p and point positions are identical to the composite figure.

Style updates (user request, shares the same palette and font-size set as
the 10 mixed-BP violin script):
  * unified palette with the 08 delta-alpha figure: HC blue #3182bd /
    disease red #de2d26; boxes translucent-filled, dark bold median line
  * scatter: small white-edged points, reduced jitter
  * disease labels in one right-aligned column, staggered vertically, thin
    grey leader lines back to each point
  * MWU result as a significance bracket between the groups (not in the
    title); title only "Weighted UniFrac"; medians + BP counts under the
    x labels
  * y axis 0-0.8

Inputs : results/single_pair/03_beta/pcoa_permanova/pcoa_site_r2_hc_vs_disease.tsv
Outputs: results/new/figures/09_site_r2_hc_vs_disease/
  pcoa_site_r2_weighted_unifrac_hc_vs_disease.{png,pdf}
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path("/home/hexincheng/projects/hxc/oral_project/oral_gut_pair")
sys_path = ROOT / "scripts" / "python"
if str(sys_path) not in __import__("sys").path:
    __import__("sys").path.insert(0, str(sys_path))

from oral_gut_pair import report_figures as rf  # noqa: E402

IN_R2 = ROOT / "results/single_pair/03_beta/pcoa_permanova/pcoa_site_r2_hc_vs_disease.tsv"
OUT_DIR = ROOT / "results/new/figures/09_site_r2_hc_vs_disease"
METRIC = "weighted_unifrac"
METRIC_TITLE = "Weighted UniFrac"
# same palette as the 08 delta-alpha figure (HC blue / disease red)
HC_COLOR = "#3182bd"
DISEASE_COLOR = "#de2d26"

# shared font-size set (shrunk with the canvas, user request 2026-09-30:
# shrink both figures; ticks 8 / label 9 / title 10 matched to figs 05/06)
FONT_TICK = 8
FONT_LABEL = 9
FONT_TITLE = 10
FONT_ANNO = 7.5

plt.rcParams.update(
    {
        "pdf.fonttype": 42,
        "font.size": 9,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
    }
)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    perm = pd.read_csv(IN_R2, sep="\t")
    sub = perm.loc[perm["metric"].astype(str).eq(METRIC)].copy()

    hc = sub.loc[sub["status_group"].eq("HC"), "r2_term"].astype(float).to_numpy()
    dz = sub.loc[sub["status_group"].eq("disease_study"), "r2_term"].astype(float).to_numpy()
    p = float(mannwhitneyu(dz, hc, alternative="two-sided").pvalue)

    rng = np.random.default_rng(2026)
    # canvas shrunk (user request): 7.0x4.2 -> 4.7x2.9
    fig, ax = plt.subplots(figsize=(4.7, 2.9))

    data = []
    for i, group in enumerate(["HC", "disease_study"]):
        vals = sub.loc[sub["status_group"].eq(group), "r2_term"].astype(float).to_numpy()
        data.append(vals)
        jitter = rng.uniform(-0.05, 0.05, len(vals))
        color = HC_COLOR if group == "HC" else DISEASE_COLOR
        ax.scatter(
            np.full(len(vals), i + 1) + jitter,
            vals,
            s=14,
            color=color,
            alpha=0.85,
            edgecolors="white",
            linewidths=0.4,
            zorder=3,
        )

    bplot = ax.boxplot(
        data,
        positions=[1, 2],
        widths=0.45,
        showfliers=False,
        patch_artist=True,
        medianprops=dict(color="#1a1a1a", linewidth=1.8),
        whiskerprops=dict(color="#555", linewidth=0.9),
        capprops=dict(color="#555", linewidth=0.9),
    )
    # box fill identical to the 10 mixed-BP figure: solid group colour at
    # alpha=0.8 with #333 edges
    for patch, color in zip(bplot["boxes"], [HC_COLOR, DISEASE_COLOR]):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)
        patch.set_edgecolor("#333")
        patch.set_linewidth(0.9)

    # disease labels: one aligned column at the right edge of the panel,
    # staggered vertically (MIN_GAP), thin grey leader lines back to the
    # corresponding points
    labs = sorted(
        [(float(r["r2_term"]), rf._short_disease_phrase(r["status_values"]))
         for _, r in sub.loc[sub["status_group"].eq("disease_study")].iterrows()]
    )
    X_COL = 2.55   # label column close to the right edge of the disease bar
    MIN_GAP = 0.05
    ys = []
    for y, _t in labs:
        if ys and y < ys[-1] + MIN_GAP:
            y = ys[-1] + MIN_GAP
        ys.append(y)
    for (r2, text), y_col in zip(labs, ys):
        ax.annotate(
            text,
            xy=(2.08, r2),
            xytext=(X_COL, y_col),
            fontsize=7,
            color="#1a1a1a",
            va="center",
            ha="left",
            annotation_clip=False,
            arrowprops=dict(arrowstyle="-", color="#999", linewidth=0.6,
                            shrinkA=0, shrinkB=2),
        )

    ax.set_xticks([1, 2])
    ax.set_xticklabels(
        [
            f"HC-only (n={len(hc)} BP)",
            f"Disease-only (n={len(dz)} BP)",
        ],
        fontsize=FONT_TICK,
    )
    ax.set_ylabel("PERMANOVA site R²", fontsize=FONT_LABEL)
    ax.set_ylim(0, 0.8)
    ax.set_xlim(0.25, 4.15)   # provisional; right edge tightened below

    # significance bracket between the two groups (replaces the in-title p)
    y_br = 0.73
    ax.plot([1, 1, 2, 2], [y_br - 0.02, y_br, y_br, y_br - 0.02],
            lw=0.8, color="black", zorder=4)
    ax.text(1.5, y_br + 0.004, f"MWU p = {p:.4f}", ha="center", va="bottom",
            fontsize=FONT_ANNO, color="black")
    ax.set_title(METRIC_TITLE, fontsize=FONT_TITLE, pad=6)
    ax.grid(axis="y", linestyle=":", linewidth=0.5, alpha=0.4)
    ax.set_axisbelow(True)

    fig.tight_layout()
    # right border hugs the labels: fit xlim to the longest disease label.
    # MUST run AFTER tight_layout — it resizes the axes, which changes how
    # many data units a fixed-pixel text spans, and the longest label
    # ("COVID conval. +1") leaked past the border when measured before.
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    anns = [t for t in ax.texts if t.get_text()]
    x_right = max(inv.transform((b[0] + b[2], 0.0))[0]
                  for b in (a.get_window_extent(r).bounds for a in anns))
    ax.set_xlim(0.25, x_right + 0.04)   # was 4.15 (large empty right margin)
    for ext in ("png", "pdf"):
        fig.savefig(
            OUT_DIR / f"pcoa_site_r2_weighted_unifrac_hc_vs_disease.{ext}",
            dpi=300,
            bbox_inches="tight",
        )
    plt.close(fig)
    print(
        f"HC n={len(hc)} BP (median {np.median(hc):.3f}), "
        f"disease n={len(dz)} BP (median {np.median(dz):.3f}), MWU p={p:.5g}"
    )
    print(f"outputs -> {OUT_DIR}")


if __name__ == "__main__":
    main()