#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_mixed_bp_unifrac_violin.py
(still this name for path stability — figure kept under
10_mixed_bp_unifrac_violin/)

Re-render of
results/single_pair/03_beta/pcoa_permanova/plots/paired_unifrac_hc_vs_disease_mixed_bp.png
(original: scripts/python/oral_gut_pair/report_figures.py ->
plot_mixed_bp_paired_distance). After an aesthetic iteration with violin
plots the user preferred the ORIGINAL boxplot + jitter drawing style, so
the group rendering is now the original code again.

Data pipeline is byte-identical to the original function:
  * paired oral-gut weighted UniFrac distances, dataset_label="filter",
    complete pairs only
  * pair status from pair_metrics_all_cohorts.tsv (single_pair set) via
    report_figures._pair_status_binary (HC / disease / Unknown)
  * mixed BioProjects = classify_bp_status_group(...) == "HC+disease"
  * Mann-Whitney U (disease vs HC) within each BP, BH-corrected across BPs,
    stars via report_figures._sig_stars
  * BP order = ascending (median_disease - median_hc)
  * same seed (2026) for jitter

Change: each HC / disease group is drawn as filled boxplot (alpha=0.8) +
jitter scatter, exactly like the original function; palette unchanged
(HC #3182bd / disease #de2d26, same as the 08 figure). The only deviations
from the original are user-requested:
  * BP order: significant BPs first (ascending q), the rest in the
    original ascending-delta order
  * legend inside the axes, top right

NOTE: x ticks/labels MUST be set after the boxes are drawn — ax.boxplot()
registers each box position as a tick (this bit fig 10 once already).

Inputs:
  results/single_pair/03_beta/distance_matrix/paired_beta_distances_all_pairs.tsv
  results/cross_comparison/01_merged_metrics/pair_metrics_all_cohorts.tsv
Outputs: results/new/figures/10_mixed_bp_unifrac_violin/
  paired_unifrac_hc_vs_disease_mixed_bp_violin.{png,pdf}
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path("/home/hexincheng/projects/hxc/oral_project/oral_gut_pair")
sys.path.insert(0, str(ROOT / "scripts" / "python"))

from oral_gut_pair import report_figures as rf  # noqa: E402

IN_DIST = ROOT / "results/single_pair/03_beta/distance_matrix/paired_beta_distances_all_pairs.tsv"
IN_META = ROOT / "results/cross_comparison/01_merged_metrics/pair_metrics_all_cohorts.tsv"
OUT_DIR = ROOT / "results/new/figures/10_mixed_bp_unifrac_violin"

METRIC = "weighted_unifrac"
# same palette as the 08 delta-alpha figure (HC blue / disease red)
PAIR_STATUS_COLORS = {"HC": "#3182bd", "disease": "#de2d26", "Unknown": "#969696"}

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

    dist = pd.read_csv(IN_DIST, sep="\t")
    dist = dist.loc[
        dist["metric"].astype(str).eq(METRIC)
        & dist["dataset_label"].astype(str).eq("filter")
        & dist["has_complete_oral_gut_pair"].astype(str).isin(["True", "true", "1"])
    ].copy()
    meta = pd.read_csv(IN_META, sep="\t")
    meta = meta.loc[
        meta["analysis_set"].astype(str).eq("single_pair"),
        ["pair_id", "bioproject", "status_cleaned"],
    ].drop_duplicates("pair_id")
    df = dist.merge(meta, on=["pair_id", "bioproject"], how="inner")
    df["pair_status"] = df["status_cleaned"].map(rf._pair_status_binary)
    df["distance"] = pd.to_numeric(df["distance"], errors="coerce")
    df = df.dropna(subset=["distance"])

    bp_group = (
        df.groupby("bioproject")["status_cleaned"]
        .apply(lambda s: rf.classify_bp_status_group(s.tolist())[1])
    )
    mixed = bp_group[bp_group.eq("HC+disease")].index.tolist()
    df = df.loc[
        df["bioproject"].isin(mixed) & df["pair_status"].isin(["HC", "disease"])
    ].copy()

    labels = {}
    for bp, bp_sub in df.groupby("bioproject"):
        disease = rf._bp_disease_label(bp_sub)
        labels[bp] = f"{bp}\n{disease}" if disease else str(bp)

    sub = df.loc[df["metric"].eq(METRIC)]
    test_rows = []
    for bp in mixed:
        bp_sub = sub.loc[sub["bioproject"].eq(bp)]
        hc = bp_sub.loc[bp_sub["pair_status"].eq("HC"), "distance"].to_numpy(dtype=float)
        dz = bp_sub.loc[bp_sub["pair_status"].eq("disease"), "distance"].to_numpy(dtype=float)
        if len(hc) == 0 or len(dz) == 0:
            continue
        try:
            p = float(mannwhitneyu(dz, hc, alternative="two-sided").pvalue)
        except ValueError:
            p = float("nan")
        test_rows.append(
            {
                "bioproject": bp,
                "n_hc": int(len(hc)),
                "n_disease": int(len(dz)),
                "median_hc": float(np.median(hc)),
                "median_disease": float(np.median(dz)),
                "mwu_p": p,
            }
        )
    tests = pd.DataFrame(test_rows)
    tests["q_BH"] = rf._bh_fdr(tests["mwu_p"].tolist())
    tests["sig_stars"] = tests["q_BH"].map(rf._sig_stars).fillna("")
    tests["delta"] = tests["median_disease"] - tests["median_hc"]
    # user ordering: significant BPs first (ascending q), the rest in the
    # original ascending-delta order
    sig = tests.loc[tests["sig_stars"].ne("")].sort_values("q_BH")
    nonsig = tests.loc[tests["sig_stars"].eq("")].sort_values("delta")
    bp_order = sig["bioproject"].tolist() + nonsig["bioproject"].tolist()

    # canvas shrunk (user request: 缩小) — ~0.62x linear scale of the
    # previous 11.5x4.8; fonts dropped one step each (title 10, ticks 6.5)
    fig, ax = plt.subplots(figsize=(max(7.0, len(bp_order) * 0.78), 3.0))
    rng = np.random.default_rng(2026)
    metric_tests = tests.set_index("bioproject")
    positions, data, colors, ns, tick_pos, tick_labels = [], [], [], [], [], []
    star_spans = []
    for i, bp in enumerate(bp_order):
        bp_sub = sub.loc[sub["bioproject"].eq(bp)]
        xs = []
        for st in ("HC", "disease"):
            vals = bp_sub.loc[bp_sub["pair_status"].eq(st), "distance"].to_numpy(dtype=float)
            x = i + (-0.18 if st == "HC" else 0.18)
            positions.append(x)
            data.append(vals)
            colors.append(PAIR_STATUS_COLORS[st])
            ns.append(len(vals))
            xs.append(x)
        tick_pos.append(i)
        tick_labels.append(labels.get(bp, bp))
        stars = str(metric_tests.loc[bp, "sig_stars"] or "")
        if stars and len(xs) == 2:
            star_spans.append((xs[0], xs[1], stars, np.concatenate(data[-2:])))

    bp_artist = ax.boxplot(
        data,
        positions=positions,
        widths=0.28,
        showfliers=False,
        patch_artist=True,
        medianprops=dict(color="#222", linewidth=1.2),
    )
    for patch, color in zip(bp_artist["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)
        patch.set_edgecolor("#333")
    for pos, vals, color in zip(positions, data, colors):
        if len(vals) == 0:
            continue
        ax.scatter(
            np.full(len(vals), pos) + rng.uniform(-0.05, 0.05, len(vals)),
            vals,
            s=8,
            alpha=0.28,
            color=color,
            linewidths=0,
            zorder=3,
        )
    ymax = max(float(np.nanmax(v)) for v in data if len(v))
    ymin = min(float(np.nanmin(v)) for v in data if len(v))
    span = ymax - ymin if ymax > ymin else 0.1
    ax.set_ylim(ymin - 0.04 * span, ymax + 0.22 * span)
    for x0, x1, stars, vals in star_spans:
        y_br = float(np.nanmax(vals)) + 0.06 * span
        ax.plot(
            [x0, x0, x1, x1],
            [y_br - 0.02 * span, y_br, y_br, y_br - 0.02 * span],
            color="#8b0000",
            linewidth=1.0,
            zorder=4,
        )
        ax.text(
            (x0 + x1) / 2,
            y_br + 0.01 * span,
            stars,
            ha="center",
            va="bottom",
            fontsize=8,
            fontweight="bold",
            color="#8b0000",
        )
    y_top = ax.get_ylim()[1]
    for pos, n in zip(positions, ns):
        ax.text(pos, y_top, f"n={n}", ha="center", va="bottom",
                fontsize=5.5, color="#444", clip_on=False)
    ax.set_ylabel("Weighted UniFrac\noral–gut distance", fontsize=9)
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_labels, rotation=0, fontsize=6)
    ax.set_xlim(-0.6, len(bp_order) - 0.4)
    ax.grid(axis="y", linestyle=":", linewidth=0.5, alpha=0.4)
    ax.set_axisbelow(True)

    from matplotlib.patches import Patch

    handles = [
        Patch(facecolor=PAIR_STATUS_COLORS["HC"], edgecolor="#333", label="HC"),
        Patch(facecolor=PAIR_STATUS_COLORS["disease"], edgecolor="#333", label="disease"),
    ]
    ax.legend(handles=handles, loc="upper right", ncol=2, frameon=False,
              fontsize=8, columnspacing=1.0, handletextpad=0.5)
    fig.suptitle(
        "Paired oral–gut weighted UniFrac within mixed BioProjects\n"
        "(Mann–Whitney within BP; BH across BioProjects; * q<0.05)",
        # y=1.08 was tuned for the ORIGINAL layout whose fig.legend filled
        # the top band; with the legend moved inside the axes that band is
        # empty, so the title sits just above the n= row instead (y in
        # figure coords; va="top" is the default)
        y=0.99,
        fontsize=10,
    )
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(
            OUT_DIR / f"paired_unifrac_hc_vs_disease_mixed_bp_violin.{ext}",
            dpi=300,
            bbox_inches="tight",
        )
    plt.close(fig)
    print(tests[["bioproject", "n_hc", "n_disease", "mwu_p", "q_BH", "sig_stars"]]
          .to_string(index=False))
    print(f"outputs -> {OUT_DIR}")


if __name__ == "__main__":
    main()