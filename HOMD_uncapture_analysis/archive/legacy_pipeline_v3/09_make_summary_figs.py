#!/usr/bin/env python3
"""
pipeline_v3 Step 09: summary figures (v2 12 equivalent, cluster units).

Figures (pipeline_v3/analysis_report/figs/):
  1. fig_funnel.png        - pipeline funnel (counts + retention)
  2. fig_novelty.png       - HOMD novelty classes + oral clade split
  3. fig_support.png       - support evidence of strict candidates
  4. fig_geography.png     - continent coverage of strict candidates
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V3 = BASE / "pipeline_v3"
OUT = PIPE_V3 / "analysis_report" / "figs"
OUT.mkdir(parents=True, exist_ok=True)

BLUE = "#4E79A7"; RED = "#E15759"; GRAY = "#BAB0AC"; TEAL = "#76B7B2"
INK = "#2B2B2B"; MUTED = "#8A8A8A"; GRID = "#E8E8E8"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10,
    "axes.edgecolor": "#BBBBBB", "axes.linewidth": 0.8,
    "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
    "text.color": INK, "axes.titlecolor": INK, "figure.dpi": 150,
})

log = lambda s: print(f"[09_figs] {s}", flush=True)


def style_ax(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#BBBBBB")
    ax.tick_params(length=3)


def fig_funnel():
    n_screened = 3396776
    n_cand = len(pd.read_csv(BASE / "pipeline_v2/01_chimera_filtered.tsv", sep="\t",
                             low_memory=False, usecols=["derep_id"]))
    n_floor = len(pd.read_csv(BASE / "pipeline_v2/04_support_floor_passed.tsv", sep="\t",
                              low_memory=False, usecols=["derep_id"]))
    s01 = json.load(open(PIPE_V3 / "summary.json"))
    s02 = json.load(open(PIPE_V3 / "02_summary.json"))
    s08 = json.load(open(PIPE_V3 / "08_supported_summary.json"))
    stages = ["ASVs screened (334 BioProjects)",
              "HOMD-uncaptured & non-chimeric",
              "Support floor (≥2 BP, ≥50 reads)",
              "ASVs placed on GG2 tree",
              "TreeCluster clusters (t=0.03,\nincl. singleton ASVs)",
              "Supported catalog\n(clean, subj>50, BP≥2, cont≥2)"]
    counts = [n_screened, n_cand, n_floor, s01["n_asv"],
              s02["n_clusters"], s08["n_supported_clean"]]

    fig, ax = plt.subplots(figsize=(9, 4.6))
    y = np.arange(len(stages))[::-1]
    colors = [GRAY, GRAY, GRAY, TEAL, TEAL, RED]
    ax.barh(y, counts, color=colors, height=0.62)
    for yi, c_, prev in zip(y, counts, [None] + counts[:-1]):
        ax.text(c_ * 1.02, yi, f"{c_:,}", va="center", fontsize=9, color=INK)
        if prev:
            ax.text(c_ * 1.02, yi + 0.34, f"{c_ / prev * 100:.1f}% retained",
                    va="center", fontsize=7.5, color=MUTED)
    ax.set_yticks(y)
    ax.set_yticklabels(stages, fontsize=9)
    ax.set_xscale("log")
    ax.set_xlim(1, max(counts) * 8)
    ax.set_xlabel("Count (log scale)")
    ax.set_title("pipeline_v3 funnel — TreeCluster t=0.03 scheme (all singleton ASVs included)",
                 fontsize=12.5, fontweight="bold", pad=10)
    style_ax(ax)
    fig.tight_layout()
    fig.savefig(OUT / "fig_funnel.png", bbox_inches="tight")
    plt.close(fig)
    log("funnel saved")


def fig_novelty():
    agg = pd.read_csv(PIPE_V3 / "02_clusters_aggregated.tsv", sep="\t", low_memory=False)
    cap = int((agg["cluster_novelty_class"] == "homd_captured").sum())
    uncap = int((agg["cluster_novelty_class"] == "homd_uncaptured").sum())
    u = agg[agg["cluster_novelty_class"] == "homd_uncaptured"]
    in_h = int(u["annotation"].map(
        lambda g: str(g).split("_")[0] in HOMD_GENUS if pd.notna(g) and str(g).strip() else False).sum())
    out_h = len(u) - in_h

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    ax = axes[0]
    bars = ax.bar(["HOMD-captured", "HOMD-uncaptured"], [cap, uncap],
                  color=[GRAY, RED], width=0.5)
    for b, v in zip(bars, [cap, uncap]):
        ax.text(b.get_x() + b.get_width() / 2, v * 1.02, f"{v:,} ({v / (cap + uncap) * 100:.1f}%)",
                ha="center", va="bottom", fontsize=9)
    ax.set_ylabel(f"Clusters ({cap + uncap:,} with ASVs)")
    ax.set_ylim(0, (cap + uncap) * 1.15)
    ax.set_title("HOMD novelty of TreeCluster clusters", fontsize=11, fontweight="bold")
    style_ax(ax)

    ax = axes[1]
    bars = ax.bar(["in HOMD genus", "outside HOMD genus"], [in_h, out_h],
                  color=[BLUE, RED], width=0.5)
    for b, v in zip(bars, [in_h, out_h]):
        ax.text(b.get_x() + b.get_width() / 2, v * 1.02, f"{v:,} ({v / len(u) * 100:.1f}%)",
                ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("HOMD-uncaptured clusters")
    ax.set_ylim(0, len(u) * 1.15)
    ax.set_title("Oral clade split of uncaptured", fontsize=11, fontweight="bold")
    style_ax(ax)
    fig.tight_layout()
    fig.savefig(OUT / "fig_novelty.png", bbox_inches="tight")
    plt.close(fig)
    log("novelty saved")


# HOMD genus set (module-level, used by fig_novelty)
HOMD_GENUS = set()
with open("/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03_full.fasta") as f:
    for line in f:
        if line.startswith(">"):
            parts = line.split()
            if len(parts) >= 3:
                HOMD_GENUS.add(parts[2].strip())


if __name__ == "__main__":
    fig_funnel()
    fig_novelty()
    log("Done!")
