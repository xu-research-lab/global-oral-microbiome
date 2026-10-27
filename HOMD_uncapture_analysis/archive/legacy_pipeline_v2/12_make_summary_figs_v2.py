#!/usr/bin/env python3
"""
pipeline_v2 Step 12: summary figures for the NEW scheme (GG2 2024.09
closed-reference annotation + edge x label phylotypes + STRICT candidates).

Figures (analysis_report/figs_v2/):
  1. fig_v2_funnel.png        — pipeline funnel (counts + retention)
  2. fig_v2_novelty.png       — HOMD novelty classes + oral clade split
  3. fig_v2_strict_support.png— support evidence of STRICT candidates
  4. fig_v2_strict_geography.png — continent coverage of STRICT candidates

Y-axis labels use species (or genus) names instead of cryptic phylotype IDs.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"
OUT = PIPE_V2 / "analysis_report" / "figs_v2"
OUT.mkdir(parents=True, exist_ok=True)

BLUE = "#4E79A7"; RED = "#E15759"; GRAY = "#BAB0AC"; TEAL = "#76B7B2"
INK = "#2B2B2B"; MUTED = "#8A8A8A"; GRID = "#E8E8E8"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10,
    "axes.edgecolor": "#BBBBBB", "axes.linewidth": 0.8,
    "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
    "text.color": INK, "axes.titlecolor": INK, "figure.dpi": 150,
})

log = lambda s: print(f"[12_figs] {s}", flush=True)


def style_ax(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#BBBBBB")
    ax.tick_params(length=3)


def load_strict():
    c = pd.read_csv(PIPE_V2 / "11_headline_candidates.tsv", sep="\t", low_memory=False)
    for col in ["rep_identity", "max_subject_prevalence", "n_bioprojects_total",
                "n_regions", "edge_num"]:
        c[col] = pd.to_numeric(c[col], errors="coerce")
    c["oral_clade"] = c["oral_clade"].astype(str).str.lower().isin(["true", "1"])
    # readable row label: species -> genus -> "Edge N (unnamed)"
    c["label"] = c["rep_species"].fillna(c["rep_genus"])
    unnamed = c["label"].isna()
    c.loc[unnamed, "label"] = "Edge " + c.loc[unnamed, "edge_num"].astype(int).astype(str) + " (unnamed)"
    return c


def fig_funnel():
    n_screened = 3396776
    n_cand = len(pd.read_csv(PIPE_V2 / "01_chimera_filtered.tsv", sep="\t",
                             low_memory=False, usecols=["derep_id"]))
    n_floor = len(pd.read_csv(PIPE_V2 / "04_support_floor_passed.tsv", sep="\t",
                              low_memory=False, usecols=["derep_id"]))
    s10 = json.load(open(PIPE_V2 / "10_summary.json"))
    s11 = json.load(open(PIPE_V2 / "11_summary.json"))
    stages = ["ASVs screened (334 BioProjects)",
              "HOMD-uncaptured & non-chimeric",
              "Support floor (≥2 BP, ≥50 reads)",
              "GG2 2024.09 annotated ASVs",
              "Phylotypes (edge × label)",
              "Supported catalog\n(clean)",
              "Headline\n(GG2≥99% & HOMD-missing)"]
    counts = [n_screened, n_cand, n_floor, s10["n_asvs_annotated"],
              s10["n_phylotypes"], s11["n_supported_clean"],
              s11["headline"]]

    fig, ax = plt.subplots(figsize=(9, 4.6))
    y = np.arange(len(stages))[::-1]
    colors = [GRAY, GRAY, GRAY, TEAL, TEAL, BLUE, RED]
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
    ax.set_title("pipeline_v2 funnel — GG2 2024.09 closed-reference scheme",
                 fontsize=12.5, fontweight="bold", pad=10)
    style_ax(ax)
    fig.tight_layout()
    fig.savefig(OUT / "fig_v2_funnel.png", bbox_inches="tight")
    plt.close(fig)
    log("funnel saved")


def fig_novelty():
    phy = pd.read_csv(PIPE_V2 / "11_phylotypes_homd_novelty.tsv", sep="\t", low_memory=False)
    phy["oral_clade"] = phy["oral_clade"].astype(str).str.lower().isin(["true", "1"])
    cap = int((phy["novelty_class"] == "homd_captured").sum())
    uncap = int((phy["novelty_class"] == "homd_uncaptured").sum())
    u = phy[phy["novelty_class"] == "homd_uncaptured"]
    in_h = int(u["oral_clade"].sum()); out_h = len(u) - in_h

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    ax = axes[0]
    bars = ax.bar(["HOMD-captured", "HOMD-uncaptured"], [cap, uncap],
                  color=[GRAY, RED], width=0.5)
    for b, v in zip(bars, [cap, uncap]):
        ax.text(b.get_x() + b.get_width() / 2, v * 1.02, f"{v:,} ({v / (cap + uncap) * 100:.1f}%)",
                ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("Phylotypes (26,316)")
    ax.set_ylim(0, (cap + uncap) * 1.15)
    ax.set_title("HOMD novelty of new phylotypes", fontsize=11, fontweight="bold")
    style_ax(ax)
    ax.grid(axis="y", color=GRID, linewidth=0.7)

    ax = axes[1]
    bars = ax.bar(["Within HOMD genus", "Outside HOMD genera"], [in_h, out_h],
                  color=[BLUE, RED], width=0.5)
    for b, v in zip(bars, [in_h, out_h]):
        ax.text(b.get_x() + b.get_width() / 2, v * 1.02, f"{v:,} ({v / len(u) * 100:.1f}%)",
                ha="center", va="bottom", fontsize=9)
    ax.set_ylabel(f"HOMD-uncaptured ({len(u):,})")
    ax.set_ylim(0, len(u) * 1.15)
    ax.set_title("Uncaptured phylotypes by oral clade", fontsize=11, fontweight="bold")
    style_ax(ax)
    ax.grid(axis="y", color=GRID, linewidth=0.7)

    fig.tight_layout()
    fig.savefig(OUT / "fig_v2_novelty.png", bbox_inches="tight")
    plt.close(fig)
    log("novelty saved")


def fig_strict_support():
    c = load_strict().sort_values("max_subject_prevalence", ascending=False)
    y = np.arange(len(c))
    labels = c["label"].astype(str) + "  (" + c["rep_identity"].round(1).astype(str) + "% GG2)"

    fig, axes = plt.subplots(1, 3, figsize=(12.5, max(5.0, len(c) * 0.34)), sharey=True)
    for ax, col, title, color, logscale in [
        (axes[0], "max_subject_prevalence", "Subjects (max)", BLUE, True),
        (axes[1], "n_bioprojects_total", "BioProjects", TEAL, False),
        (axes[2], "n_regions", "Amplicon regions", MUTED, False)]:
        vals = c[col].astype(float)
        ax.barh(y, vals, color=color, height=0.62)
        for yi, v in zip(y, vals):
            ax.text(v + (max(vals) * 0.02 if not logscale else v * 0.12), yi,
                    f"{v:g}", va="center", fontsize=8.5)
        if logscale:
            ax.set_xscale("log")
            ax.set_xlim(1, max(vals) * 20)
        else:
            ax.set_xlim(0, max(vals) * 1.22)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_title(title, fontsize=11, fontweight="bold")
        style_ax(ax)
        ax.grid(axis="x", color=GRID, linewidth=0.7)
    fig.suptitle(f"{len(c)} HEADLINE candidates — "
                 "support evidence (sorted by subjects)",
                 fontsize=12.5, fontweight="bold", y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(OUT / "fig_v2_strict_support.png", bbox_inches="tight")
    plt.close(fig)
    log("strict support saved")


def fig_strict_geography():
    c = load_strict()
    mp = pd.read_csv(PIPE_V2 / "10_asv_phylotype_map.tsv", sep="\t", low_memory=False,
                     usecols=["derep_id", "phylotype_id"])
    s04 = pd.read_csv(PIPE_V2 / "04_support_floor_passed.tsv", sep="\t", low_memory=False,
                      usecols=["derep_id", "continents"])

    def parse_set(val):
        if pd.isna(val) or val == "":
            return []
        s = str(val).strip("{}")
        for sep in ["; ", ";", ", ", ","]:
            if sep in s:
                return [x.strip().strip("'\"") for x in s.split(sep)]
        return [s.strip().strip("'\"")]

    cont_rows = []
    for pid in c["phylotype_id"]:
        members = mp.loc[mp["phylotype_id"] == pid, "derep_id"]
        all_c = set()
        for v in s04.loc[s04["derep_id"].isin(members), "continents"].dropna():
            all_c.update(parse_set(v))
        cont_rows.append(all_c)
    c = c.copy()
    c["continents"] = cont_rows
    c = c.sort_values("max_subject_prevalence", ascending=False)

    continent_colors = {"Asia": "#4E79A7", "Europe": "#59A14F", "North America": "#F28E2B",
                        "South America": "#E15759", "Africa": "#B07AA1", "Oceania": "#76B7B2"}
    fig, ax = plt.subplots(figsize=(8.8, max(4.8, len(c) * 0.34)))
    y = np.arange(len(c))
    left = np.zeros(len(c))
    for cont, color in continent_colors.items():
        vals = c["continents"].apply(lambda s: 1 if cont in s else 0).to_numpy()
        if vals.sum() == 0:
            continue
        ax.barh(y, vals, left=left, color=color, height=0.62, label=cont)
        left += vals
    ax.set_yticks(y)
    ax.set_yticklabels(c["label"], fontsize=8.5)
    ax.set_xlim(0, 6.3)
    ax.set_xticks(range(0, 7))
    ax.set_xlabel("Number of continents (presence)")
    ax.set_title(f"{len(c)} HEADLINE candidates — continent coverage", fontsize=12.5,
                 fontweight="bold", pad=10)
    ax.legend(ncol=3, fontsize=8.5, frameon=False, loc="lower right")
    style_ax(ax)
    fig.tight_layout()
    fig.savefig(OUT / "fig_v2_strict_geography.png", bbox_inches="tight")
    plt.close(fig)
    log("geography saved")


if __name__ == "__main__":
    fig_funnel()
    fig_novelty()
    fig_strict_support()
    fig_strict_geography()
    log("Done!")
