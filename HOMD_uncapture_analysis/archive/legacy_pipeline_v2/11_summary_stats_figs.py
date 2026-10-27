#!/usr/bin/env python3
"""
pipeline_v2 Step 11: Summary statistics figures — the whole phylotype landscape.

Visualizes the supplementary table for ALL phylotypes (23,088), marking the 276
strict-evidence candidates:
  A. member ASVs per phylotype
  B. contamination-flagged members per phylotype (clean vs flagged)
  C. BioProjects per phylotype
  D. amplicon regions per phylotype
  E. contamination status (clean vs flagged) + clean candidates by oral clade
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

BLUE = "#4E79A7"
RED  = "#E15759"
TEAL = "#76B7B2"
GRAY = "#BAB0AC"
GOLD = "#F28E2B"
INK  = "#2B2B2B"
GRID = "#E6E6E6"

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"
OUT = PIPE_V2 / "analysis_report" / "figs"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10,
    "axes.edgecolor": "#BABABA", "axes.linewidth": 0.8,
    "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
    "text.color": INK, "axes.titlecolor": INK, "figure.dpi": 150,
})
log = lambda s: print(f"[11_stats] {s}", flush=True)

# ---- load ----
agg = pd.read_csv(PIPE_V2 / "06b_phylotypes_aggregated.tsv", sep="\t", low_memory=False)
strict_ids = set(pd.read_csv(PIPE_V2 / "07_strict_candidates.tsv", sep="\t")["phylotype_id"])
agg["is_strict"] = agg["phylotype_id"].isin(strict_ids)
log(f"all phylotypes: {len(agg):,} | strict candidates: {int(agg['is_strict'].sum()):,}")

# ============================================================
# Panel figure A-D: support landscape
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(12, 8.5))

def dist_panel(ax, values, strict_vals, bins, xlab, title, xlim=None, logx=False):
    ax.hist(values, bins=bins, color=GRAY, alpha=0.6, edgecolor="white", linewidth=0.3)
    if len(strict_vals):
        ax.hist(strict_vals, bins=bins, color=RED, alpha=0.85, edgecolor="white", linewidth=0.3,
                label=f"strict candidates (n={len(strict_vals)})")
    if logx:
        ax.set_xscale("log")
    ax.set_xlabel(xlab); ax.set_ylabel("phylotypes")
    ax.set_title(title, fontsize=11)
    ax.grid(axis="y", color=GRID, linewidth=0.7); ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    if xlim: ax.set_xlim(xlim)
    ax.legend(frameon=False, fontsize=8.5)

# A: member ASVs
ax = axes[0, 0]
asv_bins = [1, 2, 3, 5, 10, 20, 50, 100, 200]
asv_vals = agg["n_members"].clip(upper=200)
dist_panel(ax, asv_vals, agg.loc[agg["is_strict"], "n_members"].clip(upper=200),
           asv_bins, "member ASVs per phylotype", "A. Phylotype size (ASVs)", logx=True)

# B: flagged members
ax = axes[0, 1]
flag_vals = agg["n_flagged"].clip(upper=20)
dist_panel(ax, flag_vals, agg.loc[agg["is_strict"], "n_flagged"].clip(upper=20),
           np.arange(0, 21), "contamination-flagged members per phylotype",
           "B. Contamination flag composition")

# C: BioProjects
ax = axes[1, 0]
bp_vals = agg["n_bioprojects_total"].clip(upper=50)
dist_panel(ax, bp_vals, agg.loc[agg["is_strict"], "n_bioprojects_total"].clip(upper=50),
           [2, 3, 5, 10, 20, 50], "BioProjects per phylotype", "C. Cross-BioProject prevalence",
           logx=True)

# D: regions
ax = axes[1, 1]
reg_vals = agg["n_regions"].clip(upper=10)
dist_panel(ax, reg_vals, agg.loc[agg["is_strict"], "n_regions"].clip(upper=10),
           np.arange(0.5, 11.5), "amplicon regions per phylotype", "D. Cross-region detection")

fig.suptitle(f"Phylotype landscape (n={len(agg):,}) — strict candidates highlighted (n={int(agg['is_strict'].sum()):,})",
             fontsize=13, y=1.0)
plt.tight_layout()
fig.savefig(OUT / "fig_phylotype_landscape.png", bbox_inches="tight")
plt.close(fig)
log("  saved fig_phylotype_landscape.png")

# ============================================================
# Panel E: contamination status + clean-candidate novelty split
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(11, 4))

# E1: clean vs flagged phylotypes
ax = axes[0]
cats = ["clean\n(n_flagged=0)", "flagged\n(n_flagged≥1)"]
vals = [(agg["n_flagged"] == 0).sum(), (agg["n_flagged"] > 0).sum()]
bars = ax.bar(cats, vals, color=[TEAL, GOLD], edgecolor="white", linewidth=0.6, width=0.55)
for b, v in zip(bars, vals):
    ax.text(b.get_x() + b.get_width()/2, v + 200, f"{v:,}", ha="center", fontsize=10, fontweight="bold")
ax.set_ylabel("phylotypes")
ax.set_title("E1. Contamination status of all phylotypes", fontsize=11)
ax.grid(axis="y", color=GRID, linewidth=0.7); ax.set_axisbelow(True)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax.set_ylim(0, max(vals) * 1.12)

# E2: clean candidate phylotypes by oral clade (all are homd_uncaptured;
# the informative split is whether the placement-clade genus is in HOMD)
ax = axes[1]
cc = agg[agg["clean_phylotype"] & agg["candidate_novel"]]
nov_cats = ["outside\nHOMD genus", "within\nHOMD genus"]
nov_vals = [int((~cc["oral_clade"]).sum()),
            int(cc["oral_clade"].sum())]
bars = ax.bar(nov_cats, nov_vals, color=[RED, BLUE], edgecolor="white", linewidth=0.6, width=0.55)
for b, v in zip(bars, nov_vals):
    ax.text(b.get_x() + b.get_width()/2, v + 40, f"{v:,}", ha="center", fontsize=10, fontweight="bold")
ax.set_ylabel("clean candidate phylotypes")
ax.set_title(f"E2. Clean candidates by oral clade (n={len(cc):,})", fontsize=11)
ax.grid(axis="y", color=GRID, linewidth=0.7); ax.set_axisbelow(True)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax.set_ylim(0, max(nov_vals) * 1.15)

fig.tight_layout()
fig.savefig(OUT / "fig_phylotype_contamination.png", bbox_inches="tight")
plt.close(fig)
log("  saved fig_phylotype_contamination.png")

log("\nDone!")
