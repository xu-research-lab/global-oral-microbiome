#!/usr/bin/env python3
"""
pipeline_v2 Step 11: Authenticity Tiering Figures.

Figures:
  1. fig_tiering_candidates.png  — tier distribution among candidate novel phylotypes
  2. fig_tiering_final.png       — final-novel set by novelty class (Tier1+2)
  3. fig_tiering_sensitivity.png — where the strict-pipeline 56 land in V2 tiers
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

# ---- Style: match pipeline_v2 (Tableau-10 colorblind-safe) ----
BLUE  = "#4E79A7"   # known_not_in_homd
RED   = "#E15759"   # novel_all_refs / Tier-1
TEAL  = "#76B7B2"   # Tier-2
GOLD  = "#F28E2B"   # Tier-3
GRAY  = "#BAB0AC"   # below threshold
INK   = "#2B2B2B"
GRID  = "#E6E6E6"

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"
OUT = PIPE_V2 / "analysis_report" / "figs"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.edgecolor": "#BABABA",
    "axes.linewidth": 0.8,
    "axes.labelcolor": INK,
    "xtick.color": INK,
    "ytick.color": INK,
    "text.color": INK,
    "axes.titlecolor": INK,
    "figure.dpi": 150,
})

log = lambda s: print(f"[11_tier] {s}", flush=True)

# ============================================================
# 1. Load tiered table
# ============================================================
log("Loading 06b_phylotypes_tiered.tsv...")
t = pd.read_csv(PIPE_V2 / "06b_phylotypes_tiered.tsv", sep="\t", low_memory=False)
cn = t[t["candidate_novel"]].copy()
fn = t[t["final_novel"]].copy()

log(f"  candidates: {len(cn):,} | final_novel: {len(fn):,}")

# ============================================================
# Fig 1: Tier distribution among candidates (stacked horizontal)
# ============================================================
log("Fig 1: candidate tier distribution...")

tier_order = [("Tier 1\nresidence-validated", 1, RED),
              ("Tier 2\nreproducibility-validated", 2, TEAL),
              ("Tier 3\nconcern", 3, GOLD),
              ("below threshold", 0, GRAY)]
tier_vals = [int((cn["tier"] == k).sum()) for _, k, _ in tier_order]

fig, ax = plt.subplots(figsize=(7.5, 4.2))
y = 0
for (lab, k, col), v in zip(tier_order, tier_vals):
    ax.barh(y, v, left=0, color=col, edgecolor="white", linewidth=0.6, height=0.55,
            alpha=0.92)
    ax.text(v + 60, y, f"{v:,}", va="center", fontsize=10, fontweight="bold")
    y += 1

ax.set_yticks(np.arange(len(tier_order)))
ax.set_yticklabels([l[0] for l in tier_order], fontsize=9)
ax.set_xlabel("candidate novel phylotypes")
ax.set_title(f"Authenticity tiers among candidate novel phylotypes (n={len(cn):,})",
             fontsize=12, pad=12)
ax.grid(axis="x", color=GRID, linewidth=0.7)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_xlim(0, max(tier_vals) * 1.15)
plt.tight_layout()
fig.savefig(OUT / "fig_tiering_candidates.png", bbox_inches="tight")
plt.close(fig)
log("  saved fig_tiering_candidates.png")

# ============================================================
# Fig 2: Final-novel set by novelty class + tier
# ============================================================
log("Fig 2: final-novel composition...")

# stack by tier, within novelty class
novel_classes = ["known_not_in_homd", "novel_all_refs"]
n_cols = {1: RED, 2: TEAL}
fig, ax = plt.subplots(figsize=(7.5, 4.6))
x = np.arange(len(novel_classes))
for tier, col in n_cols.items():
    vals = [int(((fn["novelty_class"] == c) & (fn["tier"] == tier)).sum())
            for c in novel_classes]
    if tier == 1:
        bottom = [0, 0]
    else:
        bottom = [int(((fn["novelty_class"] == c) & (fn["tier"] == 1)).sum())
                  for c in novel_classes]
    ax.bar(x, vals, bottom=bottom, width=0.5, color=col, edgecolor="white",
           linewidth=0.5, label=f"Tier {tier}")

for xi, c in enumerate(novel_classes):
    tot = int((fn["novelty_class"] == c).sum())
    ax.text(xi, tot + 20, f"{tot:,}", ha="center", fontsize=11, fontweight="bold")

ax.set_xticks(x)
ax.set_xticklabels(["known not in HOMD", "novel to all refs"], fontsize=9.5)
ax.set_ylabel("final-novel phylotypes")
ax.set_title(f"Final novel oral phylotypes (Tier 1+2, n={len(fn):,})",
             fontsize=12, pad=12)
ax.grid(axis="y", color=GRID, linewidth=0.7)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.legend(frameon=False, fontsize=9)
ax.set_ylim(0, max([int((fn["novelty_class"] == c).sum()) for c in novel_classes]) * 1.12)
plt.tight_layout()
fig.savefig(OUT / "fig_tiering_final.png", bbox_inches="tight")
plt.close(fig)
log("  saved fig_tiering_final.png")

# ============================================================
# Fig 3: Strict-pipeline 56 sensitivity
# ============================================================
log("Fig 3: strict-56 sensitivity...")

sens = pd.read_csv(PIPE_V2 / "06b_strict56_sensitivity.tsv", sep="\t")
cats = [("Tier 1", 1, RED), ("Tier 2", 2, TEAL), ("Tier 3", 3, GOLD),
        ("below threshold", 0, GRAY), ("not in V2 space", -1, "#8C8C8C")]
vals = [int((sens["tier"] == k).sum()) for _, k, _ in cats]
lab = [c[0] for c in cats]

fig, ax = plt.subplots(figsize=(7.5, 3.8))
bars = ax.bar(lab, vals, color=[c[2] for c in cats], edgecolor="white",
              linewidth=0.6, width=0.6)
for b, v in zip(bars, vals):
    ax.text(b.get_x() + b.get_width() / 2, v + 0.15, f"{v}", ha="center",
            fontsize=10, fontweight="bold")
ax.set_ylabel("strict-pipeline candidates")
ax.set_title("Fate of the 56 strict-pipeline candidates in the V2 tiering",
             fontsize=12, pad=12)
ax.set_ylim(0, max(vals) * 1.15)
ax.grid(axis="y", color=GRID, linewidth=0.7)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(axis="x", labelsize=8.5)
plt.tight_layout()
fig.savefig(OUT / "fig_tiering_sensitivity.png", bbox_inches="tight")
plt.close(fig)
log("  saved fig_tiering_sensitivity.png")

log("\nDone!")
