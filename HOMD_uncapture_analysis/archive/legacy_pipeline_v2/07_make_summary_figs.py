#!/usr/bin/env python3
"""
pipeline_v2 Step 07: Summary & Novelty Figures

Figures:
  1. fig_funnel_pipeline.png      — pipeline_v2 flow funnel
  2. fig_novelty_classification.png — novelty classes (all + candidate)
  3. fig_candidate_support.png    — support metrics of candidate novel phylotypes
  4. fig_candidate_geography.png  — continent profile of candidate novel phylotypes

Color: Tableau-10 (colorblind-safe) categorical palette.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from pathlib import Path

# ---- Style: colorblind-safe palette (Tableau 10) ----
BLUE   = "#4E79A7"   # within HOMD genus (oral_clade = True) / captured
RED    = "#E15759"   # outside HOMD genus / uncaptured
GRAY   = "#BAB0AC"   # homd_captured / neutral
TEAL   = "#76B7B2"   # funnel accent
INK    = "#2B2B2B"
MUTED  = "#8A8A8A"
GRID   = "#E6E6E6"

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

log = lambda s: print(f"[07_figs] {s}", flush=True)

# ============================================================
# 1. Load data
# ============================================================
log("Loading novelty-classified phylotypes...")
df = pd.read_csv(PIPE_V2 / "06_phylotypes_novelty_classified.tsv", sep="\t", low_memory=False)
cn = df[df["candidate_novel"] == True].copy()

# numeric coercion (prevalence is subject-level, 2026-08-13)
for c in ["n_asvs", "max_n_bioprojects", "max_subject_prevalence", "total_reads_sum", "n_continents"]:
    df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)
    cn[c] = pd.to_numeric(cn[c], errors="coerce").fillna(0)

# oral_clade as bool (binary novelty; informative split among candidates)
df["oral_clade"] = df["oral_clade"].astype(str).str.strip().str.lower().isin(["true", "1"])
cn["oral_clade"] = cn["oral_clade"].astype(str).str.strip().str.lower().isin(["true", "1"])

log(f"  All phylotypes: {len(df):,} | Candidate novel: {len(cn):,}")

# ============================================================
# Fig 1: Pipeline funnel
# ============================================================
log("Fig 1: pipeline funnel...")

# counts read from the actual intermediate tables (2026-08-14)
import json
with open(PIPE_V2 / "02_summary.json") as f:
    s02 = json.load(f)
n_chimera   = s02["total_candidates"]           # 01 chimera-filtered
n_firewall  = s02["clean"]                      # 02 firewall clean
n_floor     = len(pd.read_csv(PIPE_V2 / "04_support_floor_passed.tsv", sep="\t", low_memory=False,
                              usecols=["derep_id"]))
n_placed    = len(pd.read_csv(PIPE_V2 / "05_asv_phylotype_assignments.tsv", sep="\t", low_memory=False,
                              usecols=["derep_id"]))
n_phylo     = len(pd.read_csv(PIPE_V2 / "05_phylotypes_all.tsv", sep="\t", low_memory=False,
                              usecols=["phylotype_id"]))
stages = [
    ("Chimera-filtered\ncandidates",            n_chimera),
    ("Contamination\nfirewall (clean)",          n_firewall),
    ("Support floor\n≥2 BP, ≥50 reads",          n_floor),
    ("SEPP-placed\nASVs",                        n_placed),
    ("Phylotypes\n(0.03 species boundary)",      n_phylo),
    ("Candidate novel\n(≥10 subjects)",          len(cn)),
]
labels = [s[0] for s in stages]
values = [s[1] for s in stages]
w = 0.7
x = np.arange(len(stages))

fig, ax = plt.subplots(figsize=(8.5, 4.6))
# funnel: decreasing bar widths centered
widths = np.array(values) / max(values)
bars = []
for i, (lab, v) in enumerate(stages):
    b = ax.bar(i, 1, width=w * widths[i], bottom=0,
               color=TEAL if i < len(stages) - 1 else RED,
               edgecolor="white", linewidth=0.6, alpha=0.9)
    bars.append(b)
    ax.text(i, 0.5, f"{v:,}", ha="center", va="center",
            fontsize=12, fontweight="bold",
            color="white" if i == len(stages) - 1 else INK)

ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=9)
ax.set_ylim(0, 1.15)
ax.set_yticks([])
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.set_title("pipeline_v2 — ASV to candidate-novel-phylotype funnel", fontsize=13, pad=14)
plt.tight_layout()
fig.savefig(OUT / "fig_funnel_pipeline.png", bbox_inches="tight")
plt.close(fig)
log("  saved fig_funnel_pipeline.png")

# ============================================================
# Fig 2: Novelty classification
# ============================================================
log("Fig 2: novelty classification...")

# binary novelty (nearest species in HOMD?)
nov_all = df["novelty_class"].value_counts()
order = ["homd_captured", "homd_uncaptured"]
nov_all = nov_all.reindex([o for o in order if o in nov_all.index])

# among candidates (all homd_uncaptured): the informative split is oral_clade
cn_oral = cn["oral_clade"].value_counts().reindex([True, False]).fillna(0).astype(int)

fig, axes = plt.subplots(1, 2, figsize=(10, 4.6))

# Panel A: all phylotypes — binary novelty
ax = axes[0]
cols_all = [GRAY if c == "homd_captured" else RED for c in nov_all.index]
bars = ax.bar(range(len(nov_all)), nov_all.values, color=cols_all,
              edgecolor="white", linewidth=0.6, width=0.62)
for b, v in zip(bars, nov_all.values):
    ax.text(b.get_x() + b.get_width() / 2, v + 250, f"{v:,}",
            ha="center", fontsize=10, fontweight="bold")
ax.set_xticks(range(len(nov_all)))
ax.set_xticklabels(["HOMD-captured", "HOMD-uncaptured"], fontsize=8.5)
ax.set_ylabel("phylotypes")
ax.set_title(f"All phylotypes (n={len(df):,})", fontsize=11)
ax.set_ylim(0, max(nov_all.values) * 1.12)
ax.grid(axis="y", color=GRID, linewidth=0.7)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Panel B: candidate novel — by oral clade
ax = axes[1]
cols_cn = [BLUE if o else RED for o in cn_oral.index]   # True = within HOMD genus
bars = ax.bar(range(len(cn_oral)), cn_oral.values, color=cols_cn,
              edgecolor="white", linewidth=0.6, width=0.5)
for b, v in zip(bars, cn_oral.values):
    ax.text(b.get_x() + b.get_width() / 2, v + 60, f"{v:,}",
            ha="center", fontsize=10, fontweight="bold")
ax.set_xticks(range(len(cn_oral)))
ax.set_xticklabels(["within HOMD genus", "outside HOMD genus"], fontsize=8.5)
ax.set_ylabel("phylotypes")
ax.set_title(f"Candidate novel, ≥10 subjects (n={len(cn):,})", fontsize=11)
ax.set_ylim(0, max(cn_oral.values) * 1.12)
ax.grid(axis="y", color=GRID, linewidth=0.7)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.suptitle("Novelty classification of pipeline_v2 phylotypes (binary: nearest species in HOMD?)",
             fontsize=13, y=1.02)
plt.tight_layout()
fig.savefig(OUT / "fig_novelty_classification.png", bbox_inches="tight")
plt.close(fig)
log("  saved fig_novelty_classification.png")

# ============================================================
# Fig 3: Candidate support metrics
# ============================================================
log("Fig 3: candidate support metrics...")

fig, axes = plt.subplots(2, 2, figsize=(10, 7.2))

# Panel A: Bioprojects
ax = axes[0, 0]
bp_bins = [1, 2, 5, 10, 20, 50, 100]
bp_counts, bp_edges = np.histogram(cn["max_n_bioprojects"], bins=bp_bins)
bp_labels = [f"{bp_bins[i]}-{bp_bins[i+1] - 1}" for i in range(len(bp_bins) - 1)]
bars = ax.bar(bp_labels, bp_counts, color=BLUE, edgecolor="white", linewidth=0.6, width=0.62)
for b, v in zip(bars, bp_counts):
    ax.text(b.get_x() + b.get_width() / 2, v + 20, f"{v:,}", ha="center", fontsize=8.5)
ax.set_xlabel("BioProjects per phylotype")
ax.set_ylabel("phylotypes")
ax.set_title("Independent BioProjects", fontsize=11)
ax.grid(axis="y", color=GRID, linewidth=0.7)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Panel B: Prevalence (distinct subjects)
ax = axes[0, 1]
prev_bins = [10, 20, 50, 100, 200, 500, 1000, 5000]
prev_counts, prev_edges = np.histogram(cn["max_subject_prevalence"], bins=prev_bins)
prev_labels = [f"{prev_bins[i]}-{prev_bins[i+1] - 1}" if i < len(prev_bins) - 2
               else f"{prev_bins[i]}+" for i in range(len(prev_bins) - 1)]
bars = ax.bar(prev_labels, prev_counts, color=TEAL, edgecolor="white", linewidth=0.6, width=0.62)
for b, v in zip(bars, prev_counts):
    ax.text(b.get_x() + b.get_width() / 2, v + 15, f"{v:,}", ha="center", fontsize=8.5)
ax.set_xlabel("subjects (max subject prevalence)")
ax.set_ylabel("phylotypes")
ax.set_title("Subject prevalence", fontsize=11)
ax.grid(axis="y", color=GRID, linewidth=0.7)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Panel C: ASVs per phylotype
ax = axes[1, 0]
asv_bins = [1, 2, 3, 5, 10, 20, 50, 150]
asv_counts, _ = np.histogram(cn["n_asvs"], bins=asv_bins)
asv_labels = [f"{asv_bins[i]}-{asv_bins[i+1] - 1}" for i in range(len(asv_bins) - 1)]
bars = ax.bar(asv_labels, asv_counts, color=GRAY, edgecolor="white", linewidth=0.6, width=0.62)
for b, v in zip(bars, asv_counts):
    ax.text(b.get_x() + b.get_width() / 2, v + 15, f"{v:,}", ha="center", fontsize=8.5)
ax.set_xlabel("member ASVs")
ax.set_ylabel("phylotypes")
ax.set_title("Cross-region ASV merging", fontsize=11)
ax.grid(axis="y", color=GRID, linewidth=0.7)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Panel D: Continents per phylotype
ax = axes[1, 1]
cont_bins = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5]
cont_counts, _ = np.histogram(cn["n_continents"], bins=cont_bins)
cont_labels = ["1", "2", "3", "4", "5", "6"]
bars = ax.bar(cont_labels, cont_counts, color=RED, edgecolor="white", linewidth=0.6, width=0.62)
for b, v in zip(bars, cont_counts):
    ax.text(b.get_x() + b.get_width() / 2, v + 20, f"{v:,}", ha="center", fontsize=8.5)
ax.set_xlabel("continents")
ax.set_ylabel("phylotypes")
ax.set_title("Geographic spread", fontsize=11)
ax.grid(axis="y", color=GRID, linewidth=0.7)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.suptitle(f"Support metrics — candidate novel phylotypes (n={len(cn):,})",
             fontsize=13, y=1.0)
plt.tight_layout()
fig.savefig(OUT / "fig_candidate_support.png", bbox_inches="tight")
plt.close(fig)
log("  saved fig_candidate_support.png")

# ============================================================
# Fig 4: Continent profile
# ============================================================
log("Fig 4: continent profile...")

def parse_cont(val):
    if pd.isna(val) or val == "":
        return set()
    s = str(val).strip("{}")
    return set(x.strip().strip("'\"").strip() for x in s.split(";") if x.strip())

CONT_ORDER = ["Asia", "Europe", "North America", "South America", "Africa", "Oceania"]
cont_names = {"N.America": "North America", "North America": "North America",
              "S.America": "South America", "South America": "South America"}
def norm(c):
    return cont_names.get(c, c)

cont_counts = {c: 0 for c in CONT_ORDER}
for v in cn["continents"].dropna():
    for c in parse_cont(v):
        c = norm(c)
        if c in cont_counts:
            cont_counts[c] += 1

# split by oral clade (outside vs within HOMD genus)
cont_out = {c: 0 for c in CONT_ORDER}
cont_in  = {c: 0 for c in CONT_ORDER}
for _, row in cn.iterrows():
    for c in parse_cont(row["continents"]):
        c = norm(c)
        if c not in cont_counts:
            continue
        if row["oral_clade"]:
            cont_in[c] += 1
        else:
            cont_out[c] += 1

fig, ax = plt.subplots(figsize=(9, 4.6))
x = np.arange(len(CONT_ORDER))
wbar = 0.6
b1 = ax.bar(x, [cont_out[c] for c in CONT_ORDER], width=wbar, color=RED,
            edgecolor="white", linewidth=0.5, label="outside HOMD genus")
b2 = ax.bar(x, [cont_in[c] for c in CONT_ORDER], width=wbar, bottom=[cont_out[c] for c in CONT_ORDER],
            color=BLUE, edgecolor="white", linewidth=0.5, label="within HOMD genus")
for xi, c in enumerate(CONT_ORDER):
    tot = cont_out[c] + cont_in[c]
    if tot > 0:
        ax.text(xi, tot + 60, f"{tot:,}", ha="center", fontsize=9)
ax.set_xticks(x)
ax.set_xticklabels(CONT_ORDER, fontsize=9)
ax.set_ylabel("candidate phylotypes present on continent")
ax.set_ylim(0, max(cont_out[c] + cont_in[c] for c in CONT_ORDER) * 1.12)
ax.grid(axis="y", color=GRID, linewidth=0.7)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.legend(frameon=False, fontsize=9, loc="upper right")
ax.set_title(f"Continent distribution of candidate novel phylotypes (n={len(cn):,})",
             fontsize=12, pad=12)
plt.tight_layout()
fig.savefig(OUT / "fig_candidate_geography.png", bbox_inches="tight")
plt.close(fig)
log("  saved fig_candidate_geography.png")

# ============================================================
log(f"\nAll figures saved to: {OUT}")
log("Done!")
