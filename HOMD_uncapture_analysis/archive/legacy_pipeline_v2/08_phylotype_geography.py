#!/usr/bin/env python3
"""
pipeline_v2 Step 08: Phylum-level distribution of candidate novel phylotypes.

Rewritten 2026-08-13. The novelty decision is now placement-clade species vs HOMD
(no gg2_target / SILVA comparisons), so phylum is mapped from the placement-clade
GG2 genus (05c phylo_genus) through the GG2 taxonomy.tsv genus→phylum lookup.

Maps each candidate novel phylotype to phylum/class via its placement-clade genus,
then plots the distribution across the bacterial tree of life.

Input : 06_phylotypes_novelty_classified.tsv (candidate_novel subset)
Output: fig_phylotype_phylum.png + 08_candidate_novel_with_phylum.tsv
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"
OUT = PIPE_V2 / "analysis_report" / "figs"
OUT.mkdir(parents=True, exist_ok=True)
TAX = "/beegfs/dongbiao/greengene2/exported-taxonomy/taxonomy.tsv"

BLUE = "#4E79A7"
RED = "#E15759"
GRAY = "#BAB0AC"
TEAL = "#76B7B2"
INK = "#2B2B2B"
GRID = "#E6E6E6"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.edgecolor": "#BABABA",
    "axes.linewidth": 0.8,
    "axes.labelcolor": INK,
    "xtick.color": INK,
    "ytick.color": INK,
    "text.color": INK,
    "figure.dpi": 150,
})

log = lambda s: print(f"[08_phylo] {s}", flush=True)

# ============================================================
# 1. Load candidate novel phylotypes (placement-clade genus present)
# ============================================================
log("Loading candidate novel phylotypes...")
nov = pd.read_csv(PIPE_V2 / "06_phylotypes_novelty_classified.tsv", sep="\t", low_memory=False)
cn = nov[nov["candidate_novel"] == True].copy()
log(f"  {len(cn):,} candidate novel phylotypes")

# ============================================================
# 2. Load GG2 taxonomy → genus→phylum/class lookup
# ============================================================
log("Loading GG2 taxonomy.tsv...")
tax = pd.read_csv(TAX, sep="\t", low_memory=False, dtype=str)

def parse_ranks(taxon_str):
    """Parse 'd__;p__;c__;o__;f__;g__;s__' into dict of ranks."""
    out = {}
    for part in str(taxon_str).split(";"):
        part = part.strip()
        if len(part) >= 3 and part[1:3] == "__":
            out[part[0]] = part[3:]
    return out

# GTDB genus names keep strain tokens ("Eubacterium_G_180878"); base = plain stem.
genus_to_phylum = {}
genus_to_class = {}
base_to_phylum = {}
base_to_class = {}
for _, row in tax.iterrows():
    ranks = parse_ranks(row["Taxon"])
    g = ranks.get("g", "")
    if not g:
        continue
    ph = ranks.get("p", "Unknown")
    cl = ranks.get("c", "Unknown")
    genus_to_phylum.setdefault(g, ph)
    genus_to_class.setdefault(g, cl)
    base = g.split("_")[0]
    base_to_phylum.setdefault(base, ph)
    base_to_class.setdefault(base, cl)
log(f"  taxonomy entries: {len(tax):,} | genus lookups: {len(genus_to_phylum):,}")

def lookup_ranks(genus):
    """phylo_genus ('Roseibacillus_B') → (phylum, class)."""
    genus = str(genus)
    if genus in genus_to_phylum:
        return genus_to_phylum[genus], genus_to_class[genus]
    base = genus.split("_")[0]
    if base in base_to_phylum:
        return base_to_phylum[base], base_to_class[base]
    return "Unknown", "Unknown"

# ============================================================
# 3. Map phylotypes to phylum via placement-clade genus
# ============================================================
log("Mapping phylotype placement-clade genera to phylum/class...")
phylum_arr = []
class_arr = []
for _, row in cn.iterrows():
    g = row.get("phylo_genus", "")
    if not pd.notna(g) or not str(g).strip():
        g = row.get("rep_gg2_genus", "")          # fallback to rep ASV GG2 genus
    ph, cl = lookup_ranks(g)
    phylum_arr.append(ph)
    class_arr.append(cl)

cn["gg2_phylum"] = phylum_arr
cn["gg2_class"] = class_arr

n_mapped = (cn["gg2_phylum"] != "Unknown").sum()
log(f"  Mapped to phylum: {n_mapped:,} / {len(cn):,} ({n_mapped/len(cn)*100:.1f}%)")

# ============================================================
# 4. Phylum distribution
# ============================================================
log("\nPhylum distribution (top 20):")
phylum_counts = cn["gg2_phylum"].value_counts()
for ph, n in phylum_counts.head(20).items():
    log(f"  {ph}: {n:,}")

# Split by oral clade (binary novelty: all candidates are homd_uncaptured; the
# informative split is whether the placement-clade genus is in HOMD)
cn["is_oral_clade"] = cn["oral_clade"].astype(str).str.strip().str.lower().isin(["true", "1"])

# ============================================================
# Fig: Phylum distribution (split by novelty class)
# ============================================================
log("Fig: phylum distribution...")

TOP_N = 15
top_phyla = phylum_counts.head(TOP_N).index.tolist()

cn_plot = cn.copy()
cn_plot["phylum_plot"] = cn_plot["gg2_phylum"].apply(
    lambda p: p if p in top_phyla else "Other"
)
phylum_plot_order = cn_plot["phylum_plot"].value_counts().index.tolist()
phylum_plot_order = [p for p in phylum_plot_order if p != "Other"] + (["Other"] if "Other" in phylum_plot_order else [])

out_by_ph = cn_plot[~cn_plot["is_oral_clade"]]["phylum_plot"].value_counts()
in_by_ph  = cn_plot[cn_plot["is_oral_clade"]]["phylum_plot"].value_counts()

fig, ax = plt.subplots(figsize=(10, 5.6))
x = np.arange(len(phylum_plot_order))
out_vals = [out_by_ph.get(p, 0) for p in phylum_plot_order]
in_vals  = [in_by_ph.get(p, 0) for p in phylum_plot_order]

b1 = ax.barh(x, out_vals, color=RED, edgecolor="white", linewidth=0.4, height=0.6,
             label="outside HOMD genus")
b2 = ax.barh(x, in_vals, left=out_vals, color=BLUE, edgecolor="white", linewidth=0.4,
             height=0.6, label="within HOMD genus")

for xi, (nv, kv) in enumerate(zip(out_vals, in_vals)):
    tot = nv + kv
    if tot > 0:
        ax.text(tot + 20, xi, f"{tot:,}", va="center", fontsize=9, color=INK)

ax.set_yticks(x)
ax.set_yticklabels(phylum_plot_order, fontsize=9)
ax.set_xlabel("candidate novel phylotypes")
ax.set_title("Phylum distribution of candidate novel phylotypes (placement-clade GG2 taxonomy)",
             fontsize=12, pad=12)
ax.grid(axis="x", color=GRID, linewidth=0.7)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.legend(frameon=False, fontsize=9, loc="lower right")
ax.set_xlim(0, max(tot for tot in [sum(p) for p in zip(out_vals, in_vals)]) * 1.12)

plt.tight_layout()
fig.savefig(OUT / "fig_phylotype_phylum.png", bbox_inches="tight")
plt.close(fig)
log("  saved fig_phylotype_phylum.png")

# ============================================================
# 5. Save phylum-annotated table (subject-level prevalence)
# ============================================================
out_cols = ["phylotype_id", "novelty_class", "candidate_novel", "gg2_phylum", "gg2_class",
            "n_asvs", "max_n_bioprojects", "max_subject_prevalence", "total_reads_sum",
            "n_continents", "continents"]
avail = [c for c in out_cols if c in cn.columns]
cn[avail].to_csv(PIPE_V2 / "08_candidate_novel_with_phylum.tsv", sep="\t", index=False)
log(f"  Saved: pipeline_v2/08_candidate_novel_with_phylum.tsv")

log("\nDone!")
