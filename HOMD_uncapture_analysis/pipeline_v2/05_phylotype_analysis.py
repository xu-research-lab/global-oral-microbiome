#!/usr/bin/env python3
"""
pipeline_v2 Step 05: Phylotype Analysis (Tree-Edge Clustering)

Define phylotypes based on SEPP placement on GG2 backbone.
ASVs placed on the same tree edge (edge_num) → same phylotype.

This is the KEY advantage of SEPP over identity-based clustering:
phylotypes are defined by actual phylogenetic position, not arbitrary thresholds.

Inputs:
  - 04_support_floor_passed.tsv: ASVs passing support floor
  - sepp_placements_reused.tsv: reused SEPP placements
  - sepp_missing/placements_export/placements.json: new SEPP for 140 missing

Outputs:
  - 05_phylotypes_all.tsv: all phylotypes with member details
  - 05_phylotype_stats.tsv: per-phylotype statistics
  - 05_summary.json
"""
import pandas as pd
import numpy as np
import json, os, sys
from pathlib import Path
from collections import defaultdict

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"
os.makedirs(PIPE_V2, exist_ok=True)
log = lambda s: print(f"[05_phylo] {s}", flush=True)

# ============================================================
# 1. Load support floor passed ASVs
# ============================================================
log("Loading support floor passed ASVs...")
passed = pd.read_csv(PIPE_V2 / "04_support_floor_passed.tsv", sep="\t", low_memory=False)
log(f"  {len(passed):,} ASVs")

# ============================================================
# 2. Load SEPP placements (reused + new)
# ============================================================
log("\nLoading SEPP placements...")

# Reused placements
reused_path = PIPE_V2 / "sepp_placements_reused.tsv"
sepp_reused = pd.read_csv(reused_path, sep="\t", low_memory=False)
log(f"  Reused: {len(sepp_reused):,}")

# New placements for 140 missing ASVs
missing_json = PIPE_V2 / "sepp_missing" / "placements_export" / "placements.json"
new_rows = []
if missing_json.exists():
    log(f"  Loading new SEPP placements from {missing_json}...")
    with open(missing_json) as f:
        data = json.load(f)

    for entry in data['placements']:
        nm = entry.get('nm', [])
        if not nm or not nm[0]:
            continue
        query_id = nm[0][0]

        positions = entry.get('p', [])
        if not positions:
            continue

        # Best placement: highest like_weight_ratio (field index 2)
        best_idx = np.argmax([p[2] for p in positions])
        best = positions[best_idx]
        distal_length, edge_num, lwr, likelihood, pendant_length = best

        new_rows.append({
            'query_id': query_id,
            'edge_num': int(edge_num),
            'distal_length': float(distal_length),
            'like_weight_ratio': float(lwr),
            'likelihood': float(likelihood),
            'pendant_length': float(pendant_length),
        })
    log(f"  New placements: {len(new_rows)}")
else:
    log(f"  WARNING: {missing_json} not found - skipping {140 if 'missing' in dir() else '?'} ASVs")

sepp_new = pd.DataFrame(new_rows)

# Merge all placements
sepp_all = pd.concat([sepp_reused, sepp_new], ignore_index=True)
sepp_all = sepp_all.drop_duplicates(subset="query_id", keep="first")
log(f"  Merged placements: {len(sepp_all):,}")

log(f"  Median pendant_length: {sepp_all['pendant_length'].median():.4f}")
log(f"  Median like_weight_ratio: {sepp_all['like_weight_ratio'].median():.4f}")

# Compute phylotype = EDGE_{edge_num}
sepp_all["phylotype_id"] = "EDGE_" + sepp_all["edge_num"].astype(str)
log(f"  Unique phylotypes (edges): {sepp_all['phylotype_id'].nunique():,}")

# ============================================================
# 3. Merge SEPP placements with ASV metadata
# ============================================================
log("\nMerging placements with ASV metadata...")
df = passed.merge(sepp_all, left_on="derep_id", right_on="query_id", how="left")
n_unmatched = df["phylotype_id"].isna().sum()
if n_unmatched:
    log(f"  WARNING: {n_unmatched} ASVs without SEPP placement")
    missing_ids = df[df["phylotype_id"].isna()]["derep_id"].tolist()
    log(f"  Missing: {missing_ids[:10]}...")
df = df[df["phylotype_id"].notna()].copy()
log(f"  Merged: {len(df):,} ASVs with placement")

# ---- Subject-level prevalence (01b) ----
# distinct (bioproject, subject_id) individuals, not samples. 2026-08-13.
subj = pd.read_csv(PIPE_V2 / "01b_derep_subject_prevalence.tsv", sep="\t", low_memory=False)
df = df.merge(subj[["derep_id", "subject_prevalence"]], on="derep_id", how="left")
df["subject_prevalence"] = pd.to_numeric(df["subject_prevalence"], errors="coerce").fillna(0).astype(int)
log(f"  subject_prevalence merged | max={df['subject_prevalence'].max():,}")

# ============================================================
# 3b. Within-edge phylogenetic splitting (T = 0.05)
#     ASVs placed on the same SEPP edge can still be phylogenetically far apart.
#     Distance = placement-position difference: |distal_i - distal_j| + |pend_i - pend_j|.
#     (NOT pend_i + pend_j: novel taxa far from the GG2 backbone have large pendant
#      lengths, so the sum over-splits same-taxon variants. Decision 2026-08-05.)
#     Split by single-linkage connected components at T = 0.05 (≈ species boundary).
# ============================================================
log("\nSplitting phylotypes by within-edge position distance (T=0.05)...")
SPLIT_T = 0.05
df["dist"] = pd.to_numeric(df.get("distal_length"), errors="coerce").fillna(0)
df["pend"] = pd.to_numeric(df.get("pendant_length"), errors="coerce").fillna(0)

new_phy_ids = {}
n_split = 0
n_multi = 0
for edge, g in df.groupby("edge_num"):
    idxs = g.index.to_numpy()
    dists = g["dist"].to_numpy()
    pends = g["pend"].to_numpy()
    n = len(idxs)
    if n < 2:
        for i in idxs:
            new_phy_ids[i] = f"EDGE_{int(edge)}"
        continue
    n_multi += 1
    parent = list(range(n))
    def _find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for i in range(n):
        for j in range(i + 1, n):
            d = abs(dists[i] - dists[j]) + abs(pends[i] - pends[j])
            if d <= SPLIT_T:
                parent[_find(i)] = _find(j)
    comp = {}
    comp_ids = []
    for i in range(n):
        r = _find(i)
        if r not in comp:
            comp[r] = len(comp)
        comp_ids.append(comp[r])
    ncomp = len(comp)
    if ncomp > 1:
        n_split += 1
    for i, c in zip(idxs, comp_ids):
        new_phy_ids[i] = f"EDGE_{int(edge)}" if ncomp == 1 else f"EDGE_{int(edge)}_{c}"
df["phylotype_id"] = df.index.map(new_phy_ids)
n_unique_after = df["phylotype_id"].nunique()
log(f"  multi-ASV edges: {n_multi:,} | split into >1 phylotype: {n_split:,} ({n_split/max(n_multi,1)*100:.1f}%)")
log(f"  phylotypes after split: {n_unique_after:,}")

# ============================================================
# 4. Phylotype-level aggregation
# ============================================================
log("\nComputing phylotype-level statistics...")

# Parse regions and bioprojects from ASV-level data if available
def parse_set(val):
    """Parse region/bioproject list - handles string set literals"""
    if pd.isna(val) or val == "":
        return set()
    if isinstance(val, set):
        return val
    s = str(val)
    # Remove Python set notation { }
    s = s.strip("{}")
    # Handle semicolons, commas, or spaces as separators
    items = []
    for sep in ["; ", ";", ", ", ","]:
        if sep in s:
            items = [x.strip().strip("'\"") for x in s.split(sep)]
            break
    if not items:
        items = [s.strip().strip("'\"")]
    return set(i for i in items if i)

phylo = df.groupby("phylotype_id").agg(
    n_asvs=("derep_id", "count"),
    rep_derep_id=("derep_id", "first"),
    rep_edge_num=("edge_num", "first"),
    rep_pendant_length=("pendant_length", "first"),
    rep_gg2_genus=("gg2_genus", "first"),
    member_derep_ids=("derep_id", lambda x: ";".join(x)),
    member_gg2_genera=("gg2_genus", lambda x: ";".join(sorted(set(str(g) for g in x if pd.notna(g))))),
    total_reads_sum=("total_reads", "sum"),
    max_n_bioprojects=("n_bioprojects", "max"),
    avg_n_bioprojects=("n_bioprojects", "mean"),
    max_prevalence=("prevalence_samples", "max"),
    max_subject_prevalence=("subject_prevalence", "max"),
).reset_index()

# Count unique regions per phylotype
if "amplicon_regions" in df.columns:
    phylo["regions"] = phylo["phylotype_id"].map(
        df.groupby("phylotype_id")["amplicon_regions"].apply(
            lambda x: ";".join(sorted(set(
                r for val in x.dropna() for r in parse_set(val)
            )))
        )
    )

# Count unique bioprojects per phylotype
if "bioprojects" in df.columns:
    phylo["n_bioprojects_total"] = phylo["phylotype_id"].map(
        df.groupby("phylotype_id")["bioprojects"].apply(
            lambda x: len(set(
                r for val in x.dropna() for r in parse_set(val)
            ))
        )
    )

# Count unique continents per phylotype
if "continents" in df.columns:
    phylo["continents"] = phylo["phylotype_id"].map(
        df.groupby("phylotype_id")["continents"].apply(
            lambda x: ";".join(sorted(set(
                r for val in x.dropna() for r in parse_set(val)
            )))
        )
    )
    phylo["n_continents"] = phylo["phylotype_id"].map(
        df.groupby("phylotype_id")["continents"].apply(
            lambda x: len(set(
                r for val in x.dropna() for r in parse_set(val)
            ))
        )
    )

# Contamination flags per phylotype
phylo["n_contaminated_asvs"] = phylo["phylotype_id"].map(
    df.groupby("phylotype_id")["contamination_flag"].sum()
)
phylo["any_contamination"] = phylo["n_contaminated_asvs"] > 0

# Order columns
first_cols = ["phylotype_id", "rep_derep_id", "n_asvs", "member_derep_ids",
              "rep_edge_num", "rep_pendant_length",
              "rep_gg2_genus", "member_gg2_genera",
              "total_reads_sum", "max_n_bioprojects", "avg_n_bioprojects",
              "max_prevalence", "max_subject_prevalence", "n_continents", "continents",
              "n_contaminated_asvs", "any_contamination"]
extra_cols = [c for c in phylo.columns if c not in first_cols]
phylo = phylo[first_cols + extra_cols]

# ============================================================
# 5. Summary statistics
# ============================================================
log(f"\n=== Phylotype Analysis Summary ===")
log(f"Total ASVs with placement: {len(df):,}")
log(f"Total phylotypes: {len(phylo):,}")
log(f"Mean ASVs/phylotype: {phylo['n_asvs'].mean():.2f}")
log(f"Median ASVs/phylotype: {phylo['n_asvs'].median():.0f}")

# Single vs multi-ASV phylotypes
single_asv = (phylo["n_asvs"] == 1).sum()
multi_asv = (phylo["n_asvs"] > 1).sum()
log(f"Single-ASV phylotypes: {single_asv:,}")
log(f"Multi-ASV phylotypes: {multi_asv:,} (cross-region unification)")

# ASV distribution in multi-ASV phylotypes
multi_df = phylo[phylo["n_asvs"] > 1]
if len(multi_df) > 0:
    log(f"\nMulti-ASV phylotype size distribution:")
    for lo, hi, label in [(2,3,"2-3"),(3,5,"3-4"),(5,10,"5-9"),
                            (10,20,"10-19"),(20,50,"20-49"),(50,1000,"50+")]:
        n = len(multi_df[(multi_df["n_asvs"] >= lo) & (multi_df["n_asvs"] < hi)])
        if n > 0:
            log(f"  {label} ASVs: {n:,} phylotypes")

# Bioproject distribution
if "n_bioprojects_total" in phylo.columns:
    log(f"\nPhylotype bioproject distribution:")
    for lo, hi, label in [(1,2,"1"),(2,5,"2-4"),(5,10,"5-9"),
                            (10,20,"10-19"),(20,50,"20-49"),(50,1000,"50+")]:
        n = len(phylo[(phylo["n_bioprojects_total"] >= lo) & (phylo["n_bioprojects_total"] < hi)])
        log(f"  {label} BPs: {n:,} phylotypes")

# GG2 genus distribution
log(f"\nTop 20 phylotype genera:")
genus_counts = phylo["rep_gg2_genus"].value_counts()
for genus, count in genus_counts.head(20).items():
    log(f"  {genus}: {count:,} phylotypes")

# ============================================================
# 6. Save
# ============================================================
phylo_out = PIPE_V2 / "05_phylotypes_all.tsv"
phylo.to_csv(phylo_out, sep="\t", index=False)
log(f"\nSaved: {phylo_out} ({os.path.getsize(phylo_out):,} bytes)")

# Also save the merged ASV-level table with phylotype assignments
asv_out = PIPE_V2 / "05_asv_phylotype_assignments.tsv"
df_out_cols = ["derep_id", "phylotype_id", "edge_num", "pendant_length",
               "like_weight_ratio", "likelihood", "distal_length",
               "gg2_genus", "n_bioprojects", "total_reads",
               "contamination_flag", "homd_uncaptured", "is_chimera"]
avail = [c for c in df_out_cols if c in df.columns]
df[avail].to_csv(asv_out, sep="\t", index=False)
log(f"Saved: {asv_out}")

summary = {
    "total_asvs_with_placement": int(len(df)),
    "total_phylotypes": int(len(phylo)),
    "single_asv_phylotypes": int(single_asv),
    "multi_asv_phylotypes": int(multi_asv),
    "mean_asvs_per_phylotype": round(phylo["n_asvs"].mean(), 2),
    "median_asvs_per_phylotype": int(phylo["n_asvs"].median()),
    "median_pendant_length": round(sepp_all["pendant_length"].median(), 4),
    "median_lwr": round(sepp_all["like_weight_ratio"].median(), 4),
    "unique_genera": len(genus_counts),
    "top_genera": {k: int(v) for k, v in genus_counts.head(10).items()},
}
with open(PIPE_V2 / "05_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
log(f"Saved: {PIPE_V2 / '05_summary.json'}")

log("Done!")
