#!/usr/bin/env python3
"""Cross-reference updated contamination flags onto existing SEPP phylotypes.

Does NOT regenerate phylotypes or novelty (those are sequence-based).
Adds: per-phylotype contamination evidence summary from updated integrated annotations.
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path
from collections import defaultdict

PIPE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline")

# Load updated integrated annotations
integ = pd.read_csv(PIPE / "05_integrated_candidates/integrated_candidate_annotations.tsv", sep="\t", low_memory=False)
log = lambda s: print(f"[update-phylo] {s}", flush=True)

log(f"Integrated annotations: {len(integ):,} candidates")

# Load existing SEPP phylotypes (all, not just candidate_novel)
phy_all = pd.read_csv(PIPE / "06_phylotypes_sepp/phylotypes_all.tsv", sep="\t", low_memory=False)
log(f"SEPP phylotypes (all): {len(phy_all):,}")

# Load existing SEPP novelty
nov = pd.read_csv(PIPE / "07_novelty_sepp/phylotypes_novelty_classified.tsv", sep="\t", low_memory=False)
log(f"Novelty classified: {len(nov):,}")

# ============================================================
# Build derep_id -> contamination info
# ============================================================
contam_cols = [
    "silva_genus", "contamination_concern", "contamination_concern_level",
    "genus_decontam_A_n_bps", "genus_decontam_B_contaminant",
    "catalog_contaminant", "literature_classification",
    "literature_oral_core", "literature_contaminant",
    "residence_pass",
]
derep_info = integ.set_index("derep_id")[contam_cols].to_dict("index")

# ============================================================
# Aggregate contamination evidence per phylotype
# ============================================================
def aggregate_phylotype_contam(row):
    """For a phylotype row, aggregate contamination evidence across members."""
    members = str(row["member_derep_ids"]).split(";")
    members = [m.strip() for m in members if m.strip()]

    n_total = len(members)
    n_contam_concern = 0
    n_auto_remove = 0
    n_flag_review = 0
    n_oral_core = 0
    n_decontam_a2 = 0
    n_decontam_b = 0
    n_catalog = 0
    n_residence_pass = 0
    concern_levels = []
    genera = []

    for did in members:
        info = derep_info.get(did, {})
        if not info:
            continue
        if info.get("contamination_concern", False):
            n_contam_concern += 1
        lit = str(info.get("literature_classification", "unknown"))
        if lit == "auto_remove": n_auto_remove += 1
        elif lit == "flag_review": n_flag_review += 1
        elif lit == "oral_core": n_oral_core += 1
        if info.get("genus_decontam_A_n_bps", 0) >= 2: n_decontam_a2 += 1
        if info.get("genus_decontam_B_contaminant", False): n_decontam_b += 1
        if info.get("catalog_contaminant", False): n_catalog += 1
        if info.get("residence_pass", False): n_residence_pass += 1
        concern_levels.append(str(info.get("contamination_concern_level", "unknown")))
        g = str(info.get("silva_genus", "Unknown"))
        if g != "Unknown" and g != "nan":
            genera.append(g)

    # Phylotype-level contamination concern
    has_concern = n_contam_concern > 0
    # If ALL members are flagged as auto_remove, the phylotype is likely contaminant
    if n_auto_remove >= n_total * 0.5:
        phylo_concern = "likely_contaminant"
    elif n_contam_concern >= n_total * 0.5:
        phylo_concern = "contamination_suspect"
    elif n_oral_core > 0:
        phylo_concern = "oral_core_protected"
    elif has_concern:
        phylo_concern = "minor_concern"
    else:
        phylo_concern = "clean"

    return pd.Series({
        "n_members_total": n_total,
        "n_members_contamination_concern": n_contam_concern,
        "n_members_auto_remove": n_auto_remove,
        "n_members_flag_review": n_flag_review,
        "n_members_oral_core": n_oral_core,
        "n_members_decontam_a2": n_decontam_a2,
        "n_members_decontam_b": n_decontam_b,
        "n_members_catalog": n_catalog,
        "n_members_residence_pass": n_residence_pass,
        "phylotype_contamination_concern": phylo_concern,
        "member_genera": ";".join(sorted(set(genera))) if genera else "Unknown",
        "dominant_concern_level": max(set(concern_levels), key=concern_levels.count) if concern_levels else "unknown",
    })

log("Aggregating contamination per phylotype...")
contam_agg = phy_all.apply(aggregate_phylotype_contam, axis=1)
phy_updated = pd.concat([phy_all, contam_agg], axis=1)

# ============================================================
# Merge novelty classification
# ============================================================
nov_cols = ["phylotype_id", "novelty_class", "homd_best_id", "gg2_best_id", "silva_best_id"]
nov_subset = nov[nov_cols].copy()
phy_updated = phy_updated.merge(nov_subset, on="phylotype_id", how="left")

# ============================================================
# Save updated tables
# ============================================================
OUTDIR = PIPE / "06_phylotypes_sepp"
phy_updated.to_csv(OUTDIR / "phylotypes_all_updated.tsv", sep="\t", index=False)

# Filter to candidate novel only
phy_novel = phy_updated[phy_updated["candidate_novel"] == True].copy()
phy_novel.to_csv(OUTDIR / "phylotypes_candidate_novel_updated.tsv", sep="\t", index=False)

log(f"Saved: phylotypes_all_updated.tsv ({len(phy_updated):,} rows)")
log(f"Saved: phylotypes_candidate_novel_updated.tsv ({len(phy_novel):,} rows)")

# ============================================================
# Summary
# ============================================================
log(f"\n{'='*60}")
log(f"UPDATED PHYLOTYPE CONTAMINATION SUMMARY")
log(f"{'='*60}")
log(f"Total phylotypes: {len(phy_updated):,}")
log(f"Candidate novel phylotypes: {len(phy_novel):,}")

log(f"\nPhylotype contamination concern (all phylotypes):")
for level in phy_updated["phylotype_contamination_concern"].value_counts().index:
    n = (phy_updated["phylotype_contamination_concern"] == level).sum()
    log(f"  {level}: {n:,}")

log(f"\nCandidate novel phylotype contamination:")
for level in phy_novel["phylotype_contamination_concern"].value_counts().index:
    n = (phy_novel["phylotype_contamination_concern"] == level).sum()
    log(f"  {level}: {n:,}")

log(f"\nNovelty classification of candidate phylotypes:")
for cls in phy_novel["novelty_class"].value_counts().index:
    n = (phy_novel["novelty_class"] == cls).sum()
    log(f"  {cls}: {n:,}")

# Show any phylotypes flagged as likely_contaminant
flagged = phy_novel[phy_novel["phylotype_contamination_concern"] == "likely_contaminant"]
if len(flagged) > 0:
    log(f"\n  ⚠️ Candidate phylotypes flagged as likely_contaminant:")
    for _, row in flagged.iterrows():
        log(f"    {row['phylotype_id']}: genera={row['member_genera']}, "
            f"auto_remove={row['n_members_auto_remove']}/{row['n_members_total']}, "
            f"novelty={row.get('novelty_class','?')}")

phy_concern = phy_novel[phy_novel["phylotype_contamination_concern"] == "contamination_suspect"]
if len(phy_concern) > 0:
    log(f"\n  ⚠️ Candidate phylotypes flagged as contamination_suspect:")
    for _, row in phy_concern.head(10).iterrows():
        log(f"    {row['phylotype_id']}: genera={row['member_genera']}, "
            f"concern={row['n_members_contamination_concern']}/{row['n_members_total']}")

# Oral-core protected
oral_core_phy = phy_novel[phy_novel["phylotype_contamination_concern"] == "oral_core_protected"]
log(f"\n  Oral-core protected candidate phylotypes: {len(oral_core_phy)}")
for _, row in oral_core_phy.iterrows():
    log(f"    {row['phylotype_id']}: genera={row['member_genera']}")

log("\nDone!")
