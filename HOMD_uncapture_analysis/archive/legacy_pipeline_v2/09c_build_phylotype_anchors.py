#!/usr/bin/env python3
"""
pipeline_v2 Step 09c: Build phylotype anchor table for the circle tree.

Each phylotype is represented by ONE circle placed at its representative ASV tip
(rep_derep_id, always present in the tree). Circle size encodes n_members.
Species label = GG2 taxonomy species of gg2_target (fallback: genus, then phylotype id).

Output: 07_phylotype_anchors.tsv
"""
import pandas as pd
from pathlib import Path

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"
log = lambda s: print(f"[09c_anchor] {s}", flush=True)

# ---- strict candidates + novelty annotations ----
strict = pd.read_csv(PIPE_V2 / "07_strict_candidates.tsv", sep="\t")
nov = pd.read_csv(PIPE_V2 / "06_phylotypes_novelty_classified.tsv", sep="\t", low_memory=False)
st = strict.merge(nov[["phylotype_id", "rep_derep_id", "rep_gg2_genus", "oral_clade"]],
                  on="phylotype_id", how="left")
# 05e tree-distance species anchor (rep ASV's nearest-tip species at D<=0.03) is
# the preferred species source — 2026-08-14; 05c clade-majority is the fallback.
phy05 = pd.read_csv(PIPE_V2 / "05_phylotypes_all.tsv", sep="\t", low_memory=False)
st = st.merge(phy05[["phylotype_id", "anchor_species", "anchor_genus"]],
              on="phylotype_id", how="left")
ann = pd.read_csv(PIPE_V2 / "05c_phylotype_phylo_annotation.tsv", sep="\t", low_memory=False)
st = st.merge(ann[["phylotype_id", "phylo_genus", "phylo_species"]],
              on="phylotype_id", how="left")
# genus: prefer the 05d nearest-tip genus (already strain-stripped), else clade-majority
st["phylo_genus"] = st["anchor_genus"].where(st["anchor_genus"].notna(), st["phylo_genus"])
log(f"  strict phylotypes: {len(st)}")

# ---- species: 05d anchor > 05c clade-majority ----
st["species"] = st["anchor_species"].where(st["anchor_species"].notna(), st["phylo_species"])
# explicit species annotation flag (1 = s__ resolved; 0 = genus / phylotype fallback)
st["has_species"] = st["species"].notna().astype(int)

# ---- label: anchor_species > clade_species > genus > phylotype_id ----
def make_label(row):
    if pd.notna(row["anchor_species"]) and str(row["anchor_species"]).strip():
        return row["anchor_species"]
    if pd.notna(row["phylo_species"]) and str(row["phylo_species"]).strip():
        return row["phylo_species"]
    g = row["phylo_genus"] if pd.notna(row["phylo_genus"]) else row["rep_gg2_genus"]
    if pd.notna(g) and g not in ("Unknown", "Unannotated"):
        return g
    return row["phylotype_id"]

st["species_label"] = st.apply(make_label, axis=1)

# ---- anchor table (subject-level prevalence) ----
anchors = st[["phylotype_id", "rep_derep_id", "novelty_class", "species_label",
              "phylo_genus", "phylo_species", "rep_gg2_genus", "has_species",
              "oral_clade",
              "max_subject_prevalence", "max_prevalence",
              "n_bioprojects_total", "n_continents", "n_regions",
              "n_members", "n_flagged"]].copy()
anchors = anchors.rename(columns={"rep_derep_id": "anchor_tip"})
# deliver strongest candidates first by subject-level prevalence (distinct
# (bioproject, subject_id) individuals); support_score is not exported here so
# the raw table otherwise looks unsorted
anchors = anchors.sort_values("max_subject_prevalence", ascending=False)
anchors.to_csv(PIPE_V2 / "07_phylotype_anchors.tsv", sep="\t", index=False)

n_species = anchors["species_label"].notna().sum()
log(f"  saved 07_phylotype_anchors.tsv ({len(anchors)} rows)")
log(f"  species labels: {n_species}/{len(anchors)}")
log("\nSample:")
print(anchors[["phylotype_id", "anchor_tip", "species_label", "max_subject_prevalence", "n_bioprojects_total", "n_members"]].head(10).to_string(index=False))
