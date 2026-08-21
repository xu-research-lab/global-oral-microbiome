#!/usr/bin/env python3
"""
pipeline_v2 Step 06d: Apply tree labels for the phylotype-circle figure.

Label the phylotypes that are the defensible "novel to the oral catalog" set:
  - phylogenetic annotation gave a REAL genus name (not a GTDB/UBA placeholder)
  - that genus is NOT present in HOMD
  - not an obviously non-oral lineage (plant/animal pathogens etc.)

The circle fill already encodes the oral clade (red = outside HOMD genus, blue =
in HOMD genus); the label shows "Genus species" (or "Genus" when no species was
resolved). All strict candidates are homd_uncaptured, so `tree_label_novel` now
just flags the labelled (top-30) subset.

Adds two columns to 07_phylotype_anchors.tsv:
  tree_label       : "Genus species" / "Genus" / NA
  tree_label_novel : 1 if labelled (homd_uncaptured by construction)

The number of labelled phylotypes is reported so the figure stays legible.
"""
import re
from pathlib import Path

import pandas as pd

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"
HOMD_FASTA = BASE.parent.parent / "HOMD_16S_rRNA_RefSeq_V16.03_full.fasta"
log = lambda s: print(f"[06d_label] {s}", flush=True)

# ---- HOMD genus set ----
homd_genus = set()
with open(HOMD_FASTA) as f:
    for line in f:
        if line.startswith(">"):
            parts = line.split()
            if len(parts) >= 3:
                homd_genus.add(parts[2].strip())
log(f"HOMD genera: {len(homd_genus)}")

# ---- anchors (09c already carries phylo_genus / phylo_species / novelty_class;
#      merging 05c again would collide -> _x/_y suffixes under pandas 2.x) ----
anchors = pd.read_csv(PIPE_V2 / "07_phylotype_anchors.tsv", sep="\t")
m = anchors.copy()

# ---- label criteria ----
m["base"] = m["phylo_genus"].str.split("_").str[0].fillna("")
m["not_homd"] = ~m["base"].isin(homd_genus)
m["real_genus"] = m["base"].str.match(r"^[A-Z][a-zA-Z]+$") & ~m["base"].str.startswith(
    ("UBA", "UMGS", "GCA", "GWA", "CAG", "J"))
NON_ORAL = {"Phytoplasma", "Snodgrassella", "Pseudoalteromonas", "Alteromonas",
            "Glaciecola", "Planktomarina", "Roseovarius", "Vibrio", "Photobacterium"}
m["labelled"] = m["not_homd"] & m["real_genus"] & ~m["base"].isin(NON_ORAL)

# ---- keep only the top-30 strongest labels (support_score) ----
# support_score = max_subject_prevalence x n_bioprojects_total x n_regions (same
# metric used to rank the strict candidates); 2026-08-05: label top-30 only;
# 2026-08-13: subject-level prevalence
m["support_score"] = (pd.to_numeric(m["max_subject_prevalence"], errors="coerce").fillna(0)
                      * pd.to_numeric(m["n_bioprojects_total"], errors="coerce").fillna(0)
                      * pd.to_numeric(m["n_regions"], errors="coerce").fillna(0))
labelled_ids = m.loc[m["labelled"], "phylotype_id"]
top30_ids = m.loc[m["labelled"]].sort_values("support_score", ascending=False).head(30)["phylotype_id"].tolist()
m["labelled"] = m["phylotype_id"].isin(top30_ids)
log(f"  labels after top-30 filter: {int(m['labelled'].sum())}")

def make_label(r):
    if not r["labelled"]:
        return None
    # GG2 s__ is a full binomial ("Aquabacter cavernae") — use it directly
    if pd.notna(r["phylo_species"]) and r["phylo_species"]:
        return r["phylo_species"]
    return r["phylo_genus"]

m["tree_label"] = m.apply(make_label, axis=1)
m["tree_label_novel"] = ((m["labelled"]) & (m["novelty_class"] == "homd_uncaptured")).astype(int)

n_labelled = m["tree_label"].notna().sum()
n_novel = int((m["tree_label_novel"] == 1).sum())
log(f"tree labels: {n_labelled}/{len(m)} phylotypes ({n_novel} homd_uncaptured)")

# ---- save back ----
out_cols = [c for c in anchors.columns if c not in ("tree_label", "tree_label_novel")]
out = m[out_cols + ["tree_label", "tree_label_novel"]]
out.to_csv(PIPE_V2 / "07_phylotype_anchors.tsv", sep="\t", index=False)
log("saved updated 07_phylotype_anchors.tsv")

# ---- preview ----
prev = out[out["tree_label"].notna()].sort_values("max_subject_prevalence", ascending=False)
for _, r in prev.head(20).iterrows():
    mark = "  [uncaptured]" if r["tree_label_novel"] == 1 else ""
    print(f"  {r['phylotype_id']:>16} | {r['tree_label']:<38} | s{int(r['max_subject_prevalence'])}{mark}")
log("Done!")
