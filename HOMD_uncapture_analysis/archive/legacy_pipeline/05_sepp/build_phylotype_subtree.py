#!/usr/bin/env python3
"""
Rebuild phylotype subtree from full GG2+SEPP tree,
keeping actual DEREP ASV tips + sister GG2 reference tips.

Uses single-pass traversal to collect all needed tips efficiently.
Output: phylotype_subtree_corrected.nwk
"""
import pandas as pd
import sys, os, warnings
from pathlib import Path
from ete3 import Tree

warnings.filterwarnings("ignore")

PIPE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline")
ITOL = PIPE / "analysis_report/itol"
FULL_TREE = PIPE / "05_sepp/tree_export/tree.nwk"
CAND_NOVEL = PIPE / "06_phylotypes_sepp/phylotypes_candidate_novel_updated.tsv"

def log(s):
    print(f"[build-subtree] {s}", flush=True)

# ============================================================
# 1. Load candidate novel phylotypes and their member DEREP IDs
# ============================================================
log("Loading phylotype member ASVs...")
phy = pd.read_csv(CAND_NOVEL, sep="\t", low_memory=False)
log(f"  {len(phy)} candidate novel phylotypes")

# Collect all unique DEREP IDs across all 56 phylotypes
all_derep = set()
for ids in phy["member_derep_ids"].dropna():
    for did in str(ids).split(";"):
        did = did.strip()
        if did:
            all_derep.add(did)
log(f"  Total unique DEREP IDs: {len(all_derep):,}")

# DEREP → phylotype_id mapping
derep_to_phylo = {}
for _, row in phy.iterrows():
    for did in str(row["member_derep_ids"]).split(";"):
        did = did.strip()
        if did:
            derep_to_phylo[did] = row["phylotype_id"]

log("Loading full tree...")
t = Tree(str(FULL_TREE), format=1, quoted_node_names=True)
log(f"  Tree loaded: {len(t):,} tips")

# ============================================================
# 2. Single-pass: find DEREP tips of interest + collect IDs to keep
# ============================================================
log("Scanning tree for phylotype DEREP tips (single-pass)...")

# Strategy: for each leaf, if it's a DEREP in our set,
# go up one level (unnamed insertion parent) and collect all leaves
# from that parent's subtree.

tips_to_keep = set()
found_count = 0
seen_parents = set()

for leaf in t.iter_leaves():
    if leaf.name in all_derep:
        found_count += 1
        parent = leaf.up
        if parent is None:
            continue
        # The parent node is the unnamed insertion point
        # But we need to go up one more level to get the
        # context (sister GG2 refs)
        grandparent = parent.up
        if grandparent is not None:
            # Use grandparent to include more context
            if id(grandparent) not in seen_parents:
                seen_parents.add(id(grandparent))
                kids = grandparent.get_leaf_names()
                tips_to_keep.update(kids)
        else:
            if id(parent) not in seen_parents:
                seen_parents.add(id(parent))
                kids = parent.get_leaf_names()
                tips_to_keep.update(kids)

log(f"  Found DEREP in tree: {found_count}/{len(all_derep)}")
log(f"  Unique parent nodes: {len(seen_parents)}")
log(f"  Total tips to keep: {len(tips_to_keep):,}")

if found_count < len(all_derep):
    missing = all_derep - {l.name for l in t.iter_leaves()}
    log(f"  WARNING: Missing {len(missing)} DEREPs not found as tree tips")
    for m in list(missing)[:5]:
        log(f"    e.g. {m}")

# Analyze
derep_kept = {t for t in tips_to_keep if t.startswith("DEREP_")}
gg2_kept = {t for t in tips_to_keep if not t.startswith("DEREP_")}
log(f"  DEREP tips: {len(derep_kept):,}")
log(f"  GG2 ref tips: {len(gg2_kept):,}")

# ============================================================
# 3. Prune tree to only these tips
# ============================================================
log("Pruning tree...")
keep_list = list(tips_to_keep)
t.prune(keep_list, preserve_branch_length=True)
log(f"  Pruned to {len(t):,} tips")

# ============================================================
# 4. Rename DEREP tips: add phylotype ID and species annotation
# ============================================================
# We need to also label each DEREP tip with its phylotype info.
# Load phylotype_info.tsv for species annotations
phy_info = pd.read_csv(ITOL / "phylotype_info.tsv", sep="\t")
if "gg2_name" in phy_info.columns:
    phy_name_map = dict(zip(phy_info["id"], phy_info["gg2_name"]))
else:
    phy_name_map = {}

for leaf in t.get_leaves():
    if leaf.name.startswith("DEREP_"):
        phy_id = derep_to_phylo.get(leaf.name, "")
        if phy_id and phy_id in phy_name_map:
            nm = phy_name_map[phy_id]
            if pd.isna(nm) or nm == "":
                nm = ""
            # Keep DEREP name but append phylotype context
            # The visualization will use the phylotype grouping

# Actually, for ggtree, it's easier to keep DEREP names as-is
# and build a separate annotation file for colors/grouping.
# The tree should show biological positions, not renamed leaves.
log("  Keeping original tip names for ggtree annotation compatibility")

# ============================================================
# 5. Save
# ============================================================
out_nwk = ITOL / "phylotype_subtree_corrected.nwk"
t.write(outfile=str(out_nwk), format=1, quoted_node_names=True)
log(f"Saved: {out_nwk}")

# Quick stats
leaf_names = t.get_leaf_names()
n_derep = sum(1 for n in leaf_names if n.startswith("DEREP_"))
n_gg2 = sum(1 for n in leaf_names if not n.startswith("DEREP_"))
log(f"  Tips: {len(leaf_names):,} (DEREP={n_derep:,}, GG2={n_gg2:,})")

# Show a few GG2 refs that are sister to our phylotypes
gg2_samples = [n for n in leaf_names if not n.startswith("DEREP_")][:10]
log(f"  Sample GG2 sisters: {gg2_samples}")

log("Done!")
