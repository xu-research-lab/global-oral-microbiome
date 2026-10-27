#!/usr/bin/env python
"""pipeline_v3: distance from each reference sequence to its nearest HOMD reference (background distribution)"""
import time, random
from pathlib import Path
import pandas as pd
from ete3 import Tree

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V3 = BASE / "pipeline_v3"
TREE = BASE / "pipeline/05_sepp/tree_export/tree.nwk"
TAXONOMY = BASE / "pipeline_v2/gappa/gg2_2022.10_tree_taxonomy.tsv"
HOMD_FASTA = Path("/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03_full.fasta")
OUT = PIPE_V3 / "12h_genus_density"
OUT.mkdir(exist_ok=True)
log = lambda s: print(f"[12h_bg] {s}", flush=True)

tax = {}
with open(TAXONOMY) as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) == 2:
            tax[parts[0]] = parts[1]

def get_level(lineage, level):
    for part in lineage.split(";"):
        part = part.strip()
        if part.startswith(level + "__"):
            v = part[len(level)+2:].strip()
            return v if v else None
    return None

homd_species, homd_epithets = set(), {}
with open(HOMD_FASTA) as f:
    for line in f:
        if line.startswith(">"):
            parts = line[1:].split()
            if len(parts) >= 4:
                sp = f"{parts[2].strip()} {parts[3].strip()}"
                homd_species.add(sp)
                homd_epithets.setdefault(parts[3].strip(), set()).add(sp)

def renamed_in_homd(sp, homd_epithets):
    if not sp:
        return False
    parts = sp.split(" ")
    if len(parts) < 2:
        return False
    return len(homd_epithets.get(parts[-1].strip(), ())) == 1

log("Loading tree ...")
t0 = time.time()
tree = Tree(str(TREE), format=1, quoted_node_names=True)
log(f"  loaded in {time.time()-t0:.0f}s")

homd_ref_set = set()
for leaf in tree.iter_leaves():
    if leaf.name.startswith("DEREP_"):
        continue
    lin = tax.get(leaf.name)
    sp = get_level(lin, "s") if lin else None
    if sp and (sp in homd_species or renamed_in_homd(sp, homd_epithets)):
        homd_ref_set.add(leaf.name)
log(f"HOMD ref tips: {len(homd_ref_set):,}")

def nearest_ref_dp(tree, ref_set):
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            if node.name in ref_set:
                node.best_ref, node.best_dist = node.name, 0.0
            else:
                node.best_ref, node.best_dist = None, float("inf")
        else:
            best_ref, best_dist = None, float("inf")
            for c in node.children:
                if c.best_ref is not None:
                    d = c.best_dist + (c.dist or 0.0)
                    if d < best_dist:
                        best_dist, best_ref = d, c.best_ref
            node.best_ref, node.best_dist = best_ref, best_dist

def nearest_ref_for_tip(tip):
    node, child = tip.up, tip
    dist_to_node = tip.dist or 0.0
    best_ref, best_dist = None, float("inf")
    while node is not None:
        for c in node.children:
            if c is child or c.best_ref is None:
                continue
            d = dist_to_node + c.best_dist + (c.dist or 0.0)
            if d < best_dist:
                best_dist, best_ref = d, c.best_ref
        child = node
        dist_to_node += (node.dist or 0.0)
        node = node.up
    return best_ref, best_dist

log("DP nearest HOMD ref ...")
t0 = time.time()
nearest_ref_dp(tree, homd_ref_set)
log(f"  DP done in {time.time()-t0:.0f}s")

# compute the distance for every non-HOMD reference tip (background)
rows = []
t0 = time.time()
for leaf in tree.iter_leaves():
    if leaf.name.startswith("DEREP_") or leaf.name in homd_ref_set:
        continue
    r, d = nearest_ref_for_tip(leaf)
    rows.append((leaf.name, d))
log(f"  background computed in {time.time()-t0:.0f}s | n={len(rows):,}")
bg = pd.DataFrame(rows, columns=["ref_id", "nearest_homd_ref_dist"])
bg.to_csv(OUT / "background_homd_dist.tsv", sep="\t", index=False)
log("saved background_homd_dist.tsv")
