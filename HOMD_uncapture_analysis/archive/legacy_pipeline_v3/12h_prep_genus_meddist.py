#!/usr/bin/env python
"""pipeline_v3: median distance from each genus's ASVs to the nearest HOMD reference, merged into genus_stats"""
import time, math
from pathlib import Path
from collections import defaultdict
import pandas as pd
from ete3 import Tree

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V3 = BASE / "pipeline_v3"
TREE = BASE / "pipeline/05_sepp/tree_export/tree.nwk"
TAXONOMY = BASE / "pipeline_v2/gappa/gg2_2022.10_tree_taxonomy.tsv"
OUT = PIPE_V3 / "12h_genus_density"
log = lambda s: print(f"[12h_med] {s}", flush=True)

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

log("Loading tree ...")
tree = Tree(str(TREE), format=1, quoted_node_names=True)
name2node = {leaf.name: leaf for leaf in tree.iter_leaves()}

tip_genus = {}
for name, leaf in name2node.items():
    if name.startswith("DEREP_"):
        continue
    lin = tax.get(name)
    tip_genus[name] = get_level(lin, "g") if lin else None

genus_tips = defaultdict(list)
for name, g in tip_genus.items():
    if g:
        genus_tips[g].append(name)

genus_node = {}
for g, tips in genus_tips.items():
    if len(tips) == 1:
        genus_node[g] = name2node[tips[0]]
    else:
        genus_node[g] = tree.get_common_ancestor([name2node[t] for t in tips])
genus_node_ids = {id(nd): g for g, nd in genus_node.items()}

asv_genus = {}
for name, leaf in name2node.items():
    if not name.startswith("DEREP_"):
        continue
    node = leaf
    while node is not None:
        g = genus_node_ids.get(id(node))
        if g is not None:
            asv_genus[name] = g
            break
        node = node.up
log(f"ASV->genus: {len(asv_genus):,}")

# master table: asv -> nearest_homd_ref_dist
m = pd.read_csv(PIPE_V3 / "07_master_asv_table.tsv", sep="\t", usecols=["asv_id", "nearest_homd_ref_dist", "contamination_flag"])
m["genus"] = m["asv_id"].map(asv_genus)
mg = m[m["genus"].notna()].copy()
med = mg.groupby("genus")["nearest_homd_ref_dist"].median().rename("median_asv_homd_dist")
n_contam = mg.groupby("genus")["contamination_flag"].sum().rename("n_contam_asv")
med = pd.concat([med, n_contam], axis=1).reset_index()

gs = pd.read_csv(OUT / "genus_stats.tsv", sep="\t")
gs = gs.merge(med, on="genus", how="left")
gs.to_csv(OUT / "genus_stats.tsv", sep="\t", index=False)
log(f"saved genus_stats.tsv with median distance | rows={len(gs)}")
