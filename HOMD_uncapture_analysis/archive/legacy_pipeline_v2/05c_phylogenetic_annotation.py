#!/usr/bin/env python3
"""
pipeline_v2 Step 05c: Phylogenetic (placement-clade) annotation — Option A.

For each ASV placed on the GG2 SEPP tree, annotate by the taxonomy of the
reference tips in the descendant clade of its placement edge:

  1. parse the jplace reference tree (edge numbers [N] = edge labels);
  2. for each ASV's placement edge → descendant reference tips;
  3. map tips to GG2 2022.10 taxonomy (genus / species);
  4. assign the MAJORITY genus (and species) of those tips;
  5. record the placement LWR as confidence (low LWR = uncertain placement).

Unlike gappa examine assign (full placement distribution × LWR weighting), this
uses the single best placement, but it is fast, controllable, and uses the actual
tree structure (the clade the ASV is inserted into).

Outputs (pipeline_v2/):
  - 05c_asv_phylo_annotation.tsv : per-ASV annotation (all placed ASVs)
  - 05c_phylotype_phylo_annotation.tsv : aggregated per phylotype (all phylotypes)
"""
import json
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"
GAPPA = PIPE_V2 / "gappa"
log = lambda s: print(f"[05c_phylo_ann] {s}", flush=True)
sys.setrecursionlimit(1000000)

# ============================================================
# 1. Parse the jplace reference tree
# ============================================================
log("Parsing jplace reference tree...")
with open(BASE / "pipeline/05_sepp/placements_export/placements.json") as f:
    data = json.load(f)
s = data["tree"]


def parse_jplace_tree(s):
    L = len(s)
    i = 0
    nodes = []

    def mk(par, edge=None, tip=None):
        nd = {"par": par, "kids": [], "edge": edge, "tip": tip}
        nodes.append(nd)
        idx = len(nodes) - 1
        if par is not None:
            nodes[par]["kids"].append(idx)
        return idx

    stack = [mk(None)]
    while i < L:
        c = s[i]
        if c == "(":
            stack.append(mk(stack[-1]))
            i += 1
        elif c == ")":
            closed = stack.pop()
            i += 1
            if i < L and s[i] == ":":
                while i < L and s[i] not in "[],);":
                    i += 1
            if i < L and s[i] == "[":
                k = s.index("]", i)
                nodes[closed]["edge"] = int(s[i + 1:k])
                i = k + 1
        elif c == ",":
            i += 1
        elif c == ";":
            break
        else:
            par = stack[-1]
            j = i
            while j < L and s[j] not in ":,()[];":
                j += 1
            nd = mk(par, tip=s[i:j])
            i = j
            if i < L and s[i] == ":":
                while i < L and s[i] not in "[],);":
                    i += 1
            if i < L and s[i] == "[":
                k = s.index("]", i)
                nodes[nd]["edge"] = int(s[i + 1:k])
                i = k + 1
    return nodes


nodes = parse_jplace_tree(s)
n_tips = sum(1 for n in nodes if n["tip"] is not None)
log(f"  nodes: {len(nodes):,} | tips: {n_tips:,}")

edge_to_node = {}
for idx, n in enumerate(nodes):
    if n["edge"] is not None:
        edge_to_node[n["edge"]] = idx


def desc_tips(i):
    """descendant tip names of node i (iterative post-order to avoid recursion limits)"""
    n = nodes[i]
    if n["tip"]:
        return [n["tip"]]
    out = []
    st = list(n["kids"])
    while st:
        k = st.pop()
        kn = nodes[k]
        if kn["tip"]:
            out.append(kn["tip"])
        else:
            st.extend(kn["kids"])
    return out


# ============================================================
# 2. Taxonomy (GG2 2022.10, tree-tip IDs) → genus / species
# ============================================================
log("Loading GG2 2022.10 taxonomy...")
tax = {}
with open(GAPPA / "gg2_2022.10_tree_taxonomy.tsv") as f:
    for line in f:
        parts = line.rstrip().split("\t")
        if len(parts) < 2:
            continue
        genus = species = None
        for tok in parts[1].split(";"):
            tok = tok.strip()
            if tok.startswith("g__") and tok != "g__":
                genus = tok[3:]
            elif tok.startswith("s__") and tok != "s__":
                species = tok[3:]
        tax[parts[0]] = (genus, species)
log(f"  tips with taxonomy: {len(tax):,}")

# ============================================================
# 3. ASVs to annotate — ALL placed ASVs (support-floor set, step 05)
#    (2026-08-13: expanded from strict-225 members so that every
#     phylotype gets a placement-clade species for novelty calling)
# ============================================================
log("Loading placed ASVs (step 05 output)...")
asv = pd.read_csv(PIPE_V2 / "05_asv_phylotype_assignments.tsv", sep="\t", low_memory=False)
asv["edge_num"] = pd.to_numeric(asv["edge_num"], errors="coerce")
asv["like_weight_ratio"] = pd.to_numeric(asv["like_weight_ratio"], errors="coerce").fillna(0)
log(f"  ASVs with placement: {len(asv):,}")

# ============================================================
# 4. Annotate each ASV by its placement clade
# ============================================================
log("Annotating by placement clade...")
rows = []
n_edge_missing = 0
for _, r in asv.iterrows():
    did = r["derep_id"]
    edge = int(r["edge_num"]) if pd.notna(r["edge_num"]) else None
    if edge is None or edge not in edge_to_node:
        n_edge_missing += 1
        rows.append({"derep_id": did, "phylotype_id": r.get("phylotype_id", ""),
                     "edge_num": edge, "n_desc_tips": 0, "n_tax": 0,
                     "genus": None, "genus_frac": 0.0, "species": None, "species_frac": 0.0,
                     "lwr": r["like_weight_ratio"], "pend": r.get("pendant_length", None)})
        continue
    tips = desc_tips(edge_to_node[edge])
    genera = [tax[t][0] for t in tips if t in tax and tax[t][0]]
    species = [tax[t][1] for t in tips if t in tax and tax[t][1]]
    n_tax = len([t for t in tips if t in tax])
    genus = None
    gfrac = 0.0
    if genera:
        c = Counter(genera)
        genus, gcnt = c.most_common(1)[0]
        gfrac = gcnt / len(genera)
    sp = None
    sfrac = 0.0
    if species:
        c = Counter(species)
        sp, scnt = c.most_common(1)[0]
        sfrac = scnt / len(species)
    rows.append({"derep_id": did, "phylotype_id": r.get("phylotype_id", ""),
                 "edge_num": edge, "n_desc_tips": len(tips), "n_tax": n_tax,
                 "genus": genus, "genus_frac": gfrac, "species": sp, "species_frac": sfrac,
                 "lwr": r["like_weight_ratio"], "pend": r.get("pendant_length", None)})

ann = pd.DataFrame(rows)
log(f"  annotated: {len(ann):,} ASVs | edges not found: {n_edge_missing}")
log(f"  with genus: {ann['genus'].notna().sum():,} | median LWR: {ann['lwr'].median():.3f}")
ann.to_csv(PIPE_V2 / "05c_asv_phylo_annotation.tsv", sep="\t", index=False)
log(f"  saved 05c_asv_phylo_annotation.tsv")

# ============================================================
# 5. Aggregate to phylotype (majority genus across members)
# ============================================================
log("Aggregating to phylotype...")
phy_rows = []
for pid, g in ann.groupby("phylotype_id"):
    genus_ok = g["genus"].dropna()
    if len(genus_ok):
        c = Counter(genus_ok)
        g_agg, g_cnt = c.most_common(1)[0]
        g_frac = g_cnt / len(g)
    else:
        g_agg, g_frac = None, 0.0
    sp_ok = g["species"].dropna()
    if len(sp_ok):
        c = Counter(sp_ok)
        s_agg, s_cnt = c.most_common(1)[0]
        s_frac = s_cnt / len(g)
    else:
        s_agg, s_frac = None, 0.0
    phy_rows.append({"phylotype_id": pid, "n_members_annotated": len(g),
                     "phylo_genus": g_agg, "phylo_genus_frac": g_frac,
                     "phylo_species": s_agg, "phylo_species_frac": s_frac,
                     "median_lwr": g["lwr"].median(),
                     "min_genus_frac_member": genus_ok.map(lambda x: x).shape[0] if False else None})
phy = pd.DataFrame(phy_rows)
phy.to_csv(PIPE_V2 / "05c_phylotype_phylo_annotation.tsv", sep="\t", index=False)
log(f"  saved 05c_phylotype_phylo_annotation.tsv ({len(phy):,} phylotypes)")
log(f"  phylotypes with phylo_genus: {phy['phylo_genus'].notna().sum():,}")
log("\nDone!")
