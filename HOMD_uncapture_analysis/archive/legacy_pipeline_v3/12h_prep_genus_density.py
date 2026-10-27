#!/usr/bin/env python
"""pipeline_v3: genus-level ASV insertion density aggregation and genus tree construction (input for the 12h visualisation)"""
import time, math
from pathlib import Path
from collections import defaultdict
import pandas as pd
from ete3 import Tree

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V3 = BASE / "pipeline_v3"
TREE = BASE / "pipeline/05_sepp/tree_export/tree.nwk"
TAXONOMY = BASE / "pipeline_v2/gappa/gg2_2022.10_tree_taxonomy.tsv"
HOMD_FASTA = Path("/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03_full.fasta")
OUT = PIPE_V3 / "12h_genus_density"
OUT.mkdir(exist_ok=True)
log = lambda s: print(f"[12h_prep] {s}", flush=True)

# ---------- 1. taxonomy ----------
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

# ---------- 2. HOMD catalog ----------
homd_species, homd_epithets, homd_genus = set(), {}, set()
with open(HOMD_FASTA) as f:
    for line in f:
        if line.startswith(">"):
            parts = line[1:].split()
            if len(parts) >= 3:
                homd_genus.add(parts[2].strip())
                if len(parts) >= 4:
                    sp = f"{parts[2].strip()} {parts[3].strip()}"
                    homd_species.add(sp)
                    homd_epithets.setdefault(parts[3].strip(), set()).add(sp)
log(f"HOMD species: {len(homd_species)} | genera: {len(homd_genus)}")

def renamed_in_homd(sp, homd_epithets):
    if not sp:
        return False
    parts = sp.split(" ")
    if len(parts) < 2:
        return False
    return len(homd_epithets.get(parts[-1].strip(), ())) == 1

# ---------- 3. load tree ----------
log("Loading tree ...")
t0 = time.time()
tree = Tree(str(TREE), format=1, quoted_node_names=True)
name2node = {leaf.name: leaf for leaf in tree.iter_leaves()}
log(f"  loaded {len(name2node):,} tips in {time.time()-t0:.0f}s")

# ---------- 4. per-tip genus / homd ----------
tip_genus = {}; tip_phylum = {}; tip_homd = {}
n_asv = 0
for name, leaf in name2node.items():
    if name.startswith("DEREP_"):
        n_asv += 1
        continue
    lin = tax.get(name)
    tip_genus[name] = get_level(lin, "g") if lin else None
    tip_phylum[name] = get_level(lin, "p") if lin else None
    sp = get_level(lin, "s") if lin else None
    tip_homd[name] = bool(sp and (sp in homd_species or renamed_in_homd(sp, homd_epithets)))
log(f"ref tips: {len(tip_genus):,} | ASV tips: {n_asv:,} | HOMD ref tips: {sum(tip_homd.values()):,}")

# ---------- 5. genus MRCA nodes (fast) ----------
genus_tips = defaultdict(list)
for name, g in tip_genus.items():
    if g:
        genus_tips[g].append(name)
log(f"genera with refs: {len(genus_tips):,}")

genus_node = {}
t0 = time.time()
for g, tips in genus_tips.items():
    if len(tips) == 1:
        genus_node[g] = name2node[tips[0]]
    else:
        genus_node[g] = tree.get_common_ancestor([name2node[t] for t in tips])
log(f"  genus MRCAs computed in {time.time()-t0:.0f}s")

genus_node_ids = {id(nd): g for g, nd in genus_node.items()}

# ---------- 6. assign ASVs to deepest genus node ----------
asv_genus = {}
asv_unclassified = []
t0 = time.time()
for name, leaf in name2node.items():
    if not name.startswith("DEREP_"):
        continue
    node = leaf
    found = None
    while node is not None:
        g = genus_node_ids.get(id(node))
        if g is not None:
            found = g
            break
        node = node.up
    if found is not None:
        asv_genus[name] = found
    else:
        asv_unclassified.append(name)
log(f"  ASV->genus assigned in {time.time()-t0:.0f}s | unclassified: {len(asv_unclassified):,}")

# ---------- 7. aggregate per genus ----------
genus_stats = []
for g, node in genus_node.items():
    n_ref = 0; n_homd = 0
    for leaf in node.iter_leaves():
        if leaf.name.startswith("DEREP_"):
            continue
        n_ref += 1
        if tip_homd.get(leaf.name):
            n_homd += 1
    n_asv_g = sum(1 for a, ag in asv_genus.items() if ag == g)
    genus_stats.append({
        "genus": g,
        "phylum": next((tip_phylum[t] for t in genus_tips[g] if tip_phylum.get(t)), None),
        "n_ref": n_ref, "n_homd_ref": n_homd, "n_asv": n_asv_g,
        "n_backbone": n_ref - n_homd,
    })
gs = pd.DataFrame(genus_stats)
gs["asv_density"] = gs["n_asv"] / (gs["n_ref"] + 1)
gs["log10_density"] = gs["asv_density"].apply(lambda x: math.log10(x + 1))
gs = gs.sort_values("n_asv", ascending=False)
gs.to_csv(OUT / "genus_stats.tsv", sep="\t", index=False)
log(f"genus stats: {len(gs):,} genera | ASV in genera: {gs['n_asv'].sum():,} | unclassified ASV: {len(asv_unclassified):,}")
log(f"top 10 genera by ASV count:\n{gs.head(10)[['genus','phylum','n_ref','n_homd_ref','n_asv','asv_density']].to_string(index=False)}")

# ---------- 8. collapse genus clades into single tips ----------
collapsed = set()
for g, node in sorted(genus_node.items(), key=lambda kv: kv[1].get_distance(tree), reverse=True):
    anc = node.up
    inside = False
    while anc is not None:
        if id(anc) in collapsed:
            inside = True
            break
        anc = anc.up
    if inside:
        continue
    node.name = g
    for ch in list(node.children):
        node.remove_child(ch)
    collapsed.add(id(node))

# then collapse cold clades that hold neither an ASV nor an HOMD reference (to cut rendering cost)
def collapse_cold(root):
    changed = True
    while changed:
        changed = False
        for n in list(root.traverse()):
            if n.is_leaf() or n is root:
                continue
            has_asv = any(l.name.startswith("DEREP_") for l in n.iter_leaves())
            has_homd = any(tip_homd.get(l.name) for l in n.iter_leaves())
            if not has_asv and not has_homd:
                n.name = "cold"
                for ch in list(n.children):
                    n.remove_child(ch)
                changed = True
collapse_cold(tree)
log(f"collapsed tree tips: {len(tree.get_leaf_names()):,}")

tree.write(outfile=str(OUT / "genus_tree.nwk"), format=1, quoted_node_names=True)
log("saved genus_tree.nwk")
