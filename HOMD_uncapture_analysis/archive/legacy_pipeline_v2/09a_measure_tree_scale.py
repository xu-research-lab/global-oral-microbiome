#!/usr/bin/env python3
"""Measure keep-set scale for plotting ALL candidate novel phylotypes on the tree.

Loads the 385K-tip SEPP insertion tree once, finds the keep-set for all candidate
member DEREP tips under two context strategies, and reports sizes so we can decide
the plotting strategy (plot-all vs collapse).
"""
import pandas as pd
import sys, time
from pathlib import Path
from ete3 import Tree

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"
FULL_TREE = BASE / "pipeline/05_sepp/tree_export/tree.nwk"
CAND = PIPE_V2 / "06_phylotypes_novelty_classified.tsv"

log = lambda s: print(f"[09a_scale] {s}", flush=True)

# 1. candidate members
cn = pd.read_csv(CAND, sep="\t", low_memory=False)
cn = cn[cn["candidate_novel"] == True]
members = set()
for m in cn["member_derep_ids"].dropna():
    members.update(x for x in str(m).split(";") if x)
log(f"candidate phylotypes: {len(cn):,} | member DEREP ASVs: {len(members):,}")

# 2. load tree
t0 = time.time()
log(f"loading tree: {FULL_TREE}")
t = Tree(str(FULL_TREE), format=1, quoted_node_names=True)
log(f"  loaded {len(t):,} tips in {time.time()-t0:.0f}s")

found = 0
keep_grand = set()   # grandparent context (grandparent's full leaf set)
keep_parent = set()  # parent context (parent's full leaf set)
seen_g = set(); seen_p = set()
for leaf in t.iter_leaves():
    if leaf.name in members:
        found += 1
        par = leaf.up
        if par is not None:
            gpar = par.up
            # parent context
            if id(par) not in seen_p:
                seen_p.add(id(par)); keep_parent.update(par.get_leaf_names())
            # grandparent context
            if gpar is not None:
                if id(gpar) not in seen_g:
                    seen_g.add(id(gpar)); keep_grand.update(gpar.get_leaf_names())
            else:
                if id(par) not in seen_g:
                    seen_g.add(id(par)); keep_grand.update(par.get_leaf_names())

log(f"found DEREP in tree: {found}/{len(members)}")
log(f"parent-context keep-set: {len(keep_parent):,} tips")
log(f"grandparent-context keep-set: {len(keep_grand):,} tips")
log(f"  of which DEREP (query): {len(keep_grand & members):,} | GG2 refs: {len(keep_grand - members):,}")
log(f"unique grandparent context nodes: {len(seen_g):,}")
log(f"elapsed: {time.time()-t0:.0f}s")
