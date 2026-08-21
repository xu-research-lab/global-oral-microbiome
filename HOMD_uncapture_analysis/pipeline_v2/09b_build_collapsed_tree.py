#!/usr/bin/env python3
"""
pipeline_v2 Step 09b: Build ALL-candidate collapsed tree.

Candidate novel phylotypes span ~99.7% of the 385K-tip SEPP insertion tree, so a
"local context" subtree is impossible. Instead we collapse: every clade that contains
no candidate member tip is replaced by a single tip `CLADE_<n>` (a triangle in the
R figure, n = number of collapsed reference tips). All candidate member DEREP tips are
retained, so the full candidate distribution stays visible in real phylogenetic position.

Outputs (into pipeline_v2/analysis_report/):
  - tree_allcandidates_collapsed.nwk   : collapsed tree (candidate tips + CLADE_<n> tips)
  - derep_to_phylotype.tsv             : candidate tip → phylotype + annotation rings
"""
import pandas as pd
import numpy as np
import sys, os, time
from pathlib import Path
from ete3 import Tree

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"
REPORT = PIPE_V2 / "analysis_report"
REPORT.mkdir(parents=True, exist_ok=True)

FULL_TREE = BASE / "pipeline/05_sepp/tree_export/tree.nwk"
CAND = PIPE_V2 / "06_phylotypes_novelty_classified.tsv"
RES_EVAL = BASE / "pipeline/04_oral_vs_blank/abundance_long.tsv"
RES_PASS = BASE / "pipeline/04_oral_vs_blank/oral_vs_blank_enrichment.tsv"
FLAG_TSV = PIPE_V2 / "02_contamination_flagged.tsv"

log = lambda s: print(f"[09b_collapse] {s}", flush=True)

# ============================================================
# 1. Candidate phylotypes + per-member annotations
# ============================================================
log("Loading candidate phylotypes...")
cn = pd.read_csv(CAND, sep="\t", low_memory=False)
cn = cn[cn["candidate_novel"] == True].copy()
log(f"  {len(cn):,} candidate phylotypes")

# residence pass / evaluated sets
res_pass = set(pd.read_csv(RES_PASS, sep="\t")["feature_id"].astype(str))
eval_set = set(pd.read_csv(RES_EVAL, sep="\t")["feature_id"].astype(str))
# flagged set (for composition; currently all phylotype members are clean by construction)
fl = pd.read_csv(FLAG_TSV, sep="\t", low_memory=False, dtype=str)
flagged_set = set(fl.loc[fl["contamination_flag"].astype(str).str.lower().isin(["true","1","t"]), "derep_id"].astype(str))

# one tip per phylotype (the representative ASV) — avoids over-weighting multi-ASV phylotypes
candidate_tips = set()
tip_info = {}
for _, row in cn.iterrows():
    rep = str(row["rep_derep_id"]).strip()
    if rep:
        candidate_tips.add(rep)
        tip_info[rep] = {
            "phylotype_id": row["phylotype_id"],
            "novelty": row["novelty_class"],
            "n_members": row["n_asvs"],
            "subjects": row["max_prevalence"],
            "bps": row["n_bioprojects_total"],
            "continents": row["n_continents"],
            "total_reads": row["total_reads_sum"],
            "member_derep_ids": row["member_derep_ids"],
        }
# per-phylotype composition from ALL members (for annotation rings)
for rep, info in tip_info.items():
    members = [m for m in str(info["member_derep_ids"]).split(";") if m]
    info["n_flagged"] = sum(1 for m in members if m in flagged_set)
    info["n_residence_pass"] = sum(1 for m in members if m in res_pass)
    info["n_tested"] = sum(1 for m in members if m in eval_set)
    info["n_untested"] = len(members) - info["n_tested"]
    info.pop("member_derep_ids")
log(f"  candidate phylotype tips (rep ASV): {len(candidate_tips):,}")

# ============================================================
# 2. Load full tree
# ============================================================
t0 = time.time()
log(f"Loading full tree ({FULL_TREE})...")
t = Tree(str(FULL_TREE), format=1, quoted_node_names=True)
log(f"  {len(t):,} tips in {time.time()-t0:.0f}s")

# ============================================================
# 3. Per-node candidate counts (bottom-up)
# ============================================================
log("Computing per-node candidate counts...")
n_cand = {}
n_leaves = {}
for node in t.traverse("postorder"):
    nid = id(node)
    if node.is_leaf():
        n_cand[nid] = 1 if node.name in candidate_tips else 0
        n_leaves[nid] = 1
    else:
        c_cand = sum(n_cand[id(c)] for c in node.children)
        c_leaf = sum(n_leaves[id(c)] for c in node.children)
        n_cand[nid] = c_cand
        n_leaves[nid] = c_leaf

found = n_cand[id(t)]
log(f"  candidate tips found in tree: {found:,}/{len(candidate_tips):,}")

# ============================================================
# 4. Emit collapsed newick (post-order)
# ============================================================
log("Emitting collapsed newick...")

# iterative post-order emit (tree is too deep for recursion)
frag = {}
for node in t.traverse("postorder"):
    nid = id(node)
    dist = f":{node.dist:.6f}"
    if node.is_leaf():
        if node.name in candidate_tips:
            frag[nid] = node.name + dist
        else:
            frag[nid] = f"CLADE_1{dist}"
    elif n_cand[nid] == 0:
        # collapse: single tip representing n_leaves references
        frag[nid] = f"CLADE_{n_leaves[nid]}{dist}"
    else:
        kids = ",".join(frag[id(c)] for c in node.children)
        frag[nid] = f"({kids}){dist}"

nwk = frag[id(t)].rstrip(":0123456789.") + ";"
out_nwk = REPORT / "tree_allcandidates_collapsed.nwk"
with open(out_nwk, "w") as f:
    f.write(nwk + "\n")
log(f"  saved: {out_nwk} ({os.path.getsize(out_nwk)/1e6:.1f} MB)")
log(f"  CLADE_ tips: {nwk.count('CLADE_')}")

# ============================================================
# 5. Save derep → phylotype annotation table
# ============================================================
# one row per phylotype (rep tip); composition already computed per phylotype
rows = [{"tip": did, **info} for did, info in tip_info.items()]
d2p = pd.DataFrame(rows)
d2p_out = REPORT / "derep_to_phylotype_allcand.tsv"
d2p.to_csv(d2p_out, sep="\t", index=False)
log(f"  saved: {d2p_out} ({len(d2p):,} rows)")

log(f"\nelapsed {time.time()-t0:.0f}s — Done!")
