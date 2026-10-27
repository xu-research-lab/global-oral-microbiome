#!/usr/bin/env python3
"""
pipeline_v3 Step 06: Build ALL-candidate collapsed tree (v2 09b equivalent).

Candidate novel clusters span most of the 385K-tip SEPP insertion tree, so a
"local context" subtree is impossible. Instead we collapse: every clade that
contains no candidate member tip is replaced by a single tip `CLADE_<n>`
(a triangle in the R figure, n = number of collapsed reference tips). All
candidate rep ASV tips are retained, so the full candidate distribution stays
visible in real phylogenetic position.

Outputs (into pipeline_v3/analysis_report/):
  - tree_allcandidates_collapsed.nwk
  - derep_to_cluster_allcand.tsv : candidate rep tip -> cluster + annotation
"""
import os
import time
from pathlib import Path

import pandas as pd
from ete3 import Tree

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V3 = BASE / "pipeline_v3"
REPORT = PIPE_V3 / "analysis_report"
REPORT.mkdir(parents=True, exist_ok=True)

FULL_TREE = BASE / "pipeline/05_sepp/tree_export/tree.nwk"
CAND = PIPE_V3 / "02_clusters_aggregated.tsv"

log = lambda s: print(f"[06_tree] {s}", flush=True)


def main():
    log("Loading candidate novel clusters ...")
    cn = pd.read_csv(CAND, sep="\t", low_memory=False)
    cn = cn[cn["candidate_novel"] == True]
    log(f"  candidate novel clusters: {len(cn):,}")

    candidate_tips = set()
    tip_info = {}
    for _, row in cn.iterrows():
        rep = str(row["rep_derep_id"]).strip()
        if rep:
            candidate_tips.add(rep)
            tip_info[rep] = {
                "cluster_id": row["cluster_id"],
                "annotation": row["annotation"],
                "annotation_level": row["annotation_level"],
                "n_members": row["n_asvs"],
                "subjects": row["max_subject_prevalence"],
                "bps": row["n_bioprojects_total"],
                "continents": row["n_continents"],
                "total_reads": row["total_reads_sum"],
                "member_derep_ids": row["member_derep_ids"],
            }
    log(f"  candidate rep tips: {len(candidate_tips):,}")

    t0 = time.time()
    log(f"Loading full tree ({FULL_TREE})...")
    t = Tree(str(FULL_TREE), format=1, quoted_node_names=True)
    log(f"  {len(t):,} tips in {time.time()-t0:.0f}s")

    log("Computing per-node candidate counts...")
    n_cand = {}
    n_leaves = {}
    for node in t.traverse("postorder"):
        nid = id(node)
        if node.is_leaf():
            n_cand[nid] = 1 if node.name in candidate_tips else 0
            n_leaves[nid] = 1
        else:
            n_cand[nid] = sum(n_cand[id(c)] for c in node.children)
            n_leaves[nid] = sum(n_leaves[id(c)] for c in node.children)

    found = n_cand[id(t)]
    log(f"  candidate tips found in tree: {found:,}/{len(candidate_tips):,}")

    log("Emitting collapsed newick...")
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

    rows = [{"tip": did, **info} for did, info in tip_info.items()]
    d2p = pd.DataFrame(rows)
    d2p_out = REPORT / "derep_to_cluster_allcand.tsv"
    d2p.to_csv(d2p_out, sep="\t", index=False)
    log(f"  saved: {d2p_out} ({len(d2p):,} rows)")

    log(f"\nelapsed {time.time()-t0:.0f}s - Done!")


if __name__ == "__main__":
    main()
