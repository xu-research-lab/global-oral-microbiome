#!/usr/bin/env python3
"""
pipeline_v2 Step 09: Build phylotype subtree for ggtree (FAST).

Two phases:
  1. ete3: find the tip set to keep — phylotype member DEREP tips +
     sister GG2 reference tips (grandparent context). Single pass, fast.
  2. Streaming newick parser: single pass over the 20 MB newick string,
     dropping any subtree that contains no kept tip. Builds the pruned
     tree directly — no in-memory prune of the 385 K-tip tree.

Outputs (into pipeline_v2/analysis_report/):
  - phylotype_subtree_top{n}.nwk   : corrected subtree
  - derep_to_phylotype.tsv         : DEREP → phylotype + support metadata
"""
import pandas as pd
import numpy as np
import sys, os, re, warnings
from pathlib import Path
from ete3 import Tree

warnings.filterwarnings("ignore")

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"
REPORT = PIPE_V2 / "analysis_report"
REPORT.mkdir(parents=True, exist_ok=True)

FULL_TREE = BASE / "pipeline/05_sepp/tree_export/tree.nwk"
CAND_NOVEL = PIPE_V2 / "06_candidate_novel_phylotypes.tsv"
TOP_N = 150

log = lambda s: print(f"[09_tree] {s}", flush=True)

# ============================================================
# 1. Load candidate novel phylotypes, rank by total reads
# ============================================================
log("Loading candidate novel phylotypes...")
phy = pd.read_csv(CAND_NOVEL, sep="\t", low_memory=False)
for c in ["total_reads_sum", "max_prevalence", "max_n_bioprojects"]:
    phy[c] = pd.to_numeric(phy[c], errors="coerce").fillna(0)

phy = phy.sort_values("total_reads_sum", ascending=False).reset_index(drop=True)
top = phy.head(TOP_N).copy()
top["support_rank"] = range(1, len(top) + 1)
log(f"  Top-{TOP_N} phylotypes selected (by total_reads_sum)")

all_derep = set()
derep_to_info = {}
for _, row in top.iterrows():
    for did in str(row["member_derep_ids"]).split(";"):
        did = did.strip()
        if did:
            all_derep.add(did)
            derep_to_info[did] = {
                "phylotype": row["phylotype_id"],
                "novelty": row["novelty_class"],
                "n_asvs": row["n_asvs"],
                "reads": row["total_reads_sum"],
                "bps": row["max_n_bioprojects"],
                "subjects": row["max_prevalence"],
                "rank": row["support_rank"],
            }
log(f"  Member DEREP ASVs: {len(all_derep):,}")

# ============================================================
# 2. Phase 1 — ete3 single pass: find keep set
# ============================================================
log("Phase 1: finding tip set to keep (ete3, single pass)...")
t = Tree(str(FULL_TREE), format=1, quoted_node_names=True)
log(f"  Full tree tips: {len(t):,}")

keep_tips = set()
found = 0
seen_anchors = set()
for leaf in t.iter_leaves():
    if leaf.name in all_derep:
        found += 1
        parent = leaf.up
        if parent is None:
            continue
        grandparent = parent.up
        anchor = grandparent if grandparent is not None else parent
        if id(anchor) not in seen_anchors:
            seen_anchors.add(id(anchor))
            keep_tips.update(anchor.get_leaf_names())

log(f"  Found DEREP in tree: {found}/{len(all_derep)}")
log(f"  Unique context subtrees: {len(seen_anchors)}")
log(f"  Tips to keep: {len(keep_tips):,}")
del t  # free memory before streaming

missing = all_derep - keep_tips
if missing:
    log(f"  WARNING: {len(missing)} DEREPs missing from keep set")

# ============================================================
# 3. Phase 2 — streaming newick parser that prunes in one pass
# ============================================================
log("Phase 2: streaming prune (single pass over newick)...")

def read_token(s, i):
    """Read one node token starting at i: optional quoted name + optional :length.
    Returns (full_token, clean_name, new_i). full_token is the complete newick
    fragment (name + optional ':length'), clean_name is the name without quotes
    (used for the keep-set test)."""
    n = len(s)
    if s[i] == "'":
        j = i + 1
        while j < n and s[j] != "'":
            j += 1
        clean = s[i + 1:j]
        j += 1
        start = j
        while j < n and s[j] not in ",);":
            j += 1
        return s[start:j], clean, j
    else:
        j = i
        while j < n and s[j] not in "():,;":
            j += 1
        clean = s[i:j]
        start = i
        if j < n and s[j] == ":":
            j += 1
            while j < n and s[j] not in ",);":
                j += 1
            return s[start:j], clean, j
        return s[start:j], clean, j

def parse(s, i, keep):
    """Parse subtree at i; return (newick_fragment_or_None, new_i)."""
    n = len(s)
    if s[i] == "(":
        i += 1
        frags = []
        while True:
            frag, i = parse(s, i, keep)
            if frag is not None:
                frags.append(frag)
            c = s[i]
            if c == ",":
                i += 1
            elif c == ")":
                i += 1
                break
        token, clean, i = read_token(s, i)
        if frags:
            return "(" + ",".join(frags) + ")" + token, i
        return None, i
    else:
        token, clean, i = read_token(s, i)
        if clean in keep:
            return token if token else clean, i
        return None, i

with open(FULL_TREE) as f:
    content = f.read().strip()

out_frag, _ = parse(content, 0, keep_tips)
if out_frag is None:
    log("  ERROR: nothing kept")
    sys.exit(1)

# Normalize: ensure trailing semicolon
pruned = out_frag.rstrip(";") + ";"

out_nwk = REPORT / f"phylotype_subtree_top{TOP_N}.nwk"
with open(out_nwk, "w") as f:
    f.write(pruned + "\n")
log(f"  Saved: {out_nwk} ({os.path.getsize(out_nwk):,} bytes)")

# ============================================================
# 4. Save derep → phylotype mapping
# ============================================================
rows = [{"tip": did, **info} for did, info in derep_to_info.items()]
d2p = pd.DataFrame(rows)
d2p_out = REPORT / "derep_to_phylotype.tsv"
d2p.to_csv(d2p_out, sep="\t", index=False)
log(f"  Saved: {d2p_out} ({len(d2p):,} rows)")

# ============================================================
# Summary
# ============================================================
log("\n=== Top-10 phylotypes ===")
for _, r in top.head(10).iterrows():
    log(f"  {r['phylotype_id']} | {r['novelty_class']} | ASVs={r['n_asvs']} | "
        f"BPs={r['max_n_bioprojects']:.0f} | subj={r['max_prevalence']:.0f} | "
        f"reads={r['total_reads_sum']:,.0f}")

log("\nDone!")
