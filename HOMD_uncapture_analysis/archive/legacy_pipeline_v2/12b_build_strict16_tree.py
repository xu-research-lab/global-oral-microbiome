#!/usr/bin/env python3
"""
pipeline_v2 Step 12b: build the phylotype subtree for the 16 STRICT candidates
(new scheme: GG2 2024.09 annotation, edge x label phylotypes, step 11).

Same algorithm as 09_build_strict_tree.py:
  1. ete3: keep DEREP members + sister GG2 reference tips (grandparent context);
  2. streaming newick parser over the 20 MB SEPP tree.

Outputs (analysis_report/):
  - phylotype_subtree_strict16.nwk
  - derep_to_phylotype_strict16.tsv
"""
import os
import sys
from pathlib import Path

import pandas as pd
from ete3 import Tree

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"
REPORT = PIPE_V2 / "analysis_report"
FULL_TREE = BASE / "pipeline/05_sepp/tree_export/tree.nwk"
STRICT = Path(sys.argv[1]) if len(sys.argv) > 1 else PIPE_V2 / "11_headline_candidates.tsv"

log = lambda s: print(f"[12b_tree] {s}", flush=True)

log("Loading strict candidates (step 11)...")
strict = pd.read_csv(STRICT, sep="\t", low_memory=False)
log(f"  {len(strict)} strict phylotypes")
PREFIX = sys.argv[2] if len(sys.argv) > 2 else "headline"
SUFFIX = f"{PREFIX}{len(strict)}"

all_derep = set()
derep_to_info = {}
for _, row in strict.iterrows():
    for did in str(row["member_derep_ids"]).split(";"):
        did = did.strip()
        if did:
            all_derep.add(did)
            derep_to_info[did] = {
                "phylotype": row["phylotype_id"],
                "novelty": "homd_uncaptured",
                "n_asvs": row["n_asvs"],
                "subjects": row["max_subject_prevalence"],
                "bps": row["n_bioprojects_total"],
                "regions": row["n_regions"],
                "rep_species": row.get("rep_species", ""),
                "rep_genus": row.get("rep_genus", ""),
                "identity": row.get("rep_identity", ""),
            }
log(f"  Member DEREP ASVs: {len(all_derep):,}")

# ---------------- Phase 1: ete3 keep set ----------------
log("Phase 1: finding keep-set (ete3)...")
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
log(f"  Unique context subtrees: {len(seen_anchors)} | tips to keep: {len(keep_tips):,}")
del t

# ---------------- Phase 2: streaming prune ----------------
log("Phase 2: streaming prune...")


def read_token(s, i):
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
    j = i
    while j < n and s[j] not in "():,;":
        j += 1
    clean = s[i:j]
    if j < n and s[j] == ":":
        j += 1
        while j < n and s[j] not in ",);":
            j += 1
        return s[i:j], clean, j
    return s[i:j], clean, j


def parse(s, i, keep):
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
        token, _, i = read_token(s, i)
        if frags:
            return "(" + ",".join(frags) + ")" + token, i
        return None, i
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

pruned = out_frag.rstrip(";") + ";"
out_nwk = REPORT / f"phylotype_subtree_{SUFFIX}.nwk"
with open(out_nwk, "w") as f:
    f.write(pruned + "\n")
log(f"  Saved: {out_nwk} ({os.path.getsize(out_nwk):,} bytes)")

# ---------------- derep -> phylotype mapping ----------------
rows = [{"tip": did, **info} for did, info in derep_to_info.items()]
d2p = pd.DataFrame(rows)
out_tsv = REPORT / f"derep_to_phylotype_{SUFFIX}.tsv"
d2p.to_csv(out_tsv, sep="\t", index=False)
log(f"  Saved: {out_tsv} ({len(d2p):,} rows)")
log("Done!")
