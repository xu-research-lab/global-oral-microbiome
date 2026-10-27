#!/usr/bin/env python3
"""
pipeline_v2 Step 09: Build phylotype subtree for the 276 strict candidate phylotypes.

Two phases:
  1. ete3: find the tip set to keep — strict candidate member DEREP tips +
     sister GG2 reference tips (grandparent context). Single pass, fast.
  2. Streaming newick parser: single pass over the 20 MB newick, dropping any
     subtree that contains no kept tip.

Outputs (into pipeline_v2/analysis_report/):
  - phylotype_subtree_{SUFFIX}.nwk
  - derep_to_phylotype_{SUFFIX}.tsv
"""
import pandas as pd
import sys, os
from pathlib import Path
from ete3 import Tree

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"
REPORT = PIPE_V2 / "analysis_report"
REPORT.mkdir(parents=True, exist_ok=True)

FULL_TREE = BASE / "pipeline/05_sepp/tree_export/tree.nwk"
STRICT = PIPE_V2 / "07_strict_candidates.tsv"

# optional top-N (by support_score) for a focused tree
TOP_N = int(sys.argv[1]) if len(sys.argv) > 1 else None

log = lambda s: print(f"[09_tree] {s}", flush=True)

# ============================================================
# 1. Strict candidate phylotypes + members
# ============================================================
log("Loading strict candidate phylotypes...")
strict = pd.read_csv(STRICT, sep="\t", low_memory=False)
# merge member_derep_ids from the novelty-classified table (not in aggregated table)
nov = pd.read_csv(PIPE_V2 / "06_phylotypes_novelty_classified.tsv", sep="\t", low_memory=False,
                  usecols=["phylotype_id", "member_derep_ids"])
strict = strict.merge(nov, on="phylotype_id", how="left")
if TOP_N:
    strict = strict.sort_values("support_score", ascending=False).head(TOP_N)
SUFFIX = f"top{TOP_N}" if TOP_N else f"strict{len(strict)}"
log(f"  {len(strict):,} strict phylotypes (suffix: {SUFFIX})")

all_derep = set()
derep_to_info = {}
for _, row in strict.iterrows():
    for did in str(row["member_derep_ids"]).split(";"):
        did = did.strip()
        if did:
            all_derep.add(did)
            derep_to_info[did] = {
                "phylotype": row["phylotype_id"],
                "novelty": row["novelty_class"],
                "n_asvs": row["n_members"],
                "subjects": row["max_subject_prevalence"],
                "bps": row["n_bioprojects_total"],
                "continents": row["n_continents"],
                "regions": row["n_regions"],
                "n_flagged": row["n_flagged"],
            }
log(f"  Member DEREP ASVs: {len(all_derep):,}")

# ============================================================
# 2. Phase 1 — ete3 grandparent-context keep set
# ============================================================
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
log(f"  Unique context subtrees: {len(seen_anchors)}")
log(f"  Tips to keep: {len(keep_tips):,}")
del t

missing = all_derep - keep_tips
if missing:
    log(f"  WARNING: {len(missing)} DEREPs missing from keep set")

# ============================================================
# 3. Phase 2 — streaming prune
# ============================================================
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

pruned = out_frag.rstrip(";") + ";"
out_nwk = REPORT / f"phylotype_subtree_{SUFFIX}.nwk"
with open(out_nwk, "w") as f:
    f.write(pruned + "\n")
log(f"  Saved: {out_nwk} ({os.path.getsize(out_nwk):,} bytes)")

# ============================================================
# 4. derep → phylotype mapping
# ============================================================
rows = [{"tip": did, **info} for did, info in derep_to_info.items()]
d2p = pd.DataFrame(rows)
# merge genus for the tree labels
nov_genus = pd.read_csv(PIPE_V2 / "06_phylotypes_novelty_classified.tsv", sep="\t", low_memory=False,
                        usecols=["phylotype_id", "rep_gg2_genus"])
d2p = d2p.merge(nov_genus, left_on="phylotype", right_on="phylotype_id", how="left").drop(columns=["phylotype_id"])
d2p_out = REPORT / f"derep_to_phylotype_{SUFFIX}.tsv"
d2p.to_csv(d2p_out, sep="\t", index=False)
log(f"  Saved: {d2p_out} ({len(d2p):,} rows)")

log("\nDone!")
