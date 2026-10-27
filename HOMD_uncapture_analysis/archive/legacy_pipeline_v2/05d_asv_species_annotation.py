#!/usr/bin/env python3
"""
pipeline_v2 Step 05d: Per-ASV species / rank annotation by tree distance.

For EACH placed ASV, measure the genetic distance (substitutions per site on the
GG2 2022.10 16S-based reference tree) to its nearest reference tip, then annotate:

  species  : D <= 0.03 AND the nearest tip has a real (non-placeholder) binomial
             -> that species  (0.03 ~ 97% identity, the conventional species boundary)
             otherwise -> "novel" (candidate new species)
  genus    : D <= 0.05 (~95%)  -> nearest tip's genus
  family   : D <= 0.10 (~90%)
  order    : D <= 0.15 (~85%)
  class    : D <= 0.20 (~80%)
  phylum   : D <= 0.25 (~75%)
  else     : domain ("Bacteria" / "Archaea") only

D is computed from the SEPP placement geometry:
  insertion point P on edge A->B at `distal_length` dl from A, branch length L:
    D = pendant + min( (L - dl) + sub_min[B],   # nearest tip in B's subtree
                       dl        + out_min[B] ) # nearest tip outside B's subtree
where sub_min/out_min are precomputed per-node nearest-tip distances
(two-pass post-order + reroot).

Outputs (pipeline_v2/):
  - 05d_asv_species_annotation.tsv : per-ASV annotation (all placed ASVs)
"""
import json
import math
import re
import sys
from pathlib import Path

import pandas as pd

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"
GAPPA = PIPE_V2 / "gappa"
log = lambda s: print(f"[05d_asv_ann] {s}", flush=True)
sys.setrecursionlimit(1000000)

# conventional 16S p-distance thresholds for the rank cascade
RANK_CUTOFFS = [("species", 0.03), ("genus", 0.05), ("family", 0.10),
                ("order", 0.15), ("class", 0.20), ("phylum", 0.25)]
GENUS_PLACEHOLDER = re.compile(r"^(UBA|UMGS|GCA|GWA|CAG|SURF|JAAX|UQr|sp)\d")


def clean_species(s):
    """GG2 s__ ('Treponema_B_986142 vincentii') → 'Genus species', or None.

    Placeholder species (spNNN, e.g. 'Dialister sp001553355') are NOT names →
    return None (treated as unresolved)."""
    if not isinstance(s, str) or not s:
        return None
    parts = s.split(" ")
    if len(parts) < 2:
        return None
    genus = parts[0].split("_")[0]
    sp = " ".join(parts[1:]).strip()
    if re.fullmatch(r"sp\d+", sp):
        return None
    return f"{genus} {sp}".strip()


def clean_genus(s):
    """GG2 g__ ('Leptotrichia_A_993758') → 'Leptotrichia' (strip GTDB strain token)."""
    if not isinstance(s, str) or not s:
        return None
    return s.split("_")[0]


# ============================================================
# 1. Parse the jplace reference tree, KEEPING branch lengths
#    and edge labels. Unlike 05c, named internal nodes
#    ( ')...):NAME:len[edge]' ) are handled so deep edges are kept.
# ============================================================
log("Parsing jplace reference tree (keeping branch lengths)...")
with open(BASE / "pipeline/05_sepp/placements_export/placements.json") as f:
    data = json.load(f)
s = data["tree"]
L = len(s)

# per-node arrays: par, kids, edge label, tip name, incoming branch length
par = []
kids = []
edge = []
tip = []
elen = []


def mk(par_i, edg=None, t=None, ln=None):
    idx = len(par)
    par.append(par_i)
    kids.append([])
    edge.append(edg)
    tip.append(t)
    elen.append(ln)
    if par_i is not None:
        kids[par_i].append(idx)
    return idx


i = 0
root = mk(None)
stack = [root]
while i < L:
    c = s[i]
    if c == "(":
        stack.append(mk(stack[-1]))
        i += 1
    elif c == ")":
        closed = stack.pop()
        i += 1
        # optional internal node name (e.g. UQrYOlnDN...__Bacteria)
        if i < L and s[i] not in ":[],();":
            while i < L and s[i] not in ":[],();":
                i += 1
        if i < L and s[i] == ":":
            j = i + 1
            while j < L and s[j] not in "[],();":
                j += 1
            elen[closed] = float(s[i + 1:j])
            i = j
        if i < L and s[i] == "[":
            k = s.index("]", i)
            edge[closed] = int(s[i + 1:k])
            i = k + 1
    elif c == ",":
        i += 1
    elif c == ";":
        break
    else:
        j = i
        while j < L and s[j] not in ":,()[];":
            j += 1
        nd = mk(stack[-1], t=s[i:j])
        i = j
        if i < L and s[i] == ":":
            j = i + 1
            while j < L and s[j] not in "[],();":
                j += 1
            elen[nd] = float(s[i + 1:j])
            i = j
        if i < L and s[i] == "[":
            k = s.index("]", i)
            edge[nd] = int(s[i + 1:k])
            i = k + 1

# The top-of-tree wrapper node (between the virtual root and the Bacteria node)
# is closed with no ':len[edge]' annotation → its edge length is undefined. It
# carries no edge label, so no placement ever lands on it; treating its missing
# length as 0.0 is exact (distances through it are unchanged).
elen = [0.0 if x is None else x for x in elen]

n = len(par)
n_tips = sum(1 for t in tip if t is not None)
log(f"  nodes: {n:,} | tips: {n_tips:,}")

edge_to_node = {}
for idx, e in enumerate(edge):
    if e is not None:
        edge_to_node[e] = idx
log(f"  edges with node: {len(edge_to_node):,}")

# ============================================================
# 2. Per-node nearest-tip distances (sub_min post-order,
#    out_min pre-order / reroot).  (dist, tip) pairs; inf = none.
# ============================================================
INF = math.inf
sub_d = [INF] * n
sub_t = [None] * n
out_d = [INF] * n
out_t = [None] * n

order = []
st = [root]
while st:
    u = st.pop()
    order.append(u)
    st.extend(kids[u])
# order is parent-before-children; reversed gives children-before-parent

for u in reversed(order):
    if tip[u] is not None:
        sub_d[u], sub_t[u] = 0.0, tip[u]
    else:
        bd, bt = INF, None
        for c in kids[u]:
            d = elen[c] + sub_d[c]
            if d < bd:
                bd, bt = d, sub_t[c]
        sub_d[u], sub_t[u] = bd, bt

# out_min[root] stays INF (nothing outside the root's subtree)
for u in order:
    k = kids[u]
    if not k:
        continue
    # best / second-best child contribution for the sibling term
    if len(k) == 1:
        c0 = k[0]
        best_d, best_c = INF, None
        sec_d = INF
    else:
        contribs = [(elen[c] + sub_d[c], c) for c in k]
        contribs.sort(key=lambda x: x[0])  # ascending
        best_d, best_c = contribs[0]
        sec_d = contribs[1][0]
    for c in k:
        sib_d = best_d if best_c != c else sec_d
        cand_d = out_d[u] if out_d[u] < sib_d else sib_d
        cand_t = out_t[u] if out_d[u] <= sib_d else None
        if out_d[u] == INF and sib_d == INF:
            cand_d, cand_t = INF, None
        elif out_d[u] < sib_d:
            cand_d, cand_t = out_d[u], out_t[u]
        else:
            cand_d, cand_t = sib_d, None  # tip resolved below
            # find which sibling contributed sib_d
            for s2 in k:
                if s2 == c:
                    continue
                if abs((elen[s2] + sub_d[s2]) - sib_d) < 1e-12:
                    cand_t = sub_t[s2]
                    break
        if cand_d != INF:
            out_d[c] = elen[c] + cand_d
            out_t[c] = cand_t
log("  sub_min / out_min computed")

# ============================================================
# 3. Taxonomy (GG2 2022.10, tree-tip IDs) → full lineage
# ============================================================
log("Loading GG2 2022.10 taxonomy...")
RANKS = {"d": "domain", "p": "phylum", "c": "class", "o": "order",
         "f": "family", "g": "genus", "s": "species"}
tax = {}
with open(GAPPA / "gg2_2022.10_tree_taxonomy.tsv") as f:
    for line in f:
        parts = line.rstrip().split("\t")
        if len(parts) < 2:
            continue
        lineage = {}
        for tok in parts[1].split(";"):
            tok = tok.strip()
            if len(tok) >= 4 and tok[1:3] == "__" and tok[3:]:
                lineage[RANKS.get(tok[0])] = tok[3:]
        tax[parts[0]] = lineage
log(f"  tips with taxonomy: {len(tax):,}")


def resolve_rank(D, lineage):
    """(rank, label) — finest reliable label for distance D on the tree."""
    sp_clean = clean_species(lineage.get("species"))
    genus = clean_genus(lineage.get("genus"))
    for rank, cutoff in RANK_CUTOFFS:
        if D > cutoff:
            continue
        if rank == "species":
            if sp_clean:
                return rank, sp_clean
            continue  # species within boundary but no real name -> coarser
        if rank == "genus":
            if genus:
                return rank, genus
            continue
        if lineage.get(rank):
            return rank, lineage[rank]
    if lineage.get("domain"):
        return "domain", lineage["domain"]
    return "unclassified", None


# ============================================================
# 4. Per-ASV annotation
# ============================================================
log("Loading placed ASVs (step 05 output)...")
asv = pd.read_csv(PIPE_V2 / "05_asv_phylotype_assignments.tsv", sep="\t", low_memory=False)
log(f"  ASVs with placement: {len(asv):,}")

rows = []
n_edge_missing = 0
n_no_tip = 0
for _, r in asv.iterrows():
    edge_num = r["edge_num"]
    if pd.isna(edge_num) or int(edge_num) not in edge_to_node:
        n_edge_missing += 1
        rows.append({"derep_id": r["derep_id"], "phylotype_id": r.get("phylotype_id", ""),
                     "edge_num": None, "distal_length": None, "pendant_length": None,
                     "nearest_tip": None, "dist_to_species": None,
                     "species_annotation": None, "species_novel": None,
                     "rank_annotation": "unknown", "rank_label": None})
        continue
    B = edge_to_node[int(edge_num)]
    L_len = elen[B] if elen[B] is not None else INF
    dl = r["distal_length"]
    p = r["pendant_length"]
    if pd.isna(dl) or pd.isna(p) or L_len == INF:
        n_edge_missing += 1
        rows.append({"derep_id": r["derep_id"], "phylotype_id": r.get("phylotype_id", ""),
                     "edge_num": edge_num, "distal_length": dl, "pendant_length": p,
                     "nearest_tip": None, "dist_to_species": None,
                     "species_annotation": None, "species_novel": None,
                     "rank_annotation": "unknown", "rank_label": None})
        continue
    d_child = max(0.0, L_len - dl) + sub_d[B]
    d_out = dl + out_d[B]
    if d_child <= d_out:
        D_to_tip, nt = d_child, sub_t[B]
    else:
        D_to_tip, nt = d_out, out_t[B]
    if nt is None or D_to_tip == INF:
        n_no_tip += 1
        rows.append({"derep_id": r["derep_id"], "phylotype_id": r.get("phylotype_id", ""),
                     "edge_num": edge_num, "distal_length": dl, "pendant_length": p,
                     "nearest_tip": None, "dist_to_species": None,
                     "species_annotation": None, "species_novel": None,
                     "rank_annotation": "unknown", "rank_label": None})
        continue
    D = p + D_to_tip
    lineage = tax.get(nt, {})
    nearest_species = lineage.get("species")
    sp_clean = clean_species(nearest_species)
    if D <= 0.03 and sp_clean:
        sp_ann, sp_novel = sp_clean, False
    else:
        sp_ann, sp_novel = "novel", True
    rank_ann, rank_label = resolve_rank(D, lineage)
    rows.append({"derep_id": r["derep_id"], "phylotype_id": r.get("phylotype_id", ""),
                 "edge_num": edge_num, "distal_length": dl, "pendant_length": p,
                 "nearest_tip": nt, "dist_to_species": round(D, 6),
                 "nearest_genus": clean_genus(lineage.get("genus")),
                 "nearest_species": nearest_species,
                 "species_annotation": sp_ann, "species_novel": sp_novel,
                 "rank_annotation": rank_ann, "rank_label": rank_label})

ann = pd.DataFrame(rows)
ann.to_csv(PIPE_V2 / "05d_asv_species_annotation.tsv", sep="\t", index=False)
log(f"  saved 05d_asv_species_annotation.tsv ({len(ann):,} rows)")
log(f"  edges missing: {n_edge_missing} | no reachable tip: {n_no_tip}")

ok = ann["species_annotation"].notna()
sp_hit = ok & (ann["species_novel"] == False)
log(f"  annotated to species (D<=0.03, real binomial): {sp_hit.sum():,} "
    f"({sp_hit.mean()*100:.1f}%)")
log(f"  novel / not species: {(ok & ann['species_novel']).sum():,}")
log("\nD distribution (annotated ASVs):")
print(ann.loc[ok, "dist_to_species"].describe(percentiles=[.5, .75, .9, .99]).to_string())
log("\nrank_annotation counts:")
print(ann["rank_annotation"].value_counts().to_string())
log("\nDone!")
