#!/usr/bin/env python3
"""
pipeline_v2 Step 05e: Phylotype re-clustering at the 0.03 species boundary.

Replaces the 05 within-edge split (|distal_i-distal_j|+|pend_i-pend_j| <= 0.05,
which under-estimates divergence) with a species-boundary-consistent scheme:

  species-level ASVs (D <= 0.03 & real binomial, from 05d):
    phylotype key = the raw GTDB species token (e.g. 'Treponema_D succinifaciens')
    -> one phylotype per species, MERGED ACROSS EDGES
       (fixes the old behaviour that split one species over many edges;
        e.g. Anaerostipes hadrus spans 65 edges).
  novel ASVs (D > 0.03):
    within-edge single-linkage at corrected metric
       d = |distal_i - distal_j| + pend_i + pend_j  <= 0.03
       (BOTH pendants included -> consistent with the species boundary;
        the old 05 metric dropped pend_i+pend_j).

phylotype_id scheme:
  species: raw token spaces -> underscores (e.g. 'Treponema_D_succinifaciens')
  novel:   'NOVEL_EDGE_<edge>' / 'NOVEL_EDGE_<edge>_<c>' (split components)

Representative ASV (rep) = member with max subject_prevalence
(tie-break: total_reads), NOT .first() as in 05.

Outputs (overwrites the 05 outputs so the downstream chain picks up the new IDs):
  - 05_asv_phylotype_assignments.tsv  (phylotype_id column replaced, + anchor cols)
  - 05_phylotypes_all.tsv             (05 aggregation + species anchor columns)
  - 05e_branch_check.tsv              (QC: multiple phylotypes per branch)
  - 05e_branch_over_split.tsv         (flagged pairs with rep distance <= 0.03)
  - 05e_summary.json

The branch check answers "does one branch produce multiple phylotypes?":
  * edge-level: every placement edge with >=2 phylotypes -> pairwise rep distance
    within the edge  = |dl_i - dl_j| + pend_i + pend_j  (exact same-edge tree dist)
  * genus-level (cross-edge): every pair of SPECIES phylotypes sharing a clean
    genus -> full tree distance between the two rep insertion points (via LCA on
    the GG2 tree). Flags strain-split species that 05d keeps apart
    (e.g. Fusobacterium nucleatum_A/_B/_C), i.e. real over-splitting.
"""
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"
GAPPA = PIPE_V2 / "gappa"
log = lambda s: print(f"[05e_phylo] {s}", flush=True)
sys.setrecursionlimit(1000000)

SPLIT_T = 0.03          # species boundary
SPECIES_COLS = ["derep_id", "dist_to_species", "nearest_genus", "nearest_species",
                "species_annotation", "species_novel", "rank_annotation", "rank_label"]


# ============================================================
# 1b. GG2 reference tree: parser + node distances + binary-lifting LCA.
#     Reuses the 05d jplace parser (edge label -> child node, branch lengths).
#     Needed only for the cross-edge branch check (species pairs in a genus).
# ============================================================
def _build_tree():
    log("Parsing jplace reference tree (for branch check LCA)...")
    with open(BASE / "pipeline/05_sepp/placements_export/placements.json") as f:
        data = json.load(f)
    s = data["tree"]
    L = len(s)

    par, kids, edge, tip, elen = [], [], [], [], []

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

    elen = [0.0 if x is None else x for x in elen]
    n = len(par)
    par[0] = 0                      # root's parent = itself (for LCA lifting)
    edge_to_node = {e: idx for idx, e in enumerate(edge) if e is not None}

    # depth[] (branch-length) and level[] (edge count): children always have a
    # higher index than their parent (stack construction order), so a single
    # pass in index order is topological.
    depth = [0.0] * n
    level = [0] * n
    for u in range(1, n):
        depth[u] = depth[par[u]] + elen[u]
        level[u] = level[par[u]] + 1

    # binary lifting (up[k][v] = 2^k-th ancestor; root is its own ancestor)
    LOG = n.bit_length()
    up = np.zeros((LOG, n), dtype=np.int32)
    up[0] = par
    for k in range(1, LOG):
        up[k] = up[k - 1][up[k - 1]]

    def lca(u, v):
        if level[u] < level[v]:
            u, v = v, u
        diff = level[u] - level[v]
        bit = 0
        while diff:
            if diff & 1:
                u = up[bit][u]
            diff >>= 1
            bit += 1
        if u == v:
            return u
        for k in range(LOG - 1, -1, -1):
            if up[k][u] != up[k][v]:
                u = up[k][u]
                v = up[k][v]
        return up[0][u]

    def node_dist(u, v):
        if u == v:
            return 0.0
        w = lca(u, v)
        return depth[u] + depth[v] - 2.0 * depth[w]

    def point_distance(self, e1, dl1, pend1, e2, dl2, pend2):
        """Full tree distance between two ASV tips inserted at (edge, dl) + pend."""
        if e1 == e2:
            return pend1 + abs(dl1 - dl2) + pend2
        B1, B2 = edge_to_node[e1], edge_to_node[e2]
        A1, A2 = par[B1], par[B2]
        L1, L2 = elen[B1], elen[B2]
        d_insert = min(dl1 + node_dist(A1, A2) + dl2,
                       dl1 + node_dist(A1, B2) + (L2 - dl2),
                       (L1 - dl1) + node_dist(B1, A2) + dl2,
                       (L1 - dl1) + node_dist(B1, B2) + (L2 - dl2))
        return pend1 + d_insert + pend2

    log(f"  nodes: {n:,} | edges labeled: {len(edge_to_node):,}")
    return type("Tree", (), {"point_distance": point_distance})()


# ============================================================
# 2. Load & merge ASV-level data
# ============================================================
log("Loading ASV-level data...")
asv = pd.read_csv(PIPE_V2 / "05_asv_phylotype_assignments.tsv", sep="\t", low_memory=False)
log(f"  placements (from 05): {len(asv):,}")

passed = pd.read_csv(PIPE_V2 / "04_support_floor_passed.tsv", sep="\t", low_memory=False,
                     usecols=["derep_id", "prevalence_samples", "bioprojects",
                              "amplicon_regions", "continents"])
df = asv.merge(passed, on="derep_id", how="left")

subj = pd.read_csv(PIPE_V2 / "01b_derep_subject_prevalence.tsv", sep="\t", low_memory=False,
                   usecols=["derep_id", "subject_prevalence"])
df = df.merge(subj, on="derep_id", how="left")
df["subject_prevalence"] = pd.to_numeric(df["subject_prevalence"], errors="coerce").fillna(0).astype(int)

ann = pd.read_csv(PIPE_V2 / "05d_asv_species_annotation.tsv", sep="\t", low_memory=False,
                  usecols=SPECIES_COLS)
df = df.merge(ann, on="derep_id", how="left")
log(f"  merged: {len(df):,} ASVs | species-level: {(df['species_novel'] == False).sum():,} | "
    f"novel: {(df['species_novel'] != False).sum():,}")

# ============================================================
# 2. Clustering at the 0.03 species boundary
# ============================================================
df["dist"] = pd.to_numeric(df["distal_length"], errors="coerce").fillna(0)
df["pend"] = pd.to_numeric(df["pendant_length"], errors="coerce").fillna(0)

is_species = (df["species_novel"] == False) & df["nearest_species"].notna()
df["species_key"] = np.where(is_species, df["nearest_species"].str.replace(" ", "_", regex=False), np.nan)

log(f"Clustering: species keys (cross-edge) + novel within-edge {SPLIT_T}...")

# -- species phylotypes: one key per raw species token, cross-edge --
new_phy_ids = {}
for key, g in df.loc[is_species].groupby("species_key"):
    for i in g.index:
        new_phy_ids[i] = key
n_species_phy = df.loc[is_species, "species_key"].nunique()

# -- novel phylotypes: within-edge single-linkage at corrected metric --
novel = df.loc[~is_species]
n_split_edges = 0
n_multi_edges = 0
for edge, g in novel.groupby("edge_num"):
    idxs = g.index.to_numpy()
    dists = g["dist"].to_numpy()
    pends = g["pend"].to_numpy()
    n = len(idxs)
    if n < 2:
        for i in idxs:
            new_phy_ids[i] = f"NOVEL_EDGE_{int(edge)}"
        continue
    n_multi_edges += 1
    parent = list(range(n))

    def _find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(n):
        for j in range(i + 1, n):
            d = abs(dists[i] - dists[j]) + pends[i] + pends[j]
            if d <= SPLIT_T:
                parent[_find(i)] = _find(j)
    comp = {}
    comp_ids = []
    for i in range(n):
        r = _find(i)
        if r not in comp:
            comp[r] = len(comp)
        comp_ids.append(comp[r])
    ncomp = len(comp)
    if ncomp > 1:
        n_split_edges += 1
    for i, c in zip(idxs, comp_ids):
        new_phy_ids[i] = f"NOVEL_EDGE_{int(edge)}" if ncomp == 1 else f"NOVEL_EDGE_{int(edge)}_{c}"

# guard: any ASV without an id (should not happen; all have placement + D)
missing = [i for i in df.index if i not in new_phy_ids]
for k, i in enumerate(missing):
    new_phy_ids[i] = f"NOVEL_UNPLACED_{k}"
df["phylotype_id"] = df.index.map(new_phy_ids)
n_novel_phy = df.loc[~is_species, "phylotype_id"].nunique()
log(f"  species phylotypes: {n_species_phy:,} | novel phylotypes: {n_novel_phy:,} | "
    f"multi-ASV edges split: {n_split_edges:,}")

# ============================================================
# 3. Phylotype-level aggregation (rep = max subject_prevalence, tie total_reads)
# ============================================================
log("Aggregating phylotype-level statistics...")
df = df.sort_values(["subject_prevalence", "total_reads"], ascending=False)

def parse_set(val):
    if pd.isna(val) or val == "":
        return set()
    if isinstance(val, set):
        return val
    s = str(val).strip("{}")
    items = []
    for sep in ["; ", ";", ", ", ","]:
        if sep in s:
            items = [x.strip().strip("'\"") for x in s.split(sep)]
            break
    if not items:
        items = [s.strip().strip("'\"")]
    return set(i for i in items if i)

phylo = df.groupby("phylotype_id").agg(
    n_asvs=("derep_id", "count"),
    rep_derep_id=("derep_id", "first"),
    rep_edge_num=("edge_num", "first"),
    rep_distal_length=("distal_length", "first"),
    rep_pendant_length=("pendant_length", "first"),
    rep_gg2_genus=("gg2_genus", "first"),
    member_derep_ids=("derep_id", lambda x: ";".join(x)),
    member_gg2_genera=("gg2_genus", lambda x: ";".join(sorted(set(str(g) for g in x if pd.notna(g))))),
    total_reads_sum=("total_reads", "sum"),
    max_n_bioprojects=("n_bioprojects", "max"),
    avg_n_bioprojects=("n_bioprojects", "mean"),
    max_prevalence=("prevalence_samples", "max"),
    max_subject_prevalence=("subject_prevalence", "max"),
    species_annotation=("species_annotation", "first"),
    species_novel=("species_novel", "first"),
    anchor_species=("species_annotation", "first"),
    anchor_genus=("nearest_genus", "first"),
    rank_annotation=("rank_annotation", "first"),
    rank_label=("rank_label", "first"),
).reset_index()

# anchor_species only meaningful for species-level phylotypes (novel -> blank)
phylo["anchor_species"] = np.where(phylo["species_novel"] == False, phylo["anchor_species"], "")

if "amplicon_regions" in df.columns:
    phylo["regions"] = phylo["phylotype_id"].map(
        df.groupby("phylotype_id")["amplicon_regions"].apply(
            lambda x: ";".join(sorted(set(r for val in x.dropna() for r in parse_set(val)))))
    )
if "bioprojects" in df.columns:
    phylo["n_bioprojects_total"] = phylo["phylotype_id"].map(
        df.groupby("phylotype_id")["bioprojects"].apply(
            lambda x: len(set(r for val in x.dropna() for r in parse_set(val)))))
if "continents" in df.columns:
    phylo["continents"] = phylo["phylotype_id"].map(
        df.groupby("phylotype_id")["continents"].apply(
            lambda x: ";".join(sorted(set(r for val in x.dropna() for r in parse_set(val)))))
    )
    phylo["n_continents"] = phylo["phylotype_id"].map(
        df.groupby("phylotype_id")["continents"].apply(
            lambda x: len(set(r for val in x.dropna() for r in parse_set(val)))))
phylo["n_contaminated_asvs"] = phylo["phylotype_id"].map(
    df.groupby("phylotype_id")["contamination_flag"].sum())
phylo["any_contamination"] = phylo["n_contaminated_asvs"] > 0

first_cols = ["phylotype_id", "rep_derep_id", "n_asvs", "member_derep_ids",
              "rep_edge_num", "rep_distal_length", "rep_pendant_length",
              "rep_gg2_genus", "member_gg2_genera",
              "total_reads_sum", "max_n_bioprojects", "avg_n_bioprojects",
              "max_prevalence", "max_subject_prevalence", "n_continents", "continents",
              "n_contaminated_asvs", "any_contamination",
              "species_annotation", "species_novel", "anchor_species", "anchor_genus",
              "rank_annotation", "rank_label"]
extra_cols = [c for c in phylo.columns if c not in first_cols]
phylo = phylo[first_cols + extra_cols]

# ============================================================
# 4. Branch check: multiple phylotypes per branch
# ============================================================
log("Branch check: multiple phylotypes per branch...")

# --- 4a. within-edge pairwise distances (exact same-edge tree distance) ---
edge_rows = []
over_pairs = []
for edge, g in phylo.groupby("rep_edge_num"):
    if len(g) < 2:
        continue
    dls = g["rep_distal_length"].to_numpy(float)
    pds = g["rep_pendant_length"].to_numpy(float)
    pid = g["phylotype_id"].to_numpy()
    dists = []
    flags = 0
    for i in range(len(pid)):
        for j in range(i + 1, len(pid)):
            d = abs(dls[i] - dls[j]) + pds[i] + pds[j]
            dists.append(d)
            if d <= SPLIT_T:
                flags += 1
                over_pairs.append((pid[i], pid[j], f"edge:{int(edge)}", "edge", round(d, 6)))
    edge_rows.append({"branch_id": int(edge), "branch_type": "edge",
                      "n_phylotypes": len(g), "n_pairs": len(dists),
                      "min_pairwise_dist": round(min(dists), 6) if dists else None,
                      "max_pairwise_dist": round(max(dists), 6) if dists else None,
                      "n_over_split": flags})

# --- 4b. cross-edge species pairs in the same clean genus (via LCA tree distance) ---
sp_phy = phylo[phylo["species_novel"] == False].copy()
sp_phy = sp_phy[sp_phy["anchor_genus"].notna() & (sp_phy["anchor_genus"] != "")]
if len(sp_phy) > 1 and len(sp_phy[["anchor_genus"]].drop_duplicates()) < len(sp_phy):
    tree = _build_tree()  # defined below
    genus_rows = defaultdict(lambda: {"branch_id": None, "branch_type": "genus",
                                      "n_phylotypes": 0, "n_pairs": 0,
                                      "min_pairwise_dist": None, "max_pairwise_dist": None,
                                      "n_over_split": 0})
    for genus, g in sp_phy.groupby("anchor_genus"):
        if len(g) < 2:
            continue
        pid = g["phylotype_id"].to_numpy()
        eds = g["rep_edge_num"].to_numpy()
        dls = g["rep_distal_length"].to_numpy(float)
        pds = g["rep_pendant_length"].to_numpy(float)
        row = genus_rows[genus]
        row["branch_id"] = genus
        row["n_phylotypes"] = len(g)
        dists = []
        for i in range(len(pid)):
            for j in range(i + 1, len(pid)):
                if eds[i] == eds[j]:
                    d = abs(dls[i] - dls[j]) + pds[i] + pds[j]
                else:
                    d = tree.point_distance(eds[i], dls[i], pds[i], eds[j], dls[j], pds[j])
                dists.append(d)
                if d <= SPLIT_T:
                    row["n_over_split"] += 1
                    over_pairs.append((pid[i], pid[j], genus, "genus", round(d, 6)))
        row["n_pairs"] = len(dists)
        row["min_pairwise_dist"] = round(min(dists), 6) if dists else None
        row["max_pairwise_dist"] = round(max(dists), 6) if dists else None
    genus_rows = list(genus_rows.values())
else:
    genus_rows = []

branch_check = pd.DataFrame(edge_rows + genus_rows)
branch_check = branch_check[["branch_id", "branch_type", "n_phylotypes", "n_pairs",
                             "min_pairwise_dist", "max_pairwise_dist", "n_over_split"]]
branch_check = branch_check.sort_values(["n_phylotypes", "n_over_split"], ascending=False)
over_df = pd.DataFrame(over_pairs, columns=["phylotype_a", "phylotype_b", "branch",
                                            "branch_type", "pairwise_dist"])
over_df = over_df.sort_values("pairwise_dist")

# ============================================================
# 5. Summary + save
# ============================================================
single_asv = (phylo["n_asvs"] == 1).sum()
multi_asv = (phylo["n_asvs"] > 1).sum()
log(f"\n=== Phylotype re-clustering summary ===")
log(f"Total phylotypes: {len(phylo):,} (species {n_species_phy:,} + novel {n_novel_phy:,})")
log(f"Single-ASV: {single_asv:,} | Multi-ASV: {multi_asv:,}")
log(f"Mean ASVs/phylotype: {phylo['n_asvs'].mean():.2f}")
log(f"Max multi-ASV size: {phylo['n_asvs'].max()}")
log(f"Over-split pairs (rep distance <= {SPLIT_T}): {len(over_df):,}")

phylo.to_csv(PIPE_V2 / "05_phylotypes_all.tsv", sep="\t", index=False)
log(f"Saved: 05_phylotypes_all.tsv ({len(phylo):,} rows)")

out_cols = ["derep_id", "phylotype_id", "edge_num", "distal_length", "pendant_length",
            "like_weight_ratio", "likelihood", "gg2_genus", "n_bioprojects", "total_reads",
            "contamination_flag", "homd_uncaptured", "is_chimera",
            "subject_prevalence", "dist_to_species", "nearest_genus", "nearest_species",
            "species_annotation", "species_novel"]
avail = [c for c in out_cols if c in df.columns]
df[avail].sort_values("derep_id").to_csv(PIPE_V2 / "05_asv_phylotype_assignments.tsv",
                                         sep="\t", index=False)
log("Saved: 05_asv_phylotype_assignments.tsv (overwritten, new phylotype_id)")

branch_check.to_csv(PIPE_V2 / "05e_branch_check.tsv", sep="\t", index=False)
log(f"Saved: 05e_branch_check.tsv ({len(branch_check):,} branches)")
over_df.to_csv(PIPE_V2 / "05e_branch_over_split.tsv", sep="\t", index=False)
log(f"Saved: 05e_branch_over_split.tsv ({len(over_df):,} pairs)")

summary = {"total_asvs_with_placement": int(len(df)),
           "total_phylotypes": int(len(phylo)),
           "species_phylotypes": int(n_species_phy),
           "novel_phylotypes": int(n_novel_phy),
           "single_asv_phylotypes": int(single_asv),
           "multi_asv_phylotypes": int(multi_asv),
           "max_multi_asv_size": int(phylo["n_asvs"].max()),
           "over_split_pairs_le_003": int(len(over_df)),
           "n_edges_with_multi_phylotype": int((branch_check["branch_type"] == "edge").sum()),
           "n_genera_with_multi_species_phylotype": int((branch_check["branch_type"] == "genus").sum())}
with open(PIPE_V2 / "05e_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
log("Saved: 05e_summary.json")
log("Done!")
