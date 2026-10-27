#!/usr/bin/env python3
"""
Step 10 - TreeCluster OTU-like clustering and tree-based ASV annotation.

Goal: group the SEPP-inserted ASVs into OTU-like clusters on the GG2 2022.10
tree (keeping the backbone reference tips as cluster members).

Annotation logic (user, 2026-09-01; HARD rule, no distance gate):
  1. Cluster containing reference tips: annotation = LCA of the taxonomy of the
     reference tips INSIDE the cluster (cluster clade root's ref summary).
     Per-ASV annotation_distance = minimum tree distance from the ASV to a
     reference tip in the same cluster (by max_clade t=0.03 construction this
     distance is always <= 0.03).
  2. Cluster with NO reference inside: annotation = LCA of the references under
     the FIRST ancestor (one level up and upwards) of the cluster clade root
     that contains reference tips (the first reference-bearing ancestor above
     the cluster clade). Per-ASV
     annotation_distance = tree distance from the ASV tip to that ancestor
     node. The ancestor LCA is accepted AT ITS OWN RANK regardless of how large
     the distance is (d>0.03 or d<0.03 is not considered).
  3. Singleton ASV (no cluster partner): annotation = LCA of the references
     under the NEAREST ancestor of the tip that contains reference tips
     (the nearest reference-bearing ancestor); it may NOT borrow the reference
     strains of any existing
     cluster. Same rule as (2): no distance gate, the ancestor LCA is accepted
     at its own rank; if the walk-up reaches the root with no reference-bearing
     ancestor, the ASV is left unannotated (level 'none').
  4. Annotation is completed FIRST; the HOMD relationship is computed AFTER
     (annotation_in_homd: species-level annotation exactly in the HOMD catalog
     or a unique GTDB->HOMD renamed epithet; annotation_genus_in_homd: first
     token of the annotation in the HOMD genus set).

Nearest-reference columns (nearest_ref*, nearest_homd_ref*, covered) are KEPT
as informational legacy columns only; they no longer drive classification.

Outputs (10_clustering/results/):
  - treecluster/10_clusters.tsv : TreeCluster raw output (all tips)
  - 10_asv_clusters.tsv         : per-ASV cluster + tree-based annotation
  - 10_clusters.tsv             : per-cluster summary (ref species lists)
  - 10_summary.json
"""
import argparse
import gzip
import json
import re
import subprocess
import sys
import time
from collections import Counter, deque
from pathlib import Path

import pandas as pd
from treeswift import read_tree_newick

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE = BASE / "pipeline"
TREECLUSTER = "/home/huangxc/.local/bin/TreeCluster.py"
TREE = PIPE / "inputs" / "tree.nwk"
TAXONOMY = PIPE / "inputs" / "gg2_2022.10_tree_taxonomy.tsv.gz"
HOMD_FASTA = Path("/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03_full.fasta")

log = lambda s: print(f"[10_clustering] {s}", flush=True)

RANK_ORDER = ["d", "p", "c", "o", "f", "g", "s"]
LEVEL_NAMES = {"d": "domain", "p": "phylum", "c": "class", "o": "order",
               "f": "family", "g": "genus", "s": "species"}


def load_homd_catalog(path):
    homd_species, homd_epithets, homd_genus = set(), {}, set()
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                parts = line.split()
                if len(parts) >= 3:
                    homd_genus.add(parts[2].strip())
                if len(parts) >= 4:
                    sp = f"{parts[2].strip()} {parts[3].strip()}"
                    homd_species.add(sp)
                    homd_epithets.setdefault(parts[3].strip(), set()).add(sp)
    return homd_species, homd_epithets, homd_genus


def get_tax_level(tax_str, level):
    if not tax_str:
        return None
    for part in tax_str.split("; "):
        if part.startswith(level + "__"):
            v = part[len(level) + 2:].strip()
            return v if v else None
    return None


def clean_species(s):
    """GG2 s__ ('Treponema_B_986142 vincentii') -> 'Genus species', or None.

    Placeholder species (spNNN) are not names -> None (unresolved)."""
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
    """GTDB g__ ('Leptotrichia_A_993758') -> 'Leptotrichia' (strip strain token)."""
    if not isinstance(s, str) or not s:
        return None
    return s.split("_")[0]


def renamed_in_homd(sp, homd_epithets):
    """GTDB genus != HOMD genus but the epithet matches exactly one HOMD species."""
    if not isinstance(sp, str):
        return False
    parts = sp.split(" ")
    if len(parts) < 2:
        return False
    return len(homd_epithets.get(parts[-1].strip(), ())) == 1


def is_placeholder_species(s):
    """GTDB placeholder species epithet (e.g. 'Dialister sp001553355')."""
    if not isinstance(s, str) or not s:
        return False
    return bool(re.fullmatch(r"sp\d+", s.split(" ")[-1]))


# ---------------------------------------------------------------- summaries
# A "tax summary" of a node = (stopped, labels[7]) where labels[i] is the
# single taxonomy label shared by ALL references in the node's subtree at rank
# i (RANK_ORDER), '' where not claimed; `stopped` = first rank index that is
# not claimed (7 = species fully resolved). Missing labels act as wildcards
# (no constraint), EXCEPT a missing or placeholder species which blocks the
# species rank for the whole subtree (conservative: never over-claim species).

def leaf_tax_summary(tax_str):
    labels = [get_tax_level(tax_str, lv) for lv in RANK_ORDER]
    if labels[5]:
        labels[5] = clean_genus(labels[5])
    sp = labels[6]
    if not sp or is_placeholder_species(sp):
        labels[6] = ""
        stopped = 6
    else:
        cleaned = clean_species(sp)
        labels[6] = cleaned if cleaned else ""
        stopped = 7 if labels[6] else 6
    return (stopped, labels)


def merge_summaries(children):
    """LCA summary of a node from its children summaries (children with refs)."""
    if not children:
        return None
    labels = [""] * 7
    stopped = 7
    n_children = len(children)
    for i in range(7):
        cands = [c for c in children if c[0] > i]
        if len(cands) != n_children:
            stopped = i
            break
        vals = {c[1][i] for c in cands}
        vals.discard("")
        if len(vals) != 1:
            stopped = i
            break
        labels[i] = vals.pop()
    return (stopped, labels)


def summary_annotation(summary):
    """(label, level_name) from a tax summary; ('', 'none') if nothing claimed."""
    if summary is None or summary[0] == 0:
        return "", "none"
    i = summary[0] - 1
    return summary[1][i], LEVEL_NAMES[RANK_ORDER[i]]


# ---------------------------------------------------------------- tree utils

def nearest_ref_dp(tree, ref_set):
    """Bottom-up: best reference tip (label + dist) inside each subtree."""
    for node in tree.traverse_postorder():
        if node.is_leaf():
            if node.label in ref_set:
                node.best_ref, node.best_dist = node.label, 0.0
            else:
                node.best_ref, node.best_dist = None, float("inf")
        else:
            best_ref, best_dist = None, float("inf")
            for c in node.children:
                if c.best_ref is not None:
                    d = c.best_dist + (c.edge_length or 0.0)
                    if d < best_dist:
                        best_dist, best_ref = d, c.best_ref
            node.best_ref, node.best_dist = best_ref, best_dist


def nearest_ref_for_tip(tip):
    """Walk up from a tip, combining best refs of sibling subtrees (exact)."""
    node, child = tip.parent, tip
    dist_to_node = tip.edge_length or 0.0
    best_ref, best_dist = None, float("inf")
    while node is not None:
        for c in node.children:
            if c is child or c.best_ref is None:
                continue
            d = dist_to_node + c.best_dist + (c.edge_length or 0.0)
            if d < best_dist:
                best_dist, best_ref = d, c.best_ref
        child = node
        dist_to_node += (node.edge_length or 0.0)
        node = node.parent
    return best_ref, best_dist


def cluster_root_node(leaf_map, tips):
    """MRCA of a list of tip labels by lock-step ancestor counting."""
    if len(tips) == 1:
        return leaf_map[tips[0]]
    count = Counter()
    q = deque(leaf_map[t] for t in tips)
    n = len(tips)
    while q:
        cur = q.popleft()
        count[cur] += 1
        if count[cur] == n:
            return cur
        if cur.parent is not None:
            q.append(cur.parent)
    raise RuntimeError("no MRCA (should never happen)")


def attach_local_refs(root, cluster_ref_set):
    """Local DP: nearest ref of THIS cluster inside the cluster subtree."""
    for node in root.traverse_postorder():
        if node.is_leaf():
            if node.label in cluster_ref_set:
                node.lref, node.ldist = node.label, 0.0
            else:
                node.lref, node.ldist = None, float("inf")
        else:
            best, bd = None, float("inf")
            for c in node.children:
                if c.lref is not None:
                    d = c.ldist + (c.edge_length or 0.0)
                    if d < bd:
                        bd, best = d, c.lref
            node.lref, node.ldist = best, bd


def min_dist_to_cluster_ref(tip, root):
    """Exact min tree distance from an ASV tip to a ref of the same cluster."""
    node, child = tip.parent, tip
    d0 = tip.edge_length or 0.0
    best, bd = None, float("inf")
    while node is not None and node is not root:
        for c in node.children:
            if c is child or c.lref is None:
                continue
            d = d0 + c.ldist + (c.edge_length or 0.0)
            if d < bd:
                bd, best = d, c.lref
        child = node
        d0 += (node.edge_length or 0.0)
        node = node.parent
    if node is not None:
        for c in root.children:
            if c is child or c.lref is None:
                continue
            d = d0 + c.ldist + (c.edge_length or 0.0)
            if d < bd:
                bd, best = d, c.lref
    return (bd if best is not None else float("inf"))


def dist_to_ancestor(tip, anc):
    d = 0.0
    node = tip
    while node is not anc and node is not None:
        d += (node.edge_length or 0.0)
        node = node.parent
    return d


def first_ref_ancestor(node):
    """First node at-or-above `node` whose subtree contains >=1 reference tip.
    Returns (anchor_node, hops_up)."""
    hops = 0
    cur = node
    while cur is not None and not cur.has_ref:
        cur = cur.parent
        hops += 1
    return cur, hops


def ancestor_annotation(anc_tsum, basis, anchor_str, dist):
    """Annotate from an ancestor's reference LCA at its own rank (HARD rule).

    dist = tree distance from the ASV tip up to the ancestor node whose
    reference LCA is used; it is RECORDED ONLY, never used to gate the rank:
    the ancestor LCA is accepted at its own rank (species / genus / ...)
    regardless of whether dist > 0.03 or < 0.03. If the ancestor has no usable
    labels the ASV is left unannotated (level 'none').
    """
    label, level = summary_annotation(anc_tsum)
    return label, level, basis, round(float(dist), 6), anchor_str


def lineage_from_summary(summary, max_rank_idx=6):
    """Full d__..s__ lineage string from a tax summary (blank below max_rank_idx)."""
    if summary is None:
        return ""
    labels = [str(x) if x is not None else "" for x in summary[1][:max_rank_idx + 1]]
    labels += [""] * (7 - len(labels))
    return "; ".join(labels)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tree", default=str(TREE))
    ap.add_argument("--threshold", type=float, default=0.03)
    ap.add_argument("--method", default="max_clade")
    ap.add_argument("--treecluster", default=TREECLUSTER)
    ap.add_argument("--treecluster-output", default=None,
                    help="Reuse an existing TreeCluster output instead of re-running")
    ap.add_argument("--taxonomy", default=str(TAXONOMY))
    ap.add_argument("--homd-fasta", default=str(HOMD_FASTA))
    ap.add_argument("--outdir", default=str(PIPE / "10_clustering" / "results"))
    args = ap.parse_args()

    outdir = Path(args.outdir)
    tc_dir = outdir / "treecluster"
    tc_dir.mkdir(parents=True, exist_ok=True)
    tc_out = Path(args.treecluster_output) if args.treecluster_output else tc_dir / "10_clusters.tsv"

    # 1. TreeCluster
    if not tc_out.exists():
        log(f"Running TreeCluster {args.method} t={args.threshold} ...")
        t0 = time.time()
        subprocess.run([args.treecluster, "-i", args.tree, "-t", str(args.threshold),
                        "-m", args.method, "-o", str(tc_out)], check=True)
        log(f"  TreeCluster done in {time.time()-t0:.0f}s")
    else:
        log(f"Reusing TreeCluster output: {tc_out}")

    cl = pd.read_csv(tc_out, sep="\t", dtype={"SequenceName": str, "ClusterNumber": int})
    cl["is_asv"] = cl["SequenceName"].str.startswith("DEREP_")
    log(f"TreeCluster: {len(cl):,} tips | {cl['ClusterNumber'].nunique():,} clusters")

    # 2. Taxonomy + HOMD catalog
    tax = {}
    opener = gzip.open if str(args.taxonomy).endswith(".gz") else open
    with opener(args.taxonomy, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2:
                tax[parts[0]] = parts[1]
    homd_species, homd_epithets, homd_genus = load_homd_catalog(args.homd_fasta)
    log(f"taxonomy tips: {len(tax):,} | HOMD species: {len(homd_species):,} | "
        f"HOMD genera: {len(homd_genus):,}")

    # 3. Load tree
    log(f"Loading tree: {args.tree}")
    t0 = time.time()
    tree = read_tree_newick(args.tree)
    n_tips = sum(1 for _ in tree.traverse_leaves())
    log(f"  loaded {n_tips:,} tips in {time.time()-t0:.0f}s")

    ref_set = set()
    homd_ref_set = set()
    for node in tree.traverse_leaves():
        if node.label is None:
            continue
        if node.label.startswith("DEREP_"):
            continue
        ref_set.add(node.label)
        sp = clean_species(get_tax_level(tax.get(node.label), "s"))
        if sp is not None and (sp in homd_species or renamed_in_homd(sp, homd_epithets)):
            homd_ref_set.add(node.label)
    log(f"ref tips: {len(ref_set):,} | ref tips w/ HOMD species: {len(homd_ref_set):,}")

    # 4. DP + walk-up: nearest ref (any), then nearest HOMD-species ref (legacy info)
    asv_tips = [n for n in tree.traverse_leaves() if n.label and n.label.startswith("DEREP_")]
    log(f"ASV tips: {len(asv_tips):,}")

    log("DP: nearest reference ...")
    t0 = time.time()
    nearest_ref_dp(tree, ref_set)
    log(f"  done in {time.time()-t0:.0f}s")
    log("Walking up for nearest reference ...")
    t0 = time.time()
    rows = []
    for tip in asv_tips:
        r, d = nearest_ref_for_tip(tip)
        rows.append((tip.label, r, d))
    log(f"  done in {time.time()-t0:.0f}s")
    asv_nn = pd.DataFrame(rows, columns=["asv_id", "nearest_ref", "nearest_ref_dist"])
    asv_nn["nearest_ref_dist"] = asv_nn["nearest_ref_dist"].round(6)

    log("DP: nearest HOMD-species reference ...")
    t0 = time.time()
    nearest_ref_dp(tree, homd_ref_set)
    log(f"  done in {time.time()-t0:.0f}s")
    log("Walking up for nearest HOMD-species ref ...")
    t0 = time.time()
    rows_h = []
    for tip in asv_tips:
        rh, dh = nearest_ref_for_tip(tip)
        rows_h.append((tip.label, rh, dh))
    log(f"  done in {time.time()-t0:.0f}s")
    asv_h = pd.DataFrame(rows_h, columns=["asv_id", "nearest_homd_ref", "nearest_homd_ref_dist"])
    asv_h["nearest_homd_ref_dist"] = asv_h["nearest_homd_ref_dist"].round(6)

    # 5. Node summaries: LCA of all refs below each node + has_ref
    log("Computing per-node reference summaries (bottom-up LCA) ...")
    t0 = time.time()
    nid = 0
    for node in tree.traverse_postorder():
        node.nid = nid
        nid += 1
        if node.is_leaf():
            if node.label in ref_set:
                node.has_ref, node.nrefs = True, 1
                node.tsum = leaf_tax_summary(tax.get(node.label, ""))
            else:
                node.has_ref, node.nrefs = False, 0
                node.tsum = None
        else:
            has = False
            nrefs = 0
            tsums = []
            for c in node.children:
                nrefs += c.nrefs
                if c.has_ref:
                    has = True
                    tsums.append(c.tsum)
            node.has_ref, node.nrefs = has, nrefs
            node.tsum = merge_summaries(tsums) if has else None
    log(f"  done in {time.time()-t0:.0f}s  ({nid:,} nodes)")

    # 6. Merge with clusters + taxonomy (+ legacy nearest-ref / HOMD columns)
    asv = cl[cl["is_asv"]].rename(columns={"SequenceName": "asv_id", "ClusterNumber": "cluster_id"})
    asv = asv.merge(asv_nn, on="asv_id", how="left").merge(asv_h, on="asv_id", how="left")
    asv["nearest_ref_tax"] = asv["nearest_ref"].map(tax)
    asv["nearest_ref_species"] = asv["nearest_ref_tax"].map(lambda t: clean_species(get_tax_level(t, "s")))
    asv["nearest_ref_genus"] = asv["nearest_ref_tax"].map(lambda t: get_tax_level(t, "g"))
    asv["nearest_homd_ref_tax"] = asv["nearest_homd_ref"].map(tax)
    asv["nearest_homd_ref_species"] = asv["nearest_homd_ref_tax"].map(lambda t: clean_species(get_tax_level(t, "s")))
    asv["nearest_ref_in_homd"] = asv["nearest_ref_species"].map(lambda s: s in homd_species if s else False)
    asv["nearest_ref_renamed"] = asv["nearest_ref_species"].map(
        lambda s: renamed_in_homd(s, homd_epithets) if s else False)
    asv["covered"] = asv["nearest_homd_ref_dist"] <= args.threshold
    asv["annotation_resolution"] = asv["nearest_ref_tax"].map(
        lambda t: "species" if get_tax_level(t, "s") else
                  ("genus" if get_tax_level(t, "g") else
                   ("family" if get_tax_level(t, "f") else
                    ("order" if get_tax_level(t, "o") else
                     ("class" if get_tax_level(t, "c") else
                      ("phylum" if get_tax_level(t, "p") else "none"))))))

    # cluster-level info
    cl_size = cl.groupby("ClusterNumber").size().rename("cluster_size")
    n_ref = cl[~cl["is_asv"]].groupby("ClusterNumber").size().rename("n_refs_in_cluster")
    asv = asv.merge(cl_size, left_on="cluster_id", right_index=True, how="left")
    asv = asv.merge(n_ref, left_on="cluster_id", right_index=True, how="left")
    asv["n_refs_in_cluster"] = asv["n_refs_in_cluster"].fillna(0).astype(int)
    asv.loc[asv["cluster_id"] == -1, "cluster_size"] = 1
    asv.loc[asv["cluster_id"] == -1, "n_refs_in_cluster"] = 0

    # any cluster ref species in HOMD (info only; same concept as v2)
    refs = cl[(~cl["is_asv"]) & (cl["ClusterNumber"] >= 0)].copy()
    refs["species"] = refs["SequenceName"].map(lambda n: clean_species(get_tax_level(tax.get(n), "s")))
    refs["sp_in_homd"] = refs["species"].map(lambda s: s in homd_species if s else False)
    refs["sp_renamed"] = refs["species"].map(lambda s: renamed_in_homd(s, homd_epithets) if s else False)
    cl_homd = refs.groupby("ClusterNumber").agg(
        any_in_homd=("sp_in_homd", "any"),
        any_renamed=("sp_renamed", "any"),
        ref_species_list=("species", lambda s: ";".join(sorted(set(x for x in s.dropna())))),
        homd_species_list=("species", lambda s: ";".join(sorted(set(
            x for x in s.dropna() if x in homd_species or renamed_in_homd(x, homd_epithets))))),
    ).reset_index()
    asv = asv.merge(cl_homd, left_on="cluster_id", right_on="ClusterNumber", how="left")
    asv["any_cluster_ref_species_in_homd"] = (asv["cluster_id"] >= 0) & (
        asv["any_in_homd"].fillna(False) | asv["any_renamed"].fillna(False))

    # ---- 7. Tree-based annotation (user logic) -----------------------------
    log("Tree-based annotation (cluster ref LCA / common ancestor / nearest ancestor) ...")
    t0 = time.time()
    leaf_map = {n.label: n for n in tree.traverse_leaves() if n.label}
    cluster_map = {}
    for r in cl.itertuples(index=False):
        cluster_map.setdefault(r.ClusterNumber, []).append(r.SequenceName)

    n_ref_ok = n_noref = n_sing = 0
    ann = {}
    for cid, tips in cluster_map.items():
        asv_members = [t for t in tips if t.startswith("DEREP_")]
        if not asv_members:
            continue
        ref_members = [t for t in tips if not t.startswith("DEREP_")]
        if cid == -1:
            # TreeCluster bucket -1: every unclustered ASV is its own singleton
            for t in asv_members:
                tip_node = leaf_map[t]
                anc, hops = first_ref_ancestor(tip_node)
                if anc is None:
                    ann[t] = ("", "none", "singleton_ancestor", float("nan"), "no-ref-ancestor", "")
                    continue
                d = dist_to_ancestor(tip_node, anc)
                a5 = ancestor_annotation(
                    anc.tsum, "singleton_ancestor",
                    f"ancestor:nid{anc.nid}:{anc.nrefs}refs:{hops}up",
                    d)
                lineage = lineage_from_summary(anc.tsum)
                ann[t] = a5 + (lineage,)
                n_sing += 1
            continue
        root = cluster_root_node(leaf_map, tips)
        if ref_members:
            ref_set_this = set(ref_members)
            if root.nrefs != len(ref_set_this):
                # safety fallback: use the cluster refs' own taxonomy directly
                tstrs = [tax.get(x) for x in ref_members if tax.get(x)]
                summary = leaf_tax_summary(tstrs[0]) if len(tstrs) == 1 else merge_summaries(
                    [leaf_tax_summary(t) for t in tstrs]) if tstrs else None
            else:
                summary = root.tsum
            attach_local_refs(root, ref_set_this)
            label, level = summary_annotation(summary)
            lineage = lineage_from_summary(summary)
            for t in asv_members:
                tip_node = leaf_map[t]
                d = min_dist_to_cluster_ref(tip_node, root)
                ann[t] = (label, level, "cluster_refs", round(d, 6),
                          f"cluster:nid{root.nid}:{len(ref_set_this)}refs", lineage)
                n_ref_ok += 1
        else:
            anc, hops = first_ref_ancestor(root)
            if anc is None:
                for t in asv_members:
                    ann[t] = ("", "none", "cluster_ancestor", float("nan"), "no-ref-ancestor", "")
                continue
            for t in asv_members:
                tip_node = leaf_map[t]
                d = dist_to_ancestor(tip_node, anc)
                a5 = ancestor_annotation(
                    anc.tsum, "cluster_ancestor",
                    f"ancestor:nid{anc.nid}:{anc.nrefs}refs:{hops}up",
                    d)
                lineage = lineage_from_summary(anc.tsum)
                ann[t] = a5 + (lineage,)
                n_noref += 1

    ann_df = pd.DataFrame.from_dict(ann, orient="index",
                                    columns=["annotation", "annotation_level",
                                             "annotation_basis", "annotation_distance",
                                             "annotation_anchor", "annotation_lineage"])
    ann_df.index.name = "asv_id"
    asv = asv.merge(ann_df.reset_index(), on="asv_id", how="left")
    asv["annotation"] = asv["annotation"].fillna("")
    asv["annotation_level"] = asv["annotation_level"].fillna("none")
    asv["annotation_basis"] = asv["annotation_basis"].fillna("none")
    # HOMD relationship computed AFTER annotation
    asv["annotation_in_homd"] = asv.apply(
        lambda r: bool(r["annotation_level"] == "species" and r["annotation"]
                       and (r["annotation"] in homd_species
                            or renamed_in_homd(r["annotation"], homd_epithets))), axis=1)
    asv["annotation_genus_in_homd"] = asv["annotation"].map(
        lambda a: bool(str(a)) and str(a).split(" ")[0] in homd_genus)
    log(f"  annotation done in {time.time()-t0:.0f}s  "
        f"(by ref LCA: {n_ref_ok:,} | common-ancestor: {n_noref:,} | singleton: {n_sing:,})")

    # 8. Outputs
    out_cols = ["asv_id", "cluster_id", "cluster_size", "n_refs_in_cluster",
                "annotation", "annotation_level", "annotation_basis",
                "annotation_distance", "annotation_anchor", "annotation_lineage",
                "annotation_in_homd", "annotation_genus_in_homd",
                "nearest_ref", "nearest_ref_species", "nearest_ref_genus", "nearest_ref_dist",
                "nearest_ref_in_homd", "nearest_ref_renamed",
                "nearest_homd_ref", "nearest_homd_ref_species", "nearest_homd_ref_dist",
                "covered", "any_cluster_ref_species_in_homd", "annotation_resolution"]
    asv[out_cols].to_csv(outdir / "10_asv_clusters.tsv", sep="\t", index=False)

    cl_out = cl.groupby("ClusterNumber").size().rename("n_total").reset_index()
    cl_out = cl_out[cl_out["ClusterNumber"] >= 0].rename(columns={"ClusterNumber": "cluster_id"})
    cl_out["n_asv"] = cl[cl["is_asv"]].groupby("ClusterNumber").size().reindex(
        cl_out["cluster_id"]).fillna(0).astype(int).values
    cl_out["n_refs_in_cluster"] = n_ref.reindex(cl_out["cluster_id"]).fillna(0).astype(int).values
    cl_out = cl_out.merge(cl_homd, left_on="cluster_id", right_on="ClusterNumber", how="left")
    cl_out["ref_species_list"] = cl_out["ref_species_list"].fillna("")
    cl_out["homd_species_list"] = cl_out["homd_species_list"].fillna("")
    cl_out["any_in_homd"] = cl_out["any_in_homd"].fillna(False)
    cl_out["any_renamed"] = cl_out["any_renamed"].fillna(False)
    cl_out = cl_out[["cluster_id", "n_total", "n_asv", "n_refs_in_cluster",
                     "ref_species_list", "homd_species_list", "any_in_homd", "any_renamed"]]
    cl_out.to_csv(outdir / "10_clusters.tsv", sep="\t", index=False)

    summary = {
        "threshold": args.threshold,
        "method": args.method,
        "annotation_logic": "HARD rule (no distance gate): cluster refs -> LCA of "
                            "refs inside cluster; ref-less cluster -> LCA under first "
                            "ancestor with refs (common ancestor one level up); "
                            "singleton -> nearest ancestor with refs (never borrows "
                            "reference strains of other clusters); ancestor-based "
                            "annotations accepted at their own rank regardless of "
                            "ASV-to-ancestor distance (d>0.03 or d<0.03); "
                            "HOMD relation computed AFTER annotation",
        "n_asv": int(len(asv)),
        "n_clusters": int(cl["ClusterNumber"].nunique()),
        "n_asv_singletons": int((asv["cluster_id"] == -1).sum()),
        "n_asv_annotation_in_homd": int(asv["annotation_in_homd"].sum()),
        "n_asv_annotation_genus_in_homd": int(asv["annotation_genus_in_homd"].sum()),
        "n_asv_covered_legacy": int(asv["covered"].sum()),
        "annotation_level_counts": asv["annotation_level"].value_counts().to_dict(),
        "annotation_basis_counts": asv["annotation_basis"].value_counts().to_dict(),
    }
    with open(outdir / "10_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    log(f"summary: {json.dumps(summary)}")
    log(f"wrote {outdir/'10_asv_clusters.tsv'}, {outdir/'10_clusters.tsv'}, {outdir/'10_summary.json'}")


if __name__ == "__main__":
    main()
