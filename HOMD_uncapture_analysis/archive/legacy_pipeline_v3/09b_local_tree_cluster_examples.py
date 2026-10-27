#!/usr/bin/env python
"""Draw local phylogenetic subtrees for clusters annotated at species level
from an ancestor's reference LCA (HARD rule, distance recorded only).

ASV tips (DEREP_*) = red ; reference tips = blue.
For each cluster: root of displayed subtree = first ancestor (>= cluster MRCA)
that contains >=1 reference tip (the node whose ref-LCA was used as annotation).
"""
import os, sys, math
import pandas as pd
from collections import Counter
from ete3 import Tree, TreeStyle, TextFace, NodeStyle

PIPE = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
TREE = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline/05_sepp/tree_export/tree.nwk"
TAX  = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v2/gappa/gg2_2022.10_tree_taxonomy.tsv"
ASV  = os.path.join(PIPE, "asv_clusters.tsv")
OUT  = os.path.join(PIPE, "analysis_report", "figs")
os.makedirs(OUT, exist_ok=True)

CLUSTERS = [36553, 29573, 14392, 9690]

ASV_COLOR = "#c0392b"   # red
REF_COLOR = "#2471a3"   # blue

# ---------------------------------------------------------------- taxonomy
tax = {}
for line in open(TAX):
    p = line.rstrip("\n").split("\t")
    if len(p) == 2:
        tax[p[0]] = p[1]

def species_of(lab):
    t = tax.get(lab, "")
    if not t:
        return ""
    for part in t.split(";"):
        part = part.strip()
        if part.startswith("s__") and part[3:]:
            return part[3:].replace("_", " ")
    for part in t.split(";"):
        part = part.strip()
        if part.startswith("g__") and part[3:]:
            return part[3:].replace("_", " ") + " sp."
    return ""

def short_label(lab):
    """Short readable id for a reference tip."""
    return lab

# ------------------------------------------------------------------ asv table
a = pd.read_csv(ASV, sep="\t", low_memory=False)
a = a[a["cluster_id"].isin(CLUSTERS)]
info = {}
for _, r in a.iterrows():
    cid = int(r["cluster_id"])
    info.setdefault(cid, {"members": [], "annotation": r["annotation"],
                          "level": r["annotation_level"], "dist": r["annotation_distance"],
                          "anchor": r["annotation_anchor"]})
    info[cid]["members"].append(r["asv_id"])
for cid in CLUSTERS:
    info[cid]["members"] = sorted(set(info[cid]["members"]))

# -------------------------------------------------------------------- tree
print("loading tree ...", flush=True)
t = Tree(TREE, format=1, quoted_node_names=True)
leaves = {n.name: n for n in t.iter_leaves()}
print("tips:", len(leaves), flush=True)

def is_ref(name):
    return not name.startswith("DEREP_")

def edge_dist(node, anc):
    d, cur = 0.0, node
    while cur is not None and cur is not anc:
        d += (cur.dist or 0.0)
        cur = cur.up
    return d

rows = []
for cid in CLUSTERS:
    mem = info[cid]["members"]
    nodes = [leaves[m] for m in mem]
    mrca = t.get_common_ancestor(nodes)
    # first ancestor with >=1 ref tip
    anc = mrca
    hops = 0
    while anc is not None:
        refs = [x.name for x in anc.iter_leaves() if is_ref(x.name)]
        if refs:
            break
        anc = anc.up
        hops += 1
    if anc is None:
        print(f"cluster {cid}: NO ref ancestor found!", flush=True)
        continue
    refs = [x.name for x in anc.iter_leaves() if is_ref(x.name)]
    asvs = [x.name for x in anc.iter_leaves() if not is_ref(x.name)]
    d_tip = [round(edge_dist(leaves[m], anc), 6) for m in mem]
    d_mrca = round(edge_dist(mrca, anc), 6)
    # distance from mrca down to asv tips (cluster internal depth)
    d_in = [round(edge_dist(leaves[m], mrca), 6) for m in mem]
    sps = sorted({species_of(x) for x in refs})
    print(f"cluster {cid}: members={mem} n_tips={len(asvs)} n_refs={len(refs)} "
          f"hops={hops} d_tip={d_tip} d_mrca_to_anc={d_mrca} ref_species={sps}", flush=True)
    rows.append((cid, mem, mrca, anc))

# ------------------------------------------------------------------ drawing
def build_subtree(cid, mem, anc):
    nwk = anc.write(format=1, quoted_node_names=True)
    sub = Tree(nwk, format=1, quoted_node_names=True)
    s_leaves = {n.name: n for n in sub.iter_leaves()}
    members_sub = [s_leaves[m] for m in mem]
    cr = sub.get_common_ancestor(members_sub)
    cr_leaves = {x.name for x in cr.iter_leaves()}
    if cr_leaves != set(mem):
        print(f"  !! cluster {cid}: mrca descendants != members -> "
              f"{sorted(cr_leaves)}", flush=True)
    else:
        print(f"  ok cluster {cid}: mrca clade == exactly {len(mem)} members", flush=True)
    return sub, cr

def draw(cid, mem, anc, dist_tsv):
    sub, cr = build_subtree(cid, mem, anc)
    ts = TreeStyle()
    ts.mode = "r"
    ts.show_leaf_name = False
    ts.show_scale = True
    ts.show_branch_length = False
    ts.branch_vertical_margin = 6
    ts.min_leaf_separation = 10

    ann = info[cid]["annotation"]
    title = TextFace(f"cluster {cid}  |  {ann}  (ancestor-LCA, species; d={dist_tsv})",
                     fsize=15, bold=True)
    ts.title.add_face(title, column=0)

    for leaf in sub.iter_leaves():
        if leaf.name.startswith("DEREP_"):
            if leaf.name in mem:
                col = ASV_COLOR          # cluster member
                size, fsize, bold = 10, 10, True
                lab = leaf.name
            else:
                col = "#e67e22"          # other ASV in same ancestor subtree
                size, fsize, bold = 6, 8, False
                lab = leaf.name
        else:
            col = REF_COLOR
            size, fsize, bold = 7, 9, False
            sp = species_of(leaf.name)
            lab = leaf.name if not sp else f"{sp} [{leaf.name}]"
        ns = NodeStyle()
        ns["fgcolor"] = col
        ns["size"] = size
        ns["hz_line_color"] = col
        ns["vt_line_color"] = col
        leaf.set_style(ns)
        f = TextFace(lab, fsize=fsize, fgcolor=col, bold=bold)
        leaf.add_face(f, column=0, position="branch-right")

    # mark cluster clade
    ns = NodeStyle()
    ns["fgcolor"] = ASV_COLOR
    ns["size"] = 14
    ns["shape"] = "square"
    ns["bgcolor"] = "#fadbd8"
    cr.set_style(ns)
    cr.add_face(TextFace(f"cluster {cid}: {len(mem)} ASVs, "
                         "no reference inside (d to ref ancestor = "
                         f"{dist_tsv})", fsize=10, fgcolor=ASV_COLOR),
                column=1, position="branch-right")

    # legend (left of root)
    for i, (txt, fg) in enumerate([
        ("inserted ASV (cluster member)", ASV_COLOR),
        ("reference sequence", REF_COLOR),
        ("other ASV under same ancestor", "#e67e22"),
    ]):
        lf = TextFace(txt, fsize=10, fgcolor=fg, bold=True)
        sub.add_face(lf, column=-1 - i, position="branch-right")

    png = os.path.join(OUT, f"fig_local_tree_cluster_{cid}.png")
    pdf = os.path.join(OUT, f"fig_local_tree_cluster_{cid}.pdf")
    sub.render(png, tree_style=ts, w=1200, units="px", dpi=150)
    sub.render(pdf, tree_style=ts, w=220, units="mm", dpi=300)
    print("wrote", png, pdf, flush=True)

for cid, mem, mrca, anc in rows:
    draw(cid, mem, anc, info[cid]["dist"])
