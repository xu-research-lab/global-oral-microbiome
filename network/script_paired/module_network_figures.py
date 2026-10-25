#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
module_network_figures.py — Module-structured network figure for the four
HC/PD/DC/OSCC pooled genus backbones, styled after
result_disease/module_networks/module_network_healthy_disease_combined_roles.png
(scripts_disease/module_network_combined_healthy_disease.py):

  - modules laid out in rows of AT MOST FIVE panels (the four groups ->
    1 row x 4 columns), grey facet-strip title band above each panel
  - within each module the nodes sit on a CIRCLE, ordered clockwise by mean
    relative abundance (most abundant at 12 o'clock)
  - role marking (Guimera-Amaral cartography): module hub = max
    within-module z (yellow outline), main connector = max participation P
    (green outline), both on one node = purple; in role mode the degree hubs
    lose their black outline
  - node sizes share ONE global mean-abundance basis so panels are comparable
  - one shared dashed-box legend BELOW the panels: correlation colours,
    cartography roles, mean relative abundance circles

Two elements are intentionally absent (no longer needed, per decision):
  - the module-strength background shading (white -> orange hull fill)
  - the intra-module edge-density badge text ("density=0.xx") and the
    "Module Strength" gradient legend section

Inputs (outputs of this pipeline, nothing re-pooled):
  result_paired/backbone_v2/<grp>_genus_uw/backbone_edges.csv
  result_paired/modularity/modules_<grp>.csv
  result_paired/biom_paired/<grp>_genus/*.tsv            (genus counts)

Outputs (under result_paired/module_networks/):
  module_networks_1x4.png / .pdf       the 1x4 combined deliverable
  module_network_<grp>_genus.png       standalone panel per group (legend incl.)

Usage (from network/):
  python script_hc_pd_dc_oscc/module_network_figures.py
"""
from __future__ import annotations

import argparse
import glob
import itertools
import os

import numpy as np
import pandas as pd
import networkx as nx

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.text as mtext
import matplotlib.transforms as mtransforms
from matplotlib.patches import Rectangle

GROUPS = ["HC", "PD", "DC", "OSCC"]
GRP_LABEL = {"HC": "Healthy control", "PD": "Periodontal disease",
             "DC": "Dental caries", "OSCC": "OSCC"}

ROOT_DIR = "result_paired"
BACKBONE_DIR = "backbone_v2"
MODULARITY_DIR = "modularity"
BIOM_DIR = "biom_paired"

MOD_COLORS = ["#1b6ca8", "#d95f02", "#1a9e77", "#7570b3", "#e6ab02",
              "#e15759", "#4e79a7", "#59a14f", "#b07aa1", "#ff9da7",
              "#9c755f", "#bab0ac", "#86bcb6", "#d37295", "#499894"]

POS_C, NEG_C = "#d64541", "#3b74c4"
# cartography role colours (user choice): module hub = yellow, main
# connector = green, hub+connector stays purple
ROLE_RED, ROLE_BLUE, ROLE_PURPLE = "#ff8c00", "#2ca02c", "#7b2fa0"

# highlighted key genera (user selection, unified emphasis ring + bold
# enlarged label): node FILL and the yellow/green/purple role outlines stay
# as-is; a magenta ring is drawn UNDER the node (outside the role outline).
HL_C = "#b3006b"
HIGHLIGHT = {
    "PD": ["Fretibacterium", "Treponema_B_986141", "Treponema_B_986142",
           "Treponema_C", "Treponema_D", "Filifactor", "Streptococcus",
           "Granulicatella", "Porphyromonas_A", "Tannerella"],
    "OSCC": ["Neisseria_563205", "Haemophilus_D_735815", "Parvimonas"],
    "DC": ["Veillonella_A", "Selenomonas_A_42762"],
}
# HC shows every highlighted genus of the other panels *if present* there
HIGHLIGHT["HC"] = sorted(set().union(*HIGHLIGHT.values()))
# how far (points) hub/connector labels sit OUTSIDE each node — approx node
# radius + small gap, so the aligned text end starts at the node rim
OUTSIDE_OFFSET_PT = 3.0
# Uniform font scale for every fontsize (labels, badges, title band,
# legends). The canvas stays the same size, so >1 makes text read
# proportionally larger when the figure is scaled to page width.
# DEFAULT_FS = 1.15 is the user-approved standard (2026-10-10); the old
# 1.0 layout is still reachable via NET_FONT_SCALE=1.0. Any value other
# than DEFAULT_FS makes all output names carry a _fs{S} suffix so the
# canonical deliverable names are never silently overwritten.
DEFAULT_FS = 1.15
_FS_ENV = os.environ.get("NET_FONT_SCALE")
FS = float(_FS_ENV) if _FS_ENV else DEFAULT_FS
# NET_ROT_DEG -> rotate every network layout by this angle (degrees,
# counter-clockwise, around each panel's centroid). DEFAULT_ROT = 45 is the
# user-approved standard for the 1x4 paired deliverable (2026-10-10); the
# make_paired_2x2.py driver pins fig.ROT_DEG = 0 to keep the 2x2 layout
# unrotated. Drivers may also patch fig.FS / fig.ROT_DEG directly.
DEFAULT_ROT = 45.0
_ROT_ENV = os.environ.get("NET_ROT_DEG")
ROT_DEG = float(_ROT_ENV) if _ROT_ENV else DEFAULT_ROT
# NET_ROT_PER_SITE -> per-panel rotation overrides, "HC:0,OSCC:30" syntax;
# sites not listed fall back to ROT_DEG. Canonical mixed layout (user-chosen,
# 2026-10-10): the 1x4's HC panel keeps the UNSHARED rot0 look of the 2x2's
# HC panel, the other three panels use the 45° rotation.
DEFAULT_ROT_SITES = {"HC": 0.0}
_ROT_SITE_ENV = os.environ.get("NET_ROT_PER_SITE")
ROT_SITES = {}
if _ROT_SITE_ENV:
    for kv in _ROT_SITE_ENV.split(","):
        k, v = kv.split(":")
        ROT_SITES[k.strip().upper()] = float(v)
else:
    ROT_SITES = dict(DEFAULT_ROT_SITES)
# output-name suffix: ONLY explicitly env-overridden knobs get one, so
# canonical no-env runs always write the canonical names and cannot clobber
# them with trial variants.
SUFFIX = ("" if _FS_ENV is None else f"_fs{FS:g}") + \
         ("" if _ROT_ENV is None else f"_r{ROT_DEG:g}") + \
         ("" if os.environ.get("NET_NO_REPEL") is None else "_nr") + \
         ("" if _ROT_SITE_ENV is None else "_rotmix")
# NET_NO_REPEL -> skip the label-overlap resolver in EVERY panel (labels hug
# their nodes in the original radial placement, like the PD panel always
# does; the frame-widen pass still runs so nothing pokes outside). Trial
# knob; PD behaves the same either way.
NO_REPEL = os.environ.get("NET_NO_REPEL") == "1"
# node FILL colours mirroring circle_networks_pan_1x4.png: healthy panel all
# green; disease panels coloured by enrichment direction from
# data/pan_diease_genera.csv (+ disease-enriched red, - health-enriched blue,
# blank non-enriched grey)
HC_GREEN = "#4a9c5d"
ENR_UP, ENR_DOWN, ENR_NON = "#b2182b", "#2166ac", "#5a5a5a"
DISEASE_OF = {"PD": "Periodontal Disease", "DC": "Dental Caries",
              "OSCC": "OSCC"}
ENRICH_CSV = "data/pan_diease_genera.csv"


# ──────────────────────────────────────────────────────────────────────
# Data
# ──────────────────────────────────────────────────────────────────────

def load_mean_abundance(grp, root):
    """Mean per-sample relative abundance of every genus across ALL units of
    the group (the pooled size basis for this group's figure)."""
    sums, n_samp = {}, 0
    pattern = os.path.join(root, BIOM_DIR, f"{grp}_genus", "*.tsv")
    for fn in sorted(glob.glob(pattern)):
        df = pd.read_csv(fn, sep="\t", index_col=0)
        obs = df.to_numpy(dtype=float)         # genus x sample (counts)
        s = obs.sum(axis=0)
        s[s == 0] = np.nan
        rel = obs / s                          # per-sample relative abundance
        for i, tax in enumerate(df.index):
            v = rel[i]
            v = v[np.isfinite(v)]
            if len(v):
                sums[str(tax)] = sums.get(str(tax), 0.0) + v.sum()
        n_samp += obs.shape[1]
    if not sums or not n_samp:
        return {}
    return {g: v / n_samp for g, v in sums.items()}


def assemble(grp, root):
    """Graph + module assignment for one group network."""
    d = os.path.join(root, BACKBONE_DIR, f"{grp}_genus_uw")
    edges = pd.read_csv(os.path.join(d, "backbone_edges.csv"))
    mods = pd.read_csv(os.path.join(root, MODULARITY_DIR, f"modules_{grp}.csv"),
                       dtype={"taxon_id": str})
    mod_of = dict(zip(mods.taxon_id, mods.module))
    abu = load_mean_abundance(grp, root)

    G = nx.Graph()
    for _, e in edges.iterrows():
        G.add_edge(str(e.taxon_i), str(e.taxon_j),
                   r=float(e.pooled_obs), w=abs(float(e.pooled_obs)))
    for n in G.nodes():
        G.nodes[n]["abu"] = abu.get(n, np.nan)
        G.nodes[n]["mod"] = mod_of.get(n, 0)
        G.nodes[n]["rank"] = G.nodes[n]["mod"]

    # hub of each module: highest-degree member
    hubs = {}
    for n, d0 in G.nodes(data=True):
        m = d0["rank"]
        if m not in hubs or G.degree(n) > G.degree(hubs[m]):
            hubs[m] = n
    G.graph["hubs"] = hubs
    G.graph["site"] = grp
    return G


def load_enrichment(root):
    """Read data/pan_diease_genera.csv -> {group: {genus: '+'|'-'}}.
    Mirrors circle_network_figures_pan.py; node names in this pipeline are
    already bare genus strings (no g__ prefix), matching the CSV's Genus
    column."""
    path = os.path.normpath(os.path.join(root, os.pardir, ENRICH_CSV))
    df = pd.read_csv(path)
    out = {}
    for grp, label in DISEASE_OF.items():
        sub = df[df["Disease"] == label]
        out[grp] = {str(g): str(d) for g, d in zip(sub["Genus"], sub["Direction"])}
    return out


# ──────────────────────────────────────────────────────────────────────
# Guimerà-Amaral cartography (roles), as in the reference combined figure
# ──────────────────────────────────────────────────────────────────────

def cartography_df(G):
    """z_in (within-module degree z-score, sd=0 modules -> 0) and the
    participation coefficient P = 1 - sum_s (k_is / k_i)^2."""
    rows = []
    for m in sorted({d0["rank"] for _, d0 in G.nodes(data=True)}):
        members = [n for n in G.nodes() if G.nodes[n]["rank"] == m]
        kin = {n: sum(1 for v in G[n] if G.nodes[v]["rank"] == m)
               for n in members}
        arr = np.array([kin[n] for n in members], float)
        sd = float(arr.std())
        for n in members:
            z = (kin[n] - arr.mean()) / sd if sd > 1e-12 else 0.0
            kout = G.degree(n) - kin[n]
            rows.append(dict(taxon=n, module=m, k_in=kin[n], z_in=z,
                             k_out=kout, degree=kin[n] + kout))
    df = pd.DataFrame(rows)

    def _P(name, degree):
        if degree == 0:
            return 0.0
        shares = {}
        for v in G[name]:
            s = G.nodes[v]["rank"]
            shares[s] = shares.get(s, 0) + 1
        return 1.0 - sum((c / degree) ** 2 for c in shares.values())

    df["P"] = [_P(n, d) for n, d in zip(df.taxon, df.degree)]

    def _role(r):
        if r.z_in >= 2.5:
            if r.P < 0.3:
                return "provincial hub"
            return "connector hub" if r.P < 0.75 else "kinless hub"
        if r.P < 0.05:
            return "ultra-peripheral"
        return "peripheral" if r.P < 0.62 else "satellite connector"
    df["role"] = df.apply(_role, axis=1)
    return df


def cartography_role_nodes(G):
    """node -> outline/label colour. module hub = max within-module z (only when
    k_in varies, yellow), main connector = max P (only when cross-module links
    exist, green); hub and connector on the same node -> purple."""
    dfc = cartography_df(G)
    role_nodes = {}
    for m, sub in dfc.groupby("module"):
        sd = float(sub.k_in.std(ddof=0))
        hub = (sub.sort_values(["z_in", "degree"]).iloc[-1]
               if sd > 1e-12 else None)
        con = sub.sort_values(["P", "z_in"]).iloc[-1]
        con = con if con.P > 1e-12 else None
        if hub is not None and con is not None and \
                hub["taxon"] == con["taxon"]:
            role_nodes[hub["taxon"]] = ROLE_PURPLE
        else:
            if hub is not None:
                role_nodes[hub["taxon"]] = ROLE_RED
            if con is not None:
                role_nodes[con["taxon"]] = ROLE_BLUE
    return role_nodes


# ──────────────────────────────────────────────────────────────────────
# Fig 5a style — module-clustered network panel (no strength shading)
# ──────────────────────────────────────────────────────────────────────

def layout_modules(G, seed=7):
    """Circle layout: each module's nodes sit on a circle, ordered clockwise
    by mean relative abundance (most abundant at 12 o'clock). Module circles
    are placed on a ring around the origin, sized by member count, with a
    separation pass so circles never overlap. Deterministic."""
    from itertools import combinations
    ranks = sorted({d0["rank"] for _, d0 in G.nodes(data=True)})
    sizes = {m: sum(1 for _, d0 in G.nodes(data=True) if d0["rank"] == m)
             for m in ranks}
    smax = max(sizes.values())
    rad = {m: 0.55 + 0.65 * np.sqrt(sizes[m] / smax) for m in ranks}

    # module centres on a ring, largest module at 12 o'clock
    n_mod = len(ranks)
    order = sorted(ranks, key=lambda m: (-sizes[m], m))
    mpos = {}
    for k, m in enumerate(order):
        ang = np.pi / 2 - 2 * np.pi * k / n_mod
        mpos[m] = np.array([np.cos(ang), np.sin(ang)]) * 4.5

    # separation pass: push circle centres apart until r_i + r_j + margin
    margin = 0.55
    for _ in range(600):
        moved = False
        for a, b in combinations(ranks, 2):
            d = mpos[a] - mpos[b]
            dist = np.hypot(*d)
            need = rad[a] + rad[b] + margin
            if dist < need:
                push = (need - dist) / 2
                u_ = d / (dist or 1e-9)
                mpos[a] += u_ * push
                mpos[b] -= u_ * push
                moved = True
        if not moved:
            break

    # intra-module: nodes on the circle, clockwise by abundance
    pos = {}
    for m in ranks:
        members = [n for n, d0 in G.nodes(data=True) if d0["rank"] == m]
        abu = {n: (G.nodes[n]["abu"] if np.isfinite(G.nodes[n]["abu"]) else 0.0)
               for n in members}
        # sort by abundance descending, then clockwise from 12 o'clock
        members.sort(key=lambda n: (-abu[n], n))
        k = len(members)
        for i, n in enumerate(members):
            ang = np.pi / 2 - 2 * np.pi * i / k   # clockwise
            pos[n] = (mpos[m][0] + rad[m] * np.cos(ang),
                      mpos[m][1] + rad[m] * np.sin(ang))
    return pos, ranks


def draw_network_panel(ax, G, seed=7, smax=None, title=True, enr=None):
    """Network panel WITHOUT module-strength shading and WITHOUT the density
    badge; nodes are role-marked when G.graph['role_nodes'] is set."""
    pos, ranks = layout_modules(G, seed)
    rot = ROT_SITES.get(G.graph.get("site"), ROT_DEG)
    if rot:
        # rotate the whole layout around its centroid: edges, node circles,
        # badges and radial labels all read `pos`, so every element turns
        # together and label angles stay radial
        th = np.radians(rot)
        Rm = np.array([[np.cos(th), -np.sin(th)],
                       [np.sin(th), np.cos(th)]])
        xy = np.array(list(pos.values()))
        c0 = xy.mean(axis=0)
        rot_xy = (xy - c0) @ Rm.T + c0
        pos = {n: (rot_xy[i, 0], rot_xy[i, 1]) for i, n in enumerate(pos)}

    # edges: red positive / blue negative
    for u, v, d0 in G.edges(data=True):
        r = d0["r"]
        c = POS_C if r >= 0 else NEG_C
        ax.plot(*zip(pos[u], pos[v]), color=c, lw=0.25 + 1.6 * min(abs(r), 1),
                alpha=0.30, zorder=1, solid_capstyle="round")

    # nodes: size = mean relative abundance (smax can be a shared global basis)
    abu = np.array([d0["abu"] for _, d0 in G.nodes(data=True)])
    fin = abu[np.isfinite(abu)]
    med = float(np.median(fin)) if fin.size else 1.0
    abu_f = {n: (d0["abu"] if np.isfinite(d0["abu"]) else med)
             for n, d0 in G.nodes(data=True)}
    smax = max(abu_f.values()) if smax is None else smax
    G.graph["abu_smax"] = smax
    G.graph["abu_values"] = np.array(
        sorted(v for v in abu_f.values() if np.isfinite(v)))
    hub_nodes = set(G.graph["hubs"].values())
    # in role-marking mode the degree hubs lose their black outline
    role_nodes = G.graph.get("role_nodes") or {}
    if role_nodes:
        hub_nodes = set()
    for n, d0 in G.nodes(data=True):
        s = 8 + 60 * np.sqrt(abu_f[n] / smax)
        if n in role_nodes:            # cartography role marking (yellow/green/purple)
            ec, lw = role_nodes[n], 0.6
        else:
            hub = n in hub_nodes
            ec, lw = ("#111" if hub else "white"), (1.2 if hub else 0.7)
        # node fill: healthy panel all green; disease panels ALL coloured by
        # enrichment direction (red +, blue -, grey otherwise — including
        # genera absent from the CSV, treated as non-enriched)
        if enr is None:
            fill = HC_GREEN
        else:
            d = enr.get(n)
            fill = (ENR_UP if d == "+" else
                    ENR_DOWN if d == "-" else ENR_NON)
        # emphasis style per hl_mode: "font" = magenta bold LABEL only
        # (no ring); "plain" = no magenta anywhere, key genera appear only
        # as plain black labels (their names go in the figure caption)
        ax.scatter(*pos[n], s=s, color=fill,
                   edgecolor=ec, linewidth=lw, zorder=3)

    # module badges: M# only, centred inside the module circle
    import matplotlib.patheffects as pe
    stroke = [pe.withStroke(linewidth=1.8, foreground="white")]
    all_xy = np.array(list(pos.values()))
    span = max(all_xy.max(axis=0) - all_xy.min(axis=0))
    for m in ranks:
        members = [n for n, d0 in G.nodes(data=True) if d0["rank"] == m]
        cx = float(np.mean([pos[n][0] for n in members]))
        cy = float(np.mean([pos[n][1] for n in members]))
        ax.text(cx, cy, f"M{m}", ha="center",
                va="center", fontsize=8 * FS, fontweight="bold",
                color="#000000",
                zorder=6, path_effects=stroke)

    # internal padding between the network and the panel frame — generous,
    # especially vertically, so labels never hug the border (user request)
    ax.set_xlim(all_xy[:, 0].min() - 0.14 * span,
                all_xy[:, 0].max() + 0.14 * span)
    ax.set_ylim(all_xy[:, 1].min() - 0.22 * span,
                all_xy[:, 1].max() + 0.22 * span)

    ax.set_aspect("equal")
    place_hub_labels(ax, G, pos)
    G.graph["pos"] = pos

    ax.set_facecolor("white")
    ax.set_xticks([]); ax.set_yticks([])
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_color("#9a9a9a"); sp.set_linewidth(0.8)
    ax.grid(True, color="#e3e3e2", lw=0.55)
    ax.set_axisbelow(True)
    if title:
        ax.set_title(GRP_LABEL.get(G.graph["site"], G.graph["site"]),
                     fontsize=11 * FS)


def place_hub_labels(ax, G, pos):
    """Label each module hub (or the cartography roles in role-marking mode)
    with its genus name OUTSIDE the module circle, placed along the true
    screen-space radial direction (centre -> node); white halo keeps it
    legible. Labels are ordered for the downstream greedy placement: by module
    number ascending (M1 -> M2 -> ...), and within a module clockwise from
    12 o'clock (screen-space polar angle from +y)."""
    import matplotlib.patheffects as pe
    stroke = [pe.withStroke(linewidth=1.0, foreground="white")]
    role_nodes = G.graph.get("role_nodes")
    if role_nodes:
        items = list(role_nodes.items())          # (node, role colour)
    else:
        items = [(hn, "#222") for hn in G.graph["hubs"].values()]
    labelled = {hn for hn, _ in items}
    # highlighted key genera get an (extra) label; style depends on hl_mode:
    #   "font" -> magenta bold enlarged label (role-labelled keep role colour)
    #   "ring" -> plain black normal label; the halo ring carries the emphasis
    # (black labels for the pan-disease CSV genera were removed per decision)
    hl_nodes = set(G.graph.get("highlight") or ())
    font_mode = G.graph.get("hl_mode", "font") == "font"
    for n in hl_nodes - labelled:
        # font version: bold enlarged but BLACK (magenta retired per user);
        # plain version: same black, normal size handled below
        items.append((n, "#000000"))
        labelled.add(n)
    # module centre per module (screen px) for polar-angle sort
    mod_centre = {}
    for m in {G.nodes[n]["rank"] for n, _ in items}:
        members = [n for n, d0 in G.nodes(data=True) if d0["rank"] == m]
        mod_centre[m] = (float(np.mean([pos[n][0] for n in members])),
                         float(np.mean([pos[n][1] for n in members])))
    items.sort(key=lambda it: (
        G.nodes[it[0]]["rank"],
        # clockwise angle from +y (12 o'clock): 0 at top, +pi/2 at right
        # screen-py grows downward, so use -(dy) -> atan2(dx, -dy) gives
        # 0 at top, increasing clockwise
        (lambda c, p: np.arctan2(p[0] - c[0], -(p[1] - c[1])))(
            mod_centre[G.nodes[it[0]]["rank"]], pos[it[0]])))
    labels = []
    for hn, c in items:
        m = G.nodes[hn]["rank"]
        cx, cy = mod_centre[m]
        # radial direction in screen space (compensates for non-square data
        # unit scaling so the offset truly points outward from the module)
        disp = ax.transData.transform([(cx, cy), pos[hn]])
        sdx, sdy = disp[1] - disp[0]
        snrm = np.hypot(sdx, sdy) or 1.0
        ux, uy = sdx / snrm, sdy / snrm
        # rotation along the radial direction (angle in screen space);
        # right-side nodes: left end just outside the node rim (text grows
        # outward); left-side nodes: right end at the rim (+180° so it reads
        # outward). OUTSIDE_OFFSET_PT ≈ node radius + small gap in points.
        ang = np.degrees(np.arctan2(uy, ux))
        if ux < 0:
            rot, ha = ang + 180, "right"
        else:
            rot, ha = ang, "left"
        big = hn in hl_nodes and font_mode   # highlighted -> bold + enlarged
        # plain version: hub/connector labels KEEP their role colours
        # (user decision) — only the magenta emphasis is retired
        ann = ax.annotate(
            hn.replace("g__", ""), pos[hn],
            xytext=(ux * OUTSIDE_OFFSET_PT, uy * OUTSIDE_OFFSET_PT),
            textcoords="offset points",
            ha=ha, va="center", rotation=rot,
            rotation_mode="anchor",
            fontsize=(5.2 if big else 4) * FS,
            fontweight="bold" if big else "normal",
            color=c, zorder=7, path_effects=stroke)
        # radial unit vector (display space), reused by the overlap resolver:
        # labels can additionally be pushed FARTHER OUT along it
        ann.label_dir = (ux, uy)
        labels.append(ann)
    G.graph["label_anns"] = labels
    G.graph["badge_texts"] = [t for t in ax.texts if t not in set(labels)]


def resolve_label_overlaps(fig, graphs, widen=True, clamp=False):
    """Keep text off text after the final rendered layout is fixed.

    Labels (hub/connector/enrichment) are anchored to their node at a fixed
    radial offset; instead of translating them (which detaches the label
    from its node), each label is ROTATED about its node anchor until it no
    longer overlaps neighbours — the text stays on the circle rim, just at
    a slightly different angle. Module badges (M#) are translated freely.
    With widen=True the axis limits are expanded so no text pokes outside
    its panel. With clamp=True the rotation search additionally REJECTS
    angles whose text would sit outside the panel frame — used for the
    final pass (after equalize_panel_limits re-set the limits), so the text
    can never escape the panel even though the limits no longer move."""
    import math
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    px = 72.0 / fig.dpi      # window-extent px per offset point

    def tbox(t):
        return mtext.Text.get_window_extent(t, renderer)
    for G in graphs:
        anns = G.graph.get("label_anns", [])
        static = G.graph.get("badge_texts", [])
        ax = anns[0].axes if anns else None
        if ax is None:
            continue
        abox = ax.get_window_extent()  # panel frame, display coords
        inv = ax.transData.inverted()

        m_frame = 2.0                       # keep-out margin from the frame (px)

        # PD panel: NO label avoidance at all (user choice) — text stays
        # right next to its node, in the original radial placement. Only the
        # panel limits still expand (widen pass) so nothing pokes outside
        # the frame. NET_NO_REPEL extends the same treatment to every panel.
        if G.graph.get("site") == "PD" or NO_REPEL:
            if widen:
                x0, x1 = ax.get_xlim()
                y0, y1 = ax.get_ylim()
                for t in anns + static:
                    bb = tbox(t)
                    (bx0, by0), (bx1, by1) = ax.transData.inverted().transform(
                        [[bb.x0, bb.y0], [bb.x1, bb.y1]])
                    x0, x1 = min(x0, bx0), max(x1, bx1)
                    y0, y1 = min(y0, by0), max(y1, by1)
                ax.set_xlim(x0, x1)
                ax.set_ylim(y0, y1)
            continue

        def out_of_frame(b):
            return (b.x0 < abox.x0 + m_frame or b.x1 > abox.x1 - m_frame or
                    b.y0 < abox.y0 + m_frame or b.y1 > abox.y1 - m_frame)

        def shift(t, dxp, dyp):
            x, y = t.get_position()
            cxp, cyp = ax.transData.transform((x, y))
            t.set_position(tuple(inv.transform([(cxp + dxp,
                                                 cyp + dyp)])[0]))

        # per-label base rotation (degrees) and node anchor (display px)
        base_rot = {}
        anchor_px = {}
        for a in anns:
            base_rot[id(a)] = a.get_rotation()
            anchor_px[id(a)] = ax.transData.transform(a.xy)

        # upright envelope: effective base rotation folded into (-180, 180],
        # then into (-90, 90] if it points "down" (text mirrored). Deltas are
        # only allowed that keep base_eff + delta inside [-90, +90], so no
        # label can ever end up upside-down after avoidance.
        ROT_LIM = 90.0

        def base_eff_of(a):
            r = ((base_rot[id(a)] + 180.0) % 360.0) - 180.0
            if r > ROT_LIM:
                r -= 360.0
            return r

        base_eff = {id(a): base_eff_of(a) for a in anns}

        def upright_cands(a, cands):
            return [d for d in cands
                    if abs(base_eff[id(a)] + d) <= ROT_LIM + 1e-9]

        # candidate radial-offset multipliers (multiples of the original
        # OUTSIDE_OFFSET_PT): crowded rims escape by pushing labels farther
        # out along their radial direction, in addition to rotating. A fine
        # ladder (many steps) lets several labels share one inter-module gap
        # at staggered distances instead of all-or-nothing escapes.
        OFF_MULTS = (1.0, 1.6, 2.3, 3.1, 4.0, 5.2, 7.0, 9.5)

        def set_label_rot(a, delta_deg):
            a.set_rotation(base_rot[id(a)] + delta_deg)

        def set_label_off(a, mult):
            d = getattr(a, "label_dir", None)
            if d is not None:
                ux, uy = d
                a.set_position((ux * OUTSIDE_OFFSET_PT * mult,
                                uy * OUTSIDE_OFFSET_PT * mult))

        def rotate_relax(rounds, max_deg=55.0, step=2.5):
            """Greedy sequential placement: labels are processed in order; each
            earlier label is FROZEN once placed, and each later label is rotated
            (trying both directions) AND/OR pushed farther out radially until
            it no longer overlaps any already-placed label, then frozen itself.
            Avoids the many-neighbour cancellation that stranded mid-sequence
            labels. A single renderer snapshot (drawn once) is reused for all
            probes — moving one label doesn't move any other element, so no
            re-render needed."""
            placed = []
            placed_boxes = []
            n_step = int(round(max_deg / step))
            # candidate angles: 0, then +step, +2step, ... +max, then
            # -step, -2step, ... -max. Sequential labels thus fan out in the
            # SAME direction with increasing angle (comb-like), which is what
            # actually separates adjacent labels on a circle. Combined with
            # the radial-offset multipliers; smallest angle + offset wins.
            cands = ([0.0]
                     + [step * i for i in range(1, n_step + 1)]
                     + [-step * i for i in range(1, n_step + 1)])
            # per-label combos (upright-filtered) built inside the loop
            combos_of = {id(a): [(d, m) for d in upright_cands(a, cands)
                                 for m in OFF_MULTS] for a in anns}
            for a in anns:
                # probe EVERY (angle, offset) candidate with the cached
                # renderer; score (escapes frame, #overlaps, |delta| + 0.5*
                # offset step) so a candidate that stays inside always beats
                # one that pokes out, rotation is preferred over pushing, and
                # close-to-the-rim wins otherwise
                best = None
                for d, m in combos_of[id(a)]:
                    set_label_rot(a, d)
                    set_label_off(a, m)
                    ba = tbox(a)
                    viol = 1 if (clamp and out_of_frame(ba)) else 0
                    n = sum(1 for pb in placed_boxes if ba.overlaps(pb))
                    key = (viol, n, abs(d) + 0.5 * OFF_MULTS.index(m))
                    if best is None or key < best[0]:
                        best = (key, d, m)
                    if viol == 0 and n == 0:
                        break
                set_label_rot(a, best[1])
                set_label_off(a, best[2])
                a.off_mult = best[2]
                placed.append(a)
                # cache this label's box at its chosen (angle, offset)
                placed_boxes.append(tbox(a))

        # badges still translated freely; labels only rotated
        def badge_relax(rounds, cap=60.0):
            bitems = [("badge", t) for t in static]
            if not bitems:
                return
            # fontscale runs only: longer labels can run straight across a
            # badge, so badges are ALSO repelled from overlapping genus
            # labels (FS == 1.0 keeps the frozen legacy layout untouched)
            pushers = list(anns) if FS != 1.0 else []
            start = {}
            for kind, obj in bitems:
                bb = tbox(obj)
                start[id(obj)] = (0.5 * (bb.x0 + bb.x1),
                                  0.5 * (bb.y0 + bb.y1))

            def used_of(obj):
                bb = tbox(obj)
                return np.hypot(0.5 * (bb.x0 + bb.x1) - start[id(obj)][0],
                                0.5 * (bb.y0 + bb.y1) - start[id(obj)][1])
            for _ in range(rounds):
                moved_any = False
                for _, obj in bitems:
                    ba = tbox(obj)
                    fx = fy = 0.0
                    for other in bitems:
                        if other[1] is obj:
                            continue
                        sb = tbox(other[1])
                        if not ba.overlaps(sb):
                            continue
                        cx = 0.5 * (sb.x0 + sb.x1)
                        cy = 0.5 * (sb.y0 + sb.y1)
                        bx = 0.5 * (ba.x0 + ba.x1)
                        by = 0.5 * (ba.y0 + ba.y1)
                        vx, vy = bx - cx, by - cy
                        if vx == 0 and vy == 0:
                            vx, vy = 1.0, 0.0
                        nrm = np.hypot(vx, vy)
                        pen = max(min(ba.x1 - sb.x0, sb.x1 - ba.x0),
                                  min(ba.y1 - sb.y0, sb.y1 - ba.y0))
                        fx += vx / nrm * (pen + 6)
                        fy += vy / nrm * (pen + 6)
                    for src in pushers:
                        sa = tbox(src)
                        if not ba.overlaps(sa):
                            continue
                        cx = 0.5 * (sa.x0 + sa.x1)
                        cy = 0.5 * (sa.y0 + sa.y1)
                        bx = 0.5 * (ba.x0 + ba.x1)
                        by = 0.5 * (ba.y0 + ba.y1)
                        vx, vy = bx - cx, by - cy
                        if vx == 0 and vy == 0:
                            vx, vy = 1.0, 0.0
                        nrm = np.hypot(vx, vy)
                        pen = max(min(ba.x1 - sa.x0, sa.x1 - ba.x0),
                                  min(ba.y1 - sa.y0, sa.y1 - ba.y0))
                        fx += vx / nrm * (pen + 6)
                        fy += vy / nrm * (pen + 6)
                    if fx or fy:
                        room = cap - used_of(obj)
                        if room <= 0:
                            continue
                        if np.hypot(fx, fy) > room:
                            sc = room / np.hypot(fx, fy)
                            fx, fy = fx * sc, fy * sc
                        shift(obj, fx, fy)
                        moved_any = True
                if not moved_any:
                    break

        rotate_relax(60, max_deg=90.0, step=2.5)

        def resolve_residual(rounds=20, max_deg=90.0, step=2.5):
            """Second pass: find label pairs still overlapping after greedy
            placement (typically cross-module). For each stuck label, keep its
            first overlapping partner fixed and search (rotation, offset) for
            the least-colliding spot: first a move that clears partner AND
            everyone else, otherwise the one that strictly reduces the total
            overlap count (no exact clear exists)."""
            n_step = int(round(max_deg / step))
            cands = ([0.0]
                     + [step * i for i in range(1, n_step + 1)]
                     + [-step * i for i in range(1, n_step + 1)])
            cands_of = {id(a): upright_cands(a, cands) for a in anns}
            others_cache = {id(a): [o for o in anns if o is not a]
                            for a in anns}
            for _ in range(rounds):
                fig.canvas.draw()
                boxes = {id(a): tbox(a) for a in anns}
                moved = False
                for i, a in enumerate(anns):
                    ba = boxes[id(a)]
                    # is a overlapping anything?
                    partner = None
                    for j, b in enumerate(anns):
                        if j == i:
                            continue
                        if ba.overlaps(boxes[id(b)]):
                            partner = b
                            break
                    if partner is None:
                        continue
                    # score of the CURRENT spot, so a candidate must strictly
                    # beat it (the search may end local-min-bound otherwise)
                    base = base_rot[id(a)]
                    cur_delta = a.get_rotation() - base
                    cur_mult = getattr(a, "off_mult", 1.0)
                    oth_a = others_cache[id(a)]
                    pb = boxes[id(partner)]

                    def spot(box):
                        hit = 1 if box.overlaps(pb) else 0
                        n_new = sum(1 for o in oth_a
                                    if box.overlaps(boxes[id(o)]) and
                                    o is not partner)
                        return (1 if (clamp and out_of_frame(box)) else 0,
                                hit, n_new)
                    cur_key = spot(ba) + (0.0, 0)
                    pbest = None
                    for d in cands_of[id(a)]:
                        for m in OFF_MULTS:
                            if abs(d - cur_delta) < 1e-9 and m == cur_mult:
                                continue
                            a.set_rotation(base + d)
                            set_label_off(a, m)
                            na = tbox(a)
                            k3 = spot(na)
                            key = k3 + (abs(d), OFF_MULTS.index(m))
                            if pbest is None or key < pbest[0]:
                                pbest = (key, (d, m))
                            if k3 == (0, 0, 0):
                                break
                        else:
                            continue
                        break
                    # restore current placement, then move only if the best
                    # candidate strictly improves on the current spot
                    a.set_rotation(base + cur_delta)
                    set_label_off(a, cur_mult)
                    if pbest is not None and pbest[0][:3] < cur_key[:3]:
                        a.set_rotation(base + pbest[1][0])
                        set_label_off(a, pbest[1][1])
                        a.off_mult = pbest[1][1]
                        moved = True
                if not moved:
                    break
                fig.canvas.draw()

        resolve_residual(rounds=20)
        badge_relax(80)
        for _ in range(3):
            moved = False
            for t in static:
                bb = tbox(t)
                dxp = dyp = 0.0
                if bb.x0 < abox.x0 + 3:
                    dxp = abox.x0 + 3 - bb.x0
                elif bb.x1 > abox.x1 - 3:
                    dxp = abox.x1 - 3 - bb.x1
                if bb.y0 < abox.y0 + 3:
                    dyp = abox.y0 + 3 - bb.y0
                elif bb.y1 > abox.y1 - 3:
                    dyp = abox.y1 - 3 - bb.y1
                if dxp or dyp:
                    shift(t, dxp, dyp)
                    moved = True
            if not moved:
                break
            badge_relax(40)
        if os.environ.get("NET_DEBUG"):
            left = []
            for a in anns:
                ba = tbox(a)
                for b in anns:
                    if b is not a and ba.overlaps(tbox(b)):
                        left.append(("label", "label"))
                        break
            print(f"[debug] residual label overlaps (widen={widen}):",
                  len(left))
        if not widen:
            continue
        x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
        for t in ax.texts:
            bb = tbox(t)
            (bx0, by0), (bx1, by1) = ax.transData.inverted().transform(
                [[bb.x0, bb.y0], [bb.x1, bb.y1]])
            x0, x1 = min(x0, bx0), max(x1, bx1)
            y0, y1 = min(y0, by0), max(y1, by1)
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)


# ──────────────────────────────────────────────────────────────────────
# Facet-strip titles + shared bottom legend (reference-figure machinery)
# ──────────────────────────────────────────────────────────────────────

def add_title_band(fig, ax, text):
    """Grey facet-strip title band above one panel."""
    x0, y0, x1, y1 = ax.get_subplotspec().get_position(fig).extents
    bh = 0.024
    fig.add_artist(Rectangle((x0, y1 + 0.004), x1 - x0, bh,
                             facecolor="#dcdcdc", edgecolor="none",
                             transform=fig.transFigure, zorder=10))
    fig.text((x0 + x1) / 2, y1 + 0.004 + bh / 2, text, ha="center",
             va="center", fontsize=10.5, fontweight="bold", zorder=11)


def equalize_panel_limits(fig, axes):
    """Give every panel the same drawn box size: expand each panel's data
    limits so its aspect matches its subplot cell."""
    for ax in axes:
        sp = ax.get_subplotspec().get_position(fig)
        fw, fh = fig.get_size_inches()
        r = (sp.width * fw) / (sp.height * fh)
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        xs, ys = x1 - x0, y1 - y0
        if xs / ys < r:
            need = ys * r
            c = 0.5 * (x0 + x1)
            ax.set_xlim(c - need / 2, c + need / 2)
        else:
            need = xs / r
            c = 0.5 * (y0 + y1)
            ax.set_ylim(c - need / 2, c + need / 2)


def add_shared_legend(fig, data, slot, width=0.50):
    """Compact dashed-box legend BELOW the panels, mirroring the layout of
    circle_networks_pan_1x4.png: four measured, tightly-spaced sections —
    Correlation, Role, Enrichment Status, Mean Relative Abundance (global
    basis). Section widths are measured from the renderer so inter-section
    gaps stay small."""
    leg = fig.add_subplot(slot)
    pos = leg.get_position()
    leg.set_position([0.5 - width / 2, pos.y0, width, pos.height])
    leg.set_xlim(0, 1); leg.set_ylim(0, 1)
    leg.axis("off")
    leg.set_zorder(0)
    fs = 8 * FS
    yhead, ybody = 0.78, 0.34
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    dpi = fig.dpi
    fig_w_in = fig.get_size_inches()[0]

    def tw(s, size, bold=False):
        t = leg.text(0, 0, s, fontsize=size,
                     fontweight="bold" if bold else "normal")
        w = t.get_window_extent(r).width / dpi
        t.remove()
        return w

    def dia(s):
        return 2.0 * np.sqrt(s / np.pi) / 72.0

    # section content widths (inches)
    LINE = 0.11
    GAP, PAD = 0.10, 0.12
    # Correlation Type: Negative/Positive stacked vertically -> width is the
    # single widest line entry, not the sum of both
    w1 = max(tw("Correlation Type", fs, True),
             LINE + 0.05 + max(tw("Negative", fs - 1.5),
                                tw("Positive", fs - 1.5)))
    w2 = max(tw("Role", fs, True),
             dia(64) + 0.05 + max(tw("Module hub", fs - 1.5),
                                  tw("Main connector", fs - 1.5),
                                  tw("Hub + connector", fs - 1.5)))
    enr_labels = ["Disease-enriched", "Health-enriched", "Non-enriched"]
    w3 = max(tw("Enrichment Status", fs, True),
             dia(70) + 0.06 + max(tw(l, fs - 1.5) for l in enr_labels))
    allv = np.concatenate([G.graph["abu_values"] for G in data.values()])
    smax = next(iter(data.values())).graph["abu_smax"]
    v_hi = float(np.max(allv))
    v_lo = round(float(np.median(allv)), 3)
    s_hi = 8 + 60 * np.sqrt(v_hi / smax)
    s_lo = 8 + 60 * np.sqrt(v_lo / smax)
    row4 = (dia(s_hi) + 0.05 + tw(f"{v_hi:.0%}", fs - 1.5) + 0.16 +
            dia(s_lo) + 0.05 + tw(f"{v_lo:.1%}", fs - 1.5))
    w4 = max(tw("Mean Relative Abundance", fs, True), row4)

    total = PAD + w1 + GAP + w2 + GAP + w3 + GAP + w4 + PAD
    width_frac = total / fig_w_in
    pos = leg.get_position()
    leg.set_position([0.5 - width_frac / 2, pos.y0, width_frac, pos.height])

    def ax(xin):
        return xin / total

    leg.add_patch(Rectangle((0.004, 0.03), 0.988, 0.94, fill=False,
                            ec="#444", lw=0.9))
    yh, yb = 0.78, 0.34

    # ── Correlation Type (Negative / Positive stacked vertically) ──
    x = PAD
    leg.text(ax(x), yh, "Correlation Type", fontsize=fs, fontweight="bold",
             va="center")
    leg.plot([ax(x), ax(x + LINE)], [0.55, 0.55], color=NEG_C, lw=1.8,
             solid_capstyle="round")
    leg.text(ax(x + LINE + 0.05), 0.55, "Negative", fontsize=fs - 1.5,
             va="center")
    leg.plot([ax(x), ax(x + LINE)], [0.36, 0.36], color=POS_C, lw=1.8,
             solid_capstyle="round")
    leg.text(ax(x + LINE + 0.05), 0.36, "Positive", fontsize=fs - 1.5,
             va="center")

    # ── Role ──
    x = PAD + w1 + GAP
    leg.text(ax(x), yh, "Role", fontsize=fs, fontweight="bold", va="center")
    role_labels = ["Module hub", "Main connector", "Hub + connector"]
    role_cols = [ROLE_RED, ROLE_BLUE, ROLE_PURPLE]
    # "Highlighted genus" legend row removed (font mode black per user);
    # rows stay at the top three slot heights
    for yy, col, lab in zip((0.60, 0.44, 0.28),
                            role_cols, role_labels):
        leg.scatter([ax(x + dia(64) / 2)], [yy], s=64,
                    facecolor="#dddddd", edgecolor=col,
                    linewidth=1.2, zorder=5)
        leg.text(ax(x + dia(64) + 0.05), yy, lab, fontsize=fs - 1.5,
                 va="center")

    # ── Enrichment Status ──
    x = PAD + w1 + GAP + w2 + GAP
    leg.text(ax(x), yh, "Enrichment Status", fontsize=fs,
             fontweight="bold", va="center")
    for yy, col, lab in zip((0.55, 0.36, 0.17),
                            (ENR_UP, ENR_DOWN, ENR_NON), enr_labels):
        leg.scatter([ax(x + dia(70) / 2)], [yy], s=70, facecolor=col,
                    edgecolor="none", zorder=5)
        leg.text(ax(x + dia(70) + 0.06), yy, lab, fontsize=fs - 1.5,
                 va="center")

    # ── Mean Relative Abundance ──
    x = PAD + w1 + GAP + w2 + GAP + w3 + GAP
    leg.text(ax(x), yh, "Mean Relative Abundance", fontsize=fs,
             fontweight="bold", va="center")
    cx = x + dia(s_hi) / 2
    leg.scatter([ax(cx)], [yb], s=s_hi, color="black", zorder=5)
    leg.text(ax(x + dia(s_hi) + 0.05), yb, f"{v_hi:.0%}",
             fontsize=fs - 1.5, va="center")
    x3 = x + dia(s_hi) + 0.05 + tw(f"{v_hi:.0%}", fs - 1.5) + 0.16
    leg.scatter([ax(x3 + dia(s_lo) / 2)], [yb], s=s_lo, color="black",
                zorder=5)
    leg.text(ax(x3 + dia(s_lo) + 0.05), yb, f"{v_lo:.1%}",
             fontsize=fs - 1.5, va="center")
    return leg


def snug_shared_legend(fig, leg, panel_axes, gap=0.006):
    """Move the legend axes directly below the bottom row of panels."""
    fig.canvas.draw()
    bottom = min(ax.get_position().y0 for ax in panel_axes)
    pos = leg.get_position()
    leg.set_position([pos.x0, bottom - gap - pos.height,
                      pos.width, pos.height])


# ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", default=".")
    ap.add_argument("--groups", default=",".join(GROUPS),
                    help="comma-separated subset of groups to draw "
                         "(default: all four)")
    ap.add_argument("--layouts", default="1x4",
                    help="additional/different combined-figure arrangement(s) "
                         "as comma-separated RxC specs, e.g. '2x2' or "
                         "'1x4,2x2'; the output name follows "
                         "module_networks_{R}x{C}_{mode}.png (default 1x4)")
    args = ap.parse_args()
    groups = [g.strip().upper() for g in args.groups.split(",") if g.strip()]
    bad = [g for g in groups if g not in GROUPS]
    if bad:
        raise SystemExit(f"unknown groups: {bad} (choose from {GROUPS})")
    layouts = [s.strip().lower() for s in args.layouts.split(",") if s.strip()]
    root = os.path.join(args.root, ROOT_DIR)
    out_base = os.path.join(root, "module_networks")
    os.makedirs(out_base, exist_ok=True)

    graphs = {}
    for grp in groups:
        mods_path = os.path.join(root, MODULARITY_DIR, f"modules_{grp}.csv")
        edges_path = os.path.join(root, BACKBONE_DIR, f"{grp}_genus_uw",
                                  "backbone_edges.csv")
        if not (os.path.exists(mods_path) and os.path.exists(edges_path)):
            print(f"[net] {grp}: missing backbone/modules output, skipped")
            continue
        G = assemble(grp, root)
        G.graph["highlight"] = [n for n in HIGHLIGHT.get(grp, []) if n in G]
        ranks = sorted({d0["rank"] for _, d0 in G.nodes(data=True)})
        print(f"[net] {grp}: {G.number_of_nodes()} nodes, "
              f"{G.number_of_edges()} edges, {len(ranks)} modules")
        graphs[grp] = G

    if not graphs:
        print("[net] nothing to draw")
        return

    # cartography roles for every panel
    for grp, G in graphs.items():
        G.graph["role_nodes"] = cartography_role_nodes(G)
        print(f"[roles] {grp}: " + ", ".join(
            f"{n.replace('g__', '')}({c[1:]})"
            for n, c in G.graph["role_nodes"].items()))
        print(f"[hl] {grp}: {len(G.graph['highlight'])} highlighted -> "
              + ", ".join(G.graph["highlight"]))

    # one global abundance basis so node sizes are comparable across panels
    gmax = max(d0["abu"] for G in graphs.values() for _, d0 in G.nodes(data=True)
               if np.isfinite(d0["abu"]))
    print(f"[net] global max mean rel. abundance: {gmax:.4f}")

    # enrichment directions for node fills + extra black labels (mirrors
    # circle_networks_pan_1x4.png); HC panel stays all-green (enr=None)
    enrich = load_enrichment(root)
    for grp in DISEASE_OF:
        if grp in enrich:
            e = enrich[grp]
            up = sum(1 for n in e.values() if n == "+")
            dn = sum(1 for n in e.values() if n == "-")
            print(f"[enrich] {grp}: disease+={up} health-={dn}")

    # ── two emphasis versions of the deliverable (plus singles):
    #    "font"  = key genera marked by magenta bold label only, no ring
    #    "plain" = no magenta anywhere: key genera as plain black labels,
    #              hub/connector labels keep their role colours
    order = [g for g in GROUPS if g in graphs]
    for mode in ("font", "plain"):
        for G in graphs.values():
            G.graph["hl_mode"] = mode
        print(f"[mode] {mode}")

        # combined deliverable(s): each layout spec "RxC" gives its panels-
        # per-row count; shared legend sits BELOW the bottom row
        for spec in layouts:
            make_combined(mode, order, graphs, gmax, enrich, out_base,
                          ncol=int(spec.split("x")[-1]))

        # standalone panel per group (own legend, role-marked as well)
        for grp in order:
            G = graphs[grp]
            fig, ax = plt.subplots(figsize=(9.6, 8.6))
            enr = enrich.get(grp) if grp != "HC" else None
            draw_network_panel(ax, G, enr=enr)
            _add_single_legend(ax, G)
            fig.tight_layout()
            resolve_label_overlaps(fig, [G])
            fname = f"module_network_{grp}_genus_{mode}{SUFFIX}.png"
            fig.savefig(os.path.join(out_base, fname), dpi=300)
            plt.close(fig)
            print(f"[net] {grp} -> {fname}")


def make_combined(mode, order, graphs, gmax, enrich, out_base, ncol):
    """One combined deliverable: `ncol` panels per row
    (ceil(len(order)/ncol) rows), shared legend BELOW the bottom row.
    Single-row grids keep a legend band under the panels; multi-row grids
    stack the panel rows with ZERO vertical gap between them (user choice),
    with a single legend row at the very bottom. Writes
    module_networks_{nrow}x{ncol}_{mode}.{png,pdf} into out_base."""
    nrow = int(np.ceil(len(order) / ncol))
    # multi-row grids stack the rows with ZERO vertical gap between them
    fig = plt.figure(figsize=(5.2 * ncol, 6.2 * nrow))
    if nrow == 1:
        gs = fig.add_gridspec(2, ncol, height_ratios=[1, 0.18],
                              hspace=0.22, wspace=0.05)
    else:
        gs = fig.add_gridspec(nrow + 1, ncol,
                              height_ratios=[1] * nrow + [0.18],
                              hspace=0.0, wspace=0.05)
    panel_axes = []
    for k, grp in enumerate(order):
        r, c = divmod(k, ncol)
        ax = fig.add_subplot(gs[r, c])
        enr = enrich.get(grp) if grp != "HC" else None
        draw_network_panel(ax, graphs[grp], smax=gmax, title=False, enr=enr)
        ax.set_anchor("N")
        add_title_band(fig, ax, GRP_LABEL[grp])
        panel_axes.append(ax)
    leg = add_shared_legend(fig, graphs, gs[nrow, :], width=0.50)
    resolve_label_overlaps(fig, [graphs[g] for g in order])
    equalize_panel_limits(fig, panel_axes)
    resolve_label_overlaps(fig, [graphs[g] for g in order], widen=False,
                           clamp=True)
    snug_shared_legend(fig, leg, panel_axes)
    for ext in ("png", "pdf"):      # pdf = vector version for publication
        base = f"module_networks_{nrow}x{ncol}_{mode}{SUFFIX}"
        out = os.path.join(out_base, f"{base}.{ext}")
        fig.savefig(out, dpi=300, bbox_inches="tight", pad_inches=0.10)
        print(f"[net] -> {out}")
    plt.close(fig)


def _add_single_legend(ax, G):
    """Standalone-panel legend: correlation + roles (incl. the highlighted
    key genera) + this panel's abundance. Sections are laid out on absolute
    inset-fraction rows so everything stays inside the legend frame."""
    from matplotlib.patches import Rectangle as _R
    w, h = 0.34, 0.22
    leg = ax.inset_axes([0.01, 0.02, w, h])
    leg.set_xlim(0, 1); leg.set_ylim(0, 1)
    leg.set_xticks([]); leg.set_yticks([])
    leg.set_facecolor("white"); leg.set_zorder(10)
    for sp in leg.spines.values():
        sp.set_visible(False)
    leg.add_patch(_R((0.005, 0.005), 0.99, 0.99, fill=False,
                     ec="#444", lw=0.8,
                     transform=leg.transAxes, zorder=10))
    fs = 6.3 * FS

    # ── Correlation ──
    y = 0.86
    leg.text(0.07, y, "Correlation", fontsize=fs, fontweight="bold", va="center")
    yl = 0.70
    leg.plot([0.10, 0.22], [yl, yl], color=NEG_C, lw=1.4, solid_capstyle="round")
    leg.text(0.25, yl, "Negative", fontsize=fs - 0.8, va="center")
    leg.plot([0.58, 0.70], [yl, yl], color=POS_C, lw=1.4, solid_capstyle="round")
    leg.text(0.73, yl, "Positive", fontsize=fs - 0.8, va="center")

    # ── Role (hub / connector / hub+connector; highlight row removed) ──
    y = 0.51
    leg.text(0.07, y, "Role", fontsize=fs, fontweight="bold", va="center")
    yc = 0.34
    for xc, col, lab in ((0.13, ROLE_RED, "Hub"),
                         (0.44, ROLE_BLUE, "Conn."),
                         (0.72, ROLE_PURPLE, "Hub+conn.")):
        leg.scatter([xc], [yc], s=40, facecolor="#dddddd",
                    edgecolor=col, linewidth=1.5, zorder=11)
        leg.text(xc + 0.03, yc, lab, fontsize=fs - 1.2, va="center")

    # ── Mean Relative Abundance ──
    y2 = 0.17
    leg.text(0.07, y2, "Mean Relative Abundance", fontsize=fs,
             fontweight="bold", va="center")
    smax = G.graph["abu_smax"]
    av = G.graph["abu_values"]
    v_hi = float(av.max())
    v_lo = float(np.median(av))
    yc2 = 0.045
    leg.scatter([0.15], [yc2], s=8 + 60 * np.sqrt(v_hi / smax),
                color="black", zorder=11)
    leg.text(0.32, yc2, f"{v_hi:.1%}", fontsize=fs - 0.8, va="center")
    leg.scatter([0.55], [yc2], s=8 + 60 * np.sqrt(v_lo / smax),
                color="black", zorder=11)
    leg.text(0.68, yc2, f"{v_lo:.2%}", fontsize=fs - 0.8, va="center")


if __name__ == "__main__":
    main()