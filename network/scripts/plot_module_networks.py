"""Render the 2x2 module networks with highlighted genus labels."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from network_layout import FS, layout_modules, place_hub_labels, resolve_label_overlaps, equalize_panel_limits
from network_common import (MODULE_DIR, GROUPS, GRP_LABEL,
                            ROLE_RED, ROLE_BLUE, ROLE_PURPLE)

POS_C, NEG_C = "#d64541", "#3b74c4"
HC_GREEN = "#4a9c5d"
ENR_UP, ENR_DOWN, ENR_NON = "#b2182b", "#2166ac", "#5a5a5a"
DISEASE_OF = {"PD": "Periodontal Disease", "DC": "Dental Caries", "OSCC": "OSCC"}
HIGHLIGHT = {
    "PD": ["Fretibacterium", "Treponema_B_986141", "Treponema_B_986142",
           "Treponema_C", "Treponema_D", "Filifactor", "Streptococcus",
           "Granulicatella", "Porphyromonas_A", "Tannerella"],
    "OSCC": ["Neisseria_563205", "Haemophilus_D_735815", "Parvimonas"],
    "DC": ["Veillonella_A", "Selenomonas_A_42762"],
}
HIGHLIGHT["HC"] = sorted(set().union(*HIGHLIGHT.values()))


def load_enrichment():
    df = pd.read_csv(MODULE_DIR / "data/pan_diease_genera.csv")
    return {group: dict(zip(df.loc[df.Disease == disease, "Genus"],
                            df.loc[df.Disease == disease, "Direction"]))
            for group, disease in DISEASE_OF.items()}


def draw_network_panel(ax, G, seed=7, smax=None, title=True, enr=None):
    """Network panel WITHOUT module-strength shading and WITHOUT the density
    badge; nodes are role-marked when G.graph['role_nodes'] is set."""
    pos, ranks = layout_modules(G, seed)
    for u, v, d0 in G.edges(data=True):
        r = d0["r"]
        c = POS_C if r >= 0 else NEG_C
        ax.plot(*zip(pos[u], pos[v]), color=c, lw=0.25 + 1.6 * min(abs(r), 1),
                alpha=0.30, zorder=1, solid_capstyle="round")

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
        if enr is None:
            fill = HC_GREEN
        else:
            d = enr.get(n)
            fill = (ENR_UP if d == "+" else
                    ENR_DOWN if d == "-" else ENR_NON)
        ax.scatter(*pos[n], s=s, color=fill,
                   edgecolor=ec, linewidth=lw, zorder=3)

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


def add_title_band(fig, ax, text):
    """Grey facet-strip title band above one panel."""
    x0, y0, x1, y1 = ax.get_subplotspec().get_position(fig).extents
    bh = 0.024
    fig.add_artist(Rectangle((x0, y1 + 0.004), x1 - x0, bh,
                             facecolor="#dcdcdc", edgecolor="none",
                             transform=fig.transFigure, zorder=10))
    fig.text((x0 + x1) / 2, y1 + 0.004 + bh / 2, text, ha="center",
             va="center", fontsize=10.5, fontweight="bold", zorder=11)


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

    LINE = 0.11
    GAP, PAD = 0.10, 0.12
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

    x = PAD + w1 + GAP
    leg.text(ax(x), yh, "Role", fontsize=fs, fontweight="bold", va="center")
    role_labels = ["Module hub", "Main connector", "Hub + connector"]
    role_cols = [ROLE_RED, ROLE_BLUE, ROLE_PURPLE]
    for yy, col, lab in zip((0.60, 0.44, 0.28),
                            role_cols, role_labels):
        leg.scatter([ax(x + dia(64) / 2)], [yy], s=64,
                    facecolor="#dddddd", edgecolor=col,
                    linewidth=1.2, zorder=5)
        leg.text(ax(x + dia(64) + 0.05), yy, lab, fontsize=fs - 1.5,
                 va="center")

    x = PAD + w1 + GAP + w2 + GAP
    leg.text(ax(x), yh, "Enrichment Status", fontsize=fs,
             fontweight="bold", va="center")
    for yy, col, lab in zip((0.55, 0.36, 0.17),
                            (ENR_UP, ENR_DOWN, ENR_NON), enr_labels):
        leg.scatter([ax(x + dia(70) / 2)], [yy], s=70, facecolor=col,
                    edgecolor="none", zorder=5)
        leg.text(ax(x + dia(70) + 0.06), yy, lab, fontsize=fs - 1.5,
                 va="center")

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

def make_combined(graphs, out_base):
    """Draw four unrotated panels with a shared abundance scale and legend."""
    order = GROUPS
    nrow, ncol = 2, 2
    enrich = load_enrichment()
    gmax = max(d["abu"] for G in graphs.values() for _, d in G.nodes(data=True)
               if np.isfinite(d["abu"]))
    for group, graph in graphs.items():
        graph.graph["highlight"] = [n for n in HIGHLIGHT.get(group, []) if n in graph]
    fig = plt.figure(figsize=(10.4, 12.4))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 0.18], hspace=0.0, wspace=0.05)
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
    os.makedirs(out_base, exist_ok=True)
    for ext in ("png", "pdf"):      # pdf = vector version for publication
        base = "module_networks_2x2_font"
        out = os.path.join(out_base, f"{base}.{ext}")
        fig.savefig(out, dpi=300, bbox_inches="tight", pad_inches=0.10)
        print(f"[net] -> {out}")
    plt.close(fig)
