"""Module-circle geometry and label placement shared by network panels."""
import numpy as np
import matplotlib.text as mtext

FS = 1.15
OUTSIDE_OFFSET_PT = 3.0

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

    n_mod = len(ranks)
    order = sorted(ranks, key=lambda m: (-sizes[m], m))
    mpos = {}
    for k, m in enumerate(order):
        ang = np.pi / 2 - 2 * np.pi * k / n_mod
        mpos[m] = np.array([np.cos(ang), np.sin(ang)]) * 4.5

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

    pos = {}
    for m in ranks:
        members = [n for n, d0 in G.nodes(data=True) if d0["rank"] == m]
        abu = {n: (G.nodes[n]["abu"] if np.isfinite(G.nodes[n]["abu"]) else 0.0)
               for n in members}
        members.sort(key=lambda n: (-abu[n], n))
        k = len(members)
        for i, n in enumerate(members):
            ang = np.pi / 2 - 2 * np.pi * i / k   # clockwise
            pos[n] = (mpos[m][0] + rad[m] * np.cos(ang),
                      mpos[m][1] + rad[m] * np.sin(ang))
    return pos, ranks

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
    hl_nodes = set(G.graph.get("highlight") or ())
    for n in hl_nodes - labelled:
        items.append((n, "#000000"))
        labelled.add(n)
    mod_centre = {}
    for m in {G.nodes[n]["rank"] for n, _ in items}:
        members = [n for n, d0 in G.nodes(data=True) if d0["rank"] == m]
        mod_centre[m] = (float(np.mean([pos[n][0] for n in members])),
                         float(np.mean([pos[n][1] for n in members])))
    items.sort(key=lambda it: (
        G.nodes[it[0]]["rank"],
        (lambda c, p: np.arctan2(p[0] - c[0], -(p[1] - c[1])))(
            mod_centre[G.nodes[it[0]]["rank"]], pos[it[0]])))
    labels = []
    for hn, c in items:
        m = G.nodes[hn]["rank"]
        cx, cy = mod_centre[m]
        disp = ax.transData.transform([(cx, cy), pos[hn]])
        sdx, sdy = disp[1] - disp[0]
        snrm = np.hypot(sdx, sdy) or 1.0
        ux, uy = sdx / snrm, sdy / snrm
        ang = np.degrees(np.arctan2(uy, ux))
        if ux < 0:
            rot, ha = ang + 180, "right"
        else:
            rot, ha = ang, "left"
        big = hn in hl_nodes   # highlighted -> bold + enlarged
        ann = ax.annotate(
            hn.replace("g__", ""), pos[hn],
            xytext=(ux * OUTSIDE_OFFSET_PT, uy * OUTSIDE_OFFSET_PT),
            textcoords="offset points",
            ha=ha, va="center", rotation=rot,
            rotation_mode="anchor",
            fontsize=(5.2 if big else 4) * FS,
            fontweight="bold" if big else "normal",
            color=c, zorder=7, path_effects=stroke)
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
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

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

        if G.graph.get("site") == "PD":
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

        base_rot = {}
        for a in anns:
            base_rot[id(a)] = a.get_rotation()

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

        OFF_MULTS = (1.0, 1.6, 2.3, 3.1, 4.0, 5.2, 7.0, 9.5)

        def set_label_rot(a, delta_deg):
            a.set_rotation(base_rot[id(a)] + delta_deg)

        def set_label_off(a, mult):
            d = getattr(a, "label_dir", None)
            if d is not None:
                ux, uy = d
                a.set_position((ux * OUTSIDE_OFFSET_PT * mult,
                                uy * OUTSIDE_OFFSET_PT * mult))

        def rotate_relax(max_deg=90.0, step=2.5):
            """Greedy sequential placement: labels are processed in order; each
            earlier label is FROZEN once placed, and each later label is rotated
            (trying both directions) AND/OR pushed farther out radially until
            it no longer overlaps any already-placed label, then frozen itself.
            Avoids the many-neighbour cancellation that stranded mid-sequence
            labels. A single renderer snapshot (drawn once) is reused for all
            probes — moving one label doesn't move any other element, so no
            re-render needed."""
            placed_boxes = []
            n_step = int(round(max_deg / step))
            cands = ([0.0]
                     + [step * i for i in range(1, n_step + 1)]
                     + [-step * i for i in range(1, n_step + 1)])
            combos_of = {id(a): [(d, m) for d in upright_cands(a, cands)
                                 for m in OFF_MULTS] for a in anns}
            for a in anns:
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
                placed_boxes.append(tbox(a))

        def badge_relax(rounds, cap=60.0):
            bitems = [("badge", t) for t in static]
            if not bitems:
                return
            pushers = list(anns)
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

        rotate_relax()

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
                    partner = None
                    for j, b in enumerate(anns):
                        if j == i:
                            continue
                        if ba.overlaps(boxes[id(b)]):
                            partner = b
                            break
                    if partner is None:
                        continue
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

