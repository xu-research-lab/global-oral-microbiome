#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_delta_alpha_by_bp_status.py

Re-render of results/single_pair/02_alpha/plots/delta_alpha_by_bioproject_status.png
(original: scripts/python/oral_gut_pair/report_figures.py ->
plot_alpha_by_bp_status; same pair_annot via build_alpha_status_table, so
data, BP ordering, Mann-Whitney tests and Cliff's delta are byte-identical).

Layout change requested:
  BEFORE: 3 metrics x 3 status-group columns (HC | HC+disease | disease_study)
  NOW:    3 metrics x 2 panels
    LEFT  "HC / disease": all pure-HC BioProjects merged into ONE "HC" box,
          followed by the disease_study BioProjects (one box per BP, split by
          pair_status if a BP also has Unknown pairs). Unknown-status pairs
          (status_cleaned blank, PRJNA798597 pre-eclampsia) are EXCLUDED
          from this figure.
    RIGHT "HC+disease": mixed BioProjects, HC vs disease boxes side by side,
          with the original significance highlighting (yellow band, dark-red
          box edge + bracket + BH q<0.05 stars)

Inputs :
  results/cross_comparison/01_merged_metrics/pair_metrics_all_cohorts.tsv
Outputs: results/new/figures/08_delta_alpha_by_bp_status/
  delta_alpha_by_bioproject_status.{png,pdf}
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, ttest_1samp

ROOT = Path("/home/hexincheng/projects/hxc/oral_project/oral_gut_pair")
sys.path.insert(0, str(ROOT / "scripts" / "python"))

from oral_gut_pair import report_figures as rf  # noqa: E402

IN_PAIR = ROOT / "results/cross_comparison/01_merged_metrics/pair_metrics_all_cohorts.tsv"
OUT_DIR = ROOT / "results/new/figures/08_delta_alpha_by_bp_status"

PAIR_STATUS_ORDER = rf.PAIR_STATUS_ORDER
PAIR_STATUS_COLORS = rf.PAIR_STATUS_COLORS
ALPHA_DELTA_STATUS_METRICS = rf.ALPHA_DELTA_STATUS_METRICS

# user request (2026-09-30):
#   * LEFT panel: drop the uninformative "Other" BP box (PRJDB8122, the only
#     BP whose disease_label is the catch-all "Other" status)
#   * RIGHT panel: keep "Other" under the normal HC+disease handling but
#     display its actual disease — PRJNA799190's disease side is helminth
#     infection (12 HC vs 12 helminth pairs per its pair_time_point metadata)
OTHER_DISEASE_OVERRIDE = {"Other": "Helminth infection"}
DISEASE_LABEL_DISPLAY = {
    "Parkinson's Disease": "Parkinson",
    "Helminth infection": "HI",
}

plt.rcParams.update(
    {
        "pdf.fonttype": 42,
        "font.size": 9,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
    }
)


def draw_violin_group(ax, x: float, vals: np.ndarray, color: str, sig: bool,
                      width: float = 0.42) -> None:
    """One group: KDE violin clipped to the observed range (same style as
    figs 05/06) + slim inner box; tiny/degenerate groups (n < 3 or zero
    variance) fall back to a narrow rectangle.

    user request (2026-09-30): violin plot + ALL jittered points shown
    (previously boxplot with at most 80 sampled points).
    """
    sig = bool(sig)
    edge = "#8b0000" if sig else "#333"
    lw = 1.6 if sig else 0.8
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    ok = len(vals) >= 3 and float(np.ptp(vals)) > 0
    if ok:
        try:
            kde = gaussian_kde(vals)
            y = np.linspace(vals.min(), vals.max(), 150)
            d = kde(y)
            d = d / d.max() * (width / 2.0)
            ax.fill_betweenx(y, x - d, x + d, facecolor=color, alpha=0.32,
                             edgecolor=edge, linewidth=lw, zorder=2)
        except np.linalg.LinAlgError:
            ok = False
    if not ok:
        if len(vals) >= 2:
            ys = [vals.min(), vals.max()]
        else:
            ys = [vals[0] - 0.5, vals[0] + 0.5]
        ax.fill_betweenx(ys, x - width / 4.0, x + width / 4.0,
                         facecolor=color, alpha=0.32, edgecolor=edge,
                         linewidth=lw, zorder=2)
    if len(vals):
        ax.boxplot(
            [vals], positions=[x], widths=0.12, patch_artist=True,
            showfliers=False,
            boxprops={"facecolor": "white", "edgecolor": "black",
                      "linewidth": 0.8},
            medianprops={"color": "black", "linewidth": 1.4},
            whiskerprops={"linewidth": 0.8, "color": "black"},
            capprops={"linewidth": 0.8, "color": "black"},
            zorder=3,
        )


def draw_group_boxes(
    ax,
    sig_bps: set,
    positions: list,
    data: list,
    colors: list,
    box_bps: list,
    rng,
) -> None:
    """Shared violin + jitter drawing for one panel (user request)."""
    for pos, vals, color, bp in zip(positions, data, colors, box_bps):
        if len(vals) == 0:
            continue
        draw_violin_group(ax, pos, np.asarray(vals, float), color,
                          bp in sig_bps)
    # every pair as a small jittered dot over the violin (alpha raised above
    # the violin fill so points read through — user feedback)
    for pos, vals, color in zip(positions, data, colors):
        if len(vals) == 0:
            continue
        vals = np.asarray(vals, dtype=float)
        vals = vals[np.isfinite(vals)]
        ax.scatter(
            np.full(len(vals), pos) + rng.uniform(-0.09, 0.09, len(vals)),
            vals,
            s=9,
            alpha=0.75,
            color=color,
            edgecolors="white",
            linewidths=0.3,
            zorder=3,
        )


def draw_panel(ax, sub: pd.DataFrame, col: str, rng, merged_hc: bool,
               sig_map: dict, metric_name: str, ttest_rows: list = None) -> None:
    """One panel. merged_hc=True -> LEFT (merged HC + disease BPs);
    merged_hc=False -> RIGHT (HC+disease mixed BPs, sig highlighting).
    Appends one-sample t-test rows (vs zero) to ttest_rows for the
    cross-panel BH correction done in main()."""
    if merged_hc:
        hc_sub = sub.loc[sub["status_group"].eq("HC")]
        hc_vals = (
            hc_sub.loc[hc_sub["pair_status"].eq("HC"), col]
            .astype(float)
            .dropna()
            .values
        )
        # Unknown-status pairs (blank status_cleaned) are excluded entirely
        order = (
            sub.loc[sub["status_group"].eq("disease_study")]
            .groupby("bioproject")["delta_alpha_observed"]
            .median()
            .sort_values()
            .index.tolist()
        )
        # user request: the left panel's catch-all "Other" BP box is removed
        order = [
            bp for bp in order
            if str(sub.loc[sub["bioproject"].eq(bp), "disease_label"].iloc[0] or "")
            != "Other"
        ]
        groups = []  # list of (tick_label, [(status, x, vals), ...], sig?)
        if len(hc_vals):
            groups.append(("HC\n(all HC BPs)", [("HC", 1.0, hc_vals)], False))
        x = 2.1
        for bp in order:
            bp_sub = sub.loc[sub["bioproject"].eq(bp)]
            statuses = [s for s in PAIR_STATUS_ORDER
                        if (bp_sub["pair_status"] == s).any()]
            if not statuses:
                continue
            disease_lab = str(bp_sub["disease_label"].iloc[0] or "")
            disease_lab = DISEASE_LABEL_DISPLAY.get(disease_lab, disease_lab)
            items = []
            for st in statuses:
                vals = (
                    bp_sub.loc[bp_sub["pair_status"].eq(st), col]
                    .astype(float).dropna().values
                )
                items.append((st, x, vals))
                x += 0.62
            groups.append((f"{bp}\n({disease_lab})", items, False))
            x += 0.12
    else:
        order = (
            sub.loc[sub["status_group"].eq("HC+disease")]
            .groupby("bioproject")["delta_alpha_observed"]
            .median()
            .sort_values()
            .index.tolist()
        )
        groups = []
        x = 1.0
        for bp in order:
            bp_sub = sub.loc[sub["bioproject"].eq(bp)]
            statuses = [s for s in PAIR_STATUS_ORDER
                        if (bp_sub["pair_status"] == s).any()]
            if not statuses:
                continue
            disease_lab = str(bp_sub["disease_label"].iloc[0] or "")
            disease_lab = DISEASE_LABEL_DISPLAY.get(disease_lab, disease_lab)
            # right panel: show the actual disease behind the catch-all "Other"
            disease_lab = OTHER_DISEASE_OVERRIDE.get(disease_lab, disease_lab)
            disease_lab = DISEASE_LABEL_DISPLAY.get(disease_lab, disease_lab)
            items = []
            for st in statuses:
                vals = (
                    bp_sub.loc[bp_sub["pair_status"].eq(st), col]
                    .astype(float).dropna().values
                )
                items.append((st, x, vals))
                x += 0.62
            groups.append((f"{bp}\n({disease_lab})", items, False))
            x += 0.12

    if not groups:
        ax.set_visible(False)
        return

    sig_bps = set()
    if not merged_hc:
        for bp in order:
            info = sig_map.get((bp, metric_name))
            stars = ""
            if info is not None:
                stars = info.get("sig_stars", "") or ""
                if isinstance(stars, float) and np.isnan(stars):
                    stars = ""
            if stars:
                sig_bps.add(bp)

    positions, data, colors, ns, box_bps = [], [], [], [], []
    g_names = []
    bp_pos_map = {}
    tick_pos, tick_labels, tick_bps = [], [], []
    for label, items, _ in groups:
        bp = label.split("\n")[0]
        group_xs = []
        for st, x, vals in items:
            positions.append(x)
            data.append(vals)
            colors.append(PAIR_STATUS_COLORS[st])
            g_names.append(st)
            ns.append(len(vals))
            box_bps.append(bp)
            group_xs.append(x)
        bp_pos_map[bp] = group_xs
        tick_pos.append(float(np.mean(group_xs)))
        tick_bps.append(bp)
        tick_labels.append(label)
    if tick_bps and tick_bps[0] == "HC":
        tick_bps[0] = "HC_merged"

    # annotation tiers on the ROW-level scale (the left/right panels of a
    # row share the y axis) so the n labels, vs-zero stars and MWU brackets
    # sit at identical heights in both panels
    row_vals = sub[col].astype(float).dropna().values
    if len(row_vals):
        y_hi = float(np.nanpercentile(row_vals, 99))
        y_lo = float(np.nanpercentile(row_vals, 1))
    else:
        y_hi, y_lo = 1.0, 0.0
    pad = (y_hi - y_lo) * 0.06 if y_hi > y_lo else 1.0

    # yellow background band + brackets/stars for significant mixed BPs
    if sig_bps:
        for bp in sig_bps:
            xs = bp_pos_map.get(bp, [])
            if not xs:
                continue
            ax.axvspan(min(xs) - 0.27, max(xs) + 0.27,
                       facecolor="#fff3cd", edgecolor="none",
                       alpha=0.85, zorder=0)
            if len(xs) >= 2:
                info = sig_map.get((bp, metric_name))
                stars = info.get("sig_stars", "") if info is not None else ""
                # bracket raised above the vs-zero star tier (which sits at
                # y_hi + 0.85*pad, ~0.5pad tall) so the two star sets and
                # the bracket line never collide (overlap-checked by extent)
                y_br = y_hi + pad * 2.6
                ax.plot([xs[0], xs[0], xs[-1], xs[-1]],
                        [y_br - pad * 0.12, y_br, y_br, y_br - pad * 0.12],
                        color="#8b0000", linewidth=1.3, zorder=5)
                ax.text((xs[0] + xs[-1]) / 2, y_br + pad * 0.05, stars,
                        ha="center", va="bottom", fontsize=9,
                        fontweight="bold", color="#8b0000", zorder=5)

    draw_group_boxes(ax, sig_bps, positions, data, colors, box_bps, rng)

    # n labels (tiers already on the row-level scale, see above)
    for pos, n in zip(positions, ns):
        ax.text(pos, y_hi + pad * 0.15, f"n={n}", ha="center", va="bottom",
                fontsize=5, color="#444")

    # user request: each group first tested against ZERO (two-sided
    # one-sample t-test on the paired Δ values); p values are collected and
    # BH-corrected in main() across all groups x metrics, stars drawn after
    if ttest_rows is not None:
        for pos, vals, gname, bp in zip(positions, data, g_names, box_bps):
            v = np.asarray(vals, dtype=float)
            v = v[np.isfinite(v)]
            if len(v) < 3:
                continue
            t, p = ttest_1samp(v, 0.0)
            ttest_rows.append({
                "metric": metric_name,
                "panel": "left" if merged_hc else "right",
                "bioproject": bp,
                "group": gname,
                "n": int(len(v)),
                "mean_delta": float(v.mean()),
                "t_stat": float(t),
                "df": int(len(v) - 1),
                "p_value": float(p),
                "ax": ax,
                "x": float(pos),
                "y_star": float(y_hi + pad * 0.85),
            })

    ax.axhline(0, color="#d62728", linestyle="--", linewidth=1, zorder=1)
    ax.set_xticks(tick_pos)
    # axis font sizes unified with fig 05 (x ticks 8 pt) — kept vertical
    # (user request: 标签还是竖着), so per-box width limits usable size
    tick_fs = 8.0
    tick_objs = ax.set_xticklabels(tick_labels, rotation=45, ha="right",
                                   va="top", rotation_mode="anchor",
                                   fontsize=tick_fs)
    for lab, bp in zip(tick_objs, tick_bps):
        if bp in sig_bps:
            lab.set_fontweight("bold")
            lab.set_color("#8b0000")
            lab.set_fontsize(tick_fs + 0.6)
    ax.set_xlim(0.5, positions[-1] + 0.45)
    ax.tick_params(axis="x", pad=2)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    pair_df = pd.read_csv(IN_PAIR, sep="\t")
    pair_annot, tests_df, _pooled = rf.build_alpha_status_table(pair_df)
    # drop Unknown-status pairs (blank status_cleaned, PRJNA798597) entirely
    n_unk = int(pair_annot["pair_status"].eq("Unknown").sum())
    pair_annot = pair_annot.loc[~pair_annot["pair_status"].eq("Unknown")].copy()
    print(f"dropped {n_unk} Unknown-status pairs")
    sig_map = {}
    if tests_df is not None and not tests_df.empty:
        for _, row in tests_df.iterrows():
            sig_map[(row["bioproject"], row["metric"])] = row

    n_boxes_left = 1 + pair_annot.loc[
        pair_annot["status_group"].eq("disease_study")
        & pair_annot["disease_label"].ne("Other"),   # Other BP box dropped (user request)
        "bioproject",
    ].nunique()
    n_boxes_right = 0
    for bp, bp_sub in pair_annot.loc[
        pair_annot["status_group"].eq("HC+disease")
    ].groupby("bioproject"):
        n_boxes_right += max(1, len(
            [s for s in PAIR_STATUS_ORDER if (bp_sub["pair_status"] == s).any()]
        ))
    width_ratios = [n_boxes_left, max(1, n_boxes_right)]

    # canvas shrunk (user request): ~0.75x of the previous 16.75x13.2 in
    # width, rows trimmed to 3.0 in each; fonts unchanged so labels read
    # slightly larger in the crop
    fig = plt.figure(figsize=(max(11.0, (sum(width_ratios)) * 0.4 + 2.0), 9.3))
    gs = fig.add_gridspec(3, 2, width_ratios=width_ratios,
                          wspace=0.08, hspace=0.55)

    legend_handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=PAIR_STATUS_COLORS["HC"],
                      edgecolor="#333", label="HC"),
        plt.Rectangle((0, 0), 1, 1, facecolor=PAIR_STATUS_COLORS["disease"],
                      edgecolor="#333", label="disease"),
        plt.Rectangle((0, 0), 1, 1, facecolor="#fff3cd", edgecolor="#8b0000",
                      linewidth=1.5, label="HC vs disease significant (BH q<0.05)"),
    ]

    panel_titles = ["HC (merged) / disease", "HC+disease"]

    # one-sample t-test rows accumulated across ALL panels and metrics
    ttest_rows = []

    for row_i, (metric_name, col) in enumerate(ALPHA_DELTA_STATUS_METRICS):
        axes_row = []
        for col_i, merged_hc in enumerate((True, False)):
            sharey = axes_row[0] if axes_row else None
            ax = fig.add_subplot(gs[row_i, col_i], sharey=sharey)
            axes_row.append(ax)
            rng = np.random.default_rng(2026 + row_i * 17 + col_i)
            draw_panel(ax, pair_annot, col, rng, merged_hc, sig_map,
                       metric_name, ttest_rows=ttest_rows)
            if col_i == 0:
                ax.set_ylabel(f"Δ {metric_name}\n(Gut − Oral)", fontsize=9)
            else:
                plt.setp(ax.get_yticklabels(), visible=False)
            if row_i == 0:
                ax.set_title(panel_titles[col_i], fontsize=10, pad=6)

        row_vals = pair_annot[col].astype(float).dropna().values
        if len(row_vals):
            lo = float(np.nanpercentile(row_vals, 1))
            hi = float(np.nanpercentile(row_vals, 99))
            span = hi - lo if hi > lo else 1.0
            for ax in axes_row:
                if ax.get_visible():
                    ax.set_ylim(lo - 0.08 * span, hi + 0.35 * span)

    fig.legend(handles=legend_handles, loc="upper center", ncol=3,
               frameon=False, bbox_to_anchor=(0.5, 0.998), fontsize=8)

    # ---- finalize the vs-zero t-tests: BH across all groups x metrics
    # (single family), stars above each tested group, TSV written ----
    if ttest_rows:
        qs = rf._bh_fdr([r["p_value"] for r in ttest_rows])
        for r, q in zip(ttest_rows, qs):
            r["q_BH"] = float(q)
            stars = rf._sig_stars(q)
            r["sig_stars"] = stars
            if stars:
                r["ax"].text(r["x"], r["y_star"], stars, ha="center",
                             va="bottom", fontsize=7, fontweight="bold",
                             color="#1a1a1a", zorder=6)
        out_cols = ["metric", "panel", "bioproject", "group", "n",
                    "mean_delta", "t_stat", "df", "p_value", "q_BH",
                    "sig_stars"]
        pd.DataFrame([{k: r[k] for k in out_cols} for r in ttest_rows]).to_csv(
            OUT_DIR / "delta_alpha_one_sample_t_tests.tsv", sep="\t", index=False)
        n_sig = int(sum(bool(r["sig_stars"]) for r in ttest_rows))
        print(f"one-sample t-tests vs 0: {len(ttest_rows)} groups tested, "
              f"{n_sig} with BH q<0.05 -> delta_alpha_one_sample_t_tests.tsv")

    fig.suptitle(
        "Single-pair Δ alpha by BioProject — HC merged (left), "
        "HC+disease mixed BPs (right)",
        y=1.025, fontsize=10,
    )
    fig.subplots_adjust(bottom=0.12, top=0.92)
    for ext in ("png", "pdf"):
        fig.savefig(OUT_DIR / f"delta_alpha_by_bioproject_status.{ext}",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"outputs -> {OUT_DIR}")
    print(f"left boxes={n_boxes_left}, right boxes={n_boxes_right}, "
          f"mixed BPs={pair_annot.loc[pair_annot['status_group'].eq('HC+disease'), 'bioproject'].nunique()}, "
          f"HC-merged n={int(((pair_annot['status_group'] == 'HC') & (pair_annot['pair_status'] == 'HC')).sum())}")


if __name__ == "__main__":
    main()