#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_alpha_oral_gut_spearman_by_disease.py

Re-render of results/single_pair/02_alpha/plots/alpha_oral_gut_spearman.png
(original: scripts/python/oral_gut_pair/report_figures.py ->
plot_alpha_oral_gut_spearman, same data source: delta_alpha.tsv, dataset
filter_ra, pair status from pair_metrics_all_cohorts.tsv), with changes:

  * points colored by HC and by DISEASE CLASS unified with the fig 08
    delta-alpha violin classification (per-BioProject disease label from
    report_figures._bp_disease_label): every non-HC pair takes its
    BioProject's disease label, so statuses that alternate inside one
    study are merged the same way as in the violin figure, e.g.
    PRJEB11425 ICU + "ICU discharge or on ICU day 10" -> "ICU +1" and
    PRJNA756623 -> "COVID conval. +1"; PRJNA799190's catch-all "Other"
    disease side reads "Helminth infection" like fig 08's right panel;
    HC pairs all keep the blue "HC" class (fig 08's merged-HC box)
  * axes on log10(x + 10⁻⁶) scale for the "..._by_disease" version (raw values
    pile up in the lower-left corner) and raw scale for the "..._raw"
    version; Spearman is rank-based, so rho is identical on both scales
  * TWO smoothed trend lines per panel (no straight lines):
      - HC (blue):  lowess fit over all HC pairs
      - all disease (black): lowess fit over all non-HC pairs
    each with a bootstrap 95% CI band (100 resamples, frac=0.3)

Inputs :
  results/single_pair/02_alpha/delta_alpha.tsv                      (filter_ra)
  results/cross_comparison/01_merged_metrics/pair_metrics_all_cohorts.tsv
Outputs: results/new/figures/07_alpha_spearman_by_disease/
  alpha_oral_gut_spearman_by_disease.{png,pdf}       (log10(x+10⁻⁶) scale)
  alpha_oral_gut_spearman_by_disease_raw.{png,pdf}   (raw scale)
  alpha_oral_gut_spearman_by_disease_raw_{observed,shannon,invsimpson}.{png,pdf}
      (the same raw-scale plot saved as one file per alpha-diversity metric)
  alpha_oral_gut_spearman_by_disease_log10_dropzero.{png,pdf}
      (third version, user request 2026-09-30: zero-valued samples removed
       with their paired partner first, then bare log10 — no pseudo-count)
  alpha_oral_gut_spearman_by_disease_log10_dropzero_density.{png,pdf}
      (density version on zero-filtered raw metadata values; no log transform,
       with the existing KDE and LOWESS styling)
  alpha_oral_gut_spearman_by_status.tsv   (per-metric per-status rho/p/n +
     q_value_BH_within_metric — BH-FDR within each metric across the nine
     subgroup statuses, All/HC pooled rows not in the family)
  alpha_oral_gut_spearman_by_status_log10_dropzero.tsv
     (same table layout for the zero-filtered version — family = the eight
     plotted disease classes, All/HC rows not in the family)
"""

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import gaussian_kde, spearmanr
from statsmodels.nonparametric.smoothers_lowess import lowess

ROOT = Path("/home/hexincheng/projects/hxc/oral_project/oral_gut_pair")
# shared disease-abbreviation table (user-supplied mapping, same source as
# figs 03/08/09/10)
sys.path.insert(0, str(ROOT / "scripts" / "python"))
from oral_gut_pair.report_figures import DISEASE_LABEL_SHORT as DISEASE_ABBREV  # noqa: E402
from oral_gut_pair import report_figures as rf  # noqa: E402
IN_DELTA = ROOT / "results/single_pair/02_alpha/delta_alpha.tsv"
IN_META = ROOT / "results/cross_comparison/01_merged_metrics/pair_metrics_all_cohorts.tsv"
OUT_DIR = ROOT / "results/new/figures/07_alpha_spearman_by_disease"

MIN_N = 25          # merged disease classes with fewer pairs lump into "Other"

ALPHA_METRICS = [
    ("Observed", "oral_Observed", "gut_Observed"),
    ("Shannon", "oral_Shannon", "gut_Shannon"),
    ("InvSimpson", "oral_InvSimpson", "gut_InvSimpson"),
]

# unified plot-area geometry (user request 2026-10-20): identical to figs
# 03/05/06/11 — canvas 3.40 in tall, axes 2.26 in, bottom edge 0.70 in /
# top edge 2.96 in from the canvas bottom, saved WITHOUT
# bbox_inches="tight"; the panel is an exact 2.26x2.26-in SQUARE because
# draw_panel uses aspect="equal" (matplotlib would otherwise shrink the
# spine box away from the assigned rect and break the alignment)
CANVAS_H = 3.40
AX_BOTTOM_IN = 0.70
AX_TOP_IN = 2.96
PANEL_IN = 2.26          # AX_TOP_IN - AX_BOTTOM_IN, square side
FIG_W = 4.60             # panel + right-hand legend zone
AX_RECT07 = [0.62 / FIG_W, AX_BOTTOM_IN / CANVAS_H,
             PANEL_IN / FIG_W, PANEL_IN / CANVAS_H]

HC_COLOR = "#2166AC"
OTHER_COLOR = "#8C8C8C"
HC_LINE_COLOR = HC_COLOR
DISEASE_LINE_COLOR = "#1A1A1A"
# qualitative palette for the individual diseases (tab20 + extras)
_DISEASE_COLORS = [
    "#D65F4B", "#E8A33D", "#63B6C7", "#7A5195", "#76B7B2",
    "#E15759", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7",
    "#9C755F", "#BAB0AC", "#1B9E77", "#D95F02", "#66A61E",
]

plt.rcParams.update(
    {
        "pdf.fonttype": 42,
        "font.size": 9,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
    }
)


def load_pairs() -> pd.DataFrame:
    delta = pd.read_csv(IN_DELTA, sep="\t")
    delta = delta.loc[delta["dataset_label"].eq("filter_ra")].copy()

    meta = pd.read_csv(IN_META, sep="\t")
    sp_meta = meta.loc[
        meta["analysis_set"].astype(str).eq("single_pair"),
        ["pair_id", "status_cleaned"],
    ].drop_duplicates("pair_id")
    df = delta.merge(sp_meta, on="pair_id", how="left")
    df["status_cleaned"] = df["status_cleaned"].fillna("Unknown").astype(str)
    df.loc[df["status_cleaned"].isin(["", "nan", "NA"]), "status_cleaned"] = "Unknown"

    # disease classes unified with fig 08 (user request 2026-09-30):
    # non-HC pairs are grouped by the BioProject's _bp_disease_label — the
    # exact label set used by the 08 delta-alpha violin figure — so
    # alternate time-point statuses inside one study are merged here as
    # well (e.g. PRJEB11425 -> "ICU +1", PRJNA756623 -> "COVID conval. +1")
    df["pair_status"] = df["status_cleaned"].map(rf._pair_status_binary)
    bp_labels = {}
    for bp, bp_sub in df.groupby("bioproject"):
        lab = rf._bp_disease_label(bp_sub) or "Other"
        # fig 08 right-panel override: PRJNA799190's disease side is
        # helminth infection, not the catch-all "Other"
        if bp == "PRJNA799190":
            lab = "Helminth infection"
        bp_labels[bp] = lab
    df["status_plot"] = np.where(df["pair_status"].eq("HC"), "HC",
                                 df["bioproject"].map(bp_labels))

    # user request (2026-09-30): AFTER merging into BioProject-level classes,
    # the old visibility rule still applies — only classes with >= MIN_N
    # pairs are drawn individually; the rest (small diseases, Unknown and
    # catch-all status pairs) go back into the grey "Other" bag
    counts = df["status_plot"].value_counts()
    keep = [s for s in counts.index
            if s not in ("HC", "Other") and counts[s] >= MIN_N]
    df["status_plot"] = df["status_plot"].where(
        df["status_plot"].isin(keep) | df["status_plot"].eq("HC"), "Other")
    return df


def status_order_and_colors(df: pd.DataFrame):
    counts = df["status_plot"].value_counts()
    diseases = [s for s in counts.index if s not in ("HC", "Other")]
    order = ["HC"] + sorted(diseases, key=lambda s: (-counts[s], s)) + ["Other"]
    colors = {}
    di = 0
    for s in order:
        if s == "HC":
            colors[s] = HC_COLOR
        elif s == "Other":
            colors[s] = OTHER_COLOR
        else:
            colors[s] = _DISEASE_COLORS[di % len(_DISEASE_COLORS)]
            di += 1
    return order, colors


def log10p(v: np.ndarray) -> np.ndarray:
    # user request: pseudo-count 1e-6 instead of +1 (abundance-scale offset)
    return np.log10(np.asarray(v, dtype=float) + 1e-6)


def log10_plain(v: np.ndarray) -> np.ndarray:
    # third version (user request 2026-09-30): no pseudo-count at all —
    # requires zero-valued samples/pairs to be removed first (see drop_zero)
    return np.log10(np.asarray(v, dtype=float))


def drop_zero_pairs(df: pd.DataFrame) -> pd.DataFrame:
    """User request (2026-09-30): before taking the bare log10, remove every
    PAIR in which any plotted metric has a zero value on either side (a
    zero-valued sample cannot be dropped alone — the Spearman unit is the
    oral-gut pair, so its matched partner is removed with it).  All six
    metric columns (Observed/Shannon/InvSimpson x oral/gut) must be finite
    and > 0; the same filtered set feeds all three panels so the panels stay
    comparable."""
    mask = pd.Series(True, index=df.index)
    for _, oral_col, gut_col in ALPHA_METRICS:
        for c in (oral_col, gut_col):
            v = pd.to_numeric(df[c], errors="coerce")
            mask &= v.notna() & v.gt(0)
    out = df.loc[mask].copy()
    n_drop = len(df) - len(out)
    print(f"zero-pair filter: dropped {n_drop} of {len(df)} pairs "
          f"-> {len(out)} pairs kept ({len(out)/len(df)*100:.1f}%)")
    return out


def fmt_p_display(p: float) -> str:
    """Figure-annotation p display. Underflowed p (p == 0.0 — the true
    value falls below the ~1e-308 double-precision floor; e.g. the
    All-pairs Observed p is ~1e-385 by the t-approximation) would render
    as the nonsensical '0.00e+00', so show those as an inequality bound.
    User request 2026-10-25: p=0.00e+00 is not an acceptable display."""
    if not np.isfinite(p) or p <= 0.0:
        return "<1e-300"
    return f"{p:.2e}"


def spearman_safe(x: np.ndarray, y: np.ndarray):
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = int(len(x))
    if n < 3:
        return float("nan"), float("nan"), n
    rho, p = spearmanr(x, y)
    return float(rho), float(p), n


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.clip(ranked, 0, 1)
    return q


def _lowess_curve(x: np.ndarray, y: np.ndarray, grid: np.ndarray,
                  frac: float) -> np.ndarray:
    res = lowess(y, x, frac=frac, it=0, return_sorted=True)
    return np.interp(grid, res[:, 0], res[:, 1])


def class_outer_contour(ax, xs: np.ndarray, ys: np.ndarray, color: str,
                        lim, n_grid: int = 200, frac_of_peak: float = 0.10,
                        lw: float = 1.3) -> None:
    """ONE outermost KDE contour per class (line at 10% of the class's own
    peak density).  A single coloured outline per class keeps every class
    identifiable on the density background; multiple levels per class just
    tangl up when all class clouds overlap the same region."""
    kde = gaussian_kde(np.vstack([xs, ys]))
    xg = np.linspace(lim[0], lim[1], n_grid)
    yg = np.linspace(lim[0], lim[1], n_grid)
    xx, yy = np.meshgrid(xg, yg)
    zz = kde(np.vstack([xx.ravel(), yy.ravel()])).reshape(len(yg), len(xg))
    ax.contour(xx, yy, zz, levels=[frac_of_peak * zz.max()],
               colors=[color], linewidths=lw, zorder=2.6)


def lowess_band(ax, x: np.ndarray, y: np.ndarray, color: str, rng,
                n_boot: int = 100, frac: float = 0.3,
                z_fill: float = 1.5, z_line: float = 1.6) -> None:
    """LOWESS smooth + bootstrap 95% CI band (the geom_smooth look)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 25:
        return
    grid = np.linspace(x.min(), x.max(), 200)
    fits = [_lowess_curve(x, y, grid, frac)]
    idx = rng.integers(0, len(x), size=(n_boot, len(x)))
    for b in range(n_boot):
        fits.append(_lowess_curve(x[idx[b]], y[idx[b]], grid, frac))
    fits = np.asarray(fits)
    lo, hi = np.percentile(fits, [2.5, 97.5], axis=0)
    ax.fill_between(grid, lo, hi, color=color, alpha=0.15, linewidth=0,
                    zorder=z_fill)
    ax.plot(grid, fits[0], color=color, linewidth=1.8, zorder=z_line)


def draw_panel(ax, df: pd.DataFrame, name: str, oral_col: str, gut_col: str,
               order, colors, tf=None, tf_label: str = "", rng=None,
               density: bool = False, label_metric: bool = True) -> dict:
    """tf=None -> raw scale; tf=log10p -> log10(x+10⁻⁶) scale.
    density=True (user request 2026-09-30): replace the scatter dots with
    one filled 2-D KDE layer per class (same class colours, common_norm
    off so every class keeps its own scale), everything else unchanged."""
    sub = df[[oral_col, gut_col, "status_plot"]].copy()
    for c in (oral_col, gut_col):
        sub[c] = pd.to_numeric(sub[c], errors="coerce")
    sub = sub.dropna(subset=[oral_col, gut_col])

    oral_raw = sub[oral_col].to_numpy(float)
    gut_raw = sub[gut_col].to_numpy(float)
    oral = oral_raw if tf is None else tf(oral_raw)
    gut = gut_raw if tf is None else tf(gut_raw)

    lo = float(min(oral.min(), gut.min()))
    hi = float(max(oral.max(), gut.max()))
    pad = (hi - lo) * 0.04
    lim = (lo - pad, hi + pad)

    ax.plot(lim, lim, linestyle="--", color="#bdbdbd", linewidth=1.0, zorder=1)

    is_hc = sub["status_plot"].eq("HC").to_numpy()
    if density:
        # density render with FULL class information (user request
        # 2026-09-30): the big HC class is a blue FILLED KDE (nested levels,
        # like the shaded mass everywhere the HC pairs pile up), while every
        # disease class contributes exactly ONE thick coloured outline (the
        # outermost contour of its own KDE at 5% of its own peak) — stacking
        # 9 translucent fills blurs into mud and 9 multi-level contour
        # families tangle into spaghetti, but one outer ring per class keeps
        # the disease classification readable on the density background
        ax_lim = lim
        for st in order:
            st_sub = sub.loc[sub["status_plot"].eq(st)]
            if len(st_sub) < 25:
                continue
            xs_ = st_sub[oral_col].to_numpy(float)
            ys_ = st_sub[gut_col].to_numpy(float)
            if tf is not None:
                xs_ = tf(xs_)
                ys_ = tf(ys_)
            if st == "HC":
                sns.kdeplot(x=xs_, y=ys_, ax=ax, fill=True, levels=8,
                            thresh=0.04, color=colors[st], alpha=0.25,
                            warn_singular=False, zorder=2.4)
            else:
                class_outer_contour(ax, xs_, ys_, colors[st], ax_lim)
        lowess_band(ax, oral[is_hc], gut[is_hc], HC_LINE_COLOR, rng,
                    z_fill=4.4, z_line=4.5)
        lowess_band(ax, oral[~is_hc], gut[~is_hc], DISEASE_LINE_COLOR, rng,
                    z_fill=4.4, z_line=4.5)
    else:
        # two smoothed trend lines: all HC vs all disease (non-HC)
        lowess_band(ax, oral[is_hc], gut[is_hc], HC_LINE_COLOR, rng)
        lowess_band(ax, oral[~is_hc], gut[~is_hc], DISEASE_LINE_COLOR, rng)

        for st in order:
            st_sub = sub.loc[sub["status_plot"].eq(st)]
            if st_sub.empty:
                continue
            xs = st_sub[oral_col].to_numpy(float)
            ys = st_sub[gut_col].to_numpy(float)
            if tf is not None:
                xs, ys = tf(xs), tf(ys)
            ax.scatter(xs, ys, s=9, color=colors[st], alpha=0.55,
                       linewidths=0, zorder=3)

    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal", adjustable="box")
    if label_metric:
        ax.set_xlabel(f"Oral {name}{tf_label}", fontsize=8)
        ax.set_ylabel(f"Gut {name}{tf_label}", fontsize=8)
        ax.set_title(name, fontsize=10)
    else:
        # per-metric raw figures (user request 2026-10-14): the panel
        # subtitle ("Observed" etc.) is DELETED and the metric name moves
        # into the axis labels instead — "Oral observed OTUs" (richness is
        # a count of observed OTUs), "Oral Shannon", "Oral InvSimpson"
        desc = "observed OTUs" if name == "Observed" else name
        ax.set_xlabel(f"Oral {desc}{tf_label}", fontsize=8)
        ax.set_ylabel(f"Gut {desc}{tf_label}", fontsize=8)

    rho_all, p_all, n_all = spearman_safe(oral, gut)
    rho_hc, p_hc, n_hc = spearman_safe(oral[is_hc], gut[is_hc])
    rho_di, p_di, n_di = spearman_safe(oral[~is_hc], gut[~is_hc])
    ax.text(
        0.04, 0.98,
        (f"All: ρ={rho_all:.3f}, p={fmt_p_display(p_all)}, n={n_all}\n"
         f"HC: ρ={rho_hc:.3f}, p={fmt_p_display(p_hc)}, n={n_hc}\n"
         f"Disease: ρ={rho_di:.3f}, p={fmt_p_display(p_di)}, n={n_di}"),
        transform=ax.transAxes, ha="left", va="top", fontsize=6.5,
        family="monospace",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.9,
                  edgecolor="#ccc"),
    )

    out_rows = []
    for st in order:
        st_sub = sub.loc[sub["status_plot"].eq(st)]
        xs = st_sub[oral_col].to_numpy(float)
        ys = st_sub[gut_col].to_numpy(float)
        if tf is not None:
            xs, ys = tf(xs), tf(ys)
        rho, p, n = spearman_safe(xs, ys)
        out_rows.append({"metric": name, "status": st, "n_pairs": n,
                         "spearman_rho": rho, "p_value": p})

    # --- subgroup multiplicity: BH-FDR within this metric across the nine
    # subgroup statuses (All/HC pooled rows are primary, not in the family);
    # computed for the TSV only — NOT annotated on the figure (user request)
    sub_rows = [r for r in out_rows if r["status"] != "HC"]
    qs = bh_fdr(np.array([r["p_value"] for r in sub_rows]))
    for r, q in zip(sub_rows, qs):
        r["q_value_BH_within_metric"] = q

    return {"all": {"metric": name, "status": "All", "n_pairs": n_all,
                    "spearman_rho": rho_all, "p_value": p_all},
            "rows": out_rows}


def make_figure(df: pd.DataFrame, order, colors, tf, tf_label, subtitle,
                stem: str, density: bool = False,
                parkinson_short: bool = False) -> list:
    # canvas shrunk to near-minimum (user request: as small as possible)
    # — ~3.0-in square panels with fonts/scatter dropped proportionally
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.7))
    all_rows = []
    for ax, (name, oral_col, gut_col) in zip(axes, ALPHA_METRICS):
        rng = np.random.default_rng(2026)   # same draws for every panel
        res = draw_panel(ax, df, name, oral_col, gut_col, order, colors,
                         tf=tf, tf_label=tf_label, rng=rng,
                         density=density)
        all_rows.append(res["all"])
        all_rows.extend(res["rows"])

    if density:
        n_hc = int(df["status_plot"].eq("HC").sum())
        handles = [
            plt.Rectangle((0, 0), 1, 1, facecolor=HC_COLOR, alpha=0.25,
                          edgecolor="none",
                          label=f"HC (n={n_hc}), filled density"),
        ] + [
            plt.Line2D([], [], color=colors[st], lw=1.5,
                       label=f"{DISEASE_ABBREV.get(st, st)} "
                             f"(n={int((df['status_plot'] == st).sum())}), "
                             f"density outline")
            for st in order if st not in ("HC", "Other")
        ]
        handles.append(plt.Line2D([], [], color=colors["Other"], lw=1.5,
                                  label=f"Other (n="
                                        f"{int((df['status_plot'] == 'Other').sum())}), "
                                        f"density outline"))
    else:
        # user request (2026-10-02): the raw 3-panel figure must use the same
        # disease abbreviations as the per-metric raw figures, which shorten
        # "Parkinson's Disease" to "Parkinson"; every other label is
        # unchanged, so only the raw stem is rendered with this flag on.
        label_of = dict(DISEASE_ABBREV)
        if parkinson_short:
            label_of["Parkinson's Disease"] = "Parkinson"
        handles = [
            plt.Line2D([], [], marker="o", ls="", color=colors[st],
                       label=f"{label_of.get(st, st)} "
                             f"(n={int((df['status_plot'] == st).sum())})")
            for st in order
        ]
    handles += [
        plt.Line2D([], [], color=HC_LINE_COLOR, lw=2,
                   label="LOWESS — HC"),
        plt.Line2D([], [], color=DISEASE_LINE_COLOR, lw=2,
                   label="LOWESS — all disease"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=6, frameon=False,
               fontsize=7, bbox_to_anchor=(0.5, 1.06), columnspacing=0.8,
               handletextpad=0.2)
    fig.suptitle(
        "Single-pair Oral–Gut alpha diversity Spearman correlation "
        f"({subtitle}; line = LOWESS ± bootstrap 95% CI)",
        y=1.10, fontsize=9.5,
    )
    fig.tight_layout()
    fig.savefig(OUT_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)
    return all_rows


def make_metric_figures(df: pd.DataFrame, order, colors) -> None:
    """Save each raw-scale alpha-diversity metric as its own figure."""
    metric_stems = {
        "Observed": "observed",
        "Shannon": "shannon",
        "InvSimpson": "invsimpson",
    }
    for name, oral_col, gut_col in ALPHA_METRICS:
        # explicit inches-based unified geometry (user request 2026-10-20):
        # square 2.26-in panel, bottom edge 0.70 in, top edge 2.96 in, saved
        # at 300 dpi WITHOUT bbox_inches="tight" so the panel box sits at the
        # exact same pixel rect in all of 03/05/06/07/11; legend zone to the
        # right (panel right edge 2.88 in, legend starts 2.98 in), suptitle
        # anchored at a FIXED y in the shared 0.44-in top band
        fig = plt.figure(figsize=(FIG_W, CANVAS_H))
        ax = fig.add_axes(AX_RECT07)
        rng = np.random.default_rng(2026)
        draw_panel(ax, df, name, oral_col, gut_col, order, colors,
                   tf=None, tf_label="", rng=rng, label_metric=False)

        metric_labels = {
            st: "Parkinson" if st == "Parkinson's Disease" else DISEASE_ABBREV.get(st, st)
            for st in order
        }
        handles = [
            plt.Line2D([], [], marker="o", ls="", color=colors[st],
                       label=f"{metric_labels[st]} "
                             f"(n={int((df['status_plot'] == st).sum())})")
            for st in order
        ]
        # user request (2026-10-14): drop the "LOWESS —" prefix from the two
        # smooth-line entries and rename the disease line to
        # "all diseases combined"
        handles += [
            plt.Line2D([], [], color=HC_LINE_COLOR, lw=2, label="HC"),
            plt.Line2D([], [], color=DISEASE_LINE_COLOR, lw=2,
                       label="all diseases combined"),
        ]
        # user request (2026-10-06): legend to the RIGHT of the panel,
        # stacked vertically, flush against it; panel rect is fixed (see
        # AX_RECT07) so no subplots_adjust / anchored-suptitle dance is
        # needed any more — the suptitle sits at a fixed y in the shared
        # top band, ~6pt above the panel like the 05 beta title
        fig.legend(handles=handles, loc="center left", ncol=1,
                   frameon=False, bbox_to_anchor=((PANEL_IN + 0.62 + 0.10) / FIG_W, 0.5),
                   fontsize=7, handletextpad=0.2, borderaxespad=0)
        # user request (2026-10-14): title shortened to exactly
        # "alpha diversity Spearman correlation"
        fig.suptitle(
            "alpha diversity Spearman correlation",
            x=0.5, y=(AX_TOP_IN + 0.06) / CANVAS_H,
            fontsize=10, va="bottom",
        )
        stem = f"alpha_oral_gut_spearman_by_disease_raw_{metric_stems[name]}"
        fig.savefig(OUT_DIR / f"{stem}.png", dpi=300)  # no bbox_inches="tight"
        fig.savefig(OUT_DIR / f"{stem}.pdf")
        plt.close(fig)
        print(f"wrote {stem}.png and {stem}.pdf")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = load_pairs()
    order, colors = status_order_and_colors(df)
    print("status groups (pairs):",
          {s: int((df["status_plot"] == s).sum()) for s in order})

    # log-scale version (original re-render)
    rows_log = make_figure(df, order, colors, log10p, "  [log10(x+10⁻⁶)]",
                           "log10(x+10⁻⁶) scale",
                           "alpha_oral_gut_spearman_by_disease")
    # raw-scale version (no log, diseases still separated); its legend uses
    # "Parkinson" like the per-metric raw figures (parkinson_short=True)
    rows_raw = make_figure(df, order, colors, None, "", "raw scale",
                           "alpha_oral_gut_spearman_by_disease_raw",
                           parkinson_short=True)
    make_metric_figures(df, order, colors)

    # third version (user request 2026-09-30): remove zero-valued samples
    # first — the whole PAIR goes (see drop_zero_pairs), then take a bare
    # log10 with NO pseudo-count; the two versions above stay untouched
    df_nz = drop_zero_pairs(df)
    order_nz, colors_nz = status_order_and_colors(df_nz)
    print("status groups after zero-pair filter:",
          {s: int((df_nz["status_plot"] == s).sum()) for s in order_nz})
    rows_nz = make_figure(df_nz, order_nz, colors_nz, log10_plain, "  [log10]",
                          "log10 scale, zero-value pairs excluded",
                          "alpha_oral_gut_spearman_by_disease_log10_dropzero")
    print("drop-zero version, All-pairs Spearman:")
    for r in rows_nz:
        if r["status"] == "All":
            print(f"  {r['metric']}: rho={r['spearman_rho']:.3f} "
                  f"p={r['p_value']:.2e} n={r['n_pairs']}")
    stats_nz = pd.DataFrame(rows_nz)
    stats_nz = stats_nz[["metric", "status", "n_pairs", "spearman_rho",
                         "p_value", "q_value_BH_within_metric"]]
    stats_nz.to_csv(
        OUT_DIR / "alpha_oral_gut_spearman_by_status_log10_dropzero.tsv",
        sep="\t", index=False)

    # density render: use the original zero-filtered metadata values without
    # a log transform; the KDE and LOWESS styling stays unchanged
    make_figure(df_nz, order_nz, colors_nz, None, "",
                "raw scale, zero-value pairs excluded (density)",
                "alpha_oral_gut_spearman_by_disease_log10_dropzero_density",
                density=True)

    # Spearman is rank-based: values identical on raw and log10(x+10⁻⁶) scales
    stats = pd.DataFrame(rows_raw)
    stats = stats[["metric", "status", "n_pairs", "spearman_rho", "p_value",
                   "q_value_BH_within_metric"]]
    stats.to_csv(OUT_DIR / "alpha_oral_gut_spearman_by_status.tsv", sep="\t",
                 index=False)
    print(stats.to_string(index=False))
    print(f"outputs -> {OUT_DIR}")


if __name__ == "__main__":
    main()