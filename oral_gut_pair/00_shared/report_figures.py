#!/usr/bin/env python3
"""Lightweight report figures for oral–gut pair sections (cu01-friendly)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Patch, PathPatch
from matplotlib.path import Path as MplPath
import numpy as np
import pandas as pd
from scipy import stats

ALPHA_METRICS = [
    ("Observed", "oral_Observed", "gut_Observed", "delta_Observed"),
    ("Shannon", "oral_Shannon", "gut_Shannon", "delta_Shannon"),
    ("InvSimpson", "oral_InvSimpson", "gut_InvSimpson", "delta_InvSimpson"),
]

METRIC_DISPLAY = {
    "delta_alpha_observed": "Δ Observed",
    "delta_alpha_shannon": "Δ Shannon",
    "delta_alpha_invsimpson": "Δ InvSimpson",
    "beta_braycurtis_distance": "Bray–Curtis",
    "beta_weighted_unifrac_distance": "Weighted UniFrac",
    "beta_unweighted_unifrac_distance": "Unweighted UniFrac",
    "shared_otu_count": "Shared OTU count",
    "shared_otu_count_fraction_of_gut": "Shared OTU / gut",
    "shared_gut_abundance_fraction": "Shared gut abundance",
    "sourcetracker_paired_oral_fraction": "ST2 paired oral fraction",
}

STATUS_GROUP_ORDER = ["HC", "HC+disease", "disease_study"]
STATUS_GROUP_COLORS = {
    "HC": "#3182bd",
    "HC+disease": "#e6550d",
    "disease_study": "#de2d26",
}
ANALYSIS_SET_ORDER = ["single_pair", "time_series"]
ANALYSIS_SET_COLORS = {
    "single_pair": "#3182bd",
    "time_series": "#e6550d",
}
ANALYSIS_SET_LABELS = {
    "single_pair": "Single-pair",
    "time_series": "Time-series",
}
AGE_GROUP_ORDER = [
    "0-6 months",
    "6-36 months",
    "3-12 years",
    "12-18 years",
    "18-60 years",
    ">60 years",
    "Unknown",
]
AGE_GROUP_ALIASES = {
    "＞60 years": ">60 years",
    ">60 years": ">60 years",
}
PAIR_STATUS_ORDER = ["HC", "disease", "Unknown"]
PAIR_STATUS_COLORS = {
    "HC": "#3182bd",
    "disease": "#de2d26",
    "Unknown": "#969696",
}
ALPHA_DELTA_STATUS_METRICS = [
    ("Observed", "delta_alpha_observed"),
    ("Shannon", "delta_alpha_shannon"),
    ("InvSimpson", "delta_alpha_invsimpson"),
]
# Display abbreviations for disease/status labels (user-supplied mapping,
# 2026-09-30; applied in all new-render figures via the label helpers
# _short_disease_phrase / _bp_disease_label). Statuses whose abbreviation is
# the name itself (Asthma, Celiac Disease, Depression, Diabetes, Obesity,
# Pancreatic Cancer, Parkinson's Disease) and statuses never abbreviated by
# the user (e.g. "COVID-19 discharge") are omitted and pass through verbatim;
# the two long phrase statuses keep their earlier short forms.
DISEASE_LABEL_SHORT = {
    "Acute Myeloid Leukemia": "AML",
    "Alcohol Use Disorder": "AUD",
    "Alzheimer's Disease": "AD",
    "Aphthous": "AS",
    "Atopic Dermatitis": "AtD",
    "Autoimmune Hepatitis": "AIH",
    "Axial Spondyloarthritis": "axSpA",
    "Breast Cancer": "BC",
    "Chagas Disease": "CD",
    "Chronic Obstructive Pulmonary Disease": "COPD",
    "Clostridium Difficile Infection": "CDI",
    "Colorectal Cancer": "CRC",
    "Coronavirus Disease 2019": "COVID-19",
    "Dental Caries": "DC",
    "Digestive Tract Cancer": "DTC",
    "Esophageal Cancer": "EC",
    "Gastroesophageal Reflux Disease": "GERD",
    "Hand, Foot and Mouth Disease": "HFMD",
    "Head and Neck Cancer": "HNC",
    "Helicobacter pylori Infection": "HPI",
    "Human Immunodeficiency Virus": "HIV",
    "Hypertension": "HTN",
    "Inflammatory Bowel Disease": "IBD",
    "Intensive Care Unit": "ICU",
    "Lung Cancer": "LC",
    "Non-Carious Cervical Lesions": "NCCL",
    "Nonmalignant Breast Disease": "NBD",
    "Obstructive Sleep Apnea": "OSA",
    "Oral Lichen Planus": "OLP",
    "Oral Mucositis": "OM",
    "Oral Squamous Cell Carcinoma": "OSCC",
    "Oral Submucous Fibrosis": "OSMF",
    "Oral Verrucous Hyperplasia": "OVH",
    "Periodontal Disease": "PD",
    "Post-Traumatic Stress Disorder": "PTSD",
    "Preterm": "PTB",
    "Prostate Cancer": "PC",
    "Rheumatoid Arthritis": "RA",
    "Sjögren's Syndrome": "SS",
    "Thyroid Cancer": "TC",
    "1-year COVID-19 convalescent": "COVID conval.",
    "ICU discharge or on ICU day 10": "ICU day10+",
}


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def savefig(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {path}")


def load_delta(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if "dataset_label" in df.columns:
        prefer = df.loc[df["dataset_label"].eq("filter_ra")]
        if not prefer.empty:
            return prefer.copy()
    return df


def plot_alpha_pooled(delta: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.2))
    for ax, (name, oral_col, gut_col, delta_col) in zip(axes, ALPHA_METRICS):
        oral = delta[oral_col].astype(float)
        gut = delta[gut_col].astype(float)
        d = delta[delta_col].astype(float)
        keep = oral.notna() & gut.notna()
        oral, gut, d = oral[keep], gut[keep], d[keep]
        data = [oral.values, gut.values]
        ax.boxplot(
            data,
            labels=["Oral", "Gut"],
            widths=0.55,
            showfliers=False,
            patch_artist=True,
            boxprops=dict(facecolor="#9ecae1", alpha=0.8),
            medianprops=dict(color="#08519c", linewidth=1.5),
        )
        # light paired jitter overlay of a subsample for visual density
        rng = np.random.default_rng(2026)
        n = min(400, len(oral))
        idx = rng.choice(len(oral), size=n, replace=False)
        for i, vals in enumerate([oral.values[idx], gut.values[idx]], start=1):
            ax.scatter(
                np.full(n, i) + rng.uniform(-0.08, 0.08, n),
                vals,
                s=4,
                alpha=0.12,
                color="#3182bd",
                linewidths=0,
            )
        try:
            stat, p = stats.wilcoxon(gut, oral, zero_method="wilcox", alternative="two-sided")
        except ValueError:
            p = np.nan
            stat = np.nan
        ax.axhline(np.nanmedian(oral), color="#bdbdbd", linestyle=":", linewidth=0.8)
        ax.set_title(name)
        ax.set_ylabel(name)
        ax.text(
            0.5,
            0.98,
            f"Δ median={np.nanmedian(d):.3g}\nWilcoxon p={p:.2e}",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=8,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.85, edgecolor="#ccc"),
        )
    fig.suptitle("Paired Oral vs Gut alpha diversity (filter_ra)", y=1.02, fontsize=12)
    fig.tight_layout()
    savefig(fig, out / "delta_alpha_pooled_box.png")


def plot_alpha_by_bp(delta: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(max(10, delta["bioproject"].nunique() * 0.32), 10), sharex=False)
    for ax, (name, _, _, delta_col) in zip(axes, ALPHA_METRICS):
        order = (
            delta.groupby("bioproject")[delta_col]
            .median()
            .sort_values()
            .index.tolist()
        )
        data = [delta.loc[delta["bioproject"].eq(bp), delta_col].dropna().astype(float).values for bp in order]
        ax.boxplot(
            data,
            labels=order,
            showfliers=False,
            widths=0.6,
            patch_artist=True,
            boxprops=dict(facecolor="#c6dbef", alpha=0.85),
            medianprops=dict(color="#08519c", linewidth=1.2),
        )
        ax.axhline(0, color="#d62728", linestyle="--", linewidth=1)
        ax.set_ylabel(f"Δ {name}\n(Gut − Oral)")
        ax.tick_params(axis="x", labelrotation=75, labelsize=7)
    fig.suptitle("Per-BioProject Δ alpha (Gut − Oral)", y=1.01, fontsize=12)
    fig.tight_layout()
    savefig(fig, out / "delta_alpha_by_bioproject.png")


# ---------------------------------------------------------------------------
# Single-pair Δ alpha stratified by disease status
# ---------------------------------------------------------------------------

def _cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """Cliff's δ = P(a > b) − P(a < b). Positive means a tends larger than b."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    diff = a[:, None] - b[None, :]
    return float(np.sign(diff).sum() / (len(a) * len(b)))


def _hodges_lehmann_shift(a: np.ndarray, b: np.ndarray) -> float:
    """Median of all pairwise (a − b) differences."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    return float(np.median(a[:, None] - b[None, :]))


def _bootstrap_cliffs_ci(
    a: np.ndarray,
    b: np.ndarray,
    n_boot: int = 2000,
    seed: int = 2026,
    alpha: float = 0.05,
) -> tuple[float, float]:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    deltas = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        aa = rng.choice(a, size=len(a), replace=True)
        bb = rng.choice(b, size=len(b), replace=True)
        deltas[i] = _cliffs_delta(aa, bb)
    lo, hi = np.quantile(deltas, [alpha / 2, 1 - alpha / 2])
    return float(lo), float(hi)


def _bh_fdr(pvals: list[float] | np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    q = np.full(n, np.nan)
    valid = np.isfinite(p)
    if not valid.any():
        return q
    idx = np.where(valid)[0]
    order = idx[np.argsort(p[idx])]
    ranked = p[order]
    m = len(ranked)
    adj = ranked * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    q[order] = adj
    return q


def _dersimonian_laird(
    yi: np.ndarray, sei: np.ndarray
) -> dict:
    """Random-effects pooled estimate with DerSimonian–Laird τ² and I²."""
    yi = np.asarray(yi, dtype=float)
    sei = np.asarray(sei, dtype=float)
    keep = np.isfinite(yi) & np.isfinite(sei) & (sei > 0)
    yi, sei = yi[keep], sei[keep]
    k = len(yi)
    if k == 0:
        return {
            "k": 0,
            "pooled": float("nan"),
            "se": float("nan"),
            "ci_low": float("nan"),
            "ci_high": float("nan"),
            "tau2": float("nan"),
            "Q": float("nan"),
            "I2": float("nan"),
        }
    vi = sei ** 2
    wi = 1.0 / vi
    mean_fe = float(np.sum(wi * yi) / np.sum(wi))
    Q = float(np.sum(wi * (yi - mean_fe) ** 2))
    df = k - 1
    c = float(np.sum(wi) - np.sum(wi ** 2) / np.sum(wi)) if k > 1 else float("nan")
    tau2 = max(0.0, (Q - df) / c) if k > 1 and np.isfinite(c) and c > 0 else 0.0
    wi_star = 1.0 / (vi + tau2)
    pooled = float(np.sum(wi_star * yi) / np.sum(wi_star))
    se = float(np.sqrt(1.0 / np.sum(wi_star)))
    i2 = max(0.0, (Q - df) / Q) * 100.0 if k > 1 and Q > 0 else 0.0
    return {
        "k": int(k),
        "pooled": pooled,
        "se": se,
        "ci_low": pooled - 1.96 * se,
        "ci_high": pooled + 1.96 * se,
        "tau2": float(tau2),
        "Q": Q,
        "I2": float(i2),
    }


def _pair_status_binary(status_cleaned) -> str:
    if pd.isna(status_cleaned):
        return "Unknown"
    s = normalize_status_value(status_cleaned)
    if s == "NA":
        return "Unknown"
    return "HC" if is_hc_status(s) else "disease"


def _bp_disease_label(sub: pd.DataFrame) -> str:
    """Most common non-HC status; append '+k' if multiple disease labels."""
    disease = sub.loc[sub["pair_status"].eq("disease"), "status_cleaned"]
    disease = disease.dropna().map(normalize_status_value)
    disease = disease[disease != "NA"]
    if disease.empty:
        return ""
    counts = disease.value_counts()
    top = str(counts.index[0])
    top = DISEASE_LABEL_SHORT.get(top, top)
    if len(top) > 20:   # 20 fits "Parkinson's Disease" untruncated
        top = top[:18] + "…"
    if len(counts) > 1:
        return f"{top} +{len(counts) - 1}"
    return top


def _sig_stars(q: float) -> str:
    if not np.isfinite(q):
        return ""
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return ""


def build_alpha_status_table(pair_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Annotate single_pair rows and compute within-BP HC vs disease tests."""
    sp = pair_df.loc[pair_df["analysis_set"].astype(str).eq("single_pair")].copy()
    if sp.empty:
        raise ValueError("No single_pair rows in pair metrics")

    sp["pair_status"] = sp["status_cleaned"].map(_pair_status_binary)
    bp_status = (
        sp.groupby("bioproject", sort=False)["status_cleaned"]
        .apply(lambda s: classify_bp_status_group(s.tolist()))
        .apply(pd.Series)
        .rename(columns={0: "status_rank", 1: "status_group"})
    )
    sp = sp.join(bp_status, on="bioproject")

    disease_labels = {}
    for bp, sub in sp.groupby("bioproject"):
        disease_labels[bp] = _bp_disease_label(sub)
    sp["disease_label"] = sp["bioproject"].map(disease_labels)

    # Fixed BP order within each status_group by median Δ Observed
    order_metric = "delta_alpha_observed"
    bp_order_parts: list[str] = []
    for sg in STATUS_GROUP_ORDER:
        bps = (
            sp.loc[sp["status_group"].eq(sg)]
            .groupby("bioproject")[order_metric]
            .median()
            .sort_values()
            .index.tolist()
        )
        bp_order_parts.extend(bps)
    sp["bp_order"] = sp["bioproject"].map({bp: i for i, bp in enumerate(bp_order_parts)})

    test_rows = []
    for metric_name, col in ALPHA_DELTA_STATUS_METRICS:
        for bp, sub in sp.groupby("bioproject"):
            if sub["status_group"].iloc[0] != "HC+disease":
                continue
            hc = sub.loc[sub["pair_status"].eq("HC"), col].astype(float).dropna().values
            dz = sub.loc[sub["pair_status"].eq("disease"), col].astype(float).dropna().values
            if len(hc) == 0 or len(dz) == 0:
                continue
            try:
                mwu_p = float(stats.mannwhitneyu(dz, hc, alternative="two-sided").pvalue)
            except ValueError:
                mwu_p = float("nan")
            cd = _cliffs_delta(dz, hc)
            lo, hi = _bootstrap_cliffs_ci(dz, hc, seed=abs(hash((bp, col))) % (2**31 - 1) or 2026)
            test_rows.append(
                {
                    "bioproject": bp,
                    "status_group": "HC+disease",
                    "disease_label": disease_labels.get(bp, ""),
                    "metric": metric_name,
                    "metric_col": col,
                    "n_hc": int(len(hc)),
                    "n_disease": int(len(dz)),
                    "median_hc": float(np.median(hc)),
                    "median_disease": float(np.median(dz)),
                    "hodges_lehmann_shift": _hodges_lehmann_shift(dz, hc),
                    "cliffs_delta": cd,
                    "ci_low": lo,
                    "ci_high": hi,
                    "mwu_p": mwu_p,
                }
            )
    tests_df = pd.DataFrame(test_rows)
    if not tests_df.empty:
        tests_df["q_BH"] = np.nan
        for metric_name, _ in ALPHA_DELTA_STATUS_METRICS:
            mask = tests_df["metric"].eq(metric_name)
            tests_df.loc[mask, "q_BH"] = _bh_fdr(tests_df.loc[mask, "mwu_p"].tolist())
        tests_df["sig_stars"] = tests_df["q_BH"].map(_sig_stars)

    pooled: dict = {}
    for metric_name, _ in ALPHA_DELTA_STATUS_METRICS:
        sub = tests_df.loc[tests_df["metric"].eq(metric_name)] if not tests_df.empty else tests_df
        if sub.empty:
            pooled[metric_name] = _dersimonian_laird(np.array([]), np.array([]))
            continue
        yi = sub["cliffs_delta"].to_numpy(dtype=float)
        half = (sub["ci_high"] - sub["ci_low"]).to_numpy(dtype=float) / (2 * 1.96)
        pooled[metric_name] = _dersimonian_laird(yi, half)
        pooled[metric_name]["n_bioprojects"] = int(sub["bioproject"].nunique())

    return sp, tests_df, pooled


def plot_alpha_by_bp_status(
    pair_annot: pd.DataFrame,
    tests_df: pd.DataFrame,
    out: Path,
) -> None:
    """3 metrics × 3 status_group facets; mixed BPs split HC/disease boxes."""
    # Count boxes per facet for width_ratios
    box_counts = {}
    for sg in STATUS_GROUP_ORDER:
        sub = pair_annot.loc[pair_annot["status_group"].eq(sg)]
        n_boxes = 0
        for bp, bp_sub in sub.groupby("bioproject"):
            statuses = [s for s in PAIR_STATUS_ORDER if (bp_sub["pair_status"] == s).any()]
            n_boxes += max(1, len(statuses))
        box_counts[sg] = max(1, n_boxes)
    width_ratios = [box_counts[sg] for sg in STATUS_GROUP_ORDER]

    fig = plt.figure(figsize=(max(18.5, sum(width_ratios) * 0.55 + 3), 13.2))
    gs = fig.add_gridspec(3, 3, width_ratios=width_ratios, wspace=0.12, hspace=0.42)

    # BP order fixed by Δ Observed median within each status_group
    bp_orders: dict[str, list[str]] = {}
    for sg in STATUS_GROUP_ORDER:
        sub = pair_annot.loc[pair_annot["status_group"].eq(sg)]
        bp_orders[sg] = (
            sub.groupby("bioproject")["delta_alpha_observed"]
            .median()
            .sort_values()
            .index.tolist()
        )

    sig_map = {}
    if tests_df is not None and not tests_df.empty:
        for _, row in tests_df.iterrows():
            sig_map[(row["bioproject"], row["metric"])] = row

    legend_handles = [
        Patch(facecolor=PAIR_STATUS_COLORS["HC"], edgecolor="#333", label="HC"),
        Patch(facecolor=PAIR_STATUS_COLORS["disease"], edgecolor="#333", label="disease"),
        Patch(facecolor=PAIR_STATUS_COLORS["Unknown"], edgecolor="#333", label="Unknown"),
        Patch(
            facecolor="#fff3cd",
            edgecolor="#8b0000",
            linewidth=1.5,
            label="HC vs disease significant (BH q<0.05)",
        ),
    ]

    for row_i, (metric_name, col) in enumerate(ALPHA_DELTA_STATUS_METRICS):
        axes_row = []
        for col_i, sg in enumerate(STATUS_GROUP_ORDER):
            sharey = axes_row[0] if axes_row else None
            ax = fig.add_subplot(gs[row_i, col_i], sharey=sharey)
            axes_row.append(ax)

            sub = pair_annot.loc[pair_annot["status_group"].eq(sg)]
            order = bp_orders[sg]
            positions = []
            data = []
            colors = []
            ns = []
            box_bps = []  # bioproject for each box (for sig highlight)
            tick_pos = []
            tick_labels = []
            tick_bps = []
            bp_pos_map: dict[str, list[float]] = {}
            x = 1.0
            rng = np.random.default_rng(2026 + row_i * 17 + col_i)

            for bp in order:
                bp_sub = sub.loc[sub["bioproject"].eq(bp)]
                statuses = [s for s in PAIR_STATUS_ORDER if (bp_sub["pair_status"] == s).any()]
                if not statuses:
                    continue
                disease_lab = str(bp_sub["disease_label"].iloc[0] or "")
                group_xs = []
                for st in statuses:
                    vals = (
                        bp_sub.loc[bp_sub["pair_status"].eq(st), col]
                        .astype(float)
                        .dropna()
                        .values
                    )
                    positions.append(x)
                    data.append(vals)
                    colors.append(PAIR_STATUS_COLORS[st])
                    ns.append(len(vals))
                    box_bps.append(bp)
                    group_xs.append(x)
                    x += 0.9
                bp_pos_map[bp] = group_xs
                tick_pos.append(float(np.mean(group_xs)))
                tick_bps.append(bp)
                if disease_lab and sg != "HC":
                    tick_labels.append(f"{bp}\n({disease_lab})")
                elif "Unknown" in statuses:
                    tick_labels.append(f"{bp}\n(Unknown)")
                else:
                    tick_labels.append(bp)
                x += 0.55  # gap between BPs

            if not data:
                ax.set_visible(False)
                continue

            sig_bps = set()
            if sg == "HC+disease":
                for bp in tick_bps:
                    info = sig_map.get((bp, metric_name))
                    stars = ""
                    if info is not None:
                        stars = info.get("sig_stars", "") or ""
                        if isinstance(stars, float) and np.isnan(stars):
                            stars = ""
                    if stars:
                        sig_bps.add(bp)

            # light background band behind significant BP pairs
            if sig_bps:
                for bp in sig_bps:
                    xs = bp_pos_map[bp]
                    ax.axvspan(
                        min(xs) - 0.42,
                        max(xs) + 0.42,
                        facecolor="#fff3cd",
                        edgecolor="none",
                        alpha=0.85,
                        zorder=0,
                    )

            bp_artist = ax.boxplot(
                data,
                positions=positions,
                widths=0.55,
                showfliers=False,
                patch_artist=True,
                medianprops=dict(color="#222", linewidth=1.2),
                whiskerprops=dict(color="#555"),
                capprops=dict(color="#555"),
            )
            for patch, color, bp in zip(bp_artist["boxes"], colors, box_bps):
                patch.set_facecolor(color)
                patch.set_alpha(0.78)
                if bp in sig_bps:
                    patch.set_edgecolor("#8b0000")
                    patch.set_linewidth(2.0)
                    patch.set_zorder(4)
                else:
                    patch.set_edgecolor("#333")
                    patch.set_linewidth(0.8)

            # jittered points (clipped later via ylim)
            for pos, vals, color in zip(positions, data, colors):
                if len(vals) == 0:
                    continue
                n_show = min(len(vals), 80)
                idx = rng.choice(len(vals), size=n_show, replace=False)
                ax.scatter(
                    np.full(n_show, pos) + rng.uniform(-0.12, 0.12, n_show),
                    vals[idx],
                    s=8,
                    alpha=0.28,
                    color=color,
                    linewidths=0,
                    zorder=3,
                )

            # n labels
            y_vals_all = np.concatenate([v for v in data if len(v)])
            y_hi = float(np.nanpercentile(y_vals_all, 99)) if len(y_vals_all) else 0.0
            y_lo = float(np.nanpercentile(y_vals_all, 1)) if len(y_vals_all) else 0.0
            pad = (y_hi - y_lo) * 0.06 if y_hi > y_lo else 1.0
            for pos, n in zip(positions, ns):
                ax.text(
                    pos,
                    y_hi + pad * 0.15,
                    f"n={n}",
                    ha="center",
                    va="bottom",
                    fontsize=5,
                    color="#444",
                )

            # significance brackets for mixed BPs
            if sg == "HC+disease":
                for bp, xs in bp_pos_map.items():
                    if bp not in sig_bps or len(xs) < 2:
                        continue
                    info = sig_map.get((bp, metric_name))
                    stars = info.get("sig_stars", "") if info is not None else ""
                    x0, x1 = xs[0], xs[1]
                    y_br = y_hi + pad * 0.55
                    ax.plot(
                        [x0, x0, x1, x1],
                        [y_br - pad * 0.12, y_br, y_br, y_br - pad * 0.12],
                        color="#8b0000",
                        linewidth=1.3,
                        zorder=5,
                    )
                    ax.text(
                        (x0 + x1) / 2,
                        y_br + pad * 0.05,
                        stars,
                        ha="center",
                        va="bottom",
                        fontsize=9,
                        fontweight="bold",
                        color="#8b0000",
                        zorder=5,
                    )

            ax.axhline(0, color="#d62728", linestyle="--", linewidth=1, zorder=1)
            ax.set_xticks(tick_pos)
            tick_fs = 4.2 if sg == "HC" else 4.8
            tick_objs = ax.set_xticklabels(
                tick_labels, rotation=90, ha="center", va="top", fontsize=tick_fs
            )
            for lab, bp in zip(tick_objs, tick_bps):
                if bp in sig_bps:
                    lab.set_fontweight("bold")
                    lab.set_color("#8b0000")
                    lab.set_fontsize(tick_fs + 0.6)
            if col_i == 0:
                ax.set_ylabel(f"Δ {metric_name}\n(Gut − Oral)", fontsize=9)
            else:
                plt.setp(ax.get_yticklabels(), visible=False)
            if row_i == 0:
                ax.set_title(sg, fontsize=11, pad=6)
            ax.set_xlim(0.3, positions[-1] + 0.7)
            ax.tick_params(axis="x", pad=2)

        # shared y-limits for the row (1–99 pct of all row data + margin for labels)
        row_vals = pair_annot[col].astype(float).dropna().values
        if len(row_vals):
            lo = float(np.nanpercentile(row_vals, 1))
            hi = float(np.nanpercentile(row_vals, 99))
            span = hi - lo if hi > lo else 1.0
            for ax in axes_row:
                if ax.get_visible():
                    ax.set_ylim(lo - 0.08 * span, hi + 0.35 * span)

    fig.legend(
        handles=legend_handles,
        loc="upper center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, 0.998),
        fontsize=8,
    )
    fig.suptitle(
        "Single-pair Δ alpha by BioProject, stratified by health status",
        y=1.025,
        fontsize=12,
    )
    fig.subplots_adjust(bottom=0.14, top=0.92)
    savefig(fig, out / "delta_alpha_by_bioproject_status.png")


def plot_alpha_status_forest(
    tests_df: pd.DataFrame,
    pooled: dict,
    out: Path,
) -> None:
    """Forest plot of within-BP Cliff's δ (disease − HC) for mixed BioProjects."""
    if tests_df is None or tests_df.empty:
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.text(0.5, 0.5, "No mixed HC+disease BioProjects", ha="center", va="center")
        ax.axis("off")
        savefig(fig, out / "delta_alpha_hc_vs_disease_forest.png")
        return

    # Shared BP order by Observed cliffs_delta
    obs = tests_df.loc[tests_df["metric"].eq("Observed")].copy()
    if obs.empty:
        bp_order = sorted(tests_df["bioproject"].unique())
    else:
        bp_order = obs.sort_values("cliffs_delta")["bioproject"].tolist()
    # leave room at bottom for pooled diamond
    n_bp = len(bp_order)
    fig, axes = plt.subplots(1, 3, figsize=(13.5, max(4.5, 0.42 * n_bp + 1.8)), sharey=True)

    for ax, (metric_name, _) in zip(axes, ALPHA_DELTA_STATUS_METRICS):
        sub = tests_df.loc[tests_df["metric"].eq(metric_name)].set_index("bioproject")
        ys = []
        for i, bp in enumerate(bp_order):
            y = n_bp - i  # top = first
            ys.append(y)
            if bp not in sub.index:
                continue
            row = sub.loc[bp]
            cd = float(row["cliffs_delta"])
            lo = float(row["ci_low"])
            hi = float(row["ci_high"])
            n_min = min(int(row["n_hc"]), int(row["n_disease"]))
            size = 40 + 8 * n_min
            ax.errorbar(
                cd,
                y,
                xerr=[[cd - lo], [hi - cd]],
                fmt="none",
                ecolor="#555",
                elinewidth=1.2,
                capsize=2.5,
                zorder=2,
            )
            ax.scatter(
                cd,
                y,
                s=size,
                color="#de2d26",
                edgecolors="#333",
                linewidths=0.6,
                zorder=3,
                alpha=0.9,
            )
            label = bp
            dlab = str(row.get("disease_label", "") or "")
            if dlab:
                label = f"{bp} ({dlab})"
            ax.text(
                -1.02,
                y,
                f"{label}  [{int(row['n_hc'])}/{int(row['n_disease'])}]",
                ha="right",
                va="center",
                fontsize=7,
                clip_on=False,
            )

        # pooled diamond
        info = pooled.get(metric_name, {})
        pooled_est = info.get("pooled", float("nan"))
        if np.isfinite(pooled_est):
            y0 = 0.15
            half_h = 0.28
            lo = info["ci_low"]
            hi = info["ci_high"]
            ax.fill(
                [lo, pooled_est, hi, pooled_est],
                [y0, y0 + half_h, y0, y0 - half_h],
                color="#08519c",
                alpha=0.85,
                zorder=3,
            )
            ax.text(
                0.0,
                -0.55,
                f"RE pooled δ={pooled_est:.2f} "
                f"[{info['ci_low']:.2f}, {info['ci_high']:.2f}]\n"
                f"I²={info['I2']:.0f}%  (k={info['k']})",
                ha="center",
                va="top",
                fontsize=7.5,
                color="#08519c",
            )

        ax.axvline(0, color="#333", linestyle="--", linewidth=1, zorder=1)
        ax.set_xlim(-1.15, 1.15)
        ax.set_ylim(-1.1, n_bp + 0.7)
        ax.set_yticks([])
        ax.set_xlabel("Cliff's δ (disease − HC)", fontsize=9)
        ax.set_title(f"Δ {metric_name}", fontsize=11)
        ax.spines["left"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

    fig.suptitle(
        "Within-BioProject effect of disease on oral–gut Δ alpha\n"
        "(HC+disease BioProjects only; positive δ = larger Gut−Oral gap in disease)",
        fontsize=11,
        y=1.02,
    )
    fig.tight_layout()
    savefig(fig, out / "delta_alpha_hc_vs_disease_forest.png")


def make_alpha_status_figures(root: Path) -> None:
    pair_path = root / "results" / "cross_comparison" / "01_merged_metrics" / "pair_metrics_all_cohorts.tsv"
    if not pair_path.is_file():
        raise FileNotFoundError(f"Missing merged pair metrics: {pair_path}")
    out_dir = ensure_dir(root / "results" / "single_pair" / "02_alpha")
    plots = ensure_dir(out_dir / "plots")

    pair_df = pd.read_csv(pair_path, sep="\t")
    pair_annot, tests_df, pooled = build_alpha_status_table(pair_df)

    tests_path = out_dir / "delta_alpha_status_within_bp_tests.tsv"
    out_cols = [
        "bioproject",
        "metric",
        "n_hc",
        "n_disease",
        "median_hc",
        "median_disease",
        "hodges_lehmann_shift",
        "cliffs_delta",
        "ci_low",
        "ci_high",
        "mwu_p",
        "q_BH",
        "disease_label",
        "sig_stars",
    ]
    if tests_df.empty:
        pd.DataFrame(columns=out_cols).to_csv(tests_path, sep="\t", index=False)
    else:
        out_df = tests_df[out_cols].copy()
        out_df["sig_stars"] = out_df["sig_stars"].fillna("")
        out_df.to_csv(tests_path, sep="\t", index=False)

    pooled_path = out_dir / "delta_alpha_status_pooled_cliffs.json"
    with open(pooled_path, "w") as handle:
        json.dump(pooled, handle, indent=2)

    plot_alpha_by_bp_status(pair_annot, tests_df, plots)
    plot_alpha_status_forest(tests_df, pooled, plots)
    alpha_pairs = _load_single_pair_alpha_pairs(root)
    plot_alpha_oral_gut_spearman(alpha_pairs, plots)
    plot_alpha_oral_gut_spearman_by_bp(alpha_pairs, plots)
    print(
        f"alpha status figures: {pair_annot['bioproject'].nunique()} BP, "
        f"{len(pair_annot)} pairs, {0 if tests_df.empty else tests_df['bioproject'].nunique()} mixed BP "
        f"-> {out_dir}"
    )


def _spearman_safe(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int]:
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = int(len(x))
    if n < 3:
        return float("nan"), float("nan"), n
    rho, p = stats.spearmanr(x, y)
    return float(rho), float(p), n


def _load_single_pair_alpha_pairs(root: Path) -> pd.DataFrame:
    delta_path = root / "results" / "single_pair" / "02_alpha" / "delta_alpha.tsv"
    meta_path = (
        root / "results" / "cross_comparison" / "01_merged_metrics" / "pair_metrics_all_cohorts.tsv"
    )
    delta = load_delta(delta_path)
    meta = pd.read_csv(meta_path, sep="\t")
    sp_meta = meta.loc[
        meta["analysis_set"].astype(str).eq("single_pair"),
        ["pair_id", "status_cleaned"],
    ].drop_duplicates("pair_id")
    df = delta.merge(sp_meta, on="pair_id", how="left")
    df["pair_status"] = df["status_cleaned"].map(_pair_status_binary)
    return df


def _bp_spearman_title(sub: pd.DataFrame) -> str:
    """Panel title: BioProject plus HC / disease composition."""
    bp = str(sub["bioproject"].iloc[0])
    _, status_group = classify_bp_status_group(sub["status_cleaned"].tolist())
    disease = _bp_disease_label(sub)
    has_unknown = bool(sub["pair_status"].eq("Unknown").any())
    if status_group == "HC":
        label = "HC; Unknown" if has_unknown else "HC"
    elif status_group == "disease_study":
        label = disease or "disease"
    else:
        label = f"HC + {disease}" if disease else "HC + disease"
    return f"{bp}\n{label}"


def _ordered_bioprojects(df: pd.DataFrame) -> list[str]:
    ranks = {}
    for bp, sub in df.groupby("bioproject", sort=False):
        rank, _ = classify_bp_status_group(sub["status_cleaned"].tolist())
        ranks[bp] = rank
    return sorted(df["bioproject"].astype(str).unique(), key=lambda bp: (ranks.get(bp, 9), bp))


def _scatter_oral_gut(ax, sub: pd.DataFrame, oral_col: str, gut_col: str, point_size: float) -> None:
    for st in ("Unknown", "HC", "disease"):
        st_sub = sub.loc[sub["pair_status"].eq(st)]
        if st_sub.empty:
            continue
        ax.scatter(
            st_sub[oral_col],
            st_sub[gut_col],
            s=point_size,
            color=PAIR_STATUS_COLORS[st],
            alpha=0.55 if point_size < 20 else 0.42,
            linewidths=0,
            zorder=3,
        )


def plot_alpha_oral_gut_spearman(df: pd.DataFrame, out: Path) -> None:
    """Pooled Oral vs Gut Spearman scatter, colored by HC / disease."""
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 4.8))
    status_plot_order = ["HC", "disease", "Unknown"]

    for ax, (name, oral_col, gut_col, _) in zip(axes, ALPHA_METRICS):
        sub = df[[oral_col, gut_col, "pair_status"]].copy()
        sub[oral_col] = pd.to_numeric(sub[oral_col], errors="coerce")
        sub[gut_col] = pd.to_numeric(sub[gut_col], errors="coerce")
        sub = sub.dropna(subset=[oral_col, gut_col])

        oral = sub[oral_col].to_numpy(dtype=float)
        gut = sub[gut_col].to_numpy(dtype=float)
        lo = float(np.nanmin([oral.min(), gut.min()]))
        hi = float(np.nanmax([oral.max(), gut.max()]))
        pad = (hi - lo) * 0.04 if hi > lo else 1.0
        lim = (lo - pad, hi + pad)

        rho_all, p_all, n_all = _spearman_safe(oral, gut)
        stat_lines = [f"All: ρ={rho_all:.3f}, p={p_all:.2e}, n={n_all}"]
        for st in ("HC", "disease"):
            st_sub = sub.loc[sub["pair_status"].eq(st)]
            rho, p, n = _spearman_safe(
                st_sub[oral_col].to_numpy(dtype=float),
                st_sub[gut_col].to_numpy(dtype=float),
            )
            if n:
                stat_lines.append(f"{st}: ρ={rho:.3f}, p={p:.2e}, n={n}")

        ax.plot(lim, lim, linestyle="--", color="#bdbdbd", linewidth=1.0, zorder=1)
        _scatter_oral_gut(ax, sub, oral_col, gut_col, point_size=14)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(f"Oral {name}")
        ax.set_ylabel(f"Gut {name}")
        ax.set_title(name, fontsize=11)
        ax.text(
            0.04,
            0.98,
            "\n".join(stat_lines),
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7.2,
            family="monospace",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.9, edgecolor="#ccc"),
        )

    handles = [
        Patch(facecolor=PAIR_STATUS_COLORS[st], edgecolor="#333", label=st)
        for st in status_plot_order
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 1.04),
        fontsize=9,
    )
    fig.suptitle(
        "Single-pair Oral–Gut alpha diversity Spearman correlation",
        y=1.10,
        fontsize=12,
    )
    fig.tight_layout()
    savefig(fig, out / "alpha_oral_gut_spearman.png")


def plot_alpha_oral_gut_spearman_by_bp(df: pd.DataFrame, out: Path) -> None:
    """One Oral–Gut Spearman panel per BioProject, faceted by metric."""
    bps = _ordered_bioprojects(df)
    n_bp = len(bps)
    ncols = 6
    nrows = int(np.ceil(n_bp / ncols))
    titles = {}
    for bp, sub in df.groupby("bioproject"):
        titles[str(bp)] = _bp_spearman_title(sub)

    fig, axes = plt.subplots(
        3 * nrows,
        ncols,
        figsize=(ncols * 2.45, 3 * nrows * 2.35 + 0.6),
        squeeze=False,
    )
    for metric_i, (name, oral_col, gut_col, _) in enumerate(ALPHA_METRICS):
        for i, bp in enumerate(bps):
            r = metric_i * nrows + i // ncols
            c = i % ncols
            ax = axes[r, c]
            sub = df.loc[df["bioproject"].astype(str).eq(bp), [oral_col, gut_col, "pair_status"]].copy()
            sub[oral_col] = pd.to_numeric(sub[oral_col], errors="coerce")
            sub[gut_col] = pd.to_numeric(sub[gut_col], errors="coerce")
            sub = sub.dropna(subset=[oral_col, gut_col])
            if sub.empty:
                ax.axis("off")
                continue
            oral = sub[oral_col].to_numpy(dtype=float)
            gut = sub[gut_col].to_numpy(dtype=float)
            lo = float(min(oral.min(), gut.min()))
            hi = float(max(oral.max(), gut.max()))
            pad = (hi - lo) * 0.08 if hi > lo else 1.0
            lim = (lo - pad, hi + pad)
            ax.plot(lim, lim, linestyle="--", color="#bdbdbd", linewidth=0.7, zorder=1)
            _scatter_oral_gut(ax, sub, oral_col, gut_col, point_size=10)
            ax.set_xlim(lim)
            ax.set_ylim(lim)
            rho, p, n = _spearman_safe(oral, gut)
            rho_txt = "ρ=NA" if not np.isfinite(rho) else f"ρ={rho:.2f}"
            p_txt = "" if not np.isfinite(p) else f", p={p:.1e}"
            ax.set_title(f"{titles[bp]}\n{name}  {rho_txt}{p_txt}  n={n}", fontsize=6.2, pad=2)
            ax.tick_params(labelsize=5.5)
            if c == 0:
                ax.set_ylabel(f"Gut {name}", fontsize=6.5)
            if i // ncols == nrows - 1:
                ax.set_xlabel(f"Oral {name}", fontsize=6.5)
        for i in range(n_bp, nrows * ncols):
            axes[metric_i * nrows + i // ncols, i % ncols].axis("off")

    handles = [
        Patch(facecolor=PAIR_STATUS_COLORS[st], edgecolor="#333", label=st)
        for st in ("HC", "disease", "Unknown")
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 1.005),
        fontsize=9,
    )
    fig.suptitle(
        "Per-BioProject Oral–Gut alpha Spearman correlation (single-pair)",
        y=1.012,
        fontsize=13,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.992))
    savefig(fig, out / "alpha_oral_gut_spearman_by_bioproject.png")


def plot_st2_pooled(st2: pd.DataFrame, out: Path) -> None:
    vals = st2["paired_oral_fraction"].astype(float).dropna()
    vals_plot = vals.clip(lower=1e-6)
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.hist(np.log10(vals_plot), bins=40, color="#6baed6", edgecolor="white", alpha=0.9)
    med, mean = float(vals.median()), float(vals.mean())
    ax.axvline(np.log10(max(med, 1e-6)), color="#e6550d", linestyle="--", label=f"median={med:.4g}")
    ax.axvline(np.log10(max(mean, 1e-6)), color="#31a354", linestyle=":", label=f"mean={mean:.4g}")
    ax.set_xlabel("log10(paired oral fraction)")
    ax.set_ylabel("Number of gut sinks")
    ax.set_title("SourceTracker2 paired oral fraction (pooled)")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    savefig(fig, out / "st2_paired_oral_fraction_pooled.png")


def plot_st2_by_bp(st2: pd.DataFrame, out: Path) -> None:
    order = (
        st2.groupby("bioproject")["paired_oral_fraction"]
        .median()
        .sort_values()
        .index.tolist()
    )
    data = [
        st2.loc[st2["bioproject"].eq(bp), "paired_oral_fraction"].dropna().astype(float).values
        for bp in order
    ]
    fig, ax = plt.subplots(figsize=(max(10, len(order) * 0.35), 4.8))
    ax.boxplot(
        data,
        labels=order,
        showfliers=False,
        widths=0.6,
        patch_artist=True,
        boxprops=dict(facecolor="#c7e9c0", alpha=0.85),
        medianprops=dict(color="#006d2c", linewidth=1.2),
    )
    ax.set_yscale("log")
    ax.set_ylabel("Paired oral fraction")
    ax.set_title("SourceTracker2 paired oral fraction by BioProject")
    ax.tick_params(axis="x", labelrotation=75, labelsize=7)
    fig.tight_layout()
    savefig(fig, out / "st2_paired_oral_fraction_by_bioproject.png")


def plot_dispersion(disp: pd.DataFrame, out: Path) -> None:
    # Prefer per-bioproject rows; fall back to all non-pooled if needed.
    df = disp.copy()
    if "dispersion_scope" in df.columns:
        per = df.loc[df["dispersion_scope"].astype(str).str.contains("per_bioproject", case=False, na=False)]
        if not per.empty:
            df = per
    if "test_bioproject" in df.columns:
        df = df.loc[~df["test_bioproject"].astype(str).isin(["ALL", "pooled", "nan"])].copy()
        bp_col = "test_bioproject"
    else:
        bp_col = "bioproject" if "bioproject" in df.columns else None
    if bp_col is None or df.empty:
        print("skip dispersion plot: no per-bioproject rows")
        return

    metrics = sorted(df["metric"].dropna().unique().tolist())
    fig, axes = plt.subplots(1, len(metrics), figsize=(4.2 * len(metrics), 5.2), sharey=True)
    if len(metrics) == 1:
        axes = [axes]
    for ax, metric in zip(axes, metrics):
        sub = df.loc[df["metric"].eq(metric)].copy()
        sub = sub.sort_values("median_delta_gut_minus_oral")
        y = np.arange(len(sub))
        colors = []
        for _, row in sub.iterrows():
            p = row.get("p_signflip_two_sided", np.nan)
            colors.append("#e6550d" if pd.notna(p) and p < 0.05 else "#6baed6")
        ax.scatter(sub["median_delta_gut_minus_oral"], y, c=colors, s=28, zorder=3)
        ax.axvline(0, color="#d62728", linestyle="--", linewidth=1)
        ax.set_yticks(y)
        ax.set_yticklabels(sub[bp_col].astype(str), fontsize=7)
        ax.set_xlabel("median Δ dispersion\n(Gut − Oral)")
        ax.set_title(metric)
    fig.suptitle("Paired dispersion delta by BioProject (orange = sign-flip p < 0.05)", y=1.02, fontsize=11)
    fig.tight_layout()
    savefig(fig, out / "dispersion_delta_by_metric.png")


def plot_concordance(conc: pd.DataFrame, desc: pd.DataFrame, out: Path) -> None:
    ensure_dir(out)
    order = list(METRIC_DISPLAY.keys())
    conc = conc.set_index("metric").reindex([m for m in order if m in conc["metric"].values]).reset_index()
    y = np.arange(len(conc))
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    deltas = conc["cliffs_delta_single_minus_time"].astype(float).values
    qs = conc["q_between_designs"].astype(float).values
    colors = ["#e6550d" if q < 0.05 else "#6baed6" for q in qs]
    ax.axvline(0, color="#969696", linestyle="--", linewidth=1)
    ax.hlines(y, 0, deltas, color="#bdbdbd", linewidth=1)
    ax.scatter(deltas, y, c=colors, s=55, zorder=3)
    labels = [METRIC_DISPLAY.get(m, m) for m in conc["metric"]]
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Cliff's δ (single_pair − time_series)")
    ax.set_title("Cross-design effect sizes\n(orange = BH q_between < 0.05)")
    for yi, q, d in zip(y, qs, deltas):
        ax.text(d + (0.01 if d >= 0 else -0.01), yi, f"q={q:.2g}", va="center", ha="left" if d >= 0 else "right", fontsize=7)
    fig.tight_layout()
    savefig(fig, out / "concordance_cliffs_delta_forest.png")

    # median lollipop: normalize each metric by pooled IQR of single_pair for visual scale
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    med_s = conc["median_single_pair"].astype(float).values
    med_t = conc["median_time_series"].astype(float).values
    # scale by max(|med|) per metric to put on common axis for display
    scale = np.maximum(np.maximum(np.abs(med_s), np.abs(med_t)), 1e-12)
    xs = med_s / scale
    xt = med_t / scale
    ax.hlines(y, xs, xt, color="#bdbdbd", linewidth=1.5)
    ax.scatter(xs, y, color="#3182bd", s=45, label="single_pair (scaled)", zorder=3)
    ax.scatter(xt, y, color="#e6550d", s=45, label="time_series (scaled)", zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Median / max(|median| of two designs)")
    ax.set_title("Median metric values: single_pair vs time_series")
    ax.legend(frameon=False, loc="lower right", fontsize=8)
    # annotate raw medians
    for yi, ms, mt in zip(y, med_s, med_t):
        ax.text(1.02, yi, f"{ms:.3g} | {mt:.3g}", transform=ax.get_yaxis_transform(), va="center", fontsize=7, color="#636363")
    fig.tight_layout()
    savefig(fig, out / "median_single_vs_time.png")


# ---------------------------------------------------------------------------
# Cohort overview (dataset description)
# ---------------------------------------------------------------------------

def normalize_status_value(x) -> str:
    x = str(x).strip()
    if x == "" or x.lower() in {"nan", "none", "na", "null"}:
        return "NA"
    return x


def is_hc_status(x) -> bool:
    s = str(x).strip().lower()
    return s in {
        "hc", "h", "healthy", "health", "control", "controls", "normal",
        "healthy_control", "healthy control", "non-disease", "non_disease",
        "non disease", "non-diseased", "non_diseased", "non diseased",
        "noncase", "non-case", "non_case",
    }


def classify_bp_status_group(status_values) -> tuple[int, str]:
    values = [normalize_status_value(x) for x in status_values]
    values = [x for x in values if x != "NA"]
    if not values:
        return 2, "disease_study"
    has_hc = any(is_hc_status(x) for x in values)
    has_disease = any(not is_hc_status(x) for x in values)
    if has_hc and not has_disease:
        return 0, "HC"
    if has_hc and has_disease:
        return 1, "HC+disease"
    return 2, "disease_study"


def normalize_age_group(x) -> str:
    if pd.isna(x):
        return "Unknown"
    s = str(x).strip()
    if s == "" or s.lower() in {"nan", "none", "na", "null"}:
        return "Unknown"
    return AGE_GROUP_ALIASES.get(s, s)


def normalize_continent(x) -> str:
    if pd.isna(x):
        return "Unknown"
    s = str(x).strip()
    if s == "" or s.lower() in {"nan", "none", "na", "null"}:
        return "Unknown"
    return s


def _safe_unique_join(values) -> str:
    vals = []
    for v in values:
        s = normalize_status_value(v)
        if s != "NA":
            vals.append(s)
    return ";".join(sorted(set(vals))) if vals else "NA"


def build_cohort_overview_tables(pair_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    df = pair_df.copy()
    required = {
        "bioproject",
        "pair_id",
        "analysis_set",
        "status_cleaned",
        "continent_group",
        "age_cleaned_group",
    }
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"pair table missing columns: {missing}")

    bp_rows = []
    for bp, sub in df.groupby("bioproject", sort=False):
        statuses = sub["status_cleaned"].tolist()
        status_rank, status_group = classify_bp_status_group(statuses)
        status_values = _safe_unique_join(statuses)
        analysis_sets = sorted(set(sub["analysis_set"].astype(str)))
        if len(analysis_sets) != 1:
            raise ValueError(f"BioProject {bp} spans multiple analysis_set values: {analysis_sets}")
        analysis_set = analysis_sets[0]
        continents = sorted({normalize_continent(x) for x in sub["continent_group"]})
        n_pairs = int(len(sub))
        n_subjects = int(sub["subject_id"].nunique()) if "subject_id" in sub.columns else np.nan
        bp_rows.append(
            {
                "bioproject": bp,
                "analysis_set": analysis_set,
                "status_group": status_group,
                "status_rank": status_rank,
                "status_values": status_values,
                "bp_status_label": f"{bp}({status_values})",
                "n_pairs": n_pairs,
                "n_subjects": n_subjects,
                "n_samples": n_pairs * 2,
                "continent_values": ";".join(continents),
                "primary_continent": max(
                    continents,
                    key=lambda c: int((sub["continent_group"].map(normalize_continent) == c).sum()),
                ),
            }
        )
    bp_summary = pd.DataFrame(bp_rows).sort_values(
        ["status_rank", "analysis_set", "n_pairs", "bioproject"],
        ascending=[True, True, False, True],
    ).reset_index(drop=True)

    pair_annot = df.merge(
        bp_summary[["bioproject", "status_group", "status_rank", "status_values", "bp_status_label"]],
        on="bioproject",
        how="left",
        validate="many_to_one",
    )
    pair_annot["age_group"] = pair_annot["age_cleaned_group"].map(normalize_age_group)
    pair_annot["continent"] = pair_annot["continent_group"].map(normalize_continent)

    n_pairs = int(len(pair_annot))
    n_bp = int(pair_annot["bioproject"].nunique())
    by_set = pair_annot["analysis_set"].value_counts().to_dict()
    by_status_pairs = pair_annot["status_group"].value_counts().reindex(STATUS_GROUP_ORDER).fillna(0).astype(int).to_dict()
    by_status_bp = bp_summary["status_group"].value_counts().reindex(STATUS_GROUP_ORDER).fillna(0).astype(int).to_dict()
    by_continent = pair_annot["continent"].value_counts().to_dict()
    by_age = pair_annot["age_group"].value_counts().to_dict()
    by_set_status = (
        pair_annot.groupby(["analysis_set", "status_group"], dropna=False)
        .size()
        .reset_index(name="n_pairs")
        .to_dict(orient="records")
    )

    summary = {
        "n_bioprojects": n_bp,
        "n_pairs": n_pairs,
        "n_samples": n_pairs * 2,
        "n_pairs_by_analysis_set": {k: int(by_set.get(k, 0)) for k in ANALYSIS_SET_ORDER},
        "n_bioprojects_by_analysis_set": {
            k: int((bp_summary["analysis_set"] == k).sum()) for k in ANALYSIS_SET_ORDER
        },
        "n_pairs_by_status_group": {k: int(by_status_pairs.get(k, 0)) for k in STATUS_GROUP_ORDER},
        "n_bioprojects_by_status_group": {k: int(by_status_bp.get(k, 0)) for k in STATUS_GROUP_ORDER},
        "n_pairs_by_continent": {str(k): int(v) for k, v in by_continent.items()},
        "n_pairs_by_age_group": {str(k): int(v) for k, v in by_age.items()},
        "n_pairs_by_analysis_set_status_group": by_set_status,
        "status_group_definition": {
            "HC": "BioProject contains only HC-like status_cleaned values",
            "HC+disease": "BioProject contains both HC-like and non-HC status_cleaned values",
            "disease_study": "BioProject contains only non-HC status_cleaned values",
        },
        "count_unit": "pair",
    }
    return pair_annot, bp_summary, summary


def plot_cohort_overview_composite(pair_annot: pd.DataFrame, bp_summary: pd.DataFrame, summary: dict, out: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.2))
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    # Panel A: summary counts
    ax_a.axis("off")
    lines = [
        "Oral–gut paired dataset",
        "",
        f"BioProjects: {summary['n_bioprojects']}",
        f"  Single-pair: {summary['n_bioprojects_by_analysis_set']['single_pair']}",
        f"  Time-series: {summary['n_bioprojects_by_analysis_set']['time_series']}",
        "",
        f"Pairs: {summary['n_pairs']:,}",
        f"  Single-pair: {summary['n_pairs_by_analysis_set']['single_pair']:,}",
        f"  Time-series: {summary['n_pairs_by_analysis_set']['time_series']:,}",
        "",
        f"Samples (Oral+Gut): {summary['n_samples']:,}",
        "",
        "Status groups (BioProjects)",
        f"  HC: {summary['n_bioprojects_by_status_group']['HC']}",
        f"  HC+disease: {summary['n_bioprojects_by_status_group']['HC+disease']}",
        f"  disease_study: {summary['n_bioprojects_by_status_group']['disease_study']}",
    ]
    ax_a.text(
        0.04,
        0.97,
        "\n".join(lines),
        transform=ax_a.transAxes,
        va="top",
        ha="left",
        fontsize=11,
        family="DejaVu Sans",
        linespacing=1.35,
    )
    ax_a.set_title("A. Dataset summary", loc="left", fontsize=12, fontweight="bold")

    # Panel B: analysis_set × status_group stacked bars (pair counts)
    ct = (
        pair_annot.groupby(["analysis_set", "status_group"], dropna=False)
        .size()
        .unstack(fill_value=0)
        .reindex(index=ANALYSIS_SET_ORDER, columns=STATUS_GROUP_ORDER)
        .fillna(0)
    )
    y = np.arange(len(ANALYSIS_SET_ORDER))
    left = np.zeros(len(ANALYSIS_SET_ORDER))
    for status in STATUS_GROUP_ORDER:
        vals = ct[status].astype(float).values
        ax_b.barh(
            y,
            vals,
            left=left,
            color=STATUS_GROUP_COLORS[status],
            edgecolor="white",
            linewidth=0.6,
            height=0.62,
            label=status,
        )
        for yi, v, l in zip(y, vals, left):
            if v >= max(40, 0.04 * ct.values.sum()):
                ax_b.text(l + v / 2, yi, f"{int(v)}", ha="center", va="center", fontsize=8, color="white")
        left += vals
    ax_b.set_yticks(y)
    ax_b.set_yticklabels([ANALYSIS_SET_LABELS[a] for a in ANALYSIS_SET_ORDER])
    ax_b.set_xlabel("Number of pairs")
    ax_b.set_title("B. Study design × health status", loc="left", fontsize=12, fontweight="bold")
    ax_b.legend(frameon=False, fontsize=8, title="status_group", loc="lower right")
    ax_b.set_xlim(0, left.max() * 1.05)

    # Panel C: continent × analysis_set
    cont_counts = (
        pair_annot.groupby(["continent", "analysis_set"], dropna=False)
        .size()
        .unstack(fill_value=0)
        .reindex(columns=ANALYSIS_SET_ORDER)
        .fillna(0)
    )
    cont_order = cont_counts.sum(axis=1).sort_values(ascending=False).index.tolist()
    cont_counts = cont_counts.reindex(cont_order)
    x = np.arange(len(cont_order))
    width = 0.38
    for i, aset in enumerate(ANALYSIS_SET_ORDER):
        ax_c.bar(
            x + (i - 0.5) * width,
            cont_counts[aset].values,
            width=width,
            color=ANALYSIS_SET_COLORS[aset],
            edgecolor="white",
            linewidth=0.5,
            label=ANALYSIS_SET_LABELS[aset],
        )
    ax_c.set_xticks(x)
    ax_c.set_xticklabels(cont_order, rotation=25, ha="right", fontsize=9)
    ax_c.set_ylabel("Number of pairs")
    ax_c.set_title("C. Continent × study design", loc="left", fontsize=12, fontweight="bold")
    ax_c.legend(frameon=False, fontsize=8)

    # Panel D: age × analysis_set
    age_counts = (
        pair_annot.groupby(["age_group", "analysis_set"], dropna=False)
        .size()
        .unstack(fill_value=0)
        .reindex(columns=ANALYSIS_SET_ORDER)
        .fillna(0)
    )
    age_order = [a for a in AGE_GROUP_ORDER if a in age_counts.index] + [
        a for a in age_counts.index if a not in AGE_GROUP_ORDER
    ]
    age_counts = age_counts.reindex(age_order).fillna(0)
    x = np.arange(len(age_order))
    for i, aset in enumerate(ANALYSIS_SET_ORDER):
        ax_d.bar(
            x + (i - 0.5) * width,
            age_counts[aset].values,
            width=width,
            color=ANALYSIS_SET_COLORS[aset],
            edgecolor="white",
            linewidth=0.5,
            label=ANALYSIS_SET_LABELS[aset],
        )
    ax_d.set_xticks(x)
    ax_d.set_xticklabels(age_order, rotation=30, ha="right", fontsize=8)
    ax_d.set_ylabel("Number of pairs")
    ax_d.set_title("D. Age layer × study design", loc="left", fontsize=12, fontweight="bold")
    ax_d.legend(frameon=False, fontsize=8)

    fig.suptitle(
        "Cohort overview of oral–gut paired studies (count unit = pair)",
        fontsize=13,
        y=1.01,
    )
    fig.tight_layout()
    savefig(fig, out / "cohort_overview_composite.png")


def _sankey_layer_layout(labels: list[str], totals: dict[str, float], y_max: float, gap_frac: float = 0.03):
    """Return bottom/top y for each node so stacked heights sum to y_max with gaps."""
    total = sum(totals[l] for l in labels)
    if total <= 0:
        return {l: (0.0, 0.0) for l in labels}
    gap = y_max * gap_frac
    usable = y_max - gap * max(0, len(labels) - 1)
    coords = {}
    y = 0.0
    for i, lab in enumerate(labels):
        h = usable * (totals[lab] / total)
        coords[lab] = (y, y + h)
        y += h + (gap if i < len(labels) - 1 else 0.0)
    return coords


def _draw_sankey_flow(ax, x0, y0_bot, y0_top, x1, y1_bot, y1_top, color, alpha=0.35):
    verts = [
        (x0, y0_bot),
        (x0 + (x1 - x0) * 0.45, y0_bot),
        (x0 + (x1 - x0) * 0.55, y1_bot),
        (x1, y1_bot),
        (x1, y1_top),
        (x0 + (x1 - x0) * 0.55, y1_top),
        (x0 + (x1 - x0) * 0.45, y0_top),
        (x0, y0_top),
        (x0, y0_bot),
    ]
    codes = [
        MplPath.MOVETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.LINETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CLOSEPOLY,
    ]
    path = MplPath(verts, codes)
    patch = PathPatch(path, facecolor=color, edgecolor="none", alpha=alpha, lw=0, zorder=1)
    ax.add_patch(patch)
    return patch


def plot_cohort_overview_sankey(pair_annot: pd.DataFrame, out: Path) -> None:
    layers = [
        ("analysis_set", ANALYSIS_SET_ORDER, {a: ANALYSIS_SET_COLORS[a] for a in ANALYSIS_SET_ORDER}),
        ("status_group", STATUS_GROUP_ORDER, STATUS_GROUP_COLORS),
        (
            "continent",
            pair_annot["continent"].value_counts().index.tolist(),
            None,
        ),
        (
            "age_group",
            [a for a in AGE_GROUP_ORDER if a in set(pair_annot["age_group"])]
            + [a for a in pair_annot["age_group"].unique() if a not in AGE_GROUP_ORDER],
            None,
        ),
    ]
    # continent color map
    continent_labels = layers[2][1]
    cmap_c = plt.cm.Set2(np.linspace(0, 1, max(len(continent_labels), 1)))
    continent_colors = {lab: cmap_c[i] for i, lab in enumerate(continent_labels)}
    layers[2] = (layers[2][0], continent_labels, continent_colors)

    age_labels = layers[3][1]
    cmap_a = plt.cm.YlGnBu(np.linspace(0.25, 0.9, max(len(age_labels), 1)))
    age_colors = {lab: cmap_a[i] for i, lab in enumerate(age_labels)}
    layers[3] = (layers[3][0], age_labels, age_colors)

    y_max = float(len(pair_annot))
    node_coords = []
    for col, labels, colors in layers:
        totals = pair_annot[col].value_counts().to_dict()
        bp_totals = pair_annot.groupby(col)["bioproject"].nunique().to_dict()
        labels = [l for l in labels if totals.get(l, 0) > 0]
        coords = _sankey_layer_layout(labels, totals, y_max=1.0)
        node_coords.append(
            {
                "col": col,
                "labels": labels,
                "colors": colors,
                "coords": coords,
                "totals": totals,
                "bp_totals": bp_totals,
            }
        )

    # Wider gutters so outside labels sit clear of ribbons.
    xs = [0.0, 1.55, 3.10, 4.65]
    fig, ax = plt.subplots(figsize=(18.5, 8.8))
    ax.set_xlim(-1.70, 6.35)
    ax.set_ylim(-0.05, 1.12)
    ax.axis("off")

    node_half = 0.065
    # flows between adjacent layers
    for i in range(len(layers) - 1):
        left = node_coords[i]
        right = node_coords[i + 1]
        left_col, right_col = left["col"], right["col"]
        flow = (
            pair_annot.groupby([left_col, right_col], dropna=False)
            .size()
            .reset_index(name="n")
        )
        left_cursor = {lab: left["coords"][lab][1] for lab in left["labels"]}
        right_cursor = {lab: right["coords"][lab][1] for lab in right["labels"]}
        left_rank = {lab: r for r, lab in enumerate(left["labels"])}
        right_rank = {lab: r for r, lab in enumerate(right["labels"])}
        flow["_lr"] = flow[left_col].map(left_rank)
        flow["_rr"] = flow[right_col].map(right_rank)
        flow = flow.sort_values(["_lr", "_rr"])
        scale = 1.0 / y_max
        for _, row in flow.iterrows():
            a, b, n = row[left_col], row[right_col], float(row["n"])
            if a not in left_cursor or b not in right_cursor or n <= 0:
                continue
            h = n * scale
            y0_top = left_cursor[a]
            y0_bot = y0_top - h
            y1_top = right_cursor[b]
            y1_bot = y1_top - h
            left_cursor[a] = y0_bot
            right_cursor[b] = y1_bot
            color = left["colors"].get(a, "#9ecae1")
            _draw_sankey_flow(
                ax,
                xs[i] + node_half,
                y0_bot,
                y0_top,
                xs[i + 1] - node_half,
                y1_bot,
                y1_top,
                color,
                alpha=0.38,
            )

    layer_titles = ["Study design", "Health status", "Continent", "Age layer"]
    n_layers = len(node_coords)
    label_jobs = []
    for i, node in enumerate(node_coords):
        for lab in node["labels"]:
            y0, y1 = node["coords"][lab]
            height = y1 - y0
            ax.add_patch(
                FancyBboxPatch(
                    (xs[i] - node_half, y0),
                    2 * node_half,
                    max(height, 1e-4),
                    boxstyle="round,pad=0.004,rounding_size=0.01",
                    facecolor=node["colors"].get(lab, "#bdbdbd"),
                    edgecolor="white",
                    linewidth=0.8,
                    mutation_aspect=0.5,
                    zorder=3,
                )
            )
            label = ANALYSIS_SET_LABELS.get(lab, lab)
            n_pairs = int(node["totals"].get(lab, 0))
            n_bp = int(node["bp_totals"].get(lab, 0))
            label_jobs.append((i, lab, label, n_pairs, n_bp, y0, y1, height))
        ax.text(
            xs[i],
            1.06,
            layer_titles[i],
            ha="center",
            va="bottom",
            fontsize=11,
            fontweight="bold",
            zorder=6,
        )

    # White masks: left of study-design column, right of age column.
    ax.add_patch(
        plt.Rectangle(
            (xs[0] - 1.55, -0.02),
            1.45,
            1.08,
            facecolor="white",
            edgecolor="none",
            zorder=8,
            clip_on=False,
        )
    )
    ax.add_patch(
        plt.Rectangle(
            (xs[-1] + node_half + 0.02, -0.02),
            1.55,
            1.08,
            facecolor="white",
            edgecolor="none",
            zorder=8,
            clip_on=False,
        )
    )

    # Labels outside colored bars. First column → left; all others → right.
    min_label_gap = 0.060
    label_dx = 0.11
    for i in range(n_layers):
        layer_items = [j for j in label_jobs if j[0] == i]
        centers = []
        for _i, _lab, label, n_pairs, n_bp, y0, y1, _h in layer_items:
            centers.append([(y0 + y1) / 2.0, label, n_pairs, n_bp])
        centers.sort(key=lambda t: -t[0])
        ys = [c[0] for c in centers]
        for k in range(1, len(ys)):
            ys[k] = min(ys[k], ys[k - 1] - min_label_gap)
        ys[-1] = max(ys[-1], 0.03)
        for k in range(len(ys) - 2, -1, -1):
            ys[k] = max(ys[k], ys[k + 1] + min_label_gap)
        ys[0] = min(ys[0], 1.02)
        for k in range(1, len(ys)):
            ys[k] = min(ys[k], ys[k - 1] - min_label_gap)

        for (cy0, label, n_pairs, n_bp), ty in zip(centers, ys):
            bp_unit = "BP" if n_bp == 1 else "BPs"
            text = f"{label}\n{n_pairs:,} pairs | {n_bp} {bp_unit}"
            if i == 0:
                tx, ha = xs[i] - node_half - label_dx, "right"
                x_bar_edge = xs[i] - node_half
            else:
                tx, ha = xs[i] + node_half + label_dx, "left"
                x_bar_edge = xs[i] + node_half
            ax.text(
                tx,
                ty,
                text,
                ha=ha,
                va="center",
                fontsize=7.0,
                color="#252525",
                fontweight="normal",
                zorder=10,
                clip_on=False,
                bbox=dict(
                    boxstyle="round,pad=0.16",
                    facecolor="white",
                    edgecolor="#d9d9d9",
                    linewidth=0.35,
                    alpha=0.98,
                ),
            )
            if abs(ty - cy0) > 0.012:
                ax.plot(
                    [x_bar_edge, tx],
                    [cy0, ty],
                    color="#bdbdbd",
                    linewidth=0.6,
                    zorder=9,
                    clip_on=False,
                )

    ax.set_title(
        "Cohort flow: study design → health status → continent → age\n"
        "(ribbon width = pairs; labels show pair and BioProject counts)",
        fontsize=12.5,
        pad=16,
    )
    fig.tight_layout()
    savefig(fig, out / "cohort_overview_sankey.png")


def plot_cohort_overview_by_bioproject(bp_summary: pd.DataFrame, out: Path) -> None:
    d = bp_summary.copy().sort_values(
        ["status_rank", "n_pairs", "bioproject"],
        ascending=[True, False, True],
    ).reset_index(drop=True)
    n = len(d)
    fig_h = max(8.5, 0.28 * n + 1.8)
    fig, ax = plt.subplots(figsize=(11.5, fig_h))
    y = np.arange(n)
    colors = [STATUS_GROUP_COLORS.get(s, "#bdbdbd") for s in d["status_group"]]
    hatches = ["///" if a == "time_series" else "" for a in d["analysis_set"]]
    bars = ax.barh(y, d["n_pairs"].values, color=colors, edgecolor="white", height=0.78)
    for bar, hatch in zip(bars, hatches):
        bar.set_hatch(hatch)
        if hatch:
            bar.set_edgecolor("#636363")
            bar.set_linewidth(0.6)
    labels = []
    for _, row in d.iterrows():
        short_status = row["status_values"]
        if len(short_status) > 42:
            short_status = short_status[:39] + "..."
        labels.append(f"{row['bioproject']}  ({short_status})")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlabel("Number of pairs")
    ax.set_title("BioProject sample sizes by health-status group\n(hatch = time_series)", fontsize=12)
    legend_handles = [
        Patch(facecolor=STATUS_GROUP_COLORS[s], edgecolor="white", label=s) for s in STATUS_GROUP_ORDER
    ] + [
        Patch(facecolor="#ffffff", edgecolor="#636363", hatch="///", label="time_series"),
        Patch(facecolor="#bdbdbd", edgecolor="white", label="single_pair"),
    ]
    ax.legend(handles=legend_handles, frameon=False, fontsize=8, loc="lower right")
    for yi, n_pairs in zip(y, d["n_pairs"].values):
        ax.text(n_pairs + max(d["n_pairs"]) * 0.01, yi, str(int(n_pairs)), va="center", fontsize=7, color="#525252")
    ax.set_xlim(0, d["n_pairs"].max() * 1.18)

    # status-group side labels
    for status in STATUS_GROUP_ORDER:
        idx = np.where(d["status_group"].eq(status).values)[0]
        if len(idx) == 0:
            continue
        mid = float(np.mean(idx))
        ax.text(
            1.01,
            mid,
            status,
            va="center",
            ha="left",
            fontsize=8,
            color=STATUS_GROUP_COLORS[status],
            fontweight="bold",
            transform=ax.get_yaxis_transform(),
            clip_on=False,
        )
    fig.tight_layout()
    savefig(fig, out / "cohort_overview_by_bioproject.png")


def make_cohort_overview_figures(root: Path) -> None:
    pair_path = root / "results" / "cross_comparison" / "01_merged_metrics" / "pair_metrics_all_cohorts.tsv"
    if not pair_path.is_file():
        raise FileNotFoundError(f"Missing merged pair metrics: {pair_path}")
    out = ensure_dir(root / "results" / "cross_comparison" / "00_cohort_overview")
    plots = ensure_dir(out / "plots")

    pair_df = pd.read_csv(pair_path, sep="\t")
    pair_annot, bp_summary, summary = build_cohort_overview_tables(pair_df)

    pair_out_cols = [
        c
        for c in [
            "bioproject",
            "pair_id",
            "subject_id",
            "analysis_set",
            "status_cleaned",
            "status_group",
            "status_values",
            "bp_status_label",
            "continent_group",
            "continent",
            "age_cleaned_group",
            "age_group",
            "country_cleaned",
            "sex_cleaned",
        ]
        if c in pair_annot.columns
    ]
    pair_annot[pair_out_cols].to_csv(out / "cohort_overview_by_pair.tsv", sep="\t", index=False)
    bp_summary.to_csv(out / "cohort_overview_by_bioproject.tsv", sep="\t", index=False)
    with open(out / "cohort_overview_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2)

    plot_cohort_overview_composite(pair_annot, bp_summary, summary, plots)
    plot_cohort_overview_sankey(pair_annot, plots)
    plot_cohort_overview_by_bioproject(bp_summary, plots)
    print(
        f"cohort overview: {summary['n_bioprojects']} BP, {summary['n_pairs']} pairs "
        f"-> {out}"
    )


def make_cohort_figures(root: Path, cohort: str) -> None:
    base = root / "results" / cohort
    alpha_out = ensure_dir(base / "02_alpha" / "plots")
    delta = load_delta(base / "02_alpha" / "delta_alpha.tsv")
    plot_alpha_pooled(delta, alpha_out)
    plot_alpha_by_bp(delta, alpha_out)

    st2_out = ensure_dir(base / "04_source_tracking" / "sourcetracker" / "plots")
    st2 = pd.read_csv(base / "04_source_tracking" / "sourcetracker" / "paired_oral_fraction.tsv", sep="\t")
    st2 = st2.loc[st2["task_status"].astype(str).eq("ok")].copy() if "task_status" in st2.columns else st2
    plot_st2_pooled(st2, st2_out)
    plot_st2_by_bp(st2, st2_out)

    disp_path = base / "03_beta" / "dispersion" / "dispersion_test_summary.tsv"
    if disp_path.is_file():
        disp_out = ensure_dir(base / "03_beta" / "dispersion" / "plots")
        plot_dispersion(pd.read_csv(disp_path, sep="\t"), disp_out)


def make_cross_figures(root: Path) -> None:
    conc_dir = root / "results" / "cross_comparison" / "02_concordance"
    plot_concordance(
        pd.read_csv(conc_dir / "metric_concordance.tsv", sep="\t"),
        pd.read_csv(conc_dir / "metric_descriptives.tsv", sep="\t"),
        ensure_dir(conc_dir / "plots"),
    )


def _short_disease_phrase(status_values) -> str:
    if pd.isna(status_values):
        return "disease"
    parts = [p.strip() for p in str(status_values).split(";") if p.strip()]
    disease = []
    for part in parts:
        if is_hc_status(part):
            continue
        disease.append(DISEASE_LABEL_SHORT.get(part, part))
    if not disease:
        return "HC"
    if len(disease) == 1:
        return disease[0]
    return f"{disease[0]} +{len(disease) - 1}"


def make_pcoa_status_figures(root: Path) -> None:
    """PCoA of HC-only vs disease-only BioProjects, plus site-R² comparison.

    Mixed HC+disease BioProjects are excluded. Site R² is the existing
    per-BioProject paired adonis2 (D ~ site).
    """
    perm_path = root / "results" / "single_pair" / "03_beta" / "pcoa_permanova" / "pcoa_permanova_by_bp.tsv"
    coord_path = root / "results" / "single_pair" / "03_beta" / "pcoa_permanova" / "pcoa_coordinates_per_bp.tsv"
    overview_path = root / "results" / "cross_comparison" / "00_cohort_overview" / "cohort_overview_by_bioproject.tsv"
    for path in (perm_path, coord_path, overview_path):
        if not path.is_file():
            raise FileNotFoundError(f"Missing {path}")

    overview = pd.read_csv(overview_path, sep="\t")
    overview = overview.loc[overview["analysis_set"].astype(str).eq("single_pair")].copy()
    keep = overview.loc[overview["status_group"].isin(["HC", "disease_study"])].copy()
    keep["panel_label"] = keep.apply(
        lambda r: (
            str(r["bioproject"])
            if r["status_group"] == "HC"
            else f"{r['bioproject']}\n({_short_disease_phrase(r['status_values'])})"
        ),
        axis=1,
    )

    perm = pd.read_csv(perm_path, sep="\t")
    perm = perm.loc[
        perm["analysis_scope"].astype(str).eq("per_bioproject")
        & perm["metric"].isin(["weighted_unifrac", "unweighted_unifrac"])
    ].copy()
    perm = perm.merge(
        keep[["bioproject", "status_group", "status_values", "panel_label"]],
        on="bioproject",
        how="inner",
    )

    test_rows = []
    for metric, sub in perm.groupby("metric"):
        hc = sub.loc[sub["status_group"].eq("HC"), "r2_term"].astype(float)
        dz = sub.loc[sub["status_group"].eq("disease_study"), "r2_term"].astype(float)
        try:
            p = float(stats.mannwhitneyu(dz, hc, alternative="two-sided").pvalue)
        except ValueError:
            p = float("nan")
        test_rows.append(
            {
                "metric": metric,
                "n_hc_bp": int(hc.shape[0]),
                "n_disease_bp": int(dz.shape[0]),
                "median_r2_hc": float(hc.median()),
                "median_r2_disease": float(dz.median()),
                "mwu_p": p,
            }
        )
    tests = pd.DataFrame(test_rows)
    out_dir = root / "results" / "single_pair" / "03_beta" / "pcoa_permanova"
    plots = ensure_dir(out_dir / "plots")
    detail_cols = [
        "bioproject",
        "status_group",
        "status_values",
        "metric",
        "r2_term",
        "pseudo_F",
        "p_value",
        "pc1_percent",
        "pc2_percent",
        "n_samples",
    ]
    perm[detail_cols].sort_values(["metric", "status_group", "r2_term"]).to_csv(
        out_dir / "pcoa_site_r2_hc_vs_disease.tsv", sep="\t", index=False
    )

    metric = "weighted_unifrac"
    stats_m = perm.loc[perm["metric"].eq(metric)].copy()
    coords = pd.read_csv(
        coord_path,
        sep="\t",
        usecols=["PC1", "PC2", "site_for_analysis", "bioproject", "metric", "analysis_scope"],
    )
    coords = coords.loc[
        coords["analysis_scope"].astype(str).eq("per_bioproject")
        & coords["metric"].astype(str).eq(metric)
    ].copy()

    hc_order = (
        stats_m.loc[stats_m["status_group"].eq("HC")]
        .sort_values("r2_term")["bioproject"]
        .tolist()
    )
    dz_order = (
        stats_m.loc[stats_m["status_group"].eq("disease_study")]
        .sort_values("r2_term")["bioproject"]
        .tolist()
    )
    site_colors = {"Oral": "#D65F4B", "Gut": "#63B6C7"}
    label_map = keep.set_index("bioproject")["panel_label"].to_dict()
    stat_map = stats_m.set_index("bioproject")

    n_hc_cols = 4
    n_rows = max(int(np.ceil(len(hc_order) / n_hc_cols)), len(dz_order))
    fig = plt.figure(figsize=(18.6, n_rows * 2.55 + 3.4))
    gs = fig.add_gridspec(
        n_rows + 1,
        n_hc_cols + 1,
        height_ratios=[1] * n_rows + [1.15],
        width_ratios=[1, 1, 1, 1, 1.05],
        hspace=0.62,
        wspace=0.38,
    )

    def _draw_pcoa(ax, bp: str) -> None:
        sub = coords.loc[coords["bioproject"].astype(str).eq(str(bp))]
        for site in ("Oral", "Gut"):
            s = sub.loc[sub["site_for_analysis"].astype(str).eq(site)]
            if s.empty:
                continue
            ax.scatter(
                s["PC1"],
                s["PC2"],
                s=9,
                alpha=0.65,
                c=site_colors[site],
                linewidths=0,
                zorder=3,
            )
        ax.axhline(0, color="#e0e0e0", linewidth=0.6)
        ax.axvline(0, color="#e0e0e0", linewidth=0.6)
        row = stat_map.loc[bp]
        r2 = float(row["r2_term"])
        p = float(row["p_value"])
        p_txt = "P<0.001" if p <= 0.001 else f"P={p:.3f}"
        ax.set_title(f"{label_map.get(bp, bp)}\nR²={r2:.3f}, {p_txt}", fontsize=6.4, pad=2)
        ax.tick_params(labelsize=5)
        ax.set_xlabel(f"PC1 ({float(row['pc1_percent']):.0f}%)", fontsize=6)
        ax.set_ylabel(f"PC2 ({float(row['pc2_percent']):.0f}%)", fontsize=6)

    for i, bp in enumerate(hc_order):
        ax = fig.add_subplot(gs[i // n_hc_cols, i % n_hc_cols])
        _draw_pcoa(ax, bp)
    for i in range(len(hc_order), n_rows * n_hc_cols):
        ax = fig.add_subplot(gs[i // n_hc_cols, i % n_hc_cols])
        ax.axis("off")
    for i, bp in enumerate(dz_order):
        ax = fig.add_subplot(gs[i, n_hc_cols])
        _draw_pcoa(ax, bp)
    for i in range(len(dz_order), n_rows):
        ax = fig.add_subplot(gs[i, n_hc_cols])
        ax.axis("off")

    fig.text(0.34, 0.985, "HC-only BioProjects", ha="center", va="top", fontsize=11, color="#3182bd")
    fig.text(0.86, 0.985, "Disease-only BioProjects", ha="center", va="top", fontsize=11, color="#de2d26")

    bottom = gs[n_rows, :].subgridspec(1, 2, wspace=0.28)
    metric_titles = {
        "weighted_unifrac": "Weighted UniFrac",
        "unweighted_unifrac": "Unweighted UniFrac",
    }
    rng = np.random.default_rng(2026)
    for col, metric_name in enumerate(["weighted_unifrac", "unweighted_unifrac"]):
        ax = fig.add_subplot(bottom[0, col])
        sub = perm.loc[perm["metric"].eq(metric_name)].copy()
        data = []
        for i, group in enumerate(["HC", "disease_study"]):
            vals = sub.loc[sub["status_group"].eq(group), "r2_term"].astype(float).to_numpy()
            data.append(vals)
            jitter = rng.uniform(-0.08, 0.08, len(vals))
            color = "#3182bd" if group == "HC" else "#de2d26"
            ax.scatter(np.full(len(vals), i + 1) + jitter, vals, s=28, color=color, alpha=0.85, zorder=3)
            if group == "disease_study":
                labels = sub.loc[sub["status_group"].eq(group), ["bioproject", "r2_term", "status_values"]]
                for _, row in labels.iterrows():
                    ax.annotate(
                        _short_disease_phrase(row["status_values"]),
                        (2, float(row["r2_term"])),
                        textcoords="offset points",
                        xytext=(8, 0),
                        fontsize=6.5,
                        color="#de2d26",
                        va="center",
                    )
        ax.boxplot(
            data,
            positions=[1, 2],
            widths=0.45,
            showfliers=False,
            patch_artist=True,
            boxprops=dict(facecolor="none", edgecolor="#444"),
            medianprops=dict(color="#222", linewidth=1.4),
        )
        info = tests.loc[tests["metric"].eq(metric_name)].iloc[0]
        ax.set_xticks([1, 2])
        ax.set_xticklabels([f"HC-only\n({int(info['n_hc_bp'])} BP)", f"Disease-only\n({int(info['n_disease_bp'])} BP)"])
        ax.set_ylabel("PERMANOVA site R²")
        ax.set_ylim(0, 1)
        ax.set_title(
            f"{metric_titles[metric_name]}\n"
            f"median {info['median_r2_hc']:.2f} vs {info['median_r2_disease']:.2f}, "
            f"MWU p={info['mwu_p']:.3g}",
            fontsize=9,
        )

    handles = [
        plt.Line2D([], [], marker="o", linestyle="None", color=site_colors["Oral"], label="Oral"),
        plt.Line2D([], [], marker="o", linestyle="None", color=site_colors["Gut"], label="Gut"),
    ]
    fig.legend(handles=handles, loc="upper right", frameon=False, fontsize=8)
    fig.suptitle(
        "Single-pair oral–gut structure: HC-only vs disease-only BioProjects\n"
        "(mixed HC+disease studies excluded; PCoA = weighted UniFrac)",
        fontsize=12,
        y=1.02,
    )
    savefig(fig, plots / "pcoa_hc_vs_disease_weighted_unifrac.png")
    plot_mixed_bp_paired_distance(root, plots)
    print(
        f"pcoa status figure: HC {len(hc_order)} BP, disease {len(dz_order)} BP -> {plots}"
    )


def plot_mixed_bp_paired_distance(root: Path, out: Path) -> None:
    """Within mixed BioProjects, boxplot of paired oral–gut UniFrac by HC vs disease."""
    dist_path = root / "results" / "single_pair" / "03_beta" / "distance_matrix" / "paired_beta_distances_all_pairs.tsv"
    meta_path = root / "results" / "cross_comparison" / "01_merged_metrics" / "pair_metrics_all_cohorts.tsv"
    dist = pd.read_csv(dist_path, sep="\t")
    dist = dist.loc[
        dist["metric"].astype(str).eq("weighted_unifrac")
        & dist["dataset_label"].astype(str).eq("filter")
        & dist["has_complete_oral_gut_pair"].astype(str).isin(["True", "true", "1"])
    ].copy()
    meta = pd.read_csv(meta_path, sep="\t")
    meta = meta.loc[
        meta["analysis_set"].astype(str).eq("single_pair"),
        ["pair_id", "bioproject", "status_cleaned"],
    ].drop_duplicates("pair_id")
    df = dist.merge(meta, on=["pair_id", "bioproject"], how="inner")
    df["pair_status"] = df["status_cleaned"].map(_pair_status_binary)
    df["distance"] = pd.to_numeric(df["distance"], errors="coerce")
    df = df.dropna(subset=["distance"])

    bp_group = (
        df.groupby("bioproject")["status_cleaned"]
        .apply(lambda s: classify_bp_status_group(s.tolist())[1])
    )
    mixed = bp_group[bp_group.eq("HC+disease")].index.tolist()
    df = df.loc[df["bioproject"].isin(mixed) & df["pair_status"].isin(["HC", "disease"])].copy()

    labels = {}
    for bp, sub in df.groupby("bioproject"):
        disease = _bp_disease_label(sub)
        labels[bp] = f"{bp}\n{disease}" if disease else str(bp)

    metric = "weighted_unifrac"
    test_rows = []
    sub = df.loc[df["metric"].eq(metric)]
    for bp in mixed:
        bp_sub = sub.loc[sub["bioproject"].eq(bp)]
        hc = bp_sub.loc[bp_sub["pair_status"].eq("HC"), "distance"].to_numpy(dtype=float)
        dz = bp_sub.loc[bp_sub["pair_status"].eq("disease"), "distance"].to_numpy(dtype=float)
        if len(hc) == 0 or len(dz) == 0:
            continue
        try:
            p = float(stats.mannwhitneyu(dz, hc, alternative="two-sided").pvalue)
        except ValueError:
            p = float("nan")
        test_rows.append(
            {
                "bioproject": bp,
                "disease_label": _bp_disease_label(bp_sub),
                "metric": metric,
                "n_hc": int(len(hc)),
                "n_disease": int(len(dz)),
                "median_hc": float(np.median(hc)),
                "median_disease": float(np.median(dz)),
                "mwu_p": p,
            }
        )
    tests = pd.DataFrame(test_rows)
    if not tests.empty:
        tests["q_BH"] = _bh_fdr(tests["mwu_p"].tolist())
        tests["sig_stars"] = tests["q_BH"].map(_sig_stars).fillna("")
    tests_path = (
        root
        / "results"
        / "single_pair"
        / "03_beta"
        / "pcoa_permanova"
        / "paired_unifrac_hc_vs_disease_within_bp.tsv"
    )
    tests.to_csv(tests_path, sep="\t", index=False)

    if tests.empty:
        bp_order = sorted(mixed)
    else:
        tests = tests.copy()
        tests["delta"] = tests["median_disease"] - tests["median_hc"]
        bp_order = tests.sort_values("delta")["bioproject"].tolist()

    fig, ax = plt.subplots(figsize=(max(11, len(bp_order) * 1.15), 4.8))
    rng = np.random.default_rng(2026)
    metric_tests = tests.set_index("bioproject") if not tests.empty else tests
    positions, data, colors, ns, tick_pos, tick_labels = [], [], [], [], [], []
    star_spans = []
    for i, bp in enumerate(bp_order):
        bp_sub = sub.loc[sub["bioproject"].eq(bp)]
        xs = []
        for st in ("HC", "disease"):
            vals = bp_sub.loc[bp_sub["pair_status"].eq(st), "distance"].to_numpy(dtype=float)
            x = i + (-0.18 if st == "HC" else 0.18)
            positions.append(x)
            data.append(vals)
            colors.append(PAIR_STATUS_COLORS[st])
            ns.append(len(vals))
            xs.append(x)
        tick_pos.append(i)
        tick_labels.append(labels.get(bp, bp))
        stars = ""
        if bp in getattr(metric_tests, "index", []):
            stars = str(metric_tests.loc[bp, "sig_stars"] or "")
        if stars and len(xs) == 2:
            star_spans.append((xs[0], xs[1], stars, np.concatenate(data[-2:])))

    bp_artist = ax.boxplot(
        data,
        positions=positions,
        widths=0.28,
        showfliers=False,
        patch_artist=True,
        medianprops=dict(color="#222", linewidth=1.2),
    )
    for patch, color in zip(bp_artist["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)
        patch.set_edgecolor("#333")
    for pos, vals, color in zip(positions, data, colors):
        if len(vals) == 0:
            continue
        ax.scatter(
            np.full(len(vals), pos) + rng.uniform(-0.05, 0.05, len(vals)),
            vals,
            s=8,
            alpha=0.28,
            color=color,
            linewidths=0,
            zorder=3,
        )
    ymax = max(float(np.nanmax(v)) for v in data if len(v))
    ymin = min(float(np.nanmin(v)) for v in data if len(v))
    span = ymax - ymin if ymax > ymin else 0.1
    ax.set_ylim(ymin - 0.04 * span, ymax + 0.22 * span)
    for x0, x1, stars, vals in star_spans:
        y_br = float(np.nanmax(vals)) + 0.06 * span
        ax.plot(
            [x0, x0, x1, x1],
            [y_br - 0.02 * span, y_br, y_br, y_br - 0.02 * span],
            color="#8b0000",
            linewidth=1.0,
            zorder=4,
        )
        ax.text(
            (x0 + x1) / 2,
            y_br + 0.01 * span,
            stars,
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold",
            color="#8b0000",
        )
    y_top = ax.get_ylim()[1]
    for pos, n in zip(positions, ns):
        ax.text(pos, y_top, f"n={n}", ha="center", va="bottom", fontsize=6, color="#444", clip_on=False)
    ax.set_ylabel("Weighted UniFrac\noral–gut distance")
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_labels, rotation=0, fontsize=7.5)
    ax.set_xlim(-0.6, len(bp_order) - 0.4)

    handles = [
        Patch(facecolor=PAIR_STATUS_COLORS["HC"], edgecolor="#333", label="HC"),
        Patch(facecolor=PAIR_STATUS_COLORS["disease"], edgecolor="#333", label="disease"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle(
        "Paired oral–gut weighted UniFrac within mixed BioProjects\n"
        "(Mann–Whitney within BP; BH across BioProjects; * q<0.05)",
        y=1.08,
        fontsize=12,
    )
    fig.tight_layout()
    savefig(fig, out / "paired_unifrac_hc_vs_disease_mixed_bp.png")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        default="/home/zhangcaiyi/project/oral_gut_pair",
        help="oral_gut_pair project root",
    )
    parser.add_argument(
        "--cohorts",
        default="single_pair,time_series",
        help="Comma-separated cohorts for alpha/ST2/dispersion figures",
    )
    parser.add_argument("--skip-cross", action="store_true")
    parser.add_argument(
        "--overview-only",
        action="store_true",
        help="Only build cross_comparison/00_cohort_overview tables and plots",
    )
    parser.add_argument(
        "--skip-overview",
        action="store_true",
        help="Skip cohort overview figures when running full report figures",
    )
    parser.add_argument(
        "--alpha-status-only",
        action="store_true",
        help="Only build single_pair Δ alpha disease-stratified figures and tests",
    )
    parser.add_argument(
        "--skip-alpha-status",
        action="store_true",
        help="Skip single_pair Δ alpha disease-stratified figures",
    )
    parser.add_argument(
        "--pcoa-status-only",
        action="store_true",
        help="Only build single_pair PCoA figures split by HC-only vs disease-only BioProjects",
    )
    args = parser.parse_args()
    root = Path(args.project_root)

    if args.overview_only:
        print("=== figures for cohort overview ===")
        make_cohort_overview_figures(root)
        return

    if args.alpha_status_only:
        print("=== figures for single_pair alpha status ===")
        make_alpha_status_figures(root)
        return

    if args.pcoa_status_only:
        print("=== figures for single_pair PCoA HC vs disease ===")
        make_pcoa_status_figures(root)
        return

    for cohort in [c.strip() for c in args.cohorts.split(",") if c.strip()]:
        print(f"=== figures for {cohort} ===")
        make_cohort_figures(root, cohort)
    if not args.skip_cross:
        print("=== figures for cross_comparison ===")
        make_cross_figures(root)
    if not args.skip_overview:
        print("=== figures for cohort overview ===")
        make_cohort_overview_figures(root)
    if not args.skip_alpha_status:
        print("=== figures for single_pair alpha status ===")
        make_alpha_status_figures(root)


if __name__ == "__main__":
    main()
