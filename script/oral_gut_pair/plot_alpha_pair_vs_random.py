#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_alpha_pair_vs_random.py

Alpha-diversity analogue of the beta "paired vs random" analysis, using the
SAME plotting template and the SAME random-sampling PROCEDURE as the beta run
(s04_2_beta_diversity_pair_vs_random.py). The draws are fresh (not literally
the same gut samples as the beta run), but the sampling process is identical:

  For each focal Oral sample (within a BioProject, complete pair):
    paired  = MEAN over its matched gut sample(s) of |alpha_oral - alpha_gut|
    random  = one gut drawn WITH replacement per iteration, 100 iterations,
              from the same BioProject + same stratum
              (stratify-cols = status_cleaned, the beta run's default),
              EXCLUDING the focal's own pair_id and its subject_id
              (longitudinal-safe); the 100 per-iteration values are
              aggregated by MEDIAN ("median of 100 gut")
  Test: two-sided Wilcoxon signed-rank on the matched values (identical to
  the beta pooled test), BH-corrected across the 3 alpha metrics.

  One shared rng (seed 2026, the beta run's base seed) drives all draws, so
  the three alpha metrics use the SAME random gut draws and are directly
  comparable to each other.

Metadata helpers (standardize_site / normalize_status_value /
prepare_metadata) are imported from the original s04_2 script so the
stratification and site standardization are byte-identical to beta.

Inputs:
  data/derived/single_pair/filter/filter_metadata.csv
  results/single_pair/02_alpha/alpha_diversity_all_samples.tsv
      (dataset_label=filter; per-sample Observed / Shannon / InvSimpson)

Outputs: results/new/figures/06_alpha_pair_vs_random/
  alpha_paired_vs_random_filter_observed_otus.{png,pdf}
  alpha_paired_vs_random_filter_shannon.{png,pdf}
  alpha_paired_vs_random_filter_invsimpson.{png,pdf}
  alpha_paired_vs_random_long_plot_table.tsv
  alpha_paired_vs_random_stats.tsv
"""

import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, wilcoxon

ROOT = Path("/home/hexincheng/projects/hxc/oral_project/oral_gut_pair")
S04_2_PATH = Path(
    "/home/hexincheng/projects/hxc/oral_project/oral_gut_pair_TimeSeries_11bp"
) / "scripts/s04_2_beta_diversity_pair_vs_random.py"
IN_META = ROOT / "data/derived/single_pair/filter/filter_metadata.csv"
IN_ALPHA = ROOT / "results/single_pair/02_alpha/alpha_diversity_all_samples.tsv"
OUT_DIR = ROOT / "results/new/figures/06_alpha_pair_vs_random"

N_RANDOM = 100
SEED = 2026
FOCAL_SITE = "Oral"
RANDOM_SITE = "Gut"
EXCLUDE_SAME_SUBJECT = True

METRICS = [
    {"col": "Observed", "suffix": "observed_otus", "title": "Observed OTUs"},
    {"col": "Shannon", "suffix": "shannon", "title": "Shannon index"},
    {"col": "InvSimpson", "suffix": "invsimpson", "title": "Inverse Simpson index"},
]
PAIRED_COLOR = "#D65F4B"
RANDOM_COLOR = "#63B6C7"
X_LABELS = ["Paired\n(same pair_id)", "Random-unpaired\n(median of 100 gut)"]

plt.rcParams.update(
    {
        "pdf.fonttype": 42,
        "font.size": 9,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
    }
)


def load_s04_2():
    spec = importlib.util.spec_from_file_location("s04_2_beta", S04_2_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_alpha_index() -> pd.DataFrame:
    a = pd.read_csv(
        IN_ALPHA, sep="\t",
        usecols=["dataset_label", "sample_id", "Observed", "Shannon", "InvSimpson"],
    )
    a = a.loc[a["dataset_label"] == "filter"].drop_duplicates("sample_id")
    for c in ("Observed", "Shannon", "InvSimpson"):
        a[c] = pd.to_numeric(a[c], errors="coerce")
    # zero-read samples have NaN alpha and were also absent from the beta
    # distance matrices (e.g. DRR182827, PRJDB8487, sample_total=0) — dropping
    # them here keeps the focal-unit set identical to the beta run
    n_bad = int(a[["Observed", "Shannon", "InvSimpson"]].isna().any(axis=1).sum())
    if n_bad:
        print(f"dropping {n_bad} samples with missing alpha values")
    a = a.loc[a[["Observed", "Shannon", "InvSimpson"]].notna().all(axis=1)]
    return a.set_index("sample_id")


def build_long_table(s04_2, alpha_idx: pd.DataFrame, rng) -> pd.DataFrame:
    md = s04_2.prepare_metadata(
        pd.read_csv(IN_META),
        sample_col="sample_id",
        bioproject_col="bioproject",
        site_col="site_group_cleaned",
        pair_col="pair_id",
        subject_col="subject_id",
        status_col="status_cleaned",
        disease_col="disease_study_cleaned",
        stratify_cols=["status_cleaned"],  # the beta run used this default
    )

    alpha_cols = [m["col"] for m in METRICS]
    have_alpha = set(alpha_idx.index)

    rows = []
    for bp, sub in md.groupby("bioproject_for_analysis", dropna=False):
        bp = str(bp)
        oral_df = sub.loc[sub["site_for_analysis"].eq(FOCAL_SITE), :]
        gut_df = sub.loc[sub["site_for_analysis"].eq(RANDOM_SITE), :].copy()
        if oral_df.empty or gut_df.empty:
            continue

        random_by_pair = {
            str(pid): g["sample_id_for_analysis"].astype(str).tolist()
            for pid, g in gut_df.groupby("pair_id_for_analysis", dropna=False)
        }

        # candidate pools per (bp, stratum) with alpha values available
        gut_df["has_alpha"] = gut_df["sample_id_for_analysis"].isin(have_alpha)
        gut_df = gut_df.loc[gut_df["has_alpha"]].copy()
        gut_pos = {s: i for i, s in enumerate(gut_df["sample_id_for_analysis"])}
        gut_vals = np.column_stack(
            [alpha_idx.loc[gut_df["sample_id_for_analysis"], c].to_numpy(float)
             for c in alpha_cols]
        )
        gut_pair = gut_df["pair_id_for_analysis"].astype(str).to_numpy()
        gut_subj = gut_df["subject_id_for_analysis"].astype(str).to_numpy()
        gut_stratum = gut_df["stratum_key_for_random"].astype(str).to_numpy()

        for _, focal in oral_df.iterrows():
            focal_sample = str(focal["sample_id_for_analysis"])
            if focal_sample not in have_alpha:
                continue
            pair_id = str(focal["pair_id_for_analysis"])
            stratum = str(focal["stratum_key_for_random"])
            focal_subject = str(focal.get("subject_id_for_analysis", "NA"))

            matched_guts = [s for s in random_by_pair.get(pair_id, [])
                            if s in have_alpha]
            if not matched_guts:
                continue  # incomplete pair, same criterion as beta

            focal_vals = alpha_idx.loc[focal_sample, alpha_cols].to_numpy(float)
            # (n_matched_guts, n_metrics) — same layout as gut_vals below
            paired_vals = np.array(
                [alpha_idx.loc[s, alpha_cols].to_numpy(float) for s in matched_guts]
            )
            paired = np.nanmean(np.abs(focal_vals - paired_vals), axis=0)

            pool = np.where(
                (gut_stratum == stratum)
                & (gut_pair != pair_id)
                & ((gut_subj != focal_subject)
                   | (focal_subject in {"", "NA", "nan", "None"}))
            )[0]
            pool_size = int(pool.size)
            if pool_size == 0:
                continue

            picks = rng.integers(0, pool_size, size=N_RANDOM)
            random = np.nanmedian(np.abs(focal_vals - gut_vals[pool[picks]]), axis=0)

            for j, m in enumerate(METRICS):
                rows.append({
                    "bioproject": bp,
                    "pair_id": pair_id,
                    "focal_sample_id": focal_sample,
                    "metric": m["suffix"],
                    "distance_type": "Paired",
                    "value": float(paired[j]),
                    "stratum_key": stratum,
                    "candidate_pool_size": pool_size,
                })
                rows.append({
                    "bioproject": bp,
                    "pair_id": pair_id,
                    "focal_sample_id": focal_sample,
                    "metric": m["suffix"],
                    "distance_type": "Random-unpaired",
                    "value": float(random[j]),
                    "stratum_key": stratum,
                    "candidate_pool_size": pool_size,
                })
    return pd.DataFrame(rows)


def p_to_stars(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def format_p(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    return f"{p:.2e}" if p < 1e-4 else f"{p:.4f}"


def draw_violin(ax, x: float, vals: np.ndarray, color: str, width: float = 0.6) -> None:
    """Half-free violin (KDE clamped to the observed data range) + slim inner box."""
    kde = gaussian_kde(vals)
    y = np.linspace(vals.min(), vals.max(), 200)
    d = kde(y)
    d = d / d.max() * (width / 2.0)
    ax.fill_betweenx(y, x - d, x + d, facecolor=color, edgecolor="black",
                     linewidth=0.8, alpha=0.85, zorder=2)
    ax.boxplot(
        [vals], positions=[x], widths=0.12, patch_artist=True, showfliers=False,
        boxprops={"facecolor": "white", "edgecolor": "black", "linewidth": 0.8},
        medianprops={"color": "black", "linewidth": 1.4},
        whiskerprops={"linewidth": 0.8, "color": "black"},
        capprops={"linewidth": 0.8, "color": "black"},
        zorder=3,
    )


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.clip(ranked, 0, 1)
    return q


def draw_panel(d: pd.DataFrame, metric: dict, stats_row: pd.Series) -> None:
    paired = d.loc[d["distance_type"] == "Paired", "value"].to_numpy(float)
    random = d.loc[d["distance_type"] == "Random-unpaired", "value"].to_numpy(float)

    p_text = format_p(float(stats_row["p_value"]))
    stars = str(stats_row["p_significance_label"])
    n_units = int(stats_row["n_focal_units"])
    median_diff = float(stats_row["median_difference_paired_minus_random"])

    # canvas shrunk to near-minimum (user request: as small as possible):
    # 3.0x3.0 in with fonts dropped ~1 pt; xlim tightened to (-0.42, 1.08)
    # (violins keep 0.12 data-unit side margins) so the two long x labels
    # still clear each other — with the old (-0.55, 1.3) they collided
    # below 3.6 in width
    fig, ax = plt.subplots(figsize=(3.0, 3.4))

    # --- violins (KDE clamped to the data range) with slim inner boxes ---
    draw_violin(ax, 0, paired, PAIRED_COLOR)
    draw_violin(ax, 0.75, random, RANDOM_COLOR)

    ax.set_xlim(-0.42, 1.17)

    y_min = min(paired.min(), random.min())
    y_max = max(paired.max(), random.max())
    span = y_max - y_min
    # bracket hugging the violins (user request 2026-10-07: the old
    # +0.045·span bracket under a +0.11·span head margin floated too high)
    ax.set_ylim(y_min - 0.03 * span, y_max + 0.12 * span)
    y_br = y_max + 0.028 * span
    tick = 0.016 * span
    # bracket ends back ON the violin centre lines (user request 2026-10-07:
    # 括号要在中线 — the ±0.135 widening put the ticks off the whisker lines)
    x_l = 0
    x_r = 0.75
    ax.plot([x_l, x_l, x_r, x_r], [y_br - tick, y_br, y_br, y_br - tick],
            lw=1.0, color="black", zorder=4)
    ax.text((x_l + x_r) / 2, y_br + 0.008 * span, stars, ha="center",
            va="bottom", fontsize=9.5, color="black")

    # stats line removed (user request); bracket + stars remain
    ax.set_title(metric["title"], fontsize=10, pad=6)

    ax.set_xticks([0, 0.75])
    ax.set_xticklabels(X_LABELS, fontsize=8)
    # two lines so the rotated label fits the 3.0-in canvas height
    ax.set_ylabel(f"Absolute oral−gut difference,\n{metric['title']}", fontsize=9)
    ax.grid(axis="y", linestyle=":", linewidth=0.5, alpha=0.4)
    ax.set_axisbelow(True)

    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT_DIR / f"alpha_paired_vs_random_filter_{metric['suffix']}.{ext}",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    s04_2 = load_s04_2()
    alpha_idx = load_alpha_index()
    rng = np.random.default_rng(SEED)
    long_df = build_long_table(s04_2, alpha_idx, rng)
    long_df.to_csv(OUT_DIR / "alpha_paired_vs_random_long_plot_table.tsv", sep="\t",
                   index=False)

    stat_rows = []
    for m in METRICS:
        d = long_df.loc[long_df["metric"] == m["suffix"]]
        paired = d.loc[d["distance_type"] == "Paired", "value"].to_numpy(float)
        random = d.loc[d["distance_type"] == "Random-unpaired", "value"].to_numpy(float)
        res = wilcoxon(paired, random, alternative="two-sided",
                       zero_method="wilcox", mode="auto")
        stat_rows.append({
            "metric": m["suffix"],
            "metric_display": m["title"],
            "n_focal_units": len(paired),
            "median_paired": float(np.median(paired)),
            "median_random_unpaired": float(np.median(random)),
            "median_difference_paired_minus_random": float(np.median(paired - random)),
            "frac_paired_less_random": float((paired < random).mean()),
            "statistic": float(res.statistic),
            "p_value": float(res.pvalue),
        })
    stats = pd.DataFrame(stat_rows)
    stats["p_significance_label"] = stats["p_value"].map(p_to_stars)
    stats["q_value_BH_across_metrics"] = bh_fdr(stats["p_value"].to_numpy())
    stats.to_csv(OUT_DIR / "alpha_paired_vs_random_stats.tsv", sep="\t", index=False)

    for m in METRICS:
        d = long_df.loc[long_df["metric"] == m["suffix"]]
        row = stats.loc[stats["metric"] == m["suffix"]].iloc[0]
        draw_panel(d, m, row)
        print(f"done: {m['suffix']}  P={format_p(row['p_value'])} "
              f"({row['p_significance_label']})  q={row['q_value_BH_across_metrics']:.3g}  "
              f"n={row['n_focal_units']}  "
              f"dmedian={row['median_difference_paired_minus_random']:.4g}")
    print(f"outputs -> {OUT_DIR}")


if __name__ == "__main__":
    main()