#!/usr/bin/env python3
"""
01_technical_assessment.py
===========================
Evaluate how technical factors (platform, sequencing region, primer sets)
affect red complex detection sensitivity.

This is a PREREQUISITE for all downstream analyses — without understanding
technical biases, biological conclusions may be spurious.

Analyses:
  1. Detection rate matrix: platform × sequencing_region
  2. Detection rate matrix: platform × sequencing_region × species
  3. Relative abundance distribution by platform/region (HC only)
  4. Quantify variance explained by technical vs biological factors
  5. Identify "comparable" subsets for downstream sensitivity analyses

Output:
  - technical_detection_matrix.csv
  - technical_variance_decomposition.csv
  - platform_region_rc_detection_heatmap.png
  - rc_abundance_by_platform_region.png
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

OUTDIR = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/red_complex")
DATA_PATH = OUTDIR / "red_complex_sample_abundance.tsv"

# ── Style ──────────────────────────────────────────────────────────
plt.rcParams.update({
    "figure.dpi": 150,
    "font.size": 9,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
})


def load_data():
    df = pd.read_csv(DATA_PATH, sep="\t", low_memory=False)
    # Ensure numeric columns
    for col in ["P_gingivalis", "T_denticola", "T_forsythia", "RedComplexTotal"]:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)
    for col in ["P_gingivalis_detected", "T_denticola_detected",
                "T_forsythia_detected", "RedComplex_detected", "AllThree_detected"]:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)
    return df


def detection_matrix(df, group_cols, detect_cols):
    """Build detection rate matrix for given grouping columns."""
    records = []
    for group_keys, subset in df.groupby(group_cols):
        if not isinstance(group_keys, tuple):
            group_keys = (group_keys,)
        row = dict(zip(group_cols, group_keys))
        row["n_samples"] = len(subset)
        for dc in detect_cols:
            row[f"{dc}_rate"] = subset[dc].mean()
        records.append(row)
    return pd.DataFrame(records)


def plot_detection_heatmap(mat_df, x_col, y_col, value_col,
                           title, filename, annot=True):
    """Plot a detection rate heatmap."""
    pivot = mat_df.pivot(index=y_col, columns=x_col, values=value_col)
    # Sort rows/cols by overall rate
    pivot["_avg"] = pivot.mean(axis=1)
    pivot = pivot.sort_values("_avg", ascending=False).drop(columns="_avg")

    fig, ax = plt.subplots(figsize=(max(10, len(pivot.columns) * 1.2),
                                    max(4, len(pivot) * 0.5)))
    sns.heatmap(pivot, annot=annot, fmt=".2f", cmap="RdYlGn",
                vmin=0, vmax=1, linewidths=0.5, ax=ax,
                cbar_kws={"label": "Detection Rate"})
    ax.set_title(title)
    ax.set_ylabel(y_col.replace("_", " ").title())
    ax.set_xlabel(x_col.replace("_", " ").title())
    plt.tight_layout()
    fig.savefig(OUTDIR / filename, dpi=200)
    plt.close(fig)
    print(f"  Saved: {filename}")


def plot_abundance_by_tech(df, tech_col, filename):
    """Boxplot of RC abundance by technical factor (HC only)."""
    hc = df[df["DiseaseStatus"] == "HC"].copy()

    # Focus on top categories
    top_cats = hc[tech_col].value_counts().head(10).index
    plot_df = hc[hc[tech_col].isin(top_cats)]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    species = ["P_gingivalis", "T_denticola", "T_forsythia", "RedComplexTotal"]

    for ax, sp in zip(axes.flat, species):
        # Only show detected samples for clarity
        sub = plot_df[plot_df[sp] > 0]
        if len(sub) < 10:
            ax.text(0.5, 0.5, f"Insufficient data\n(n={len(sub)})",
                    transform=ax.transAxes, ha="center")
            continue
        order = sub.groupby(tech_col)[sp].median().sort_values(ascending=False).index
        sns.boxplot(data=sub, x=tech_col, y=sp, order=order, ax=ax,
                    showfliers=False, palette="Set2")
        ax.set_yscale("log")
        ax.set_ylabel("Relative Abundance (log scale)")
        ax.set_xlabel("")
        ax.set_title(sp)
        ax.tick_params(axis="x", rotation=45)

    fig.suptitle(f"Red Complex Abundance by {tech_col} (HC, detected only)",
                 fontsize=13, y=1.01)
    plt.tight_layout()
    fig.savefig(OUTDIR / filename, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {filename}")


def variance_decomposition(df):
    """Quantify variance in RC detection explained by technical vs biological factors.

    Uses ANOVA (Type II via linear model) on logit-transformed detection.
    """
    from statsmodels.formula.api import ols
    from statsmodels.stats.anova import anova_lm

    # Use HC only to isolate technical from disease effects
    hc = df[df["DiseaseStatus"] == "HC"].copy()

    # Focus on samples with known platform, region, site, geography
    sub = hc.dropna(subset=["Platform", "SeqRegion", "Site", "Geography",
                             "Sex", "AgeGroup"]).copy()
    # Exclude "other" / small categories
    for col in ["Platform", "SeqRegion", "Site", "Geography"]:
        valid = sub[col].value_counts()
        valid = valid[valid >= 50].index
        sub = sub[sub[col].isin(valid)]

    print(f"\n  Variance decomposition on {len(sub)} HC samples")

    # Model: detection ~ platform + seq_region + site + geography + sex + age_group
    formula = ("AllThree_detected ~ C(Platform) + C(SeqRegion) + C(Site) "
               "+ C(Geography) + C(Sex) + C(AgeGroup)")
    try:
        model = ols(formula, data=sub).fit()
        anova = anova_lm(model, typ=2)
        anova["pct_variance"] = 100 * anova["sum_sq"] / anova["sum_sq"].sum()
        anova["significant"] = anova["PR(>F)"].apply(lambda p: "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "ns")))
        print(anova[["sum_sq", "pct_variance", "PR(>F)", "significant"]].to_string())
        return anova
    except Exception as e:
        print(f"  Warning: ANOVA failed: {e}")
        return None


def comparable_subsets(df):
    """Identify subsets where technical factors are comparable across groups.

    Key question: which platform+region combinations have enough samples
    in both HC and disease groups for valid comparison?
    """
    print("\n── Comparable Subsets for HC vs Disease Comparison ──")

    # Define major platform × region combinations
    df["tech_combo"] = df["Platform"].astype(str) + " | " + df["SeqRegion"].astype(str)

    combos = df["tech_combo"].value_counts()
    major_combos = combos[combos >= 500].index

    records = []
    for combo in major_combos:
        sub = df[df["tech_combo"] == combo]
        diseases = sub["DiseaseStatus"].value_counts()
        hc_n = diseases.get("HC", 0)
        # Find disease with most samples besides HC
        disease_n = {k: v for k, v in diseases.items() if k != "HC"}
        if disease_n:
            top_disease = max(disease_n, key=disease_n.get)
            top_n = disease_n[top_disease]
        else:
            top_disease = "N/A"
            top_n = 0
        records.append({
            "tech_combo": combo,
            "total_n": len(sub),
            "HC_n": hc_n,
            "top_disease": top_disease,
            "disease_n": top_n,
            "RC_detection": sub["RedComplex_detected"].mean(),
        })

    combo_df = pd.DataFrame(records).sort_values("total_n", ascending=False)
    print(combo_df.to_string(index=False))
    return combo_df


def main():
    print("Loading data...")
    df = load_data()
    print(f"  Samples: {len(df):,}")

    detect_cols = ["RedComplex_detected", "AllThree_detected",
                   "P_gingivalis_detected", "T_denticola_detected",
                   "T_forsythia_detected"]

    # ── 1. Detection rate matrix: platform × sequencing_region ─────
    print("\n" + "=" * 60)
    print("1. DETECTION MATRIX: Platform × Sequencing Region")
    print("=" * 60)

    # Use HC only to isolate technical effects from disease
    hc = df[df["DiseaseStatus"] == "HC"].copy()

    mat = detection_matrix(hc, ["Platform", "SeqRegion"], detect_cols)
    mat = mat[mat["n_samples"] >= 20]  # filter tiny groups
    mat.to_csv(OUTDIR / "technical_detection_matrix.csv", index=False)
    print(f"  Saved: technical_detection_matrix.csv ({len(mat)} combinations)")

    # Pivot and display top rows
    for dc in ["AllThree_detected_rate", "RedComplex_detected_rate"]:
        print(f"\n  Top 10 by {dc}:")
        top10 = mat.nlargest(10, dc)[["Platform", "SeqRegion", "n_samples", dc]]
        print(top10.to_string(index=False))
        print(f"\n  Bottom 10 by {dc}:")
        bot10 = mat.nsmallest(10, dc)[["Platform", "SeqRegion", "n_samples", dc]]
        print(bot10.to_string(index=False))

    # ── 1b. Platform × Region heatmaps ───────────────────────────────
    print("\n── Heatmaps ──")
    heatmap_cols = [f"{dc}_rate" for dc in detect_cols]
    for dc, hc_col in zip(detect_cols, heatmap_cols):
        if hc_col not in mat.columns:
            print(f"  Skipping {hc_col} (not in matrix)")
            continue
        plot_detection_heatmap(
            mat, "SeqRegion", "Platform", hc_col,
            f"HC {dc.replace('_detected', '').replace('_', ' ').title()} Detection Rate\nby Platform and Sequencing Region",
            f"heatmap_platform_region_{dc}.png"
        )

    # ── 2. Abundance distributions ─────────────────────────────────
    print("\n" + "=" * 60)
    print("2. ABUNDANCE BY TECHNICAL FACTORS (HC only)")
    print("=" * 60)

    plot_abundance_by_tech(hc, "Platform", "boxplot_rc_by_platform_hc.png")
    plot_abundance_by_tech(hc, "SeqRegion", "boxplot_rc_by_seqregion_hc.png")

    # ── 3. Variance decomposition ──────────────────────────────────
    print("\n" + "=" * 60)
    print("3. VARIANCE DECOMPOSITION")
    print("=" * 60)

    anova_result = variance_decomposition(df)
    if anova_result is not None:
        anova_result.to_csv(OUTDIR / "technical_variance_decomposition.csv")
        print(f"  Saved: technical_variance_decomposition.csv")

    # ── 4. Comparable subsets ──────────────────────────────────────
    print("\n" + "=" * 60)
    print("4. COMPARABLE SUBSETS FOR DISEASE ANALYSIS")
    print("=" * 60)

    combo_df = comparable_subsets(df)
    combo_df.to_csv(OUTDIR / "comparable_subsets.csv", index=False)
    print(f"  Saved: comparable_subsets.csv")

    # ── 5. Summary for downstream analysis ─────────────────────────
    print("\n" + "=" * 60)
    print("5. KEY FINDINGS FOR DOWNSTREAM ANALYSIS")
    print("=" * 60)

    # Platform warning
    low_plat = hc.groupby("Platform")["RedComplex_detected"].mean()
    low_plat = low_plat[low_plat < 0.2]
    if len(low_plat) > 0:
        print(f"\n  ⚠️  Low-detection platforms (<20%):")
        for plat, rate in low_plat.items():
            n = (hc["Platform"] == plat).sum()
            print(f"     - {plat}: {rate:.1%} (n={n})")

    # Region warning
    low_reg = hc.groupby("SeqRegion")["RedComplex_detected"].mean()
    low_reg = low_reg[low_reg < 0.2]
    if len(low_reg) > 0:
        print(f"\n  ⚠️  Low-detection regions (<20%):")
        for reg, rate in low_reg.items():
            n = (hc["SeqRegion"] == reg).sum()
            print(f"     - {reg}: {rate:.1%} (n={n})")

    # Recommendation
    print(f"\n  ✅ Recommended primary analysis subset:")
    v4_miseq = len(hc[(hc["SeqRegion"] == "V4") & (hc["Platform"] == "Illumina MiSeq")])
    v3v4_miseq = len(hc[(hc["SeqRegion"] == "V3-V4") & (hc["Platform"] == "Illumina MiSeq")])
    print(f"     - V4 + MiSeq: {v4_miseq:,} HC samples (largest homogeneous subset)")
    print(f"     - V3-V4 + MiSeq: {v3v4_miseq:,} HC samples")
    print(f"\n  ✅ All models MUST include Platform + SeqRegion as covariates")
    print(f"  ✅ Sensitivity analysis: repeat key results in V4+MiSeq subset")
    print(f"  ⚠️  Exclude 454 GS FLX Titanium + GAIIx from primary analysis")
    print(f"     (detection rates too low, likely due to read length / reference DB)")

    print("\nDone.")


if __name__ == "__main__":
    main()
