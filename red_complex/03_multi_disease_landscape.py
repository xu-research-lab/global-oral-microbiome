#!/usr/bin/env python3
"""
03_multi_disease_landscape.py (P1.1)
=====================================
Multi-disease red complex landscape analysis.

For each major disease (vs HC), estimates:
  1. Detection rate odds ratio (logistic regression)
  2. Abundance fold change (linear model on log10, among detected)
  3. All models control: Site + AgeGroup + Sex + Platform + SeqRegion

Outputs:
  - multi_disease_effects.csv : effect size table
  - multi_disease_forest_plot.png : forest plot
  - multi_disease_detection_rates.png : detection rate bar chart
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
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

EXCLUDE_PLATFORMS = ["454 GS FLX Titanium", "Illumina GAIIx"]
MIN_SAMPLES_PER_DISEASE = 50

# Diseases to include (must have >= MIN_SAMPLES_PER_DISEASE)
DISEASE_ORDER = [
    "HC",
    "Periodontal Disease",
    "Dental Caries",
    "Head and Neck Cancer",
    "OSCC",
    "Breast Cancer",
    "Colorectal Cancer",
    "HIV",
    "Lung Cancer",
    "COVID-19",
    "ICU",
    "Sjögren's Syndrome",
    "IBD",
    "Pancreatic Cancer",
    "Acute Myeloid Leukemia",
    "Aphthous",
    "RA",
    "Preterm",
    "Oral Mucositis",
    "Other",
]

plt.rcParams.update({"figure.dpi": 150, "font.size": 9})


def load_and_clean():
    df = pd.read_csv(DATA_PATH, sep="\t", low_memory=False)

    # Numeric columns
    for col in ["P_gingivalis", "T_denticola", "T_forsythia",
                "RedComplexTotal", "RedComplex_detected", "AllThree_detected",
                "P_gingivalis_detected", "T_denticola_detected",
                "T_forsythia_detected"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)

    # Exclude low-detection platforms
    df = df[~df["Platform"].isin(EXCLUDE_PLATFORMS)]

    # Drop NA in key covariates
    for col in ["Site", "AgeGroup", "Sex", "Platform", "SeqRegion", "DiseaseStatus"]:
        df = df[df[col].notna()]

    # Filter to major diseases
    diseases_present = [d for d in DISEASE_ORDER if d in df["DiseaseStatus"].values]
    df = df[df["DiseaseStatus"].isin(diseases_present)]

    # Set reference categories
    df["DiseaseStatus"] = pd.Categorical(df["DiseaseStatus"],
                                          categories=[d for d in DISEASE_ORDER
                                                     if d in diseases_present],
                                          ordered=True)
    df["AgeGroup"] = pd.Categorical(df["AgeGroup"],
                                     categories=["0-6 months", "6-36 months",
                                                 "3-12 years", "12-18 years",
                                                 "18-60 years", "＞60 years"],
                                     ordered=True)
    # Set HC as reference for disease
    df["DiseaseStatus"] = df["DiseaseStatus"].cat.reorder_categories(
        ["HC"] + [d for d in diseases_present if d != "HC"]
    )

    print(f"  Cleaned data: {len(df):,} samples")
    print(f"  Diseases: {df['DiseaseStatus'].nunique()}")
    print(f"  HC samples: {(df['DiseaseStatus']=='HC').sum():,}")

    return df


def run_disease_logistic(df, disease, outcome="RedComplex_detected"):
    """Logistic regression for disease vs HC, controlling for confounders."""
    sub = df[df["DiseaseStatus"].isin(["HC", disease])].copy()
    if len(sub) < MIN_SAMPLES_PER_DISEASE * 2:
        return None

    # Reference: HC
    sub["DiseaseStatus"] = sub["DiseaseStatus"].cat.remove_unused_categories()

    try:
        formula = (f"{outcome} ~ C(DiseaseStatus, Treatment('HC')) + "
                   f"C(Site) + C(AgeGroup) + C(Sex) + C(Platform) + C(SeqRegion)")
        model = smf.logit(formula, data=sub).fit(disp=0, maxiter=200)
        # Extract disease coefficient
        disease_param = [p for p in model.params.index if "DiseaseStatus" in p]
        if len(disease_param) == 0:
            return None
        param_name = disease_param[0]
        coef = model.params[param_name]
        ci = model.conf_int().loc[param_name]
        pval = model.pvalues[param_name]
        return {
            "disease": disease,
            "n_disease": (sub["DiseaseStatus"] == disease).sum(),
            "n_hc": (sub["DiseaseStatus"] == "HC").sum(),
            "outcome": outcome,
            "or_": np.exp(coef),
            "or_ci_low": np.exp(ci[0]),
            "or_ci_high": np.exp(ci[1]),
            "coef": coef,
            "pvalue": pval,
            "model": "logistic",
        }
    except Exception as e:
        print(f"    {disease}: logistic model failed: {e}")
        return None


def run_disease_lm(df, disease, outcome="RedComplexTotal"):
    """Linear model on log10 abundance among detected samples."""
    sub = df[(df["DiseaseStatus"].isin(["HC", disease])) &
             (df[outcome] > 0)].copy()
    if len(sub) < MIN_SAMPLES_PER_DISEASE:
        return None
    if (sub["DiseaseStatus"] == disease).sum() < 20:
        return None

    sub["DiseaseStatus"] = sub["DiseaseStatus"].cat.remove_unused_categories()
    sub["log10_rc"] = np.log10(sub[outcome] + 1e-6)

    try:
        formula = ("log10_rc ~ C(DiseaseStatus, Treatment('HC')) + "
                   "C(Site) + C(AgeGroup) + C(Sex) + C(Platform) + C(SeqRegion)")
        model = smf.ols(formula, data=sub).fit()
        disease_param = [p for p in model.params.index if "DiseaseStatus" in p]
        if len(disease_param) == 0:
            return None
        param_name = disease_param[0]
        coef = model.params[param_name]
        ci = model.conf_int().loc[param_name]
        pval = model.pvalues[param_name]
        return {
            "disease": disease,
            "n_detected_disease": (sub["DiseaseStatus"] == disease).sum(),
            "n_detected_hc": (sub["DiseaseStatus"] == "HC").sum(),
            "outcome": outcome,
            "log10_fc": coef,
            "fc": 10**coef,
            "fc_ci_low": 10**ci[0],
            "fc_ci_high": 10**ci[1],
            "pvalue": pval,
            "model": "lm_log10",
        }
    except Exception as e:
        print(f"    {disease}: lm model failed: {e}")
        return None


def compute_detection_rates(df, outcome="RedComplex_detected"):
    """Raw detection rates by disease."""
    rates = df.groupby("DiseaseStatus", observed=False).agg(
        n=("SampleID", "count"),
        detection_rate=(outcome, "mean"),
        all3_rate=("AllThree_detected", "mean"),
        pg_rate=("P_gingivalis_detected", "mean"),
        td_rate=("T_denticola_detected", "mean"),
        tf_rate=("T_forsythia_detected", "mean"),
        rc_mean=("RedComplexTotal", "mean"),
        rc_median=("RedComplexTotal", "median"),
    ).reset_index()
    return rates


def plot_forest(combined, filename):
    """Forest plot of adjusted OR (logistic) and FC (linear model)."""
    # Filter to significant results
    plot_df = combined[combined["pvalue"] < 0.05].copy()
    if len(plot_df) == 0:
        plot_df = combined.copy()

    # Create two columns: one for OR (detection), one for FC (abundance)
    or_df = plot_df[plot_df["model"] == "logistic"].copy()
    fc_df = plot_df[plot_df["model"] == "lm_log10"].copy()

    # Sort by OR
    or_df = or_df.sort_values("or_", ascending=True)
    disease_order = or_df["disease"].tolist()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, max(6, len(disease_order) * 0.4)),
                                    gridspec_kw={"width_ratios": [1, 1]})

    # Panel A: Detection Odds Ratio
    y_pos = range(len(disease_order))
    for i, disease in enumerate(disease_order):
        row = or_df[or_df["disease"] == disease].iloc[0]
        color = "#c0392b" if row["pvalue"] < 0.001 else ("#e67e22" if row["pvalue"] < 0.01 else "#95a5a6")
        ax1.errorbar(row["or_"], i,
                     xerr=[[row["or_"] - row["or_ci_low"]],
                           [row["or_ci_high"] - row["or_"]]],
                     fmt="o", color=color, capsize=3, markersize=6)
        # Label with N
        ax1.annotate(f"n={row['n_disease']:,}", (row["or_ci_high"] + 0.2, i),
                    fontsize=7, va="center")

    ax1.axvline(x=1, color="grey", linestyle="--", linewidth=0.8)
    ax1.set_yticks(list(y_pos))
    ax1.set_yticklabels(disease_order)
    ax1.set_xscale("log")
    ax1.set_xlabel("Adjusted Odds Ratio (vs HC)")
    ax1.set_title("Detection Odds Ratio\n(logistic regression)")
    ax1.grid(axis="x", alpha=0.3)

    # Panel B: Abundance Fold Change
    fc_df_ordered = fc_df[fc_df["disease"].isin(disease_order)].copy()
    fc_df_ordered["disease"] = pd.Categorical(fc_df_ordered["disease"],
                                                categories=disease_order,
                                                ordered=True)
    fc_df_ordered = fc_df_ordered.sort_values("disease")

    for i, (_, row) in enumerate(fc_df_ordered.iterrows()):
        if row["fc"] <= 0:
            continue
        color = "#c0392b" if row["pvalue"] < 0.001 else ("#e67e22" if row["pvalue"] < 0.01 else "#95a5a6")
        ax2.errorbar(row["fc"], i,
                     xerr=[[row["fc"] - row["fc_ci_low"]],
                           [row["fc_ci_high"] - row["fc"]]],
                     fmt="s", color=color, capsize=3, markersize=6)

    ax2.axvline(x=1, color="grey", linestyle="--", linewidth=0.8)
    ax2.set_yticks(list(y_pos))
    ax2.set_yticklabels(disease_order)
    ax2.set_xscale("log")
    ax2.set_xlabel("Adjusted Fold Change (vs HC, | detected)")
    ax2.set_title("Abundance Fold Change\n(linear model, among detected)")
    ax2.grid(axis="x", alpha=0.3)

    fig.suptitle("Red Complex in Human Diseases: Adjusted Effect Size vs Healthy Control",
                 fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()
    fig.savefig(OUTDIR / filename, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Saved: {filename}")


def plot_detection_rates(rates_df, filename):
    """Bar chart of detection rates by disease."""
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))

    # Sort by total RC detection rate
    rates_df = rates_df.sort_values("detection_rate", ascending=True)
    diseases = rates_df["DiseaseStatus"].tolist()
    n_list = rates_df["n"].tolist()

    # Panel A: stacked detection bar
    x = np.arange(len(diseases))
    width = 0.7

    # Plot: none detected, partial, all three
    none_rate = 1 - rates_df["detection_rate"].values
    any1_rate = rates_df["detection_rate"].values - rates_df["all3_rate"].values
    all3_rate = rates_df["all3_rate"].values

    axes[0].barh(x, none_rate, height=width, color="#ecf0f1", label="None detected")
    axes[0].barh(x, any1_rate, height=width, left=none_rate, color="#f39c12", label="1-2 members")
    axes[0].barh(x, all3_rate, height=width, left=none_rate+any1_rate, color="#c0392b", label="All three")
    axes[0].set_yticks(x)
    axes[0].set_yticklabels([f"{d} (n={n:,})" for d, n in zip(diseases, n_list)], fontsize=8)
    axes[0].set_xlabel("Proportion of samples")
    axes[0].set_title("Red Complex Detection by Disease Status")
    axes[0].legend(loc="lower right")
    axes[0].set_xlim(0, 1)

    # Panel B: species-level detection rates
    sp_colors = {"P_gingivalis": "#e74c3c", "T_denticola": "#3498db",
                 "T_forsythia": "#2ecc71"}
    for sp, color in sp_colors.items():
        col = f"{sp.lower()[:2]}_rate" if sp == "P_gingivalis" else \
              f"{sp.lower()[:2]}_rate" if sp == "T_denticola" else \
              f"{sp.lower()[:2]}_rate" if sp == "T_forsythia" else None
        # Map to actual column names
        col_map = {"P_gingivalis": "pg_rate", "T_denticola": "td_rate",
                   "T_forsythia": "tf_rate"}
        if col_map[sp] in rates_df.columns:
            axes[1].barh(x + {"P_gingivalis": 0.2, "T_denticola": 0,
                              "T_forsythia": -0.2}[sp] * width/3,
                        rates_df[col_map[sp]].values,
                        height=width/3.5,
                        color=color, alpha=0.85, label=sp)

    axes[1].set_yticks(x)
    axes[1].set_yticklabels([f"{d} (n={n:,})" for d, n in zip(diseases, n_list)], fontsize=8)
    axes[1].set_xlabel("Detection rate")
    axes[1].set_title("Individual Species Detection Rate by Disease Status")
    axes[1].legend(loc="lower right")
    axes[1].set_xlim(0, 1)

    plt.tight_layout()
    fig.savefig(OUTDIR / filename, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {filename}")


def main():
    print("Loading data...")
    df = load_and_clean()

    # ============================================================
    # 1. Detection rates (raw)
    # ============================================================
    print("\n── Raw Detection Rates by Disease ──")
    rates = compute_detection_rates(df)
    rates.to_csv(OUTDIR / "multi_disease_detection_rates.csv", index=False)
    print(rates[["DiseaseStatus", "n", "detection_rate", "all3_rate",
                  "rc_mean"]].to_string(index=False))
    print(f"  Saved: multi_disease_detection_rates.csv")

    # ============================================================
    # 2. Adjusted effect sizes (logistic + linear model)
    # ============================================================
    print("\n── Adjusted Effect Sizes (disease vs HC) ──")
    diseases = [d for d in DISEASE_ORDER
                if d != "HC" and d in df["DiseaseStatus"].values
                and (df["DiseaseStatus"] == d).sum() >= MIN_SAMPLES_PER_DISEASE]

    logistic_results = []
    lm_results = []

    for disease in diseases:
        n_disease = (df["DiseaseStatus"] == disease).sum()
        print(f"\n  {disease} (n={n_disease:,}):")

        # Logistic model
        lr = run_disease_logistic(df, disease, "RedComplex_detected")
        if lr:
            print(f"    Detection OR: {lr['or_']:.2f} [{lr['or_ci_low']:.2f}-{lr['or_ci_high']:.2f}], p={lr['pvalue']:.2e}")
            logistic_results.append(lr)

        # Linear model
        lm = run_disease_lm(df, disease, "RedComplexTotal")
        if lm:
            print(f"    Abundance FC: {lm['fc']:.2f} [{lm['fc_ci_low']:.2f}-{lm['fc_ci_high']:.2f}], p={lm['pvalue']:.2e}")
            lm_results.append(lm)

    # Combine
    combined = pd.DataFrame(logistic_results + lm_results)
    combined.to_csv(OUTDIR / "multi_disease_effects.csv", index=False)
    print(f"\n  Saved: multi_disease_effects.csv ({len(combined)} rows)")

    # ============================================================
    # 3. Visualization
    # ============================================================
    print("\n── Generating plots ──")
    plot_detection_rates(rates, "multi_disease_detection_rates.png")
    plot_forest(combined, "multi_disease_forest_plot.png")

    # ============================================================
    # 4. Summary statistics
    # ============================================================
    print("\n" + "=" * 60)
    print("TOP FINDINGS")
    print("=" * 60)

    # Highest detection odds ratio
    if logistic_results:
        top_or = pd.DataFrame(logistic_results).nsmallest(5, "pvalue")
        print("\n  Top 5 diseases by detection OR significance:")
        for _, row in top_or.iterrows():
            print(f"    {row['disease']}: OR={row['or_']:.2f} [{row['or_ci_low']:.2f}-{row['or_ci_high']:.2f}], p={row['pvalue']:.2e}")

    # Highest abundance fold change
    if lm_results:
        top_fc = pd.DataFrame(lm_results).nsmallest(5, "pvalue")
        print("\n  Top 5 diseases by abundance FC significance:")
        for _, row in top_fc.iterrows():
            print(f"    {row['disease']}: FC={row['fc']:.2f} [{row['fc_ci_low']:.2f}-{row['fc_ci_high']:.2f}], p={row['pvalue']:.2e}")

    # Diseases with LOW RC (lower than HC)
    if logistic_results:
        low_or = pd.DataFrame(logistic_results)
        low_or = low_or[low_or["or_"] < 1.0].sort_values("or_")
        if len(low_or) > 0:
            print("\n  Diseases with LOWER RC detection than HC:")
            for _, row in low_or.iterrows():
                print(f"    {row['disease']}: OR={row['or_']:.2f}, p={row['pvalue']:.2e}")

    print("\nDone.")


if __name__ == "__main__":
    main()
