#!/usr/bin/env python3
"""
Random Forest v2: Plot results.
Updated paths for v2 output.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn import metrics

BASE = Path("/beegfs/huangxiaochang/global_oral/new_analysis/randomforest_v2")
DISEASE_DIR = Path("/beegfs/huangxiaochang/global_oral/new_analysis/disease_prediction")
sns.set_theme(style="whitegrid")


def save_heatmap(df, index, columns, values, out_png, title, annot_fmt=".3f",
                 vmin=0, vmax=1, cmap="RdYlBu_r"):
    mat = df.pivot(index=index, columns=columns, values=values)
    plt.figure(figsize=(max(6, 1.3 * len(mat.columns)),
                         max(4, 0.5 * len(mat.index) + 1)))
    sns.heatmap(mat, annot=True, fmt=annot_fmt, cmap=cmap, vmin=vmin, vmax=vmax,
                linewidths=0.5, cbar_kws={"label": values})
    plt.title(title)
    plt.tight_layout()
    plt.savefig(BASE / out_png, dpi=300)
    plt.close()


def save_bar(df, x, y, hue, out_png, title, ylim=(0, 1)):
    plt.figure(figsize=(max(8, 0.7 * len(df[x].astype(str).unique()) + 2), 5))
    sns.barplot(data=df, x=x, y=y, hue=hue)
    if ylim is not None: plt.ylim(*ylim)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(BASE / out_png, dpi=300)
    plt.close()


def main():
    # Region OVA
    if (BASE / "region_ova_summary.csv").exists():
        region_summary = pd.read_csv(BASE / "region_ova_summary.csv")
        save_heatmap(region_summary, index="target_region", columns="site",
                     values="mean_auc", out_png="region_ova_heatmap.png",
                     title="Geography RF Mean AUC by Site and Target Region")
        save_bar(region_summary.sort_values(["site", "mean_auc"], ascending=[True, False]),
                 x="target_region", y="mean_auc", hue="site",
                 out_png="region_ova_barplot.png", title="Geography RF Mean AUC by Site")

        if (BASE / "region_ova_predictions.csv").exists():
            preds = pd.read_csv(BASE / "region_ova_predictions.csv")
            fold = preds.groupby(["site", "target_region", "fold"]).apply(
                lambda x: pd.Series({
                    "fold_auc": metrics.roc_auc_score(x["true_label"], x["pred_prob"])
                    if x["true_label"].nunique() > 1 else float("nan")
                })).reset_index().dropna()
            plt.figure(figsize=(11, 5))
            sns.boxplot(data=fold, x="target_region", y="fold_auc", hue="site", showfliers=False)
            plt.ylim(0, 1); plt.xticks(rotation=30, ha="right")
            plt.title("Geography RF Fold-level AUC Distribution")
            plt.tight_layout()
            plt.savefig(BASE / "region_ova_fold_auc_boxplot.png", dpi=300)
            plt.close()

    # Sex BP LOO
    if (BASE / "sex_bp_loo_summary.csv").exists():
        sex_summary = pd.read_csv(BASE / "sex_bp_loo_summary.csv")
        save_heatmap(sex_summary, index="region", columns="site",
                     values="mean_auc", out_png="sex_bp_loo_heatmap.png",
                     title="Sex RF Mean AUC by Site and Region")
        save_bar(sex_summary.sort_values(["site", "mean_auc"], ascending=[True, False]),
                 x="region", y="mean_auc", hue="site",
                 out_png="sex_bp_loo_barplot.png", title="Sex RF Mean AUC by Site")

        if (BASE / "sex_bp_loo_predictions.csv").exists():
            preds = pd.read_csv(BASE / "sex_bp_loo_predictions.csv")
            fold = preds.groupby(["site", "continent_group", "fold"]).apply(
                lambda x: pd.Series({
                    "fold_auc": metrics.roc_auc_score(x["true_label"], x["pred_prob"])
                    if x["true_label"].nunique() > 1 else float("nan")
                })).reset_index().rename(columns={"continent_group": "region"}).dropna()
            plt.figure(figsize=(11, 5))
            sns.boxplot(data=fold, x="region", y="fold_auc", hue="site", showfliers=False)
            plt.ylim(0, 1); plt.xticks(rotation=30, ha="right")
            plt.title("Sex RF Fold-level AUC Distribution")
            plt.tight_layout()
            plt.savefig(BASE / "sex_bp_loo_fold_auc_boxplot.png", dpi=300)
            plt.close()

    # Disease RF
    disease_summary_path = DISEASE_DIR / "disease_rf_summary.csv"
    if disease_summary_path.exists():
        disease_summary = pd.read_csv(disease_summary_path)
        plt.figure(figsize=(7, 4))
        sns.barplot(data=disease_summary.sort_values("mean_auc", ascending=False),
                     x="disease", y="mean_auc", color="#4C78A8")
        plt.ylim(0, 1); plt.xticks(rotation=30, ha="right")
        plt.title("Disease RF Mean AUC")
        plt.tight_layout()
        plt.savefig(DISEASE_DIR / "disease_rf_barplot.png", dpi=300)
        plt.close()

    # Region OVO
    if (BASE / "region_ovo_summary.csv").exists():
        ovo_summary = pd.read_csv(BASE / "region_ovo_summary.csv")
        for site, sdf in ovo_summary.groupby("site"):
            regions = sorted(set(sdf["region_1"]).union(set(sdf["region_2"])))
            mat = pd.DataFrame(index=regions, columns=regions, dtype=float)
            for _, row in sdf.iterrows():
                mat.loc[row["region_1"], row["region_2"]] = row["mean_auc"]
                mat.loc[row["region_2"], row["region_1"]] = row["mean_auc"]
            for r in regions: mat.loc[r, r] = 1.0
            plt.figure(figsize=(8, 6))
            sns.heatmap(mat, annot=True, fmt=".3f", cmap="RdYlBu_r", vmin=0, vmax=1)
            plt.title(f"One-vs-One Geography RF Mean AUC: {site}")
            plt.tight_layout()
            plt.savefig(BASE / f"region_ovo_heatmap_{site}.png", dpi=300)
            plt.close()

        if (BASE / "region_ovo_predictions.csv").exists():
            ovo_preds = pd.read_csv(BASE / "region_ovo_predictions.csv")
            pair_fold = ovo_preds.groupby(["site", "region_1", "region_2", "fold"]).apply(
                lambda x: pd.Series({
                    "fold_auc": metrics.roc_auc_score(x["true_label"], x["pred_prob"])
                    if x["true_label"].nunique() > 1 else float("nan")
                })).reset_index().dropna()
            pair_fold["pair"] = pair_fold["region_1"] + " vs " + pair_fold["region_2"]
            plt.figure(figsize=(14, 6))
            sns.boxplot(data=pair_fold, x="pair", y="fold_auc", hue="site", showfliers=False)
            plt.ylim(0, 1); plt.xticks(rotation=60, ha="right")
            plt.title("One-vs-One Geography RF Fold-level AUC Distribution")
            plt.tight_layout()
            plt.savefig(BASE / "region_ovo_fold_auc_boxplot.png", dpi=300)
            plt.close()

    # ========================================================================
    # Site OVA
    # ========================================================================
    site_ova_file = BASE / "site_ova_summary.csv"
    if site_ova_file.exists():
        site_summary = pd.read_csv(site_ova_file)
        # OVA has single target column, use barplot only (heatmap only makes sense for pairwise)
        plt.figure(figsize=(8, 5))
        sns.barplot(data=site_summary.sort_values("mean_auc", ascending=False),
                     x="target_site", y="mean_auc", color="#4C78A8")
        plt.ylim(0.5, 1.0); plt.title("Site One-vs-All RF Mean AUC")
        plt.tight_layout()
        plt.savefig(BASE / "site_ova_barplot.png", dpi=300)
        plt.close()

        if (BASE / "site_ova_predictions.csv").exists():
            preds = pd.read_csv(BASE / "site_ova_predictions.csv")
            fold = preds.groupby(["target_site", "fold"]).apply(
                lambda x: pd.Series({
                    "fold_auc": metrics.roc_auc_score(x["true_label"], x["pred_prob"])
                    if x["true_label"].nunique() > 1 else float("nan")
                })).reset_index().dropna()
            plt.figure(figsize=(10, 5))
            sns.boxplot(data=fold, x="target_site", y="fold_auc", showfliers=False,
                        order=site_summary.sort_values("mean_auc", ascending=False)["target_site"].tolist())
            plt.ylim(0.5, 1.0)
            plt.title("Site RF Fold-level AUC Distribution (One-vs-All)")
            plt.tight_layout()
            plt.savefig(BASE / "site_ova_fold_auc_boxplot.png", dpi=300)
            plt.close()

    # ========================================================================
    # Site OVO
    # ========================================================================
    site_ovo_file = BASE / "site_ovo_summary.csv"
    if site_ovo_file.exists():
        ovo_summary = pd.read_csv(site_ovo_file)
        sites = sorted(set(ovo_summary["site_1"]).union(set(ovo_summary["site_2"])))
        mat = pd.DataFrame(index=sites, columns=sites, dtype=float)
        for _, row in ovo_summary.iterrows():
            mat.loc[row["site_1"], row["site_2"]] = row["mean_auc"]
            mat.loc[row["site_2"], row["site_1"]] = row["mean_auc"]
        for s in sites:
            mat.loc[s, s] = 1.0
        # Mask upper triangle + diagonal for lower-triangle-only display
        mask = np.triu(np.ones_like(mat, dtype=bool), k=0)
        plt.figure(figsize=(8, 6))
        sns.heatmap(mat, annot=True, fmt=".3f", cmap="RdYlBu_r", vmin=0, vmax=1,
                    mask=mask, linewidths=0.5, cbar_kws={"label": "Mean AUC"})
        plt.title("Site One-vs-One RF Mean AUC")
        plt.tight_layout()
        plt.savefig(BASE / "site_ovo_heatmap.png", dpi=300)
        plt.close()

        if (BASE / "site_ovo_predictions.csv").exists():
            preds = pd.read_csv(BASE / "site_ovo_predictions.csv")
            preds["pair"] = preds["site_1"] + " vs " + preds["site_2"]
            fold = preds.groupby(["pair", "fold"]).apply(
                lambda x: pd.Series({
                    "fold_auc": metrics.roc_auc_score(x["true_label"], x["pred_prob"])
                    if x["true_label"].nunique() > 1 else float("nan")
                })).reset_index().dropna()
            plt.figure(figsize=(16, 6))
            sns.boxplot(data=fold, x="pair", y="fold_auc", showfliers=False)
            plt.ylim(0, 1); plt.xticks(rotation=45, ha="right")
            plt.title("Site One-vs-One RF Fold-level AUC Distribution")
            plt.tight_layout()
            plt.savefig(BASE / "site_ovo_fold_auc_boxplot.png", dpi=300)
            plt.close()

    print("Plotting done.")


if __name__ == "__main__":
    main()
