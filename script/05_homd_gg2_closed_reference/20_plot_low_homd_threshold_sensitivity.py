#!/usr/bin/env python3
"""Plot sensitivity of low-HOMD bacterial ASV sets across recurrence/prevalence thresholds."""

from __future__ import annotations

import argparse
import os

import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

DEFAULT_OUTDIR = "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/region_aware_pipeline"

sns.set_theme(style="whitegrid", context="paper", font_scale=1.12)
plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "figure.dpi": 150, "savefig.dpi": 300, "axes.spines.top": False, "axes.spines.right": False})

QC_COLORS = {
    "bacterial": "#2F7D73",
    "unassigned": "#8A8A8A",
    "chloroplast": "#8CB369",
    "mitochondria": "#D99058",
    "archaea": "#7B6EA8",
    "other_non_bacterial": "#B85C5C",
    "eukaryota": "#B85C5C",
}

REGION_COLORS = {
    "Asia": "#CC6677",
    "North America": "#4477AA",
    "Europe": "#228833",
    "Asia;Europe": "#AA3377",
    "Asia;North America": "#66CCEE",
    "Asia;Europe;North America": "#EE7733",
    "Asia;Europe;North America;Oceania": "#BBBBBB",
}


def threshold_label(row):
    return f"BP≥{int(row.min_bioprojects)}, Subj≥{int(row.min_subject_prevalence)}, Reads≥{int(row.min_total_reads)}"


def savefig(fig, path_base):
    fig.tight_layout()
    fig.savefig(path_base + ".pdf", bbox_inches="tight")
    fig.savefig(path_base + ".png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    sens_dir = os.path.join(args.outdir, "threshold_sensitivity")
    figdir = os.path.join(args.outdir, "figures")
    os.makedirs(figdir, exist_ok=True)

    summary = pd.read_csv(os.path.join(sens_dir, "low_homd_threshold_sensitivity_summary.tsv"), sep="\t")
    summary["threshold_display"] = summary.apply(threshold_label, axis=1)
    order = summary["threshold_display"].tolist()

    # QC class composition across thresholds.
    qc_rows = []
    for row in summary.itertuples(index=False):
        qc = pd.read_csv(row.qc_summary, sep="\t")
        qc["threshold_label"] = row.threshold_label
        qc["threshold_display"] = row.threshold_display
        qc_rows.append(qc)
    qc_long = pd.concat(qc_rows, ignore_index=True)
    qc_long["threshold_display"] = pd.Categorical(qc_long["threshold_display"], categories=order, ordered=True)
    qc_long.to_csv(os.path.join(sens_dir, "low_homd_threshold_sensitivity_qc_long.tsv"), sep="\t", index=False)

    # Summary figure.
    fig, axes = plt.subplots(2, 2, figsize=(13.4, 8.4))
    ax = axes[0, 0]
    count_long = summary.melt(
        id_vars=["threshold_display"],
        value_vars=["filtered_low_homd_asvs", "gg2_bacterial_asvs"],
        var_name="metric",
        value_name="count",
    )
    count_long["metric"] = count_long["metric"].map({"filtered_low_homd_asvs": "All filtered low-HOMD ASVs", "gg2_bacterial_asvs": "GG2 bacterial ASVs"})
    sns.lineplot(data=count_long, x="threshold_display", y="count", hue="metric", marker="o", linewidth=2.2, ax=ax, palette=["#444444", "#2F7D73"])
    ax.set_title("A. Retained ASV counts", loc="left", fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("ASV count")
    ax.tick_params(axis="x", rotation=35)
    ax.legend(title="", frameon=False)

    ax = axes[0, 1]
    qc_pivot = qc_long.pivot_table(index="threshold_display", columns="gg2_qc_class", values="derep_sequences", aggfunc="sum", fill_value=0)
    qc_pivot = qc_pivot.reindex(order)
    bottom = None
    x = range(len(qc_pivot))
    for cls in [c for c in ["bacterial", "unassigned", "chloroplast", "mitochondria", "archaea", "other_non_bacterial", "eukaryota"] if c in qc_pivot.columns]:
        vals = qc_pivot[cls].to_numpy()
        ax.bar(x, vals, bottom=bottom, color=QC_COLORS.get(cls, "#BBBBBB"), label=cls, width=0.72)
        bottom = vals if bottom is None else bottom + vals
    ax.set_xticks(list(x))
    ax.set_xticklabels(order, rotation=35, ha="right")
    ax.set_title("B. GG2 QC class composition", loc="left", fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("ASV count")
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")

    ax = axes[1, 0]
    burden_long = summary.melt(
        id_vars=["threshold_display"],
        value_vars=["bacterial_total_reads", "bacterial_subject_prevalence_burden"],
        var_name="metric",
        value_name="value",
    )
    burden_long["metric"] = burden_long["metric"].map({"bacterial_total_reads": "Total reads", "bacterial_subject_prevalence_burden": "Subject-prevalence burden"})
    sns.lineplot(data=burden_long, x="threshold_display", y="value", hue="metric", marker="o", linewidth=2.2, ax=ax, palette=["#9A4F2D", "#3B6EA8"])
    ax.set_title("C. Bacterial low-HOMD burden", loc="left", fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("Burden")
    ax.tick_params(axis="x", rotation=35)
    ax.legend(title="", frameon=False)
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, pos: f"{v/1e6:.0f}M" if v >= 1e6 else f"{v/1e3:.0f}K"))

    ax = axes[1, 1]
    summary["bacterial_fraction"] = summary["gg2_bacterial_asvs"] / summary["filtered_low_homd_asvs"]
    summary["unassigned_fraction"] = summary["gg2_unassigned_asvs"] / summary["filtered_low_homd_asvs"]
    frac_long = summary.melt(
        id_vars=["threshold_display"],
        value_vars=["bacterial_fraction", "unassigned_fraction"],
        var_name="metric",
        value_name="fraction",
    )
    frac_long["metric"] = frac_long["metric"].map({"bacterial_fraction": "Bacterial", "unassigned_fraction": "Unassigned"})
    sns.lineplot(data=frac_long, x="threshold_display", y="fraction", hue="metric", marker="o", linewidth=2.2, ax=ax, palette=["#2F7D73", "#777777"])
    ax.set_title("D. Fraction of filtered low-HOMD ASVs", loc="left", fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("Fraction")
    ax.set_ylim(0, 1)
    ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(1.0))
    ax.tick_params(axis="x", rotation=35)
    ax.legend(title="", frameon=False)

    savefig(fig, os.path.join(figdir, "fig12_low_homd_threshold_sensitivity_summary"))

    # Region intersection figure.
    reg_rows = []
    for row in summary.itertuples(index=False):
        reg = pd.read_csv(row.region_intersections, sep="\t")
        reg["threshold_label"] = row.threshold_label
        reg["threshold_display"] = row.threshold_display
        reg_rows.append(reg)
    reg_long = pd.concat(reg_rows, ignore_index=True)
    reg_long.to_csv(os.path.join(sens_dir, "low_homd_threshold_sensitivity_region_intersections_long.tsv"), sep="\t", index=False)

    top_regions = (
        reg_long.groupby("regions", as_index=False)["bacterial_asv_count"].max()
        .sort_values("bacterial_asv_count", ascending=False)
        .head(10)["regions"].tolist()
    )
    plot_reg = reg_long[reg_long["regions"].isin(top_regions)].copy()
    plot_reg["threshold_display"] = pd.Categorical(plot_reg["threshold_display"], categories=order, ordered=True)
    region_order = (
        plot_reg.groupby("regions")["bacterial_asv_count"].max().sort_values(ascending=False).index.tolist()
    )
    plot_reg["regions"] = pd.Categorical(plot_reg["regions"], categories=region_order, ordered=True)

    fig, ax = plt.subplots(figsize=(12.8, 5.8))
    sns.lineplot(data=plot_reg, x="threshold_display", y="bacterial_asv_count", hue="regions", marker="o", linewidth=2.0, ax=ax)
    ax.set_title("Top geographic intersections across threshold sensitivity sets", loc="left", fontweight="bold")
    ax.set_xlabel("Filtering threshold")
    ax.set_ylabel("GG2 bacterial low-HOMD ASV count")
    ax.tick_params(axis="x", rotation=35)
    ax.legend(title="Region intersection", frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    savefig(fig, os.path.join(figdir, "fig13_low_homd_threshold_sensitivity_region_intersections"))

    print(os.path.join(figdir, "fig12_low_homd_threshold_sensitivity_summary.pdf"))
    print(os.path.join(figdir, "fig13_low_homd_threshold_sensitivity_region_intersections.pdf"))


if __name__ == "__main__":
    main()
