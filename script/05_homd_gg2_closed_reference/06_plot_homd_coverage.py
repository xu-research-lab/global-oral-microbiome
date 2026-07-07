#!/usr/bin/env python3
"""Generate HOMD coverage figures from qcov sensitivity outputs."""

from __future__ import annotations

import argparse
import json
import os
import re

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

DEFAULT_OUTDIR = "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/region_aware_pipeline"

sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)
plt.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

METRIC_LABELS = {
    "derep": "Exact sequence",
    "asv": "ASV-weighted",
    "reads": "Read-weighted",
    "prevalence": "Prevalence-weighted",
}
METRIC_ORDER = ["derep", "asv", "reads", "prevalence"]
IDENTITY_ORDER = [100.0, 99.5, 99.0, 97.0]
QCOV_ORDER = [100.0, 99.0, 95.0, 80.0]
COLORS = {
    "Exact sequence": "#355C7D",
    "ASV-weighted": "#6C5B7B",
    "Read-weighted": "#C06C84",
    "Prevalence-weighted": "#F67280",
}


def log(msg):
    print(f"[plot] {msg}", flush=True)


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def capture_col(identity, qcov):
    return f"captured_id{str(identity).replace('.', 'p')}_qcov{str(qcov).replace('.', 'p')}"


def savefig(fig, path_base):
    fig.tight_layout()
    fig.savefig(path_base + ".pdf", bbox_inches="tight")
    fig.savefig(path_base + ".png", bbox_inches="tight")
    plt.close(fig)


def stats_to_long(stats):
    rows = []
    denom = stats["denominators"]
    for identity in IDENTITY_ORDER:
        for qcov in QCOV_ORDER:
            col = capture_col(identity, qcov)
            if col not in stats["capture_totals"]:
                continue
            val = stats["capture_totals"][col]
            rows.extend([
                {"identity": identity, "qcov": qcov, "metric": "derep", "fraction": val["derep"] / stats["derep_sequences"]},
                {"identity": identity, "qcov": qcov, "metric": "asv", "fraction": val["asvs"] / denom["asvs"]},
                {"identity": identity, "qcov": qcov, "metric": "reads", "fraction": val["reads"] / denom["reads"]},
                {"identity": identity, "qcov": qcov, "metric": "prevalence", "fraction": val["prevalence"] / denom["prevalence"]},
            ])
    df = pd.DataFrame(rows)
    df["metric_label"] = df["metric"].map(METRIC_LABELS)
    return df


def plot_global_sensitivity(stats, figdir):
    df = stats_to_long(stats)
    main = df[df["qcov"].eq(99.0)].copy()
    main["identity_label"] = main["identity"].map(lambda x: f">={x:g}%" if x < 100 else "100%")

    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    sns.lineplot(
        data=main,
        x="identity_label",
        y="fraction",
        hue="metric_label",
        style="metric_label",
        markers=True,
        dashes=False,
        palette=COLORS,
        linewidth=2.2,
        markersize=7,
        ax=ax,
    )
    ax.set_ylim(0, 1)
    ax.set_xlabel("HOMD identity threshold (query coverage >=99%)")
    ax.set_ylabel("Captured fraction")
    ax.set_title("Global HOMD capture across identity thresholds")
    ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(1.0))
    ax.legend(title="Weighting", frameon=False, loc="lower left")
    savefig(fig, os.path.join(figdir, "fig1_global_identity_thresholds_qcov99"))

    id99 = df[df["identity"].eq(99.0)].copy()
    id99["qcov_label"] = id99["qcov"].map(lambda x: f">={x:g}%")
    fig, ax = plt.subplots(figsize=(8.0, 4.4))
    sns.barplot(
        data=id99,
        x="qcov_label",
        y="fraction",
        hue="metric_label",
        hue_order=[METRIC_LABELS[m] for m in METRIC_ORDER],
        palette=COLORS,
        ax=ax,
    )
    ax.set_ylim(0, 1)
    ax.set_xlabel("Query coverage threshold (identity >=99%)")
    ax.set_ylabel("Captured fraction")
    ax.set_title("Capture sensitivity to query coverage threshold")
    ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(1.0))
    ax.legend(title="Weighting", frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    savefig(fig, os.path.join(figdir, "fig1b_qcov_sensitivity_id99"))


def plot_region_coverage(grouped, figdir, top_n=14):
    sub = grouped[grouped["grouping"].eq("sequencing_regions")].copy()
    frac_cols = {
        "derep": "captured_id99p0_qcov99p0_derep_fraction",
        "asv": "captured_id99p0_qcov99p0_asv_fraction",
        "reads": "captured_id99p0_qcov99p0_read_fraction",
        "prevalence": "captured_id99p0_qcov99p0_prevalence_fraction",
    }
    sub["level"] = sub["level"].fillna("NA").replace({"nan": "NA", "NaN": "NA"})
    sub = sub.sort_values("source_asvs", ascending=False).head(top_n)
    order = sub.sort_values("source_asvs", ascending=True)["level"].tolist()
    long = sub.melt(
        id_vars=["level", "source_asvs", "total_reads"],
        value_vars=list(frac_cols.values()),
        var_name="metric",
        value_name="fraction",
    )
    long["metric_label"] = long["metric"].map({v: METRIC_LABELS[k] for k, v in frac_cols.items()})

    long["level"] = pd.Categorical(long["level"], categories=order, ordered=True)
    long = long.sort_values("level")
    fig, ax = plt.subplots(figsize=(8.8, max(5.2, 0.42 * len(order))))
    sns.scatterplot(
        data=long,
        x="fraction",
        y="level",
        hue="metric_label",
        size="source_asvs",
        sizes=(35, 220),
        palette=COLORS,
        ax=ax,
    )
    ax.set_xlim(0, 1)
    ax.set_xlabel("HOMD captured fraction (identity >=99%, qcov >=99%)")
    ax.set_ylabel("Sequencing region")
    ax.set_title("Region-aware HOMD capture")
    ax.xaxis.set_major_formatter(mpl.ticker.PercentFormatter(1.0))
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    savefig(fig, os.path.join(figdir, "fig2_region_aware_capture_dotplot"))


def identity_bin(series):
    return pd.cut(
        series.fillna(0),
        bins=[-0.1, 80, 90, 97, 99, 100.0001],
        labels=["No/low hit (<80)", "80-90", "90-97", "97-99", ">=99"],
        include_lowest=True,
        right=False,
    )


def plot_reproducibility_scatter(perseq_path, figdir):
    cols = ["derep_id", "n_asvs", "n_bioprojects", "total_reads", "prevalence_samples", "identity", "query_cov"]
    df = pd.read_csv(perseq_path, sep="\t", usecols=cols)
    df["identity_bin"] = identity_bin(df["identity"])
    df["log_reads"] = np.log10(df["total_reads"] + 1)
    df["log_prevalence"] = np.log10(df["prevalence_samples"] + 1)
    df["is_candidate_like"] = (df["n_bioprojects"] >= 3) & (df["total_reads"] >= 1000) & (df["prevalence_samples"] >= 5)

    # Plot all highly recurrent/abundant points plus a deterministic sample of the rest.
    keep = df[df["is_candidate_like"]].copy()
    rest = df[~df["is_candidate_like"]].sample(n=min(120000, (~df["is_candidate_like"]).sum()), random_state=7)
    plot_df = pd.concat([rest, keep], ignore_index=True)

    palette = {
        "No/low hit (<80)": "#8C8C8C",
        "80-90": "#DDAA33",
        "90-97": "#BB5566",
        "97-99": "#004488",
        ">=99": "#228833",
    }
    fig, ax = plt.subplots(figsize=(8.2, 5.8))
    sns.scatterplot(
        data=plot_df,
        x="n_bioprojects",
        y="log_reads",
        hue="identity_bin",
        size="n_asvs",
        sizes=(8, 80),
        alpha=0.42,
        linewidth=0,
        palette=palette,
        ax=ax,
    )
    ax.axvline(3, color="black", linestyle="--", linewidth=1)
    ax.axhline(np.log10(1000 + 1), color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("Number of bioprojects sharing exact sequence")
    ax.set_ylabel("log10(total reads + 1)")
    ax.set_title("Reproducibility and abundance of ASVs by HOMD identity")
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", title="Best HOMD identity")
    savefig(fig, os.path.join(figdir, "fig3_reproducibility_abundance_homd_identity"))


def heatmap_from_grouped(grouped, grouping, figdir, top_n=14):
    sub = grouped[grouped["grouping"].eq(grouping)].copy()
    sub["level"] = sub["level"].fillna("NA").replace({"nan": "NA", "NaN": "NA"})
    sub = sub.sort_values("source_asvs", ascending=False).head(top_n)
    value_cols = [
        "captured_id100p0_qcov99p0_asv_fraction",
        "captured_id99p5_qcov99p0_asv_fraction",
        "captured_id99p0_qcov99p0_asv_fraction",
        "captured_id97p0_qcov99p0_asv_fraction",
    ]
    labels = ["100%", ">=99.5%", ">=99%", ">=97%"]
    mat = sub.set_index("level")[value_cols]
    mat.columns = labels
    mat = mat.sort_values(">=99%", ascending=True)

    fig, ax = plt.subplots(figsize=(6.4, max(4.8, 0.34 * len(mat))))
    sns.heatmap(
        mat,
        cmap="YlGnBu",
        vmin=0,
        vmax=1,
        annot=True,
        fmt=".0%",
        cbar_kws={"label": "ASV-weighted captured fraction"},
        ax=ax,
    )
    ax.set_xlabel("HOMD identity threshold (qcov >=99%)")
    ax.set_ylabel(grouping.replace("_", " "))
    ax.set_title(f"HOMD capture by {grouping.replace('_', ' ')}")
    savefig(fig, os.path.join(figdir, f"fig4_heatmap_{grouping}_identity_capture"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--figdir", default=None)
    args = parser.parse_args()

    figdir = args.figdir or os.path.join(args.outdir, "figures")
    ensure_dir(figdir)

    stats_path = os.path.join(args.outdir, "qcov_sensitivity_stats.json")
    grouped_path = os.path.join(args.outdir, "qcov_sensitivity_grouped_summary.tsv")
    perseq_path = os.path.join(args.outdir, "qcov_sensitivity_by_derep_sequence.tsv")

    with open(stats_path) as handle:
        stats = json.load(handle)
    grouped = pd.read_csv(grouped_path, sep="\t")

    log("plotting global sensitivity")
    plot_global_sensitivity(stats, figdir)
    log("plotting region-aware coverage")
    plot_region_coverage(grouped, figdir)
    log("plotting reproducibility/abundance scatter")
    plot_reproducibility_scatter(perseq_path, figdir)
    log("plotting heatmaps")
    for grouping in ["continents", "countries", "statuses"]:
        heatmap_from_grouped(grouped, grouping, figdir)
    log(f"figures written to {figdir}")


if __name__ == "__main__":
    main()
