#!/usr/bin/env python3
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FuncFormatter

BASE = Path("/beegfs/huangxiaochang/global_oral/new_analysis")
INPUT = BASE / "metadata_all_for_qc_raw.tsv"
OUTDIR = BASE / "QC"
PREFIX = "metadata_raw_qc_panels"

RAW_READS_THRESHOLD = 1000
DECHIMERA_READS_THRESHOLD = 1000
CHIMERA_THRESHOLD = 25.0
ASV_THRESHOLD = 30
LOW_ASV_RATIO_THRESHOLD = 0.50
HIGH_CHIMERA_RATIO_THRESHOLD = 0.30
OBSERVED_RETENTION_THRESHOLD = 0.50

PALETTE = [
    "#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#76B7B2", "#B07AA1",
    "#EDC948", "#9C755F", "#BAB0AC", "#1F77B4", "#FF7F0E", "#2CA02C"
]


def parse_numeric(series):
    return pd.to_numeric(series.replace({"NA": np.nan, "": np.nan}), errors="coerce")


def comma_formatter(x, _pos):
    if x >= 1:
        return f"{int(x):,}"
    return f"{x:g}"


def percent_formatter(x, _pos):
    return f"{x * 100:.0f}%"


def write_csv(path, df):
    df.to_csv(path, index=False)


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(INPUT, sep="	", dtype=str)
    df["study_id"] = df["bioproject_id"].astype(str)
    df["sample_id_norm"] = df["sample_id"].astype(str)
    df["total_reads_num"] = parse_numeric(df["total_reads"])
    df["chim_perc_num"] = parse_numeric(df["chim_perc"])
    df["asv_count_num"] = parse_numeric(df["asv_count"])
    df["observed_num"] = parse_numeric(df["Observed"])
    df["dechimera_reads"] = df["total_reads_num"] * (1.0 - df["chim_perc_num"] / 100.0)
    df.loc[df["total_reads_num"].isna() | df["chim_perc_num"].isna(), "dechimera_reads"] = np.nan
    df["low_asv_sample"] = df["asv_count_num"] <= ASV_THRESHOLD
    df["high_chimera_sample"] = df["chim_perc_num"] > CHIMERA_THRESHOLD
    df["observed_present"] = df["observed_num"].notna()

    study_metrics = df.groupby("study_id", dropna=False).agg(
        metadata_sample_count=("sample_id_norm", "size"),
        mean_raw_reads=("total_reads_num", "mean"),
        median_raw_reads=("total_reads_num", "median"),
        low_asv_sample_count=("low_asv_sample", "sum"),
        high_chimera_sample_count=("high_chimera_sample", "sum"),
        observed_present_count=("observed_present", "sum"),
    ).reset_index()
    study_metrics["low_asv_sample_ratio"] = study_metrics["low_asv_sample_count"] / study_metrics["metadata_sample_count"]
    study_metrics["high_chimera_sample_ratio"] = study_metrics["high_chimera_sample_count"] / study_metrics["metadata_sample_count"]
    study_metrics["closed_reference_retention_ratio"] = study_metrics["observed_present_count"] / study_metrics["metadata_sample_count"]
    study_metrics["low_asv_flag"] = study_metrics["low_asv_sample_ratio"] > LOW_ASV_RATIO_THRESHOLD
    study_metrics["high_chimera_flag"] = study_metrics["high_chimera_sample_ratio"] > HIGH_CHIMERA_RATIO_THRESHOLD
    study_metrics["low_observed_retention_flag"] = study_metrics["closed_reference_retention_ratio"] < OBSERVED_RETENTION_THRESHOLD
    study_metrics = study_metrics.sort_values(["metadata_sample_count", "study_id"], ascending=[False, True]).reset_index(drop=True)
    study_metrics["study_rank"] = np.arange(1, len(study_metrics) + 1)

    write_csv(OUTDIR / f"{PREFIX}_study_metrics.csv", study_metrics)
    write_csv(OUTDIR / f"{PREFIX}_low_asv_flagged_studies.csv", study_metrics[study_metrics["low_asv_flag"]])
    write_csv(OUTDIR / f"{PREFIX}_high_chimera_flagged_studies.csv", study_metrics[study_metrics["high_chimera_flag"]])
    write_csv(OUTDIR / f"{PREFIX}_low_observed_retention_studies.csv", study_metrics[study_metrics["low_observed_retention_flag"]])

    fig, axes = plt.subplots(5, 1, figsize=(30, 16), sharex=True, constrained_layout=True)
    x = np.arange(len(study_metrics))

    axes[0].bar(x, study_metrics["mean_raw_reads"], color="#6B6ECF", width=0.85)
    axes[0].axhline(RAW_READS_THRESHOLD, linestyle="--", color="#D62728", linewidth=1.0)
    axes[0].set_yscale("log")
    axes[0].yaxis.set_major_formatter(FuncFormatter(comma_formatter))
    axes[0].set_ylabel("Mean reads")
    axes[0].set_title("Mean raw reads per study (log10 scale; red line: 1,000 reads)")

    axes[1].bar(
        x,
        study_metrics["low_asv_sample_ratio"],
        color=np.where(study_metrics["low_asv_flag"], "#D62728", "#4C78A8"),
        width=0.85,
    )
    axes[1].axhline(LOW_ASV_RATIO_THRESHOLD, linestyle="--", color="#D62728", linewidth=1.0)
    axes[1].yaxis.set_major_formatter(FuncFormatter(percent_formatter))
    axes[1].set_ylabel("Sample ratio")
    axes[1].set_title("Samples with ASV count <= 30 (red: study ratio > 50%)")

    axes[2].bar(x, study_metrics["closed_reference_retention_ratio"], color="#59A14F", width=0.85)
    axes[2].axhline(OBSERVED_RETENTION_THRESHOLD, linestyle="--", color="#7F7F7F", linewidth=1.0)
    axes[2].yaxis.set_major_formatter(FuncFormatter(percent_formatter))
    axes[2].set_ylabel("Observed present")
    axes[2].set_title("Closed-reference alpha retention ratio (grey line: 50%; informational)")

    axes[3].bar(
        x,
        study_metrics["high_chimera_sample_ratio"],
        color=np.where(study_metrics["high_chimera_flag"], "#D62728", "#9C755F"),
        width=0.85,
    )
    axes[3].axhline(HIGH_CHIMERA_RATIO_THRESHOLD, linestyle="--", color="#D62728", linewidth=1.0)
    axes[3].yaxis.set_major_formatter(FuncFormatter(percent_formatter))
    axes[3].set_ylabel("Sample ratio")
    axes[3].set_title("Samples with chimera percent > 25% (red: study ratio > 30%)")

    axes[4].bar(x, study_metrics["metadata_sample_count"], color="#F28E2B", width=0.85)
    axes[4].set_yscale("log")
    axes[4].yaxis.set_major_formatter(FuncFormatter(comma_formatter))
    axes[4].set_ylabel("Samples")
    axes[4].set_title("Samples in metadata_all_for_qc_raw.tsv (log10 scale)")
    axes[4].set_xticks(x)
    axes[4].set_xticklabels(study_metrics["study_id"], rotation=90, fontsize=8)
    axes[4].set_xlabel("Study")

    for ax in axes[:-1]:
        ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    for ax in axes:
        ax.grid(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    study_png = OUTDIR / f"{PREFIX}_study_level.png"
    study_pdf = OUTDIR / f"{PREFIX}_study_level.pdf"
    fig.savefig(study_png, dpi=300)
    with PdfPages(study_pdf) as pdf:
        pdf.savefig(fig)
    plt.close(fig)

    study_order = study_metrics[["study_id", "metadata_sample_count"]].copy()
    study_order["study_rank"] = np.arange(1, len(study_order) + 1)
    study_order["study_color"] = [PALETTE[(i - 1) % len(PALETTE)] for i in study_order["study_rank"]]
    write_csv(OUTDIR / f"{PREFIX}_sample_study_order.csv", study_order)

    sample_df = df.merge(study_order, on="study_id", how="left")
    sample_df = sample_df.sort_values(["study_rank", "total_reads_num", "sample_id_norm"], ascending=[True, False, True]).reset_index(drop=True)
    sample_df["sample_rank"] = np.arange(1, len(sample_df) + 1)
    write_csv(
        OUTDIR / f"{PREFIX}_sample_ordered_values.csv",
        sample_df[[
            "study_id", "sample_id_norm", "sample_rank", "total_reads_num", "dechimera_reads",
            "chim_perc_num", "asv_count_num", "observed_num", "metadata_sample_count"
        ]].rename(columns={
            "sample_id_norm": "sample_id",
            "total_reads_num": "raw_reads",
            "chim_perc_num": "chim_perc",
            "asv_count_num": "asv_count",
            "observed_num": "Observed"
        })
    )

    threshold_summary = pd.DataFrame([
        {"metric": "raw_reads", "filter_rule": f"raw_reads < {RAW_READS_THRESHOLD}", "threshold": RAW_READS_THRESHOLD,
         "failed_samples": int((sample_df["total_reads_num"] < RAW_READS_THRESHOLD).sum()),
         "non_missing_samples": int(sample_df["total_reads_num"].notna().sum())},
        {"metric": "dechimera_reads", "filter_rule": f"dechimera_reads < {DECHIMERA_READS_THRESHOLD}", "threshold": DECHIMERA_READS_THRESHOLD,
         "failed_samples": int((sample_df["dechimera_reads"] < DECHIMERA_READS_THRESHOLD).sum()),
         "non_missing_samples": int(sample_df["dechimera_reads"].notna().sum())},
        {"metric": "chim_perc", "filter_rule": f"chim_perc > {CHIMERA_THRESHOLD}", "threshold": CHIMERA_THRESHOLD,
         "failed_samples": int((sample_df["chim_perc_num"] > CHIMERA_THRESHOLD).sum()),
         "non_missing_samples": int(sample_df["chim_perc_num"].notna().sum())},
        {"metric": "asv_count", "filter_rule": f"asv_count <= {ASV_THRESHOLD}", "threshold": ASV_THRESHOLD,
         "failed_samples": int((sample_df["asv_count_num"] <= ASV_THRESHOLD).sum()),
         "non_missing_samples": int(sample_df["asv_count_num"].notna().sum())},
        {"metric": "Observed_na", "filter_rule": "Observed is NA", "threshold": np.nan,
         "failed_samples": int(sample_df["observed_num"].isna().sum()),
         "non_missing_samples": int(len(sample_df))},
    ])
    threshold_summary["failed_ratio"] = threshold_summary["failed_samples"] / threshold_summary["non_missing_samples"]
    write_csv(OUTDIR / f"{PREFIX}_sample_threshold_summary.csv", threshold_summary)

    fig, axes = plt.subplots(3, 1, figsize=(28, 12), sharex=True, constrained_layout=True)
    x = sample_df["sample_rank"].to_numpy()
    colors = sample_df["study_color"].to_numpy()

    axes[0].bar(x, sample_df["total_reads_num"].to_numpy(), color=colors, width=1.0)
    axes[0].axhline(RAW_READS_THRESHOLD, linestyle="--", color="#D62728", linewidth=1.0)
    axes[0].axhline(DECHIMERA_READS_THRESHOLD, linestyle=":", color="#D62728", linewidth=1.0)
    axes[0].set_yscale("log")
    axes[0].yaxis.set_major_formatter(FuncFormatter(comma_formatter))
    axes[0].set_ylabel("Reads")
    axes[0].set_title("Sample raw reads by study/sample order")

    axes[1].bar(x, sample_df["chim_perc_num"].to_numpy(), color=colors, width=1.0)
    axes[1].axhline(CHIMERA_THRESHOLD, linestyle="--", color="#D62728", linewidth=1.0)
    axes[1].set_ylim(0, 100)
    axes[1].set_ylabel("Chimera %")
    axes[1].set_title("Sample chimera percentage (red line: 25%)")

    axes[2].bar(x, sample_df["asv_count_num"].to_numpy(), color=colors, width=1.0)
    axes[2].axhline(ASV_THRESHOLD, linestyle="--", color="#D62728", linewidth=1.0)
    axes[2].set_yscale("log")
    axes[2].yaxis.set_major_formatter(FuncFormatter(comma_formatter))
    axes[2].set_ylabel("ASV count")
    axes[2].set_title("Sample ASV count (red line: 30)")
    axes[2].set_xlabel("Samples ordered by study size and raw reads")

    boundaries = sample_df.groupby("study_id")["sample_rank"].agg(["min", "max"]).reset_index()
    for ax in axes:
        for _, row in boundaries.iloc[:-1].iterrows():
            ax.axvline(row["max"] + 0.5, color="grey35", linewidth=0.15, alpha=0.35)
        ax.grid(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    axes[-1].tick_params(axis="x", which="both", bottom=False, labelbottom=False)

    sample_png = OUTDIR / f"{PREFIX}_sample_level.png"
    sample_pdf = OUTDIR / f"{PREFIX}_sample_level.pdf"
    fig.savefig(sample_png, dpi=300)
    with PdfPages(sample_pdf) as pdf:
        pdf.savefig(fig)
    plt.close(fig)

    print(f"Study metrics: {OUTDIR / f'{PREFIX}_study_metrics.csv'}")
    print(f"Study plot PNG: {study_png}")
    print(f"Study plot PDF: {study_pdf}")
    print(f"Sample plot PNG: {sample_png}")
    print(f"Sample plot PDF: {sample_pdf}")
    print(f"Threshold summary: {OUTDIR / f'{PREFIX}_sample_threshold_summary.csv'}")


if __name__ == "__main__":
    main()
