#!/usr/bin/env python3
"""
Summarize per-sample read fractions mapped/unmapped by the containment-assisted
GG2 closed-reference profile and draw BioProject-level boxplots.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


DEFAULT_ROOT = "/beegfs/huangxiaochang/global_oral"
DEFAULT_METADATA = f"{DEFAULT_ROOT}/metadata_v2_merge/metadata-v2-merge-oral.csv"
DEFAULT_ASV_ROOT = "/beegfs/dongbiao/oral"
DEFAULT_PROFILE_DIR = f"{DEFAULT_ROOT}/ASV_vs_HOMD_analysis/containment_derep/unified_gg2_closed_reference_profile"
DEFAULT_OUTDIR = f"{DEFAULT_ROOT}/ASV_vs_HOMD_analysis/containment_derep/unified_gg2_closed_reference_profile/figures"


def log(message: str) -> None:
    print(f"[mapfrac] {message}", flush=True)


def retained_bioprojects(metadata: Path) -> list[str]:
    md = pd.read_csv(metadata, low_memory=False, usecols=["bioproject"])
    return sorted(md["bioproject"].dropna().astype(str).str.strip().loc[lambda x: x.ne("")].unique())


def compute_total_reads(asv_root: Path, bioprojects: list[str]) -> pd.DataFrame:
    rows = []
    for i, bp in enumerate(bioprojects, start=1):
        path = asv_root / bp / "results" / "ASVs_counts.tsv"
        if not path.exists():
            continue
        counts = pd.read_csv(path, sep="\t", index_col=0)
        numeric = counts.apply(pd.to_numeric, errors="coerce").fillna(0)
        sample_totals = numeric.sum(axis=0)
        for sample_id, total in sample_totals.items():
            rows.append(
                {
                    "bioproject": bp,
                    "sample_id": str(sample_id),
                    "sample_key": f"{bp}|{sample_id}",
                    "total_reads": float(total),
                }
            )
        if i % 25 == 0:
            log(f"computed total reads for {i}/{len(bioprojects)} BioProjects")
    return pd.DataFrame(rows)


def make_sample_fraction_table(args: argparse.Namespace) -> pd.DataFrame:
    profile_dir = Path(args.profile_dir)
    sample_summary = pd.read_csv(
        profile_dir / "all_samples_gg2_closed_reference_sample_summary.tsv",
        sep="\t",
    )
    sample_summary["sample_key"] = sample_summary["sample_key"].astype(str)
    sample_summary["closed_reference_reads"] = pd.to_numeric(
        sample_summary["closed_reference_reads"], errors="coerce"
    ).fillna(0)

    bps = retained_bioprojects(Path(args.metadata))
    totals = compute_total_reads(Path(args.asv_root), bps)
    df = totals.merge(
        sample_summary[["sample_key", "closed_reference_reads", "closed_reference_features"]],
        on="sample_key",
        how="left",
    )
    df["closed_reference_reads"] = df["closed_reference_reads"].fillna(0)
    df["closed_reference_features"] = df["closed_reference_features"].fillna(0).astype(int)
    df["unmapped_reads"] = df["total_reads"] - df["closed_reference_reads"]
    df.loc[df["unmapped_reads"] < 0, "unmapped_reads"] = 0
    df["mapped_fraction"] = np.where(df["total_reads"] > 0, df["closed_reference_reads"] / df["total_reads"], np.nan)
    df["unmapped_fraction"] = np.where(df["total_reads"] > 0, df["unmapped_reads"] / df["total_reads"], np.nan)
    return df


def plot_box_jitter(df: pd.DataFrame, outdir: Path, metric: str, ylabel: str, prefix: str) -> None:
    plot_df = df.dropna(subset=[metric]).copy()
    order = (
        plot_df.groupby("bioproject")[metric]
        .median()
        .sort_values(ascending=False)
        .index.tolist()
    )
    plot_df["bioproject"] = pd.Categorical(plot_df["bioproject"], categories=order, ordered=True)

    sns.set_theme(style="whitegrid", context="paper")
    n_projects = len(order)
    width = max(18, n_projects * 0.065)
    fig, ax = plt.subplots(figsize=(width, 6.2))
    sns.boxplot(
        data=plot_df,
        x="bioproject",
        y=metric,
        order=order,
        ax=ax,
        color="#D6EAF8",
        fliersize=0,
        linewidth=0.55,
        width=0.62,
    )
    sns.stripplot(
        data=plot_df,
        x="bioproject",
        y=metric,
        order=order,
        ax=ax,
        color="#1B4F72",
        size=1.15,
        alpha=0.42,
        jitter=0.24,
        linewidth=0,
        rasterized=True,
    )
    # Add 30% sample-level threshold
    THRESHOLD = 0.30
    ax.axhline(THRESHOLD, color="#D62728", linestyle="--", linewidth=1.0, alpha=0.7)
    n_total = len(plot_df)
    n_pass = int((plot_df[metric] >= THRESHOLD).sum())
    n_fail = n_total - n_pass
    ax.text(
        0.98, THRESHOLD + 0.02,
        f"30% sample threshold: pass {n_pass}, fail {n_fail}",
        color="#D62728", ha="right", va="bottom", fontsize=8,
        transform=ax.get_yaxis_transform(),
    )

    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("BioProject, ordered by sample-level median")
    ax.set_ylabel(ylabel)
    ax.set_title("Per-sample GG2 closed-reference read fraction by BioProject")
    ax.tick_params(axis="x", labelrotation=90, labelsize=4.5)
    ax.yaxis.set_major_formatter(lambda x, pos: f"{x * 100:.0f}%")
    fig.tight_layout()
    for ext in ["png", "pdf"]:
        fig.savefig(outdir / f"{prefix}.{ext}", dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", default=DEFAULT_METADATA)
    parser.add_argument("--asv-root", default=DEFAULT_ASV_ROOT)
    parser.add_argument("--profile-dir", default=DEFAULT_PROFILE_DIR)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    profile_dir = Path(args.profile_dir)

    df = make_sample_fraction_table(args)
    out_table = profile_dir / "all_samples_gg2_closed_reference_mapping_fraction.tsv"
    df.to_csv(out_table, sep="\t", index=False)

    project_summary = (
        df.groupby("bioproject", as_index=False)
        .agg(
            samples=("sample_key", "nunique"),
            median_mapped_fraction=("mapped_fraction", "median"),
            mean_mapped_fraction=("mapped_fraction", "mean"),
            median_unmapped_fraction=("unmapped_fraction", "median"),
            total_reads=("total_reads", "sum"),
            closed_reference_reads=("closed_reference_reads", "sum"),
            unmapped_reads=("unmapped_reads", "sum"),
        )
    )
    project_summary["overall_mapped_fraction"] = project_summary["closed_reference_reads"] / project_summary["total_reads"]
    project_summary = project_summary.sort_values("median_mapped_fraction", ascending=False)
    project_summary.to_csv(profile_dir / "all_projects_gg2_closed_reference_mapping_fraction_summary.tsv", sep="\t", index=False)

    plot_box_jitter(
        df,
        outdir,
        metric="mapped_fraction",
        ylabel="Reads mapped to GG2 closed-reference / total reads",
        prefix="fig_gg2_closed_reference_mapped_read_fraction_by_bioproject",
    )
    plot_box_jitter(
        df,
        outdir,
        metric="unmapped_fraction",
        ylabel="Reads not mapped to GG2 closed-reference / total reads",
        prefix="fig_gg2_closed_reference_unmapped_read_fraction_by_bioproject",
    )

    log(f"samples: {df.sample_key.nunique():,}")
    log(f"BioProjects: {df.bioproject.nunique():,}")
    log(f"median mapped fraction: {df.mapped_fraction.median():.4f}")
    log(f"wrote: {out_table}")
    log(f"figures: {outdir}")


if __name__ == "__main__":
    main()
