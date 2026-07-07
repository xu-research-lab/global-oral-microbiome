#!/usr/bin/env python3
"""Scatter plot: unique ASV HOMD best identity vs prevalence, colored by bioproject recurrence."""

from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

DEFAULT_OUTDIR = "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/region_aware_pipeline"

plt.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def log(msg):
    print(f"[scatter] {msg}", flush=True)


def load_plot_data(path):
    cols = ["derep_id", "identity", "prevalence_samples", "n_bioprojects", "total_reads", "n_asvs", "query_cov"]
    df = pd.read_csv(path, sep="\t", usecols=cols)
    for col in ["identity", "prevalence_samples", "n_bioprojects", "total_reads", "n_asvs", "query_cov"]:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)
    # No HOMD hit has identity=0 after fillna; keep it visible as the leftmost group.
    df["prevalence_plot"] = df["prevalence_samples"] + 1
    df["n_bioprojects_plot"] = df["n_bioprojects"].clip(lower=1)
    return df


def scatter(df, out_base, title_suffix="", max_points=None, hit_only=False, fixed_size=False):
    if hit_only:
        df = df[df["identity"] >= 80].copy()
    if max_points and len(df) > max_points:
        # Keep all recurrent/high-prevalence points and sample the dense singleton tail.
        important = df[(df["n_bioprojects"] >= 3) | (df["prevalence_samples"] >= 50) | (df["total_reads"] >= 1000)]
        rest = df.drop(index=important.index)
        n_rest = max(0, max_points - len(important))
        if len(rest) > n_rest:
            rest = rest.sample(n=n_rest, random_state=13)
        plot_df = pd.concat([important, rest], ignore_index=True)
    else:
        plot_df = df

    fig, ax = plt.subplots(figsize=(8.8, 5.8))
    norm = mpl.colors.LogNorm(vmin=1, vmax=max(2, float(df["n_bioprojects_plot"].max())))
    point_size = 8 if fixed_size else np.clip(np.sqrt(plot_df["n_asvs"].to_numpy()) * 4, 4, 45)
    sc = ax.scatter(
        plot_df["identity"],
        plot_df["prevalence_plot"],
        c=plot_df["n_bioprojects_plot"],
        s=point_size,
        cmap="viridis",
        norm=norm,
        alpha=0.35,
        linewidths=0,
        rasterized=True,
    )
    ax.axvline(99, color="#B22222", linestyle="--", linewidth=1.2, label="99% identity")
    ax.axvline(97, color="#555555", linestyle=":", linewidth=1.0, label="97% identity")
    ax.set_yscale("log")
    ax.set_xlim(80, 100.5) if hit_only else ax.set_xlim(-2, 101)
    ax.set_xlabel("Best HOMD identity (%)")
    ax.set_ylabel("ASV prevalence across samples + 1")
    ax.set_title("Unique ASV prevalence vs best HOMD identity" + title_suffix)
    cb = fig.colorbar(sc, ax=ax, pad=0.015)
    cb.set_label("Number of bioprojects sharing exact ASV")
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(out_base + ".png", bbox_inches="tight")
    fig.savefig(out_base + ".pdf", bbox_inches="tight")
    plt.close(fig)


def binned_summary(df, out_path):
    bins = [-0.1, 0.1, 80, 85, 90, 95, 97, 99, 99.5, 100.0001]
    labels = ["no_hit", "below_80", "80_85", "85_90", "90_95", "95_97", "97_99", "99_99.5", "99.5_100"]
    df = df.copy()
    df["identity_bin"] = pd.cut(df["identity"], bins=bins, labels=labels, include_lowest=True, right=False)
    summary = df.groupby("identity_bin", observed=False).agg(
        derep_sequences=("derep_id", "count"),
        source_asvs=("n_asvs", "sum"),
        total_reads=("total_reads", "sum"),
        prevalence_samples=("prevalence_samples", "sum"),
        median_prevalence=("prevalence_samples", "median"),
        max_prevalence=("prevalence_samples", "max"),
        median_bioprojects=("n_bioprojects", "median"),
        max_bioprojects=("n_bioprojects", "max"),
    ).reset_index()
    summary.to_csv(out_path, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--max-pdf-points", type=int, default=250000)
    args = parser.parse_args()

    in_path = os.path.join(args.outdir, "qcov_sensitivity_by_derep_sequence.tsv")
    figdir = os.path.join(args.outdir, "figures")
    os.makedirs(figdir, exist_ok=True)

    log("loading per-unique-ASV coverage table")
    df = load_plot_data(in_path)
    log(f"loaded {len(df):,} unique ASV/derep sequences")

    full_base = os.path.join(figdir, "fig5_unique_asv_identity_vs_prevalence_by_bioprojects_full")
    sampled_base = os.path.join(figdir, "fig5_unique_asv_identity_vs_prevalence_by_bioprojects_sampled")
    hit_full_base = os.path.join(figdir, "fig5_unique_asv_identity_vs_prevalence_by_bioprojects_hit_only_80to100_full")
    hit_sampled_base = os.path.join(figdir, "fig5_unique_asv_identity_vs_prevalence_by_bioprojects_hit_only_80to100_sampled")
    hit_fixed_base = os.path.join(figdir, "fig5_unique_asv_identity_vs_prevalence_by_bioprojects_hit_only_80to100_fixed_size")
    summary_path = os.path.join(args.outdir, "unique_asv_identity_prevalence_bins.tsv")

    log("writing binned summary")
    binned_summary(df, summary_path)
    log("plotting full rasterized scatter")
    scatter(df, full_base, title_suffix=" (all unique ASVs)", max_points=None)
    log("plotting sampled/recurrent scatter")
    scatter(df, sampled_base, title_suffix=" (recurrent/high-prevalence retained; tail sampled)", max_points=args.max_pdf_points)
    log("plotting hit-only 80-100 full scatter")
    scatter(df, hit_full_base, title_suffix=" (HOMD hit only, identity >=80%)", max_points=None, hit_only=True)
    log("plotting hit-only 80-100 sampled scatter")
    scatter(df, hit_sampled_base, title_suffix=" (HOMD hit only, identity >=80%; recurrent/high-prevalence retained)", max_points=args.max_pdf_points, hit_only=True)
    log("plotting hit-only 80-100 fixed-size scatter")
    scatter(df, hit_fixed_base, title_suffix=" (HOMD hit only, identity >=80%; fixed point size)", max_points=None, hit_only=True, fixed_size=True)
    log(f"done: {hit_fixed_base}.png/pdf")


if __name__ == "__main__":
    main()
