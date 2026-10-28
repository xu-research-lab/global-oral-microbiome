#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot Gut and Oral alpha diversity for PRJNA767939."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path("/home/hexincheng/projects/hxc/oral_project/oral_gut_pair")
ALPHA_PATH = ROOT / "results/single_pair/02_alpha/delta_alpha.tsv"
META_PATH = ROOT / "results/cross_comparison/01_merged_metrics/pair_metrics_all_cohorts.tsv"
OUT_DIR = ROOT / "results/new/figures/08_delta_alpha_by_bp_status/PRJNA767939"
BIOPROJECT = "PRJNA767939"
GROUP_ORDER = ["HC", "COVID-19"]
GROUP_LABELS = {"HC": "Healthy", "COVID-19": "COVID-19"}
GROUP_COLORS = {"HC": "#2ca25f", "COVID-19": "#de2d26"}
METRICS = [
    ("Observed OTUs", "Observed"),
    ("Shannon", "Shannon"),
    ("Simpson (InvSimpson)", "InvSimpson"),
]


def load_data() -> pd.DataFrame:
    alpha = pd.read_csv(ALPHA_PATH, sep="\t")
    alpha = alpha.loc[
        alpha["dataset_label"].eq("filter_ra")
        & alpha["bioproject"].astype(str).eq(BIOPROJECT)
    ].copy()

    meta = pd.read_csv(META_PATH, sep="\t")
    meta = meta.loc[
        meta["analysis_set"].astype(str).eq("single_pair")
        & meta["bioproject"].astype(str).eq(BIOPROJECT),
        ["pair_id", "status_cleaned"],
    ].drop_duplicates("pair_id")
    data = alpha.merge(meta, on="pair_id", how="left", validate="one_to_one")
    data = data.loc[data["status_cleaned"].isin(GROUP_ORDER)].copy()
    if data.empty:
        raise ValueError(f"No alpha data found for {BIOPROJECT}")

    for site in ("oral", "gut"):
        for _, metric in METRICS:
            col = f"{site}_{metric}"
            if col not in data.columns:
                raise ValueError(f"Missing alpha metric column: {col}")
            data[col] = pd.to_numeric(data[col], errors="coerce")
    counts = data["status_cleaned"].value_counts()
    missing_groups = [group for group in GROUP_ORDER if counts.get(group, 0) == 0]
    if missing_groups:
        raise ValueError(f"Missing groups for {BIOPROJECT}: {missing_groups}")
    return data


def plot_alpha(data: pd.DataFrame) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(8.4, 5.2))
    rng = np.random.default_rng(2026)
    panels = [
        ("Gut", "gut", "gut"),
        ("Oral", "oral", "oral"),
    ]
    panel_letters = iter("ABCDEF")

    for row, (site_label, site, x_suffix) in enumerate(panels):
        for col, (metric_label, metric) in enumerate(METRICS):
            ax = axes[row, col]
            groups = [
                data.loc[data["status_cleaned"].eq(group), f"{site}_{metric}"]
                .dropna()
                .to_numpy(dtype=float)
                for group in GROUP_ORDER
            ]
            bp = ax.boxplot(
                groups,
                positions=[0, 1],
                widths=0.48,
                patch_artist=True,
                showfliers=False,
                medianprops={"color": "#222222", "linewidth": 1.2},
                whiskerprops={"color": "#444444", "linewidth": 0.8},
                capprops={"color": "#444444", "linewidth": 0.8},
                boxprops={"edgecolor": "#444444", "linewidth": 0.8},
            )
            for patch, group in zip(bp["boxes"], GROUP_ORDER):
                patch.set_facecolor(GROUP_COLORS[group])
                patch.set_alpha(0.28)

            for pos, vals, group in zip([0, 1], groups, GROUP_ORDER):
                if len(vals):
                    ax.scatter(
                        pos + rng.uniform(-0.075, 0.075, len(vals)),
                        vals,
                        s=13,
                        alpha=0.78,
                        color=GROUP_COLORS[group],
                        edgecolors="white",
                        linewidths=0.25,
                        zorder=3,
                    )

            ax.set_xticks([0, 1])
            ax.set_xticklabels(
                [f"Healthy_{x_suffix}", f"COVID_{x_suffix}"],
                fontsize=7,
                rotation=18,
                ha="right",
            )
            ax.set_ylabel(metric_label, fontsize=8)
            ax.set_title(f"{next(panel_letters)}. {site_label} · {metric_label}", fontsize=8.5)
            ax.tick_params(axis="y", labelsize=7, length=2)
            ax.tick_params(axis="x", length=0, pad=2)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

    legend_handles = [
        plt.Line2D(
            [], [], marker="o", linestyle="None", markersize=5,
            color=GROUP_COLORS[group], label=GROUP_LABELS[group],
        )
        for group in GROUP_ORDER
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        ncol=2,
        frameon=False,
        fontsize=8,
        bbox_to_anchor=(0.5, 1.01),
    )
    fig.suptitle(f"{BIOPROJECT} alpha diversity", y=1.06, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94), w_pad=0.65, h_pad=0.8)

    stem = OUT_DIR / f"{BIOPROJECT}_alpha_diversity_6panel"
    fig.savefig(stem.with_suffix(".png"), dpi=220, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {stem}.png and {stem}.pdf")
    for group in GROUP_ORDER:
        print(f"{GROUP_LABELS[group]} n={int(data['status_cleaned'].eq(group).sum())} pairs")


def main() -> None:
    plot_alpha(load_data())


if __name__ == "__main__":
    main()
