#!/usr/bin/env python3
"""Plot BioProject-level GG2 closed-reference read retention as a barplot."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


DEFAULT_PROFILE_DIR = (
    "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/"
    "containment_derep/unified_gg2_closed_reference_profile"
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile-dir", default=DEFAULT_PROFILE_DIR)
    parser.add_argument("--threshold", type=float, default=0.3)
    args = parser.parse_args()

    profile_dir = Path(args.profile_dir)
    fig_dir = profile_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(profile_dir / "all_projects_gg2_closed_reference_mapping_fraction_summary.tsv", sep="\t")
    df["overall_mapped_fraction"] = pd.to_numeric(df["overall_mapped_fraction"], errors="coerce")
    df = df.dropna(subset=["overall_mapped_fraction"]).sort_values("overall_mapped_fraction", ascending=False).copy()
    df["retention_status"] = df["overall_mapped_fraction"].ge(args.threshold).map({True: "Retained", False: "Dropped"})

    retained = int((df["overall_mapped_fraction"] >= args.threshold).sum())
    dropped = int((df["overall_mapped_fraction"] < args.threshold).sum())

    out_summary = profile_dir / "all_projects_gg2_closed_reference_retention_threshold30_summary.tsv"
    df.to_csv(out_summary, sep="\t", index=False)

    # Also compute sample-level stats for annotation
    sample_frac = pd.read_csv(profile_dir / "all_samples_gg2_closed_reference_mapping_fraction.tsv", sep="\t")
    sample_frac = sample_frac[sample_frac["bioproject"].isin(df[df["retention_status"] == "Retained"]["bioproject"])]
    samples_total = len(sample_frac)
    samples_retained = int((sample_frac["mapped_fraction"] >= args.threshold).sum())
    samples_dropped = samples_total - samples_retained

    sns.set_theme(style="whitegrid", context="paper")
    width = max(18, len(df) * 0.055)
    fig, ax = plt.subplots(figsize=(width, 5.8))
    colors = df["retention_status"].map({"Retained": "#2E86AB", "Dropped": "#A9A9A9"}).tolist()
    ax.bar(range(len(df)), df["overall_mapped_fraction"], color=colors, width=0.86, linewidth=0)
    ax.axhline(args.threshold, color="#D62728", linestyle="--", linewidth=1.2)
    ax.text(
        len(df) - 1,
        args.threshold + 0.025,
        f"30% threshold: studies retained {retained}, dropped {dropped} | "
        f"samples retained {samples_retained}, dropped {samples_dropped}",
        color="#D62728",
        ha="right",
        va="bottom",
        fontsize=9,
    )
    ax.set_ylim(0, 1.02)
    ax.set_xlim(-0.8, len(df) - 0.2)
    ax.set_xlabel("BioProject, ordered by dataset-level read retention")
    ax.set_ylabel("Dataset-level GG2 closed-reference read retention")
    ax.set_title("BioProject-Level Read Retention After GG2 Closed-Reference Mapping")
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df["bioproject"], rotation=90, fontsize=4.5)
    ax.yaxis.set_major_formatter(lambda x, pos: f"{x * 100:.0f}%")
    fig.tight_layout()

    out_prefix = fig_dir / "fig_gg2_closed_reference_project_read_retention_barplot"
    fig.savefig(f"{out_prefix}.png", dpi=300)
    fig.savefig(f"{out_prefix}.pdf")
    plt.close(fig)

    print(f"retained_projects\t{retained}")
    print(f"dropped_projects\t{dropped}")
    print(f"summary\t{out_summary}")
    print(f"figure_png\t{out_prefix}.png")
    print(f"figure_pdf\t{out_prefix}.pdf")


if __name__ == "__main__":
    main()
