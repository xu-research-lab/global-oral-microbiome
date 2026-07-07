#!/usr/bin/env python3
"""Plot rank-decay curves for GG2 feature reads and sample prevalence."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_PROFILE_DIR = (
    "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/"
    "containment_derep/unified_gg2_closed_reference_profile"
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile-dir", default=DEFAULT_PROFILE_DIR)
    args = parser.parse_args()

    profile_dir = Path(args.profile_dir)
    in_path = profile_dir / "all_samples_gg2_closed_reference_feature_summary.tsv"
    fig_dir = profile_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(in_path, sep="\t")
    df["closed_reference_reads"] = pd.to_numeric(df["closed_reference_reads"], errors="coerce").fillna(0)
    df["sample_prevalence"] = pd.to_numeric(df["sample_prevalence"], errors="coerce").fillna(0)
    df = df.sort_values(["closed_reference_reads", "sample_prevalence"], ascending=False).reset_index(drop=True)
    df["rank"] = np.arange(1, len(df) + 1)

    out_table = profile_dir / "all_samples_gg2_closed_reference_feature_rank_decay.tsv"
    df[["rank", "gg2_target", "closed_reference_reads", "sample_prevalence"]].to_csv(out_table, sep="\t", index=False)

    rank = df["rank"].to_numpy()
    reads = df["closed_reference_reads"].to_numpy()
    prevalence = df["sample_prevalence"].to_numpy()

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.8), sharex=False)

    axes[0].plot(rank, reads, color="#1F5A7A", linewidth=1.5)
    axes[0].set_yscale("log")
    axes[0].set_xlabel("GG2 target rank by reads")
    axes[0].set_ylabel("Closed-reference reads, log scale")
    axes[0].set_title("A. Read Abundance Decay")
    axes[0].grid(True, which="major", axis="both", linewidth=0.35, alpha=0.35)
    axes[0].grid(True, which="minor", axis="y", linewidth=0.2, alpha=0.18)

    prev_df = df.sort_values(["sample_prevalence", "closed_reference_reads"], ascending=False).reset_index(drop=True)
    prev_rank = np.arange(1, len(prev_df) + 1)
    prev_values = prev_df["sample_prevalence"].to_numpy()
    axes[1].plot(prev_rank, prev_values, color="#C45A11", linewidth=1.5)
    axes[1].set_yscale("log")
    axes[1].set_xlabel("GG2 target rank by sample prevalence")
    axes[1].set_ylabel("Sample prevalence, log scale")
    axes[1].set_title("B. Sample Prevalence Decay")
    axes[1].grid(True, which="major", axis="both", linewidth=0.35, alpha=0.35)
    axes[1].grid(True, which="minor", axis="y", linewidth=0.2, alpha=0.18)

    fig.suptitle("Rank-Decay of GG2 Closed-Reference Features", y=1.02, fontsize=12)
    fig.tight_layout()
    out_prefix = fig_dir / "fig_gg2_closed_reference_feature_rank_decay_two_panel"
    fig.savefig(f"{out_prefix}.png", dpi=300)
    fig.savefig(f"{out_prefix}.pdf")
    plt.close(fig)

    print(f"features\t{len(df)}")
    print(f"table\t{out_table}")
    print(f"figure_png\t{out_prefix}.png")
    print(f"figure_pdf\t{out_prefix}.pdf")


if __name__ == "__main__":
    main()
