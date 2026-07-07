#!/usr/bin/env python3
"""Draw region UpSet plot for GG2 bacterial subset of strict low-HOMD-identity ASVs."""

from __future__ import annotations

import argparse
import json
import os
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

DEFAULT_OUTDIR = "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/region_aware_pipeline"


def split_values(x):
    if pd.isna(x):
        return []
    return [v for v in str(x).split(";") if v and v not in {"nan", "NA"}]


def plot_upset(set_counter, output_base, title, top_n=40):
    items = set_counter.most_common(top_n)
    if not items:
        raise RuntimeError("No region intersections to plot")
    combos = [k for k, v in items]
    counts = [v for k, v in items]
    regions = sorted(set().union(*[set(c) for c in combos]))

    fig = plt.figure(figsize=(max(10, top_n * 0.32), 6.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.5, 1.5], hspace=0.05)
    ax_bar = fig.add_subplot(gs[0])
    ax_mat = fig.add_subplot(gs[1], sharex=ax_bar)

    x = np.arange(len(combos))
    ax_bar.bar(x, counts, color="#2F5D62", width=0.82)
    ax_bar.set_ylabel("Bacterial ASV count")
    ax_bar.set_title(title)
    ax_bar.tick_params(axis="x", labelbottom=False)
    ax_bar.spines[["top", "right"]].set_visible(False)

    for i, combo in enumerate(combos):
        present = set(combo)
        ys = []
        for j, region in enumerate(regions):
            if region in present:
                ax_mat.scatter(i, j, s=36, color="#222222", zorder=3)
                ys.append(j)
            else:
                ax_mat.scatter(i, j, s=12, color="#D2D2D2", zorder=2)
        if len(ys) > 1:
            ax_mat.plot([i, i], [min(ys), max(ys)], color="#222222", linewidth=0.8, zorder=1)

    ax_mat.set_yticks(np.arange(len(regions)))
    ax_mat.set_yticklabels(regions)
    ax_mat.set_xlabel("Region intersection")
    ax_mat.set_xlim(-0.6, len(combos) - 0.4)
    ax_mat.set_xticks([])
    ax_mat.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(output_base + ".pdf", bbox_inches="tight")
    fig.savefig(output_base + ".png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--top-n", type=int, default=40)
    args = parser.parse_args()

    in_path = os.path.join(args.outdir, "low_identity_lt97_subject_prevalence_filtered_gg2_bacterial_only.tsv")
    out_prefix = os.path.join(args.outdir, "low_identity_lt97_subject_prevalence_filtered_gg2_bacterial_only")
    figdir = os.path.join(args.outdir, "figures")
    os.makedirs(figdir, exist_ok=True)
    fig_base = os.path.join(figdir, "fig8_low_identity_subject_prevalence_region_upset_gg2_bacterial_only")

    df = pd.read_csv(in_path, sep="\t", usecols=["derep_id", "subject_geo_regions", "subject_prevalence", "n_bioprojects", "total_reads", "gg2_identity", "gg2_phylum", "gg2_genus", "gg2_species"])
    region_counter = Counter(tuple(split_values(x)) for x in df["subject_geo_regions"])
    region_counter.pop(tuple(), None)

    inter = pd.DataFrame([
        {"regions": ";".join(k), "bacterial_asv_count": v}
        for k, v in region_counter.most_common()
    ])
    inter_path = out_prefix + ".region_intersections.tsv"
    inter.to_csv(inter_path, sep="\t", index=False)

    plot_upset(
        region_counter,
        fig_base,
        "GG2 bacterial low-HOMD-identity ASVs shared across geographic regions",
        top_n=args.top_n,
    )

    stats = {
        "input_bacterial_asvs": int(len(df)),
        "n_region_intersections": int(len(region_counter)),
        "top20_region_intersections": [
            {"regions": ";".join(k), "bacterial_asv_count": int(v)}
            for k, v in region_counter.most_common(20)
        ],
        "outputs": {
            "region_intersections": inter_path,
            "upset_pdf": fig_base + ".pdf",
            "upset_png": fig_base + ".png",
        },
    }
    stats_path = out_prefix + ".region_upset.stats.json"
    with open(stats_path, "w") as handle:
        json.dump(stats, handle, indent=2)
    print(json.dumps(stats, indent=2), flush=True)


if __name__ == "__main__":
    main()
