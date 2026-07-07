#!/usr/bin/env python3
"""Plot top reference-gap HOMD references: identity distributions plus burden heatmap."""

from __future__ import annotations

import argparse
import os
import textwrap

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec

DEFAULT_OUTDIR = "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/region_aware_pipeline"

sns.set_theme(style="whitegrid", context="paper", font_scale=1.0)
plt.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def log(msg):
    print(f"[top-gap] {msg}", flush=True)


def load_hits(path):
    cols = [
        "derep_id", "homd_target", "identity", "alnlen", "mismatches", "gaps", "qlo", "qhi", "tlo", "thi",
        "query_length", "target_length", "query_cov", "target_cov",
    ]
    hits = pd.read_csv(path, sep="\t", names=cols)
    hits = hits.sort_values(["derep_id", "identity", "query_cov", "alnlen", "homd_target"], ascending=[True, False, False, False, True])
    hits = hits.drop_duplicates("derep_id", keep="first")
    return hits


def short_label(row):
    species = str(row.homd_species) if pd.notna(row.homd_species) else "NA"
    hmt = str(row.homd_hmt) if pd.notna(row.homd_hmt) else str(row.homd_target).split("_")[0]
    ref_suffix = str(row.homd_target).split("_", 1)[1] if "_" in str(row.homd_target) else str(row.homd_target)
    label = f"{hmt} {species} {ref_suffix}"
    return "\n".join(textwrap.wrap(label, width=26, break_long_words=False))


def robust_minmax(series):
    x = pd.to_numeric(series, errors="coerce").replace([np.inf, -np.inf], np.nan)
    return x


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--top-n", type=int, default=100)
    parser.add_argument("--min-derep", type=int, default=100)
    parser.add_argument("--min-bioprojects", type=int, default=10)
    parser.add_argument("--max-points-per-reference", type=int, default=350)
    args = parser.parse_args()

    ann_path = os.path.join(args.outdir, "homd_tree_tip_annotation.tsv")
    hit_path = os.path.join(args.outdir, "homd_vsearch_best_hits.tsv")
    figdir = os.path.join(args.outdir, "figures")
    os.makedirs(figdir, exist_ok=True)
    out_base = os.path.join(figdir, f"fig7_top{args.top_n}_reference_gap_identity_distribution_heatmap")
    selected_path = os.path.join(args.outdir, f"top{args.top_n}_reference_gap_selected.tsv")
    plot_data_path = os.path.join(args.outdir, f"top{args.top_n}_reference_gap_identity_plot_data.tsv")

    log("selecting top references")
    ann = pd.read_csv(ann_path, sep="\t")
    eligible = ann[(ann["any_hit_derep_sequences"] >= args.min_derep) & (ann["n_bioprojects_any"] >= args.min_bioprojects)].copy()
    selected = eligible.sort_values("reference_gap_score", ascending=False).head(args.top_n).copy()
    selected["plot_label"] = selected.apply(short_label, axis=1)
    selected["plot_order"] = np.arange(len(selected))
    selected.to_csv(selected_path, sep="\t", index=False)
    targets = set(selected["homd_target"])

    log("loading HOMD best hits")
    hits = load_hits(hit_path)
    hits = hits[hits["homd_target"].isin(targets)].copy()
    hits = hits.merge(selected[["homd_target", "plot_label", "plot_order"]], on="homd_target", how="left")
    hits["identity"] = pd.to_numeric(hits["identity"], errors="coerce")
    hits["query_cov"] = pd.to_numeric(hits["query_cov"], errors="coerce")
    hits.to_csv(plot_data_path, sep="\t", index=False)

    log("sampling jitter points")
    sampled = []
    for _, sub in hits.groupby("homd_target", sort=False):
        if len(sub) > args.max_points_per_reference:
            # Keep identity extremes and a deterministic random middle sample.
            extremes = pd.concat([sub.nsmallest(25, "identity"), sub.nlargest(25, "identity")]).drop_duplicates("derep_id")
            rest = sub.drop(index=extremes.index, errors="ignore")
            need = max(0, args.max_points_per_reference - len(extremes))
            if len(rest) > need:
                rest = rest.sample(n=need, random_state=17)
            sub = pd.concat([extremes, rest], ignore_index=True)
        sampled.append(sub)
    jitter = pd.concat(sampled, ignore_index=True)

    order = selected.sort_values("reference_gap_score", ascending=False)["plot_label"].tolist()
    selected = selected.set_index("plot_label").loc[order].reset_index()
    hits["plot_label"] = pd.Categorical(hits["plot_label"], categories=order, ordered=True)
    jitter["plot_label"] = pd.Categorical(jitter["plot_label"], categories=order, ordered=True)

    log("building heatmap matrix")
    heat = selected.copy()
    heat_metrics = pd.DataFrame({
        "log10 unique ASVs": np.log10(heat["any_hit_derep_sequences"] + 1),
        "Bioprojects": heat["n_bioprojects_any"],
        "Median identity": heat["any_hit_median_identity"],
        "Fraction id>=99%": heat["fraction_asv_id99_qcov99"],
        "Gap score": heat["reference_gap_score"],
    }, index=heat["plot_label"])
    heat_scaled = heat_metrics.copy()
    for col in heat_scaled.columns:
        vals = robust_minmax(heat_scaled[col])
        if col == "Fraction id>=99%":
            heat_scaled[col] = vals
        else:
            lo, hi = vals.min(), vals.max()
            heat_scaled[col] = 0 if hi == lo else (vals - lo) / (hi - lo)

    log("plotting")
    fig = plt.figure(figsize=(max(24, args.top_n * 0.28), 10.5))
    gs = GridSpec(2, 1, height_ratios=[4.3, 1.25], hspace=0.05, figure=fig)
    ax = fig.add_subplot(gs[0, 0])
    ax_h = fig.add_subplot(gs[1, 0], sharex=ax)

    sns.boxplot(
        data=hits,
        x="plot_label",
        y="identity",
        order=order,
        ax=ax,
        color="#D8E3EA",
        fliersize=0,
        linewidth=0.65,
    )
    sns.stripplot(
        data=jitter,
        x="plot_label",
        y="identity",
        order=order,
        ax=ax,
        hue=(jitter["identity"] >= 99).map({True: ">=99%", False: "<99%"}),
        palette={">=99%": "#1B9E77", "<99%": "#D95F02"},
        size=1.2,
        alpha=0.35,
        jitter=0.25,
        dodge=False,
        linewidth=0,
    )
    ax.axhline(99, color="#B22222", linestyle="--", linewidth=1.0)
    ax.axhline(97, color="#666666", linestyle=":", linewidth=0.9)
    ax.set_ylim(78, 100.5)
    ax.set_ylabel("Best HOMD identity (%)")
    ax.set_xlabel("")
    ax.set_title(f"Panel A. ASV identity distributions for top {args.top_n} reference-gap HOMD sequences")
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2], title="ASV identity", frameon=False, loc="lower right")
    ax.tick_params(axis="x", labelbottom=False)

    # Secondary line on Panel A for unique ASV burden, scaled to a right axis.
    x = np.arange(len(order))
    burden = selected["any_hit_derep_sequences"].to_numpy(dtype=float)
    burden_log = np.log10(burden + 1)
    ax2 = ax.twinx()
    ax2.plot(x, burden_log, color="#2C3E50", linewidth=1.15, alpha=0.9)
    ax2.set_ylabel("log10(best-mapped unique ASVs + 1)", color="#2C3E50")
    ax2.tick_params(axis="y", colors="#2C3E50")
    ax2.grid(False)

    sns.heatmap(
        heat_scaled.T,
        ax=ax_h,
        cmap="viridis",
        cbar_kws={"label": "Scaled metric value"},
        xticklabels=order,
        yticklabels=heat_scaled.columns,
        vmin=0,
        vmax=1,
    )
    ax_h.set_title("Panel B. Reference-level burden, recurrence and gap metrics", loc="left")
    ax_h.set_xlabel("HOMD reference, ranked by reference_gap_score")
    ax_h.set_ylabel("")
    ax_h.tick_params(axis="x", labelrotation=90, labelsize=5.5)
    ax_h.tick_params(axis="y", labelrotation=0)

    fig.savefig(out_base + ".pdf", bbox_inches="tight")
    fig.savefig(out_base + ".png", bbox_inches="tight", dpi=300)
    plt.close(fig)
    log(f"wrote {out_base}.pdf/.png")
    log(f"wrote {selected_path}")
    log(f"wrote {plot_data_path}")


if __name__ == "__main__":
    main()
