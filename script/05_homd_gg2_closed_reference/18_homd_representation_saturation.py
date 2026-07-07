#!/usr/bin/env python3
"""HOMD representation saturation analyses using bioproject accumulation."""

from __future__ import annotations

import argparse
import json
import os
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

DEFAULT_OUTDIR = "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/region_aware_pipeline"

sns.set_theme(style="whitegrid", context="paper", font_scale=1.15)
plt.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

STATUS_COLORS = {
    "All ASVs": "#4A4A4A",
    "HOMD captured": "#2F7D73",
    "HOMD uncaptured": "#B44E3A",
}


def log(msg: str) -> None:
    print(f"[saturation] {msg}", flush=True)


def split_projects(x):
    if pd.isna(x):
        return []
    return [v for v in str(x).split(";") if v and v != "NA" and v != "nan"]


def captured(identity, qcov, threshold):
    identity = pd.to_numeric(identity, errors="coerce").fillna(0)
    qcov = pd.to_numeric(qcov, errors="coerce").fillna(0)
    return (identity >= threshold) & (qcov >= 99)


def make_light_cache(outdir):
    src = os.path.join(outdir, "homd_coverage_by_derep_sequence.tsv")
    cache = os.path.join(outdir, "homd_saturation_input_light.tsv")
    if os.path.exists(cache) and os.path.getsize(cache) > 0:
        log(f"using existing light cache: {cache}")
        return cache
    cols = ["derep_id", "n_bioprojects", "total_reads", "prevalence_samples", "bioprojects", "identity", "query_cov", "has_homd_hit"]
    log("building light cache from homd_coverage_by_derep_sequence.tsv")
    first = True
    for chunk in pd.read_csv(src, sep="\t", usecols=cols, chunksize=200000):
        for c in ["identity", "query_cov", "total_reads", "prevalence_samples", "n_bioprojects"]:
            chunk[c] = pd.to_numeric(chunk[c], errors="coerce").fillna(0)
        chunk.to_csv(cache, sep="\t", index=False, mode="w" if first else "a", header=first)
        first = False
        log(f"cached {'' if first else ''}{chunk.shape[0]:,} rows chunk")
    return cache


def load_dataset(input_path, robust_path=None, threshold=99):
    if robust_path is None:
        cols = ["derep_id", "n_bioprojects", "total_reads", "prevalence_samples", "bioprojects", "identity", "query_cov"]
        df = pd.read_csv(input_path, sep="\t", usecols=cols)
        df["weight"] = 1.0
        df["subject_prevalence"] = np.nan
    else:
        cols = ["derep_id", "n_bioprojects", "total_reads", "subject_prevalence", "bioprojects", "identity", "query_cov", "gg2_qc_class"]
        df = pd.read_csv(robust_path, sep="\t", usecols=cols)
        df = df[df["gg2_qc_class"].eq("bacterial")].copy()
        df["prevalence_samples"] = df["subject_prevalence"]
        df["weight"] = pd.to_numeric(df["subject_prevalence"], errors="coerce").fillna(0)
    for c in ["n_bioprojects", "total_reads", "prevalence_samples", "identity", "query_cov"]:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)
    df["captured"] = captured(df["identity"], df["query_cov"], threshold)
    df["project_list"] = df["bioprojects"].map(split_projects)
    df = df[df["project_list"].map(len) > 0].copy()
    return df


def build_project_index(df):
    projects = sorted(set(p for vals in df["project_list"] for p in vals))
    project_to_idx = defaultdict(list)
    for i, vals in enumerate(df["project_list"]):
        for p in vals:
            project_to_idx[p].append(i)
    project_to_idx = {p: np.array(v, dtype=np.int64) for p, v in project_to_idx.items()}
    return projects, project_to_idx


def accumulate_once(order, project_to_idx, captured_arr, weights=None):
    n = len(captured_arr)
    seen = np.zeros(n, dtype=bool)
    out = []
    if weights is None:
        weights = np.ones(n, dtype=float)
    captured_arr = captured_arr.astype(bool)
    for step, project in enumerate(order, start=1):
        idx = project_to_idx.get(project)
        if idx is not None and len(idx):
            seen[idx] = True
        seen_cap = seen & captured_arr
        seen_uncap = seen & (~captured_arr)
        out.append((
            step,
            float(weights[seen].sum()),
            float(weights[seen_cap].sum()),
            float(weights[seen_uncap].sum()),
        ))
    return out


def run_permutations(df, n_perm=100, seed=17, use_weights=False):
    projects, project_to_idx = build_project_index(df)
    captured_arr = df["captured"].to_numpy(bool)
    weights = df["weight"].to_numpy(float) if use_weights else None
    rng = np.random.default_rng(seed)
    rows = []
    for perm in range(n_perm):
        order = np.array(projects, dtype=object)
        rng.shuffle(order)
        acc = accumulate_once(order, project_to_idx, captured_arr, weights=weights)
        for step, all_v, cap_v, uncap_v in acc:
            rows.append({"perm": perm, "step": step, "status": "All ASVs", "value": all_v})
            rows.append({"perm": perm, "step": step, "status": "HOMD captured", "value": cap_v})
            rows.append({"perm": perm, "step": step, "status": "HOMD uncaptured", "value": uncap_v})
        if (perm + 1) % 20 == 0:
            log(f"completed {perm + 1}/{n_perm} permutations")
    return pd.DataFrame(rows), len(projects)


def summarize_curves(curves):
    return (
        curves.groupby(["step", "status"], as_index=False)["value"]
        .agg(median="median", q025=lambda x: np.quantile(x, 0.025), q975=lambda x: np.quantile(x, 0.975))
    )


def plot_three_panel(all_summary, robust_count_summary, robust_weight_summary, threshold, figdir, suffix):
    fig, axes = plt.subplots(1, 3, figsize=(15.8, 4.6), sharex=False)
    panels = [
        (all_summary, "A. All dereplicated ASVs", "Cumulative unique ASVs"),
        (robust_count_summary, "B. Robust GG2-bacterial low-HOMD ASVs", "Cumulative ASV count"),
        (robust_weight_summary, "C. Robust ASV subject-prevalence burden", "Cumulative subject prevalence"),
    ]
    for ax, (df, title, ylabel) in zip(axes, panels):
        for status in ["All ASVs", "HOMD captured", "HOMD uncaptured"]:
            sub = df[df["status"].eq(status)].sort_values("step")
            if sub.empty:
                continue
            x = sub["step"].to_numpy(float)
            med = sub["median"].to_numpy(float)
            lo = sub["q025"].to_numpy(float)
            hi = sub["q975"].to_numpy(float)
            ax.plot(x, med, color=STATUS_COLORS[status], linewidth=2.0, label=status)
            ax.fill_between(x, lo, hi, color=STATUS_COLORS[status], alpha=0.16, linewidth=0)
        ax.set_title(title, loc="left", fontweight="bold")
        ax.set_xlabel("Accumulated bioprojects")
        ax.set_ylabel(ylabel)
        ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x, pos: f"{x/1e6:.1f}M" if x >= 1e6 else (f"{x/1e3:.0f}K" if x >= 1e3 else f"{x:.0f}")))
    axes[0].legend(frameon=False, loc="upper left")
    fig.suptitle(f"HOMD representation saturation (identity >= {threshold}%, query coverage >= 99%)", y=1.03, fontweight="bold")
    fig.tight_layout()
    out = os.path.join(figdir, f"fig11_homd_representation_saturation_{suffix}")
    fig.savefig(out + ".pdf", bbox_inches="tight")
    fig.savefig(out + ".png", bbox_inches="tight", dpi=300)
    plt.close(fig)
    return out + ".pdf", out + ".png"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--n-perm", type=int, default=100)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--thresholds", default="99,97")
    args = parser.parse_args()

    figdir = os.path.join(args.outdir, "figures")
    os.makedirs(figdir, exist_ok=True)
    light_path = make_light_cache(args.outdir)
    robust_path = os.path.join(args.outdir, "low_identity_lt97_subject_prevalence_filtered_gg2_bacterial_only.tsv")

    all_stats = {}
    for threshold in [float(x) for x in args.thresholds.split(",")]:
        suffix = f"id{str(threshold).replace('.', 'p')}_qcov99"
        log(f"loading all-ASV dataset for threshold {threshold}")
        all_df = load_dataset(light_path, threshold=threshold)
        log(f"all-ASV rows: {len(all_df):,}")
        all_curves, n_projects_all = run_permutations(all_df, n_perm=args.n_perm, seed=args.seed, use_weights=False)
        all_summary = summarize_curves(all_curves)
        all_curves_path = os.path.join(args.outdir, f"homd_saturation_all_asvs_{suffix}.permutation_curves.tsv")
        all_summary_path = os.path.join(args.outdir, f"homd_saturation_all_asvs_{suffix}.summary.tsv")
        all_curves.to_csv(all_curves_path, sep="\t", index=False)
        all_summary.to_csv(all_summary_path, sep="\t", index=False)

        log(f"loading robust bacterial dataset for threshold {threshold}")
        robust_df = load_dataset(light_path, robust_path=robust_path, threshold=threshold)
        log(f"robust bacterial rows: {len(robust_df):,}")
        robust_count_curves, n_projects_robust = run_permutations(robust_df, n_perm=args.n_perm, seed=args.seed, use_weights=False)
        robust_weight_curves, _ = run_permutations(robust_df, n_perm=args.n_perm, seed=args.seed, use_weights=True)
        robust_count_summary = summarize_curves(robust_count_curves)
        robust_weight_summary = summarize_curves(robust_weight_curves)
        robust_count_summary_path = os.path.join(args.outdir, f"homd_saturation_robust_bacterial_counts_{suffix}.summary.tsv")
        robust_weight_summary_path = os.path.join(args.outdir, f"homd_saturation_robust_bacterial_subject_prevalence_{suffix}.summary.tsv")
        robust_count_summary.to_csv(robust_count_summary_path, sep="\t", index=False)
        robust_weight_summary.to_csv(robust_weight_summary_path, sep="\t", index=False)

        pdf, png = plot_three_panel(all_summary, robust_count_summary, robust_weight_summary, threshold, figdir, suffix)
        all_stats[suffix] = {
            "threshold_identity": threshold,
            "qcov_threshold": 99,
            "n_permutations": args.n_perm,
            "all_asv_rows": int(len(all_df)),
            "all_asv_bioprojects": int(n_projects_all),
            "robust_bacterial_rows": int(len(robust_df)),
            "robust_bacterial_bioprojects": int(n_projects_robust),
            "outputs": {
                "all_asv_curves": all_curves_path,
                "all_asv_summary": all_summary_path,
                "robust_bacterial_count_summary": robust_count_summary_path,
                "robust_bacterial_subject_prevalence_summary": robust_weight_summary_path,
                "figure_pdf": pdf,
                "figure_png": png,
            },
        }

    stats_path = os.path.join(args.outdir, "homd_representation_saturation_stats.json")
    with open(stats_path, "w") as handle:
        json.dump(all_stats, handle, indent=2)
    log(json.dumps(all_stats, indent=2))


if __name__ == "__main__":
    main()
