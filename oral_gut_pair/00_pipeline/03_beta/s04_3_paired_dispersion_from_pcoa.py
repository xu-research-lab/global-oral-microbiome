#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
s04_3_paired_dispersion_from_pcoa.py

Paired beta-dispersion test from existing s04 PCoA coordinates.

Purpose:
  Supplement pair-blocked PERMANOVA by testing whether Oral and Gut differ
  in multivariate dispersion under the matched Oral-Gut design.

Input:
  pcoa_coordinates_pooled.tsv
  pcoa_coordinates_per_bp.tsv
from s04_beta_diversity_pcoa_permanova_paired_adonis2.py output directory.

Method:
  1. Use all available PCoA axes PC1...PCn.
  2. Calculate distance of each sample to its site-specific centroid.
     - pooled coordinates: centroid is computed within BioProject x site
       to avoid BioProject-driven dispersion inflation.
     - per-BioProject coordinates: centroid is computed within site.
  3. For each complete Oral-Gut pair:
       delta = Gut_distance_to_centroid - Oral_distance_to_centroid
  4. Test whether mean(delta) differs from zero using paired sign-flip
     permutation across pair-level deltas.

Outputs:
  sample_dispersion_to_centroid.tsv
  paired_dispersion_delta.tsv
  dispersion_test_summary.tsv
  dispersion_summary.json
"""

import argparse
import itertools
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd


def parse_list(x):
    if x is None:
        return []
    x = str(x).strip()
    if x == "":
        return []
    return [i.strip() for i in x.split(",") if i.strip()]


def p_to_stars(p):
    try:
        p = float(p)
    except Exception:
        return "NA"
    if not np.isfinite(p):
        return "NA"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def clean_for_json(x):
    if isinstance(x, dict):
        return {k: clean_for_json(v) for k, v in x.items()}
    if isinstance(x, list):
        return [clean_for_json(v) for v in x]
    if isinstance(x, tuple):
        return [clean_for_json(v) for v in x]
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.floating):
        if np.isnan(x):
            return None
        return float(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def numeric_pc_sort_key(col):
    m = re.match(r"^PC(\d+)$", str(col))
    if m:
        return int(m.group(1))
    return 10**9


def get_pc_columns(df):
    pcs = [c for c in df.columns if re.match(r"^PC\d+$", str(c))]
    pcs = sorted(pcs, key=numeric_pc_sort_key)
    if not pcs:
        raise ValueError("No PCoA axis columns found. Expected columns like PC1, PC2, ...")
    return pcs


def read_coords(path, label):
    path = Path(path)
    if not path.exists():
        print(f"WARNING: {label} not found: {path}")
        return pd.DataFrame()

    df = pd.read_csv(path, sep="\t", dtype=str)
    if df.empty:
        print(f"WARNING: {label} is empty: {path}")
        return df

    required = [
        "sample_id_for_analysis",
        "metric",
        "site_for_analysis",
        "pair_id_for_analysis",
        "pair_block_for_permutation",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{label} missing columns: {missing}")

    if "bioproject_for_analysis" not in df.columns:
        if "bioproject" in df.columns:
            df["bioproject_for_analysis"] = df["bioproject"]
        else:
            raise ValueError(f"{label} missing bioproject_for_analysis / bioproject column")

    pc_cols = get_pc_columns(df)
    for c in pc_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.dropna(subset=pc_cols, how="all").copy()
    for c in pc_cols:
        df[c] = df[c].fillna(0.0)

    df["site_for_analysis"] = df["site_for_analysis"].astype(str).str.strip()
    df = df.loc[df["site_for_analysis"].isin(["Oral", "Gut"]), :].copy()

    return df


def add_distances_to_centroid(df, scope):
    """
    scope:
      pooled_all_bioprojects:
        centroid group = metric + bioproject_for_analysis + site
      per_bioproject:
        centroid group = metric + bioproject_for_analysis + site
        because each BP has its own PCoA coordinates.
    """
    if df.empty:
        return df

    pc_cols = get_pc_columns(df)
    out = df.copy()

    group_cols = ["metric", "bioproject_for_analysis", "site_for_analysis"]

    centers = out.groupby(group_cols, dropna=False)[pc_cols].transform("mean")
    diff = out[pc_cols].astype(float).values - centers.astype(float).values
    out["distance_to_site_centroid"] = np.sqrt(np.sum(diff * diff, axis=1))

    out["dispersion_scope"] = scope
    if scope == "pooled_all_bioprojects":
        out["centroid_definition"] = "within_BioProject_x_site_centroid_in_pooled_PCoA_space"
        out["test_bioproject"] = "ALL"
    else:
        out["centroid_definition"] = "within_site_centroid_in_per_BioProject_PCoA_space"
        out["test_bioproject"] = out["bioproject_for_analysis"].astype(str)

    return out


def build_pair_delta(sample_disp):
    if sample_disp.empty:
        return pd.DataFrame()

    rows = []

    group_cols = [
        "dispersion_scope",
        "centroid_definition",
        "metric",
        "test_bioproject",
        "bioproject_for_analysis",
        "pair_id_for_analysis",
        "pair_block_for_permutation",
    ]

    for key, sub in sample_disp.groupby(group_cols, dropna=False):
        d = dict(zip(group_cols, key))

        oral = sub.loc[sub["site_for_analysis"] == "Oral", :].copy()
        gut = sub.loc[sub["site_for_analysis"] == "Gut", :].copy()

        if oral.empty or gut.empty:
            continue

        # Strict paired data should have one Oral and one Gut. If duplicated,
        # take mean distance and retain all sample IDs for auditing.
        oral_dist = pd.to_numeric(oral["distance_to_site_centroid"], errors="coerce").mean()
        gut_dist = pd.to_numeric(gut["distance_to_site_centroid"], errors="coerce").mean()

        if not np.isfinite(oral_dist) or not np.isfinite(gut_dist):
            continue

        row = {
            **d,
            "oral_sample_ids": ";".join(sorted(oral["sample_id_for_analysis"].astype(str).unique())),
            "gut_sample_ids": ";".join(sorted(gut["sample_id_for_analysis"].astype(str).unique())),
            "n_oral_samples_in_pair": int(oral["sample_id_for_analysis"].nunique()),
            "n_gut_samples_in_pair": int(gut["sample_id_for_analysis"].nunique()),
            "oral_distance_to_centroid": float(oral_dist),
            "gut_distance_to_centroid": float(gut_dist),
            "delta_gut_minus_oral": float(gut_dist - oral_dist),
        }
        rows.append(row)

    return pd.DataFrame(rows)


def signflip_pvalue(delta, nperm=9999, seed=2026, exact_max_pairs=16):
    delta = np.asarray(delta, dtype=float)
    delta = delta[np.isfinite(delta)]
    n = len(delta)

    if n == 0:
        return np.nan, np.nan, 0, "no_pairs"

    obs = float(np.mean(delta))

    if np.allclose(delta, 0):
        return obs, 1.0, n, "all_zero_delta"

    abs_obs = abs(obs)

    # Exact enumeration for small n.
    if n <= exact_max_pairs:
        n_total = 0
        n_extreme = 0
        for signs in itertools.product([-1.0, 1.0], repeat=n):
            stat = float(np.mean(delta * np.asarray(signs)))
            n_total += 1
            if abs(stat) >= abs_obs - 1e-15:
                n_extreme += 1
        p = n_extreme / n_total
        return obs, float(p), n, f"exact_signflip_{n_total}_permutations"

    rng = np.random.default_rng(seed)
    n_extreme = 0
    done = 0
    batch = 5000

    while done < nperm:
        b = min(batch, nperm - done)
        signs = rng.choice(np.array([-1.0, 1.0]), size=(b, n), replace=True)
        stats = (signs * delta[None, :]).mean(axis=1)
        n_extreme += int((np.abs(stats) >= abs_obs - 1e-15).sum())
        done += b

    # Add-one correction.
    p = (1 + n_extreme) / (1 + nperm)
    return obs, float(p), n, f"monte_carlo_signflip_{nperm}_permutations"


def wilcoxon_pvalue(delta):
    try:
        from scipy import stats
    except Exception:
        return np.nan, "scipy_not_available"

    delta = np.asarray(delta, dtype=float)
    delta = delta[np.isfinite(delta)]

    if len(delta) == 0:
        return np.nan, "no_pairs"
    if np.allclose(delta, 0):
        return 1.0, "all_zero_delta"

    try:
        res = stats.wilcoxon(delta, zero_method="wilcox", alternative="two-sided", mode="auto")
        return float(res.pvalue), "ok"
    except Exception as e:
        return np.nan, f"wilcoxon_failed:{type(e).__name__}:{str(e)}"


def summarize_tests(delta_df, permutations, seed):
    if delta_df.empty:
        return pd.DataFrame()

    rows = []

    # Pooled rows: test_bioproject == ALL
    # Per-BP rows: test_bioproject == actual BioProject
    group_cols = [
        "dispersion_scope",
        "centroid_definition",
        "metric",
        "test_bioproject",
    ]

    for key, sub in delta_df.groupby(group_cols, dropna=False):
        d = dict(zip(group_cols, key))

        delta = pd.to_numeric(sub["delta_gut_minus_oral"], errors="coerce").dropna().values
        oral = pd.to_numeric(sub["oral_distance_to_centroid"], errors="coerce")
        gut = pd.to_numeric(sub["gut_distance_to_centroid"], errors="coerce")

        obs_mean_delta, p_sf, n_pairs, sf_note = signflip_pvalue(
            delta,
            nperm=permutations,
            seed=seed,
        )
        p_wilcox, wilcox_note = wilcoxon_pvalue(delta)

        row = {
            **d,
            "n_pairs": int(n_pairs),
            "n_oral_samples": int(sub["oral_sample_ids"].nunique()),
            "n_gut_samples": int(sub["gut_sample_ids"].nunique()),
            "mean_oral_dispersion": float(oral.mean()) if len(oral) else np.nan,
            "mean_gut_dispersion": float(gut.mean()) if len(gut) else np.nan,
            "median_oral_dispersion": float(oral.median()) if len(oral) else np.nan,
            "median_gut_dispersion": float(gut.median()) if len(gut) else np.nan,
            "mean_delta_gut_minus_oral": float(np.mean(delta)) if len(delta) else np.nan,
            "median_delta_gut_minus_oral": float(np.median(delta)) if len(delta) else np.nan,
            "sd_delta_gut_minus_oral": float(np.std(delta, ddof=1)) if len(delta) > 1 else np.nan,
            "signflip_statistic": "mean_delta_gut_minus_oral",
            "p_signflip_two_sided": p_sf,
            "p_signflip_label": p_to_stars(p_sf),
            "signflip_note": sf_note,
            "p_wilcoxon_two_sided": p_wilcox,
            "p_wilcoxon_label": p_to_stars(p_wilcox),
            "wilcoxon_note": wilcox_note,
            "permutations_requested": int(permutations),
            "seed": int(seed),
            "interpretation": "positive_delta_means_Gut_more_dispersed_than_Oral",
        }
        rows.append(row)

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Calculate paired beta-dispersion test from existing s04 PCoA coordinates."
    )

    parser.add_argument("--pcoa-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--analysis-label", default="pair_blocked_main")
    parser.add_argument("--metrics", default="braycurtis,weighted_unifrac,unweighted_unifrac")
    parser.add_argument("--permutations", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--overwrite", action="store_true")

    args = parser.parse_args()

    pcoa_dir = Path(args.pcoa_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    metrics = set(parse_list(args.metrics))

    pooled_fp = pcoa_dir / "pcoa_coordinates_pooled.tsv"
    per_bp_fp = pcoa_dir / "pcoa_coordinates_per_bp.tsv"

    outputs = {
        "sample_dispersion": out_dir / "sample_dispersion_to_centroid.tsv",
        "paired_delta": out_dir / "paired_dispersion_delta.tsv",
        "test_summary": out_dir / "dispersion_test_summary.tsv",
        "summary_json": out_dir / "dispersion_summary.json",
    }

    for p in outputs.values():
        if p.exists() and not args.overwrite:
            raise FileExistsError(f"Output exists: {p}. Use --overwrite.")

    print("Resolved paths:", flush=True)
    print(f"  pcoa_dir: {pcoa_dir}", flush=True)
    print(f"  pooled_coords: {pooled_fp}", flush=True)
    print(f"  per_bp_coords: {per_bp_fp}", flush=True)
    print(f"  out_dir: {out_dir}", flush=True)
    print(f"  analysis_label: {args.analysis_label}", flush=True)

    pooled = read_coords(pooled_fp, "pooled PCoA coordinates")
    per_bp = read_coords(per_bp_fp, "per-BioProject PCoA coordinates")

    if not pooled.empty:
        pooled = pooled.loc[pooled["metric"].isin(metrics), :].copy()
        pooled = add_distances_to_centroid(pooled, scope="pooled_all_bioprojects")

    if not per_bp.empty:
        per_bp = per_bp.loc[per_bp["metric"].isin(metrics), :].copy()
        per_bp = add_distances_to_centroid(per_bp, scope="per_bioproject")

    sample_disp = pd.concat([pooled, per_bp], axis=0, ignore_index=True)

    if sample_disp.empty:
        raise RuntimeError("No PCoA coordinate rows available after filtering metrics.")

    sample_disp["analysis_label"] = args.analysis_label

    delta_df = build_pair_delta(sample_disp)
    delta_df["analysis_label"] = args.analysis_label

    test_summary = summarize_tests(
        delta_df,
        permutations=args.permutations,
        seed=args.seed,
    )
    test_summary["analysis_label"] = args.analysis_label

    # Put analysis_label first.
    for df_name, df in [("sample_disp", sample_disp), ("delta_df", delta_df), ("test_summary", test_summary)]:
        if "analysis_label" in df.columns:
            cols = ["analysis_label"] + [c for c in df.columns if c != "analysis_label"]
            if df_name == "sample_disp":
                sample_disp = df[cols]
            elif df_name == "delta_df":
                delta_df = df[cols]
            else:
                test_summary = df[cols]

    sample_disp.to_csv(outputs["sample_dispersion"], sep="\t", index=False)
    delta_df.to_csv(outputs["paired_delta"], sep="\t", index=False)
    test_summary.to_csv(outputs["test_summary"], sep="\t", index=False)

    summary = {
        "pcoa_dir": str(pcoa_dir),
        "out_dir": str(out_dir),
        "analysis_label": args.analysis_label,
        "metrics": sorted(metrics),
        "permutations": int(args.permutations),
        "seed": int(args.seed),
        "method": {
            "dispersion_measure": "Euclidean distance to site centroid using all PC axes in existing PCoA coordinates",
            "pooled_centroid": "within BioProject x site centroid in pooled PCoA space",
            "per_bioproject_centroid": "within site centroid in per-BioProject PCoA space",
            "paired_test": "Gut-Oral distance-to-centroid delta tested by paired sign-flip permutation",
            "delta_definition": "Gut_distance_to_centroid - Oral_distance_to_centroid",
        },
        "outputs": {k: str(v) for k, v in outputs.items()},
        "n_sample_dispersion_rows": int(sample_disp.shape[0]),
        "n_paired_delta_rows": int(delta_df.shape[0]),
        "n_test_summary_rows": int(test_summary.shape[0]),
    }

    with open(outputs["summary_json"], "w") as f:
        json.dump(clean_for_json(summary), f, indent=2)

    print("\nDone.", flush=True)
    print("Main outputs:", flush=True)
    for k, v in outputs.items():
        print(f"  {k}: {v}", flush=True)

    print("\nTest summary:", flush=True)
    show_cols = [
        "analysis_label",
        "dispersion_scope",
        "metric",
        "test_bioproject",
        "n_pairs",
        "mean_oral_dispersion",
        "mean_gut_dispersion",
        "mean_delta_gut_minus_oral",
        "p_signflip_two_sided",
        "p_signflip_label",
        "p_wilcoxon_two_sided",
        "p_wilcoxon_label",
    ]
    show_cols = [c for c in show_cols if c in test_summary.columns]
    print(test_summary[show_cols].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
