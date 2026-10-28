#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
s16_beta_pair_vs_random_from_distance_matrix_sampled.py

Integrated focal-only paired-vs-random beta-diversity analysis using
precomputed distance matrices from:

  03.beta_diversity/distance_matrix/

This version integrates the logic of:
  - s13 random unpaired Bray-Curtis
  - s15 paired/random UniFrac
  - s16 focal-only paired-vs-random summary/plotting

Design:
  focal_site  = Oral
  random_site = Gut

For each Oral focal sample:
  paired:
    Oral_i vs matched Gut_i within the same BioProject and same pair_id

  random-unpaired:
    Repeat n_random times. Each iteration randomly samples one Gut_j from
    the same BioProject, same optional stratum, and different pair_id.

Distance calculation:
  This script does not recompute Bray-Curtis or UniFrac. It only looks up
  precomputed sample-sample distances from the distance matrices.

Default output:
  03.beta_diversity/pair_vs_random/
"""

import argparse
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from scipy.stats import wilcoxon


DEFAULT_ROOT = Path("/beegfs/zhangcaiyi/oral_gut_pair_TimeSeries_11bp")
DEFAULT_DISTANCE_DIR = DEFAULT_ROOT / "03.beta_diversity" / "distance_matrix"
DEFAULT_METADATA = DEFAULT_ROOT / "01.filter_data" / "filter_metadata.csv"
DEFAULT_OUT_DIR = DEFAULT_ROOT / "03.beta_diversity" / "pair_vs_random"


METRIC_SEED_OFFSET = {
    "braycurtis": 11,
    "weighted_unifrac": 23,
    "unweighted_unifrac": 37,
}


# -----------------------------
# General helpers
# -----------------------------

def parse_list(x):
    if x is None:
        return []
    x = str(x).strip()
    if x == "":
        return []
    return [i.strip() for i in x.split(",") if i.strip()]


def clean_for_filename(x):
    x = str(x)
    x = re.sub(r"[^\w.\-]+", "_", x)
    x = x.strip("_")
    return x if x else "NA"


def read_table(path):
    path = str(path)
    if path.endswith((".tsv", ".txt")):
        return pd.read_csv(path, sep="\t", dtype=str)
    return pd.read_csv(path, sep=None, engine="python", dtype=str)


def safe_unique_join(x):
    vals = []
    for v in x:
        if pd.isna(v):
            continue
        v = str(v).strip()
        if v == "" or v.lower() in {"nan", "na", "none", "null"}:
            continue
        vals.append(v)
    vals = sorted(set(vals))
    return ";".join(vals) if vals else "NA"


def safe_bp_list(x):
    vals = []
    for v in x:
        if pd.isna(v):
            continue
        v = str(v).strip()
        if v:
            vals.append(v)
    vals = sorted(set(vals))
    return ";".join(vals) if vals else "NA"


def json_safe_records(df):
    tmp = df.copy().replace({np.nan: None})
    records = tmp.to_dict(orient="records")

    def convert(v):
        if isinstance(v, np.integer):
            return int(v)
        if isinstance(v, np.floating):
            if np.isnan(v):
                return None
            return float(v)
        if isinstance(v, np.bool_):
            return bool(v)
        return v

    return [{k: convert(v) for k, v in rec.items()} for rec in records]


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


def format_p(p):
    try:
        p = float(p)
    except Exception:
        return "NA"
    if not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4f}"


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    q = np.full_like(p, np.nan, dtype=float)

    valid = np.isfinite(p)
    if valid.sum() == 0:
        return q

    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    m = len(ranked)

    adjusted = ranked * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)

    q_valid = np.empty_like(adjusted)
    q_valid[order] = adjusted
    q[valid] = q_valid
    return q


def classify_direction(diff):
    if pd.isna(diff):
        return "NA"
    diff = float(diff)
    if diff < 0:
        return "paired_less_random"
    if diff > 0:
        return "paired_greater_random"
    return "paired_equal_random"


def standardize_site(x):
    x0 = str(x).strip()
    xl = x0.lower()

    site_map = {
        "oral": "Oral",
        "mouth": "Oral",
        "saliva": "Oral",
        "tongue": "Oral",
        "buccal": "Oral",
        "dental": "Oral",
        "plaque": "Oral",
        "gut": "Gut",
        "stool": "Gut",
        "fecal": "Gut",
        "faecal": "Gut",
        "feces": "Gut",
        "faeces": "Gut",
    }

    return site_map.get(xl, x0)


def normalize_status_value(x):
    if pd.isna(x):
        return "NA"
    x = str(x).strip()
    if x == "" or x.lower() in {"nan", "na", "none", "null"}:
        return "NA"
    return x


def is_hc_status(x):
    s = normalize_status_value(x).lower()
    return s in {
        "hc",
        "h",
        "healthy",
        "health",
        "healthy control",
        "healthy_control",
        "control",
        "controls",
        "normal",
        "non-disease",
        "non_disease",
        "non disease",
        "non-diseased",
        "non_diseased",
        "non diseased",
        "noncase",
        "non-case",
        "non_case",
    }


def status_sort_key(x):
    x = normalize_status_value(x)
    if x == "NA":
        return (2, x)
    if is_hc_status(x):
        return (0, x.lower())
    return (1, x.lower())


def classify_bp_status_group(status_values):
    values = [normalize_status_value(x) for x in status_values]
    values = [x for x in values if x != "NA"]

    if not values:
        return 2, "disease_study"

    has_hc = any(is_hc_status(x) for x in values)
    has_disease = any(not is_hc_status(x) for x in values)

    if has_hc and not has_disease:
        return 0, "HC"
    if has_hc and has_disease:
        return 1, "HC+disease"
    return 2, "disease_study"


def metric_display_name(metric):
    mapping = {
        "braycurtis": "Bray-Curtis",
        "weighted_unifrac": "Weighted UniFrac",
        "unweighted_unifrac": "Unweighted UniFrac",
    }
    return mapping.get(metric, metric)


def metric_to_dataset_label(metric):
    if metric == "braycurtis":
        return "filter_ra"
    if metric in {"weighted_unifrac", "unweighted_unifrac"}:
        return "filter"
    return "unknown"


def dataset_display_name(dataset_label):
    mapping = {
        "filter_ra": "Filtered relative abundance",
        "filter": "Filtered count table",
    }
    return mapping.get(dataset_label, dataset_label)


# -----------------------------
# Input loading
# -----------------------------

def read_distance_matrix(path):
    path = Path(path)
    df = pd.read_csv(path, sep="\t", index_col=0)

    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)

    common = [sid for sid in df.index if sid in set(df.columns)]
    df = df.loc[common, common].copy()

    df = df.apply(pd.to_numeric, errors="coerce")

    if df.isna().values.any():
        n_na = int(df.isna().sum().sum())
        raise RuntimeError(f"Distance matrix contains {n_na} NA values: {path}")

    arr = df.values.astype(float)

    if not np.allclose(arr, arr.T, atol=1e-8, rtol=1e-8):
        max_diff = float(np.max(np.abs(arr - arr.T)))
        print(f"WARNING: distance matrix is not perfectly symmetric; max_diff={max_diff}. Symmetrizing: {path}")
        arr = (arr + arr.T) / 2.0

    np.fill_diagonal(arr, 0.0)
    df = pd.DataFrame(arr, index=df.index, columns=df.columns)

    return df


def prepare_metadata(
    metadata,
    sample_col,
    bioproject_col,
    site_col,
    pair_col,
    subject_col,
    status_col,
    disease_col,
    stratify_cols,
):
    md = metadata.copy()

    required = [sample_col, bioproject_col, site_col, pair_col]
    missing = [c for c in required if c not in md.columns]
    if missing:
        raise ValueError(
            "Metadata missing required columns:\n  "
            + "\n  ".join(missing)
            + "\n\nAvailable columns:\n  "
            + "\n  ".join(md.columns)
        )

    for c in [subject_col, status_col, disease_col] + stratify_cols:
        if c not in md.columns:
            md[c] = "NA"

    for c in [sample_col, bioproject_col, site_col, pair_col, subject_col, status_col, disease_col] + stratify_cols:
        md[c] = md[c].astype(str).str.strip()

    md = md.drop_duplicates(subset=[sample_col]).copy()

    md["sample_id_for_analysis"] = md[sample_col]
    md["bioproject_for_analysis"] = md[bioproject_col]
    md["site_for_analysis"] = md[site_col].apply(standardize_site)
    md["pair_id_for_analysis"] = md[pair_col]
    md["subject_id_for_analysis"] = md[subject_col]
    md["status_for_analysis"] = md[status_col].map(normalize_status_value)
    md["disease_for_analysis"] = md[disease_col].map(normalize_status_value)

    md = md.loc[md["site_for_analysis"].isin(["Oral", "Gut"]), :].copy()

    if stratify_cols:
        md["stratum_key_for_random"] = md.apply(
            lambda row: "||".join([normalize_status_value(row.get(c, "NA")) for c in stratify_cols]),
            axis=1,
        )
    else:
        md["stratum_key_for_random"] = "ALL"

    return md


def build_bp_order(md):
    rows = []

    for bp, sub in md.groupby("bioproject_for_analysis", dropna=False):
        statuses = (
            sub["status_for_analysis"]
            .dropna()
            .astype(str)
            .map(normalize_status_value)
            .unique()
            .tolist()
        )

        statuses_clean = [s for s in statuses if s != "NA"]
        status_text = "NA" if not statuses_clean else "+".join(sorted(statuses_clean, key=status_sort_key))

        status_rank, status_group = classify_bp_status_group(statuses)

        rows.append({
            "bioproject": str(bp),
            "bp_status_label": f"{bp}({status_text})",
            "status_values": status_text,
            "status_rank": int(status_rank),
            "status_group_for_order": status_group,
        })

    out = pd.DataFrame(rows).sort_values(["status_rank", "bioproject"], ascending=[True, True])
    out["bp_plot_order"] = np.arange(out.shape[0], dtype=int)
    return out


def add_bp_annotations(df, bp_order):
    out = df.copy()
    if out.empty:
        return out

    out = out.merge(
        bp_order[["bioproject", "bp_status_label", "status_values", "status_group_for_order", "bp_plot_order"]],
        on="bioproject",
        how="left",
    )

    out["bp_status_label"] = out["bp_status_label"].fillna(out["bioproject"])
    out["status_values"] = out["status_values"].fillna("NA")
    out["status_group_for_order"] = out["status_group_for_order"].fillna("unknown")
    out["bp_plot_order"] = pd.to_numeric(out["bp_plot_order"], errors="coerce")

    missing = out["bp_plot_order"].isna()
    if missing.any():
        extra_bps = sorted(out.loc[missing, "bioproject"].dropna().unique())
        extra_map = {bp: 10000 + i for i, bp in enumerate(extra_bps)}
        out.loc[missing, "bp_plot_order"] = out.loc[missing, "bioproject"].map(extra_map)

    out["bp_plot_order"] = out["bp_plot_order"].astype(int)
    return out


# -----------------------------
# Paired-vs-random construction
# -----------------------------

def build_random_iterations_for_metric(
    dm_df,
    md,
    metric,
    n_random,
    seed,
    focal_site="Oral",
    random_site="Gut",
    exclude_same_subject=True,
):
    dataset_label = metric_to_dataset_label(metric)
    rng = np.random.default_rng(seed + METRIC_SEED_OFFSET.get(metric, 0))

    sample_ids_in_dm = set(dm_df.index)
    md2 = md.loc[md["sample_id_for_analysis"].isin(sample_ids_in_dm), :].copy()

    rows = []
    skipped_rows = []

    for bp, sub in md2.groupby("bioproject_for_analysis", dropna=False):
        bp = str(bp)
        sub = sub.copy()

        focal_df = sub.loc[sub["site_for_analysis"].eq(focal_site), :].copy()
        random_df = sub.loc[sub["site_for_analysis"].eq(random_site), :].copy()

        if focal_df.empty or random_df.empty:
            skipped_rows.append({
                "metric": metric,
                "bioproject": bp,
                "reason": "missing_focal_or_random_site",
                "n_focal_site": int(focal_df.shape[0]),
                "n_random_site": int(random_df.shape[0]),
            })
            continue

        random_by_pair = {
            str(pair_id): g["sample_id_for_analysis"].astype(str).tolist()
            for pair_id, g in random_df.groupby("pair_id_for_analysis", dropna=False)
        }

        candidate_index = {
            stratum: g.copy()
            for stratum, g in random_df.groupby("stratum_key_for_random", dropna=False)
        }

        for _, focal in focal_df.iterrows():
            focal_sample = str(focal["sample_id_for_analysis"])
            pair_id = str(focal["pair_id_for_analysis"])
            stratum = str(focal["stratum_key_for_random"])
            focal_subject = str(focal.get("subject_id_for_analysis", "NA"))

            if focal_sample not in dm_df.index:
                skipped_rows.append({
                    "metric": metric,
                    "bioproject": bp,
                    "reason": "focal_not_in_distance_matrix",
                    "focal_sample_id": focal_sample,
                    "pair_id": pair_id,
                })
                continue

            matched_random_samples = [
                sid for sid in random_by_pair.get(pair_id, [])
                if sid in dm_df.columns
            ]

            if len(matched_random_samples) == 0:
                skipped_rows.append({
                    "metric": metric,
                    "bioproject": bp,
                    "reason": "no_matched_random_site_same_pair_id",
                    "focal_sample_id": focal_sample,
                    "pair_id": pair_id,
                })
                continue

            candidates = candidate_index.get(stratum, pd.DataFrame()).copy()
            if candidates.empty:
                skipped_rows.append({
                    "metric": metric,
                    "bioproject": bp,
                    "reason": "no_candidate_same_stratum",
                    "focal_sample_id": focal_sample,
                    "pair_id": pair_id,
                    "stratum_key": stratum,
                })
                continue

            candidates = candidates.loc[candidates["pair_id_for_analysis"].astype(str) != pair_id, :].copy()
            if exclude_same_subject and focal_subject not in {"", "NA", "nan", "None"}:
                candidates = candidates.loc[
                    candidates["subject_id_for_analysis"].astype(str) != focal_subject, :
                ].copy()
            candidates = candidates.loc[candidates["sample_id_for_analysis"].isin(dm_df.columns), :].copy()

            candidate_pool_size = int(candidates.shape[0])
            if candidate_pool_size == 0:
                skipped_rows.append({
                    "metric": metric,
                    "bioproject": bp,
                    "reason": "no_candidate_after_excluding_same_pair_or_subject",
                    "focal_sample_id": focal_sample,
                    "pair_id": pair_id,
                    "stratum_key": stratum,
                })
                continue

            paired_vals = dm_df.loc[focal_sample, matched_random_samples].astype(float).values
            paired_vals = paired_vals[np.isfinite(paired_vals)]

            if len(paired_vals) == 0:
                skipped_rows.append({
                    "metric": metric,
                    "bioproject": bp,
                    "reason": "nonfinite_paired_distance",
                    "focal_sample_id": focal_sample,
                    "pair_id": pair_id,
                })
                continue

            paired_distance = float(np.mean(paired_vals))
            candidate_rows = candidates.reset_index(drop=True)

            for r in range(n_random):
                pick_idx = int(rng.integers(0, candidate_pool_size))
                picked = candidate_rows.iloc[pick_idx]

                random_sample = str(picked["sample_id_for_analysis"])
                random_pair_id = str(picked["pair_id_for_analysis"])
                random_distance = float(dm_df.loc[focal_sample, random_sample])

                rows.append({
                    "bioproject": bp,
                    "dataset_label": dataset_label,
                    "metric": metric,
                    "random_iter": int(r + 1),
                    "pair_id": pair_id,
                    "focal_sample_id": focal_sample,
                    "focal_site": focal_site,
                    "random_sample_id": random_sample,
                    "random_site": random_site,
                    "random_pair_id": random_pair_id,
                    "same_bioproject": True,
                    "same_status_or_stratum": bool(stratum != "ALL"),
                    "stratum_key": stratum,
                    "candidate_pool_size": candidate_pool_size,
                    "random_note": "ok",
                    "paired_distance": paired_distance,
                    "paired_distance_median_if_multiple": float(np.median(paired_vals)),
                    "paired_distance_min_if_multiple": float(np.min(paired_vals)),
                    "paired_distance_max_if_multiple": float(np.max(paired_vals)),
                    "n_matched_random_site_samples": int(len(paired_vals)),
                    "matched_random_site_sample_ids": ";".join(sorted(set(matched_random_samples))),
                    "random_unpaired_distance_iter": random_distance,
                    "subject_id": str(focal.get("subject_id_for_analysis", "NA")),
                    "focal_status_cleaned": str(focal.get("status_for_analysis", "NA")),
                    "focal_disease_study_cleaned": str(focal.get("disease_for_analysis", "NA")),
                })

    random_iter_df = pd.DataFrame(rows)
    skipped_df = pd.DataFrame(skipped_rows)
    return random_iter_df, skipped_df


def aggregate_random_iterations(random_iter_df):
    if random_iter_df.empty:
        return pd.DataFrame()

    group_cols = [
        "bioproject",
        "dataset_label",
        "metric",
        "pair_id",
        "focal_sample_id",
        "focal_site",
        "random_site",
        "subject_id",
        "focal_status_cleaned",
        "focal_disease_study_cleaned",
        "stratum_key",
    ]

    matched = (
        random_iter_df.groupby(group_cols, dropna=False)
        .agg(
            paired_distance=("paired_distance", "first"),
            paired_distance_median_if_multiple=("paired_distance_median_if_multiple", "first"),
            paired_distance_min_if_multiple=("paired_distance_min_if_multiple", "first"),
            paired_distance_max_if_multiple=("paired_distance_max_if_multiple", "first"),
            n_matched_random_site_samples=("n_matched_random_site_samples", "first"),
            matched_random_site_sample_ids=("matched_random_site_sample_ids", "first"),
            random_unpaired_distance=("random_unpaired_distance_iter", "median"),
            random_unpaired_distance_mean=("random_unpaired_distance_iter", "mean"),
            random_unpaired_distance_min=("random_unpaired_distance_iter", "min"),
            random_unpaired_distance_max=("random_unpaired_distance_iter", "max"),
            random_unpaired_distance_sd=("random_unpaired_distance_iter", "std"),
            n_random_distances=("random_unpaired_distance_iter", "size"),
            n_unique_random_samples=("random_sample_id", "nunique"),
            n_unique_random_pair_ids=("random_pair_id", "nunique"),
            random_sample_ids=("random_sample_id", lambda x: ";".join(sorted(set(map(str, x))))),
            random_pair_ids=("random_pair_id", lambda x: ";".join(sorted(set(map(str, x))))),
            candidate_pool_size_median=("candidate_pool_size", "median"),
        )
        .reset_index()
    )

    matched["random_unpaired_distance_sd"] = matched["random_unpaired_distance_sd"].fillna(0.0)
    matched["difference_paired_minus_random"] = matched["paired_distance"] - matched["random_unpaired_distance"]
    matched["paired_less_than_random"] = matched["paired_distance"] < matched["random_unpaired_distance"]
    matched["unit_direction"] = matched["difference_paired_minus_random"].apply(classify_direction)

    return matched


# -----------------------------
# Statistics and summaries
# -----------------------------

def wilcoxon_paired_test(paired_values, random_values):
    x = pd.to_numeric(pd.Series(paired_values), errors="coerce")
    y = pd.to_numeric(pd.Series(random_values), errors="coerce")

    valid = x.notna() & y.notna()
    x = x.loc[valid].astype(float).values
    y = y.loc[valid].astype(float).values
    n = len(x)

    if n == 0:
        return {
            "n_focal_units": 0,
            "statistic": np.nan,
            "p_value": np.nan,
            "p_significance_label": "NA",
            "median_paired": np.nan,
            "median_random_unpaired": np.nan,
            "median_difference_paired_minus_random": np.nan,
            "mean_difference_paired_minus_random": np.nan,
            "n_paired_less_random": 0,
            "n_paired_greater_random": 0,
            "frac_paired_less_random": np.nan,
            "test_note": "not_tested_no_data",
        }

    diff = x - y

    base = {
        "n_focal_units": int(n),
        "median_paired": float(np.nanmedian(x)),
        "median_random_unpaired": float(np.nanmedian(y)),
        "median_difference_paired_minus_random": float(np.nanmedian(diff)),
        "mean_difference_paired_minus_random": float(np.nanmean(diff)),
        "n_paired_less_random": int((x < y).sum()),
        "n_paired_greater_random": int((x > y).sum()),
        "frac_paired_less_random": float((x < y).mean()),
    }

    if n < 3:
        base.update({
            "statistic": np.nan,
            "p_value": np.nan,
            "p_significance_label": "NA",
            "test_note": "not_tested_n_less_than_3",
        })
        return base

    if np.allclose(diff, 0):
        base.update({
            "statistic": 0.0,
            "p_value": 1.0,
            "p_significance_label": "ns",
            "test_note": "all_differences_zero",
        })
        return base

    try:
        res = wilcoxon(
            x,
            y,
            alternative="two-sided",
            zero_method="wilcox",
            correction=False,
            mode="auto",
        )
        p = float(res.pvalue)
        base.update({
            "statistic": float(res.statistic),
            "p_value": p,
            "p_significance_label": p_to_stars(p),
            "test_note": "ok",
        })
        return base
    except Exception as e:
        base.update({
            "statistic": np.nan,
            "p_value": np.nan,
            "p_significance_label": "NA",
            "test_note": f"wilcoxon_failed:{type(e).__name__}:{str(e)}",
        })
        return base


def calculate_stats(matched):
    rows = []

    for (dataset_label, metric, bp), sub in matched.groupby(
        ["dataset_label", "metric", "bioproject"],
        dropna=False,
    ):
        res = wilcoxon_paired_test(sub["paired_distance"], sub["random_unpaired_distance"])
        row = {
            "scope": "per_bioproject",
            "dataset_label": dataset_label,
            "metric": metric,
            "bioproject": bp,
        }
        row.update(res)
        rows.append(row)

    for (dataset_label, metric), sub in matched.groupby(
        ["dataset_label", "metric"],
        dropna=False,
    ):
        res = wilcoxon_paired_test(sub["paired_distance"], sub["random_unpaired_distance"])
        row = {
            "scope": "pooled_all_bioprojects",
            "dataset_label": dataset_label,
            "metric": metric,
            "bioproject": "ALL",
        }
        row.update(res)
        rows.append(row)

    stat_df = pd.DataFrame(rows)
    stat_df["q_value_BH_within_dataset_metric_scope"] = np.nan

    for (dataset_label, metric, scope), idx in stat_df.groupby(
        ["dataset_label", "metric", "scope"],
        dropna=False,
    ).groups.items():
        idx = list(idx)
        stat_df.loc[idx, "q_value_BH_within_dataset_metric_scope"] = bh_fdr(
            stat_df.loc[idx, "p_value"].values
        )

    stat_df["direction_by_median_diff"] = stat_df["median_difference_paired_minus_random"].apply(classify_direction)
    return stat_df


def build_long_table(matched):
    rows = []

    for _, row in matched.iterrows():
        base = {
            "bioproject": row["bioproject"],
            "dataset_label": row["dataset_label"],
            "metric": row["metric"],
            "pair_id": row["pair_id"],
            "focal_sample_id": row["focal_sample_id"],
            "subject_id": row.get("subject_id", "NA"),
            "focal_status_cleaned": row.get("focal_status_cleaned", "NA"),
            "focal_disease_study_cleaned": row.get("focal_disease_study_cleaned", "NA"),
            "bp_plot_order": row.get("bp_plot_order", np.nan),
            "bp_status_label": row.get("bp_status_label", row["bioproject"]),
            "status_group_for_order": row.get("status_group_for_order", "unknown"),
        }

        paired = base.copy()
        paired["distance_type"] = "Paired"
        paired["distance"] = row["paired_distance"]
        rows.append(paired)

        random = base.copy()
        random["distance_type"] = "Random-unpaired"
        random["distance"] = row["random_unpaired_distance"]
        rows.append(random)

    return pd.DataFrame(rows)


def summarize_bp_level(bp_stats_annot):
    d = bp_stats_annot.copy()

    return (
        d.groupby(
            ["dataset_label", "metric", "direction_by_median_diff", "status_group_for_order"],
            dropna=False,
        )
        .agg(
            n_bioprojects=("bioproject", "nunique"),
            total_focal_units=("n_focal_units", "sum"),
            median_focal_units_per_bp=("n_focal_units", "median"),
            status_values=("status_values", safe_unique_join),
            bioprojects=("bioproject", safe_bp_list),
        )
        .reset_index()
    )


def summarize_unit_level(unit_annot):
    d = unit_annot.copy()

    by_status = (
        d.groupby(
            [
                "dataset_label",
                "metric",
                "unit_direction",
                "focal_status_cleaned",
                "focal_disease_study_cleaned",
                "status_group_for_order",
            ],
            dropna=False,
        )
        .agg(
            n_focal_units=("focal_sample_id", "size"),
            n_bioprojects=("bioproject", "nunique"),
            n_pair_ids=("pair_id", "nunique"),
            median_paired_distance=("paired_distance", "median"),
            median_random_unpaired_distance=("random_unpaired_distance", "median"),
            median_difference_paired_minus_random=("difference_paired_minus_random", "median"),
            mean_difference_paired_minus_random=("difference_paired_minus_random", "mean"),
            bioprojects=("bioproject", safe_bp_list),
        )
        .reset_index()
    )

    by_bp_status = (
        d.groupby(
            [
                "dataset_label",
                "metric",
                "bioproject",
                "unit_direction",
                "focal_status_cleaned",
                "focal_disease_study_cleaned",
            ],
            dropna=False,
        )
        .agg(
            n_focal_units=("focal_sample_id", "size"),
            n_pair_ids=("pair_id", "nunique"),
            median_paired_distance=("paired_distance", "median"),
            median_random_unpaired_distance=("random_unpaired_distance", "median"),
            median_difference_paired_minus_random=("difference_paired_minus_random", "median"),
            mean_difference_paired_minus_random=("difference_paired_minus_random", "mean"),
        )
        .reset_index()
    )

    return by_status, by_bp_status


# -----------------------------
# Plotting
# -----------------------------

def plot_per_bp_boxplot(
    long_df,
    stat_df,
    dataset_label,
    metric,
    out_png,
    width_per_bp=0.72,
    height=6.4,
    x_rotation=30,
    point_alpha=0.35,
):
    d = long_df.loc[
        (long_df["dataset_label"] == dataset_label)
        & (long_df["metric"] == metric),
        :,
    ].copy()

    if d.empty:
        print(f"WARNING: no data for {dataset_label} / {metric}")
        return None

    bp_info = (
        d[["bioproject", "bp_plot_order", "bp_status_label"]]
        .drop_duplicates()
        .sort_values(["bp_plot_order", "bioproject"])
        .copy()
    )

    bp_list = bp_info["bioproject"].tolist()
    bp_labels = bp_info["bp_status_label"].tolist()
    n_bp = len(bp_list)

    fig_width = max(10, width_per_bp * n_bp)
    fig, ax = plt.subplots(figsize=(fig_width, height))

    paired_color = "#D65F4B"
    random_color = "#63B6C7"
    rng = np.random.default_rng(2026)

    box_data = []
    positions = []
    colors = []

    for i, bp in enumerate(bp_list):
        base_x = i
        sub = d.loc[d["bioproject"] == bp, :].copy()

        paired = sub.loc[sub["distance_type"] == "Paired", "distance"].astype(float).values
        random = sub.loc[sub["distance_type"] == "Random-unpaired", "distance"].astype(float).values

        pos_paired = base_x - 0.18
        pos_random = base_x + 0.18

        if len(paired) > 0:
            box_data.append(paired)
            positions.append(pos_paired)
            colors.append(paired_color)
            xj = pos_paired + rng.normal(0, 0.025, size=len(paired))
            ax.scatter(xj, paired, s=8, color=paired_color, alpha=point_alpha, edgecolors="none", zorder=3)

        if len(random) > 0:
            box_data.append(random)
            positions.append(pos_random)
            colors.append(random_color)
            xj = pos_random + rng.normal(0, 0.025, size=len(random))
            ax.scatter(xj, random, s=8, color=random_color, alpha=point_alpha, edgecolors="none", zorder=3)

    if box_data:
        box = ax.boxplot(
            box_data,
            positions=positions,
            widths=0.28,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "black", "linewidth": 1.0},
            boxprops={"linewidth": 0.8, "color": "black"},
            whiskerprops={"linewidth": 0.8, "color": "black"},
            capprops={"linewidth": 0.8, "color": "black"},
        )
        for patch, color in zip(box["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.85)

    stat_sub = stat_df.loc[
        (stat_df["scope"] == "per_bioproject")
        & (stat_df["dataset_label"] == dataset_label)
        & (stat_df["metric"] == metric),
        :,
    ].copy()

    y_min, y_max = ax.get_ylim()
    y_range = y_max - y_min if y_max > y_min else 1.0

    for i, bp in enumerate(bp_list):
        row_df = stat_sub.loc[stat_sub["bioproject"] == bp, :]
        if row_df.empty:
            continue
        row = row_df.iloc[0]
        stars = str(row.get("p_significance_label", "NA"))
        sub_y = d.loc[d["bioproject"] == bp, "distance"].astype(float)
        if sub_y.empty:
            continue
        local_y = float(sub_y.max()) + 0.035 * y_range
        ax.text(i, local_y, stars, ha="center", va="bottom", fontsize=9, fontweight="bold" if stars not in {"ns", "NA"} else "normal")

    ax.set_xticks(np.arange(n_bp))
    ax.set_xticklabels(bp_labels, rotation=x_rotation, ha="right", fontsize=8)
    ax.set_ylabel("Distance")
    ax.set_xlabel("BioProject(status_cleaned)")
    ax.set_title(
        f"{metric_display_name(metric)}: Paired vs random-unpaired\n"
        f"{dataset_display_name(dataset_label)}; focal Oral comparison",
        fontsize=12,
    )
    ax.grid(axis="y", linestyle=":", linewidth=0.5, alpha=0.55)
    ax.legend(
        handles=[
            Patch(facecolor=paired_color, edgecolor="black", label="Paired"),
            Patch(facecolor=random_color, edgecolor="black", label="Random-unpaired"),
        ],
        frameon=False,
        loc="best",
    )

    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)
    return str(out_png)


def plot_pooled_boxplot(long_df, stat_df, dataset_label, metric, out_png, point_alpha=0.35):
    d = long_df.loc[
        (long_df["dataset_label"] == dataset_label)
        & (long_df["metric"] == metric),
        :,
    ].copy()

    if d.empty:
        print(f"WARNING: no pooled data for {dataset_label} / {metric}")
        return None

    paired = d.loc[d["distance_type"] == "Paired", "distance"].astype(float).values
    random = d.loc[d["distance_type"] == "Random-unpaired", "distance"].astype(float).values

    stat_row = stat_df.loc[
        (stat_df["scope"] == "pooled_all_bioprojects")
        & (stat_df["dataset_label"] == dataset_label)
        & (stat_df["metric"] == metric),
        :,
    ]

    if stat_row.empty:
        p_text = "NA"
        stars = "NA"
        n_units = len(paired)
        median_diff = np.nan
    else:
        row = stat_row.iloc[0]
        p_text = format_p(row["p_value"])
        stars = str(row["p_significance_label"])
        n_units = int(row["n_focal_units"])
        median_diff = float(row["median_difference_paired_minus_random"])

    paired_color = "#D65F4B"
    random_color = "#63B6C7"

    fig, ax = plt.subplots(figsize=(5.2, 6.2))
    box = ax.boxplot(
        [paired, random],
        positions=[0, 1],
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.0},
        boxprops={"linewidth": 0.8, "color": "black"},
        whiskerprops={"linewidth": 0.8, "color": "black"},
        capprops={"linewidth": 0.8, "color": "black"},
    )

    for patch, color in zip(box["boxes"], [paired_color, random_color]):
        patch.set_facecolor(color)
        patch.set_alpha(0.85)

    rng = np.random.default_rng(2026)
    xj = 0 + rng.normal(0, 0.045, size=len(paired))
    ax.scatter(xj, paired, s=8, color=paired_color, alpha=point_alpha, edgecolors="none", zorder=3)
    xj = 1 + rng.normal(0, 0.045, size=len(random))
    ax.scatter(xj, random, s=8, color=random_color, alpha=point_alpha, edgecolors="none", zorder=3)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Paired\nsame pair_id", "Random-unpaired\nmedian of 100 Gut"], fontsize=9)
    ax.set_ylabel("Distance")
    ax.set_title(
        f"{metric_display_name(metric)} distance\n"
        f"{dataset_display_name(dataset_label)}; focal Oral\n"
        f"Wilcoxon P={p_text} ({stars}); n={n_units}; median diff={median_diff:.4g}",
        fontsize=11,
    )
    ax.grid(axis="y", linestyle=":", linewidth=0.5, alpha=0.55)

    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)
    return str(out_png)


# -----------------------------
# Main
# -----------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Integrated focal-only paired-vs-random beta diversity from precomputed distance matrices."
    )

    parser.add_argument("--distance-dir", default=str(DEFAULT_DISTANCE_DIR))
    parser.add_argument("--metadata", default=str(DEFAULT_METADATA))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))

    parser.add_argument("--braycurtis-matrix", default=None)
    parser.add_argument("--weighted-unifrac-matrix", default=None)
    parser.add_argument("--unweighted-unifrac-matrix", default=None)

    parser.add_argument("--sample-col", default="sample_id")
    parser.add_argument("--bioproject-col", default="bioproject")
    parser.add_argument("--site-col", default="site_group_cleaned")
    parser.add_argument("--pair-col", default="pair_id")
    parser.add_argument("--subject-col", default="subject_id")
    parser.add_argument("--status-col", default="status_cleaned")
    parser.add_argument("--disease-col", default="disease_study_cleaned")

    parser.add_argument("--stratify-cols", default="status_cleaned")
    parser.add_argument("--n-random", type=int, default=100)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument(
        "--keep-same-subject",
        action="store_true",
        help="If set, random-unpaired candidates may include other timepoints of "
             "the SAME subject. Default excludes same-subject samples (correct for "
             "longitudinal/time-series data where a subject has multiple pair_ids).",
    )

    parser.add_argument("--metrics", default="weighted_unifrac,unweighted_unifrac")
    parser.add_argument("--dataset-labels", default="filter")

    parser.add_argument("--focal-site", default="Oral", choices=["Oral"])
    parser.add_argument("--random-site", default="Gut", choices=["Gut"])

    parser.add_argument("--width-per-bp", type=float, default=0.72)
    parser.add_argument("--per-bp-height", type=float, default=6.4)
    parser.add_argument("--x-rotation", type=float, default=30)
    parser.add_argument("--point-alpha", type=float, default=0.35)

    parser.add_argument("--overwrite", action="store_true")

    args = parser.parse_args()

    if args.n_random < 1:
        raise ValueError("--n-random must be >= 1")

    distance_dir = Path(args.distance_dir).resolve()
    metadata_fp = Path(args.metadata).resolve()
    out_dir = Path(args.out_dir).resolve()
    plot_dir = out_dir / "plots"

    out_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    metrics = parse_list(args.metrics)
    dataset_labels = parse_list(args.dataset_labels)
    stratify_cols = parse_list(args.stratify_cols)

    supported_metrics = {"braycurtis", "weighted_unifrac", "unweighted_unifrac"}
    bad = sorted(set(metrics) - supported_metrics)
    if bad:
        raise ValueError(f"Unsupported metrics: {bad}. Supported: {sorted(supported_metrics)}")

    matrix_paths = {
        "braycurtis": Path(args.braycurtis_matrix).resolve() if args.braycurtis_matrix else distance_dir / "braycurtis_filter_ra_distance_matrix.tsv",
        "weighted_unifrac": Path(args.weighted_unifrac_matrix).resolve() if args.weighted_unifrac_matrix else distance_dir / "weighted_unifrac_filter_distance_matrix.tsv",
        "unweighted_unifrac": Path(args.unweighted_unifrac_matrix).resolve() if args.unweighted_unifrac_matrix else distance_dir / "unweighted_unifrac_filter_distance_matrix.tsv",
    }

    outputs = {
        "random_iter_table": out_dir / "beta_pair_vs_random_random_iterations.tsv",
        "matched_table": out_dir / "beta_pair_vs_random_matched_pair_table.tsv",
        "long_table": out_dir / "beta_pair_vs_random_long_plot_table.tsv",
        "stats_table": out_dir / "beta_pair_vs_random_wilcoxon_stats.tsv",
        "bp_direction_detail": out_dir / "beta_direction_by_bp_detail.tsv",
        "bp_direction_summary": out_dir / "beta_direction_by_bp_summary.tsv",
        "unit_direction_detail": out_dir / "beta_direction_by_oral_unit_detail.tsv",
        "unit_direction_summary_by_status": out_dir / "beta_direction_by_oral_unit_summary_by_status.tsv",
        "unit_direction_summary_by_bp_status": out_dir / "beta_direction_by_oral_unit_summary_by_bp_status.tsv",
        "skipped_table": out_dir / "beta_pair_vs_random_skipped_units.tsv",
        "bp_order": out_dir / "beta_bioproject_order_by_status.tsv",
        "summary_json": out_dir / "beta_pair_vs_random_summary.json",
    }

    for p in outputs.values():
        if p.exists() and not args.overwrite:
            raise FileExistsError(f"Output exists: {p}. Use --overwrite.")

    print("Resolved paths:", flush=True)
    print(f"  distance_dir: {distance_dir}", flush=True)
    print(f"  metadata:     {metadata_fp}", flush=True)
    print(f"  out_dir:      {out_dir}", flush=True)

    if not metadata_fp.exists():
        raise FileNotFoundError(f"Metadata not found: {metadata_fp}")

    for metric in metrics:
        if not matrix_paths[metric].exists():
            raise FileNotFoundError(f"Distance matrix for {metric} not found: {matrix_paths[metric]}")

    print("[1/7] Reading metadata ...", flush=True)
    raw_md = read_table(metadata_fp)
    md = prepare_metadata(
        metadata=raw_md,
        sample_col=args.sample_col,
        bioproject_col=args.bioproject_col,
        site_col=args.site_col,
        pair_col=args.pair_col,
        subject_col=args.subject_col,
        status_col=args.status_col,
        disease_col=args.disease_col,
        stratify_cols=stratify_cols,
    )

    bp_order = build_bp_order(md)
    bp_order.to_csv(outputs["bp_order"], sep="\t", index=False)

    print("[2/7] Reading distance matrices and sampling random-unpaired pairs ...", flush=True)
    random_iter_tables = []
    skipped_tables = []
    metric_input_summaries = {}

    for metric in metrics:
        print(f"  Metric: {metric}", flush=True)
        dm_df = read_distance_matrix(matrix_paths[metric])

        random_iter_metric, skipped_metric = build_random_iterations_for_metric(
            dm_df=dm_df,
            md=md,
            metric=metric,
            n_random=args.n_random,
            seed=args.seed,
            focal_site=args.focal_site,
            random_site=args.random_site,
            exclude_same_subject=not args.keep_same_subject,
        )

        random_iter_tables.append(random_iter_metric)
        skipped_tables.append(skipped_metric)

        metric_input_summaries[metric] = {
            "distance_matrix": str(matrix_paths[metric]),
            "dataset_label": metric_to_dataset_label(metric),
            "n_samples_in_distance_matrix": int(dm_df.shape[0]),
            "n_random_iteration_rows": int(random_iter_metric.shape[0]),
            "n_skipped_units": int(skipped_metric.shape[0]) if not skipped_metric.empty else 0,
        }

        del dm_df

    random_iter = pd.concat(random_iter_tables, axis=0, ignore_index=True) if random_iter_tables else pd.DataFrame()
    skipped = pd.concat(skipped_tables, axis=0, ignore_index=True) if skipped_tables else pd.DataFrame()

    if random_iter.empty:
        raise RuntimeError("No random iteration rows were produced.")

    random_iter = random_iter.loc[
        random_iter["dataset_label"].isin(dataset_labels)
        & random_iter["metric"].isin(metrics),
        :,
    ].copy()

    if random_iter.empty:
        raise RuntimeError("Random iteration table is empty after dataset/metric filtering.")

    random_iter = add_bp_annotations(random_iter, bp_order)

    print("[3/7] Aggregating random iterations to focal-unit table ...", flush=True)
    matched = aggregate_random_iterations(random_iter)
    matched = add_bp_annotations(matched, bp_order)

    print("[4/7] Calculating Wilcoxon statistics and summaries ...", flush=True)
    stat_df = calculate_stats(matched)
    stat_df = add_bp_annotations(stat_df, bp_order)

    long_df = build_long_table(matched)
    bp_stats = stat_df.loc[stat_df["scope"].eq("per_bioproject"), :].copy()
    bp_direction_summary = summarize_bp_level(bp_stats)
    unit_summary_by_status, unit_summary_by_bp_status = summarize_unit_level(matched)

    print("[5/7] Plotting ...", flush=True)
    output_plots = []

    for dataset_label in dataset_labels:
        for metric in metrics:
            d_check = long_df.loc[
                (long_df["dataset_label"] == dataset_label)
                & (long_df["metric"] == metric),
                :,
            ]
            if d_check.empty:
                continue

            per_bp_png = plot_dir / f"beta_per_bp_pair_vs_random_{clean_for_filename(dataset_label)}_{clean_for_filename(metric)}.png"
            pooled_png = plot_dir / f"beta_pooled_pair_vs_random_{clean_for_filename(dataset_label)}_{clean_for_filename(metric)}.png"

            for p in [per_bp_png, pooled_png]:
                if p.exists() and not args.overwrite:
                    raise FileExistsError(f"Output exists: {p}. Use --overwrite.")

            print(f"  Plotting {dataset_label} / {metric}", flush=True)

            out1 = plot_per_bp_boxplot(
                long_df=long_df,
                stat_df=stat_df,
                dataset_label=dataset_label,
                metric=metric,
                out_png=per_bp_png,
                width_per_bp=args.width_per_bp,
                height=args.per_bp_height,
                x_rotation=args.x_rotation,
                point_alpha=args.point_alpha,
            )
            if out1 is not None:
                output_plots.append(out1)

            out2 = plot_pooled_boxplot(
                long_df=long_df,
                stat_df=stat_df,
                dataset_label=dataset_label,
                metric=metric,
                out_png=pooled_png,
                point_alpha=args.point_alpha,
            )
            if out2 is not None:
                output_plots.append(out2)

    print("[6/7] Writing output tables ...", flush=True)
    random_iter.to_csv(outputs["random_iter_table"], sep="\t", index=False)
    matched.to_csv(outputs["matched_table"], sep="\t", index=False)
    long_df.to_csv(outputs["long_table"], sep="\t", index=False)
    stat_df.to_csv(outputs["stats_table"], sep="\t", index=False)
    bp_stats.to_csv(outputs["bp_direction_detail"], sep="\t", index=False)
    bp_direction_summary.to_csv(outputs["bp_direction_summary"], sep="\t", index=False)
    matched.to_csv(outputs["unit_direction_detail"], sep="\t", index=False)
    unit_summary_by_status.to_csv(outputs["unit_direction_summary_by_status"], sep="\t", index=False)
    unit_summary_by_bp_status.to_csv(outputs["unit_direction_summary_by_bp_status"], sep="\t", index=False)

    if skipped.empty:
        pd.DataFrame(columns=["metric", "bioproject", "reason"]).to_csv(outputs["skipped_table"], sep="\t", index=False)
    else:
        skipped.to_csv(outputs["skipped_table"], sep="\t", index=False)

    print("[7/7] Writing summary JSON ...", flush=True)

    quick_pooled = stat_df.loc[stat_df["scope"].eq("pooled_all_bioprojects"), :].copy()

    quick_bp_counts = (
        bp_stats.groupby(["dataset_label", "metric", "direction_by_median_diff"], dropna=False)
        .agg(
            n_bioprojects=("bioproject", "nunique"),
            total_focal_units=("n_focal_units", "sum"),
            bioprojects=("bioproject", safe_bp_list),
        )
        .reset_index()
    )

    quick_unit_counts = (
        matched.groupby(["dataset_label", "metric", "unit_direction"], dropna=False)
        .agg(
            n_focal_units=("focal_sample_id", "size"),
            n_bioprojects=("bioproject", "nunique"),
            n_pair_ids=("pair_id", "nunique"),
        )
        .reset_index()
    )

    with open(outputs["summary_json"], "w") as f:
        json.dump(
            {
                "inputs": {
                    "distance_dir": os.path.abspath(str(distance_dir)),
                    "distance_matrices": {m: str(matrix_paths[m]) for m in metrics},
                    "metadata": os.path.abspath(str(metadata_fp)),
                },
                "out_dir": os.path.abspath(str(out_dir)),
                "dataset_labels": dataset_labels,
                "metrics": metrics,
                "stratify_cols": stratify_cols,
                "n_random": int(args.n_random),
                "seed": int(args.seed),
                "exclude_same_subject_in_random": bool(not args.keep_same_subject),
                "metric_input_summaries": metric_input_summaries,
                "focal_only_definition": {
                    "focal_site": args.focal_site,
                    "random_site": args.random_site,
                    "paired": "Oral_i vs matched Gut_i within the same BioProject and pair_id.",
                    "random_unpaired": (
                        "For each Oral_i, sample one Gut_j n_random times from the same "
                        "BioProject, same stratum if stratified, and different pair_id"
                        + ("; also excluding other timepoints of the SAME subject (longitudinal-safe)."
                           if not args.keep_same_subject else "; same-subject timepoints allowed.")
                    ),
                    "distance_source": "Precomputed distance matrices; no Bray-Curtis or UniFrac distances are recomputed in this script.",
                },
                "statistical_test": {
                    "name": "Wilcoxon signed-rank test",
                    "paired_unit": "BioProject × pair_id × focal_sample_id",
                    "comparison": "paired_distance vs median random_unpaired_distance across n_random sampled Gut candidates",
                    "direction": {
                        "paired_less_random": "paired_distance - random_unpaired_distance < 0",
                        "paired_greater_random": "paired_distance - random_unpaired_distance > 0",
                    },
                },
                "outputs": {k: str(v) for k, v in outputs.items()},
                "plots": output_plots,
                "n_rows_random_iter": int(random_iter.shape[0]),
                "n_rows_matched": int(matched.shape[0]),
                "n_rows_long": int(long_df.shape[0]),
                "n_rows_skipped": int(skipped.shape[0]) if not skipped.empty else 0,
                "quick_pooled_stats": json_safe_records(quick_pooled),
                "quick_bp_counts": json_safe_records(quick_bp_counts),
                "quick_unit_counts": json_safe_records(quick_unit_counts),
            },
            f,
            indent=2,
        )

    print("\nDone.")
    print("Main outputs:")
    for k, v in outputs.items():
        print(f"  {k}: {v}")

    print("\nPlots:")
    for p in output_plots:
        print(f"  {p}")

    print("\nQuick pooled stats:")
    cols = [
        "dataset_label",
        "metric",
        "n_focal_units",
        "median_paired",
        "median_random_unpaired",
        "median_difference_paired_minus_random",
        "frac_paired_less_random",
        "p_value",
        "p_significance_label",
        "q_value_BH_within_dataset_metric_scope",
    ]
    cols = [c for c in cols if c in quick_pooled.columns]
    print(quick_pooled[cols].to_string(index=False))

    print("\nQuick BP-level direction counts:")
    print(quick_bp_counts.to_string(index=False))


if __name__ == "__main__":
    main()

