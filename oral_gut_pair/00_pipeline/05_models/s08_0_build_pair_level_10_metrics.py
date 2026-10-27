#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
s08_0_build_pair_level_10_metrics.py

Build pair-level LME input table by merging scattered outputs:

  1. Alpha delta:
       02.alpha_diversity/delta_alpha.tsv

  2. Paired beta distances:
       03.beta_diversity/distance_matrix/paired_beta_distances_all_pairs.tsv

  3. Shared OTU:
       04.sourcetracking/share_otu/shared_otu_oral_to_gut_per_sample.tsv

  4. SourceTracker2:
       04.sourcetracking/sourcetracker/paired_oral_fraction.tsv

  5. Pair-level metadata / sampling time:
       00.data/pair_sampling_time_from_pair_time_point.tsv
       00.data/metadata_time_series_11bp_with_sampling_time.csv
       or 01.filter_data/filter_metadata.csv

Output:
  05.lme_time_series_full_model/lme_input_10_metrics.tsv
  05.lme_time_series_full_model/lme_input_metric_availability.tsv
  05.lme_time_series_full_model/lme_input_merge_summary.json
"""

import argparse
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd


OUTCOME_COLUMNS = [
    "delta_alpha_observed",
    "delta_alpha_shannon",
    "delta_alpha_invsimpson",
    "beta_braycurtis_distance",
    "beta_weighted_unifrac_distance",
    "beta_unweighted_unifrac_distance",
    "shared_otu_count",
    "shared_otu_count_fraction_of_gut",
    "shared_gut_abundance_fraction",
    "sourcetracker_paired_oral_fraction",
]


def clean_str(x):
    if pd.isna(x):
        return np.nan
    s = str(x).strip()
    if s == "" or s.lower() in {"nan", "na", "none", "null"}:
        return np.nan
    return s


def read_table(path):
    path = Path(path)
    if not path.exists():
        print(f"WARNING: missing input file: {path}", flush=True)
        return pd.DataFrame()

    if str(path).endswith((".tsv", ".txt")):
        return pd.read_csv(path, sep="\t", dtype=str)
    return pd.read_csv(path, dtype=str)


def normalize_keys(df):
    if df.empty:
        return df

    # Flexible rename
    ren = {}
    for c in df.columns:
        lc = c.lower()
        if lc in {"bp", "bioproject_id", "pair_bioproject"} and "bioproject" not in df.columns:
            ren[c] = "bioproject"
        if lc in {"pair", "pairid", "pair_id_for_analysis"} and "pair_id" not in df.columns:
            ren[c] = "pair_id"

    if ren:
        df = df.rename(columns=ren)

    if "bioproject" not in df.columns or "pair_id" not in df.columns:
        raise ValueError(
            "Input table must contain bioproject and pair_id columns. "
            f"Available columns: {list(df.columns)}"
        )

    df["bioproject"] = df["bioproject"].map(clean_str)
    df["pair_id"] = df["pair_id"].map(clean_str)

    return df


def first_nonmissing_unique(x):
    vals = []
    for v in x:
        v = clean_str(v)
        if pd.notna(v):
            vals.append(str(v))
    vals = sorted(set(vals))
    if len(vals) == 0:
        return np.nan
    if len(vals) == 1:
        return vals[0]
    return ";".join(vals)


def make_pair_metadata(sample_metadata_fp, pair_sampling_fp):
    """
    Collapse sample-level metadata to one row per BioProject x pair_id,
    then merge pair-level sampling time table.
    """
    md = read_table(sample_metadata_fp)
    if md.empty:
        raise FileNotFoundError(f"Metadata is empty or missing: {sample_metadata_fp}")

    md = normalize_keys(md)

    # Keep useful covariates if present.
    preferred_cols = [
        "bioproject",
        "pair_id",
        "subject_id",
        "pair_time_point",
        "status_cleaned",
        "disease_study_cleaned",
        "age_cleaned",
        "age_cleaned_group",
        "bmi_cleaned",
        "bmi_cleaned_group",
        "sex_cleaned",
        "platform_cleaned",
        "sequencing_region_cleaned",
        "country_cleaned",
        "continent_group",
        "site_group_cleaned",
        "sample_id",
    ]
    keep_cols = [c for c in preferred_cols if c in md.columns]

    md_small = md[keep_cols].copy()

    # Collapse to pair-level.
    agg = {}
    for c in md_small.columns:
        if c in {"bioproject", "pair_id"}:
            continue
        agg[c] = first_nonmissing_unique

    pair_md = (
        md_small.groupby(["bioproject", "pair_id"], dropna=False)
        .agg(agg)
        .reset_index()
    )

    # Sampling-time table from s07.
    st = read_table(pair_sampling_fp)
    if not st.empty:
        st = normalize_keys(st)
        st_cols = [
            "bioproject",
            "pair_id",
            "subject_id",
            "pair_time_point",
            "sampling_time_numeric",
            "sampling_time_parse_type",
            "sampling_time_scale_type",
            "sampling_time_z_by_bp",
            "sampling_time_rank_by_bp",
            "sampling_time_rank_z_by_bp",
            "subject_block",
            "pair_block",
        ]
        st_cols = [c for c in st_cols if c in st.columns]
        st = st[st_cols].drop_duplicates(subset=["bioproject", "pair_id"]).copy()

        pair_md = pair_md.merge(
            st,
            on=["bioproject", "pair_id"],
            how="left",
            suffixes=("", "_from_sampling_time"),
        )

        # Prefer sampling-time subject_id/pair_time_point if original missing.
        for c in ["subject_id", "pair_time_point"]:
            alt = c + "_from_sampling_time"
            if alt in pair_md.columns:
                pair_md[c] = pair_md[c].where(pair_md[c].notna(), pair_md[alt])
                pair_md = pair_md.drop(columns=[alt])

    # Add blocks if missing.
    if "subject_block" not in pair_md.columns:
        pair_md["subject_block"] = (
            pair_md["bioproject"].astype(str) + "__" + pair_md.get("subject_id", "").astype(str)
        )
    if "pair_block" not in pair_md.columns:
        pair_md["pair_block"] = (
            pair_md["bioproject"].astype(str) + "__" + pair_md["pair_id"].astype(str)
        )

    return pair_md


def load_alpha_delta(alpha_fp):
    df = read_table(alpha_fp)
    if df.empty:
        return pd.DataFrame(columns=["bioproject", "pair_id"])

    df = normalize_keys(df)

    # Keep only complete pairs when column exists.
    if "has_complete_oral_gut_pair" in df.columns:
        keep = df["has_complete_oral_gut_pair"].astype(str).str.lower().isin({"true", "1", "yes", "y"})
        if keep.any():
            df = df.loc[keep].copy()

    out = df[["bioproject", "pair_id"]].drop_duplicates().copy()

    # Long format option.
    metric_col = None
    for c in ["metric", "alpha_metric", "diversity_metric"]:
        if c in df.columns:
            metric_col = c
            break

    value_col = None
    for c in ["delta_alpha", "delta", "delta_value", "gut_minus_oral", "delta_gut_minus_oral"]:
        if c in df.columns:
            value_col = c
            break

    if metric_col and value_col:
        tmp = df[["bioproject", "pair_id", metric_col, value_col]].copy()
        tmp[value_col] = pd.to_numeric(tmp[value_col], errors="coerce")
        tmp["metric_clean"] = tmp[metric_col].astype(str).str.lower()
        metric_map = {
            "observed": "delta_alpha_observed",
            "observed_features": "delta_alpha_observed",
            "shannon": "delta_alpha_shannon",
            "invsimpson": "delta_alpha_invsimpson",
            "inverse_simpson": "delta_alpha_invsimpson",
            "simpson": "delta_alpha_invsimpson",
        }
        tmp["outcome"] = tmp["metric_clean"].map(metric_map)
        tmp = tmp.loc[tmp["outcome"].notna()].copy()
        wide = tmp.pivot_table(
            index=["bioproject", "pair_id"],
            columns="outcome",
            values=value_col,
            aggfunc="first",
        ).reset_index()
        return wide

    # Wide format option.
    patterns = {
        "delta_alpha_observed": [
            "delta_alpha_observed", "delta_observed", "observed_delta",
            "delta_Observed", "Observed_delta", "gut_minus_oral_Observed",
        ],
        "delta_alpha_shannon": [
            "delta_alpha_shannon", "delta_shannon", "shannon_delta",
            "delta_Shannon", "Shannon_delta", "gut_minus_oral_Shannon",
        ],
        "delta_alpha_invsimpson": [
            "delta_alpha_invsimpson", "delta_invsimpson", "delta_invsimpson",
            "delta_inverse_simpson", "invsimpson_delta", "InvSimpson_delta",
            "delta_InvSimpson", "gut_minus_oral_InvSimpson",
        ],
    }

    for out_col, candidates in patterns.items():
        found = None
        for cand in candidates:
            if cand in df.columns:
                found = cand
                break
        if found:
            sub = df[["bioproject", "pair_id", found]].drop_duplicates(subset=["bioproject", "pair_id"])
            sub = sub.rename(columns={found: out_col})
            sub[out_col] = pd.to_numeric(sub[out_col], errors="coerce")
            out = out.merge(sub, on=["bioproject", "pair_id"], how="left")

    return out


def load_beta_distances(beta_fp):
    df = read_table(beta_fp)
    if df.empty:
        return pd.DataFrame(columns=["bioproject", "pair_id"])

    df = normalize_keys(df)

    # Keep complete pairs if available.
    if "has_complete_oral_gut_pair" in df.columns:
        keep = df["has_complete_oral_gut_pair"].astype(str).str.lower().isin({"true", "1", "yes", "y"})
        if keep.any():
            df = df.loc[keep].copy()

    metric_col = None
    for c in ["metric", "distance_metric", "beta_metric"]:
        if c in df.columns:
            metric_col = c
            break

    value_col = None
    for c in [
        "paired_distance", "beta_distance", "distance",
        "oral_gut_distance", "paired_beta_distance",
    ]:
        if c in df.columns:
            value_col = c
            break

    if metric_col and value_col:
        tmp = df[["bioproject", "pair_id", metric_col, value_col]].copy()
        tmp[value_col] = pd.to_numeric(tmp[value_col], errors="coerce")
        tmp["metric_clean"] = tmp[metric_col].astype(str).str.lower()

        metric_map = {
            "braycurtis": "beta_braycurtis_distance",
            "bray_curtis": "beta_braycurtis_distance",
            "weighted_unifrac": "beta_weighted_unifrac_distance",
            "unweighted_unifrac": "beta_unweighted_unifrac_distance",
        }
        tmp["outcome"] = tmp["metric_clean"].map(metric_map)
        tmp = tmp.loc[tmp["outcome"].notna()].copy()

        wide = tmp.pivot_table(
            index=["bioproject", "pair_id"],
            columns="outcome",
            values=value_col,
            aggfunc="first",
        ).reset_index()
        return wide

    # Wide fallback.
    out = df[["bioproject", "pair_id"]].drop_duplicates().copy()
    patterns = {
        "beta_braycurtis_distance": ["beta_braycurtis_distance", "braycurtis_distance", "braycurtis"],
        "beta_weighted_unifrac_distance": ["beta_weighted_unifrac_distance", "weighted_unifrac_distance", "weighted_unifrac"],
        "beta_unweighted_unifrac_distance": ["beta_unweighted_unifrac_distance", "unweighted_unifrac_distance", "unweighted_unifrac"],
    }

    for out_col, candidates in patterns.items():
        found = None
        for cand in candidates:
            if cand in df.columns:
                found = cand
                break
        if found:
            sub = df[["bioproject", "pair_id", found]].drop_duplicates(subset=["bioproject", "pair_id"])
            sub = sub.rename(columns={found: out_col})
            sub[out_col] = pd.to_numeric(sub[out_col], errors="coerce")
            out = out.merge(sub, on=["bioproject", "pair_id"], how="left")

    return out


def load_shared_otu(shared_fp):
    df = read_table(shared_fp)
    if df.empty:
        return pd.DataFrame(columns=["bioproject", "pair_id"])

    df = normalize_keys(df)

    keep_cols = ["bioproject", "pair_id"]
    candidates = [
        "shared_otu_count",
        "shared_otu_count_fraction_of_gut",
        "shared_gut_abundance_fraction",
    ]
    found = [c for c in candidates if c in df.columns]
    keep_cols += found

    out = df[keep_cols].copy()
    out = out.drop_duplicates(subset=["bioproject", "pair_id"])

    for c in found:
        out[c] = pd.to_numeric(out[c], errors="coerce")

    return out


def load_sourcetracker(st2_fp):
    df = read_table(st2_fp)
    if df.empty:
        return pd.DataFrame(columns=["bioproject", "pair_id"])

    df = normalize_keys(df)

    value_col = None
    for c in ["paired_oral_fraction", "source_paired_oral_fraction", "oral_fraction"]:
        if c in df.columns:
            value_col = c
            break

    if value_col is None:
        print(f"WARNING: no paired oral fraction column found in {st2_fp}", flush=True)
        return df[["bioproject", "pair_id"]].drop_duplicates().copy()

    out = df[["bioproject", "pair_id", value_col]].copy()
    out = out.rename(columns={value_col: "sourcetracker_paired_oral_fraction"})
    out = out.drop_duplicates(subset=["bioproject", "pair_id"])
    out["sourcetracker_paired_oral_fraction"] = pd.to_numeric(
        out["sourcetracker_paired_oral_fraction"],
        errors="coerce",
    )

    return out


def main():
    parser = argparse.ArgumentParser(
        description="Build pair-level 10-metric input table for time-series LME."
    )

    parser.add_argument("--metadata", default="00.data/metadata_time_series_11bp_with_sampling_time.csv")
    parser.add_argument("--pair-sampling-time", default="00.data/pair_sampling_time_from_pair_time_point.tsv")
    parser.add_argument("--alpha-delta", default="02.alpha_diversity/delta_alpha.tsv")
    parser.add_argument("--beta-paired", default="03.beta_diversity/distance_matrix/paired_beta_distances_all_pairs.tsv")
    parser.add_argument("--shared-otu", default="04.sourcetracking/share_otu/shared_otu_oral_to_gut_per_sample.tsv")
    parser.add_argument("--sourcetracker2", default="04.sourcetracking/sourcetracker/paired_oral_fraction.tsv")
    parser.add_argument("--out-dir", default="05.lme_time_series_full_model")
    parser.add_argument("--out-name", default="lme_input_10_metrics.tsv")
    parser.add_argument("--overwrite", action="store_true")

    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    out_fp = out_dir / args.out_name
    availability_fp = out_dir / "lme_input_metric_availability.tsv"
    summary_fp = out_dir / "lme_input_merge_summary.json"

    for p in [out_fp, availability_fp, summary_fp]:
        if p.exists() and not args.overwrite:
            raise FileExistsError(f"Output exists: {p}. Use --overwrite.")

    print("Resolved inputs:", flush=True)
    for name, fp in [
        ("metadata", args.metadata),
        ("pair_sampling_time", args.pair_sampling_time),
        ("alpha_delta", args.alpha_delta),
        ("beta_paired", args.beta_paired),
        ("shared_otu", args.shared_otu),
        ("sourcetracker2", args.sourcetracker2),
    ]:
        print(f"  {name}: {Path(fp).resolve()}", flush=True)

    pair_md = make_pair_metadata(args.metadata, args.pair_sampling_time)
    alpha = load_alpha_delta(args.alpha_delta)
    beta = load_beta_distances(args.beta_paired)
    shared = load_shared_otu(args.shared_otu)
    st2 = load_sourcetracker(args.sourcetracker2)

    merged = pair_md.copy()

    source_counts = {
        "pair_metadata_rows": int(pair_md.shape[0]),
        "alpha_rows": int(alpha.shape[0]),
        "beta_rows": int(beta.shape[0]),
        "shared_otu_rows": int(shared.shape[0]),
        "sourcetracker2_rows": int(st2.shape[0]),
    }

    for label, tab in [
        ("alpha", alpha),
        ("beta", beta),
        ("shared_otu", shared),
        ("sourcetracker2", st2),
    ]:
        if tab.empty:
            print(f"WARNING: {label} table is empty; skipping merge.", flush=True)
            continue

        tab = normalize_keys(tab)
        before_cols = set(merged.columns)
        merged = merged.merge(tab, on=["bioproject", "pair_id"], how="left")
        added = [c for c in merged.columns if c not in before_cols]
        print(f"Merged {label}; added columns: {added}", flush=True)

    # Ensure all outcome columns exist.
    for c in OUTCOME_COLUMNS:
        if c not in merged.columns:
            merged[c] = np.nan
        merged[c] = pd.to_numeric(merged[c], errors="coerce")

    # Useful flag: all 10 outcomes available.
    merged["n_available_outcomes"] = merged[OUTCOME_COLUMNS].notna().sum(axis=1)
    merged["has_all_10_outcomes"] = merged["n_available_outcomes"].eq(len(OUTCOME_COLUMNS))

    # Sort.
    sort_cols = [c for c in ["bioproject", "subject_id", "pair_id"] if c in merged.columns]
    merged = merged.sort_values(sort_cols).copy()

    # Availability.
    rows = []
    for c in OUTCOME_COLUMNS:
        rows.append({
            "outcome": c,
            "n_nonmissing": int(merged[c].notna().sum()),
            "n_missing": int(merged[c].isna().sum()),
            "n_total_pairs": int(merged.shape[0]),
            "fraction_nonmissing": float(merged[c].notna().mean()),
        })

    availability = pd.DataFrame(rows)

    # By BioProject availability.
    by_bp_rows = []
    for bp, sub in merged.groupby("bioproject", dropna=False):
        row = {
            "bioproject": bp,
            "n_pairs": int(sub.shape[0]),
            "n_pairs_all_10_outcomes": int(sub["has_all_10_outcomes"].sum()),
        }
        for c in OUTCOME_COLUMNS:
            row[f"{c}__n_nonmissing"] = int(sub[c].notna().sum())
        by_bp_rows.append(row)

    by_bp = pd.DataFrame(by_bp_rows)
    by_bp_fp = out_dir / "lme_input_metric_availability_by_bioproject.tsv"

    merged.to_csv(out_fp, sep="\t", index=False)
    availability.to_csv(availability_fp, sep="\t", index=False)
    by_bp.to_csv(by_bp_fp, sep="\t", index=False)

    summary = {
        "inputs": {
            "metadata": os.path.abspath(args.metadata),
            "pair_sampling_time": os.path.abspath(args.pair_sampling_time),
            "alpha_delta": os.path.abspath(args.alpha_delta),
            "beta_paired": os.path.abspath(args.beta_paired),
            "shared_otu": os.path.abspath(args.shared_otu),
            "sourcetracker2": os.path.abspath(args.sourcetracker2),
        },
        "outputs": {
            "lme_input": os.path.abspath(out_fp),
            "availability": os.path.abspath(availability_fp),
            "availability_by_bioproject": os.path.abspath(by_bp_fp),
        },
        "source_counts": source_counts,
        "n_output_pairs": int(merged.shape[0]),
        "n_pairs_all_10_outcomes": int(merged["has_all_10_outcomes"].sum()),
        "outcome_columns": OUTCOME_COLUMNS,
        "availability": availability.to_dict(orient="records"),
    }

    with open(summary_fp, "w") as f:
        json.dump(summary, f, indent=2)

    print("\nDone.", flush=True)
    print(f"Output: {out_fp}", flush=True)
    print(f"Rows: {merged.shape[0]}", flush=True)
    print(f"Pairs with all 10 outcomes: {int(merged['has_all_10_outcomes'].sum())}", flush=True)
    print("\nOutcome availability:", flush=True)
    print(availability.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
