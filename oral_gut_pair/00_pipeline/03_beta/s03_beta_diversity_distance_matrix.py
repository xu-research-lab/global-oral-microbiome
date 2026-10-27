#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
s05_beta_distance_matrix_integrated.py

Integrated paired beta-diversity distance workflow for the new 32 non-time-series
oral-gut paired analysis.

This script integrates the old paired Bray-Curtis idea and the tree/feature
checking needed for UniFrac.

Default inputs:
  filter metadata:
    /beegfs/zhangcaiyi/oral_gut_pair_noTimeSeries_32bp/01.filter_data/filter_metadata.csv

  filtered count BIOM:
    /beegfs/zhangcaiyi/oral_gut_pair_noTimeSeries_32bp/01.filter_data/filter.biom

  filtered relative-abundance BIOM:
    /beegfs/zhangcaiyi/oral_gut_pair_noTimeSeries_32bp/01.filter_data/filter_ra.biom

  tree:
    /beegfs/zhangcaiyi/oral_gut_pair_noTimeSeries_32bp/00.data/phylogeny_32bp_from_ra_filtered_ra.nwk

Default outputs:
  /beegfs/zhangcaiyi/oral_gut_pair_noTimeSeries_32bp/03.beta_diversity/distance_matrix/
    braycurtis_filter_ra_distance_matrix.tsv
    weighted_unifrac_filter_distance_matrix.tsv
    unweighted_unifrac_filter_distance_matrix.tsv
    paired_beta_distances_all_pairs.tsv
    paired_beta_distances_pairwise_comparisons.tsv
    beta_distance_summary_by_bp.tsv
    beta_distance_run_summary.json
    feature_tree_match_summary.tsv
    feature_tree_match_summary.json
    features_in_biom_and_tree.txt
    features_missing_from_tree.txt
    samples_dropped_unifrac_zero_depth_after_tree_filter.txt

Distances:
  1. Bray-Curtis:
       calculated from filter_ra.biom

  2. Weighted UniFrac:
       calculated from filter.biom after restricting features to tree tips.
       Uses scikit-bio beta_diversity(metric="weighted_unifrac", normalized=True).

  3. Unweighted UniFrac:
       calculated from filter.biom after restricting features to tree tips.

Pairing:
  For each bioproject × pair_id, all Oral × Gut combinations are extracted
  from each distance matrix. The main pair-level distance is the mean distance
  if multiple Oral/Gut samples exist for the same pair_id.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_ROOT = Path("/home/zhangcaiyi/project/oral_gut_pair_TimeSeries_11bp")
DEFAULT_FILTER_DIR = DEFAULT_ROOT / "01.filter_data"
DEFAULT_TREE = DEFAULT_ROOT / "00.data/phylogeny_time_series_11bp.nwk"
DEFAULT_OUT_DIR = DEFAULT_ROOT / "03.beta_diversity" / "distance_matrix"


def read_metadata(path):
    path = str(path)
    if path.endswith((".tsv", ".txt")):
        return pd.read_csv(path, sep="\t", dtype=str)
    return pd.read_csv(path, sep=None, engine="python", dtype=str)


def load_biom_table(path):
    try:
        from biom import load_table
    except ImportError as e:
        raise ImportError(
            "Cannot import biom-format. Install it first:\n"
            "  conda install -c bioconda biom-format\n"
            "or:\n"
            "  pip install biom-format"
        ) from e

    return load_table(str(path))


def load_tree(path):
    try:
        from skbio import TreeNode
    except ImportError as e:
        raise ImportError(
            "Cannot import scikit-bio. Install it first:\n"
            "  conda install -c conda-forge scikit-bio\n"
            "or:\n"
            "  pip install scikit-bio"
        ) from e

    return TreeNode.read(str(path))


def beta_diversity(*args, **kwargs):
    try:
        from skbio.diversity import beta_diversity as _beta_diversity
    except ImportError as e:
        raise ImportError(
            "Cannot import skbio.diversity.beta_diversity. "
            "Install scikit-bio in the active environment."
        ) from e

    return _beta_diversity(*args, **kwargs)


def write_list(items, path):
    with open(path, "w") as f:
        for x in items:
            f.write(str(x) + "\n")


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


def validate_columns(df, cols, label):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"{label} missing required columns:\n  "
            + "\n  ".join(missing)
            + "\n\nAvailable columns:\n  "
            + "\n  ".join(df.columns)
        )


def json_safe_records(df):
    tmp = df.copy()
    tmp = tmp.replace({np.nan: None})
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


def summarize_numeric(values, prefix):
    x = pd.to_numeric(pd.Series(values), errors="coerce").dropna()

    if x.empty:
        return {
            f"{prefix}_min": np.nan,
            f"{prefix}_q25": np.nan,
            f"{prefix}_median": np.nan,
            f"{prefix}_mean": np.nan,
            f"{prefix}_q75": np.nan,
            f"{prefix}_max": np.nan,
        }

    return {
        f"{prefix}_min": float(x.min()),
        f"{prefix}_q25": float(x.quantile(0.25)),
        f"{prefix}_median": float(x.median()),
        f"{prefix}_mean": float(x.mean()),
        f"{prefix}_q75": float(x.quantile(0.75)),
        f"{prefix}_max": float(x.max()),
    }


def check_missing_branch_lengths(tree):
    missing = []
    for node in tree.postorder():
        if node.is_root():
            continue
        if node.length is None:
            missing.append(node.name)
    return missing


def table_to_sample_by_feature_matrix(table, feature_ids=None, sample_ids=None):
    """
    Convert BIOM table to sample x feature dense matrix.

    If feature_ids is provided, table is filtered and ordered by those features.
    If sample_ids is provided, table is filtered and ordered by those samples.
    """
    t = table

    if sample_ids is not None:
        t = t.filter(sample_ids, axis="sample", inplace=False)

    if feature_ids is not None:
        t = t.filter(feature_ids, axis="observation", inplace=False)

    mat = t.matrix_data.T.toarray().astype(np.float64)
    samples = [str(x) for x in t.ids(axis="sample")]
    features = [str(x) for x in t.ids(axis="observation")]

    return mat, samples, features


def write_distance_matrix(distance_matrix, out_path):
    """
    distance_matrix is a skbio DistanceMatrix.
    """
    df = distance_matrix.to_data_frame()
    df.to_csv(out_path, sep="\t", index=True)


def try_midpoint_root(tree):
    """
    Try to midpoint-root a scikit-bio TreeNode.
    If rooting is unavailable or fails, return the original tree.
    """
    try:
        return tree.root_at_midpoint()
    except Exception:
        return tree


def calculate_braycurtis_distance_matrix(ra_biom_path, metadata_sample_order=None):
    table = load_biom_table(ra_biom_path)

    biom_sample_ids = [str(x) for x in table.ids(axis="sample")]

    if metadata_sample_order is None:
        sample_ids = biom_sample_ids
    else:
        sample_set = set(biom_sample_ids)
        sample_ids = [s for s in metadata_sample_order if s in sample_set]

    if not sample_ids:
        raise RuntimeError("No samples available for Bray-Curtis distance calculation.")

    data, sample_ids, feature_ids = table_to_sample_by_feature_matrix(
        table,
        feature_ids=None,
        sample_ids=sample_ids,
    )

    dm = beta_diversity(
        metric="braycurtis",
        counts=data,
        ids=sample_ids,
        validate=False,
    )

    summary = {
        "metric": "braycurtis",
        "input_biom": str(ra_biom_path),
        "input_type": "relative_abundance",
        "n_samples": int(len(sample_ids)),
        "n_features": int(len(feature_ids)),
        "n_zero_sum_samples": int((data.sum(axis=1) <= 0).sum()),
        "note": "Bray-Curtis calculated from filter_ra.biom.",
    }

    return dm, summary


def prepare_unifrac_inputs(count_biom_path, tree_path, metadata_sample_order=None, midpoint_root_tree=False):
    table = load_biom_table(count_biom_path)
    tree = load_tree(tree_path)

    biom_sample_ids = [str(x) for x in table.ids(axis="sample")]
    biom_feature_ids = [str(x) for x in table.ids(axis="observation")]

    if metadata_sample_order is None:
        sample_ids = biom_sample_ids
    else:
        sample_set = set(biom_sample_ids)
        sample_ids = [s for s in metadata_sample_order if s in sample_set]

    if not sample_ids:
        raise RuntimeError("No samples available for UniFrac distance calculation.")

    tree_tips = {str(t.name) for t in tree.tips() if t.name is not None}

    features_in_tree = sorted(set(biom_feature_ids) & tree_tips)
    features_missing_from_tree = sorted(set(biom_feature_ids) - tree_tips)
    tree_tips_not_in_biom = sorted(tree_tips - set(biom_feature_ids))

    if not features_in_tree:
        raise RuntimeError("No BIOM features are present in the tree tips.")

    # Keep BIOM feature order for the count matrix, but only features found in tree.
    feature_ids_for_matrix = [fid for fid in biom_feature_ids if fid in set(features_in_tree)]

    data, sample_ids, feature_ids_for_matrix = table_to_sample_by_feature_matrix(
        table,
        feature_ids=feature_ids_for_matrix,
        sample_ids=sample_ids,
    )

    sample_sums_after_tree_filter = data.sum(axis=1)
    nonzero_sample_mask = sample_sums_after_tree_filter > 0

    dropped_zero_depth_samples = [
        sid for sid, keep in zip(sample_ids, nonzero_sample_mask)
        if not keep
    ]

    if not np.any(nonzero_sample_mask):
        raise RuntimeError("All samples have zero counts after restricting features to tree tips.")

    if dropped_zero_depth_samples:
        data = data[nonzero_sample_mask, :]
        sample_ids = [sid for sid, keep in zip(sample_ids, nonzero_sample_mask) if keep]

    subtree = tree.shear(feature_ids_for_matrix)
    subtree.prune()

    if midpoint_root_tree:
        subtree = try_midpoint_root(subtree)

    missing_branch_lengths = check_missing_branch_lengths(subtree)

    feature_summary = {
        "input_count_biom": str(count_biom_path),
        "input_tree": str(tree_path),
        "n_biom_features": int(len(biom_feature_ids)),
        "n_tree_tips": int(len(tree_tips)),
        "n_features_in_biom_and_tree": int(len(feature_ids_for_matrix)),
        "n_features_missing_from_tree": int(len(features_missing_from_tree)),
        "n_tree_tips_not_in_biom": int(len(tree_tips_not_in_biom)),
        "n_samples_before_tree_feature_filter": int(len(nonzero_sample_mask)),
        "n_samples_after_tree_feature_filter": int(len(sample_ids)),
        "n_samples_dropped_zero_depth_after_tree_filter": int(len(dropped_zero_depth_samples)),
        "n_missing_branch_lengths_non_root_after_pruning": int(len(missing_branch_lengths)),
        "midpoint_root_tree": bool(midpoint_root_tree),
    }

    return (
        data,
        sample_ids,
        feature_ids_for_matrix,
        subtree,
        features_in_tree,
        features_missing_from_tree,
        tree_tips_not_in_biom,
        dropped_zero_depth_samples,
        feature_summary,
    )


def calculate_unifrac_distance_matrix(metric, data, sample_ids, feature_ids, tree, normalized_weighted=True):
    if metric == "weighted_unifrac":
        dm = beta_diversity(
            metric="weighted_unifrac",
            counts=data,
            ids=sample_ids,
            taxa=feature_ids,
            tree=tree,
            normalized=normalized_weighted,
            validate=True,
        )
    elif metric == "unweighted_unifrac":
        dm = beta_diversity(
            metric="unweighted_unifrac",
            counts=data,
            ids=sample_ids,
            taxa=feature_ids,
            tree=tree,
            validate=True,
        )
    else:
        raise ValueError(f"Unsupported UniFrac metric: {metric}")

    return dm


def extract_paired_distances_from_dm(
    dm,
    metadata_df,
    metric_name,
    dataset_label,
    sample_col,
    bioproject_col,
    site_col,
    pair_col,
):
    sample_ids_in_dm = list(dm.ids)
    sample_set = set(sample_ids_in_dm)
    sample_to_index = {sid: i for i, sid in enumerate(sample_ids_in_dm)}

    md = metadata_df.copy()
    md = md.loc[md[sample_col].isin(sample_set), :].copy()
    md[site_col] = md[site_col].apply(standardize_site)
    md = md.loc[md[site_col].isin(["Oral", "Gut"]), :].copy()
    md = md.drop_duplicates(subset=[sample_col]).copy()

    pair_rows = []
    comparison_rows = []

    # DistanceMatrix supports dm[sample_a, sample_b], but indexing the numpy
    # matrix is much faster for many pairs.
    dm_data = dm.data

    for (bp, pair_id), sub in md.groupby([bioproject_col, pair_col], sort=True):
        oral_samples = (
            sub.loc[sub[site_col] == "Oral", sample_col]
            .dropna()
            .astype(str)
            .tolist()
        )

        gut_samples = (
            sub.loc[sub[site_col] == "Gut", sample_col]
            .dropna()
            .astype(str)
            .tolist()
        )

        distances = []

        for oral_sample in oral_samples:
            for gut_sample in gut_samples:
                if oral_sample not in sample_to_index or gut_sample not in sample_to_index:
                    continue

                d = float(dm_data[sample_to_index[oral_sample], sample_to_index[gut_sample]])
                distances.append(d)

                comparison_rows.append({
                    "metric": metric_name,
                    "dataset_label": dataset_label,
                    "bioproject": bp,
                    "pair_id": pair_id,
                    "oral_sample_id": oral_sample,
                    "gut_sample_id": gut_sample,
                    "distance": d,
                })

        dist_series = pd.to_numeric(pd.Series(distances), errors="coerce").dropna()

        if dist_series.empty:
            dist_mean = np.nan
            dist_median = np.nan
            dist_min = np.nan
            dist_max = np.nan
        else:
            dist_mean = float(dist_series.mean())
            dist_median = float(dist_series.median())
            dist_min = float(dist_series.min())
            dist_max = float(dist_series.max())

        pair_rows.append({
            "metric": metric_name,
            "dataset_label": dataset_label,
            "bioproject": bp,
            "pair_id": pair_id,
            "n_oral_samples": int(len(oral_samples)),
            "n_gut_samples": int(len(gut_samples)),
            "n_oral_gut_comparisons": int(len(distances)),
            "has_complete_oral_gut_pair": bool(len(distances) > 0),
            "distance": dist_mean,
            "distance_median_if_multiple": dist_median,
            "distance_min_if_multiple": dist_min,
            "distance_max_if_multiple": dist_max,
            "oral_sample_ids": ";".join(oral_samples),
            "gut_sample_ids": ";".join(gut_samples),
        })

    return pd.DataFrame(pair_rows), pd.DataFrame(comparison_rows)


def summarize_pair_distances(pair_df):
    rows = []

    if pair_df.empty:
        return pd.DataFrame()

    for (metric, dataset_label, bp), sub in pair_df.groupby(["metric", "dataset_label", "bioproject"], sort=True):
        complete = sub.loc[sub["has_complete_oral_gut_pair"], :].copy()
        numeric_summary = summarize_numeric(complete["distance"], "distance")

        row = {
            "metric": metric,
            "dataset_label": dataset_label,
            "bioproject": bp,
            "n_pair_ids_total": int(sub["pair_id"].nunique()),
            "n_complete_oral_gut_pairs": int(complete["pair_id"].nunique()),
            "n_pair_ids_missing_oral": int(((sub["n_oral_samples"] == 0) & (sub["n_gut_samples"] > 0)).sum()),
            "n_pair_ids_missing_gut": int(((sub["n_gut_samples"] == 0) & (sub["n_oral_samples"] > 0)).sum()),
            "n_pairwise_oral_gut_comparisons": int(complete["n_oral_gut_comparisons"].sum()),
        }
        row.update(numeric_summary)
        rows.append(row)

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Calculate Bray-Curtis and UniFrac distance matrices, then extract paired Oral-Gut distances."
        )
    )

    parser.add_argument("--filter-dir", default=str(DEFAULT_FILTER_DIR))
    parser.add_argument("--metadata", default=None, help="Default: <filter-dir>/filter_metadata.csv")
    parser.add_argument("--count-biom", default=None, help="Default: <filter-dir>/filter.biom")
    parser.add_argument("--ra-biom", default=None, help="Default: <filter-dir>/filter_ra.biom")
    parser.add_argument("--tree", default=str(DEFAULT_TREE))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))

    parser.add_argument("--sample-col", default="sample_id")
    parser.add_argument("--bioproject-col", default="bioproject")
    parser.add_argument("--site-col", default="site_group_cleaned")
    parser.add_argument("--pair-col", default="pair_id")

    parser.add_argument(
        "--metrics",
        default="braycurtis,weighted_unifrac,unweighted_unifrac",
        help="Comma-separated metrics to calculate. Default: braycurtis,weighted_unifrac,unweighted_unifrac",
    )

    parser.add_argument(
        "--no-normalized-weighted-unifrac",
        action="store_true",
        help="If set, weighted UniFrac is not normalized. Default: normalized weighted UniFrac.",
    )

    parser.add_argument(
        "--midpoint-root-tree",
        action="store_true",
        help=(
            "Try midpoint-rooting the pruned tree before UniFrac. "
            "Use this only if scikit-bio complains that the tree is not rooted."
        ),
    )

    parser.add_argument("--overwrite", action="store_true")

    args = parser.parse_args()

    filter_dir = Path(args.filter_dir).resolve()
    metadata_fp = Path(args.metadata).resolve() if args.metadata else filter_dir / "filter_metadata.csv"
    count_biom_fp = Path(args.count_biom).resolve() if args.count_biom else filter_dir / "filter.biom"
    ra_biom_fp = Path(args.ra_biom).resolve() if args.ra_biom else filter_dir / "filter_ra.biom"
    tree_fp = Path(args.tree).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    metrics = [m.strip() for m in str(args.metrics).split(",") if m.strip()]

    supported = {"braycurtis", "weighted_unifrac", "unweighted_unifrac"}
    unsupported = sorted(set(metrics) - supported)
    if unsupported:
        raise ValueError(f"Unsupported metrics: {unsupported}. Supported: {sorted(supported)}")

    output_paths = {
        "braycurtis": out_dir / "braycurtis_filter_ra_distance_matrix.tsv",
        "weighted_unifrac": out_dir / "weighted_unifrac_filter_distance_matrix.tsv",
        "unweighted_unifrac": out_dir / "unweighted_unifrac_filter_distance_matrix.tsv",
        "all_pairs": out_dir / "paired_beta_distances_all_pairs.tsv",
        "pairwise": out_dir / "paired_beta_distances_pairwise_comparisons.tsv",
        "summary_by_bp": out_dir / "beta_distance_summary_by_bp.tsv",
        "run_summary": out_dir / "beta_distance_run_summary.json",
        "feature_tree_summary_tsv": out_dir / "feature_tree_match_summary.tsv",
        "feature_tree_summary_json": out_dir / "feature_tree_match_summary.json",
        "features_in_biom_and_tree": out_dir / "features_in_biom_and_tree.txt",
        "features_missing_from_tree": out_dir / "features_missing_from_tree.txt",
        "samples_dropped_unifrac": out_dir / "samples_dropped_unifrac_zero_depth_after_tree_filter.txt",
    }

    for key, p in output_paths.items():
        if key in {"braycurtis", "weighted_unifrac", "unweighted_unifrac"} and key not in metrics:
            continue
        if p.exists() and not args.overwrite:
            raise FileExistsError(f"Output exists: {p}. Use --overwrite.")

    print("Resolved paths:", flush=True)
    print(f"  metadata:   {metadata_fp}", flush=True)
    print(f"  count_biom: {count_biom_fp}", flush=True)
    print(f"  ra_biom:    {ra_biom_fp}", flush=True)
    print(f"  tree:       {tree_fp}", flush=True)
    print(f"  out_dir:    {out_dir}", flush=True)

    if not metadata_fp.exists():
        raise FileNotFoundError(f"Metadata not found: {metadata_fp}")
    if "braycurtis" in metrics and not ra_biom_fp.exists():
        raise FileNotFoundError(f"RA BIOM not found: {ra_biom_fp}")
    if any(m in metrics for m in ["weighted_unifrac", "unweighted_unifrac"]):
        if not count_biom_fp.exists():
            raise FileNotFoundError(f"Count BIOM not found: {count_biom_fp}")
        if not tree_fp.exists():
            raise FileNotFoundError(f"Tree not found: {tree_fp}")

    print("[1/5] Reading metadata ...", flush=True)
    metadata_df = read_metadata(metadata_fp)
    validate_columns(
        metadata_df,
        [args.sample_col, args.bioproject_col, args.site_col, args.pair_col],
        label="metadata",
    )

    for col in [args.sample_col, args.bioproject_col, args.site_col, args.pair_col]:
        metadata_df[col] = metadata_df[col].astype(str).str.strip()

    metadata_df = metadata_df.drop_duplicates(subset=[args.sample_col]).copy()
    metadata_sample_order = metadata_df[args.sample_col].tolist()

    all_pair_dfs = []
    all_comparison_dfs = []
    metric_summaries = {}
    feature_tree_summary = None

    print("[2/5] Calculating distance matrices ...", flush=True)

    if "braycurtis" in metrics:
        print("  Bray-Curtis from filter_ra.biom ...", flush=True)
        dm_bc, summary_bc = calculate_braycurtis_distance_matrix(
            ra_biom_fp,
            metadata_sample_order=metadata_sample_order,
        )
        write_distance_matrix(dm_bc, output_paths["braycurtis"])

        pair_df, comparison_df = extract_paired_distances_from_dm(
            dm_bc,
            metadata_df=metadata_df,
            metric_name="braycurtis",
            dataset_label="filter_ra",
            sample_col=args.sample_col,
            bioproject_col=args.bioproject_col,
            site_col=args.site_col,
            pair_col=args.pair_col,
        )

        all_pair_dfs.append(pair_df)
        all_comparison_dfs.append(comparison_df)
        metric_summaries["braycurtis"] = summary_bc

    need_unifrac = any(m in metrics for m in ["weighted_unifrac", "unweighted_unifrac"])

    if need_unifrac:
        print("  Preparing UniFrac feature/tree matched table ...", flush=True)

        (
            unifrac_data,
            unifrac_sample_ids,
            unifrac_feature_ids,
            unifrac_tree,
            features_in_tree,
            features_missing_from_tree,
            tree_tips_not_in_biom,
            dropped_zero_depth_samples,
            feature_tree_summary,
        ) = prepare_unifrac_inputs(
            count_biom_fp,
            tree_fp,
            metadata_sample_order=metadata_sample_order,
            midpoint_root_tree=args.midpoint_root_tree,
        )

        pd.DataFrame([feature_tree_summary]).to_csv(
            output_paths["feature_tree_summary_tsv"],
            sep="\t",
            index=False,
        )

        with open(output_paths["feature_tree_summary_json"], "w") as f:
            json.dump(feature_tree_summary, f, indent=2)

        write_list(unifrac_feature_ids, output_paths["features_in_biom_and_tree"])
        write_list(features_missing_from_tree, output_paths["features_missing_from_tree"])
        write_list(dropped_zero_depth_samples, output_paths["samples_dropped_unifrac"])

        if feature_tree_summary["n_missing_branch_lengths_non_root_after_pruning"] > 0:
            raise RuntimeError(
                "The pruned tree contains non-root nodes without branch lengths. "
                "UniFrac may be invalid. Please fix branch lengths before running UniFrac."
            )

        if "weighted_unifrac" in metrics:
            print("  Weighted UniFrac from filter.biom + tree ...", flush=True)
            dm_wu = calculate_unifrac_distance_matrix(
                metric="weighted_unifrac",
                data=unifrac_data,
                sample_ids=unifrac_sample_ids,
                feature_ids=unifrac_feature_ids,
                tree=unifrac_tree,
                normalized_weighted=(not args.no_normalized_weighted_unifrac),
            )
            write_distance_matrix(dm_wu, output_paths["weighted_unifrac"])

            pair_df, comparison_df = extract_paired_distances_from_dm(
                dm_wu,
                metadata_df=metadata_df,
                metric_name="weighted_unifrac",
                dataset_label="filter",
                sample_col=args.sample_col,
                bioproject_col=args.bioproject_col,
                site_col=args.site_col,
                pair_col=args.pair_col,
            )

            all_pair_dfs.append(pair_df)
            all_comparison_dfs.append(comparison_df)

            metric_summaries["weighted_unifrac"] = {
                "metric": "weighted_unifrac",
                "input_biom": str(count_biom_fp),
                "input_type": "count",
                "input_tree": str(tree_fp),
                "normalized": bool(not args.no_normalized_weighted_unifrac),
                "n_samples": int(len(unifrac_sample_ids)),
                "n_features_in_biom_and_tree": int(len(unifrac_feature_ids)),
            }

        if "unweighted_unifrac" in metrics:
            print("  Unweighted UniFrac from filter.biom + tree ...", flush=True)
            dm_uwu = calculate_unifrac_distance_matrix(
                metric="unweighted_unifrac",
                data=unifrac_data,
                sample_ids=unifrac_sample_ids,
                feature_ids=unifrac_feature_ids,
                tree=unifrac_tree,
                normalized_weighted=True,
            )
            write_distance_matrix(dm_uwu, output_paths["unweighted_unifrac"])

            pair_df, comparison_df = extract_paired_distances_from_dm(
                dm_uwu,
                metadata_df=metadata_df,
                metric_name="unweighted_unifrac",
                dataset_label="filter",
                sample_col=args.sample_col,
                bioproject_col=args.bioproject_col,
                site_col=args.site_col,
                pair_col=args.pair_col,
            )

            all_pair_dfs.append(pair_df)
            all_comparison_dfs.append(comparison_df)

            metric_summaries["unweighted_unifrac"] = {
                "metric": "unweighted_unifrac",
                "input_biom": str(count_biom_fp),
                "input_type": "count_presence_absence",
                "input_tree": str(tree_fp),
                "n_samples": int(len(unifrac_sample_ids)),
                "n_features_in_biom_and_tree": int(len(unifrac_feature_ids)),
            }

    print("[3/5] Combining paired distances ...", flush=True)

    if all_pair_dfs:
        all_pairs = pd.concat(all_pair_dfs, axis=0, ignore_index=True)
    else:
        all_pairs = pd.DataFrame()

    if all_comparison_dfs:
        all_pairwise = pd.concat(all_comparison_dfs, axis=0, ignore_index=True)
    else:
        all_pairwise = pd.DataFrame()

    all_pairs.to_csv(output_paths["all_pairs"], sep="\t", index=False)
    all_pairwise.to_csv(output_paths["pairwise"], sep="\t", index=False)

    print("[4/5] Summarizing by BioProject ...", flush=True)

    summary_by_bp = summarize_pair_distances(all_pairs)
    summary_by_bp.to_csv(output_paths["summary_by_bp"], sep="\t", index=False)

    print("[5/5] Writing run summary ...", flush=True)

    global_summary = {
        "input_filter_dir": str(filter_dir),
        "input_metadata": str(metadata_fp),
        "input_count_biom": str(count_biom_fp),
        "input_ra_biom": str(ra_biom_fp),
        "input_tree": str(tree_fp),
        "out_dir": str(out_dir),
        "metrics_requested": metrics,
        "metric_summaries": metric_summaries,
        "feature_tree_summary": feature_tree_summary,
        "pairing": {
            "pair_column": args.pair_col,
            "site_column": args.site_col,
            "site_comparison": "Oral vs Gut within pair_id",
            "multiple_samples_per_pair_handling": (
                "All Oral x Gut combinations are extracted. "
                "The pair-level distance is the mean across combinations."
            ),
        },
        "output_files": {k: str(v) for k, v in output_paths.items()},
        "n_rows_all_pairs": int(all_pairs.shape[0]),
        "n_rows_pairwise_comparisons": int(all_pairwise.shape[0]),
        "n_rows_summary_by_bp": int(summary_by_bp.shape[0]),
    }

    with open(output_paths["run_summary"], "w") as f:
        json.dump(global_summary, f, indent=2)

    print("\nDone.", flush=True)
    print("Main outputs:", flush=True)
    for key in ["braycurtis", "weighted_unifrac", "unweighted_unifrac", "all_pairs", "pairwise", "summary_by_bp", "run_summary"]:
        if key in {"braycurtis", "weighted_unifrac", "unweighted_unifrac"} and key not in metrics:
            continue
        print(f"  {key}: {output_paths[key]}", flush=True)


if __name__ == "__main__":
    main()
