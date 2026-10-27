#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
s02_alpha_diversity.py

Alpha diversity workflow for 11 time-series oral-gut paired BioProjects.

Main tasks:
  1. Read filtered count BIOM and filtered relative abundance BIOM.
  2. Calculate alpha diversity for each sample:
       Observed
       Shannon
       InvSimpson
  3. Merge metadata.
  4. Calculate paired delta alpha:
       delta_alpha = Gut alpha - Oral alpha

For time-series data, delta alpha is calculated by default within:
  bioproject + pair_id + pair_time_point

If multiple Oral or Gut samples exist within the same group, values are averaged
within site before calculating Gut - Oral.

No LME models are fitted in this script.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_ROOT = Path("/home/zhangcaiyi/project/oral_gut_pair_TimeSeries_11bp")
DEFAULT_FILTER_DIR = DEFAULT_ROOT / "00.data/01.filter_data"
DEFAULT_OUT_DIR = DEFAULT_ROOT / "02.alpha_diversity"


def parse_list(x):
    if x is None:
        return []
    return [i.strip() for i in str(x).split(",") if i.strip()]


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
            "Cannot import biom-format. Activate/install biom-format first:\n"
            "  conda install -c bioconda biom-format\n"
            "or\n"
            "  pip install biom-format"
        ) from e

    return load_table(str(path))


def validate_columns(df, cols, label):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"{label} missing required columns:\n  "
            + "\n  ".join(missing)
            + "\n\nAvailable columns:\n  "
            + "\n  ".join(df.columns)
        )


def json_safe(obj):
    if isinstance(obj, dict):
        return {k: json_safe(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [json_safe(v) for v in obj]
    if isinstance(obj, tuple):
        return [json_safe(v) for v in obj]
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        if np.isnan(obj):
            return None
        return float(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    if pd.isna(obj) if not isinstance(obj, (list, dict, tuple)) else False:
        return None
    return obj


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
        "supragingival plaque": "Oral",
        "subgingival plaque": "Oral",
        "gut": "Gut",
        "stool": "Gut",
        "fecal": "Gut",
        "faecal": "Gut",
        "feces": "Gut",
        "faeces": "Gut",
    }

    return site_map.get(xl, x0)


def join_unique_ids(x):
    vals = sorted(set(str(i) for i in x if pd.notna(i)))
    return ";".join(vals)


def calculate_alpha_from_biom(table):
    """
    Calculate alpha diversity from a BIOM table.

    BIOM matrix is feature x sample.
    Values are internally normalized per sample, so this works for both
    count BIOM and relative abundance BIOM.
    """
    mat = table.matrix_data.tocsc().astype(np.float64)

    sample_ids = [str(x) for x in table.ids(axis="sample")]
    n_features = len(table.ids(axis="observation"))

    sample_sums = np.asarray(mat.sum(axis=0)).ravel().astype(np.float64)

    rows = []

    for j, sid in enumerate(sample_ids):
        start = mat.indptr[j]
        end = mat.indptr[j + 1]

        values = mat.data[start:end]
        values = values[values > 0]
        sample_sum = float(sample_sums[j])

        observed = int(values.size)

        if sample_sum <= 0 or values.size == 0:
            shannon = np.nan
            simpson_d = np.nan
            invsimpson = np.nan
            exp_shannon = np.nan
        else:
            p = values / sample_sum

            # Defensive renormalization for RA BIOM floating-point safety.
            psum = p.sum()
            if psum > 0:
                p = p / psum

            shannon = float(-np.sum(p * np.log(p)))
            simpson_d = float(np.sum(p * p))
            invsimpson = float(1.0 / simpson_d) if simpson_d > 0 else np.nan
            exp_shannon = float(np.exp(shannon))

        rows.append({
            "sample_id": sid,
            "sample_total": sample_sum,
            "n_features_in_table": int(n_features),
            "Observed": observed,
            "Shannon": shannon,
            "InvSimpson": invsimpson,
            "Simpson_D": simpson_d,
            "ExpShannon": exp_shannon,
        })

    return pd.DataFrame(rows)


def summarize_alpha(df):
    metrics = ["Observed", "Shannon", "InvSimpson"]

    out = {
        "n_rows": int(df.shape[0]),
        "n_samples": int(df["sample_id"].nunique()),
        "sample_total_min": float(df["sample_total"].min()) if df.shape[0] else None,
        "sample_total_median": float(df["sample_total"].median()) if df.shape[0] else None,
        "sample_total_mean": float(df["sample_total"].mean()) if df.shape[0] else None,
        "sample_total_max": float(df["sample_total"].max()) if df.shape[0] else None,
    }

    for metric in metrics:
        x = pd.to_numeric(df[metric], errors="coerce").dropna()

        if x.empty:
            for stat in ["min", "q25", "median", "mean", "q75", "max"]:
                out[f"{metric}_{stat}"] = None
        else:
            out[f"{metric}_min"] = float(x.min())
            out[f"{metric}_q25"] = float(x.quantile(0.25))
            out[f"{metric}_median"] = float(x.median())
            out[f"{metric}_mean"] = float(x.mean())
            out[f"{metric}_q75"] = float(x.quantile(0.75))
            out[f"{metric}_max"] = float(x.max())

    return out


def calculate_delta_alpha(
    alpha,
    dataset_label,
    metrics,
    group_cols,
    sample_col,
    site_col,
):
    """
    Calculate delta alpha = Gut - Oral.

    group_cols usually:
      bioproject,pair_id,pair_time_point

    If pair_time_point is not desired:
      bioproject,pair_id
    """
    df = alpha.loc[alpha["dataset_label"].astype(str) == dataset_label, :].copy()

    if df.empty:
        raise RuntimeError(f"No alpha rows found for dataset_label = {dataset_label}")

    required = ["dataset_label", sample_col, site_col] + group_cols + metrics
    validate_columns(df, required, "alpha table for delta alpha")

    df[site_col] = df[site_col].apply(standardize_site)
    df = df.loc[df[site_col].isin(["Oral", "Gut"]), :].copy()

    for c in [sample_col, site_col] + group_cols:
        df[c] = df[c].astype(str).str.strip()

    for m in metrics:
        df[m] = pd.to_numeric(df[m], errors="coerce")

    rows = []

    for key, sub in df.groupby(group_cols, sort=True, dropna=False):
        if not isinstance(key, tuple):
            key = (key,)

        oral = sub.loc[sub[site_col] == "Oral", :].copy()
        gut = sub.loc[sub[site_col] == "Gut", :].copy()

        n_oral = oral[sample_col].nunique()
        n_gut = gut[sample_col].nunique()

        row = {
            "dataset_label": dataset_label,
            "n_oral_samples": int(n_oral),
            "n_gut_samples": int(n_gut),
            "has_complete_oral_gut_pair": bool(n_oral > 0 and n_gut > 0),
            "oral_sample_ids": join_unique_ids(oral[sample_col]),
            "gut_sample_ids": join_unique_ids(gut[sample_col]),
        }

        for c, v in zip(group_cols, key):
            row[c] = v

        for m in metrics:
            oral_value = oral[m].mean(skipna=True) if not oral.empty else np.nan
            gut_value = gut[m].mean(skipna=True) if not gut.empty else np.nan

            row[f"oral_{m}"] = oral_value
            row[f"gut_{m}"] = gut_value
            row[f"delta_{m}"] = (
                gut_value - oral_value
                if pd.notna(gut_value) and pd.notna(oral_value)
                else np.nan
            )

        rows.append(row)

    out = pd.DataFrame(rows)

    front_cols = [
        "dataset_label",
    ] + group_cols + [
        "n_oral_samples",
        "n_gut_samples",
        "has_complete_oral_gut_pair",
        "oral_sample_ids",
        "gut_sample_ids",
    ]

    metric_cols = []
    for m in metrics:
        metric_cols.extend([f"oral_{m}", f"gut_{m}", f"delta_{m}"])

    out = out[front_cols + metric_cols].copy()

    sort_cols = ["has_complete_oral_gut_pair"] + group_cols
    out = out.sort_values(
        sort_cols,
        ascending=[False] + [True] * len(group_cols),
    ).reset_index(drop=True)

    return out


def make_summary_by_bp_site(alpha, dataset_label, metrics, bioproject_col, site_col):
    d = alpha.loc[alpha["dataset_label"].astype(str) == dataset_label, :].copy()
    d[site_col] = d[site_col].apply(standardize_site)

    rows = []

    for (bp, site), sub in d.groupby([bioproject_col, site_col], sort=True, dropna=False):
        row = {
            "dataset_label": dataset_label,
            "bioproject": bp,
            "site": site,
            "n_samples": int(sub["sample_id"].nunique()),
        }

        for m in metrics:
            x = pd.to_numeric(sub[m], errors="coerce").dropna()
            if x.empty:
                row[f"{m}_median"] = np.nan
                row[f"{m}_mean"] = np.nan
                row[f"{m}_q25"] = np.nan
                row[f"{m}_q75"] = np.nan
            else:
                row[f"{m}_median"] = float(x.median())
                row[f"{m}_mean"] = float(x.mean())
                row[f"{m}_q25"] = float(x.quantile(0.25))
                row[f"{m}_q75"] = float(x.quantile(0.75))

        rows.append(row)

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Calculate alpha diversity and paired delta alpha for time-series 11bp data. No LME."
    )

    parser.add_argument("--filter-dir", default=str(DEFAULT_FILTER_DIR))
    parser.add_argument("--count-biom", default=None, help="Default: <filter-dir>/filter.biom")
    parser.add_argument("--ra-biom", default=None, help="Default: <filter-dir>/filter_ra.biom")
    parser.add_argument("--metadata", default=None, help="Default: <filter-dir>/filter_metadata.csv")
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))

    parser.add_argument("--sample-col", default="sample_id")
    parser.add_argument("--bioproject-col", default="bioproject")
    parser.add_argument("--site-col", default="site_group_cleaned")
    parser.add_argument("--pair-col", default="pair_id")
    parser.add_argument("--time-col", default="pair_time_point")

    parser.add_argument("--calculate-dataset-labels", default="filter,filter_ra")
    parser.add_argument("--delta-dataset-label", default="filter_ra")
    parser.add_argument("--metrics", default="Observed,Shannon,InvSimpson")

    parser.add_argument(
        "--delta-group-cols",
        default="bioproject,pair_id,pair_time_point",
        help=(
            "Grouping columns for delta alpha. "
            "Default: bioproject,pair_id,pair_time_point. "
            "For non-time-series use bioproject,pair_id."
        ),
    )

    parser.add_argument("--overwrite", action="store_true")

    args = parser.parse_args()

    filter_dir = Path(args.filter_dir).resolve()
    count_biom = Path(args.count_biom).resolve() if args.count_biom else filter_dir / "filter.biom"
    ra_biom = Path(args.ra_biom).resolve() if args.ra_biom else filter_dir / "filter_ra.biom"
    metadata_fp = Path(args.metadata).resolve() if args.metadata else filter_dir / "filter_metadata.csv"
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    out_alpha_tsv = out_dir / "alpha_diversity_all_samples.tsv"
    out_dataset_summary_tsv = out_dir / "alpha_diversity_dataset_summary.tsv"
    out_bp_site_summary_tsv = out_dir / "alpha_diversity_by_bioproject_site.tsv"
    out_delta_tsv = out_dir / "delta_alpha.tsv"
    out_delta_summary_tsv = out_dir / "delta_alpha_summary_by_bioproject.tsv"
    out_summary_json = out_dir / "alpha_diversity_summary.json"

    outputs = [
        out_alpha_tsv,
        out_dataset_summary_tsv,
        out_bp_site_summary_tsv,
        out_delta_tsv,
        out_delta_summary_tsv,
        out_summary_json,
    ]

    if not args.overwrite:
        existing = [p for p in outputs if p.exists()]
        if existing:
            raise FileExistsError(f"Output exists: {existing[0]}. Use --overwrite.")

    print("Resolved paths:", flush=True)
    print(f"  count_biom: {count_biom}", flush=True)
    print(f"  ra_biom:    {ra_biom}", flush=True)
    print(f"  metadata:   {metadata_fp}", flush=True)
    print(f"  out_dir:    {out_dir}", flush=True)

    for p, label in [(count_biom, "count BIOM"), (ra_biom, "RA BIOM"), (metadata_fp, "metadata")]:
        if not p.exists():
            raise FileNotFoundError(f"{label} not found: {p}")

    metrics = parse_list(args.metrics)
    calculate_labels = parse_list(args.calculate_dataset_labels)
    delta_group_cols = parse_list(args.delta_group_cols)

    biom_inputs = {
        "filter": count_biom,
        "filter_ra": ra_biom,
    }

    bad_labels = [x for x in calculate_labels if x not in biom_inputs]
    if bad_labels:
        raise ValueError(f"Unsupported dataset labels: {bad_labels}. Supported: {list(biom_inputs)}")

    if args.delta_dataset_label not in biom_inputs:
        raise ValueError(f"Unsupported delta dataset label: {args.delta_dataset_label}")

    print("[1/5] Reading metadata ...", flush=True)
    md = read_metadata(metadata_fp)

    required_meta_cols = list(set(
        [args.sample_col, args.bioproject_col, args.site_col, args.pair_col]
        + delta_group_cols
    ))

    validate_columns(md, required_meta_cols, "metadata")

    for c in required_meta_cols:
        md[c] = md[c].astype(str).str.strip()

    duplicated_samples = (
        md.loc[md[args.sample_col].duplicated(keep=False), args.sample_col]
        .sort_values()
        .unique()
        .tolist()
    )

    if duplicated_samples:
        dup_fp = out_dir / "duplicated_sample_ids_in_metadata.txt"
        with open(dup_fp, "w") as f:
            for x in duplicated_samples:
                f.write(str(x) + "\n")
        raise RuntimeError(f"Duplicated sample IDs in metadata: {len(duplicated_samples)}. See {dup_fp}")

    md = md.rename(columns={args.sample_col: "sample_id"})

    # Update group column names if sample_col was renamed.
    sample_col = "sample_id"

    print("[2/5] Calculating alpha diversity ...", flush=True)

    all_alpha = []
    summary_rows = []

    for dataset_label in calculate_labels:
        biom_path = biom_inputs[dataset_label]
        print(f"  Dataset={dataset_label}: {biom_path}", flush=True)

        table = load_biom_table(biom_path)
        alpha_df = calculate_alpha_from_biom(table)

        alpha_df.insert(0, "dataset_label", dataset_label)
        alpha_df.insert(1, "input_biom", str(biom_path))

        dataset_summary = summarize_alpha(alpha_df)
        dataset_summary.update({
            "dataset_label": dataset_label,
            "input_biom": str(biom_path),
            "n_features_in_table": int(alpha_df["n_features_in_table"].iloc[0]) if not alpha_df.empty else 0,
        })

        summary_rows.append(dataset_summary)
        all_alpha.append(alpha_df)

    alpha = pd.concat(all_alpha, ignore_index=True)

    print("[3/5] Merging metadata ...", flush=True)

    alpha = alpha.merge(md, on="sample_id", how="left", validate="many_to_one")

    missing_md = int(alpha[args.bioproject_col].isna().sum())
    if missing_md > 0:
        raise RuntimeError(f"{missing_md} alpha rows did not match metadata by sample_id.")

    for c in [args.site_col]:
        alpha[c] = alpha[c].apply(standardize_site)

    alpha.to_csv(out_alpha_tsv, sep="\t", index=False)

    dataset_summary_df = pd.DataFrame(summary_rows)
    dataset_summary_df.to_csv(out_dataset_summary_tsv, sep="\t", index=False)

    bp_site_summary = []
    for dataset_label in calculate_labels:
        bp_site_summary.append(
            make_summary_by_bp_site(
                alpha,
                dataset_label=dataset_label,
                metrics=metrics,
                bioproject_col=args.bioproject_col,
                site_col=args.site_col,
            )
        )

    bp_site_summary = pd.concat(bp_site_summary, ignore_index=True)
    bp_site_summary.to_csv(out_bp_site_summary_tsv, sep="\t", index=False)

    print("[4/5] Calculating delta alpha ...", flush=True)

    # Keep only group columns that actually exist after metadata merge.
    missing_delta_group_cols = [c for c in delta_group_cols if c not in alpha.columns]
    if missing_delta_group_cols:
        raise ValueError(
            "Delta group columns missing from alpha table:\n  "
            + "\n  ".join(missing_delta_group_cols)
        )

    delta = calculate_delta_alpha(
        alpha=alpha,
        dataset_label=args.delta_dataset_label,
        metrics=metrics,
        group_cols=delta_group_cols,
        sample_col=sample_col,
        site_col=args.site_col,
    )

    delta.to_csv(out_delta_tsv, sep="\t", index=False)

    delta_summary_rows = []
    for bp, sub in delta.groupby(args.bioproject_col, sort=True):
        row = {
            "dataset_label": args.delta_dataset_label,
            "bioproject": bp,
            "n_delta_rows": int(sub.shape[0]),
            "n_complete_oral_gut_rows": int(sub["has_complete_oral_gut_pair"].sum()),
            "n_incomplete_rows": int((~sub["has_complete_oral_gut_pair"]).sum()),
        }

        for m in metrics:
            x = pd.to_numeric(
                sub.loc[sub["has_complete_oral_gut_pair"], f"delta_{m}"],
                errors="coerce",
            ).dropna()

            if x.empty:
                row[f"delta_{m}_median"] = np.nan
                row[f"delta_{m}_mean"] = np.nan
                row[f"delta_{m}_q25"] = np.nan
                row[f"delta_{m}_q75"] = np.nan
            else:
                row[f"delta_{m}_median"] = float(x.median())
                row[f"delta_{m}_mean"] = float(x.mean())
                row[f"delta_{m}_q25"] = float(x.quantile(0.25))
                row[f"delta_{m}_q75"] = float(x.quantile(0.75))

        delta_summary_rows.append(row)

    delta_summary = pd.DataFrame(delta_summary_rows)
    delta_summary.to_csv(out_delta_summary_tsv, sep="\t", index=False)

    print("[5/5] Writing run summary ...", flush=True)

    summary = {
        "filter_dir": str(filter_dir),
        "count_biom": str(count_biom),
        "ra_biom": str(ra_biom),
        "metadata": str(metadata_fp),
        "out_dir": str(out_dir),
        "outputs": {
            "alpha_diversity_all_samples": str(out_alpha_tsv),
            "alpha_diversity_dataset_summary": str(out_dataset_summary_tsv),
            "alpha_diversity_by_bioproject_site": str(out_bp_site_summary_tsv),
            "delta_alpha": str(out_delta_tsv),
            "delta_alpha_summary_by_bioproject": str(out_delta_summary_tsv),
        },
        "dataset_labels_calculated": calculate_labels,
        "delta_dataset_label": args.delta_dataset_label,
        "metrics": metrics,
        "delta_definition": "delta_alpha = Gut alpha - Oral alpha",
        "delta_group_cols": delta_group_cols,
        "notes": [
            "No LME models are fitted in this script.",
            "Alpha diversity is calculated using internally sample-wise normalized values.",
            "For delta alpha, multiple Oral or Gut samples within the same group are averaged before Gut - Oral.",
            "For time-series data, default delta grouping includes pair_time_point to avoid averaging across time points.",
        ],
        "metadata_columns": {
            "sample_col": "sample_id",
            "bioproject_col": args.bioproject_col,
            "site_col": args.site_col,
            "pair_col": args.pair_col,
            "time_col": args.time_col,
        },
        "n_alpha_rows": int(alpha.shape[0]),
        "n_unique_samples": int(alpha["sample_id"].nunique()),
        "n_bioprojects": int(alpha[args.bioproject_col].nunique()),
        "n_delta_rows": int(delta.shape[0]),
        "n_delta_complete_oral_gut_rows": int(delta["has_complete_oral_gut_pair"].sum()),
    }

    with open(out_summary_json, "w") as f:
        json.dump(json_safe(summary), f, indent=2)

    print("\nDone.", flush=True)
    print("Main outputs:", flush=True)
    print(f"  alpha table:       {out_alpha_tsv}", flush=True)
    print(f"  dataset summary:   {out_dataset_summary_tsv}", flush=True)
    print(f"  BP/site summary:   {out_bp_site_summary_tsv}", flush=True)
    print(f"  delta alpha:       {out_delta_tsv}", flush=True)
    print(f"  delta summary:     {out_delta_summary_tsv}", flush=True)
    print(f"  run summary JSON:  {out_summary_json}", flush=True)


if __name__ == "__main__":
    main()
