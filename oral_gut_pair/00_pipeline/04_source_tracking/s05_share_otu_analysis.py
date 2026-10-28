#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
s05_share_otu_analysis.py

Shared OTU oral-to-gut analysis for the new 32 non-time-series workflow.

Input:
  01.filter_data/filter.biom
  01.filter_data/filter_metadata.csv

Output:
  04.sourcetracking/share_otu/
    shared_otu_oral_to_gut_per_sample.tsv
    shared_otu_oral_to_gut_summary_by_bp.tsv
    shared_otu_oral_to_gut_summary_pooled.tsv
    shared_otu_oral_to_gut_long_plot_table.tsv
    shared_otu_oral_to_gut_summary.json
    plots/

Definition:
  For each Gut sample:
    paired oral OTU set = union of OTUs detected in Oral sample(s)
                          from the same bioproject and same pair_id.

  shared_otu_count:
    Number of gut-detected OTUs also detected in paired Oral sample(s).

  shared_otu_count_fraction_of_gut:
    shared_otu_count / number of OTUs detected in the Gut sample.

  shared_gut_abundance_fraction:
    Sum of Gut counts for shared OTUs / total Gut counts.

Important:
  This is OTU sharing / overlap, not proof of oral-to-gut transmission.
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


DEFAULT_ROOT = Path("/beegfs/zhangcaiyi/oral_gut_pair_noTimeSeries_32bp")
DEFAULT_METADATA = DEFAULT_ROOT / "01.filter_data" / "filter_metadata.csv"
DEFAULT_BIOM = DEFAULT_ROOT / "01.filter_data" / "filter.biom"
DEFAULT_OUT_DIR = DEFAULT_ROOT / "04.sourcetracking" / "share_otu"


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
        if v == "" or v.lower() == "nan":
            continue
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
    x = str(x).strip()
    if x == "" or x.lower() in {"nan", "none", "na", "null"}:
        return "NA"
    return x


def is_hc_status(x):
    s = str(x).strip().lower()
    return s in {
        "hc", "h", "healthy", "health", "control", "controls", "normal",
        "healthy_control", "healthy control", "non-disease", "non_disease",
        "non disease", "non-diseased", "non_diseased", "non diseased",
        "noncase", "non-case", "non_case",
    }


def status_sort_key(x):
    x = normalize_status_value(x)
    if x == "NA":
        return (2, x)
    if is_hc_status(x):
        return (0, x)
    return (1, x)


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


def build_bp_order_from_metadata(metadata, bioproject_col, status_col):
    rows = []

    for bp, sub in metadata.groupby(bioproject_col, dropna=False):
        if status_col in sub.columns:
            statuses = (
                sub[status_col]
                .dropna()
                .astype(str)
                .map(normalize_status_value)
                .unique()
                .tolist()
            )
        else:
            statuses = ["NA"]

        statuses_clean = [s for s in statuses if s != "NA"]
        status_text = "NA" if not statuses_clean else "+".join(sorted(statuses_clean, key=status_sort_key))
        status_rank, status_group = classify_bp_status_group(statuses)

        rows.append(
            {
                "bioproject": bp,
                "bp_status_label": f"{bp}({status_text})",
                "status_rank": status_rank,
                "status_group_for_order": status_group,
                "status_values": status_text,
            }
        )

    out = pd.DataFrame(rows).sort_values(
        ["status_rank", "bioproject"],
        ascending=[True, True],
    )
    out["bp_plot_order"] = np.arange(out.shape[0], dtype=int)
    return out


def validate_columns(df, required, label):
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"{label} missing required columns:\n  "
            + "\n  ".join(missing)
            + "\n\nAvailable columns:\n  "
            + "\n  ".join(df.columns)
        )


def load_biom_sparse_sample_by_feature(path):
    try:
        from biom import load_table
    except ImportError as e:
        raise ImportError(
            "Cannot import biom-format. Install it first:\n"
            "  conda install -c bioconda biom-format\n"
            "or:\n"
            "  pip install biom-format"
        ) from e

    table = load_table(str(path))

    # BIOM is feature x sample; transpose to sample x feature CSR.
    X = table.matrix_data.T.tocsr().astype(np.float64)
    sample_ids = [str(x) for x in table.ids(axis="sample")]
    feature_ids = [str(x) for x in table.ids(axis="observation")]

    return X, sample_ids, feature_ids


def prepare_metadata(metadata, sample_col, bioproject_col, site_col, pair_col, subject_col, status_col, disease_col):
    md = metadata.copy()

    required = [sample_col, bioproject_col, site_col, pair_col]
    validate_columns(md, required, "metadata")

    for c in [subject_col, status_col, disease_col]:
        if c not in md.columns:
            md[c] = "NA"

    for c in [sample_col, bioproject_col, site_col, pair_col, subject_col, status_col, disease_col]:
        md[c] = md[c].astype(str).str.strip()

    md["_site_std"] = md[site_col].apply(standardize_site)
    md = md.loc[md["_site_std"].isin(["Oral", "Gut"]), :].copy()
    md = md.drop_duplicates(subset=[sample_col]).copy()

    return md


def get_oral_presence_union(X, oral_indices, presence_threshold):
    if len(oral_indices) == 0:
        return np.zeros(X.shape[1], dtype=bool)

    oral_mat = X[oral_indices, :]

    if presence_threshold <= 0:
        presence = np.asarray(oral_mat.getnnz(axis=0)).ravel() > 0
    else:
        presence = np.asarray((oral_mat > presence_threshold).getnnz(axis=0)).ravel() > 0

    return presence.astype(bool)


def calculate_shared_otu(
    X,
    sample_ids,
    feature_ids,
    metadata,
    dataset_label,
    input_biom,
    sample_col,
    bioproject_col,
    pair_col,
    subject_col,
    status_col,
    disease_col,
    presence_threshold,
):
    sample_to_index = {sid: i for i, sid in enumerate(sample_ids)}
    sample_set = set(sample_ids)

    md = metadata.copy()
    md = md.loc[md[sample_col].isin(sample_set), :].copy()
    md = md.drop_duplicates(subset=[sample_col]).copy()

    md["_sample_index"] = md[sample_col].map(sample_to_index)
    md = md.dropna(subset=["_sample_index"]).copy()
    md["_sample_index"] = md["_sample_index"].astype(int)

    rows = []

    for (bp, pair_id), sub in md.groupby([bioproject_col, pair_col], sort=True):
        oral_sub = sub.loc[sub["_site_std"] == "Oral", :].copy()
        gut_sub = sub.loc[sub["_site_std"] == "Gut", :].copy()

        if oral_sub.empty or gut_sub.empty:
            continue

        oral_indices = oral_sub["_sample_index"].astype(int).tolist()
        oral_ids = oral_sub[sample_col].astype(str).tolist()

        oral_presence = get_oral_presence_union(
            X=X,
            oral_indices=oral_indices,
            presence_threshold=presence_threshold,
        )
        paired_oral_otu_count_union = int(oral_presence.sum())

        for _, gut_row in gut_sub.iterrows():
            gut_id = str(gut_row[sample_col])
            gut_idx = int(gut_row["_sample_index"])

            gut_vector = X.getrow(gut_idx)
            gut_total_counts = float(gut_vector.sum())

            if gut_total_counts <= 0:
                continue

            indices = gut_vector.indices
            data = gut_vector.data.astype(float)

            if presence_threshold <= 0:
                present_mask = data > 0
            else:
                present_mask = data > presence_threshold

            gut_present_indices = indices[present_mask]
            gut_present_counts = data[present_mask]
            gut_otu_count = int(len(gut_present_indices))

            if gut_otu_count > 0:
                shared_feature_mask = oral_presence[gut_present_indices]
                shared_otu_count = int(shared_feature_mask.sum())
                shared_gut_counts = float(gut_present_counts[shared_feature_mask].sum())
                shared_otu_count_fraction_of_gut = shared_otu_count / gut_otu_count
            else:
                shared_otu_count = 0
                shared_gut_counts = 0.0
                shared_otu_count_fraction_of_gut = np.nan

            shared_gut_abundance_fraction = shared_gut_counts / gut_total_counts

            rows.append(
                {
                    "bioproject": bp,
                    "dataset_label": dataset_label,
                    "pair_id": pair_id,
                    "gut_sample_id": gut_id,
                    "oral_sample_ids": ";".join(sorted(set(oral_ids))),
                    "n_oral_samples_in_pair": int(len(set(oral_ids))),
                    "n_gut_samples_in_pair": int(gut_sub[sample_col].nunique()),
                    "subject_id": gut_row.get(subject_col, "NA"),
                    "gut_status_cleaned": gut_row.get(status_col, "NA"),
                    "gut_disease_study_cleaned": gut_row.get(disease_col, "NA"),
                    "gut_total_counts": gut_total_counts,
                    "gut_otu_count": gut_otu_count,
                    "paired_oral_otu_count_union": paired_oral_otu_count_union,
                    "shared_otu_count": shared_otu_count,
                    "shared_otu_count_fraction_of_gut": shared_otu_count_fraction_of_gut,
                    "shared_gut_counts": shared_gut_counts,
                    "shared_gut_abundance_fraction": shared_gut_abundance_fraction,
                    "presence_threshold": float(presence_threshold),
                    "input_biom": str(input_biom),
                    "n_features_in_biom": int(len(feature_ids)),
                    "n_samples_in_biom": int(len(sample_ids)),
                }
            )

    return pd.DataFrame(rows)


def annotate_bp_order(df, bp_order):
    """
    Add BioProject plotting order/status labels.

    Robust to input tables that already contain columns such as status_values,
    bp_status_label, status_group_for_order, or bp_plot_order. This avoids
    pandas merge suffixes like status_values_x/status_values_y.
    """
    out = df.copy()

    if "bioproject" not in out.columns:
        raise ValueError("Input table to annotate_bp_order() missing column: bioproject")

    if bp_order is None or bp_order.empty:
        bp_order2 = pd.DataFrame({"bioproject": sorted(out["bioproject"].dropna().astype(str).unique())})
        bp_order2["bp_status_label"] = bp_order2["bioproject"]
        bp_order2["status_group_for_order"] = "NA"
        bp_order2["status_values"] = "NA"
        bp_order2["bp_plot_order"] = range(bp_order2.shape[0])
    else:
        keep_cols = [
            "bioproject",
            "bp_status_label",
            "status_group_for_order",
            "status_values",
            "bp_plot_order",
        ]
        keep_cols = [c for c in keep_cols if c in bp_order.columns]
        bp_order2 = bp_order[keep_cols].drop_duplicates(subset=["bioproject"]).copy()

    # Avoid duplicated columns after merge, e.g. status_values_x/status_values_y.
    merge_cols = [
        "bp_status_label",
        "status_group_for_order",
        "status_values",
        "bp_plot_order",
    ]
    drop_cols = [c for c in merge_cols if c in out.columns and c in bp_order2.columns]
    if drop_cols:
        out = out.drop(columns=drop_cols)

    out = out.merge(bp_order2, on="bioproject", how="left")

    default_map = {
        bp: i for i, bp in enumerate(sorted(out["bioproject"].dropna().astype(str).unique()))
    }

    if "bp_plot_order" not in out.columns:
        out["bp_plot_order"] = out["bioproject"].map(default_map)
    else:
        out["bp_plot_order"] = pd.to_numeric(out["bp_plot_order"], errors="coerce")
        out["bp_plot_order"] = out["bp_plot_order"].fillna(out["bioproject"].map(default_map))

    out["bp_plot_order"] = out["bp_plot_order"].astype(int)

    if "bp_status_label" not in out.columns:
        out["bp_status_label"] = out["bioproject"]
    else:
        out["bp_status_label"] = out["bp_status_label"].fillna(out["bioproject"])

    if "status_group_for_order" not in out.columns:
        out["status_group_for_order"] = "NA"
    else:
        out["status_group_for_order"] = out["status_group_for_order"].fillna("NA")

    if "status_values" not in out.columns:
        out["status_values"] = "NA"
    else:
        out["status_values"] = out["status_values"].fillna("NA")

    return out


def summarize_by_bp(df):
    if df.empty:
        return pd.DataFrame()

    d = df.copy()
    numeric_cols = [
        "gut_total_counts",
        "gut_otu_count",
        "paired_oral_otu_count_union",
        "shared_otu_count",
        "shared_otu_count_fraction_of_gut",
        "shared_gut_counts",
        "shared_gut_abundance_fraction",
    ]

    for c in numeric_cols:
        d[c] = pd.to_numeric(d[c], errors="coerce")

    out = (
        d.groupby(["dataset_label", "bioproject"], dropna=False)
        .agg(
            n_gut_samples=("gut_sample_id", "nunique"),
            n_pair_ids=("pair_id", "nunique"),
            n_total_rows=("gut_sample_id", "size"),
            median_gut_otu_count=("gut_otu_count", "median"),
            mean_gut_otu_count=("gut_otu_count", "mean"),
            median_paired_oral_otu_count_union=("paired_oral_otu_count_union", "median"),
            mean_paired_oral_otu_count_union=("paired_oral_otu_count_union", "mean"),
            median_shared_otu_count=("shared_otu_count", "median"),
            mean_shared_otu_count=("shared_otu_count", "mean"),
            median_shared_otu_count_fraction_of_gut=("shared_otu_count_fraction_of_gut", "median"),
            mean_shared_otu_count_fraction_of_gut=("shared_otu_count_fraction_of_gut", "mean"),
            median_shared_gut_abundance_fraction=("shared_gut_abundance_fraction", "median"),
            mean_shared_gut_abundance_fraction=("shared_gut_abundance_fraction", "mean"),
            status_values=("gut_status_cleaned", safe_unique_join),
            disease_study_values=("gut_disease_study_cleaned", safe_unique_join),
        )
        .reset_index()
    )

    return out


def summarize_pooled(df):
    if df.empty:
        return pd.DataFrame()

    d = df.copy()
    numeric_cols = [
        "gut_otu_count",
        "paired_oral_otu_count_union",
        "shared_otu_count",
        "shared_otu_count_fraction_of_gut",
        "shared_gut_abundance_fraction",
    ]

    for c in numeric_cols:
        d[c] = pd.to_numeric(d[c], errors="coerce")

    out = (
        d.groupby(["dataset_label"], dropna=False)
        .agg(
            n_bioprojects=("bioproject", "nunique"),
            n_gut_samples=("gut_sample_id", "nunique"),
            n_pair_ids=("pair_id", "nunique"),
            median_gut_otu_count=("gut_otu_count", "median"),
            median_paired_oral_otu_count_union=("paired_oral_otu_count_union", "median"),
            median_shared_otu_count=("shared_otu_count", "median"),
            mean_shared_otu_count=("shared_otu_count", "mean"),
            median_shared_otu_count_fraction_of_gut=("shared_otu_count_fraction_of_gut", "median"),
            mean_shared_otu_count_fraction_of_gut=("shared_otu_count_fraction_of_gut", "mean"),
            median_shared_gut_abundance_fraction=("shared_gut_abundance_fraction", "median"),
            mean_shared_gut_abundance_fraction=("shared_gut_abundance_fraction", "mean"),
        )
        .reset_index()
    )

    return out


def build_long_plot_table(df):
    rows = []

    metrics = [
        ("shared_otu_count", "Shared OTU count", "count"),
        ("shared_otu_count_fraction_of_gut", "Shared OTU count fraction of gut", "fraction"),
        ("shared_gut_abundance_fraction", "Shared gut abundance fraction", "fraction"),
    ]

    for _, r in df.iterrows():
        for col, label, value_type in metrics:
            rows.append(
                {
                    "bioproject": r["bioproject"],
                    "dataset_label": r["dataset_label"],
                    "pair_id": r["pair_id"],
                    "gut_sample_id": r["gut_sample_id"],
                    "metric": col,
                    "metric_label": label,
                    "value_type": value_type,
                    "value": r[col],
                    "bp_plot_order": r.get("bp_plot_order", np.nan),
                    "bp_status_label": r.get("bp_status_label", r["bioproject"]),
                    "status_group_for_order": r.get("status_group_for_order", "NA"),
                    "gut_status_cleaned": r.get("gut_status_cleaned", "NA"),
                    "gut_disease_study_cleaned": r.get("gut_disease_study_cleaned", "NA"),
                }
            )

    return pd.DataFrame(rows)


def plot_per_bp_boxplot(long_df, dataset_label, metric, out_png, width_per_bp=0.72, height=6.2, x_rotation=30, point_alpha=0.35):
    d = long_df.loc[(long_df["dataset_label"] == dataset_label) & (long_df["metric"] == metric), :].copy()

    if d.empty:
        print(f"WARNING: no data for plot {dataset_label} / {metric}")
        return None

    d["value"] = pd.to_numeric(d["value"], errors="coerce")
    d = d.dropna(subset=["value"]).copy()

    if d.empty:
        print(f"WARNING: no numeric values for plot {dataset_label} / {metric}")
        return None

    bp_info = (
        d[["bioproject", "bp_plot_order", "bp_status_label"]]
        .drop_duplicates()
        .sort_values(["bp_plot_order", "bioproject"])
    )

    bps = bp_info["bioproject"].tolist()
    labels = bp_info["bp_status_label"].tolist()

    data = []
    positions = []
    rng = np.random.default_rng(2026)

    fig_width = max(10, width_per_bp * len(bps))
    fig, ax = plt.subplots(figsize=(fig_width, height))

    for i, bp in enumerate(bps):
        vals = d.loc[d["bioproject"] == bp, "value"].astype(float).values
        if len(vals) == 0:
            continue

        data.append(vals)
        positions.append(i)

        xj = i + rng.normal(0, 0.06, size=len(vals))
        ax.scatter(xj, vals, s=10, alpha=point_alpha, color="black", edgecolors="none", zorder=3)

    box = ax.boxplot(
        data,
        positions=positions,
        widths=0.45,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.1},
        boxprops={"linewidth": 0.8, "color": "black"},
        whiskerprops={"linewidth": 0.8, "color": "black"},
        capprops={"linewidth": 0.8, "color": "black"},
    )

    for patch in box["boxes"]:
        patch.set_facecolor("#7AA6DC")
        patch.set_alpha(0.85)

    ax.set_xticks(np.arange(len(bps)))
    ax.set_xticklabels(labels, rotation=x_rotation, ha="right", fontsize=8)

    metric_label = d["metric_label"].iloc[0]
    ax.set_ylabel(metric_label)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("BioProject(status_cleaned)")
    ax.set_title(f"{metric_label} by BioProject\nDataset: {dataset_label}; unit = one Gut sample", fontsize=12)
    ax.grid(axis="y", linestyle=":", linewidth=0.5, alpha=0.55)

    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)

    return str(out_png)


def plot_pooled_boxplot(long_df, dataset_label, metric, out_png, point_alpha=0.35):
    d = long_df.loc[(long_df["dataset_label"] == dataset_label) & (long_df["metric"] == metric), :].copy()

    if d.empty:
        return None

    d["value"] = pd.to_numeric(d["value"], errors="coerce")
    d = d.dropna(subset=["value"]).copy()

    if d.empty:
        return None

    vals = d["value"].astype(float).values
    metric_label = d["metric_label"].iloc[0]

    fig, ax = plt.subplots(figsize=(4.8, 5.8))

    box = ax.boxplot(
        [vals],
        positions=[0],
        widths=0.45,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.1},
        boxprops={"linewidth": 0.8, "color": "black"},
        whiskerprops={"linewidth": 0.8, "color": "black"},
        capprops={"linewidth": 0.8, "color": "black"},
    )

    for patch in box["boxes"]:
        patch.set_facecolor("#7AA6DC")
        patch.set_alpha(0.85)

    rng = np.random.default_rng(2026)
    xj = rng.normal(0, 0.045, size=len(vals))
    ax.scatter(xj, vals, s=9, alpha=point_alpha, color="black", edgecolors="none", zorder=3)

    ax.set_xticks([0])
    ax.set_xticklabels([dataset_label])
    ax.set_ylabel(metric_label)
    ax.set_title(f"Pooled {metric_label}\n{dataset_label}; n gut samples = {len(vals)}", fontsize=11)
    ax.set_ylim(bottom=0)
    ax.grid(axis="y", linestyle=":", linewidth=0.5, alpha=0.55)

    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)

    return str(out_png)


def main():
    parser = argparse.ArgumentParser(description="Calculate paired oral-to-gut shared OTU overlap from global filter.biom.")

    parser.add_argument("--metadata", default=str(DEFAULT_METADATA))
    parser.add_argument("--biom", default=str(DEFAULT_BIOM))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))

    parser.add_argument("--sample-col", default="sample_id")
    parser.add_argument("--bioproject-col", default="bioproject")
    parser.add_argument("--site-col", default="site_group_cleaned")
    parser.add_argument("--pair-col", default="pair_id")
    parser.add_argument("--subject-col", default="subject_id")
    parser.add_argument("--status-col", default="status_cleaned")
    parser.add_argument("--disease-col", default="disease_study_cleaned")

    parser.add_argument("--dataset-label", default="filter")
    parser.add_argument("--presence-threshold", type=float, default=0.0, help="Feature is present if count > threshold. Default: 0.")

    parser.add_argument("--width-per-bp", type=float, default=0.72)
    parser.add_argument("--plot-height", type=float, default=6.2)
    parser.add_argument("--x-rotation", type=float, default=30)
    parser.add_argument("--point-alpha", type=float, default=0.35)

    parser.add_argument("--overwrite", action="store_true")

    args = parser.parse_args()

    metadata_fp = Path(args.metadata)
    biom_fp = Path(args.biom)
    out_dir = Path(args.out_dir)
    plot_dir = out_dir / "plots"

    out_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    outputs = {
        "per_sample": out_dir / "shared_otu_oral_to_gut_per_sample.tsv",
        "summary_by_bp": out_dir / "shared_otu_oral_to_gut_summary_by_bp.tsv",
        "summary_pooled": out_dir / "shared_otu_oral_to_gut_summary_pooled.tsv",
        "long_plot_table": out_dir / "shared_otu_oral_to_gut_long_plot_table.tsv",
        "summary_json": out_dir / "shared_otu_oral_to_gut_summary.json",
    }

    for p in outputs.values():
        if p.exists() and not args.overwrite:
            raise FileExistsError(f"Output exists: {p}. Use --overwrite.")

    print("Resolved paths:", flush=True)
    print(f"  metadata: {metadata_fp}", flush=True)
    print(f"  biom:     {biom_fp}", flush=True)
    print(f"  out_dir:  {out_dir}", flush=True)
    print(f"  plots:    {plot_dir}", flush=True)

    if not metadata_fp.exists():
        raise FileNotFoundError(f"Metadata not found: {metadata_fp}")
    if not biom_fp.exists():
        raise FileNotFoundError(f"BIOM not found: {biom_fp}")

    print("[1/6] Reading metadata ...", flush=True)
    metadata_raw = read_table(metadata_fp)

    metadata = prepare_metadata(
        metadata=metadata_raw,
        sample_col=args.sample_col,
        bioproject_col=args.bioproject_col,
        site_col=args.site_col,
        pair_col=args.pair_col,
        subject_col=args.subject_col,
        status_col=args.status_col,
        disease_col=args.disease_col,
    )

    bp_order = build_bp_order_from_metadata(metadata, bioproject_col=args.bioproject_col, status_col=args.status_col)

    print("[2/6] Reading BIOM ...", flush=True)
    X, sample_ids, feature_ids = load_biom_sparse_sample_by_feature(biom_fp)

    print(f"  BIOM samples: {len(sample_ids)}", flush=True)
    print(f"  BIOM features: {len(feature_ids)}", flush=True)
    print(f"  BIOM nonzero entries: {X.nnz}", flush=True)

    print("[3/6] Calculating shared OTU overlap ...", flush=True)
    per_sample = calculate_shared_otu(
        X=X,
        sample_ids=sample_ids,
        feature_ids=feature_ids,
        metadata=metadata,
        dataset_label=args.dataset_label,
        input_biom=biom_fp,
        sample_col=args.sample_col,
        bioproject_col=args.bioproject_col,
        pair_col=args.pair_col,
        subject_col=args.subject_col,
        status_col=args.status_col,
        disease_col=args.disease_col,
        presence_threshold=args.presence_threshold,
    )

    if per_sample.empty:
        raise RuntimeError("No shared OTU results were produced. Check pair_id/site metadata and sample overlap with BIOM.")

    print("[4/6] Annotating and summarizing ...", flush=True)
    per_sample = annotate_bp_order(per_sample, bp_order)

    summary_by_bp = summarize_by_bp(per_sample)
    summary_by_bp = annotate_bp_order(summary_by_bp, bp_order)

    summary_pooled = summarize_pooled(per_sample)
    long_plot = build_long_plot_table(per_sample)

    print("[5/6] Plotting ...", flush=True)
    plot_files = []
    plot_metrics = [
        "shared_otu_count",
        "shared_otu_count_fraction_of_gut",
        "shared_gut_abundance_fraction",
    ]

    for metric in plot_metrics:
        per_bp_png = plot_dir / f"per_bp_{clean_for_filename(metric)}_{clean_for_filename(args.dataset_label)}.png"
        pooled_png = plot_dir / f"pooled_{clean_for_filename(metric)}_{clean_for_filename(args.dataset_label)}.png"

        if per_bp_png.exists() and not args.overwrite:
            raise FileExistsError(f"Output exists: {per_bp_png}. Use --overwrite.")
        if pooled_png.exists() and not args.overwrite:
            raise FileExistsError(f"Output exists: {pooled_png}. Use --overwrite.")

        p1 = plot_per_bp_boxplot(
            long_df=long_plot,
            dataset_label=args.dataset_label,
            metric=metric,
            out_png=per_bp_png,
            width_per_bp=args.width_per_bp,
            height=args.plot_height,
            x_rotation=args.x_rotation,
            point_alpha=args.point_alpha,
        )
        if p1 is not None:
            plot_files.append(p1)

        p2 = plot_pooled_boxplot(
            long_df=long_plot,
            dataset_label=args.dataset_label,
            metric=metric,
            out_png=pooled_png,
            point_alpha=args.point_alpha,
        )
        if p2 is not None:
            plot_files.append(p2)

    print("[6/6] Writing outputs ...", flush=True)
    per_sample.to_csv(outputs["per_sample"], sep="\t", index=False)
    summary_by_bp.to_csv(outputs["summary_by_bp"], sep="\t", index=False)
    summary_pooled.to_csv(outputs["summary_pooled"], sep="\t", index=False)
    long_plot.to_csv(outputs["long_plot_table"], sep="\t", index=False)

    with open(outputs["summary_json"], "w") as f:
        json.dump(
            {
                "inputs": {
                    "metadata": os.path.abspath(str(metadata_fp)),
                    "biom": os.path.abspath(str(biom_fp)),
                },
                "out_dir": os.path.abspath(str(out_dir)),
                "plot_dir": os.path.abspath(str(plot_dir)),
                "dataset_label": args.dataset_label,
                "presence_threshold": args.presence_threshold,
                "definition": {
                    "paired_oral_otu_set": "Union of OTUs with count > presence_threshold across paired Oral sample(s) from the same BioProject and pair_id.",
                    "shared_otu_count": "Number of Gut OTUs also present in paired Oral OTU set.",
                    "shared_otu_count_fraction_of_gut": "shared_otu_count divided by number of OTUs detected in the Gut sample.",
                    "shared_gut_abundance_fraction": "Sum of Gut counts for shared OTUs divided by total Gut counts.",
                    "important_note": "This is shared OTU overlap, not proof of oral origin.",
                },
                "outputs": {k: str(v) for k, v in outputs.items()},
                "plots": plot_files,
                "n_rows_per_sample": int(per_sample.shape[0]),
                "n_bioprojects": int(per_sample["bioproject"].nunique()),
                "n_gut_samples": int(per_sample["gut_sample_id"].nunique()),
                "summary_pooled": json_safe_records(summary_pooled),
            },
            f,
            indent=2,
        )

    print("\nDone.", flush=True)
    print("Main outputs:", flush=True)
    for k, v in outputs.items():
        print(f"  {k}: {v}", flush=True)

    print("\nPlots:", flush=True)
    for p in plot_files:
        print(f"  {p}", flush=True)

    print("\nPooled summary:", flush=True)
    print(summary_pooled.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
