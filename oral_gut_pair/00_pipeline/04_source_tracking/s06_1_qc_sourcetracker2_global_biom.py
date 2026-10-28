#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
s06_1_qc_sourcetracker2_global_biom.py

Read-depth QC before SourceTracker2 using the current global filtered BIOM table.

Default inputs:
  01.filter_data/filter.biom
  01.filter_data/filter_metadata.csv

Default output directory:
  04.sourcetracking/sourcetracker

Filtering rule:
  Pair-level filtering. For each BioProject x pair_id, if any Oral/Gut sample
  in that pair has read_count < min_reads, remove the whole pair. Otherwise,
  keep the whole pair if it has at least one Oral and one Gut sample.

Main outputs:
  sample_read_depth_all.tsv
  pair_depth_qc_all.tsv
  sourcetracker2_depth_qc_summary_by_bp.tsv
  sourcetracker2_depth_qc_summary_pooled.tsv
  samples_keep_for_sourcetracker2.tsv
  tasks_keep_for_sourcetracker2.tsv
  sourcetracker2_depth_qc_summary.json
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_ROOT = Path("/beegfs/zhangcaiyi/oral_gut_pair_noTimeSeries_32bp")
DEFAULT_BIOM = DEFAULT_ROOT / "01.filter_data" / "filter.biom"
DEFAULT_METADATA = DEFAULT_ROOT / "01.filter_data" / "filter_metadata.csv"
DEFAULT_OUT_DIR = DEFAULT_ROOT / "04.sourcetracking" / "sourcetracker"


def read_table(path):
    path = str(path)
    if path.endswith((".tsv", ".txt")):
        return pd.read_csv(path, sep="\t", dtype=str)
    return pd.read_csv(path, sep=None, engine="python", dtype=str)


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


def safe_join(x):
    vals = []
    for v in x:
        if pd.isna(v):
            continue
        v = str(v).strip()
        if v:
            vals.append(v)
    vals = sorted(set(vals))
    return ";".join(vals) if vals else ""


def validate_columns(df, required, label):
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"{label} missing required columns:\n  "
            + "\n  ".join(missing)
            + "\n\nAvailable columns:\n  "
            + "\n  ".join(df.columns)
        )


def load_biom_sample_depth(path):
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
    sample_ids = [str(x) for x in table.ids(axis="sample")]
    mat = table.matrix_data.tocsr()  # observation x sample
    sample_depth = np.asarray(mat.sum(axis=0)).ravel()

    return pd.DataFrame({"sample_id": sample_ids, "read_count": sample_depth.astype(float)})


def prepare_metadata(metadata, sample_col, bioproject_col, site_col, pair_col, subject_col, status_col, disease_col):
    md = metadata.copy()

    required = [sample_col, bioproject_col, site_col, pair_col]
    validate_columns(md, required, "metadata")

    for c in [subject_col, status_col, disease_col]:
        if c not in md.columns:
            md[c] = "NA"

    for c in [sample_col, bioproject_col, site_col, pair_col, subject_col, status_col, disease_col]:
        md[c] = md[c].astype(str).str.strip()

    md = md.rename(
        columns={
            sample_col: "sample_id",
            bioproject_col: "bioproject",
            site_col: "site_original",
            pair_col: "pair_id",
            subject_col: "subject_id",
            status_col: "status_cleaned",
            disease_col: "disease_study_cleaned",
        }
    )

    md["site"] = md["site_original"].apply(standardize_site)
    md = md.loc[md["site"].isin(["Oral", "Gut"]), :].copy()
    md = md.drop_duplicates(subset=["sample_id"]).copy()

    return md


def build_pair_qc(sample_depth, min_reads):
    d = sample_depth.copy()
    d["read_count"] = pd.to_numeric(d["read_count"], errors="coerce")
    d["low_depth"] = d["read_count"] < float(min_reads)

    pair_rows = []

    for (bp, pair_id), sub in d.groupby(["bioproject", "pair_id"], dropna=False):
        n_samples = int(sub.shape[0])
        n_oral = int((sub["site"] == "Oral").sum())
        n_gut = int((sub["site"] == "Gut").sum())
        n_low = int(sub["low_depth"].sum())

        has_low = n_low > 0
        has_complete_oral_gut = n_oral > 0 and n_gut > 0
        keep_pair = bool((not has_low) and has_complete_oral_gut)

        if keep_pair:
            reason = "keep"
        elif has_low:
            reason = "has_low_depth_sample"
        else:
            reason = "missing_oral_or_gut"

        pair_rows.append(
            {
                "dataset_label": "filter",
                "bioproject": bp,
                "pair_id": pair_id,
                "keep_pair_for_sourcetracker2": keep_pair,
                "remove_pair_reason": reason,
                "n_samples_in_pair": n_samples,
                "n_oral_samples": n_oral,
                "n_gut_samples": n_gut,
                "n_low_depth_samples": n_low,
                "min_read_count_in_pair": float(sub["read_count"].min()),
                "median_read_count_in_pair": float(sub["read_count"].median()),
                "max_read_count_in_pair": float(sub["read_count"].max()),
                "low_depth_sample_ids": safe_join(sub.loc[sub["low_depth"], "sample_id"]),
                "oral_sample_ids": safe_join(sub.loc[sub["site"] == "Oral", "sample_id"]),
                "gut_sample_ids": safe_join(sub.loc[sub["site"] == "Gut", "sample_id"]),
                "status_values": safe_join(sub["status_cleaned"]),
                "disease_study_values": safe_join(sub["disease_study_cleaned"]),
            }
        )

    pair_qc = pd.DataFrame(pair_rows)

    d = d.merge(
        pair_qc[["bioproject", "pair_id", "keep_pair_for_sourcetracker2", "remove_pair_reason"]],
        on=["bioproject", "pair_id"],
        how="left",
    )
    d["dataset_label"] = "filter"
    d["keep_sample_for_sourcetracker2"] = d["keep_pair_for_sourcetracker2"]

    return d, pair_qc


def summarize_by_bp(sample_depth, pair_qc):
    s = sample_depth.copy()
    p = pair_qc.copy()

    sample_summary = (
        s.groupby(["dataset_label", "bioproject"], dropna=False)
        .agg(
            n_samples_total=("sample_id", "nunique"),
            n_oral_samples_total=("site", lambda x: int((x == "Oral").sum())),
            n_gut_samples_total=("site", lambda x: int((x == "Gut").sum())),
            n_samples_low_depth=("low_depth", "sum"),
            n_samples_keep=("keep_sample_for_sourcetracker2", "sum"),
            median_read_count=("read_count", "median"),
            min_read_count=("read_count", "min"),
        )
        .reset_index()
    )

    pair_summary = (
        p.groupby(["dataset_label", "bioproject"], dropna=False)
        .agg(
            n_pairs_total=("pair_id", "nunique"),
            n_pairs_keep=("keep_pair_for_sourcetracker2", "sum"),
            n_pairs_removed=("keep_pair_for_sourcetracker2", lambda x: int((~x.astype(bool)).sum())),
            n_pairs_removed_low_depth=("remove_pair_reason", lambda x: int((x == "has_low_depth_sample").sum())),
            n_pairs_removed_missing_oral_or_gut=("remove_pair_reason", lambda x: int((x == "missing_oral_or_gut").sum())),
            n_gut_sink_tasks_keep=("n_gut_samples", lambda x: int(x[p.loc[x.index, "keep_pair_for_sourcetracker2"].astype(bool)].sum())),
        )
        .reset_index()
    )

    out = sample_summary.merge(pair_summary, on=["dataset_label", "bioproject"], how="outer")

    count_cols = [
        "n_samples_total",
        "n_oral_samples_total",
        "n_gut_samples_total",
        "n_samples_low_depth",
        "n_samples_keep",
        "n_pairs_total",
        "n_pairs_keep",
        "n_pairs_removed",
        "n_pairs_removed_low_depth",
        "n_pairs_removed_missing_oral_or_gut",
        "n_gut_sink_tasks_keep",
    ]
    for c in count_cols:
        out[c] = pd.to_numeric(out[c], errors="coerce").fillna(0).astype(int)

    out["fraction_samples_keep"] = out["n_samples_keep"] / out["n_samples_total"].replace(0, np.nan)
    out["fraction_pairs_keep"] = out["n_pairs_keep"] / out["n_pairs_total"].replace(0, np.nan)

    return out


def summarize_pooled(sample_depth, pair_qc):
    s = sample_depth.copy()
    p = pair_qc.copy()

    n_samples_total = int(s["sample_id"].nunique())
    n_samples_keep = int(s.loc[s["keep_sample_for_sourcetracker2"].astype(bool), "sample_id"].nunique())
    n_pairs_total = int(p["pair_id"].nunique())
    n_pairs_keep = int(p["keep_pair_for_sourcetracker2"].sum())
    n_pairs_removed = int((~p["keep_pair_for_sourcetracker2"].astype(bool)).sum())
    n_gut_sink_tasks_keep = int(p.loc[p["keep_pair_for_sourcetracker2"].astype(bool), "n_gut_samples"].sum())

    return pd.DataFrame(
        [
            {
                "dataset_label": "filter",
                "n_bioprojects": int(s["bioproject"].nunique()),
                "n_samples_total": n_samples_total,
                "n_oral_samples_total": int((s["site"] == "Oral").sum()),
                "n_gut_samples_total": int((s["site"] == "Gut").sum()),
                "n_samples_low_depth": int(s["low_depth"].sum()),
                "n_samples_keep": n_samples_keep,
                "n_samples_removed": n_samples_total - n_samples_keep,
                "fraction_samples_keep": n_samples_keep / n_samples_total if n_samples_total else np.nan,
                "n_pairs_total": n_pairs_total,
                "n_pairs_keep": n_pairs_keep,
                "n_pairs_removed": n_pairs_removed,
                "fraction_pairs_keep": n_pairs_keep / n_pairs_total if n_pairs_total else np.nan,
                "n_gut_sink_tasks_keep": n_gut_sink_tasks_keep,
                "min_read_count": float(s["read_count"].min()),
                "q01_read_count": float(s["read_count"].quantile(0.01)),
                "q05_read_count": float(s["read_count"].quantile(0.05)),
                "median_read_count": float(s["read_count"].median()),
            }
        ]
    )


def build_keep_tasks(pair_qc):
    keep = pair_qc.loc[pair_qc["keep_pair_for_sourcetracker2"].astype(bool), :].copy()
    rows = []

    for _, r in keep.iterrows():
        gut_ids = [x for x in str(r["gut_sample_ids"]).split(";") if x]
        oral_ids = [x for x in str(r["oral_sample_ids"]).split(";") if x]

        for gut_id in gut_ids:
            rows.append(
                {
                    "dataset_label": "filter",
                    "bioproject": r["bioproject"],
                    "pair_id": r["pair_id"],
                    "gut_sample_id": gut_id,
                    "oral_sample_ids": ";".join(oral_ids),
                    "n_oral_sources": len(oral_ids),
                    "filter_rule": "pair_level_remove_if_any_sample_read_count_lt_min_reads",
                }
            )

    return pd.DataFrame(rows)


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


def main():
    parser = argparse.ArgumentParser(description="Pair-level read-depth QC before SourceTracker2 from global filter.biom.")

    parser.add_argument("--biom", default=str(DEFAULT_BIOM))
    parser.add_argument("--metadata", default=str(DEFAULT_METADATA))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--min-reads", type=float, default=1000)

    parser.add_argument("--sample-col", default="sample_id")
    parser.add_argument("--bioproject-col", default="bioproject")
    parser.add_argument("--site-col", default="site_group_cleaned")
    parser.add_argument("--pair-col", default="pair_id")
    parser.add_argument("--subject-col", default="subject_id")
    parser.add_argument("--status-col", default="status_cleaned")
    parser.add_argument("--disease-col", default="disease_study_cleaned")

    parser.add_argument("--overwrite", action="store_true")

    args = parser.parse_args()

    biom_fp = Path(args.biom).resolve(strict=False)
    metadata_fp = Path(args.metadata).resolve(strict=False)
    out_dir = Path(args.out_dir).resolve(strict=False)
    out_dir.mkdir(parents=True, exist_ok=True)

    outputs = {
        "sample_depth": out_dir / "sample_read_depth_all.tsv",
        "pair_qc": out_dir / "pair_depth_qc_all.tsv",
        "summary_by_bp": out_dir / "sourcetracker2_depth_qc_summary_by_bp.tsv",
        "summary_pooled": out_dir / "sourcetracker2_depth_qc_summary_pooled.tsv",
        "samples_keep": out_dir / "samples_keep_for_sourcetracker2.tsv",
        "tasks_keep": out_dir / "tasks_keep_for_sourcetracker2.tsv",
        "summary_json": out_dir / "sourcetracker2_depth_qc_summary.json",
    }

    for p in outputs.values():
        if p.exists() and not args.overwrite:
            raise FileExistsError(f"Output exists: {p}. Use --overwrite.")

    print("Resolved paths:", flush=True)
    print(f"  biom:     {biom_fp}", flush=True)
    print(f"  metadata: {metadata_fp}", flush=True)
    print(f"  out_dir:  {out_dir}", flush=True)

    if not biom_fp.exists():
        raise FileNotFoundError(f"BIOM not found: {biom_fp}")
    if not metadata_fp.exists():
        raise FileNotFoundError(f"Metadata not found: {metadata_fp}")

    print("[1/5] Reading metadata ...", flush=True)
    metadata = read_table(metadata_fp)
    md = prepare_metadata(
        metadata=metadata,
        sample_col=args.sample_col,
        bioproject_col=args.bioproject_col,
        site_col=args.site_col,
        pair_col=args.pair_col,
        subject_col=args.subject_col,
        status_col=args.status_col,
        disease_col=args.disease_col,
    )

    print("[2/5] Counting reads per sample from BIOM ...", flush=True)
    depth = load_biom_sample_depth(biom_fp)

    print("[3/5] Merging metadata and applying pair-level depth filter ...", flush=True)
    sample_depth0 = depth.merge(md, on="sample_id", how="inner")

    if sample_depth0.empty:
        raise RuntimeError("No metadata samples matched BIOM sample IDs.")

    sample_depth, pair_qc = build_pair_qc(sample_depth0, min_reads=args.min_reads)

    print("[4/5] Summarizing and building ST2 task table ...", flush=True)
    summary_by_bp = summarize_by_bp(sample_depth, pair_qc)
    summary_pooled = summarize_pooled(sample_depth, pair_qc)
    samples_keep = sample_depth.loc[sample_depth["keep_sample_for_sourcetracker2"].astype(bool), :].copy()
    tasks_keep = build_keep_tasks(pair_qc)

    print("[5/5] Writing outputs ...", flush=True)
    sample_depth.to_csv(outputs["sample_depth"], sep="\t", index=False)
    pair_qc.to_csv(outputs["pair_qc"], sep="\t", index=False)
    summary_by_bp.to_csv(outputs["summary_by_bp"], sep="\t", index=False)
    summary_pooled.to_csv(outputs["summary_pooled"], sep="\t", index=False)
    samples_keep.to_csv(outputs["samples_keep"], sep="\t", index=False)
    tasks_keep.to_csv(outputs["tasks_keep"], sep="\t", index=False)

    with open(outputs["summary_json"], "w") as f:
        json.dump(
            {
                "inputs": {"biom": str(biom_fp), "metadata": str(metadata_fp)},
                "out_dir": str(out_dir),
                "dataset_label": "filter",
                "min_reads": args.min_reads,
                "filter_rule": (
                    "Pair-level filtering: if any Oral/Gut sample in a BioProject x pair_id "
                    "has read_count < min_reads, the whole pair_id is removed."
                ),
                "outputs": {k: str(v) for k, v in outputs.items()},
                "summary_pooled": json_safe_records(summary_pooled),
            },
            f,
            indent=2,
        )

    print("\nDone.", flush=True)
    print("Main outputs:", flush=True)
    for k, v in outputs.items():
        print(f"  {k}: {v}", flush=True)
    print("\nPooled summary:", flush=True)
    print(summary_pooled.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
