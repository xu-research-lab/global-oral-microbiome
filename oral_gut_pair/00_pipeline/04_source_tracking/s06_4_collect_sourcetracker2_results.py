#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
s06_4_collect_sourcetracker2_results.py

Collect SourceTracker2 task outputs from st2_tasks_manifest.tsv.

Expected input:
  04.sourcetracking/sourcetracker/st2_tasks_manifest.tsv

Expected ST2 output file per task:
  <st2_out_dir>/mixing_proportions.txt

Outputs:
  paired_oral_fraction.tsv
  sourcetracker2_task_status.tsv
  sourcetracker2_summary_by_bp.tsv
  sourcetracker2_summary_pooled.tsv
  sourcetracker2_collect_summary.json
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_ROOT = Path("/beegfs/zhangcaiyi/oral_gut_pair_noTimeSeries_32bp")
DEFAULT_OUT_DIR = DEFAULT_ROOT / "04.sourcetracking" / "sourcetracker"
DEFAULT_MANIFEST = DEFAULT_OUT_DIR / "st2_tasks_manifest.tsv"


def read_table(path):
    path = str(path)
    if path.endswith((".tsv", ".txt")):
        return pd.read_csv(path, sep="\t", dtype=str)
    return pd.read_csv(path, sep=None, engine="python", dtype=str)


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


def read_mixing_proportions(fp, gut_sample_id, source_env_name="paired_oral"):
    """
    Robustly read SourceTracker2 mixing_proportions.txt.

    Most ST2 outputs are one row per sink and one column per source category,
    often including paired_oral and Unknown. This function supports either:
      - explicit first sample-id column, or
      - sample IDs in the index-like first unnamed column, or
      - a single-row table without an obvious sample column.
    """
    df0 = pd.read_csv(fp, sep="\t")

    if df0.empty:
        return {}, "empty_mixing_proportions"

    df = df0.copy()

    possible_id_cols = ["#SampleID", "SampleID", "sample_id", "sample", "Sink", "sink", "SinkID", "sink_id"]
    id_col = None
    for c in possible_id_cols:
        if c in df.columns:
            id_col = c
            break

    if id_col is None:
        first_col = df.columns[0]
        first_numeric = pd.to_numeric(df[first_col], errors="coerce").notna().mean()
        if first_numeric < 0.5:
            id_col = first_col

    if id_col is not None:
        df[id_col] = df[id_col].astype(str)
        if str(gut_sample_id) in set(df[id_col]):
            row = df.loc[df[id_col].astype(str).eq(str(gut_sample_id)), :].iloc[0].copy()
        elif df.shape[0] == 1:
            row = df.iloc[0].copy()
        else:
            return {}, f"gut_sample_not_found_in_mixing_proportions: {gut_sample_id}"
        row = row.drop(labels=[id_col], errors="ignore")
    else:
        if df.shape[0] != 1:
            return {}, "cannot_identify_sink_row"
        row = df.iloc[0].copy()

    props = {}
    for k, v in row.items():
        k2 = str(k).strip()
        if k2 == "" or k2.lower().startswith("unnamed"):
            continue
        val = pd.to_numeric(pd.Series([v]), errors="coerce").iloc[0]
        if pd.notna(val):
            props[k2] = float(val)

    if not props:
        return {}, "no_numeric_proportion_columns"

    return props, "ok"


def summarize_by_bp(df):
    if df.empty:
        return pd.DataFrame()

    d = df.copy()
    for c in ["paired_oral_fraction", "unknown_fraction"]:
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")

    out = (
        d.groupby(["dataset_label", "bioproject"], dropna=False)
        .agg(
            n_tasks=("task_id", "size"),
            n_tasks_ok=("task_status", lambda x: int((x == "ok").sum())),
            n_gut_samples=("gut_sample_id", "nunique"),
            n_pair_ids=("pair_id", "nunique"),
            median_paired_oral_fraction=("paired_oral_fraction", "median"),
            mean_paired_oral_fraction=("paired_oral_fraction", "mean"),
            median_unknown_fraction=("unknown_fraction", "median"),
            mean_unknown_fraction=("unknown_fraction", "mean"),
        )
        .reset_index()
    )

    out["fraction_tasks_ok"] = out["n_tasks_ok"] / out["n_tasks"].replace(0, np.nan)
    return out


def summarize_pooled(df):
    if df.empty:
        return pd.DataFrame()

    d = df.copy()
    for c in ["paired_oral_fraction", "unknown_fraction"]:
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")

    out = (
        d.groupby(["dataset_label"], dropna=False)
        .agg(
            n_bioprojects=("bioproject", "nunique"),
            n_tasks=("task_id", "size"),
            n_tasks_ok=("task_status", lambda x: int((x == "ok").sum())),
            n_gut_samples=("gut_sample_id", "nunique"),
            n_pair_ids=("pair_id", "nunique"),
            median_paired_oral_fraction=("paired_oral_fraction", "median"),
            mean_paired_oral_fraction=("paired_oral_fraction", "mean"),
            median_unknown_fraction=("unknown_fraction", "median"),
            mean_unknown_fraction=("unknown_fraction", "mean"),
        )
        .reset_index()
    )
    out["fraction_tasks_ok"] = out["n_tasks_ok"] / out["n_tasks"].replace(0, np.nan)
    return out


def main():
    parser = argparse.ArgumentParser(description="Collect SourceTracker2 paired oral fraction results.")
    parser.add_argument("--manifest", default=str(DEFAULT_MANIFEST))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--source-env-name", default="paired_oral")
    parser.add_argument("--unknown-column", default="Unknown")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    manifest_fp = Path(args.manifest).resolve(strict=False)
    out_dir = Path(args.out_dir).resolve(strict=False)
    out_dir.mkdir(parents=True, exist_ok=True)

    outputs = {
        "paired_oral_fraction": out_dir / "paired_oral_fraction.tsv",
        "task_status": out_dir / "sourcetracker2_task_status.tsv",
        "summary_by_bp": out_dir / "sourcetracker2_summary_by_bp.tsv",
        "summary_pooled": out_dir / "sourcetracker2_summary_pooled.tsv",
        "summary_json": out_dir / "sourcetracker2_collect_summary.json",
    }

    for p in outputs.values():
        if p.exists() and not args.overwrite:
            raise FileExistsError(f"Output exists: {p}. Use --overwrite.")

    print("Resolved paths:", flush=True)
    print(f"  manifest: {manifest_fp}", flush=True)
    print(f"  out_dir:  {out_dir}", flush=True)

    if not manifest_fp.exists():
        raise FileNotFoundError(f"Manifest not found: {manifest_fp}")

    manifest = read_table(manifest_fp)
    required = ["task_id", "dataset_label", "bioproject", "pair_id", "gut_sample_id", "oral_sample_ids", "st2_out_dir"]
    missing = [c for c in required if c not in manifest.columns]
    if missing:
        raise ValueError(
            "Manifest missing required columns:\n  "
            + "\n  ".join(missing)
            + "\n\nAvailable columns:\n  "
            + "\n  ".join(manifest.columns)
        )

    rows = []

    for _, r in manifest.iterrows():
        st2_out_dir = Path(str(r["st2_out_dir"]))
        mix_fp = st2_out_dir / "mixing_proportions.txt"

        base = r.to_dict()
        base["st2_out_dir"] = str(st2_out_dir)
        base["mixing_proportions"] = str(mix_fp)

        if not mix_fp.exists():
            base.update({
                "task_status": "missing_mixing_proportions",
                "task_status_detail": "mixing_proportions.txt_not_found",
                "paired_oral_fraction": np.nan,
                "unknown_fraction": np.nan,
            })
            rows.append(base)
            continue

        props, status = read_mixing_proportions(
            mix_fp,
            gut_sample_id=str(r["gut_sample_id"]),
            source_env_name=args.source_env_name,
        )

        base["task_status"] = status
        base["task_status_detail"] = status

        paired_val = props.get(args.source_env_name, np.nan)
        if pd.isna(paired_val):
            # Case-insensitive fallback.
            for k, v in props.items():
                if str(k).lower() == args.source_env_name.lower():
                    paired_val = v
                    break

        unknown_val = props.get(args.unknown_column, np.nan)
        if pd.isna(unknown_val):
            for k, v in props.items():
                if str(k).lower() == args.unknown_column.lower():
                    unknown_val = v
                    break

        base["paired_oral_fraction"] = paired_val
        base["unknown_fraction"] = unknown_val

        for k, v in props.items():
            safe_k = str(k).replace(" ", "_").replace("/", "_")
            base[f"proportion__{safe_k}"] = v

        rows.append(base)

    result = pd.DataFrame(rows)

    # If a file was parsed successfully but the paired_oral column was missing,
    # keep task_status as parsed status but make the missing value obvious.
    result["paired_oral_fraction"] = pd.to_numeric(result["paired_oral_fraction"], errors="coerce")
    result["unknown_fraction"] = pd.to_numeric(result["unknown_fraction"], errors="coerce")

    summary_by_bp = summarize_by_bp(result)
    summary_pooled = summarize_pooled(result)

    result.to_csv(outputs["paired_oral_fraction"], sep="\t", index=False)
    result.to_csv(outputs["task_status"], sep="\t", index=False)
    summary_by_bp.to_csv(outputs["summary_by_bp"], sep="\t", index=False)
    summary_pooled.to_csv(outputs["summary_pooled"], sep="\t", index=False)

    with open(outputs["summary_json"], "w") as f:
        json.dump(
            {
                "manifest": str(manifest_fp),
                "out_dir": str(out_dir),
                "source_env_name": args.source_env_name,
                "unknown_column": args.unknown_column,
                "outputs": {k: str(v) for k, v in outputs.items()},
                "n_tasks_in_manifest": int(manifest.shape[0]),
                "n_rows_collected": int(result.shape[0]),
                "n_tasks_ok": int((result["task_status"] == "ok").sum()),
                "n_tasks_missing_mixing_proportions": int((result["task_status"] == "missing_mixing_proportions").sum()),
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

    print("\nTask status counts:", flush=True)
    print(result["task_status"].value_counts(dropna=False).to_string(), flush=True)


if __name__ == "__main__":
    main()
