#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
s06_2_prepare_sourcetracker2_global_biom.py

Prepare per-Gut-sample SourceTracker2 inputs from the retained task table
produced by s06_1_qc_sourcetracker2_global_biom.py.

Design:
  For each retained Gut sample:
    source = paired Oral sample(s) from same BioProject and same pair_id
    sink   = this Gut sample only

Default inputs:
  04.sourcetracking/sourcetracker/tasks_keep_for_sourcetracker2.tsv
  01.filter_data/filter.biom

Default output directory:
  04.sourcetracking/sourcetracker

Outputs:
  st2_tasks_manifest.tsv
  st2_prepare_skipped_tasks.tsv
  st2_prepare_summary.json
  tasks/filter/<BioProject>/task_<task_id>/input.biom
  tasks/filter/<BioProject>/task_<task_id>/map.tsv
  tasks/filter/<BioProject>/task_<task_id>/st2_out/

Map columns:
  #SampleID
  SourceSink: source / sink
  Env: paired_oral for source, gut_sink for sink

Important:
  Manifest paths are absolute to make Slurm array execution robust.
"""

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_ROOT = Path("/beegfs/zhangcaiyi/oral_gut_pair_noTimeSeries_32bp")
DEFAULT_BIOM = DEFAULT_ROOT / "01.filter_data" / "filter.biom"
DEFAULT_OUT_DIR = DEFAULT_ROOT / "04.sourcetracking" / "sourcetracker"
DEFAULT_TASKS_KEEP = DEFAULT_OUT_DIR / "tasks_keep_for_sourcetracker2.tsv"


def clean_for_path(x):
    x = str(x)
    x = re.sub(r"[^\w.\-]+", "_", x)
    x = x.strip("_")
    return x if x else "NA"


def read_table(path):
    path = str(path)
    if path.endswith((".tsv", ".txt")):
        return pd.read_csv(path, sep="\t", dtype=str)
    return pd.read_csv(path, sep=None, engine="python", dtype=str)


def parse_oral_sample_ids(x):
    if pd.isna(x):
        return []
    return [i.strip() for i in str(x).split(";") if i.strip()]


def validate_columns(df, required, label):
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"{label} missing required columns:\n  "
            + "\n  ".join(missing)
            + "\n\nAvailable columns:\n  "
            + "\n  ".join(df.columns)
        )


def load_biom(path):
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


def write_biom(table, out_fp):
    try:
        from biom.util import biom_open
    except ImportError as e:
        raise ImportError("Cannot import biom.util.biom_open.") from e

    with biom_open(str(out_fp), "w") as f:
        table.to_hdf5(f, generated_by="s06_2_prepare_sourcetracker2_global_biom.py")


def subset_biom_for_task(table, sample_ids):
    keep = set(map(str, sample_ids))

    sub = table.filter(
        lambda vals, sid, md: str(sid) in keep,
        axis="sample",
        inplace=False,
    )

    sub = sub.filter(
        lambda vals, oid, md: np.asarray(vals).sum() > 0,
        axis="observation",
        inplace=False,
    )

    return sub


def main():
    parser = argparse.ArgumentParser(description="Prepare paired-oral-only SourceTracker2 tasks from global filter.biom.")

    parser.add_argument("--tasks-keep", default=str(DEFAULT_TASKS_KEEP))
    parser.add_argument("--biom", default=str(DEFAULT_BIOM))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))

    parser.add_argument("--source-env-name", default="paired_oral")
    parser.add_argument("--sink-env-name", default="gut_sink")
    parser.add_argument("--overwrite", action="store_true")

    args = parser.parse_args()

    tasks_keep_fp = Path(args.tasks_keep).resolve(strict=False)
    biom_fp = Path(args.biom).resolve(strict=False)
    out_dir = Path(args.out_dir).resolve(strict=False)
    task_root = out_dir / "tasks"

    out_dir.mkdir(parents=True, exist_ok=True)
    task_root.mkdir(parents=True, exist_ok=True)

    manifest_fp = out_dir / "st2_tasks_manifest.tsv"
    skipped_fp = out_dir / "st2_prepare_skipped_tasks.tsv"
    summary_fp = out_dir / "st2_prepare_summary.json"

    for p in [manifest_fp, skipped_fp, summary_fp]:
        if p.exists() and not args.overwrite:
            raise FileExistsError(f"Output exists: {p}. Use --overwrite.")

    print("Resolved paths:", flush=True)
    print(f"  tasks_keep: {tasks_keep_fp}", flush=True)
    print(f"  biom:       {biom_fp}", flush=True)
    print(f"  out_dir:    {out_dir}", flush=True)
    print(f"  manifest:   {manifest_fp}", flush=True)

    if not tasks_keep_fp.exists():
        raise FileNotFoundError(f"tasks_keep not found: {tasks_keep_fp}")
    if not biom_fp.exists():
        raise FileNotFoundError(f"BIOM not found: {biom_fp}")

    print("[1/4] Reading retained task table ...", flush=True)
    tasks = read_table(tasks_keep_fp)

    required_cols = ["dataset_label", "bioproject", "pair_id", "gut_sample_id", "oral_sample_ids", "n_oral_sources"]
    validate_columns(tasks, required_cols, "tasks_keep")

    tasks = tasks.loc[tasks["dataset_label"].astype(str).eq("filter"), :].copy()
    if tasks.empty:
        raise RuntimeError("No tasks with dataset_label == filter.")

    print(f"Retained ST2 gut-sink tasks: {tasks.shape[0]}", flush=True)

    print("[2/4] Loading global BIOM once ...", flush=True)
    table = load_biom(biom_fp)
    table_sample_ids = set(map(str, table.ids(axis="sample")))

    print(f"  BIOM samples: {len(table_sample_ids)}", flush=True)
    print(f"  BIOM features: {len(table.ids(axis='observation'))}", flush=True)

    print("[3/4] Preparing task input BIOM/map files ...", flush=True)

    manifest_rows = []
    skipped_rows = []
    task_id = 0

    for _, r in tasks.iterrows():
        dataset_label = str(r["dataset_label"])
        bp = str(r["bioproject"])
        pair_id = str(r["pair_id"])
        gut_sample_id = str(r["gut_sample_id"])
        oral_ids = parse_oral_sample_ids(r["oral_sample_ids"])
        oral_ids = [x for x in oral_ids if x in table_sample_ids]

        if gut_sample_id not in table_sample_ids:
            skipped_rows.append({
                "dataset_label": dataset_label,
                "bioproject": bp,
                "pair_id": pair_id,
                "gut_sample_id": gut_sample_id,
                "reason": "gut_sample_not_in_biom",
            })
            continue

        if len(oral_ids) == 0:
            skipped_rows.append({
                "dataset_label": dataset_label,
                "bioproject": bp,
                "pair_id": pair_id,
                "gut_sample_id": gut_sample_id,
                "reason": "no_oral_source_in_biom",
            })
            continue

        selected_samples = list(dict.fromkeys(oral_ids + [gut_sample_id]))
        sub_table = subset_biom_for_task(table, selected_samples)

        if len(sub_table.ids(axis="observation")) == 0:
            skipped_rows.append({
                "dataset_label": dataset_label,
                "bioproject": bp,
                "pair_id": pair_id,
                "gut_sample_id": gut_sample_id,
                "reason": "no_features_after_subsetting",
            })
            continue

        task_id += 1

        task_dir = (
            task_root
            / clean_for_path(dataset_label)
            / clean_for_path(bp)
            / f"task_{task_id:06d}"
        ).resolve(strict=False)
        task_dir.mkdir(parents=True, exist_ok=True)

        input_biom = (task_dir / "input.biom").resolve(strict=False)
        map_tsv = (task_dir / "map.tsv").resolve(strict=False)
        st2_out_dir = (task_dir / "st2_out").resolve(strict=False)

        if input_biom.exists() and not args.overwrite:
            raise FileExistsError(f"Output exists: {input_biom}. Use --overwrite.")
        if map_tsv.exists() and not args.overwrite:
            raise FileExistsError(f"Output exists: {map_tsv}. Use --overwrite.")

        write_biom(sub_table, input_biom)

        map_rows = []
        for oral_id in oral_ids:
            map_rows.append({
                "#SampleID": oral_id,
                "SourceSink": "source",
                "Env": args.source_env_name,
                "BioProject": bp,
                "pair_id": pair_id,
                "site": "Oral",
            })

        map_rows.append({
            "#SampleID": gut_sample_id,
            "SourceSink": "sink",
            "Env": args.sink_env_name,
            "BioProject": bp,
            "pair_id": pair_id,
            "site": "Gut",
        })

        pd.DataFrame(map_rows).to_csv(map_tsv, sep="\t", index=False)

        manifest_rows.append({
            "task_id": task_id,
            "dataset_label": dataset_label,
            "bioproject": bp,
            "pair_id": pair_id,
            "gut_sample_id": gut_sample_id,
            "oral_sample_ids": ";".join(oral_ids),
            "n_oral_sources": len(oral_ids),
            "input_biom": str(input_biom),
            "map_tsv": str(map_tsv),
            "st2_out_dir": str(st2_out_dir),
            "source_env_name": args.source_env_name,
            "sink_env_name": args.sink_env_name,
            "original_biom": str(biom_fp),
            "n_features_in_task": len(sub_table.ids(axis="observation")),
            "n_samples_in_task": len(sub_table.ids(axis="sample")),
        })

    print("[4/4] Writing manifest and summary ...", flush=True)

    manifest = pd.DataFrame(manifest_rows)
    skipped = pd.DataFrame(skipped_rows)

    if manifest.empty:
        raise RuntimeError("No ST2 tasks were prepared.")

    manifest.to_csv(manifest_fp, sep="\t", index=False)
    skipped.to_csv(skipped_fp, sep="\t", index=False)

    summary = {
        "tasks_keep": str(tasks_keep_fp),
        "biom": str(biom_fp),
        "out_dir": str(out_dir),
        "dataset_label": "filter",
        "definition": {
            "source": "paired Oral sample(s) from same BioProject and same pair_id",
            "sink": "one retained Gut sample",
            "source_env_name": args.source_env_name,
            "sink_env_name": args.sink_env_name,
        },
        "n_tasks_requested": int(tasks.shape[0]),
        "n_tasks_prepared": int(manifest.shape[0]),
        "n_tasks_skipped": int(skipped.shape[0]),
        "manifest": str(manifest_fp.resolve(strict=False)),
        "skipped_tasks": str(skipped_fp.resolve(strict=False)),
        "summary": str(summary_fp.resolve(strict=False)),
    }

    with open(summary_fp, "w") as f:
        json.dump(summary, f, indent=2)

    print("\nDone.", flush=True)
    print(f"Manifest: {manifest_fp}", flush=True)
    print(f"Skipped tasks: {skipped_fp}", flush=True)
    print(f"Summary: {summary_fp}", flush=True)
    print(f"Prepared tasks: {manifest.shape[0]}", flush=True)
    print(f"Skipped tasks: {skipped.shape[0]}", flush=True)

    print("\nPath check in manifest:", flush=True)
    for col in ["input_biom", "map_tsv", "st2_out_dir", "original_biom"]:
        n_abs = manifest[col].map(lambda x: Path(x).is_absolute()).sum()
        print(f"  {col}: absolute paths = {n_abs}/{manifest.shape[0]}", flush=True)


if __name__ == "__main__":
    main()
