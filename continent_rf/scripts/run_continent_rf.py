#!/usr/bin/env python3
"""Train oral/gut continent RF; choose feature level, task and CV grouping."""

import argparse
from dataclasses import asdict
from datetime import datetime
import json
import os
from pathlib import Path

import pandas as pd

from rf_cv import CVConfig, output_columns, run_cv
from rf_data import ORAL_SITES, REGIONS, load_dataset, load_genus_map

MODULE_DIR = Path(__file__).resolve().parents[1]
ROOT = Path(os.environ.get("GLOBAL_ORAL_ROOT") or MODULE_DIR.parent)


def positive_int(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return number


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--level", choices=["genus", "feature"], default="genus")
    parser.add_argument("--tasks", nargs="+", choices=["ovo", "ova"], default=["ovo"])
    parser.add_argument("--datasets", nargs="+", choices=["oral", "gut"], default=["oral", "gut"])
    parser.add_argument("--split-mode", choices=["project", "unit"],
                        default=os.environ.get("GLOBAL_ORAL_SPLIT_MODE", "project"))
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--taxonomy", type=Path,
                        default=Path(os.environ.get("GLOBAL_ORAL_TAXONOMY") or
                                     (ROOT / "taxonomy.tsv" if (ROOT / "taxonomy.tsv").exists() else
                                      "/beegfs/huangxiaochang/global_oral/taxonomy.tsv")))
    for name, relative in {
        "oral-metadata": "data/metadata-and-biom-v2-qc/metadata_oral_qc_dedup.tsv",
        "oral-biom": "data/metadata-and-biom-v2-qc/biom_oral_qc_dedup.biom",
        "gut-metadata": "data/metadata_gut_v2.txt",
        "gut-biom": "data/table_filter_reads5000.biom",
    }.items():
        parser.add_argument(f"--{name}", type=Path, default=ROOT / relative)
    parser.add_argument("--repeats", type=positive_int, default=5)
    parser.add_argument("--folds-max", type=positive_int, default=5)
    parser.add_argument("--trees", type=positive_int, default=200)
    parser.add_argument("--jobs", type=positive_int, default=8)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--min-study-samples", type=positive_int, default=20)
    args = parser.parse_args(argv)
    if args.folds_max < 2:
        parser.error("--folds-max must be at least 2")
    if args.split_mode not in {"project", "unit"}:
        parser.error("GLOBAL_ORAL_SPLIT_MODE must be project or unit")
    args.datasets = list(dict.fromkeys(args.datasets))
    args.tasks = list(dict.fromkeys(args.tasks))
    if args.output_dir is None:
        suffix = "_unitsplit" if args.split_mode == "unit" else ""
        args.output_dir = MODULE_DIR / f"result_{args.level}_{args.repeats}rep{suffix}"
    required = [getattr(args, f"{dataset}_{kind}") for dataset in args.datasets
                for kind in ["metadata", "biom"]]
    if args.level == "genus":
        required.append(args.taxonomy)
    for path in required:
        if not path.is_file():
            parser.error(f"Input file does not exist: {path}")
    return args


def main(argv=None):
    args = parse_args(argv)
    config = CVConfig(**{name: getattr(args, name) for name in
                         ["folds_max", "repeats", "seed", "trees", "jobs"]})
    genus_map = load_genus_map(args.taxonomy) if args.level == "genus" else None
    provenance = {
        **{key: str(value) if isinstance(value, Path) else value
           for key, value in vars(args).items()},
        **asdict(config), "generated": datetime.now().isoformat(),
        "script": Path(__file__).name, "regions": REGIONS, "oral_sites": ORAL_SITES,
        "split_unit": "project" if args.split_mode == "project" else "region||project",
        "normalization": "unit-sum per sample after genus aggregation" if genus_map is not None
                         else "unit-sum per sample",
    }
    for dataset in args.datasets:
        meta, matrix = load_dataset(
            dataset, getattr(args, f"{dataset}_metadata"), getattr(args, f"{dataset}_biom"),
            split_mode=args.split_mode, min_study_samples=args.min_study_samples,
            genus_map=genus_map,
        )
        print(f"{dataset} matrix: {matrix.shape}; split={args.split_mode}", flush=True)
        out = args.output_dir / dataset
        out.mkdir(parents=True, exist_ok=True)
        for task in args.tasks:
            columns = output_columns(task)
            prediction_path = out / f"region_{task}_predictions.tsv"
            pd.DataFrame(columns=columns["predictions"]).to_csv(prediction_path, sep="\t", index=False)
            folds, summaries = run_cv(meta, matrix, dataset, task, prediction_path, config)
            for name, rows in [("fold_metrics", folds), ("summary", summaries)]:
                pd.DataFrame(rows, columns=columns[name]).to_csv(
                    out / f"region_{task}_{name}.tsv", sep="\t", index=False)
        (out / "run_config.json").write_text(json.dumps(
            {**provenance, "datasets": [dataset]}, indent=2, sort_keys=True) + "\n")
    (args.output_dir / "run_config.json").write_text(json.dumps(
        provenance, indent=2, sort_keys=True) + "\n")
    print(f"Outputs written to {args.output_dir}")


if __name__ == "__main__":
    main()
