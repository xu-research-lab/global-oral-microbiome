#!/usr/bin/env python3
"""Train HC age random forests with project or sample cross-validation."""

import argparse
from dataclasses import asdict
import json
from pathlib import Path

from scipy.stats import spearmanr

from age_common import (GENUS_TABLE_PATH, METADATA_PATH, ORAL_SITES, RESULT,
                        SCHEMES, load_metadata, prediction_files)
from age_model import ModelConfig, fit_fold, load_genus_table, split_pool


def positive_int(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return number


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scheme", choices=SCHEMES, default="project")
    parser.add_argument("--sites", choices=ORAL_SITES, nargs="+", default=ORAL_SITES)
    parser.add_argument("--metadata", type=Path, default=METADATA_PATH)
    parser.add_argument("--genus-table", type=Path, default=GENUS_TABLE_PATH)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--trees", type=positive_int, default=500)
    parser.add_argument("--folds", type=positive_int, default=5)
    parser.add_argument("--jobs", type=positive_int, default=16)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--min-study-samples", type=positive_int, default=11)
    args = parser.parse_args(argv)
    if args.folds < 2:
        parser.error("--folds must be at least 2")
    args.sites = list(dict.fromkeys(args.sites))
    config = ModelConfig(**{key: getattr(args, key) for key in ["trees", "folds", "seed", "jobs"]})
    directory, token = SCHEMES[args.scheme]
    output = args.output_dir or RESULT / directory
    metadata = load_metadata(args.metadata)
    metadata = metadata[metadata.site_cleaned.isin(args.sites)
                        & metadata.age_cleaned.between(1, 80)].copy()
    features, sample_ids, matrix = load_genus_table(set(metadata.sample_id), args.genus_table)
    metadata = metadata[metadata.sample_id.isin(sample_ids)]
    sid_to_age = metadata.set_index("sample_id")["age_cleaned"]
    print(f"Genus table: {len(features)} features, {len(sample_ids)} samples", flush=True)
    output.mkdir(parents=True, exist_ok=True)
    completed = {}
    for site in args.sites:
        md_site = metadata[metadata.site_cleaned == site]
        counts = md_site.groupby("bioproject_id").sample_id.count()
        pool = md_site[md_site.bioproject_id.isin(counts[counts >= args.min_study_samples].index)]
        predictions = []
        for train, test in split_pool(pool, args.scheme, config):
            frame = fit_fold(matrix, sample_ids, train, test, sid_to_age, config)
            if frame is None:
                continue
            predictions.append(frame)
            print(f"[{site}] fold {len(predictions)}: n_test={len(frame)} "
                  f"spearman={spearmanr(frame.true_age, frame.pred_age)[0]:+.3f}", flush=True)
        # Finish fitting the site before replacing its prediction tables.
        for path in prediction_files(output, token, site):
            path.unlink()
        for fold, frame in enumerate(predictions, 1):
            frame.to_csv(output / f"age_predict_test_{fold}_{site}_{token}_range_1_80.csv", index=False)
        completed[site] = len(predictions)
        if not predictions:
            print(f"SKIP {site}: insufficient groups or valid HC training folds", flush=True)
    provenance = {**{key: str(value) if isinstance(value, Path) else value
                      for key, value in vars(args).items()},
                  **asdict(config), "output_dir": str(output), "completed_folds": completed,
                  "age_range": [1, 80], "training_status": "HC", "transform": "log10(relative_abundance + 1e-6)"}
    (output / "run_config.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(f"Outputs written to {output}")


if __name__ == "__main__":
    main()
