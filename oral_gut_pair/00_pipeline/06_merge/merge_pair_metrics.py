#!/usr/bin/env python3
"""Merge harmonized single-pair and time-series 10-metric tables."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


METRICS = [
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


def load(path: Path, analysis_set: str) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t")
    required = {"bioproject", "pair_id", *METRICS}
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(f"{path} missing columns: {missing}")
    if table.duplicated(["bioproject", "pair_id"]).any():
        raise ValueError(f"{path} has duplicate BioProject/pair_id keys")
    table = table.copy()
    table["analysis_set"] = analysis_set
    table["is_time_series"] = analysis_set == "time_series"
    return table


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--single", required=True)
    parser.add_argument("--time-series", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--availability", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    single_path = Path(args.single).resolve()
    time_path = Path(args.time_series).resolve()
    output = Path(args.output).resolve()
    availability_path = Path(args.availability).resolve()
    summary_path = Path(args.summary).resolve()
    for path in (output, availability_path, summary_path):
        path.parent.mkdir(parents=True, exist_ok=True)

    single = load(single_path, "single_pair")
    time = load(time_path, "time_series")
    merged = pd.concat([single, time], ignore_index=True, sort=False)
    if merged.duplicated(["bioproject", "pair_id"]).any():
        raise ValueError("Merged table has duplicate BioProject/pair_id keys")
    merged["n_available_metrics"] = merged[METRICS].notna().sum(axis=1)
    merged["has_all_10_metrics"] = merged["n_available_metrics"].eq(10)
    merged = merged.sort_values(["analysis_set", "bioproject", "pair_id"])

    rows = []
    for analysis_set, subset in merged.groupby("analysis_set"):
        for metric in METRICS:
            rows.append(
                {
                    "analysis_set": analysis_set,
                    "metric": metric,
                    "n_total_pairs": int(subset.shape[0]),
                    "n_nonmissing": int(subset[metric].notna().sum()),
                    "fraction_nonmissing": float(subset[metric].notna().mean()),
                }
            )
    availability = pd.DataFrame(rows)

    merged.to_csv(output, sep="\t", index=False)
    availability.to_csv(availability_path, sep="\t", index=False)
    summary = {
        "inputs": {
            "single_pair": str(single_path),
            "time_series": str(time_path),
        },
        "output": str(output),
        "n_rows": int(merged.shape[0]),
        "n_pairs": int(merged[["bioproject", "pair_id"]].drop_duplicates().shape[0]),
        "n_bioprojects": int(merged["bioproject"].nunique()),
        "n_pairs_by_analysis_set": merged["analysis_set"].value_counts().to_dict(),
        "n_pairs_all_10_metrics": int(merged["has_all_10_metrics"].sum()),
        "metrics": METRICS,
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
