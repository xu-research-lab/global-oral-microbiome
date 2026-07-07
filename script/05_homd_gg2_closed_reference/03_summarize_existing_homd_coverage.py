#!/usr/bin/env python3
"""Chunked summaries for completed HOMD per-sequence coverage table."""

from __future__ import annotations

import argparse
import csv
import json
import os
from collections import defaultdict

import pandas as pd

DEFAULT_OUTDIR = "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/region_aware_pipeline"

GROUPING_COLS = [
    "sequencing_regions",
    "platforms",
    "continents",
    "countries",
    "sites",
    "disease_studies",
    "statuses",
]
CAPTURE_COLS = ["captured_100p0", "captured_99p5", "captured_99p0", "captured_97p0"]
NUMERIC_COLS = ["n_asvs", "n_bioprojects", "total_reads", "prevalence_samples"]


def split_values(value):
    if pd.isna(value):
        return ["NA"]
    vals = [v for v in str(value).split(";") if v and v != "nan"]
    return vals or ["NA"]


def empty_group_stats():
    stats = {
        "derep_sequences": 0,
        "source_asvs": 0.0,
        "bioproject_sequence_links": 0.0,
        "total_reads": 0.0,
        "prevalence_samples": 0.0,
    }
    for col in CAPTURE_COLS:
        stats[f"{col}_derep"] = 0
        stats[f"{col}_asvs"] = 0.0
        stats[f"{col}_reads"] = 0.0
        stats[f"{col}_prevalence"] = 0.0
    return stats


def update_group_stats(stats, row):
    stats["derep_sequences"] += 1
    stats["source_asvs"] += float(row.n_asvs)
    stats["bioproject_sequence_links"] += float(row.n_bioprojects)
    stats["total_reads"] += float(row.total_reads)
    stats["prevalence_samples"] += float(row.prevalence_samples)
    for col in CAPTURE_COLS:
        if bool(getattr(row, col)):
            stats[f"{col}_derep"] += 1
            stats[f"{col}_asvs"] += float(row.n_asvs)
            stats[f"{col}_reads"] += float(row.total_reads)
            stats[f"{col}_prevalence"] += float(row.prevalence_samples)


def write_grouped_summary(group_stats, out_path):
    fields = [
        "grouping",
        "level",
        "derep_sequences",
        "source_asvs",
        "bioproject_sequence_links",
        "total_reads",
        "prevalence_samples",
    ]
    for col in CAPTURE_COLS:
        suffix = col.replace("captured_", "")
        fields.extend([
            f"{suffix}_derep_fraction",
            f"{suffix}_asv_fraction",
            f"{suffix}_read_fraction",
            f"{suffix}_prevalence_fraction",
        ])
    with open(out_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for (grouping, level), stats in sorted(group_stats.items()):
            row = {
                "grouping": grouping,
                "level": level,
                "derep_sequences": int(stats["derep_sequences"]),
                "source_asvs": int(stats["source_asvs"]),
                "bioproject_sequence_links": int(stats["bioproject_sequence_links"]),
                "total_reads": stats["total_reads"],
                "prevalence_samples": int(stats["prevalence_samples"]),
            }
            for col in CAPTURE_COLS:
                suffix = col.replace("captured_", "")
                row[f"{suffix}_derep_fraction"] = stats[f"{col}_derep"] / stats["derep_sequences"] if stats["derep_sequences"] else 0
                row[f"{suffix}_asv_fraction"] = stats[f"{col}_asvs"] / stats["source_asvs"] if stats["source_asvs"] else 0
                row[f"{suffix}_read_fraction"] = stats[f"{col}_reads"] / stats["total_reads"] if stats["total_reads"] else 0
                row[f"{suffix}_prevalence_fraction"] = stats[f"{col}_prevalence"] / stats["prevalence_samples"] if stats["prevalence_samples"] else 0
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--chunksize", type=int, default=100000)
    parser.add_argument("--primary-capture", default="captured_99p0")
    parser.add_argument("--min-projects", type=int, default=3)
    parser.add_argument("--min-reads", type=float, default=1000.0)
    parser.add_argument("--min-prevalence", type=int, default=5)
    args = parser.parse_args()

    in_path = os.path.join(args.outdir, "homd_coverage_by_derep_sequence.tsv")
    grouped_path = os.path.join(args.outdir, "homd_coverage_grouped_summary.tsv")
    candidates_path = os.path.join(args.outdir, "reproducible_abundant_not_high_identity_homd.tsv")
    stats_path = os.path.join(args.outdir, "homd_summary_stats.json")

    usecols = [
        "derep_id", "sequence", "length", "n_asvs", "n_bioprojects", "total_reads", "prevalence_samples",
        "representative_asv_id", "bioprojects", "identity", "query_cov", "has_homd_hit",
    ] + GROUPING_COLS + CAPTURE_COLS

    group_stats = defaultdict(empty_group_stats)
    total_rows = 0
    any_hit = 0
    candidate_rows = 0
    first_candidate_chunk = True
    capture_totals = {col: {"derep": 0, "asvs": 0.0, "reads": 0.0, "prevalence": 0.0} for col in CAPTURE_COLS}
    denominators = {"asvs": 0.0, "reads": 0.0, "prevalence": 0.0}

    for chunk in pd.read_csv(in_path, sep="\t", usecols=usecols, chunksize=args.chunksize):
        for col in NUMERIC_COLS:
            chunk[col] = pd.to_numeric(chunk[col], errors="coerce").fillna(0)
        for col in CAPTURE_COLS:
            chunk[col] = chunk[col].astype(bool)
        chunk["has_homd_hit"] = chunk["has_homd_hit"].astype(bool)

        total_rows += len(chunk)
        any_hit += int(chunk["has_homd_hit"].sum())
        denominators["asvs"] += float(chunk["n_asvs"].sum())
        denominators["reads"] += float(chunk["total_reads"].sum())
        denominators["prevalence"] += float(chunk["prevalence_samples"].sum())
        for col in CAPTURE_COLS:
            captured = chunk[col]
            capture_totals[col]["derep"] += int(captured.sum())
            capture_totals[col]["asvs"] += float(chunk.loc[captured, "n_asvs"].sum())
            capture_totals[col]["reads"] += float(chunk.loc[captured, "total_reads"].sum())
            capture_totals[col]["prevalence"] += float(chunk.loc[captured, "prevalence_samples"].sum())

        candidate_mask = (
            ~chunk[args.primary_capture]
            & (chunk["n_bioprojects"] >= args.min_projects)
            & (chunk["total_reads"] >= args.min_reads)
            & (chunk["prevalence_samples"] >= args.min_prevalence)
        )
        candidates = chunk.loc[candidate_mask].sort_values(
            ["n_bioprojects", "total_reads", "prevalence_samples"], ascending=[False, False, False]
        )
        candidate_rows += len(candidates)
        candidates.to_csv(
            candidates_path,
            sep="\t",
            index=False,
            mode="w" if first_candidate_chunk else "a",
            header=first_candidate_chunk,
        )
        first_candidate_chunk = False

        for row in chunk.itertuples(index=False):
            for grouping in GROUPING_COLS:
                for level in split_values(getattr(row, grouping)):
                    update_group_stats(group_stats[(grouping, level)], row)

        print(f"[summary] processed {total_rows:,} derep rows", flush=True)

    if first_candidate_chunk:
        pd.DataFrame(columns=usecols).to_csv(candidates_path, sep="\t", index=False)

    write_grouped_summary(group_stats, grouped_path)

    stats = {
        "derep_sequences": total_rows,
        "derep_with_any_homd_hit": any_hit,
        "capture_totals": capture_totals,
        "denominators": denominators,
        "primary_capture": args.primary_capture,
        "candidate_filters": {
            "min_projects": args.min_projects,
            "min_reads": args.min_reads,
            "min_prevalence": args.min_prevalence,
        },
        "reproducible_candidates": candidate_rows,
        "outputs": {
            "coverage_by_derep_sequence": in_path,
            "grouped_summary": grouped_path,
            "reproducible_candidates": candidates_path,
        },
    }
    with open(stats_path, "w") as handle:
        json.dump(stats, handle, indent=2)
    print(f"[summary] complete: {total_rows:,} derep rows; candidates={candidate_rows:,}", flush=True)


if __name__ == "__main__":
    main()
