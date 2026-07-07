#!/usr/bin/env python3
"""Recompute HOMD capture using stricter query coverage thresholds without rerunning vsearch."""

from __future__ import annotations

import argparse
import csv
import json
import os
from collections import defaultdict

import pandas as pd

DEFAULT_OUTDIR = "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/region_aware_pipeline"
IDENTITY_THRESHOLDS = [100.0, 99.5, 99.0, 97.0]
QCOV_THRESHOLDS = [100.0, 99.0, 95.0, 80.0]
GROUPING_COLS = [
    "sequencing_regions",
    "platforms",
    "continents",
    "countries",
    "sites",
    "disease_studies",
    "statuses",
]
NUMERIC_COLS = ["n_asvs", "n_bioprojects", "total_reads", "prevalence_samples"]


def log(msg):
    print(f"[qcov] {msg}", flush=True)


def split_values(value):
    if pd.isna(value):
        return ["NA"]
    vals = [v for v in str(value).split(";") if v and v != "nan"]
    return vals or ["NA"]


def capture_col(identity, qcov):
    return f"captured_id{str(identity).replace('.', 'p')}_qcov{str(qcov).replace('.', 'p')}"


def empty_stats(capture_cols):
    stats = {
        "derep_sequences": 0,
        "source_asvs": 0.0,
        "bioproject_sequence_links": 0.0,
        "total_reads": 0.0,
        "prevalence_samples": 0.0,
    }
    for col in capture_cols:
        stats[f"{col}_derep"] = 0
        stats[f"{col}_asvs"] = 0.0
        stats[f"{col}_reads"] = 0.0
        stats[f"{col}_prevalence"] = 0.0
    return stats


def update_stats(stats, row, capture_cols):
    stats["derep_sequences"] += 1
    stats["source_asvs"] += float(row.n_asvs)
    stats["bioproject_sequence_links"] += float(row.n_bioprojects)
    stats["total_reads"] += float(row.total_reads)
    stats["prevalence_samples"] += float(row.prevalence_samples)
    for col in capture_cols:
        if bool(getattr(row, col)):
            stats[f"{col}_derep"] += 1
            stats[f"{col}_asvs"] += float(row.n_asvs)
            stats[f"{col}_reads"] += float(row.total_reads)
            stats[f"{col}_prevalence"] += float(row.prevalence_samples)


def write_grouped(group_stats, capture_cols, out_path):
    fields = [
        "grouping", "level", "derep_sequences", "source_asvs", "bioproject_sequence_links",
        "total_reads", "prevalence_samples",
    ]
    for col in capture_cols:
        fields.extend([
            f"{col}_derep_fraction",
            f"{col}_asv_fraction",
            f"{col}_read_fraction",
            f"{col}_prevalence_fraction",
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
            for col in capture_cols:
                row[f"{col}_derep_fraction"] = stats[f"{col}_derep"] / stats["derep_sequences"] if stats["derep_sequences"] else 0
                row[f"{col}_asv_fraction"] = stats[f"{col}_asvs"] / stats["source_asvs"] if stats["source_asvs"] else 0
                row[f"{col}_read_fraction"] = stats[f"{col}_reads"] / stats["total_reads"] if stats["total_reads"] else 0
                row[f"{col}_prevalence_fraction"] = stats[f"{col}_prevalence"] / stats["prevalence_samples"] if stats["prevalence_samples"] else 0
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--chunksize", type=int, default=100000)
    parser.add_argument("--primary-identity", type=float, default=99.0)
    parser.add_argument("--primary-qcov", type=float, default=99.0)
    parser.add_argument("--min-projects", type=int, default=3)
    parser.add_argument("--min-reads", type=float, default=1000.0)
    parser.add_argument("--min-prevalence", type=int, default=5)
    args = parser.parse_args()

    in_path = os.path.join(args.outdir, "homd_coverage_by_derep_sequence.tsv")
    out_prefix = os.path.join(args.outdir, "qcov_sensitivity")
    perseq_path = out_prefix + "_by_derep_sequence.tsv"
    grouped_path = out_prefix + "_grouped_summary.tsv"
    candidates_path = out_prefix + "_reproducible_abundant_not_high_identity_homd.tsv"
    stats_path = out_prefix + "_stats.json"

    capture_cols = [capture_col(i, q) for i in IDENTITY_THRESHOLDS for q in QCOV_THRESHOLDS]
    primary_col = capture_col(args.primary_identity, args.primary_qcov)

    usecols = [
        "derep_id", "sequence", "length", "n_asvs", "n_bioprojects", "total_reads", "prevalence_samples",
        "representative_asv_id", "bioprojects", "identity", "query_cov", "has_homd_hit",
    ] + GROUPING_COLS

    group_stats = defaultdict(lambda: empty_stats(capture_cols))
    total_rows = 0
    any_hit = 0
    candidate_rows = 0
    first_perseq = True
    first_candidates = True
    capture_totals = {col: {"derep": 0, "asvs": 0.0, "reads": 0.0, "prevalence": 0.0} for col in capture_cols}
    denominators = {"asvs": 0.0, "reads": 0.0, "prevalence": 0.0}

    for chunk in pd.read_csv(in_path, sep="\t", usecols=usecols, chunksize=args.chunksize):
        for col in NUMERIC_COLS + ["identity", "query_cov"]:
            chunk[col] = pd.to_numeric(chunk[col], errors="coerce").fillna(0)
        chunk["has_homd_hit"] = chunk["has_homd_hit"].astype(bool)

        for identity in IDENTITY_THRESHOLDS:
            for qcov in QCOV_THRESHOLDS:
                col = capture_col(identity, qcov)
                chunk[col] = (chunk["identity"] >= identity) & (chunk["query_cov"] >= qcov)

        total_rows += len(chunk)
        any_hit += int(chunk["has_homd_hit"].sum())
        denominators["asvs"] += float(chunk["n_asvs"].sum())
        denominators["reads"] += float(chunk["total_reads"].sum())
        denominators["prevalence"] += float(chunk["prevalence_samples"].sum())

        for col in capture_cols:
            captured = chunk[col]
            capture_totals[col]["derep"] += int(captured.sum())
            capture_totals[col]["asvs"] += float(chunk.loc[captured, "n_asvs"].sum())
            capture_totals[col]["reads"] += float(chunk.loc[captured, "total_reads"].sum())
            capture_totals[col]["prevalence"] += float(chunk.loc[captured, "prevalence_samples"].sum())

        perseq_cols = usecols + capture_cols
        chunk[perseq_cols].to_csv(
            perseq_path, sep="\t", index=False, mode="w" if first_perseq else "a", header=first_perseq
        )
        first_perseq = False

        candidate_mask = (
            ~chunk[primary_col]
            & (chunk["n_bioprojects"] >= args.min_projects)
            & (chunk["total_reads"] >= args.min_reads)
            & (chunk["prevalence_samples"] >= args.min_prevalence)
        )
        candidates = chunk.loc[candidate_mask, perseq_cols].sort_values(
            ["n_bioprojects", "total_reads", "prevalence_samples"], ascending=[False, False, False]
        )
        candidate_rows += len(candidates)
        candidates.to_csv(
            candidates_path, sep="\t", index=False, mode="w" if first_candidates else "a", header=first_candidates
        )
        first_candidates = False

        for row in chunk.itertuples(index=False):
            for grouping in GROUPING_COLS:
                for level in split_values(getattr(row, grouping)):
                    update_stats(group_stats[(grouping, level)], row, capture_cols)
        log(f"processed {total_rows:,} derep rows")

    if first_candidates:
        pd.DataFrame(columns=usecols + capture_cols).to_csv(candidates_path, sep="\t", index=False)

    write_grouped(group_stats, capture_cols, grouped_path)
    stats = {
        "derep_sequences": total_rows,
        "derep_with_any_homd_hit": any_hit,
        "identity_thresholds": IDENTITY_THRESHOLDS,
        "query_coverage_thresholds": QCOV_THRESHOLDS,
        "primary_capture": primary_col,
        "denominators": denominators,
        "capture_totals": capture_totals,
        "candidate_filters": {
            "min_projects": args.min_projects,
            "min_reads": args.min_reads,
            "min_prevalence": args.min_prevalence,
        },
        "reproducible_candidates": candidate_rows,
        "outputs": {
            "per_sequence": perseq_path,
            "grouped_summary": grouped_path,
            "reproducible_candidates": candidates_path,
        },
    }
    with open(stats_path, "w") as handle:
        json.dump(stats, handle, indent=2)
    log(f"complete: {total_rows:,} derep rows; candidates={candidate_rows:,}; primary={primary_col}")


if __name__ == "__main__":
    main()
