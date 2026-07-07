#!/usr/bin/env python3
"""
Align exact-dereplicated oral ASVs to HOMD and summarize region-aware coverage.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
from typing import Dict, Iterable, List, Tuple

import pandas as pd


DEFAULT_OUTDIR = "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/region_aware_pipeline"
DEFAULT_HOMD = "/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03_full.fasta"
DEFAULT_VSEARCH = "/home/xuhuan/miniconda3/envs/vsearch_env/bin/vsearch"


def log(message: str) -> None:
    print(f"[homd] {message}", flush=True)


def run_vsearch(args: argparse.Namespace, derep_fasta: str, userout: str) -> None:
    if os.path.exists(userout) and os.path.getsize(userout) > 0 and not args.force:
        log(f"using existing HOMD alignment: {userout}")
        return
    cmd = [
        args.vsearch,
        "--usearch_global",
        derep_fasta,
        "--db",
        args.homd,
        "--id",
        str(args.min_id),
        "--maxaccepts",
        "1",
        "--maxrejects",
        "0",
        "--strand",
        "both",
        "--userout",
        userout,
        "--userfields",
        "query+target+id+alnlen+mism+gaps+qlo+qhi+tlo+thi+ql+tl+qcov+tcov",
        "--threads",
        str(args.threads),
    ]
    if args.quiet:
        cmd.append("--quiet")
    log("running HOMD alignment")
    log(" ".join(cmd))
    subprocess.run(cmd, check=True)


def split_values(value: str) -> List[str]:
    if pd.isna(value):
        return ["NA"]
    values = [v for v in str(value).split(";") if v and v != "nan"]
    return values or ["NA"]


def load_homd_hits(path: str) -> pd.DataFrame:
    columns = [
        "derep_id",
        "homd_target",
        "identity",
        "alnlen",
        "mismatches",
        "gaps",
        "qlo",
        "qhi",
        "tlo",
        "thi",
        "query_length",
        "target_length",
        "query_cov",
        "target_cov",
    ]
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame(columns=columns)
    hits = pd.read_csv(path, sep="\t", names=columns)
    # vsearch may emit multiple equal best hits for a query; coverage needs one row per derep sequence.
    hits = hits.sort_values(["derep_id", "identity", "query_cov", "alnlen"], ascending=[True, False, False, False])
    return hits.drop_duplicates(subset=["derep_id"], keep="first")


def add_capture_flags(df: pd.DataFrame, thresholds: Iterable[float], min_qcov: float) -> pd.DataFrame:
    df = df.copy()
    df["has_homd_hit"] = df["identity"].notna()
    df["high_qcov"] = df["query_cov"].fillna(0) >= min_qcov
    for threshold in thresholds:
        label = f"captured_{str(threshold).replace('.', 'p')}"
        df[label] = (df["identity"].fillna(0) >= threshold) & df["high_qcov"]
    return df


def weighted_summary(df: pd.DataFrame, group_col: str, capture_cols: List[str]) -> pd.DataFrame:
    rows = []
    expanded = []
    for _, row in df.iterrows():
        for value in split_values(row[group_col]):
            item = row.to_dict()
            item[group_col] = value
            expanded.append(item)
    if not expanded:
        return pd.DataFrame()
    ex = pd.DataFrame(expanded)
    for value, sub in ex.groupby(group_col, dropna=False, sort=True):
        base = {
            "grouping": group_col,
            "level": value,
            "derep_sequences": int(len(sub)),
            "source_asvs": int(sub["n_asvs"].sum()),
            "bioproject_sequence_links": int(sub["n_bioprojects"].sum()),
            "total_reads": float(sub["total_reads"].sum()),
            "prevalence_samples": int(sub["prevalence_samples"].sum()),
        }
        for col in capture_cols:
            captured = sub[sub[col]]
            suffix = col.replace("captured_", "")
            base[f"{suffix}_derep_fraction"] = float(len(captured) / len(sub)) if len(sub) else 0.0
            base[f"{suffix}_asv_fraction"] = float(captured["n_asvs"].sum() / sub["n_asvs"].sum()) if sub["n_asvs"].sum() else 0.0
            base[f"{suffix}_read_fraction"] = (
                float(captured["total_reads"].sum() / sub["total_reads"].sum()) if sub["total_reads"].sum() else 0.0
            )
            base[f"{suffix}_prevalence_fraction"] = (
                float(captured["prevalence_samples"].sum() / sub["prevalence_samples"].sum())
                if sub["prevalence_samples"].sum()
                else 0.0
            )
        rows.append(base)
    return pd.DataFrame(rows)


def reproducible_uncovered(df: pd.DataFrame, capture_col: str, args: argparse.Namespace) -> pd.DataFrame:
    candidates = df[
        (~df[capture_col])
        & (df["n_bioprojects"] >= args.min_projects)
        & (df["total_reads"] >= args.min_reads)
        & (df["prevalence_samples"] >= args.min_prevalence)
    ].copy()
    return candidates.sort_values(
        ["n_bioprojects", "total_reads", "prevalence_samples"], ascending=[False, False, False]
    )


def summarize(args: argparse.Namespace) -> None:
    derep_table = os.path.join(args.outdir, "exact_derep_summary.tsv")
    derep_fasta = os.path.join(args.outdir, "exact_derep_representatives.fa")
    if not os.path.exists(derep_table) or not os.path.exists(derep_fasta):
        raise FileNotFoundError("Run 01_collect_asvs_and_abundance.py first.")

    userout = os.path.join(args.outdir, "homd_vsearch_best_hits.tsv")
    run_vsearch(args, derep_fasta, userout)

    derep = pd.read_csv(derep_table, sep="\t")
    hits = load_homd_hits(userout)
    merged = derep.merge(hits, on="derep_id", how="left")
    merged = add_capture_flags(merged, args.identity_thresholds, args.min_qcov)

    merged_path = os.path.join(args.outdir, "homd_coverage_by_derep_sequence.tsv")
    merged.to_csv(merged_path, sep="\t", index=False)

    capture_cols = [f"captured_{str(t).replace('.', 'p')}" for t in args.identity_thresholds]
    grouping_cols = [
        "sequencing_regions",
        "platforms",
        "continents",
        "countries",
        "sites",
        "disease_studies",
        "statuses",
    ]
    summary_parts = [weighted_summary(merged, col, capture_cols) for col in grouping_cols]
    summary = pd.concat([part for part in summary_parts if not part.empty], ignore_index=True)
    summary_path = os.path.join(args.outdir, "homd_coverage_grouped_summary.tsv")
    summary.to_csv(summary_path, sep="\t", index=False)

    primary_capture = f"captured_{str(args.primary_identity).replace('.', 'p')}"
    if primary_capture not in merged.columns:
        raise ValueError(f"primary identity {args.primary_identity} is not in --identity-thresholds")
    candidates = reproducible_uncovered(merged, primary_capture, args)
    candidates_path = os.path.join(args.outdir, "reproducible_abundant_not_high_identity_homd.tsv")
    candidates.to_csv(candidates_path, sep="\t", index=False)

    stats = {
        "derep_sequences": int(len(merged)),
        "derep_with_any_homd_hit": int(merged["has_homd_hit"].sum()),
        "min_query_coverage_for_capture": args.min_qcov,
        "identity_thresholds": args.identity_thresholds,
        "primary_identity": args.primary_identity,
        "reproducible_candidate_filters": {
            "min_projects": args.min_projects,
            "min_reads": args.min_reads,
            "min_prevalence": args.min_prevalence,
        },
        "reproducible_candidates": int(len(candidates)),
        "outputs": {
            "best_hits": userout,
            "coverage_by_derep_sequence": merged_path,
            "grouped_summary": summary_path,
            "reproducible_candidates": candidates_path,
        },
    }
    stats_path = os.path.join(args.outdir, "homd_summary_stats.json")
    with open(stats_path, "w") as handle:
        json.dump(stats, handle, indent=2)

    log(f"complete: {len(merged):,} exact-derep sequences summarized")
    log(f"summary: {summary_path}")
    log(f"candidate table: {candidates_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--homd", default=DEFAULT_HOMD)
    parser.add_argument("--vsearch", default=DEFAULT_VSEARCH)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--min-id", type=float, default=0.8, help="Loose reporting threshold for vsearch.")
    parser.add_argument("--identity-thresholds", type=float, nargs="+", default=[100.0, 99.5, 99.0, 97.0])
    parser.add_argument("--primary-identity", type=float, default=99.0)
    parser.add_argument("--min-qcov", type=float, default=80.0)
    parser.add_argument("--min-projects", type=int, default=3)
    parser.add_argument("--min-reads", type=float, default=1000.0)
    parser.add_argument("--min-prevalence", type=int, default=5)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--quiet", action="store_true", default=True)
    return parser.parse_args()


if __name__ == "__main__":
    summarize(parse_args())
