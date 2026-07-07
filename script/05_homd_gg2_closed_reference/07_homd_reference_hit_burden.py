#!/usr/bin/env python3
"""Summarize observed ASV burden on each HOMD 16S reference sequence."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
from collections import defaultdict

import pandas as pd

DEFAULT_OUTDIR = "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/region_aware_pipeline"
DEFAULT_HOMD = "/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03_full.fasta"
IDENTITY_THRESHOLDS = [100.0, 99.5, 99.0, 97.0, 95.0, 90.0, 85.0, 80.0]
QCOV = 99.0


def log(msg):
    print(f"[homd-target] {msg}", flush=True)


def threshold_label(threshold):
    return str(threshold).replace(".", "p")


def parse_homd_fasta(path):
    rows = []
    current_id = None
    current_desc = None
    seq_len = 0
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    rows.append(parse_header(current_id, current_desc, seq_len))
                header = line[1:]
                parts = header.split(maxsplit=1)
                current_id = parts[0]
                current_desc = parts[1] if len(parts) > 1 else ""
                seq_len = 0
            else:
                seq_len += len(line.strip())
        if current_id is not None:
            rows.append(parse_header(current_id, current_desc, seq_len))
    return pd.DataFrame(rows)


def parse_header(ref_id, desc, seq_len):
    parts = desc.split(maxsplit=1)
    hmt = parts[0] if parts and parts[0].startswith("HMT-") else "NA"
    species = parts[1] if len(parts) > 1 else desc
    return {
        "homd_target": ref_id,
        "homd_hmt": hmt,
        "homd_species": species if species else "NA",
        "homd_header_description": desc,
        "homd_ref_length": seq_len,
    }


def load_hits(path):
    cols = [
        "derep_id", "homd_target", "identity", "alnlen", "mismatches", "gaps", "qlo", "qhi", "tlo", "thi",
        "query_length", "target_length", "query_cov", "target_cov",
    ]
    hits = pd.read_csv(path, sep="\t", names=cols)
    # One target per query; ties are sorted deterministically and collapsed.
    hits = hits.sort_values(["derep_id", "identity", "query_cov", "alnlen", "homd_target"], ascending=[True, False, False, False, True])
    return hits.drop_duplicates("derep_id", keep="first")


def split_values(value):
    if pd.isna(value):
        return []
    return [v for v in str(value).split(";") if v and v != "nan"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--homd", default=DEFAULT_HOMD)
    args = parser.parse_args()

    derep_path = os.path.join(args.outdir, "exact_derep_summary.tsv")
    hit_path = os.path.join(args.outdir, "homd_vsearch_best_hits.tsv")
    out_ref = os.path.join(args.outdir, "homd_reference_hit_burden_qcov99.tsv")
    out_hmt = os.path.join(args.outdir, "homd_hmt_hit_burden_qcov99.tsv")
    out_region = os.path.join(args.outdir, "homd_reference_region_hit_burden_qcov99.tsv")
    stats_path = os.path.join(args.outdir, "homd_reference_hit_burden_stats.json")

    log("loading HOMD references")
    refs = parse_homd_fasta(args.homd)
    log(f"HOMD references: {len(refs):,}")

    log("loading derep summary")
    derep_cols = ["derep_id", "n_asvs", "n_bioprojects", "total_reads", "prevalence_samples", "sequencing_regions", "countries", "continents", "statuses", "sites"]
    derep = pd.read_csv(derep_path, sep="\t", usecols=derep_cols)

    log("loading and deduplicating HOMD hits")
    hits = load_hits(hit_path)
    merged = hits.merge(derep, on="derep_id", how="left")
    for col in ["n_asvs", "n_bioprojects", "total_reads", "prevalence_samples", "identity", "query_cov"]:
        merged[col] = pd.to_numeric(merged[col], errors="coerce").fillna(0)

    rows = []
    region_rows = []
    for target, sub in merged.groupby("homd_target", sort=False):
        base = {
            "homd_target": target,
            "any_hit_derep_sequences": len(sub),
            "any_hit_source_asvs": int(sub["n_asvs"].sum()),
            "any_hit_total_reads": float(sub["total_reads"].sum()),
            "any_hit_prevalence_samples": int(sub["prevalence_samples"].sum()),
            "any_hit_max_identity": float(sub["identity"].max()),
            "any_hit_median_identity": float(sub["identity"].median()),
            "any_hit_max_qcov": float(sub["query_cov"].max()),
        }
        for threshold in IDENTITY_THRESHOLDS:
            label = threshold_label(threshold)
            cap = sub[(sub["identity"] >= threshold) & (sub["query_cov"] >= QCOV)]
            base[f"id{label}_qcov99_derep_sequences"] = len(cap)
            base[f"id{label}_qcov99_source_asvs"] = int(cap["n_asvs"].sum())
            base[f"id{label}_qcov99_total_reads"] = float(cap["total_reads"].sum())
            base[f"id{label}_qcov99_prevalence_samples"] = int(cap["prevalence_samples"].sum())
            base[f"id{label}_qcov99_n_bioprojects_sum"] = int(cap["n_bioprojects"].sum())
        rows.append(base)

        for _, row in sub.iterrows():
            for region in split_values(row["sequencing_regions"]):
                rec = {"homd_target": target, "sequencing_region": region, "identity": row["identity"], "query_cov": row["query_cov"],
                       "n_asvs": row["n_asvs"], "total_reads": row["total_reads"], "prevalence_samples": row["prevalence_samples"], "n_bioprojects": row["n_bioprojects"]}
                region_rows.append(rec)

    burden = refs.merge(pd.DataFrame(rows), on="homd_target", how="left")
    count_cols = [c for c in burden.columns if c.startswith("any_hit_") or re.match(r"id\d+p?\d*_qcov99_", c)]
    for col in count_cols:
        burden[col] = burden[col].fillna(0)
    burden = burden.sort_values("id99p0_qcov99_total_reads", ascending=False)
    burden.to_csv(out_ref, sep="\t", index=False)

    hmt = burden.groupby(["homd_hmt", "homd_species"], dropna=False).agg(
        homd_reference_sequences=("homd_target", "count"),
        any_hit_reference_sequences=("any_hit_derep_sequences", lambda x: int((x > 0).sum())),
        id99_qcov99_reference_sequences=("id99p0_qcov99_derep_sequences", lambda x: int((x > 0).sum())),
        any_hit_derep_sequences=("any_hit_derep_sequences", "sum"),
        id100_qcov99_derep_sequences=("id100p0_qcov99_derep_sequences", "sum"),
        id99p5_qcov99_derep_sequences=("id99p5_qcov99_derep_sequences", "sum"),
        id99_qcov99_derep_sequences=("id99p0_qcov99_derep_sequences", "sum"),
        id97_qcov99_derep_sequences=("id97p0_qcov99_derep_sequences", "sum"),
        id95_qcov99_derep_sequences=("id95p0_qcov99_derep_sequences", "sum"),
        id90_qcov99_derep_sequences=("id90p0_qcov99_derep_sequences", "sum"),
        id85_qcov99_derep_sequences=("id85p0_qcov99_derep_sequences", "sum"),
        id80_qcov99_derep_sequences=("id80p0_qcov99_derep_sequences", "sum"),
        id99_qcov99_source_asvs=("id99p0_qcov99_source_asvs", "sum"),
        id99_qcov99_total_reads=("id99p0_qcov99_total_reads", "sum"),
        id99_qcov99_prevalence_samples=("id99p0_qcov99_prevalence_samples", "sum"),
    ).reset_index().sort_values("id99_qcov99_total_reads", ascending=False)
    hmt.to_csv(out_hmt, sep="\t", index=False)

    if region_rows:
        r = pd.DataFrame(region_rows)
        region_summary_parts = []
        for threshold in IDENTITY_THRESHOLDS:
            label = threshold_label(threshold)
            cap = r[(r["identity"] >= threshold) & (r["query_cov"] >= QCOV)]
            tmp = cap.groupby(["homd_target", "sequencing_region"], dropna=False).agg(
                derep_sequences=("homd_target", "count"),
                source_asvs=("n_asvs", "sum"),
                total_reads=("total_reads", "sum"),
                prevalence_samples=("prevalence_samples", "sum"),
                n_bioprojects_sum=("n_bioprojects", "sum"),
            ).reset_index()
            tmp.insert(2, "identity_threshold", threshold)
            region_summary_parts.append(tmp)
        region_summary = pd.concat(region_summary_parts, ignore_index=True)
        region_summary = region_summary.merge(refs, on="homd_target", how="left")
        region_summary = region_summary.sort_values(["identity_threshold", "total_reads"], ascending=[False, False])
        region_summary.to_csv(out_region, sep="\t", index=False)

    stats = {
        "homd_reference_sequences": int(len(refs)),
        "references_with_any_hit": int((burden["any_hit_derep_sequences"] > 0).sum()),
        "references_with_id99_qcov99_hit": int((burden["id99p0_qcov99_derep_sequences"] > 0).sum()),
        "references_with_id100_qcov99_hit": int((burden["id100p0_qcov99_derep_sequences"] > 0).sum()),
        "outputs": {"reference_burden": out_ref, "hmt_burden": out_hmt, "reference_region_burden": out_region},
    }
    with open(stats_path, "w") as handle:
        json.dump(stats, handle, indent=2)
    log(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
