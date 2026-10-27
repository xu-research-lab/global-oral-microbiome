#!/usr/bin/env python3
"""
Step 01c — Collect ASVs and abundance for the eHOMD-uncaptured analysis.

Key features of this collection:
- Per-amplicon-region stratification (V3-V4, V4, V1-V3, etc.)
- No prevalence pre-filter (keep the rare tail)
- Output to pipeline/01_data_assembly/results/

Inputs:
- metadata-v2-merge-oral.csv (updated oral sample table)
- /beegfs/dongbiao/oral/*/results/ASVs.fa + ASVs_counts.tsv

Outputs:
- 01_derep_summary.tsv          : exact-derep sequence-level summary
- 01_asv_abundance_metadata.tsv : per-ASV abundance + metadata + amplicon_region
- 01_derep_representatives.fa    : unique sequences
- 01_project_input_summary.tsv   : per-project processing summary
- 01_collect_stats.json          : run statistics
- by_region/01_derep_summary_*.tsv : the same summary within each amplicon region
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple

import pandas as pd

DEFAULT_METADATA = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/data/metadata-v2-merge-oral.csv"
DEFAULT_ASV_ROOT = "/beegfs/dongbiao/oral"
DEFAULT_OUTDIR = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline/01_data_assembly/results"


@dataclass
class AsvRecord:
    renamed_id: str
    bioproject: str
    original_id: str
    sequence: str
    length: int
    total_reads: float
    prevalence_samples: int
    n_samples_project: int
    amplicon_region: str  # NEW: per-region tracking


def log(message: str) -> None:
    print(f"[collect] {message}", flush=True)


def read_fasta(path: str) -> Iterable[Tuple[str, str]]:
    header = None
    chunks: List[str] = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks).upper()
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            yield header, "".join(chunks).upper()


def write_fasta_record(handle, header: str, sequence: str, width: int = 80) -> None:
    handle.write(f">{header}\n")
    for start in range(0, len(sequence), width):
        handle.write(sequence[start : start + width] + "\n")


def summarize_metadata(metadata_csv: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load metadata and build per-project summary."""
    df = pd.read_csv(metadata_csv, low_memory=False)
    required = {"bioproject", "sample_id"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"metadata is missing required columns: {sorted(missing)}")

    df["bioproject"] = df["bioproject"].astype(str).str.strip()
    df = df[df["bioproject"].ne("") & df["bioproject"].ne("nan")].copy()

    # Determine amplicon region per sample
    region_col = "sequencing_region_cleaned"
    if region_col not in df.columns:
        df[region_col] = "NA"
    df[region_col] = df[region_col].fillna("NA").astype(str)

    columns = [
        "bioproject",
        region_col,
        "platform_cleaned",
        "country_cleaned",
        "continent_group",
        "site_cleaned",
        "disease_study_cleaned",
        "status_cleaned",
    ]
    for col in columns:
        if col not in df.columns:
            df[col] = "NA"
        df[col] = df[col].fillna("NA").astype(str).replace({"": "NA", "nan": "NA"})

    project_rows = []
    for bp, sub in df.groupby("bioproject", sort=True):
        row = {"bioproject": bp, "metadata_samples": int(sub["sample_id"].nunique())}
        for col in columns[1:]:
            values = sorted(v for v in sub[col].dropna().unique() if v and v != "nan")
            row[col] = ";".join(values) if values else "NA"
        project_rows.append(row)
    project_meta = pd.DataFrame(project_rows)
    return df, project_meta


def load_project_biom(biom_path: str) -> Tuple[Dict[str, Tuple[float, int]], int, List[str]]:
    """Load DADA2.biom -> {sequence: (total_reads, prevalence)}, sample_count, sample_ids."""
    from biom import load_table
    import numpy as np
    table = load_table(biom_path)
    obs_ids = list(map(str, table.ids(axis='observation')))
    sample_ids = list(map(str, table.ids(axis='sample')))
    sample_count = len(sample_ids)
    mat = table.matrix_data.tocsr()
    row_sums = np.asarray(mat.sum(axis=1)).ravel()
    prevalence = np.diff(mat.indptr)
    abundance: Dict[str, Tuple[float, int]] = {}
    for i, obs_id in enumerate(obs_ids):
        if row_sums[i] > 0:
            abundance[obs_id] = (float(row_sums[i]), int(prevalence[i]))
    return abundance, sample_count, sample_ids


def collect(args: argparse.Namespace) -> None:
    os.makedirs(args.outdir, exist_ok=True)
    md_df, project_meta = summarize_metadata(args.metadata)
    bioprojects = sorted(project_meta["bioproject"].unique())
    project_meta_map = project_meta.set_index("bioproject").to_dict("index")

    # Build sample -> region lookup
    sample_region = dict(zip(md_df["sample_id"].astype(str),
                             md_df.get("sequencing_region_cleaned", pd.Series(["NA"] * len(md_df)))))
    # Also build project -> dominant region
    bp_region = {}
    for bp, sub in md_df.groupby("bioproject"):
        regions = sub.get("sequencing_region_cleaned", pd.Series(["NA"] * len(sub)))
        dominant = regions.value_counts().index[0] if len(regions) > 0 else "NA"
        bp_region[bp] = str(dominant)

    log(f"retained bioprojects: {len(bioprojects)}")

    # Output paths
    asv_table_path = os.path.join(args.outdir, "01_asv_abundance_metadata.tsv")
    derep_fasta = os.path.join(args.outdir, "01_derep_representatives.fa")
    derep_table_path = os.path.join(args.outdir, "01_derep_summary.tsv")
    project_summary_path = os.path.join(args.outdir, "01_project_input_summary.tsv")
    stats_path = os.path.join(args.outdir, "01_collect_stats.json")

    seq_to_records: Dict[str, List[AsvRecord]] = defaultdict(list)
    # Also track per-region separately
    seq_to_records_by_region: Dict[str, Dict[str, List[AsvRecord]]] = defaultdict(lambda: defaultdict(list))

    project_summaries = []
    total_asvs = 0
    total_reads = 0.0
    missing_files = []

    asv_fields = [
        "renamed_asv_id", "bioproject", "original_asv_id", "sequence", "length",
        "total_reads", "prevalence_samples", "n_samples_project",
        "amplicon_region",  # NEW
        "metadata_samples", "sequencing_region_cleaned", "platform_cleaned",
        "country_cleaned", "continent_group", "site_cleaned",
        "disease_study_cleaned", "status_cleaned",
    ]

    with open(asv_table_path, "w", newline="") as table_out:
        writer = csv.DictWriter(table_out, fieldnames=asv_fields, delimiter="\t")
        writer.writeheader()

        for index, bp in enumerate(bioprojects, start=1):
            biom_path = os.path.join(args.asv_root, bp, "results", "DADA2.biom")
            if not os.path.exists(biom_path):
                missing_files.append({
                    "bioproject": bp,
                    "biom_exists": False,
                })
                continue

            abundance, sample_count, biom_sample_ids = load_project_biom(biom_path)
            region = bp_region.get(bp, "NA")
            project_meta_values = project_meta_map[bp]
            project_asvs = 0
            project_reads = 0.0

            for seq_idx, (sequence, (total_reads_asv, prevalence)) in enumerate(
                sorted(abundance.items(), key=lambda x: -x[1][0]), start=1
            ):
                if total_reads_asv <= 0:
                    continue
                original_id = f"ASV_{seq_idx}"
                renamed_id = f"{bp}|{original_id}"
                record = AsvRecord(
                    renamed_id=renamed_id, bioproject=bp,
                    original_id=original_id, sequence=sequence,
                    length=len(sequence), total_reads=total_reads_asv,
                    prevalence_samples=prevalence,
                    n_samples_project=sample_count,
                    amplicon_region=region,
                )
                seq_to_records[sequence].append(record)
                seq_to_records_by_region[region][sequence].append(record)

                writer.writerow({
                    "renamed_asv_id": renamed_id, "bioproject": bp,
                    "original_asv_id": original_id, "sequence": sequence,
                    "length": len(sequence),
                    "total_reads": f"{total_reads_asv:.6g}",
                    "prevalence_samples": prevalence,
                    "n_samples_project": sample_count,
                    "amplicon_region": region,
                    **{k: project_meta_values.get(k, "NA") for k in asv_fields[11:]},
                })

                total_asvs += 1
                project_asvs += 1
                total_reads += total_reads_asv
                project_reads += total_reads_asv

            project_summaries.append({
                "bioproject": bp, "amplicon_region": region,
                "asvs_in_fasta": project_asvs,
                "count_table_features": len(abundance),
                "count_table_samples": sample_count,
                "total_reads": f"{project_reads:.6g}",
            })
            if index % 50 == 0:
                log(f"processed {index}/{len(bioprojects)} projects; ASVs={total_asvs:,}")

    for missing in missing_files:
        project_summaries.append({
            "bioproject": missing["bioproject"],
            "amplicon_region": bp_region.get(missing["bioproject"], "NA"),
            "asvs_in_fasta": 0, "count_table_features": 0,
            "count_table_samples": 0, "total_reads": 0,
        })

    # ---- Write exact-derep summary ----
    derep_fields = [
        "derep_id", "sequence", "length", "n_asvs", "n_bioprojects",
        "total_reads", "prevalence_samples", "representative_asv_id",
        "bioprojects", "amplicon_regions",  # NEW
        "sequencing_regions", "platforms", "countries", "continents",
        "sites", "disease_studies", "statuses",
    ]

    with open(derep_fasta, "w") as fasta_out, open(derep_table_path, "w", newline="") as table_out:
        writer = csv.DictWriter(table_out, fieldnames=derep_fields, delimiter="\t")
        writer.writeheader()
        for idx, (sequence, records) in enumerate(
            sorted(seq_to_records.items(),
                   key=lambda item: (-sum(r.total_reads for r in item[1]), item[0])),
            start=1,
        ):
            derep_id = f"DEREP_{idx}"
            rep = max(records, key=lambda r: (r.total_reads, r.prevalence_samples, r.renamed_id))
            write_fasta_record(fasta_out, derep_id, sequence)
            meta_rows = [project_meta_map.get(r.bioproject, {}) for r in records]

            def joined_values(column: str) -> str:
                values = set()
                for meta in meta_rows:
                    for value in str(meta.get(column, "")).split(";"):
                        if value and value != "nan" and value != "NA":
                            values.add(value)
                return ";".join(sorted(values)) if values else "NA"

            # Unique amplicon regions for this derep sequence
            regions = sorted({r.amplicon_region for r in records if r.amplicon_region != "NA"})

            writer.writerow({
                "derep_id": derep_id, "sequence": sequence,
                "length": len(sequence), "n_asvs": len(records),
                "n_bioprojects": len({r.bioproject for r in records}),
                "total_reads": f"{sum(r.total_reads for r in records):.6g}",
                "prevalence_samples": sum(r.prevalence_samples for r in records),
                "representative_asv_id": rep.renamed_id,
                "bioprojects": ";".join(sorted({r.bioproject for r in records})),
                "amplicon_regions": ";".join(regions) if regions else "NA",
                "sequencing_regions": joined_values("sequencing_region_cleaned"),
                "platforms": joined_values("platform_cleaned"),
                "countries": joined_values("country_cleaned"),
                "continents": joined_values("continent_group"),
                "sites": joined_values("site_cleaned"),
                "disease_studies": joined_values("disease_study_cleaned"),
                "statuses": joined_values("status_cleaned"),
            })

    # ---- Per-region derep summaries ----
    for region, region_seqs in seq_to_records_by_region.items():
        if not region_seqs:
            continue
        region_clean = region.replace("/", "_").replace(" ", "_")
        region_table = os.path.join(args.outdir, f"01_derep_summary_{region_clean}.tsv")
        with open(region_table, "w", newline="") as table_out:
            writer = csv.DictWriter(table_out, fieldnames=derep_fields, delimiter="\t")
            writer.writeheader()
            for idx, (sequence, records) in enumerate(
                sorted(region_seqs.items(),
                       key=lambda item: (-sum(r.total_reads for r in item[1]), item[0])),
                start=1,
            ):
                derep_id = f"DEREP_{idx}"
                rep = max(records, key=lambda r: (r.total_reads, r.prevalence_samples, r.renamed_id))
                meta_rows_r = [project_meta_map.get(r.bioproject, {}) for r in records]
                writer.writerow({
                    "derep_id": derep_id, "sequence": sequence,
                    "length": len(sequence), "n_asvs": len(records),
                    "n_bioprojects": len({r.bioproject for r in records}),
                    "total_reads": f"{sum(r.total_reads for r in records):.6g}",
                    "prevalence_samples": sum(r.prevalence_samples for r in records),
                    "representative_asv_id": rep.renamed_id,
                    "bioprojects": ";".join(sorted({r.bioproject for r in records})),
                    "amplicon_regions": region,
                    "sequencing_regions": joined_values("sequencing_region_cleaned"),
                    "platforms": joined_values("platform_cleaned"),
                    "countries": joined_values("country_cleaned"),
                    "continents": joined_values("continent_group"),
                    "sites": joined_values("site_cleaned"),
                    "disease_studies": joined_values("disease_study_cleaned"),
                    "statuses": joined_values("status_cleaned"),
                })

    # ---- Stats ----
    pd.DataFrame(project_summaries).sort_values("bioproject").to_csv(
        project_summary_path, sep="\t", index=False)

    region_counts = {}
    for region, region_seqs in seq_to_records_by_region.items():
        total_asvs_region = sum(len(records) for records in region_seqs.values())
        region_counts[region] = {"asvs": total_asvs_region, "derep_sequences": len(region_seqs)}

    stats = {
        "metadata_bioprojects": len(bioprojects),
        "projects_with_missing_inputs": len(missing_files),
        "total_asvs": total_asvs,
        "exact_derep_sequences": len(seq_to_records),
        "total_reads_matched": total_reads,
        "per_region": region_counts,
        "outputs": {
            "asv_table": asv_table_path,
            "derep_fasta": derep_fasta,
            "derep_table": derep_table_path,
            "project_summary": project_summary_path,
        },
    }
    with open(stats_path, "w") as handle:
        json.dump(stats, handle, indent=2)

    log(f"complete: ASVs={total_asvs:,}; exact derep sequences={len(seq_to_records):,}")
    log(f"regions: {list(region_counts.keys())}")
    for region, counts in sorted(region_counts.items()):
        log(f"  {region}: {counts['asvs']:,} ASVs, {counts['derep_sequences']:,} unique sequences")
    log(f"stats: {stats_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", default=DEFAULT_METADATA)
    parser.add_argument("--asv-root", default=DEFAULT_ASV_ROOT)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    return parser.parse_args()


if __name__ == "__main__":
    collect(parse_args())
