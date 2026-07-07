#!/usr/bin/env python3
"""
Collect retained oral ASVs and abundance evidence for region-aware HOMD coverage.

Outputs are ASV-level and exact-sequence-level tables. No mito/chloro filtering and
no 99% clustering are performed here, so the primary HOMD analysis remains at the
original ASV resolution while still being computationally efficient after exact
dereplication.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple

import pandas as pd


DEFAULT_METADATA = "/beegfs/huangxiaochang/global_oral/metadata_v2_merge/metadata-v2-merge-oral.csv"
DEFAULT_ASV_ROOT = "/beegfs/dongbiao/oral"
DEFAULT_OUTDIR = "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/region_aware_pipeline"


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
    df = pd.read_csv(metadata_csv, low_memory=False)
    required = {"bioproject", "sample_id"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"metadata is missing required columns: {sorted(missing)}")

    df["bioproject"] = df["bioproject"].astype(str).str.strip()
    df = df[df["bioproject"].ne("") & df["bioproject"].ne("nan")].copy()

    columns = [
        "bioproject",
        "sequencing_region_cleaned",
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


def load_project_abundance(counts_path: str) -> Tuple[Dict[str, Tuple[float, int]], int]:
    """Load results/ASVs_counts.tsv as ASV_ID -> (total reads, sample prevalence)."""
    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    counts.index = counts.index.astype(str)
    sample_count = counts.shape[1]
    abundance: Dict[str, Tuple[float, int]] = {}
    for asv_id, row in counts.iterrows():
        numeric = pd.to_numeric(row, errors="coerce").fillna(0)
        abundance[asv_id] = (float(numeric.sum()), int((numeric > 0).sum()))
    return abundance, sample_count


def collect(args: argparse.Namespace) -> None:
    os.makedirs(args.outdir, exist_ok=True)
    _, project_meta = summarize_metadata(args.metadata)
    bioprojects = sorted(project_meta["bioproject"].unique())
    project_meta_map = project_meta.set_index("bioproject").to_dict("index")
    log(f"retained bioprojects in metadata: {len(bioprojects)}")

    all_asv_fasta = os.path.join(args.outdir, "all_asvs_renamed.fa")
    asv_table_path = os.path.join(args.outdir, "asv_abundance_metadata.tsv")
    derep_fasta = os.path.join(args.outdir, "exact_derep_representatives.fa")
    derep_table_path = os.path.join(args.outdir, "exact_derep_summary.tsv")
    project_summary_path = os.path.join(args.outdir, "project_input_summary.tsv")
    stats_path = os.path.join(args.outdir, "collect_stats.json")

    seq_to_records: Dict[str, List[AsvRecord]] = defaultdict(list)
    project_summaries = []
    total_asvs = 0
    total_reads = 0.0
    missing_files = []
    counts_unmatched_total = 0

    asv_fields = [
        "renamed_asv_id",
        "bioproject",
        "original_asv_id",
        "sequence",
        "length",
        "total_reads",
        "prevalence_samples",
        "n_samples_project",
        "metadata_samples",
        "sequencing_region_cleaned",
        "platform_cleaned",
        "country_cleaned",
        "continent_group",
        "site_cleaned",
        "disease_study_cleaned",
        "status_cleaned",
    ]

    with open(all_asv_fasta, "w") as fasta_out, open(asv_table_path, "w", newline="") as table_out:
        writer = csv.DictWriter(table_out, fieldnames=asv_fields, delimiter="\t")
        writer.writeheader()

        for index, bp in enumerate(bioprojects, start=1):
            asv_path = os.path.join(args.asv_root, bp, "results", "ASVs.fa")
            counts_path = os.path.join(args.asv_root, bp, "results", "ASVs_counts.tsv")
            if not os.path.exists(asv_path) or not os.path.exists(counts_path):
                missing_files.append(
                    {
                        "bioproject": bp,
                        "asv_fasta_exists": os.path.exists(asv_path),
                        "asvs_counts_exists": os.path.exists(counts_path),
                    }
                )
                continue

            abundance, sample_count = load_project_abundance(counts_path)
            seen_asv_ids = set()
            project_asvs = 0
            project_reads = 0.0
            project_prevalence = 0
            project_meta_values = project_meta_map[bp]

            for original_id, sequence in read_fasta(asv_path):
                seen_asv_ids.add(original_id)
                total_reads_asv, prevalence = abundance.get(original_id, (0.0, 0))
                renamed_id = f"{bp}|{original_id}"
                record = AsvRecord(
                    renamed_id=renamed_id,
                    bioproject=bp,
                    original_id=original_id,
                    sequence=sequence,
                    length=len(sequence),
                    total_reads=total_reads_asv,
                    prevalence_samples=prevalence,
                    n_samples_project=sample_count,
                )
                seq_to_records[sequence].append(record)
                write_fasta_record(fasta_out, renamed_id, sequence)

                row = {
                    "renamed_asv_id": renamed_id,
                    "bioproject": bp,
                    "original_asv_id": original_id,
                    "sequence": sequence,
                    "length": len(sequence),
                    "total_reads": f"{total_reads_asv:.6g}",
                    "prevalence_samples": prevalence,
                    "n_samples_project": sample_count,
                }
                row.update(project_meta_values)
                writer.writerow(row)

                total_asvs += 1
                project_asvs += 1
                total_reads += total_reads_asv
                project_reads += total_reads_asv
                project_prevalence += prevalence

            counts_unmatched = len(set(abundance).difference(seen_asv_ids))
            counts_unmatched_total += counts_unmatched
            project_summaries.append(
                {
                    "bioproject": bp,
                    "asv_fasta_exists": True,
                    "asvs_counts_exists": True,
                    "asvs_in_fasta": project_asvs,
                    "count_table_features": len(abundance),
                    "count_table_features_not_in_fasta": counts_unmatched,
                    "count_table_samples": sample_count,
                    "total_reads_in_matched_asvs": f"{project_reads:.6g}",
                    "summed_asv_sample_prevalence": project_prevalence,
                }
            )
            if index % 25 == 0:
                log(f"processed {index}/{len(bioprojects)} projects; ASVs={total_asvs:,}")

    for missing in missing_files:
        project_summaries.append(
            {
                "bioproject": missing["bioproject"],
                "asv_fasta_exists": missing["asv_fasta_exists"],
                "asvs_counts_exists": missing["asvs_counts_exists"],
                "asvs_in_fasta": 0,
                "count_table_features": 0,
                "count_table_features_not_in_fasta": 0,
                "count_table_samples": 0,
                "total_reads_in_matched_asvs": 0,
                "summed_asv_sample_prevalence": 0,
            }
        )

    derep_fields = [
        "derep_id",
        "sequence",
        "length",
        "n_asvs",
        "n_bioprojects",
        "total_reads",
        "prevalence_samples",
        "representative_asv_id",
        "bioprojects",
        "sequencing_regions",
        "platforms",
        "countries",
        "continents",
        "sites",
        "disease_studies",
        "statuses",
    ]
    with open(derep_fasta, "w") as fasta_out, open(derep_table_path, "w", newline="") as table_out:
        writer = csv.DictWriter(table_out, fieldnames=derep_fields, delimiter="\t")
        writer.writeheader()
        for idx, (sequence, records) in enumerate(
            sorted(seq_to_records.items(), key=lambda item: (-sum(r.total_reads for r in item[1]), item[0])),
            start=1,
        ):
            derep_id = f"DEREP_{idx}"
            rep = max(records, key=lambda r: (r.total_reads, r.prevalence_samples, r.renamed_id))
            write_fasta_record(fasta_out, derep_id, sequence)
            meta_rows = [project_meta_map[r.bioproject] for r in records]

            def joined_project_values(column: str) -> str:
                values = set()
                for meta in meta_rows:
                    for value in str(meta[column]).split(";"):
                        if value and value != "nan":
                            values.add(value)
                return ";".join(sorted(values)) if values else "NA"

            writer.writerow(
                {
                    "derep_id": derep_id,
                    "sequence": sequence,
                    "length": len(sequence),
                    "n_asvs": len(records),
                    "n_bioprojects": len({r.bioproject for r in records}),
                    "total_reads": f"{sum(r.total_reads for r in records):.6g}",
                    "prevalence_samples": sum(r.prevalence_samples for r in records),
                    "representative_asv_id": rep.renamed_id,
                    "bioprojects": ";".join(sorted({r.bioproject for r in records})),
                    "sequencing_regions": joined_project_values("sequencing_region_cleaned"),
                    "platforms": joined_project_values("platform_cleaned"),
                    "countries": joined_project_values("country_cleaned"),
                    "continents": joined_project_values("continent_group"),
                    "sites": joined_project_values("site_cleaned"),
                    "disease_studies": joined_project_values("disease_study_cleaned"),
                    "statuses": joined_project_values("status_cleaned"),
                }
            )

    pd.DataFrame(project_summaries).sort_values("bioproject").to_csv(project_summary_path, sep="\t", index=False)
    stats = {
        "metadata_bioprojects": len(bioprojects),
        "projects_with_missing_inputs": len(missing_files),
        "total_asvs": total_asvs,
        "exact_derep_sequences": len(seq_to_records),
        "total_reads_matched_to_fasta_asvs": total_reads,
        "count_table_features_not_in_fasta_total": counts_unmatched_total,
        "outputs": {
            "all_asv_fasta": all_asv_fasta,
            "asv_table": asv_table_path,
            "derep_fasta": derep_fasta,
            "derep_table": derep_table_path,
            "project_summary": project_summary_path,
        },
    }
    with open(stats_path, "w") as handle:
        json.dump(stats, handle, indent=2)

    log(f"complete: ASVs={total_asvs:,}; exact derep sequences={len(seq_to_records):,}")
    log(f"stats: {stats_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", default=DEFAULT_METADATA)
    parser.add_argument("--asv-root", default=DEFAULT_ASV_ROOT)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    return parser.parse_args()


if __name__ == "__main__":
    collect(parse_args())
