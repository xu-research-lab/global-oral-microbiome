#!/usr/bin/env python3
import csv
import re
from collections import defaultdict
from pathlib import Path


ROOT = Path("/beegfs/huangxiaochang/global_oral/new_analysis")
ORAL_BASE = Path("/beegfs/dongbiao/oral")

RAW_READS_PATH = ROOT / "raw_reads_stats.tsv"
SAMPLE_TABLE_PATH = ROOT / "sample_reads_asv_chim_perc.tsv"

STUDY_CONFIG = {
    "OEP000837": {"paired": False},
    "PRJEB12953": {"paired": False},
    "PRJEB48620": {"paired": False},
    "PRJEB51725": {"paired": True},
    "PRJNA306560": {"paired": True},
    "PRJNA330533": {"paired": False},
    "PRJNA419259": {"paired": True},
    "PRJNA436183": {"paired": False},
    "PRJNA598080": {"paired": True},
    "PRJNA780671": {"paired": True},
    "qiita_11259": {"paired": False},
    "qiita_11404": {"paired": False},
}


def normalize_sample(sample_id: str) -> str:
    sample_id = sample_id.strip()
    sample_id = re.sub(r"\.gz$", "", sample_id, flags=re.I)
    sample_id = re.sub(r"\.(fastq|fq|fasta|fa)$", "", sample_id, flags=re.I)
    sample_id = re.sub(r"(_R?[12])$", "", sample_id, flags=re.I)
    return sample_id


def parse_summary_reads(study_id: str, paired: bool) -> dict:
    summary_path = ORAL_BASE / study_id / "results" / "summary.tsv"
    reads = {}
    with summary_path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader, None)
        for row in reader:
            if not row or not row[0]:
                continue
            sample_id = normalize_sample(row[0])
            input_reads = int(float(row[1]))
            total_reads = input_reads * 2 if paired else input_reads
            reads[sample_id] = total_reads
    return reads


def load_backfill_reads() -> dict:
    backfill = {}
    for study_id, config in STUDY_CONFIG.items():
        sample_reads = parse_summary_reads(study_id, config["paired"])
        for sample_id, total_reads in sample_reads.items():
            backfill[(study_id, sample_id)] = total_reads
    return backfill


def rewrite_raw_reads(backfill_reads: dict) -> dict:
    existing = {}
    with RAW_READS_PATH.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            existing[(row["Project"], row["Sample"])] = row["Total_Reads"]

    for (study_id, sample_id), total_reads in backfill_reads.items():
        existing[(study_id, sample_id)] = str(total_reads)

    ordered_rows = sorted(existing.items(), key=lambda item: (item[0][0], item[0][1]))
    with RAW_READS_PATH.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Project", "Sample", "Total_Reads"])
        for (study_id, sample_id), total_reads in ordered_rows:
            writer.writerow([study_id, sample_id, total_reads])

    study_counts = defaultdict(int)
    for study_id, _sample_id in backfill_reads:
        study_counts[study_id] += 1
    return dict(study_counts)


def rewrite_sample_table(backfill_reads: dict) -> tuple:
    rows = []
    updated = defaultdict(int)
    unresolved = defaultdict(int)

    with SAMPLE_TABLE_PATH.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames
        for row in reader:
            study_id = row["bioproject_id"]
            sample_id = normalize_sample(row["sample_id"])
            key = (study_id, sample_id)
            if key in backfill_reads:
                row["total_reads"] = str(backfill_reads[key])
                updated[study_id] += 1
            elif study_id in STUDY_CONFIG and row["total_reads"].strip().upper() == "NA":
                unresolved[study_id] += 1
            rows.append(row)

    with SAMPLE_TABLE_PATH.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    return dict(updated), dict(unresolved)


def main() -> None:
    backfill_reads = load_backfill_reads()
    raw_counts = rewrite_raw_reads(backfill_reads)
    updated_counts, unresolved_counts = rewrite_sample_table(backfill_reads)

    print("Backfilled raw_reads_stats.tsv:")
    for study_id in sorted(raw_counts):
        print(f"  {study_id}: {raw_counts[study_id]} samples")

    print("Updated sample_reads_asv_chim_perc.tsv:")
    for study_id in sorted(updated_counts):
        unresolved = unresolved_counts.get(study_id, 0)
        print(f"  {study_id}: updated={updated_counts[study_id]}, unresolved={unresolved}")


if __name__ == "__main__":
    main()
