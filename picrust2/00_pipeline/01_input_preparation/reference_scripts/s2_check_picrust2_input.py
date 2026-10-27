#!/usr/bin/env python3
import argparse
import gzip
import json
import math
import sys
from pathlib import Path
from collections import Counter

import numpy as np
from biom import load_table


IUPAC_DNA = set("ACGTUNRYWSKMBDHV.-")


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def read_fasta(fasta_path):
    ids = []
    seqs = {}
    current_id = None
    chunks = []

    with open_text(fasta_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue

            if line.startswith(">"):
                if current_id is not None:
                    seqs.setdefault(current_id, []).append("".join(chunks).upper())

                # PICRUSt2 / FASTA parser 通常取 header 第一个空格前的内容作为 ID
                current_id = line[1:].strip().split()[0]
                ids.append(current_id)
                chunks = []
            else:
                chunks.append(line.strip())

        if current_id is not None:
            seqs.setdefault(current_id, []).append("".join(chunks).upper())

    return ids, seqs


def write_list(path, values):
    with open(path, "w") as w:
        for x in values:
            w.write(str(x) + "\n")


def quantiles(values):
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return {}

    return {
        "min": float(np.min(values)),
        "q25": float(np.quantile(values, 0.25)),
        "median": float(np.median(values)),
        "mean": float(np.mean(values)),
        "q75": float(np.quantile(values, 0.75)),
        "max": float(np.max(values)),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Check whether DADA2 BIOM table and representative FASTA are suitable for PICRUSt2."
    )
    parser.add_argument("--biom", required=True, help="DADA2 feature table in BIOM format.")
    parser.add_argument("--fasta", required=True, help="DADA2 ASV representative sequences FASTA.")
    parser.add_argument("--outdir", required=True, help="Output directory for QC reports.")
    parser.add_argument("--no-fail", action="store_true", help="Do not exit non-zero even if fatal checks fail.")
    args = parser.parse_args()

    biom_path = Path(args.biom)
    fasta_path = Path(args.fasta)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fatal_errors = []
    warnings = []

    # ------------------------------------------------------------
    # 1. File existence
    # ------------------------------------------------------------
    if not biom_path.exists() or biom_path.stat().st_size == 0:
        fatal_errors.append(f"BIOM file missing or empty: {biom_path}")

    if not fasta_path.exists() or fasta_path.stat().st_size == 0:
        fatal_errors.append(f"FASTA file missing or empty: {fasta_path}")

    if fatal_errors:
        summary = {
            "status": "FAIL",
            "fatal_errors": fatal_errors,
            "warnings": warnings,
        }
        (outdir / "picrust2_input_check_summary.json").write_text(
            json.dumps(summary, indent=2) + "\n"
        )
        print(json.dumps(summary, indent=2))
        sys.exit(2)

    # ------------------------------------------------------------
    # 2. Load BIOM
    # ------------------------------------------------------------
    try:
        table = load_table(str(biom_path))
    except Exception as e:
        fatal_errors.append(f"Failed to load BIOM file: {e}")
        summary = {
            "status": "FAIL",
            "fatal_errors": fatal_errors,
            "warnings": warnings,
        }
        (outdir / "picrust2_input_check_summary.json").write_text(
            json.dumps(summary, indent=2) + "\n"
        )
        print(json.dumps(summary, indent=2))
        sys.exit(2)

    feature_ids = list(table.ids(axis="observation"))
    sample_ids = list(table.ids(axis="sample"))

    feature_id_counter = Counter(feature_ids)
    duplicated_feature_ids = sorted([k for k, v in feature_id_counter.items() if v > 1])

    if duplicated_feature_ids:
        fatal_errors.append(f"Duplicated feature IDs in BIOM: {len(duplicated_feature_ids)}")

    # ------------------------------------------------------------
    # 3. Load FASTA
    # ------------------------------------------------------------
    fasta_ids, fasta_seqs_dict = read_fasta(fasta_path)
    fasta_id_counter = Counter(fasta_ids)
    duplicated_fasta_ids = sorted([k for k, v in fasta_id_counter.items() if v > 1])

    if duplicated_fasta_ids:
        fatal_errors.append(f"Duplicated sequence IDs in FASTA: {len(duplicated_fasta_ids)}")

    fasta_unique_ids = set(fasta_ids)
    biom_feature_set = set(feature_ids)

    table_not_in_fasta = sorted(biom_feature_set - fasta_unique_ids)
    fasta_not_in_table = sorted(fasta_unique_ids - biom_feature_set)
    shared_ids = sorted(biom_feature_set & fasta_unique_ids)

    if len(shared_ids) == 0:
        fatal_errors.append("No shared feature IDs between BIOM table and FASTA.")

    if table_not_in_fasta:
        fatal_errors.append(
            f"{len(table_not_in_fasta)} BIOM feature IDs are missing from FASTA. "
            "PICRUSt2 will fail or drop these features."
        )

    if fasta_not_in_table:
        warnings.append(
            f"{len(fasta_not_in_table)} FASTA sequence IDs are not present in BIOM table. "
            "This is usually not fatal, but the FASTA can be filtered to the table."
        )

    # ------------------------------------------------------------
    # 4. Sequence QC
    # ------------------------------------------------------------
    seq_lengths = []
    empty_seq_ids = []
    non_iupac_ids = []

    for fid, seq_list in fasta_seqs_dict.items():
        # 如果重复 ID，前面已经记录 fatal，这里仍然合并检查
        seq = "".join(seq_list).upper()
        seq_lengths.append(len(seq))

        if len(seq) == 0:
            empty_seq_ids.append(fid)

        bad_chars = set(seq) - IUPAC_DNA
        if bad_chars:
            non_iupac_ids.append(fid)

    if empty_seq_ids:
        fatal_errors.append(f"{len(empty_seq_ids)} FASTA records have empty sequences.")

    if non_iupac_ids:
        warnings.append(
            f"{len(non_iupac_ids)} FASTA records contain non-IUPAC DNA characters."
        )

    if seq_lengths:
        seq_len_stats = quantiles(seq_lengths)
    else:
        seq_len_stats = {}

    # 对 16S amplicon 来说太短通常可疑，但不直接 fatal
    very_short_ids = [
        fid for fid, seq_list in fasta_seqs_dict.items()
        if len("".join(seq_list)) < 100
    ]
    if very_short_ids:
        warnings.append(
            f"{len(very_short_ids)} FASTA sequences are shorter than 100 bp. "
            "Check whether these are real ASVs or artifacts."
        )

    # ------------------------------------------------------------
    # 5. Count table QC
    # ------------------------------------------------------------
    matrix = table.matrix_data
    values = matrix.data

    if values.size == 0:
        fatal_errors.append("BIOM table has no non-zero counts.")

    non_integer_values = 0
    negative_values = 0

    if values.size > 0:
        non_integer_values = int(np.sum(~np.isclose(values, np.round(values))))
        negative_values = int(np.sum(values < 0))

    if non_integer_values > 0:
        warnings.append(
            f"{non_integer_values} non-integer abundance values detected. "
            "PICRUSt2 expects counts, not relative abundance."
        )

    if negative_values > 0:
        fatal_errors.append(f"{negative_values} negative values detected in BIOM table.")

    sample_sums = np.asarray(table.sum(axis="sample"), dtype=float)
    feature_sums = np.asarray(table.sum(axis="observation"), dtype=float)

    empty_samples = [
        sid for sid, s in zip(sample_ids, sample_sums)
        if s <= 0
    ]
    empty_features = [
        fid for fid, s in zip(feature_ids, feature_sums)
        if s <= 0
    ]

    if empty_samples:
        fatal_errors.append(f"{len(empty_samples)} samples have zero total reads.")

    if empty_features:
        warnings.append(f"{len(empty_features)} features have zero total reads.")

    # 判断是否像 relative abundance table
    looks_like_ra = False
    if len(sample_sums) > 0:
        n_close_to_1 = int(np.sum((sample_sums > 0.99) & (sample_sums < 1.01)))
        frac_close_to_1 = n_close_to_1 / len(sample_sums)
        if frac_close_to_1 > 0.8:
            looks_like_ra = True
            fatal_errors.append(
                "Table looks like relative abundance: most sample sums are close to 1. "
                "PICRUSt2 should use raw count table, not RA table."
            )

    # ------------------------------------------------------------
    # 6. Write detailed lists
    # ------------------------------------------------------------
    write_list(outdir / "biom_feature_ids.txt", feature_ids)
    write_list(outdir / "fasta_ids.txt", fasta_ids)
    write_list(outdir / "shared_ids.txt", shared_ids)
    write_list(outdir / "table_not_in_fasta.txt", table_not_in_fasta)
    write_list(outdir / "fasta_not_in_table.txt", fasta_not_in_table)
    write_list(outdir / "duplicated_fasta_ids.txt", duplicated_fasta_ids)
    write_list(outdir / "duplicated_biom_feature_ids.txt", duplicated_feature_ids)
    write_list(outdir / "empty_sample_ids.txt", empty_samples)
    write_list(outdir / "empty_feature_ids.txt", empty_features)
    write_list(outdir / "empty_fasta_sequence_ids.txt", empty_seq_ids)
    write_list(outdir / "non_iupac_fasta_ids.txt", non_iupac_ids)
    write_list(outdir / "very_short_fasta_ids_lt100bp.txt", very_short_ids)

    # ------------------------------------------------------------
    # 7. Summary
    # ------------------------------------------------------------
    status = "PASS" if not fatal_errors else "FAIL"

    summary = {
        "status": status,
        "biom": str(biom_path),
        "fasta": str(fasta_path),

        "n_samples": len(sample_ids),
        "n_biom_features": len(feature_ids),
        "n_fasta_records": len(fasta_ids),
        "n_unique_fasta_ids": len(fasta_unique_ids),
        "n_shared_feature_ids": len(shared_ids),
        "n_table_not_in_fasta": len(table_not_in_fasta),
        "n_fasta_not_in_table": len(fasta_not_in_table),

        "total_reads": float(np.sum(sample_sums)) if len(sample_sums) > 0 else 0.0,
        "sample_read_count_stats": quantiles(sample_sums),
        "feature_total_count_stats": quantiles(feature_sums),
        "sequence_length_stats": seq_len_stats,

        "n_empty_samples": len(empty_samples),
        "n_empty_features": len(empty_features),
        "n_duplicate_fasta_ids": len(duplicated_fasta_ids),
        "n_duplicate_biom_feature_ids": len(duplicated_feature_ids),
        "n_empty_fasta_sequences": len(empty_seq_ids),
        "n_non_iupac_fasta_sequences": len(non_iupac_ids),
        "n_very_short_sequences_lt100bp": len(very_short_ids),

        "n_non_integer_abundance_values": non_integer_values,
        "n_negative_abundance_values": negative_values,
        "looks_like_relative_abundance": looks_like_ra,

        "fatal_errors": fatal_errors,
        "warnings": warnings,
    }

    with open(outdir / "picrust2_input_check_summary.json", "w") as w:
        json.dump(summary, w, indent=2)
        w.write("\n")

    # 简洁 TSV summary，方便批量合并
    with open(outdir / "picrust2_input_check_summary.tsv", "w") as w:
        keys = [
            "status",
            "n_samples",
            "n_biom_features",
            "n_fasta_records",
            "n_shared_feature_ids",
            "n_table_not_in_fasta",
            "n_fasta_not_in_table",
            "total_reads",
            "n_empty_samples",
            "n_duplicate_fasta_ids",
            "n_non_integer_abundance_values",
            "looks_like_relative_abundance",
        ]
        w.write("\t".join(keys) + "\n")
        w.write("\t".join(str(summary[k]) for k in keys) + "\n")

    print(json.dumps(summary, indent=2))

    if fatal_errors and not args.no_fail:
        sys.exit(2)


if __name__ == "__main__":
    main()
