#!/usr/bin/env python3
"""
Filter BIOM files: keep only samples present in metadata CSV,
then remove ASVs with prevalence = 0 (zero reads across all samples),
then remove any samples that become empty (zero total reads) as a result.

Usage:
    python filter_biom_by_metadata.py \
        --metadata /home/xuhuan/oral/metadata-v2-merge-oral.csv \
        --input-dir /home/xuhuan/oral/data_biom \
        --output-dir /home/xuhuan/oral/data_biom_oral
"""

import argparse
import os
import sys
from collections import defaultdict

import pandas as pd
from biom import load_table
from biom.util import biom_open


def main():
    parser = argparse.ArgumentParser(description="Filter BIOM files by metadata and remove zero-prevalence samples")
    parser.add_argument("--metadata", required=True, help="Path to metadata CSV file")
    parser.add_argument("--input-dir", required=True, help="Directory containing input BIOM files")
    parser.add_argument("--output-dir", required=True, help="Directory for output BIOM files")
    args = parser.parse_args()

    # Load metadata: build dict {bioproject: set(sample_ids)}
    print("Loading metadata...")
    meta = pd.read_csv(args.metadata, low_memory=False)
    project_samples = defaultdict(set)
    for _, row in meta.iterrows():
        project_samples[row["bioproject"]].add(row["sample_id"])
    print(f"  Loaded {len(meta)} sample entries across {len(project_samples)} projects")

    # Process each biom file
    biom_files = sorted([f for f in os.listdir(args.input_dir) if f.endswith(".biom")])
    print(f"Found {len(biom_files)} BIOM files to process")

    os.makedirs(args.output_dir, exist_ok=True)

    total_input_samples = 0
    total_input_obs = 0
    total_kept_samples = 0
    total_kept_obs = 0
    total_removed_no_metadata = 0
    total_removed_zero_prev_asv = 0
    total_removed_empty_samples = 0
    skipped_no_matching = 0

    for i, fname in enumerate(biom_files, 1):
        project = fname.replace(".biom", "")
        input_path = os.path.join(args.input_dir, fname)
        output_path = os.path.join(args.output_dir, fname)

        # Get valid sample IDs for this project from metadata
        valid_samples = project_samples.get(project, set())

        # Load biom table
        table = load_table(input_path)
        biom_samples = set(table.ids(axis="sample"))
        n_obs = len(table.ids(axis="observation"))
        total_input_samples += len(biom_samples)
        total_input_obs += n_obs

        # Step 1: Keep only samples in metadata
        samples_in_metadata = biom_samples & valid_samples
        removed_no_metadata = biom_samples - valid_samples
        total_removed_no_metadata += len(removed_no_metadata)

        if not samples_in_metadata:
            skipped_no_matching += 1
            if i <= 10 or i == len(biom_files):
                print(f"  [{i}/{len(biom_files)}] {project}: no matching samples, skipped")
            continue

        # Filter to metadata samples
        table = table.filter(ids_to_keep=samples_in_metadata, axis="sample", inplace=False)

        # Step 2: Remove ASVs with prevalence = 0 (zero reads across all samples)
        obs_totals = table.sum(axis="observation")
        obs_ids = table.ids(axis="observation")
        zero_prev_asvs = set()
        for oid, total in zip(obs_ids, obs_totals):
            if total == 0:
                zero_prev_asvs.add(oid)
        total_removed_zero_prev_asv += len(zero_prev_asvs)

        if zero_prev_asvs:
            table = table.filter(ids_to_keep=set(obs_ids) - zero_prev_asvs, axis="observation", inplace=False)

        # Step 3: Remove samples that became empty (zero total reads) after ASV removal
        sample_totals = table.sum(axis="sample")
        sample_ids = table.ids(axis="sample")
        empty_samples = set()
        for sid, total in zip(sample_ids, sample_totals):
            if total == 0:
                empty_samples.add(sid)
        total_removed_empty_samples += len(empty_samples)

        if empty_samples:
            table = table.filter(ids_to_keep=sample_ids - empty_samples, axis="sample", inplace=False)
            # Also remove any ASVs that became zero after removing empty samples
            obs_totals2 = table.sum(axis="observation")
            obs_ids2 = table.ids(axis="observation")
            zero_obs2 = set()
            for oid, total in zip(obs_ids2, obs_totals2):
                if total == 0:
                    zero_obs2.add(oid)
            if zero_obs2:
                table = table.filter(ids_to_keep=set(obs_ids2) - zero_obs2, axis="observation", inplace=False)

        kept_samples = len(table.ids(axis="sample"))
        kept_obs = len(table.ids(axis="observation"))
        total_kept_samples += kept_samples
        total_kept_obs += kept_obs

        # Write output
        with biom_open(output_path, "w") as f:
            table.to_hdf5(f, generated_by="filter_biom_by_metadata.py")

        if i <= 10 or i % 50 == 0 or i == len(biom_files):
            print(f"  [{i}/{len(biom_files)}] {project}: {len(biom_samples)} samples {n_obs} ASVs -> "
                  f"{kept_samples} samples {kept_obs} ASVs "
                  f"(-{len(removed_no_metadata)} no meta, -{len(zero_prev_asvs)} zero-prev ASV, "
                  f"-{len(empty_samples)} empty sample)")

    summary = (
        "\n=== Summary ===\n"
        f"  Total input samples:           {total_input_samples}\n"
        f"  Total input ASVs:              {total_input_obs}\n"
        f"  Removed (no metadata):         {total_removed_no_metadata}\n"
        f"  Removed (zero-prev ASV):       {total_removed_zero_prev_asv}\n"
        f"  Removed (empty sample):        {total_removed_empty_samples}\n"
        f"  Total kept samples:            {total_kept_samples}\n"
        f"  Total kept ASVs:               {total_kept_obs}\n"
        f"  Projects with no match:        {skipped_no_matching}\n"
        f"  Output files written:          {len(biom_files) - skipped_no_matching}\n"
    )
    print(summary)

    # Write summary to log file
    log_path = os.path.join(args.output_dir, "filter_biom_summary.log")
    with open(log_path, "w") as f:
        f.write(summary)
    print(f"Summary log written to: {log_path}")


if __name__ == "__main__":
    main()
