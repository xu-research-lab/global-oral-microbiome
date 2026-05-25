#!/usr/bin/env python3
"""
Remove singleton ASVs (total reads = 1 across all samples)
from each project's BIOM file, generate per-project singleton statistics,
and write filtered BIOM files to output directory.

Usage:
    python remove_singleton_asv.py \
        --input-dir /home/xuhuan/oral/data_biom_oral \
        --output-dir /home/xuhuan/oral/data_biom_singleton
"""

import argparse
import os

import numpy as np
import pandas as pd
from biom import load_table
from biom.util import biom_open


def main():
    parser = argparse.ArgumentParser(description="Remove singleton ASVs and generate statistics")
    parser.add_argument("--input-dir", required=True, help="Directory containing input BIOM files")
    parser.add_argument("--output-dir", required=True, help="Directory for output BIOM files (singleton removed)")
    args = parser.parse_args()

    biom_files = sorted([f for f in os.listdir(args.input_dir) if f.endswith(".biom")])
    print(f"Found {len(biom_files)} BIOM files to process")

    os.makedirs(args.output_dir, exist_ok=True)

    # Per-project stats
    stats = []

    total_input_samples = 0
    total_input_obs = 0
    total_singleton_asv = 0
    total_kept_samples = 0
    total_kept_obs = 0
    total_removed_empty_samples = 0

    for i, fname in enumerate(biom_files, 1):
        project = fname.replace(".biom", "")
        input_path = os.path.join(args.input_dir, fname)
        output_path = os.path.join(args.output_dir, fname)

        # 1. Load table
        table = load_table(input_path)
        n_samples = len(table.ids(axis="sample"))
        n_obs = len(table.ids(axis="observation"))
        total_input_samples += n_samples
        total_input_obs += n_obs

        # 2. Compute ASV total reads using built-in sparse sum (Avoids OOM)
        asv_total_reads = table.sum(axis="observation")
        obs_ids = np.array(table.ids(axis="observation")) # Convert to numpy array for boolean indexing

        # 3. Identify singleton ASVs (total reads = 1)
        singleton_mask = asv_total_reads == 1
        n_singleton = int(singleton_mask.sum())
        total_singleton_asv += n_singleton

        # Compute total reads of singleton ASVs (Each singleton has exactly 1 read, so total reads = count)
        singleton_reads = n_singleton

        # 4. Remove singleton ASVs
        if n_singleton > 0:
            keep_obs = set(obs_ids[~singleton_mask])
            table = table.filter(ids_to_keep=keep_obs, axis="observation", inplace=False)

        # 5. Remove samples that became empty after singleton removal
        sample_totals = table.sum(axis="sample")
        sample_ids = np.array(table.ids(axis="sample"))
        empty_mask = sample_totals == 0
        empty_samples = set(sample_ids[empty_mask])
        
        total_removed_empty_samples += len(empty_samples)

        if empty_samples:
            keep_samples = set(sample_ids[~empty_mask])
            table = table.filter(ids_to_keep=keep_samples, axis="sample", inplace=False)
            
            # Remove any ASVs that became zero after removing empty samples
            obs_totals2 = table.sum(axis="observation")
            obs_ids2 = np.array(table.ids(axis="observation"))
            zero_obs_mask = obs_totals2 == 0
            
            if zero_obs_mask.any():
                keep_obs2 = set(obs_ids2[~zero_obs_mask])
                table = table.filter(ids_to_keep=keep_obs2, axis="observation", inplace=False)

        kept_samples = len(table.ids(axis="sample"))
        kept_obs = len(table.ids(axis="observation"))
        total_kept_samples += kept_samples
        total_kept_obs += kept_obs

        # Per-project stats
        stats.append({
            "project": project,
            "input_samples": n_samples,
            "input_asvs": n_obs,
            "singleton_asvs": n_singleton,
            "singleton_reads": singleton_reads,
            "singleton_asv_pct": round(n_singleton / n_obs * 100, 2) if n_obs > 0 else 0,
            "empty_samples_after_removal": len(empty_samples),
            "output_samples": kept_samples,
            "output_asvs": kept_obs,
        })

        # 6. Write output
        with biom_open(output_path, "w") as f:
            table.to_hdf5(f, generated_by="remove_singleton_asv.py")

        if i <= 10 or i % 50 == 0 or i == len(biom_files):
            print(f"  [{i}/{len(biom_files)}] {project}: {n_samples} samples {n_obs} ASVs -> "
                  f"{kept_samples} samples {kept_obs} ASVs "
                  f"(-{n_singleton} singleton ASV, -{len(empty_samples)} empty sample)")

    # 7. Write per-project stats CSV
    stats_df = pd.DataFrame(stats)
    stats_path = os.path.join(args.output_dir, "singleton_stats.csv")
    stats_df.to_csv(stats_path, index=False)
    print(f"\nPer-project stats written to: {stats_path}")

    # 8. Summary
    summary = (
        "\n=== Summary ===\n"
        f"  Total input samples:              {total_input_samples}\n"
        f"  Total input ASVs:                 {total_input_obs}\n"
        f"  Removed singleton ASVs:           {total_singleton_asv}\n"
        f"  Removed empty samples:            {total_removed_empty_samples}\n"
        f"  Total kept samples:               {total_kept_samples}\n"
        f"  Total kept ASVs:                  {total_kept_obs}\n"
        f"  Output files written:             {len(biom_files)}\n"
    )
    print(summary)

    # Write summary log
    log_path = os.path.join(args.output_dir, "remove_singleton_summary.log")
    with open(log_path, "w") as f:
        f.write(summary)
    print(f"Summary log written to: {log_path}")

if __name__ == "__main__":
    main()