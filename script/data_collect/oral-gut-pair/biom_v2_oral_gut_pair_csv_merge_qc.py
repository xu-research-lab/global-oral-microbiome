#!/usr/bin/env python3
import os
import pandas as pd
import h5py
from biom import load_table

# ----------------------- 1. Configure Paths -----------------------
# Input CSV paths
csv_paths = [
    "/home/yanglun/YSD/global_oral/global_oral_v2/scripts/filter/filter-pair-0602/qc-pair-0603/results/metadata_oral_gut_pair_v2_merge_qc.csv",
    "/home/yanglun/YSD/global_oral/global_oral_v2/scripts/filter/filter-pair-0602/qc-pair-0603/results/metadata_pair_v2_additional_oral_site_merge_qc.csv",
    "/home/yanglun/YSD/global_oral/global_oral_v2/scripts/filter/filter-pair-0602/qc-pair-0603/results/metadata_pair_v2_status_unknown_merge_qc.csv"
]

# Input/Output BIOM paths
input_biom_path = "/home/yanglun/YSD/global_oral/global_oral_v2/metadata-v2/biom_all_qc.biom"
output_biom_path = "/home/yanglun/YSD/global_oral/global_oral_v2/scripts/filter/filter-pair-0602/qc-pair-0603/results/biom_oral_gut_pair_v2_qc.biom"

# ----------------------- 2. Extract Target Sample IDs -----------------------
print("Step 1: Extracting target sample IDs from CSV files...")
target_sample_ids = set()

for csv_file in csv_paths:
    if os.path.exists(csv_file):
        # Read sample_id as string to avoid type coercion issues
        df = pd.read_csv(csv_file, usecols=['sample_id'], dtype=str, low_memory=False)
        # Clean the strings: drop NaNs and strip whitespace
        cleaned_ids = df['sample_id'].dropna().str.strip().unique()
        target_sample_ids.update(cleaned_ids)
        print(f"  -> Extracted {len(cleaned_ids)} unique IDs from {os.path.basename(csv_file)}")
    else:
        print(f"  -> WARNING: File not found: {csv_file}")

print(f"\nTotal unique target sample IDs combined: {len(target_sample_ids)}")

if len(target_sample_ids) == 0:
    raise ValueError("Error: No sample IDs were extracted. Cannot filter BIOM file.")

# ----------------------- 3. Load and Filter BIOM Table -----------------------
print(f"\nStep 2: Loading input BIOM file (This may take a moment for large files)...")
table = load_table(input_biom_path)
print(f"  -> Original BIOM shape: {table.shape[1]} samples, {table.shape[0]} observations")

# Verify intersection to avoid filtering errors
biom_ids = set(table.ids(axis='sample'))
ids_to_keep = target_sample_ids.intersection(biom_ids)
print(f"  -> {len(ids_to_keep)} out of {len(target_sample_ids)} target IDs actually exist in the BIOM file.")

# Filter function for BIOM table
def keep_samples(values, id_, metadata):
    return id_ in ids_to_keep

print("\nStep 3: Filtering BIOM table...")
# filter(..., inplace=False) returns a new filtered table object containing ONLY target samples
filtered_table = table.filter(keep_samples, axis='sample', inplace=False)

# [MODIFIED] Do NOT remove empty observations to preserve the original feature space
print("  -> Retaining all original observations (even if their counts are zero in these samples)...")
# filtered_table.remove_empty(axis='observation', inplace=True)  <-- 这一行已被移除

print(f"  -> Filtered BIOM shape: {filtered_table.shape[1]} samples, {filtered_table.shape[0]} observations")

# ----------------------- 4. Export Filtered BIOM -----------------------
print("\nStep 4: Saving filtered BIOM file to HDF5 format...")
with h5py.File(output_biom_path, 'w') as f:
    filtered_table.to_hdf5(f, "Filtered by global_oral_qc python script (All observations kept)")

print(f"Success! Filtered BIOM saved to: {output_biom_path}")