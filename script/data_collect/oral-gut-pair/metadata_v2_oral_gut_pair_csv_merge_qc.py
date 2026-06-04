#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
from biom import load_table

# ----------------------- Configure Paths -----------------------
out_dir = "/home/yanglun/YSD/global_oral/global_oral_v2/scripts/filter/filter-pair-0602/qc-pair-0603/results/"
os.makedirs(out_dir, exist_ok=True)

ref_file     = "/home/yanglun/YSD/global_oral/global_oral_v2/metadata-v2/metadata_oral_qc.tsv"
ref_all_file = "/home/yanglun/YSD/global_oral/global_oral_v2/metadata-v2/metadata_all_qc.tsv"
biom_file    = "/home/yanglun/YSD/global_oral/global_oral_v2/metadata-v2/biom_all_qc.biom"

files = {
    "pair": "/home/yanglun/YSD/global_oral/metadata_v2_merge/oral-gut-pair/260602-update/metadata_oral_gut_pair_v2_merge.csv",
    "add":  "/home/yanglun/YSD/global_oral/metadata_v2_merge/oral-gut-pair/260602-update/metadata_pair_v2_additional_oral_site_merge.csv",
    "unk":  "/home/yanglun/YSD/global_oral/metadata_v2_merge/oral-gut-pair/260602-update/metadata_pair_v2_status_unknown_merge.csv"
}

# ------------------ Helper Functions ------------------
def clean_series(series):
    """Cast to string, strip whitespace, and handle NaNs."""
    return series.astype(str).str.strip().replace({'nan': np.nan, 'None': np.nan, '': np.nan})

def load_and_clean_df(filepath, sep=','):
    """Safely load dataframe and clean key columns."""
    df = pd.read_csv(filepath, sep=sep, low_memory=False, dtype=str)
    
    if 'sample_id' in df.columns:
        df['sample_id'] = clean_series(df['sample_id'])
    if 'pair_id' in df.columns:
        df['pair_id'] = clean_series(df['pair_id'])
    if 'site_group_cleaned' in df.columns:
        # Capitalize the first letter to prevent "oral" vs "Oral" mismatch issues
        df['site_group_cleaned'] = clean_series(df['site_group_cleaned']).str.capitalize()
        
    return df

# ------------------ Load Reference Files ------------------
print("0. Loading Reference Files...")
df_ref = load_and_clean_df(ref_file, sep='\t')
valid_ids = set(df_ref['sample_id'].dropna().unique())

df_ref_all = load_and_clean_df(ref_all_file, sep='\t')
valid_gut_ids = set(df_ref_all.loc[df_ref_all['site_group_cleaned'] == 'Gut', 'sample_id'].dropna().unique())

print(f"   -> Loaded {len(valid_ids)} valid Oral IDs and {len(valid_gut_ids)} valid Gut IDs.")

# ------------------ Define Processing Functions ------------------
def annotate_df(df, valid_ids_set):
    # Use native Pandas assignment to prevent Numpy from casting NaNs to string "nan"
    df['filter_reason'] = np.nan
    df.loc[~df['sample_id'].isin(valid_ids_set), 'filter_reason'] = "QC-filtered"
    return df

def write_annotation(df, base_name):
    df.to_csv(os.path.join(out_dir, f"{base_name}_qc_annotation.csv"), index=False)

def filter_by_gut(df, valid_gut_ids_set):
    invalid_mask = (df['site_group_cleaned'] == 'Gut') & (~df['sample_id'].isin(valid_gut_ids_set))
    invalid_pair_ids = set(df.loc[invalid_mask, 'pair_id'].dropna().unique())
    return df[~df['pair_id'].isin(invalid_pair_ids)].copy()

def keep_pair_ids_with_two_rows(df):
    pair_counts = df['pair_id'].value_counts()
    valid_pair_ids = set(pair_counts[pair_counts == 2].index)
    return df[df['pair_id'].isin(valid_pair_ids)].copy()

# ------------------ Main Execution Logic ------------------
print("\n1. Loading Target CSV files...")
df_pair = load_and_clean_df(files["pair"])
df_add  = load_and_clean_df(files["add"])
df_unk  = load_and_clean_df(files["unk"])
print(f"   -> Initial rows: pair={len(df_pair)}, add={len(df_add)}, unk={len(df_unk)}")

print("\n2. Annotating and Rescuing...")
df_pair = annotate_df(df_pair, valid_ids)
df_add  = annotate_df(df_add, valid_ids)
df_unk  = annotate_df(df_unk, valid_ids)

# Rescue logic: If the Oral sample in the same pair is retained, clear its filter_reason
rescued_pair_ids = set(df_pair[(df_pair['filter_reason'].isna()) & (df_pair['site_group_cleaned'] == 'Oral')]['pair_id'].dropna().unique())
df_pair.loc[df_pair['pair_id'].isin(rescued_pair_ids), 'filter_reason'] = np.nan

if 'status_cleaned' in df_unk.columns:
    df_unk.loc[df_unk['status_cleaned'].isna(), 'filter_reason'] = np.nan

print("\n3. Writing annotation files...")
write_annotation(df_pair, "metadata_oral_gut_pair_v2_merge")
write_annotation(df_add,  "metadata_pair_v2_additional_oral_site_merge")
write_annotation(df_unk,  "metadata_pair_v2_status_unknown_merge")

print("\n4. Generating clean versions and secondary filtering by Gut...")
df_pair_clean = df_pair[df_pair['filter_reason'].isna()].drop(columns=['filter_reason']).copy()
df_add_clean  = df_add[df_add['filter_reason'].isna()].drop(columns=['filter_reason']).copy()
df_unk_clean  = df_unk[df_unk['filter_reason'].isna()].drop(columns=['filter_reason']).copy()

print(f"   -> Rows after valid_ids cleanup: pair={len(df_pair_clean)}, add={len(df_add_clean)}, unk={len(df_unk_clean)}")

df_pair_clean = filter_by_gut(df_pair_clean, valid_gut_ids)
df_add_clean  = filter_by_gut(df_add_clean, valid_gut_ids)
df_unk_clean  = filter_by_gut(df_unk_clean, valid_gut_ids)
print(f"   -> Rows after Gut validation: pair={len(df_pair_clean)}, add={len(df_add_clean)}, unk={len(df_unk_clean)}")

print("\n5. Extracting Sample IDs directly from BIOM file...")
table = load_table(biom_file)
biom_sample_ids = set([str(x).strip() for x in table.ids(axis='sample')])
print(f"   -> Extracted {len(biom_sample_ids)} valid Sample IDs from BIOM.")

print("\n6. Applying BIOM sample ID filters and final row count logic...")
df_pair_final = df_pair_clean[df_pair_clean['sample_id'].isin(biom_sample_ids)].copy()
df_add_final  = df_add_clean[df_add_clean['sample_id'].isin(biom_sample_ids)].copy()
df_unk_final  = df_unk_clean[df_unk_clean['sample_id'].isin(biom_sample_ids)].copy()

df_pair_final = keep_pair_ids_with_two_rows(df_pair_final)
df_unk_final  = keep_pair_ids_with_two_rows(df_unk_final)

print(f"   -> FINAL ROWS (Ready to save): pair={len(df_pair_final)}, add={len(df_add_final)}, unk={len(df_unk_final)}")

print("\n7. Writing final QC output files...")
df_pair_final.to_csv(os.path.join(out_dir, "metadata_oral_gut_pair_v2_merge_qc.csv"), index=False)
df_add_final.to_csv(os.path.join(out_dir, "metadata_pair_v2_additional_oral_site_merge_qc.csv"), index=False)
df_unk_final.to_csv(os.path.join(out_dir, "metadata_pair_v2_status_unknown_merge_qc.csv"), index=False)

print("\n8. Generating Bioproject list...")
if 'pair_bioproject' in df_pair_final.columns:
    unique_bioprojects = df_pair_final['pair_bioproject'].dropna().unique()
    bioproject_list_path = os.path.join(out_dir, "project_oral_gut_pair_list_qc.txt")
    with open(bioproject_list_path, 'w') as f:
        for bp in unique_bioprojects:
            f.write(f"{bp}\n")

print("\nAll tasks completed successfully!")