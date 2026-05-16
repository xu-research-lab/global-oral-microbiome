#!/usr/bin/env python3
"""
Filter metadata to retain only samples present in a BIOM table,
then save the result as a QIIME2-compatible TSV metadata file.

Usage:
    python filter_metadata.py <metadata_path> <biom_path> <output_path>

Arguments:
    metadata_path  Path to the original metadata file (CSV or TSV)
    biom_path      Path to the BIOM table
    output_path    Path for the filtered metadata output (TSV)
"""

import sys
import argparse
import pandas as pd
import biom

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Intersect metadata samples with BIOM table samples "
            "and export a QIIME2-compatible metadata TSV."
        )
    )
    parser.add_argument("metadata_path", help="Path to the original metadata CSV/TSV file")
    parser.add_argument("biom_path", help="Path to the BIOM table file")
    parser.add_argument("output_path", help="Path for the filtered metadata output (TSV)")
    return parser.parse_args()

args = parse_args()
table = biom.load_table(args.biom_path)
sid = table.ids(axis="sample")
metadata = pd.read_csv(args.metadata_path, dtype={"sample_id": str})
metadata = metadata.loc[metadata.site_group.values == "oral"]
sid = np.intersect1d(sid, metadata.sample_id.values)

### save oral biom
table.filter(sid, axis="sample")
table.remove_empty()
with biom.util.biom_open(f"{args.output_path}/DADA2.biom", 'w') as f:
    table.to_hdf5(f, generated_by="fid_replacement")
    
### oral biom fasta file
seqs = table.ids(axis='observation')
with open(f"{args.output_path}/rep-seqs.fasta", 'w') as f:
    for seq in seqs:
        f.write(f'>{seq}\n{seq}\n')

### save as QIIME2-compatible TSV
sid_set = set(sid)
metadata = metadata[metadata["sample_id"].isin(sid_set)]
metadata = metadata[["sample_id", "bioproject", 'subject_id', 'age',
       'age_info', 'bmi', 'sex', 'platform', 'sequencing_region',
       'pcr_primers', 'country', 'site', 'site_group', 'disease_study',
       'status', 'title', 'public_time', 'matched_samples', 'is_time_series']]
metadata = metadata.set_index("#SampleID")
metadata.to_csv(f"{args.output_path}/metadata.tsv", sep="\t")
