#!/usr/bin/env python3

import argparse
import pandas as pd
import biom
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Filter merged metadata by project, keep oral samples, "
            "intersect with BIOM table samples, and export a QIIME2-compatible metadata TSV."
        )
    )
    parser.add_argument("metadata_path", help="Path to the merged metadata CSV/TSV file")
    parser.add_argument("biom_path", help="Path to the BIOM table file")
    parser.add_argument("output_path", help="Path for the filtered metadata output directory")
    parser.add_argument("project_id", help="Project/study ID used to subset merged metadata")
    return parser.parse_args()


args = parse_args()
table = biom.load_table(args.biom_path)
sid = table.ids(axis="sample")

metadata = pd.read_csv(
    args.metadata_path,
    dtype={"sample_id": str, "bioproject": str},
    low_memory=False,
)
metadata = metadata.loc[metadata["bioproject"].astype(str) == str(args.project_id)]
metadata = metadata.loc[metadata["site_group_cleaned"] == "Oral"]
sid = np.intersect1d(sid, metadata["sample_id"].astype(str).values)

table.filter(sid, axis="sample", inplace=True)
table.remove_empty(inplace=True)

with biom.util.biom_open(f"{args.output_path}/DADA2.biom", "w") as f:
    table.to_hdf5(f, generated_by="fid_replacement")

seqs = table.ids(axis="observation")
with open(f"{args.output_path}/rep-seqs.fasta", "w") as f:
    for seq in seqs:
        f.write(f">{seq}\n{seq}\n")

sid_set = set(sid)
metadata = metadata[metadata["sample_id"].isin(sid_set)]
metadata = metadata[
    [
        "sample_id",
        "bioproject",
        "subject_id",
        "age",
        "age_info",
        "bmi",
        "sex",
        "platform",
        "sequencing_region",
        "pcr_primers",
        "country",
        "site",
        "site_group_cleaned",
        "disease_study",
        "status",
        "title",
        "public_time",
        "matched_samples",
        "is_time_series",
    ]
]
metadata = metadata.rename(columns={"sample_id": "#SampleID"})
metadata = metadata.set_index("#SampleID")
metadata.to_csv(f"{args.output_path}/metadata.tsv", sep="\t")
