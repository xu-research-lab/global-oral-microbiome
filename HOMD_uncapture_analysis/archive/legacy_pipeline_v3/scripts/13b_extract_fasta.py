#!/usr/bin/env python
"""Extract derep_id + sequence from 02_contamination_flagged.tsv to FASTA."""
import pandas as pd

SRC = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v2/02_contamination_flagged.tsv"
OUT = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3/13b_uncaptured_asvs.fasta"

df = pd.read_csv(SRC, sep="\t", usecols=["derep_id", "sequence"], low_memory=False)
with open(OUT, "w") as f:
    for i, row in df.iterrows():
        f.write(f">{row['derep_id']}\n{row['sequence']}\n")
print(f"Wrote {len(df):,} sequences to {OUT}")
