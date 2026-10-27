#!/usr/bin/env python3
"""
Step 05 — Organelle removal, part 1: write the candidate FASTA.

All HOMD-uncaptured, non-chimeric candidates are screened against the GG2
organelle reference.  Every candidate is written, so the FASTA is the same
sequence set, in the same order, as the one the contamination firewall
screens later.

Input : 04_candidate_selection/results/04_candidates.tsv
Output: 05_organelle_removal/results/05_candidates.fasta
"""
from pathlib import Path

import pandas as pd

PIPE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline")
STEP = PIPE / "05_organelle_removal"
RESULTS = STEP / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

SRC = PIPE / "04_candidate_selection" / "results" / "04_candidates.tsv"
OUT = RESULTS / "05_candidates.fasta"

df = pd.read_csv(SRC, sep="\t", usecols=["derep_id", "sequence"], low_memory=False)
with open(OUT, "w") as f:
    for row in df.itertuples(index=False):
        f.write(f">{row.derep_id}\n{row.sequence}\n")
print(f"Wrote {len(df):,} sequences to {OUT}")
