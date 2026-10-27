#!/usr/bin/env python3
"""
Step 04 — Candidate selection: HOMD-uncaptured and non-chimeric.

Reads the chimera-annotated dereplication table (step 03) and keeps the
sequences that are both absent from eHOMD (captured = >=97 % identity and
>=99 % query coverage) and not called chimeric by uchime_ref.

Input : 03_chimera_removal/results/03_derep_chimera_annotated.tsv
Output: 04_candidate_selection/results/04_candidates.tsv
"""
import os
from pathlib import Path

import pandas as pd

PIPE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline")
STEP = PIPE / "04_candidate_selection"
RESULTS = STEP / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

log = lambda s: print(f"[04_candidates] {s}", flush=True)

CHIMERA_TSV = PIPE / "03_chimera_removal" / "results" / "03_derep_chimera_annotated.tsv"

# ---- Load chimera-annotated dereplication table ----
log(f"Loading {CHIMERA_TSV}")
df = pd.read_csv(CHIMERA_TSV, sep="\t", low_memory=False)
log(f"  Total derep sequences: {len(df):,}")

for col in ["identity", "query_cov", "n_bioprojects", "n_asvs", "total_reads", "prevalence_samples"]:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)

if "homd_uncaptured" in df.columns:
    df["homd_uncaptured"] = df["homd_uncaptured"].astype(str).str.lower().map(
        {"true": True, "false": False}).fillna(False)

if "is_chimera" in df.columns:
    df["is_chimera"] = df["is_chimera"].astype(bool)

# ---- Summary ----
log("=== Chimera removal summary ===")
log(f"Total derep: {len(df):,}")

n_uncaptured = df["homd_uncaptured"].sum() if "homd_uncaptured" in df.columns else 0
log(f"HOMD-uncaptured: {n_uncaptured:,} ({n_uncaptured/len(df)*100:.1f}%)")

n_chimera = df["is_chimera"].sum() if "is_chimera" in df.columns else 0
log(f"Chimeric: {n_chimera:,} ({n_chimera/len(df)*100:.1f}%)")

mask_uncaptured = df["homd_uncaptured"] if "homd_uncaptured" in df.columns else True
mask_nonchim = ~df["is_chimera"] if "is_chimera" in df.columns else True
candidates = df[mask_uncaptured & mask_nonchim].copy()
log(f"Candidates (HOMD-uncaptured AND non-chimeric): {len(candidates):,}")

# ---- Save ----
outpath = RESULTS / "04_candidates.tsv"
candidates.to_csv(outpath, sep="\t", index=False)
log(f"Saved: {outpath} ({os.path.getsize(outpath):,} bytes, {len(candidates):,} rows)")
log("Done!")
