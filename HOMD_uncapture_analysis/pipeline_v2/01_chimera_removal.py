#!/usr/bin/env python3
"""
pipeline_v2 Step 01: Chimera removal (before support floor).
Reuses existing vsearch --uchime_ref results from the original pipeline
(exact_derep_chimera_annotated.tsv), which ran on ALL 2,844,174 derep sequences
against HOMD + GG2 + SILVA combined reference.

Output: 01_chimera_filtered.tsv — non-chimeric derep sequences
"""
import pandas as pd
import os
from pathlib import Path

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
OUTDIR = BASE / "pipeline_v2"
os.makedirs(OUTDIR, exist_ok=True)

log = lambda s: print(f"[01_chimera] {s}", flush=True)

# ---- Load existing chimera results (no support floor filter) ----
CHIMERA_TSV = BASE / "chimera_removal/exact_derep_chimera_annotated.tsv"

log(f"Loading {CHIMERA_TSV}")
df = pd.read_csv(CHIMERA_TSV, sep="\t", low_memory=False)
log(f"  Total derep sequences: {len(df):,}")

# Convert key columns
for col in ["identity", "query_cov", "n_bioprojects", "n_asvs", "total_reads", "prevalence_samples"]:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)

if "homd_uncaptured" in df.columns:
    df["homd_uncaptured"] = df["homd_uncaptured"].astype(str).str.lower().map(
        {"true": True, "false": False}).fillna(False)

if "is_chimera" in df.columns:
    df["is_chimera"] = df["is_chimera"].astype(bool)

# ---- Summary ----
log(f"\n=== Chimera Removal Summary ===")
log(f"Total derep: {len(df):,}")

n_uncaptured = df["homd_uncaptured"].sum() if "homd_uncaptured" in df.columns else 0
log(f"HOMD-uncaptured: {n_uncaptured:,} ({n_uncaptured/len(df)*100:.1f}%)")

n_chimera = df["is_chimera"].sum() if "is_chimera" in df.columns else 0
log(f"Chimeric: {n_chimera:,} ({n_chimera/len(df)*100:.1f}%)")

# Non-chimeric + HOMD-uncaptured = our candidates
mask_uncaptured = df["homd_uncaptured"] if "homd_uncaptured" in df.columns else True
mask_nonchim = ~df["is_chimera"] if "is_chimera" in df.columns else True
candidates = df[mask_uncaptured & mask_nonchim].copy()
log(f"\nCandidates (HOMD-uncaptured AND non-chimeric): {len(candidates):,}")

# ---- Save ----
outpath = OUTDIR / "01_chimera_filtered.tsv"
candidates.to_csv(outpath, sep="\t", index=False)
log(f"Saved: {outpath} ({os.path.getsize(outpath):,} bytes, {len(candidates):,} rows)")
log("Done!")
