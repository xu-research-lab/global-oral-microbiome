#!/usr/bin/env python3
"""
pipeline_v2 Step 04: Support Floor Filtering & SEPP Placement Reuse

1. Filter clean candidates with ≥2 BP, ≥50 reads
2. Generate FASTA for SEPP
3. Match with original pipeline SEPP placements (reuse where possible)
"""
import pandas as pd
import numpy as np
import os, sys, json
from pathlib import Path

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"
PIPE_ORIG = BASE / "pipeline"
os.makedirs(PIPE_V2, exist_ok=True)
log = lambda s: print(f"[04_support] {s}", flush=True)

# ============================================================
# 1. Load updated contamination results
# ============================================================
log("Loading updated contamination-flagged candidates...")
df = pd.read_csv(PIPE_V2 / "02_contamination_flagged.tsv", sep="\t", low_memory=False)
log(f"  Total candidates: {len(df):,}")

# Ensure numeric
for c in ["n_bioprojects", "total_reads"]:
    df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)

# ============================================================
# 2. Apply support floor: ≥2 BP, ≥50 reads (support ONLY — no contamination filter)
#    Flagged ASVs that reach the floor are RETAINED with their annotations;
#    flag adjudication happens at phylotype aggregation.
# ============================================================
log("\nApplying support floor: ≥2 BP, ≥50 reads (support only)...")

bp_pass = df["n_bioprojects"] >= 2
reads_pass = df["total_reads"] >= 50

passed = bp_pass & reads_pass
n_pass = passed.sum()
log(f"  ≥2 BP: {bp_pass.sum():,}")
log(f"  ≥50 reads: {reads_pass.sum():,}")
log(f"  ≥2 BP + ≥50 reads: {n_pass:,}")

passed_df = df[passed].copy()

# Select output columns (subset of original columns)
# residence columns dropped (2026-08-05): oral-vs-blank residence test removed
out_cols = ["derep_id", "sequence", "n_bioprojects", "total_reads", "prevalence_samples",
            "homd_uncaptured", "is_chimera", "gg2_genus", "gg2_name", "gg2_identity",
            "contamination_flag", "flag_type",
            "bioprojects", "amplicon_regions", "sequencing_regions",
            "continents", "countries"]
avail_cols = [c for c in out_cols if c in passed_df.columns]
passed_df = passed_df[avail_cols]

log(f"  Support floor passed: {len(passed_df):,} ASVs")

# ============================================================
# 3. Save support floor passed table
# ============================================================
out_tsv = PIPE_V2 / "04_support_floor_passed.tsv"
passed_df.to_csv(out_tsv, sep="\t", index=False)
log(f"\nSaved: {out_tsv} ({os.path.getsize(out_tsv):,} bytes)")

# ============================================================
# 4. Generate FASTA for SEPP
# ============================================================
log("\nWriting FASTA for SEPP...")
fa_path = PIPE_V2 / "candidates_for_sepp.fasta"
with open(fa_path, "w") as f:
    for _, row in passed_df.iterrows():
        f.write(f">{row['derep_id']}\n{row['sequence']}\n")
log(f"  Sequences: {len(passed_df):,}")
log(f"  Saved: {fa_path}")

# ============================================================
# 5. Reuse original SEPP placements
# ============================================================
log("\nMatching with original SEPP placements...")
orig_sepp_path = PIPE_ORIG / "05_sepp" / "sepp_placements_parsed.tsv"

if orig_sepp_path.exists():
    orig_sepp = pd.read_csv(orig_sepp_path, sep="\t", low_memory=False)
    log(f"  Original placements: {len(orig_sepp):,}")

    # Normalize column names
    orig_sepp.columns = [c.strip() for c in orig_sepp.columns]

    # Check if query_id column exists
    if "query_id" not in orig_sepp.columns and "derep_id" in orig_sepp.columns:
        orig_sepp["query_id"] = orig_sepp["derep_id"]

    pass_ids = set(passed_df["derep_id"])
    orig_ids = set(orig_sepp["query_id"].astype(str))

    matched = pass_ids & orig_ids
    missing = pass_ids - orig_ids

    log(f"  Reused placements: {len(matched):,} / {len(pass_ids):,} ({len(matched)/max(len(pass_ids),1)*100:.1f}%)")
    log(f"  Missing (need new SEPP): {len(missing):,}")

    if missing:
        log(f"\n  Missing ASV IDs: {sorted(list(missing))[:20]}...")
        missing_fa = PIPE_V2 / "candidates_missing_sepp.fasta"
        with open(missing_fa, "w") as f:
            for _, row in passed_df[passed_df["derep_id"].isin(missing)].iterrows():
                f.write(f">{row['derep_id']}\n{row['sequence']}\n")
        log(f"  Missing FASTA saved: {missing_fa}")

    # Save reused placements
    reused = orig_sepp[orig_sepp["query_id"].astype(str).isin(pass_ids)].copy()
    reused_out = PIPE_V2 / "sepp_placements_reused.tsv"

    # Keep relevant columns (pendant_length needed for phylotype splitting, step 05b)
    keep_cols = ["query_id", "edge_num", "distal_length", "like_weight_ratio", "likelihood", "pendant_length"]
    keep_cols = [c for c in keep_cols if c in reused.columns]
    reused = reused[keep_cols]
    reused.to_csv(reused_out, sep="\t", index=False)
    log(f"  Reused placements saved: {reused_out} ({len(reused):,} rows)")
else:
    log(f"  WARNING: Original SEPP placements not found at {orig_sepp_path}")
    missing = set(passed_df["derep_id"])
    log(f"  All {len(missing):,} ASVs need new SEPP placement")

# ============================================================
# 6. Summary
# ============================================================
log(f"\n=== Support Floor Summary ===")
log(f"Total candidates (after firewall annotation): {len(df):,}")
log(f"≥2 BP + ≥50 reads (support only): {n_pass:,}")
log(f"  clean: {(passed_df['contamination_flag']==False).sum():,} | flagged: {passed_df['contamination_flag'].sum():,}")
log(f"SEPP placements reused: {len(matched) if 'matched' in dir() else 0:,}")
log(f"SEPP placements missing: {len(missing):,}")

summary = {
    "total_candidates": int(len(df)),
    "support_floor_pass": int(n_pass),
    "support_floor_clean": int((passed_df["contamination_flag"] == False).sum()),  # noqa: E712
    "support_floor_flagged": int(passed_df["contamination_flag"].sum()),
    "support_floor_criteria": {"min_bioprojects": 2, "min_reads": 50},
    "sepp_reused": int(len(matched)) if "matched" in dir() else 0,
    "sepp_missing": int(len(missing)),
    "sepp_reuse_pct": round(len(matched)/max(n_pass, 1)*100, 2) if "matched" in dir() else 0,
}
with open(PIPE_V2 / "04_summary.json", "w") as f:
    json.dump(summary, f, indent=2)

log(f"Saved summary: {PIPE_V2 / '04_summary.json'}")
log("Done!")
