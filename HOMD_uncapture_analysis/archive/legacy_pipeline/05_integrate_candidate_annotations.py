#!/usr/bin/env python3
"""
Build integrated ASV annotation table merging all evidence lines.

Inputs (all from the pipeline):
  1. pipeline/02_chimera/candidates_all_nonchimeric.tsv  (base: HOMD + chimera)
  2. pipeline/03_decontam/decontam_per_bioproject.tsv     (decontam per-study flags)
  3. pipeline/03_decontam/decontam_contaminant_catalog.tsv (cross-lab contaminants)
  4. pipeline/04_oral_vs_blank/oral_vs_blank_enrichment.tsv (residence test)
  5. contaminant_catalog/contaminant_asv_catalog.tsv       (control-derived catalog)
  6. script/.../contaminant_reference_genera.tsv            (literature blacklist)

Output:
  pipeline/05_integrated_candidates/integrated_candidate_annotations.tsv

Gate columns (each is a FLAG, not an auto-delete):
  - homd_uncaptured           : HOMD identity < 97% or qcov < 99% or no hit
  - is_chimera                : flagged by vsearch --uchime_ref
  - decontam_flagged          : flagged by decontam in >=1 study
  - decontam_cross_lab        : flagged in >=2 independent studies (high confidence)
  - catalog_contaminant       : sequence matches contaminant ASV catalog
  - catalog_contaminant_genus : genus is known contaminant per catalog
  - literature_contaminant    : genus in literature blacklist (flag_review or auto_remove)
  - literature_oral_core      : genus is protected oral core (never remove)
  - residence_pass            : oral >> blank enrichment (from step 04)
  - contamination_concern     : ANY contamination flag is raised (summary column)
"""

import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path
from collections import defaultdict

# ============================================================
# Config
# ============================================================
BASE_DIR = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPELINE_DIR = BASE_DIR / "pipeline"
SCRIPT_DIR = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/script/05_homd_gg2_closed_reference")

CANDIDATES = PIPELINE_DIR / "02_chimera/candidates_all_nonchimeric.tsv"
DECONTAM_PER_BP = PIPELINE_DIR / "03_decontam/decontam_per_bioproject.tsv"
DECONTAM_CATALOG = PIPELINE_DIR / "03_decontam/decontam_contaminant_catalog.tsv"
ORAL_VS_BLANK = PIPELINE_DIR / "04_oral_vs_blank/oral_vs_blank_enrichment.tsv"
CONTAM_CATALOG = BASE_DIR / "contaminant_catalog/contaminant_asv_catalog.tsv"
LIT_BLACKLIST = SCRIPT_DIR / "contaminant_reference_genera.tsv"

OUTDIR = PIPELINE_DIR / "05_integrated_candidates"
os.makedirs(OUTDIR, exist_ok=True)


def log(msg: str) -> None:
    print(f"[integrate] {msg}", flush=True)


# ============================================================
# 1. Load base candidate table
# ============================================================
log("Loading candidate table...")
cand = pd.read_csv(CANDIDATES, sep="\t", low_memory=False)
log(f"  {len(cand):,} non-chimeric candidates")

# Build sequence -> derep_id index
seq_to_derep = dict(zip(cand["sequence"], cand["derep_id"]))
log(f"  sequence index: {len(seq_to_derep):,} unique sequences")

# ============================================================
# 2. Load decontam results (per-study + cross-lab catalog)
# ============================================================
log("Loading decontam results...")

decontam_bp = pd.read_csv(DECONTAM_PER_BP, sep="\t", low_memory=False)
decontam_flagged_per_bp = defaultdict(set)  # derep_id -> set of bioprojects
decontam_flagged_count = defaultdict(int)   # derep_id -> n studies flagged
for _, row in decontam_bp.iterrows():
    if row["contaminant"] and not (isinstance(row["contaminant"], float) and np.isnan(row["contaminant"])):
        decontam_flagged_per_bp[row["feature_id"]].add(row["bioproject"])
        decontam_flagged_count[row["feature_id"]] += 1

decontam_catalog = pd.read_csv(DECONTAM_CATALOG, sep="\t", low_memory=False)
cross_lab_set = set(decontam_catalog["feature_id"].tolist())

n_flagged = len(decontam_flagged_count)
log(f"  decontam-flagged features: {n_flagged} (cross-lab: {len(cross_lab_set)})")

# Add to candidates
cand["decontam_n_studies_flagged"] = cand["derep_id"].map(
    lambda x: decontam_flagged_count.get(x, 0)
)
cand["decontam_cross_lab"] = cand["derep_id"].isin(cross_lab_set)
cand["decontam_flagged_bps"] = cand["derep_id"].map(
    lambda x: ";".join(sorted(decontam_flagged_per_bp.get(x, set()))) or ""
)

# ============================================================
# 3. Load oral-vs-blank enrichment
# ============================================================
log("Loading oral-vs-blank enrichment...")
ovb = pd.read_csv(ORAL_VS_BLANK, sep="\t", low_memory=False)
ovb_map = ovb.set_index("feature_id").to_dict("index")

for col in ["oral_prevalence", "blank_prevalence", "oral_mean_relab",
            "blank_mean_relab", "oral_over_blank_prev_ratio", "residence_pass"]:
    cand[col] = cand["derep_id"].map(lambda x: ovb_map.get(x, {}).get(col, np.nan))

n_residence = cand["residence_pass"].sum()
log(f"  residence_pass: {n_residence} features")

# ============================================================
# 4. Load contaminant ASV catalog (from script 09, built from control DADA2)
# ============================================================
log("Loading contaminant ASV catalog...")
cat = pd.read_csv(CONTAM_CATALOG, sep="\t", low_memory=False)
# Match by sequence -> derep_id
cat["derep_id"] = cat["sequence"].map(seq_to_derep)
cat_matched = cat[cat["derep_id"].notna()].copy()

# Build contamination per derep_id (aggregate multiple catalog ASVs)
cat_agg = cat_matched.groupby("derep_id").agg(
    catalog_n_asvs=("asv_id", "count"),
    catalog_max_confidence=("contamination_confidence", lambda x: sorted(x, key=lambda v: {"high": 3, "medium": 2, "low": 1}.get(v, 0), reverse=True)[0]),
    catalog_contaminant_genera=("genus", lambda x: ";".join(sorted(set(x)))),
    catalog_max_n_bps=("n_bioprojects", "max"),
    catalog_max_n_samples=("n_control_samples", "max"),
).reset_index()

for col in cat_agg.columns:
    if col == "derep_id":
        continue
    cand = cand.merge(cat_agg[["derep_id", col]], on="derep_id", how="left")

cand["catalog_n_asvs"] = cand["catalog_n_asvs"].fillna(0).astype(int)
cand["catalog_contaminant"] = cand["catalog_n_asvs"] > 0

n_catalog = cand["catalog_contaminant"].sum()
log(f"  derep features in contaminant catalog: {n_catalog}")

# ============================================================
# 5. Load literature-based contaminant genus blacklist
# ============================================================
log("Loading literature blacklist...")
lit = pd.read_csv(LIT_BLACKLIST, sep="\t", comment="#", low_memory=False)
lit_auto_remove = set(lit[lit["classification"] == "auto_remove"]["taxon"].str.lower())
lit_flag_review = set(lit[lit["classification"] == "flag_review"]["taxon"].str.lower())
lit_oral_core = set(lit[lit["classification"] == "oral_core"]["taxon"].str.lower())

log(f"  auto_remove genera: {len(lit_auto_remove)}")
log(f"  flag_review genera: {len(lit_flag_review)}")
log(f"  oral_core genera: {len(lit_oral_core)}")

# For genus-level matching, we need taxonomy. Use the catalog's genus info.
# Build derep_id -> best genus from catalog
derep_genus_map = {}
for _, row in cat_matched.iterrows():
    did = row["derep_id"]
    genus = str(row.get("genus", "")).strip()
    if genus and genus not in {"Unclassified", "Unknown", ""}:
        if did not in derep_genus_map:
            derep_genus_map[did] = genus

# Add genus column from catalog (where available)
cand["genus_from_catalog"] = cand["derep_id"].map(derep_genus_map)

# Check against literature blacklist
def literature_flag(genus):
    """Return literature classification for a genus."""
    if pd.isna(genus) or not genus:
        return "unknown"
    g = str(genus).lower().strip()
    if g in lit_oral_core:
        return "oral_core"
    if g in lit_auto_remove:
        return "auto_remove"
    if g in lit_flag_review:
        return "flag_review"
    return "unknown"

cand["literature_classification"] = cand["genus_from_catalog"].map(literature_flag)
cand["literature_oral_core"] = cand["literature_classification"] == "oral_core"
cand["literature_contaminant"] = cand["literature_classification"].isin(["auto_remove", "flag_review"])

n_lit_contam = cand["literature_contaminant"].sum()
n_lit_oral = cand["literature_oral_core"].sum()
log(f"  literature contaminant genus: {n_lit_contam}")
log(f"  literature oral_core genus: {n_lit_oral}")

# ============================================================
# 6. Build contamination concern summary
# ============================================================
log("Building contamination concern summary...")

# Contamination concern = ANY flag is raised
cand["contamination_concern"] = (
    (cand["decontam_n_studies_flagged"] > 0)
    | cand["catalog_contaminant"]
    | cand["literature_contaminant"]
    | (cand["residence_pass"] == False)  # failing residence is also a concern
)

# But oral_core is protected - even if flagged, note it separately
cand["contamination_concern_level"] = "clean"
cand.loc[cand["decontam_cross_lab"], "contamination_concern_level"] = "cross_lab_contaminant"
cand.loc[
    (cand["contamination_concern_level"] == "clean")
    & (cand["decontam_n_studies_flagged"] == 1),
    "contamination_concern_level"
] = "decontam_single_study"
cand.loc[
    (cand["contamination_concern_level"] == "clean")
    & cand["catalog_contaminant"],
    "contamination_concern_level"
] = "in_contaminant_catalog"
cand.loc[
    (cand["contamination_concern_level"] == "clean")
    & cand["literature_contaminant"],
    "contamination_concern_level"
] = "literature_contaminant_genus"
# Residence fail
cand.loc[
    (cand["contamination_concern_level"] == "clean")
    & (cand["residence_pass"] == False)
    & cand["residence_pass"].notna(),
    "contamination_concern_level"
] = "residence_fail"

# Oral_core override: if flagged but is oral_core, note it
cand["contamination_note"] = ""
cand.loc[
    (cand["contamination_concern_level"] != "clean")
    & cand["literature_oral_core"],
    "contamination_note"
] = "FLAGGED_BUT_ORAL_CORE_PROTECTED"

log(f"Contamination concern breakdown:")
for level in ["clean", "cross_lab_contaminant", "decontam_single_study",
              "in_contaminant_catalog", "literature_contaminant_genus", "residence_fail"]:
    n = (cand["contamination_concern_level"] == level).sum()
    log(f"  {level}: {n:,}")

# ============================================================
# 7. Final column selection & ordering
# ============================================================
log("Building final output...")

# Define output columns in logical groups
base_cols = [
    "derep_id", "sequence", "length",
    "n_asvs", "n_bioprojects", "total_reads", "prevalence_samples",
    "subject_prevalence", "support_tier",
    "amplicon_regions", "sequencing_regions", "continents", "sites",
]

homd_cols = [
    "has_homd_hit", "homd_target", "identity", "query_cov",
    "homd_uncaptured", "uncaptured_reason",
]

chimera_cols = [
    "is_chimera",
]

decontam_cols = [
    "decontam_n_studies_flagged", "decontam_cross_lab", "decontam_flagged_bps",
]

residence_cols = [
    "oral_prevalence", "blank_prevalence", "oral_mean_relab",
    "blank_mean_relab", "oral_over_blank_prev_ratio", "residence_pass",
]

catalog_cols = [
    "catalog_n_asvs", "catalog_contaminant", "catalog_max_confidence",
    "catalog_contaminant_genera", "catalog_max_n_bps", "catalog_max_n_samples",
    "genus_from_catalog",
]

literature_cols = [
    "literature_classification", "literature_oral_core", "literature_contaminant",
]

summary_cols = [
    "contamination_concern", "contamination_concern_level", "contamination_note",
]

# Extra metadata columns from original
extra_cols = [
    "representative_asv_id", "bioprojects", "platforms", "countries",
    "disease_studies", "statuses", "target_cov", "alnlen",
]

all_output_cols = (
    base_cols + homd_cols + chimera_cols + decontam_cols + residence_cols
    + catalog_cols + literature_cols + summary_cols + extra_cols
)

# Ensure all columns exist, fill missing
for col in all_output_cols:
    if col not in cand.columns:
        cand[col] = np.nan

out_df = cand[all_output_cols].copy()

# Convert boolean columns for readability
bool_cols = [
    "has_homd_hit", "homd_uncaptured", "is_chimera",
    "decontam_cross_lab", "residence_pass",
    "catalog_contaminant", "literature_oral_core", "literature_contaminant",
    "contamination_concern",
]
for col in bool_cols:
    if col in out_df.columns:
        out_df[col] = out_df[col].astype(str).str.lower()

# ============================================================
# 8. Save
# ============================================================
out_path = OUTDIR / "integrated_candidate_annotations.tsv"
out_df.to_csv(out_path, sep="\t", index=False)
log(f"\nSaved: {out_path} ({os.path.getsize(out_path):,} bytes)")
log(f"  {len(out_df):,} rows × {len(out_df.columns)} columns")

# ============================================================
# 9. Summary stats
# ============================================================
log(f"\n=== INTEGRATED TABLE SUMMARY ===")
log(f"Total candidates (non-chimeric HOMD-uncaptured): {len(out_df):,}")
log(f"Support tiers:")
for tier in ["high", "medium", "low_rare_tail"]:
    n = (out_df["support_tier"] == tier).sum()
    log(f"  {tier}: {n:,} ({n/len(out_df)*100:.1f}%)")

log(f"\nContamination flags:")
log(f"  decontam cross-lab (≥2 studies):    {(out_df['decontam_cross_lab'] == 'true').sum():,}")
log(f"  decontam single-study flagged:       {(out_df['decontam_n_studies_flagged'] == 1).sum():,}")
log(f"  in contaminant catalog:              {(out_df['catalog_contaminant'] == 'true').sum():,}")
log(f"  literature contaminant genus:         {(out_df['literature_contaminant'] == 'true').sum():,}")
log(f"  residence_fail:                       {(out_df['residence_pass'] == 'false').sum():,}")
log(f"  ANY contamination concern:            {(out_df['contamination_concern'] == 'true').sum():,}")
log(f"  literature oral_core (protected):     {(out_df['literature_oral_core'] == 'true').sum():,}")
log(f"  FLAGGED BUT ORAL_CORE:                {(out_df['contamination_note'] == 'FLAGGED_BUT_ORAL_CORE_PROTECTED').sum():,}")

# Per-region summary
if "amplicon_regions" in out_df.columns:
    log(f"\nPer-region breakdown:")
    for reg in ["V3-V4", "V4", "V1-V3", "V4-V5", "V6-V8"]:
        sub = out_df[out_df["amplicon_regions"] == reg]
        if len(sub) > 0:
            n_clean = (sub["contamination_concern_level"] == "clean").sum()
            n_decontam = (sub["decontam_cross_lab"] == "true").sum()
            log(f"  {reg}: {len(sub):,} candidates, clean={n_clean:,}, cross_lab_contam={n_decontam}")

log(f"\nOutput: {OUTDIR}/")
for f in sorted(OUTDIR.iterdir()):
    if f.is_file():
        sz = os.path.getsize(f)
        log(f"  {f.name} ({sz:,} bytes)")

log("Done!")
