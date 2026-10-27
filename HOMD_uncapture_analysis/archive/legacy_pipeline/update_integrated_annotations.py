#!/usr/bin/env python3
"""Update integrated candidate annotations with genus-level decontam + full controls.

Key changes from original 05_integrate_candidate_annotations.py:
1. Genus-level decontam (03_decontam_genus/) replaces ASV-level (03_decontam/)
2. Uses all 1,469 controls (no quality filtering)
3. SILVA genus annotation via sequence matching
4. Literature reference table matching by genus
5. Recomputes contamination_concern with updated evidence
"""

import pandas as pd
import numpy as np
import json, os
from pathlib import Path
from collections import defaultdict

BASE_DIR = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE = BASE_DIR / "pipeline"
SCRIPT_DIR = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/script/05_homd_gg2_closed_reference")

CANDIDATES = PIPE / "02_chimera/candidates_all_nonchimeric.tsv"
GENUS_DECONTAM_A = PIPE / "03_decontam_genus/decontam_A_within_study.tsv"
GENUS_DECONTAM_B = PIPE / "03_decontam_genus/decontam_B_pooled.tsv"
ORAL_VS_BLANK = PIPE / "04_oral_vs_blank/oral_vs_blank_enrichment.tsv"
CONTAM_CATALOG = BASE_DIR / "contaminant_catalog/contaminant_asv_catalog.tsv"
LIT_BLACKLIST = SCRIPT_DIR / "contaminant_reference_genera.tsv"
SCREENED = PIPE / "05_contamination_flag/candidates_screened.tsv"  # from flag_contaminants.py

OUTDIR = PIPE / "05_integrated_candidates"
os.makedirs(OUTDIR, exist_ok=True)

def log(msg): print(f"[update-integrate] {msg}", flush=True)

# ============================================================
# 1. Load base candidates
# ============================================================
log("Loading candidate table...")
cand = pd.read_csv(CANDIDATES, sep="\t", low_memory=False)
log(f"  {len(cand):,} non-chimeric candidates")

# Build sequence index
seq_to_derep = dict(zip(cand["sequence"], cand["derep_id"]))

# ============================================================
# 2. Load genus-level decontam (replaces ASV-level)
# ============================================================
log("\nLoading genus-level decontam...")

da = pd.read_csv(GENUS_DECONTAM_A, sep="\t", low_memory=False)
db = pd.read_csv(GENUS_DECONTAM_B, sep="\t", low_memory=False)

# A: within-study -> n_BPs flagged per genus
agg_a = da.groupby("feature_id")["contaminant"].sum()
genus_n_bps = agg_a.to_dict()  # genus -> n BPs flagged

# B: pooled
# Decontam R output has TRUE/FALSE strings
db["is_contam"] = db["contaminant"].astype(str).str.upper() == "TRUE"
genus_pooled_contam = set(db[db["is_contam"]]["feature_id"].str.lower())
genus_pooled_pval = dict(zip(db["feature_id"].str.lower(), db["p"]))
genus_pooled_ctrl_prev = dict(zip(db["feature_id"].str.lower(), db["ctrl_prevalence"]))
genus_pooled_oral_prev = dict(zip(db["feature_id"].str.lower(), db["oral_prevalence"]))

log(f"  Within-study: {len(genus_n_bps):,} genera flagged")
log(f"  Pooled: {len(genus_pooled_contam):,} genera flagged (p<0.5)")

# ============================================================
# 3. Build SILVA genus annotation for candidates
# ============================================================
log("\nBuilding SILVA genus annotations...")

# Load SILVA taxonomy
with open(f"{BASE_DIR}/control_quality_viz/temp/asv_mapping.json") as f:
    asv_seqs = json.load(f)['asv_id_to_seq']

silva_genus_map = {}
with open("/beegfs/dongbiao/silva_138.2/SILVA_138.2_SSURef_NR99_tax_silva.fasta") as f:
    for line in f:
        if line.startswith(">"):
            parts = line.strip().lstrip(">").split(" ", 1)
            if len(parts) > 1:
                ranks = parts[1].split(";")
                genus = ranks[5].strip() if len(ranks) > 5 else ""
            else: genus = ""
            if genus and genus != "Unknown":
                silva_genus_map[parts[0]] = genus

cols = ['query','target','id','alnlen','mism','gaps','qlo','qhi','tlo','thi','ql','tl','qcov','tcov']
hits = pd.read_csv(f"{BASE_DIR}/control_quality_viz/temp/silva_hits.txt", sep="\t", names=cols)
hits['true_id'] = (hits['alnlen'] - hits['mism'] - hits['gaps']) / hits['ql'] * 100
hits_best = hits[hits['true_id'] >= 80].drop_duplicates(subset='query', keep='first')

seq_to_genus = {}
for _, row in hits_best.iterrows():
    seq = asv_seqs.get(row['query'])
    if seq and row['target'] in silva_genus_map:
        seq_to_genus[seq] = silva_genus_map[row['target']]

cand["silva_genus"] = [seq_to_genus.get(seq, "Unknown") for seq in cand["sequence"]]
n_annotated = (cand["silva_genus"] != "Unknown").sum()
log(f"  SILVA-annotated: {n_annotated:,}/{len(cand):,} ({n_annotated/len(cand)*100:.1f}%)")

# ============================================================
# 4. Add genus-level decontam flags
# ============================================================
log("\nApplying genus-level decontam flags...")

cand["genus_decontam_A_n_bps"] = cand["silva_genus"].str.lower().map(
    lambda g: genus_n_bps.get(g, 0) if isinstance(g, str) else 0
)
cand["genus_decontam_B_contaminant"] = cand["silva_genus"].str.lower().isin(genus_pooled_contam)
cand["genus_decontam_B_pval"] = cand["silva_genus"].str.lower().map(
    lambda g: genus_pooled_pval.get(g, np.nan) if isinstance(g, str) else np.nan
)
cand["genus_decontam_B_ctrl_prev"] = cand["silva_genus"].str.lower().map(
    lambda g: genus_pooled_ctrl_prev.get(g, np.nan) if isinstance(g, str) else np.nan
)
cand["genus_decontam_B_oral_prev"] = cand["silva_genus"].str.lower().map(
    lambda g: genus_pooled_oral_prev.get(g, np.nan) if isinstance(g, str) else np.nan
)

n_a2 = (cand["genus_decontam_A_n_bps"] >= 2).sum()
n_b = cand["genus_decontam_B_contaminant"].sum()
log(f"  Genus decontam A (>=2 BP): {n_a2:,}")
log(f"  Genus decontam B (pooled): {n_b:,}")

# ============================================================
# 5. Load literature reference table
# ============================================================
log("\nLoading literature reference table...")
lit = pd.read_csv(LIT_BLACKLIST, sep="\t", comment="#", low_memory=False)
lit_auto_remove = set(lit[lit["classification"] == "auto_remove"]["taxon"].str.lower())
lit_flag_review = set(lit[lit["classification"] == "flag_review"]["taxon"].str.lower())
lit_oral_core = set(lit[lit["classification"] == "oral_core"]["taxon"].str.lower())
log(f"  auto_remove: {len(lit_auto_remove)}, flag_review: {len(lit_flag_review)}, oral_core: {len(lit_oral_core)}")

def classify_literature(genus):
    if pd.isna(genus) or genus == "Unknown":
        return "unknown"
    g = str(genus).lower().strip()
    if g in lit_oral_core: return "oral_core"
    if g in lit_auto_remove: return "auto_remove"
    if g in lit_flag_review: return "flag_review"
    return "unknown"

cand["literature_classification"] = cand["silva_genus"].map(classify_literature)
cand["literature_oral_core"] = cand["literature_classification"] == "oral_core"
cand["literature_contaminant"] = cand["literature_classification"].isin(["auto_remove", "flag_review"])

n_lit_oral = cand["literature_oral_core"].sum()
n_lit_contam = cand["literature_contaminant"].sum()
log(f"  literature oral_core: {n_lit_oral:,}")
log(f"  literature contaminant: {n_lit_contam:,}")

# ============================================================
# 6. Load oral-vs-blank enrichment
# ============================================================
log("\nLoading oral-vs-blank enrichment...")
ovb = pd.read_csv(ORAL_VS_BLANK, sep="\t", low_memory=False)
ovb_map = ovb.set_index("feature_id").to_dict("index")

for col in ["oral_prevalence", "blank_prevalence", "oral_mean_relab",
            "blank_mean_relab", "oral_over_blank_prev_ratio", "residence_pass"]:
    cand[col] = cand["derep_id"].map(lambda x: ovb_map.get(x, {}).get(col, np.nan))

n_residence = cand["residence_pass"].sum()
log(f"  residence_pass: {n_residence}")

# ============================================================
# 7. Load contaminant ASV catalog
# ============================================================
log("\nLoading contaminant ASV catalog...")
cat = pd.read_csv(CONTAM_CATALOG, sep="\t", low_memory=False)
cat["derep_id"] = cat["sequence"].map(seq_to_derep)
cat_matched = cat[cat["derep_id"].notna()].copy()

cat_agg = cat_matched.groupby("derep_id").agg(
    catalog_n_asvs=("asv_id", "count"),
    catalog_max_confidence=("contamination_confidence",
        lambda x: sorted(x, key=lambda v: {"high":3,"medium":2,"low":1}.get(v,0), reverse=True)[0]),
    catalog_max_n_bps=("n_bioprojects", "max"),
).reset_index()

for col in cat_agg.columns:
    if col != "derep_id":
        cand = cand.merge(cat_agg[["derep_id", col]], on="derep_id", how="left")

cand["catalog_n_asvs"] = cand["catalog_n_asvs"].fillna(0).astype(int)
cand["catalog_contaminant"] = cand["catalog_n_asvs"] > 0
n_catalog = cand["catalog_contaminant"].sum()
log(f"  in contaminant catalog: {n_catalog:,}")

# ============================================================
# 8. Build contamination concern summary
# ============================================================
log("\nBuilding contamination concern summary...")

# Contamination concern = ANY evidence line
cand["contamination_concern"] = (
    (cand["genus_decontam_A_n_bps"] >= 2)        # cross-lab decontam (genus-level)
    | cand["genus_decontam_B_contaminant"]         # pooled decontam
    | cand["catalog_contaminant"]                  # in cross-lab catalog
    | cand["literature_contaminant"]               # literature blacklist
    | (cand["residence_pass"] == False)            # failing residence
)

# Hierarchical concern level
cand["contamination_concern_level"] = "clean"

# Highest concern: literature auto_remove
cand.loc[cand["literature_classification"] == "auto_remove", "contamination_concern_level"] = "literature_auto_remove"

# Cross-lab decontam
cand.loc[
    (cand["contamination_concern_level"] == "clean") & (cand["genus_decontam_A_n_bps"] >= 2),
    "contamination_concern_level"
] = "genus_decontam_cross_lab"

# Pooled decontam significant
cand.loc[
    (cand["contamination_concern_level"] == "clean") & cand["genus_decontam_B_contaminant"],
    "contamination_concern_level"
] = "genus_decontam_pooled"

# Literature flag_review
cand.loc[
    (cand["contamination_concern_level"] == "clean") & (cand["literature_classification"] == "flag_review"),
    "contamination_concern_level"
] = "literature_flag_review"

# In contaminant catalog
cand.loc[
    (cand["contamination_concern_level"] == "clean") & cand["catalog_contaminant"],
    "contamination_concern_level"
] = "in_contaminant_catalog"

# Residence fail
cand.loc[
    (cand["contamination_concern_level"] == "clean")
    & (cand["residence_pass"] == False) & cand["residence_pass"].notna(),
    "contamination_concern_level"
] = "residence_fail"

# Oral_core protection note
cand["contamination_note"] = ""
cand.loc[
    (cand["contamination_concern_level"] != "clean") & cand["literature_oral_core"],
    "contamination_note"
] = "FLAGGED_BUT_ORAL_CORE_PROTECTED"

log("Contamination concern breakdown:")
for level in cand["contamination_concern_level"].value_counts().index:
    n = (cand["contamination_concern_level"] == level).sum()
    log(f"  {level}: {n:,} ({n/len(cand)*100:.1f}%)")

# ============================================================
# 9. Save
# ============================================================
log("\nSaving integrated annotations...")

# Select and order columns
output_cols = [
    # Base
    "derep_id", "sequence", "length",
    "n_asvs", "n_bioprojects", "total_reads", "prevalence_samples",
    "subject_prevalence", "support_tier",
    "amplicon_regions", "sequencing_regions", "continents", "sites",
    # HOMD
    "has_homd_hit", "homd_target", "identity", "query_cov",
    "homd_uncaptured", "uncaptured_reason",
    # Chimera
    "is_chimera",
    # SILVA
    "silva_genus",
    # Genus-level decontam
    "genus_decontam_A_n_bps", "genus_decontam_B_contaminant",
    "genus_decontam_B_pval", "genus_decontam_B_ctrl_prev", "genus_decontam_B_oral_prev",
    # Oral-vs-blank residence
    "oral_prevalence", "blank_prevalence", "oral_mean_relab",
    "blank_mean_relab", "oral_over_blank_prev_ratio", "residence_pass",
    # Contaminant catalog
    "catalog_n_asvs", "catalog_contaminant", "catalog_max_confidence", "catalog_max_n_bps",
    # Literature
    "literature_classification", "literature_oral_core", "literature_contaminant",
    # Summary
    "contamination_concern", "contamination_concern_level", "contamination_note",
]

for col in output_cols:
    if col not in cand.columns:
        cand[col] = np.nan

out_df = cand[output_cols].copy()

# Boolean -> string
for col in ["has_homd_hit", "homd_uncaptured", "is_chimera", "residence_pass",
            "catalog_contaminant", "literature_oral_core", "literature_contaminant",
            "contamination_concern", "genus_decontam_B_contaminant"]:
    if col in out_df.columns:
        out_df[col] = out_df[col].astype(str).str.lower()

out_path = OUTDIR / "integrated_candidate_annotations.tsv"
out_df.to_csv(out_path, sep="\t", index=False)
log(f"Saved: {out_path} ({os.path.getsize(out_path):,} bytes)")
log(f"  {len(out_df):,} rows × {len(out_df.columns)} columns")

# ============================================================
# 10. Summary
# ============================================================
log(f"\n{'='*60}")
log(f"UPDATED INTEGRATED TABLE SUMMARY")
log(f"{'='*60}")
log(f"Total candidates: {len(out_df):,}")
log(f"Support tiers:")
for tier in ["high", "medium", "low_rare_tail"]:
    n = (out_df["support_tier"] == tier).sum()
    log(f"  {tier}: {n:,}")

log(f"\nContamination evidence:")
log(f"  genus decontam A (>=2 BP):  {(out_df['genus_decontam_A_n_bps'] >= 2).sum():,}")
log(f"  genus decontam B (pooled):   {(out_df['genus_decontam_B_contaminant'] == 'true').sum():,}")
log(f"  contaminant catalog:         {(out_df['catalog_contaminant'] == 'true').sum():,}")
log(f"  literature auto_remove:      {(out_df['literature_classification'] == 'auto_remove').sum():,}")
log(f"  literature flag_review:      {(out_df['literature_classification'] == 'flag_review').sum():,}")
log(f"  literature oral_core:        {(out_df['literature_classification'] == 'oral_core').sum():,}")
log(f"  residence fail:              {(out_df['residence_pass'] == 'false').sum():,}")
log(f"  ANY concern:                 {(out_df['contamination_concern'] == 'true').sum():,}")
log(f"  FLAGGED BUT ORAL_CORE:       {(out_df['contamination_note'] == 'FLAGGED_BUT_ORAL_CORE_PROTECTED').sum():,}")

# Comparison with old version
log(f"\nUpdated key metrics vs old report:")
log(f"  Candidates: 54,254 (old: 54,290)")
n_clean = (out_df["contamination_concern_level"] == "clean").sum()
log(f"  Clean: {n_clean:,} (old: 47,203)")
log(f"  Contamination-flagged: {len(out_df) - n_clean:,} (old: 7,087)")

log(f"\nDone!")
