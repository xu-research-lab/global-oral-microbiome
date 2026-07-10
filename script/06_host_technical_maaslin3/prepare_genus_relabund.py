#!/usr/bin/env python3
"""
Prepare genus-level relative abundance table from BIOM file for MaAsLin3 v2.
Uses new metadata with native column names (no renaming).
Filters genera to those that can pass prevalence filters in at least one analysis group.
"""
from pathlib import Path
import biom
import numpy as np
import pandas as pd
from scipy import sparse

BASE = Path("/beegfs/huangxiaochang/global_oral/new_analysis")
BIOM_PATH = BASE / "biom_oral_qc_dedup.biom"
META_PATH = BASE / "metadata_oral_qc_dedup.tsv"
TAXONOMY_PATH = BASE.parent / "taxonomy.tsv"
OUT_DIR = BASE / "maaslin3_v2"
OUT_PATH = OUT_DIR / "genus_table_relabund.tsv"

MIN_PREVALENCE = 0.10
MIN_GROUP_N = 50

# Preferred diseases for disease analysis (from old code)
PREFERRED_DISEASES = [
    "Dental Caries",
    "Head and Neck Cancer",
    "HIV",
    "Oral Cancer",
    "OSCC",
    "Periodontal Disease",
    "RA",
    "Benign Tumor",
    "Breast Cancer",
    "CRC",
    "IBD",
    "Lung Cancer",
]

# Site name mapping (same as old code)
SITE_MAP = {
    "saliva": "Saliva",
    "tongue": "Tongue",
    "tooth": "Tooth",
    "oral swab": "Oral_swab",
    "oral mucosa": "Oral_Mucosa",
    "gingiva": "Gingiva",
    "other": "Others",
}

# Status mapping
STATUS_MAP = {
    "HC": "Healthy Control",
    "Healthy": "Healthy Control",
}


def extract_genus(taxon: str, fallback: str) -> str:
    """Extract genus name from GTDB-style taxonomy string."""
    for part in str(taxon).split(";"):
        part = part.strip()
        if part.startswith("g__"):
            genus = part[3:].strip()
            return genus if genus else fallback
    return fallback


def mark_prevalent(binary_mat, sample_idx, keep):
    """Mark genera with prevalence >= MIN_PREVALENCE in given sample subset."""
    if len(sample_idx) == 0:
        return
    prev = np.asarray(binary_mat[:, sample_idx].getnnz(axis=1)).ravel() / float(len(sample_idx))
    keep |= prev >= MIN_PREVALENCE


def load_and_prepare_metadata():
    """Load metadata with native column names and apply mappings."""
    meta = pd.read_csv(META_PATH, sep="\t", low_memory=False)

    # Clean string columns
    str_cols = ["sample_id", "bioproject_id", "site_cleaned", "continent_group",
                "status_cleaned", "sex_cleaned", "age_cleaned_group",
                "sequencing_region_cleaned", "disease_study_cleaned"]
    for col in str_cols:
        if col in meta.columns:
            meta[col] = meta[col].astype(str).str.strip()
            meta.loc[meta[col].isin(["", "nan", "None", "NA"]), col] = pd.NA

    # Apply site mapping
    if "site_cleaned" in meta.columns:
        meta["site_cleaned"] = meta["site_cleaned"].replace(SITE_MAP)

    # Apply status mapping
    if "status_cleaned" in meta.columns:
        meta["status_cleaned"] = meta["status_cleaned"].replace(STATUS_MAP)

    # Replace full-width ＞ with >
    if "age_cleaned_group" in meta.columns:
        meta["age_cleaned_group"] = meta["age_cleaned_group"].str.replace("＞", ">", regex=False)

    return meta


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load taxonomy
    print(f"Loading taxonomy: {TAXONOMY_PATH}")
    tax = pd.read_csv(TAXONOMY_PATH, sep="\t")
    tax = tax.rename(columns={"Feature ID": "FeatureID", "Taxon": "Taxon"})
    tax_map = dict(zip(tax["FeatureID"].astype(str), tax["Taxon"].astype(str)))

    # Load BIOM
    print(f"Loading BIOM: {BIOM_PATH}")
    table = biom.load_table(str(BIOM_PATH))
    obs_ids = np.array(table.ids(axis="observation"), dtype=str)
    sample_ids = np.array(table.ids(axis="sample"), dtype=str)
    mat = table.matrix_data.tocsr()

    # Aggregate to genus level
    genera = np.array([extract_genus(tax_map.get(obs, ""), obs) for obs in obs_ids], dtype=object)
    unique_genera, inverse = np.unique(genera, return_inverse=True)
    group = sparse.csr_matrix(
        (np.ones(len(obs_ids), dtype=np.float64), (inverse, np.arange(len(obs_ids)))),
        shape=(len(unique_genera), len(obs_ids)),
    )
    genus_counts = (group @ mat).tocsr()
    binary = genus_counts.copy()
    binary.data = np.ones_like(binary.data)

    # Load and prepare metadata
    meta = load_and_prepare_metadata()
    meta = meta.dropna(subset=["sample_id"]).copy()
    meta["sample_id"] = meta["sample_id"].astype(str)

    # Map sample IDs to matrix column indices
    sample_to_idx = {sid: i for i, sid in enumerate(sample_ids)}
    meta = meta[meta["sample_id"].isin(sample_to_idx)].copy()
    meta["_idx"] = meta["sample_id"].map(sample_to_idx).astype(int)

    print(f"Samples with metadata and abundance: {len(meta)}")

    keep = np.zeros(len(unique_genera), dtype=bool)

    # 1. Preserve genera prevalent in any Site group (for site-stratified analysis)
    for site, sub in meta.groupby("site_cleaned", dropna=True):
        if len(sub) >= MIN_GROUP_N:
            print(f"  Site '{site}': {len(sub)} samples")
            mark_prevalent(binary, sub["_idx"].to_numpy(), keep)

    # 2. Preserve genera prevalent in any disease analysis group
    disease_meta = meta[
        meta["status_cleaned"].notna()
        & ~meta["status_cleaned"].isin(["Other", "Healthy Control"])
        & meta["bioproject_id"].notna()
    ].copy()
    for disease in PREFERRED_DISEASES:
        if disease not in set(disease_meta["status_cleaned"]):
            continue
        disease_bp = disease_meta.loc[
            disease_meta["status_cleaned"] == disease, "bioproject_id"
        ].unique()
        sub = disease_meta[
            disease_meta["status_cleaned"].isin(["Healthy Control", disease])
            & disease_meta["bioproject_id"].isin(disease_bp)
        ].copy()
        if sub.empty:
            continue
        bp_tbl = pd.crosstab(sub["bioproject_id"], sub["status_cleaned"])
        if disease not in bp_tbl.columns or "Healthy Control" not in bp_tbl.columns:
            continue
        valid_bp = bp_tbl.index[(bp_tbl[disease] > 0) & (bp_tbl["Healthy Control"] > 0)]
        sub = sub[sub["bioproject_id"].isin(valid_bp)]
        if (sub["status_cleaned"].eq(disease).sum() >= 30
                and sub["status_cleaned"].eq("Healthy Control").sum() >= 30
                and sub["bioproject_id"].nunique() >= 2):
            print(f"  Disease '{disease}': {len(sub)} samples")
            mark_prevalent(binary, sub["_idx"].to_numpy(), keep)

    # 3. Unified region model: per-site prevalence across all samples
    for site, sub in meta.groupby("site_cleaned", dropna=True):
        if len(sub) >= MIN_GROUP_N:
            enough_regions = (
                sub["sequencing_region_cleaned"].dropna().nunique() >= 2
            )
            if enough_regions:
                mark_prevalent(binary, sub["_idx"].to_numpy(), keep)

    kept = np.where(keep)[0]
    print(f"\nKeeping {len(kept)} / {len(unique_genera)} genera that can pass prevalence filters")

    genus_counts = genus_counts[kept, :]
    kept_genera = unique_genera[kept]

    # Convert to relative abundance
    sample_totals = np.asarray(genus_counts.sum(axis=0)).ravel()
    sample_totals[sample_totals == 0] = 1.0
    genus_rel = (genus_counts @ sparse.diags(1.0 / sample_totals)).toarray()

    # Write output
    print(f"Writing genus relative abundance: {OUT_PATH}")
    with OUT_PATH.open("w") as fh:
        fh.write("Genus\t" + "\t".join(sample_ids) + "\n")
        for genus, row in zip(kept_genera, genus_rel):
            vals = "\t".join(f"{x:.8g}" for x in row)
            fh.write(f"{genus}\t{vals}\n")

    print(f"Done: {len(kept_genera)} genera x {len(sample_ids)} samples")


if __name__ == "__main__":
    main()
