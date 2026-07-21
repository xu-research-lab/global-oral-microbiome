#!/usr/bin/env python3
"""
00_extract_red_complex.py
=========================
Species-level extraction of the red complex (P. gingivalis, T. denticola,
T. forsythia) from the global oral microbiome dataset using GTDB taxonomy.

Key corrections over the previous script (new_analysis/red_complex/):
  1. Species-level extraction using GTDB taxonomy (not genus-level)
  2. Correct disease status column (status_cleaned)
  3. Proper column mapping for site, geography, etc.
  4. Detection flags and co-occurrence combinations

Input:
  - biom_oral_qc.biom (ASV-level count table)
  - taxonomy.tsv (GTDB species annotation)
  - metadata_oral_qc.tsv (sample metadata)

Output:
  - red_complex_sample_abundance.tsv : merged abundance + metadata
"""

import biom
import numpy as np
import pandas as pd
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

# Paths
BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome")
DATA_DIR = BASE / "data" / "metadata-and-biom-v2-qc"
TAXONOMY_PATH = Path("/beegfs/huangxiaochang/global_oral/taxonomy.tsv")
OUTDIR = BASE / "red_complex"
OUTDIR.mkdir(parents=True, exist_ok=True)

# ── Red Complex Species Mapping (GTDB taxonomy) ─────────────────────
# These are the strict species-level markers identified from taxonomy.tsv
RED_COMPLEX_SPECIES = {
    "P_gingivalis": {
        "label": "Porphyromonas gingivalis",
        "gtdb_marker": "s__Porphyromonas_A gingivalis",
        "genus_gtb": "g__Porphyromonas_A",
    },
    "T_denticola": {
        "label": "Treponema denticola",
        "gtdb_marker": "s__Treponema_B_986141 denticola",
        "genus_gtb": "g__Treponema_B_986141",
    },
    "T_forsythia": {
        "label": "Tannerella forsythia",
        "gtdb_marker": "s__Tannerella forsythia",
        "genus_gtb": "g__Tannerella",
    },
}

# ── Column name mapping (consistent with the project convention) ───
COLUMN_RENAME = {
    "sample_id": "SampleID",
    "bioproject_id": "BioProject",
    "subject_id": "SubjectID",
    "site_cleaned": "Site",
    "continent_group": "Geography",
    "status_cleaned": "DiseaseStatus",
    "sex_cleaned": "Sex",
    "age_cleaned_group": "AgeGroup",
    "age_cleaned": "Age",
    "bmi_cleaned": "BMI",
    "platform_cleaned": "Platform",
    "sequencing_region_cleaned": "SeqRegion",
    "is_time_series": "IsTimeSeries",
    "country_cleaned": "Country",
    "Observed": "Alpha_Observed",
    "Shannon": "Alpha_Shannon",
    "Simpson": "Alpha_Simpson",
    "InvSimpson": "Alpha_InvSimpson",
    "Chao1": "Alpha_Chao1",
    "chim_perc": "ChimeraPct",
    "asv_count": "ASVCount",
    "total_reads": "TotalReads",
}


def load_taxonomy():
    """Load GTDB taxonomy, return dict: FeatureID -> Taxon string."""
    tax = pd.read_csv(TAXONOMY_PATH, sep="\t")
    tax = tax.rename(columns={"Feature ID": "FeatureID", "Taxon": "Taxon"})
    tax["FeatureID"] = tax["FeatureID"].astype(str)
    return tax.set_index("FeatureID")


def find_rc_asvs(obs_ids, tax_df):
    """Find ASV IDs belonging to each red complex species using GTDB markers.

    Returns dict: species_key -> list of FeatureIDs
    """
    rc_asvs = {}
    for key, info in RED_COMPLEX_SPECIES.items():
        marker = info["gtdb_marker"]
        matched = []
        for obs in obs_ids:
            obs_str = str(obs)
            if obs_str in tax_df.index:
                tax_str = str(tax_df.at[obs_str, "Taxon"])
                if marker in tax_str:
                    matched.append(obs_str)
        rc_asvs[key] = matched
        print(f"  {info['label']}: {len(matched)} ASVs matched (marker: '{marker}')")
    return rc_asvs


def extract_rc_abundance(table, rc_asvs, sample_ids):
    """Extract relative abundance for each red complex species per sample.

    Returns DataFrame with columns: SampleID, P_gingivalis, T_denticola,
    T_forsythia, RedComplexTotal + detection flags.
    """
    sample_totals = np.asarray(table.sum(axis="sample")).ravel()

    abund = {}
    for key, asv_list in rc_asvs.items():
        if not asv_list:
            abund[key] = np.zeros(len(sample_ids))
            continue
        sub = table.filter(asv_list, axis="observation", inplace=False)
        counts = np.asarray(sub.matrix_data.sum(axis=0)).ravel()
        abund[key] = counts / sample_totals

    df = pd.DataFrame({
        "SampleID": sample_ids,
        "P_gingivalis": abund.get("P_gingivalis", np.zeros(len(sample_ids))),
        "T_denticola": abund.get("T_denticola", np.zeros(len(sample_ids))),
        "T_forsythia": abund.get("T_forsythia", np.zeros(len(sample_ids))),
    })

    # Total red complex abundance
    df["RedComplexTotal"] = df[["P_gingivalis", "T_denticola", "T_forsythia"]].sum(axis=1)

    # Detection flags (binary)
    for sp in ["P_gingivalis", "T_denticola", "T_forsythia"]:
        df[f"{sp}_detected"] = (df[sp] > 0).astype(int)
    df["RedComplex_detected"] = (df["RedComplexTotal"] > 0).astype(int)

    # Count how many of the three are present
    df["RC_member_count"] = (
        df["P_gingivalis_detected"]
        + df["T_denticola_detected"]
        + df["T_forsythia_detected"]
    )

    # All three present flag
    df["AllThree_detected"] = (df["RC_member_count"] == 3).astype(int)

    # 8-combination classification
    def combo_label(row):
        pg = "+" if row["P_gingivalis_detected"] else "-"
        td = "+" if row["T_denticola_detected"] else "-"
        tf = "+" if row["T_forsythia_detected"] else "-"
        return f"Pg{pg}_Td{td}_Tf{tf}"

    df["RC_combination"] = df.apply(combo_label, axis=1)

    return df


def merge_with_metadata(rc_df, meta_path):
    """Merge red complex abundance with sample metadata."""
    meta = pd.read_csv(meta_path, sep="\t", low_memory=False)

    # Rename columns for consistency
    meta = meta.rename(columns={
        k: v for k, v in COLUMN_RENAME.items() if k in meta.columns
    })

    # Ensure string types for merge keys
    meta["SampleID"] = meta["SampleID"].astype(str).str.strip()
    rc_df["SampleID"] = rc_df["SampleID"].astype(str).str.strip()

    # Merge
    merged = meta.merge(rc_df, on="SampleID", how="inner")

    # Fill NA for samples not in RC data (shouldn't happen, but safe)
    for col in ["P_gingivalis", "T_denticola", "T_forsythia", "RedComplexTotal"]:
        if col in merged.columns:
            merged[col] = merged[col].fillna(0)

    print(f"  Metadata samples: {len(meta)}")
    print(f"  BIOM samples: {len(rc_df)}")
    print(f"  Merged samples: {len(merged)}")

    return merged


def print_summary(merged):
    """Print summary statistics."""
    print("\n" + "=" * 60)
    print("RED COMPLEX EXTRACTION SUMMARY")
    print("=" * 60)

    print(f"\nTotal samples: {len(merged)}")
    print(f"BioProjects: {merged['BioProject'].nunique()}")
    print(f"Subjects: {merged['SubjectID'].nunique()}")

    print(f"\n── Detection Rates ──")
    for sp in ["P_gingivalis", "T_denticola", "T_forsythia"]:
        if f"{sp}_detected" in merged.columns:
            rate = merged[f"{sp}_detected"].mean()
            print(f"  {sp}: {rate:.3f} ({rate*100:.1f}%)")
    print(f"  Any red complex: {merged['RedComplex_detected'].mean():.3f}")
    print(f"  All three: {merged['AllThree_detected'].mean():.3f}")

    print(f"\n── Co-occurrence Combinations ──")
    combo_order = [
        "Pg+_Td+_Tf+", "Pg+_Td+_Tf-", "Pg+_Td-_Tf+", "Pg-_Td+_Tf+",
        "Pg+_Td-_Tf-", "Pg-_Td+_Tf-", "Pg-_Td-_Tf+", "Pg-_Td-_Tf-",
    ]
    for combo in combo_order:
        if combo in merged["RC_combination"].values:
            n = (merged["RC_combination"] == combo).sum()
            pct = 100 * n / len(merged)
            print(f"  {combo}: {n:6d} ({pct:.1f}%)")

    print(f"\n── Abundance (non-zero only) ──")
    for sp in ["P_gingivalis", "T_denticola", "T_forsythia", "RedComplexTotal"]:
        nonzero = merged[sp][merged[sp] > 0]
        if len(nonzero) > 0:
            print(f"  {sp}: mean={nonzero.mean():.6f}, median={nonzero.median():.6f}, "
                  f"max={nonzero.max():.6f}")


def main():
    print("Loading BIOM table...")
    table = biom.load_table(str(DATA_DIR / "biom_oral_qc.biom"))
    obs_ids = np.array(table.ids(axis="observation"))
    sample_ids = np.array(table.ids(axis="sample"))
    print(f"  Observations: {len(obs_ids):,}")
    print(f"  Samples: {len(sample_ids):,}")
    print(f"  Total reads: {table.matrix_data.sum():,.0f}")

    print("\nLoading taxonomy and finding red complex ASVs...")
    tax_df = load_taxonomy()
    print(f"  Taxonomy entries: {len(tax_df):,}")
    rc_asvs = find_rc_asvs(obs_ids, tax_df)

    print("\nExtracting red complex abundances...")
    rc_df = extract_rc_abundance(table, rc_asvs, sample_ids)

    print("\nMerging with metadata...")
    merged = merge_with_metadata(rc_df, DATA_DIR / "metadata_oral_qc.tsv")

    print_summary(merged)

    # Save
    output_path = OUTDIR / "red_complex_sample_abundance.tsv"
    merged.to_csv(output_path, sep="\t", index=False)
    print(f"\nSaved: {output_path}")
    print(f"  Rows: {len(merged):,}")
    print(f"  Columns: {len(merged.columns)}")
    print("Done.")


if __name__ == "__main__":
    main()
