# Metadata and BIOM Tables: V4 QC Release

This directory contains the quality-controlled metadata tables and closed-reference feature tables used for the Global Oral Microbiome resource.

The QC workflow starts from the DADA2-processed sample metadata (samples with valid `sample_cr_mapping_rate`), removes samples with 0% CR mapping rate (absent from the closed-reference BIOM), applies dataset-level filters using oral samples only, and then exports matched metadata and feature tables.

## QC Overview (v4)

The merged raw metadata table was generated from the curated metadata base table and QC metrics. Only samples with valid GG2 closed-reference mapping data (DADA2-processed) and CR mapping rate >0 are included as the base sample set.

**Current filtering rules (dataset-level only):**

- All thresholds are evaluated using **oral samples only** (`site_group='oral'`).
- **CR=0 filter**: Remove samples with closed-reference mapping rate == 0 (absent from the closed-reference BIOM table).
- **Chimera filter**: Remove studies where >30% of oral samples have chimera >25%.
- **CR mapping filter**: Remove studies where oral aggregate GG2 closed-reference mapping rate <70%.
- **No per-sample filter flags** applied downstream. ASV-based filtering removed.

## Key Changes from v3

| Item | v3 | v4 |
|---|---|---|
| Base sample set | Full metadata (116,643) | DADA2 + CR>0 (114,959) |
| Filter level | Dataset + Sample | **Dataset only** |
| ASV filter | Dataset-level (oral >50% ASV≤30) | **Removed** |
| CR threshold | 70% CR mapping rate | **70%** CR mapping rate |
| Sample-level filters | ASV, chimera, CR | **CR=0 removal only** |
| Metadata ↔ BIOM alignment | Mismatched | **Perfectly aligned** |

## Files

### Metadata Tables

| File | Samples | Studies | Subjects | Description |
|---|---:|---:|---:|---|
| `metadata_all_for_qc_raw.tsv` | 114,959 | 334 | 45,619 | Raw merged metadata (DADA2 + CR>0). Preserves the authoritative metadata base and adds QC metrics. |
| `metadata_all_qc.tsv` | 108,691 | 288 | 43,506 | QC-filtered metadata for all body sites. |
| `metadata_oral_qc.tsv` | 80,788 | 288 | 42,224 | QC-filtered oral subset before subject/site deduplication. |
| `metadata_oral_qc_dedup.tsv` | 51,186 | 288 | 42,224 | QC-filtered oral subset deduplicated to one sample per study × `subject_id` x `site_cleaned`. |

### BIOM Feature Tables

| File | Features | Samples | Matches Metadata |
|---|---:|---:|---|
| `biom_all_qc.biom` | 61,235 | 108,691 | `metadata_all_qc.tsv` ✓ |
| `biom_oral_qc.biom` | 61,235 | 80,788 | `metadata_oral_qc.tsv` ✓ |
| `biom_oral_qc_dedup.biom` | 61,235 | 51,186 | `metadata_oral_qc_dedup.tsv` ✓ |

### QIIME 2 Artifacts

| File | Type | Description |
|---|---|---|
| `biom_all_qc.qza` | `FeatureTable[Frequency]` | QZA version of `biom_all_qc.biom`. |
| `genus_table_all_qc.qza` | `FeatureTable[Frequency]` | Genus-level collapsed table. *Requires GG2 taxonomy file.* |
| `genus_table_relabund_all_qc.qza` | `FeatureTable[RelativeFrequency]` | Relative abundance genus table. *Requires genus_table_all_qc.qza.* |

### Auxiliary Files

| File | Description |
|---|---|
| `QC/filter_funnel_summary.tsv` | Per-step filter funnel with sample/study counts. |
| `QC/removed_studies_by_step.tsv` | List of studies removed at each filter step. |
| `data_qc_pipeline.md` | QC pipeline report (Markdown) with embedded figures. |

### QC Scripts

| Script | Purpose |
|---|---|
| `merge_all_for_qc_raw.py` | Merge metadata sources → filter DADA2 + CR>0 → compute study-level QC flags. |
| `generate_metadata_tables.py` | Apply dataset-level filters → output metadata + BIOM tables. |
| `generate_qc_report.py` | Generate Markdown QC report with embedded figures. |
| `generate_qza.py` | Convert BIOM to QZA (requires qiime2-2022.2). |

## Recommended Usage

Use matched metadata and feature-table pairs together:

- All QC-passed samples: `metadata_all_qc.tsv` + `biom_all_qc.biom` (or `.qza`).
- Oral QC-passed samples: `metadata_oral_qc.tsv` + `biom_oral_qc.biom`.
- Deduplicated oral samples: `metadata_oral_qc_dedup.tsv` + `biom_oral_qc_dedup.biom`.

For analyses sensitive to repeated sampling, use the deduplicated oral pair. For longitudinal, multi-site, or repeated-measure analyses, use the non-deduplicated metadata and explicitly model `subject_id`.
