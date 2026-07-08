# Metadata and BIOM Tables: V2 QC Release

This directory contains the quality-controlled metadata tables and closed-reference feature tables used for the Global Oral Microbiome resource. The metadata files are tab-delimited text files (`.tsv`), the feature tables are BIOM v2 files (`.biom`), and the QIIME 2 artifacts are provided as `.qza` files.

The QC workflow starts from the merged sample metadata, adds per-sample QC metrics, applies study-level and sample-level filters, and then exports matched metadata and feature tables.

## QC Overview

The merged raw metadata table was generated from the curated metadata base table and QC metrics. The metadata base table defines the authoritative sample, subject, and BioProject set.

Current filtering rules:

- Study-level filters are evaluated using oral samples only.
- Remove studies where more than 30% of oral samples have chimera rate >25%.
- Remove studies where more than 50% of oral samples have ASV count <=30.
- Remove studies where aggregate GG2 closed-reference mapping rate is <30%.
- Remove samples with chimera rate >25% or missing chimera data.
- Remove samples with ASV count <=30 or missing ASV count.
- Remove samples with GG2 closed-reference mapping rate <30% or missing mapping-rate data.
- The oral deduplicated table keeps one sample per `subject_id` x `site_cleaned` combination.

Closed-reference feature tables were generated against the GG2 v2 reference profile. The sample-level closed-reference mapping rate is recorded in the metadata column `sample_cr_mapping_rate`; the study-level rate is recorded in `study_cr_mapping_rate`.

## Files

### Metadata Tables

| File | Samples | Studies | Subjects | Description |
|---|---:|---:|---:|---|
| `metadata_all_for_qc_raw.tsv` | 116,643 | 334 | 45,880 | Raw merged metadata before QC filtering. This table preserves the complete authoritative metadata base and adds QC metrics such as chimera percentage, ASV count, raw reads, closed-reference mapping rate, and alpha diversity metrics. |
| `metadata_all_qc.tsv` | 93,647 | 293 | 42,643 | QC-filtered metadata for all body sites. This table matches `biom_all_qc.biom` and `biom_all_qc.qza`. |
| `metadata_oral_qc.tsv` | 73,340 | 293 | 41,047 | QC-filtered oral subset before subject/site deduplication. This table matches `biom_oral_qc.biom`. |
| `metadata_oral_qc_dedup.tsv` | 44,446 | 286 | 41,047 | QC-filtered oral subset deduplicated to one sample per `subject_id` x `site_cleaned`. This table matches `biom_oral_qc_dedup.biom`. |

Important metadata columns include:

- `bioproject_id`: BioProject or study accession.
- `sample_id`: sample accession or sample identifier.
- `subject_id`: subject identifier; placed as the third column.
- `chim_perc`: sample chimera percentage.
- `asv_count`: number of ASVs detected in the sample.
- `total_reads`: sample sequencing read count used for QC summaries.
- `sample_cr_mapping_rate`: sample-level GG2 closed-reference mapping rate.
- `study_cr_mapping_rate`: study-level GG2 closed-reference mapping rate.
- `Observed`, `Shannon`, `Simpson`, `InvSimpson`: alpha diversity metrics.
- `site_cleaned` and `site_group_cleaned`: harmonized body-site annotations.
- `platform_cleaned` and `sequencing_region_cleaned`: harmonized sequencing platform and 16S region annotations.

### BIOM Feature Tables

| File | Features | Samples | Total Count | Matches Metadata | Description |
|---|---:|---:|---:|---|---|
| `biom_all_qc.biom` | 61,235 | 93,647 | 5,585,338,653 | `metadata_all_qc.tsv` | Closed-reference GG2 feature table for all QC-passed samples across all body sites. |
| `biom_oral_qc.biom` | 61,235 | 73,340 | 4,417,933,287 | `metadata_oral_qc.tsv` | Closed-reference GG2 feature table for QC-passed oral samples. |
| `biom_oral_qc_dedup.biom` | 61,235 | 44,446 | 2,772,425,709 | `metadata_oral_qc_dedup.tsv` | Closed-reference GG2 feature table for the deduplicated oral sample set. |

### QIIME 2 Artifacts

| File | QIIME 2 Type | Description |
|---|---|---|
| `biom_all_qc.qza` | `FeatureTable[Frequency]` | QIIME 2 artifact version of the all-site QC-passed closed-reference feature table. It corresponds to `biom_all_qc.biom` and `metadata_all_qc.tsv`. |
| `genus_table_all_qc.qza` | `FeatureTable[Frequency]` | Genus-level abundance table generated from the QC-passed all-site feature table. Values are counts/frequencies. |
| `genus_table_relabund_all_qc.qza` | `FeatureTable[RelativeFrequency]` | Genus-level relative-abundance table generated from the QC-passed all-site feature table. Values are relative frequencies. |

## Recommended Usage

Use matched metadata and feature-table pairs together:

- All QC-passed samples: `metadata_all_qc.tsv` + `biom_all_qc.biom` or `biom_all_qc.qza`.
- Oral QC-passed samples: `metadata_oral_qc.tsv` + `biom_oral_qc.biom`.
- Deduplicated oral samples: `metadata_oral_qc_dedup.tsv` + `biom_oral_qc_dedup.biom`.
- Genus-level all-site analysis: `metadata_all_qc.tsv` + `genus_table_all_qc.qza` or `genus_table_relabund_all_qc.qza`.

For analyses sensitive to repeated sampling, use the deduplicated oral pair. For longitudinal, multi-site, or repeated-measure analyses, use the non-deduplicated metadata and explicitly model `subject_id`.
