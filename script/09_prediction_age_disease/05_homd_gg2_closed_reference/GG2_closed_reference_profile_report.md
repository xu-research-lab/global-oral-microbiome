# Containment-Assisted GG2 Closed-Reference Profile Report

Generated on 2026-06-27. Updated after mitochondria & chloroplast removal (g__GWA2-37-10 eukaryotic sequences retained), study-level threshold at 30%, and sample-level filtering added.

## Executive Summary

We built a unified closed-reference profile for oral 16S ASVs using Greengenes2 (GG2) after a containment-assisted dereplication step. The workflow started from 3,469,672 project-local ASV features, collapsed identical sequences to 2,904,097 exact-dereplicated ASVs, and then collapsed shorter sequences fully contained in longer sequences to 2,006,569 longest representatives. Using a strict closed-reference criterion of 99% identity and 100% query coverage, 535,078 containment representatives mapped to GG2.

After removing sequences matching mitochondrial and chloroplast backbone references (1,377 containment representatives → 1,430 GG2 targets removed), the final profile contained 122,046 observed GG2 reference features across 116,391 samples. g__GWA2-37-10 eukaryotic sequences were retained.

At the read-abundance level, the median per-sample mapped read fraction was 63.3%, and the overall mapped read fraction across all samples was 61.6%.

**Filtering thresholds:** We applied a 30% dataset-level read retention threshold for study filtering and a 30% mapped read fraction threshold for sample-level filtering. After study-level filtering, 285 BioProjects were retained and 49 were below threshold. Within retained studies, 96,234 samples were retained and 12,081 were dropped at the 30% sample threshold.

## Methods

### Input ASVs and abundance profiles

Per-project DADA2 ASV FASTA files and `ASVs_counts.tsv` abundance tables were used. Because ASV IDs were only project-local, each ASV was tracked as `<BioProject>|<ASV_ID>`. The count table `ASVs_counts.tsv` was used as the abundance source for profile construction.

### Exact dereplication

Identical ASV sequences were collapsed exactly. This step retained the original ASV-to-sequence relationship and produced exact dereplicated ASV IDs (`DEREP_x`). No 97% or 99% de novo clustering was used for this stage.

### Containment-assisted representative selection

To reduce redundant closed-reference searches across different amplicon lengths, exact-dereplicated sequences were further collapsed using VSEARCH. A shorter exact-derep sequence was assigned to a longer representative only if it was fully contained in the longer sequence at 100% identity and 100% query coverage. This step produced longest containment representatives while preserving a mapping from each exact-derep ID to its representative.

### GG2 closed-reference mapping

The longest containment representatives were mapped to GG2 using VSEARCH `--usearch_global` with:

- identity threshold: 99%
- query coverage threshold: 100%
- strand: both
- maxaccepts: 1
- default maxrejects behavior

For profile construction, representative-level GG2 assignments were propagated back to exact-derep ASVs and then to the original project-local ASV features. Counts from `ASVs_counts.tsv` were aggregated to GG2 target IDs for every sample.

### Mitochondria and chloroplast sequence removal

GG2 backbone sequences annotated as mitochondria (`f__Mitochondria`) or chloroplast (`c__Chloroplast`) were identified from the GG2 taxonomy (6,238 backbone sequences). Containment representatives matching these backbone sequences at 100% identity and 100% query coverage (1,377 representatives) were traced to 1,430 GG2 reference IDs via the containment assignment. These reference IDs were removed from the profile, resulting in 57,062 profile rows dropped.

g__GWA2-37-10 (Patescibacteria) sequences were **not** removed and are retained in the final profile.

### Profile QC and filtering

Dataset-level read retention was computed as closed-reference mapped reads divided by total input ASV reads. Two levels of filtering were applied:

- **Study-level**: BioProjects with overall mapped read fraction < 30% were dropped.
- **Sample-level**: Within retained studies, individual samples with mapped read fraction < 30% were dropped.

## ASV and Feature Count Summary

| Stage | Count | Note |
| --- | --- | --- |
| Original project-local ASV features | 3,469,672 | ASV IDs renamed as <BioProject>\|<ASV_ID>. |
| Exact dereplicated ASV sequences | 2,904,097 | Identical sequences collapsed exactly. |
| Longest containment representatives | 2,006,569 | 897,528 shorter exact-derep sequences collapsed into longer representatives. |
| GG2-mapped containment representatives | 535,078 | 99% identity and 100% query coverage closed-reference mapping. |
| After mitochondria & chloroplast removal | 533,701 | 1,377 suspicious containment representatives removed. |
| Original ASV features assigned to GG2 | 1,225,866 | Assignments propagated through exact derep and containment mapping. |
| Observed GG2 reference features in final profile | 122,046 | Unique GG2 target IDs with nonzero abundance across samples (after removal). |

## Profile QC Summary

| Metric | Value |
| --- | --- |
| Retained BioProjects in metadata | 334 |
| BioProjects with nonzero closed-reference profile | 332 |
| BioProjects retained at >=30% dataset-level read retention | 285 |
| BioProjects below 30% dataset-level read retention | 49 |
| Samples with closed-reference reads | 119,181 |
| Median sample-level mapped read fraction | 0.633 |
| Mean sample-level mapped read fraction | 0.587 |
| Overall mapped read fraction across samples | 0.616 |

## Sample-Level Filtering After Study 30% Retention

The table below shows, after retaining studies with >=30% overall read retention (285 studies, 108,315 samples), how many samples would be retained or dropped at different sample-level mapped read fraction thresholds.

| Sample Threshold | Total Samples | Total Retained | Total Dropped | Oral Retained | Oral Dropped |
| --- | --- | --- | --- | --- | --- |
| 30% | 108,315 | 96,234 | 12,081 | 63,776 | 8,544 |
| 40% | 108,315 | 88,399 | 19,916 | 57,838 | 14,482 |
| 50% | 108,315 | 77,794 | 30,521 | 50,214 | 22,106 |

*Oral sites include: saliva, tooth, oral mucosa, tongue, gingiva, oral swab. Samples from "other" sites are included in the total count but not in the oral count.*

## Figures

![Figure 1](/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/containment_derep/unified_gg2_closed_reference_profile/figures/fig_gg2_closed_reference_mapped_read_fraction_by_bioproject.png)

**Interpretation.** Each point is one sample. Boxes summarize sample-level GG2 mapped read fractions within each BioProject. BioProjects are ordered by the median mapped fraction. The red dashed line marks the 30% sample-level threshold. The broad spread shows that closed-reference retention varies strongly by dataset.

![Figure 2](/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/containment_derep/unified_gg2_closed_reference_profile/figures/fig_gg2_closed_reference_project_read_retention_barplot.png)

**Interpretation.** Each bar is one BioProject-level retention value, computed from total mapped reads divided by total input reads for that dataset. The red dashed line marks the 30% retention threshold; 285 datasets are retained and 49 are below the threshold. Sample-level filtering counts (within retained studies at 30%) are also annotated.

![Figure 3](/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/containment_derep/unified_gg2_closed_reference_profile/figures/fig_gg2_closed_reference_feature_rank_decay_two_panel.png)

**Interpretation.** The left panel ranks GG2 targets by read abundance, while the right panel ranks them by sample prevalence. Both axes are log-scaled, showing a small set of highly abundant/prevalent references and a long tail of rare features.

## Example Summary Tables

### Top 10 BioProjects by dataset-level read retention

| bioproject | samples | overall_mapped_fraction | retention_status |
| --- | --- | --- | --- |
| PRJEB51725 | 27 | 0.946 | Retained |
| PRJNA598825 | 30 | 0.938 | Retained |
| PRJNA622300 | 189 | 0.924 | Retained |
| PRJEB11529 | 287 | 0.921 | Retained |
| PRJEB36081 | 94 | 0.900 | Retained |
| PRJNA598080 | 26 | 0.892 | Retained |
| PRJEB29801 | 47 | 0.889 | Retained |
| PRJNA862307 | 48 | 0.888 | Retained |
| PRJNA770171 | 40 | 0.887 | Retained |
| PRJNA664107 | 46 | 0.886 | Retained |

### Bottom 10 BioProjects by dataset-level read retention

| bioproject | samples | overall_mapped_fraction | retention_status |
| --- | --- | --- | --- |
| PRJNA747639 | 149 | 0.031 | Dropped |
| PRJNA255332 | 24 | 0.020 | Dropped |
| PRJDB5944 | 75 | 0.013 | Dropped |
| PRJEB5069 | 13 | 0.012 | Dropped |
| PRJNA278076 | 140 | 0.002 | Dropped |
| PRJNA712952 | 108 | 0.001 | Dropped |
| PRJNA807496 | 52 | 0.001 | Dropped |
| PRJDB3211 | 135 | 0.001 | Dropped |
| PRJEB12953 | 228 | 0.000 | Dropped |
| PRJNA292354 | 40 | 0.000 | Dropped |

### Top 10 GG2 targets by closed-reference reads

| gg2_target | closed_reference_reads | sample_prevalence |
| --- | --- | --- |
| RS-GCF-003942595.1-NZ-RJNL01000019.1 | 205303833 | 24401 |
| GB-GCA-900767615.1-CAAFZU010000193.1 | 107046132 | 45767 |
| RS-GCF-001810545.1-NZ-KV801801.1 | 106045940 | 28272 |
| CP002122 | 97536385 | 40712 |
| RS-GCF-018619375.1-NZ-JAHHFT010000043.1 | 68614059 | 33613 |
| MJ006-2-barcode51-umi102309bins-ubs-5 | 67087967 | 8846 |
| 05084014 | 53326295 | 7091 |
| GB-GCA-900539315.1-UQBL01000026.1 | 50923731 | 49502 |
| RS-GCF-001064485.1-NZ-JVLW01000023.1 | 48059161 | 30076 |
| RS-GCF-013415365.1-NZ-JACBYB010000008.1 | 44796463 | 26030 |

## Key Output Files

- Sparse unified GG2 profile (filtered): `/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/containment_derep/unified_gg2_closed_reference_profile/all_samples_gg2_closed_reference_profile.long.tsv`
- Original ASV to GG2 assignment table: `/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/containment_derep/unified_gg2_closed_reference_profile/original_asv_to_gg2_closed_reference_assignment.tsv`
- Sample-level mapping fraction table: `/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/containment_derep/unified_gg2_closed_reference_profile/all_samples_gg2_closed_reference_mapping_fraction.tsv`
- Project-level retention summary (30% threshold): `/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/containment_derep/unified_gg2_closed_reference_profile/all_projects_gg2_closed_reference_retention_threshold30_summary.tsv`
- GG2 feature summary: `/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/containment_derep/unified_gg2_closed_reference_profile/all_samples_gg2_closed_reference_feature_summary.tsv`
- Containment mapping table: `/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/containment_derep/containment_derep_mapping.tsv`

## Interpretation and Caveats

The strict 99% identity plus 100% query coverage criterion provides a conservative, consistent GG2 reference-ID profile.

Mitochondria (`f__Mitochondria`, 1,189 backbone sequences) and chloroplast (`c__Chloroplast`, 5,049 backbone sequences) were identified from the GG2 backbone taxonomy and removed from the profile. This removed 1,430 GG2 reference IDs (1,377 containment representatives), accounting for 57,062 profile rows (0.68%). g__GWA2-37-10 (Patescibacteria) sequences were retained.

The study-level retention threshold was set at 30% (previously 50%). A parallel sample-level filter at 30% mapped read fraction is recommended. For downstream disease or age prediction, model validation should account for BioProject, sequencing region, platform, and geography.
