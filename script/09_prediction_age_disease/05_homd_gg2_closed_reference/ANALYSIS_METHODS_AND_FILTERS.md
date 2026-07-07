# Analysis methods and filtering criteria

This document summarizes the final analysis choices used in `HOMD_ASV_representation_report.html`.

## 1. Inputs

- Retained dataset list and metadata: `/beegfs/huangxiaochang/global_oral/metadata_v2_merge/metadata-v2-merge-oral.csv`
- Per-project ASV FASTA: `/beegfs/dongbiao/oral/<BioProject>/results/ASVs.fa`
- Per-project feature abundance table: `/beegfs/dongbiao/oral/<BioProject>/results/ASVs_counts.tsv`
- HOMD 16S database: `/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03_full.fasta`
- Greengenes2 sequences: `/beegfs/dongbiao/greengene2/2024.09.seqs.fna`
- Greengenes2 taxonomy: `/beegfs/dongbiao/greengene2/exported-taxonomy/taxonomy.tsv`
- VSEARCH: `/home/xuhuan/miniconda3/envs/vsearch_env/bin/vsearch`

## 2. ASV collection and exact dereplication

- DADA2 ASV IDs were project-local, so each ASV was renamed as `<BioProject>|<ASV_ID>`.
- Abundance was taken from `ASVs_counts.tsv`.
- Per-ASV total reads were calculated by summing counts across samples.
- Per-ASV sample prevalence was calculated as the number of samples with nonzero counts.
- Identical ASV sequences were collapsed by exact sequence identity; no 97% or 99% clustering was performed.
- Final collection scale:
  - BioProjects: 334
  - Original ASVs: 3,469,672
  - Exact dereplicated sequences: 2,904,097
  - FASTA-matched reads: 6,720,571,645
  - Count-table features missing from FASTA: 0

## 3. HOMD mapping and capture thresholds

- Dereplicated ASVs were mapped to HOMD using `vsearch --usearch_global`.
- For each dereplicated ASV, the best HOMD hit was chosen by highest identity, then highest query coverage, then longest alignment length.
- Primary capture threshold: identity >=99% and query coverage >=99%.
- Supplementary relaxed threshold: identity >=97% and query coverage >=99%.
- Uncaptured ASV: no HOMD hit or failure to meet the selected identity/query-coverage threshold.

## 4. Global capture weighting schemes

- Exact sequence: each dereplicated sequence has weight 1.
- ASV-weighted: each dereplicated sequence is weighted by the number of original source ASV features.
- Read-weighted: each dereplicated sequence is weighted by total reads.
- Prevalence-weighted: each dereplicated sequence is weighted by sample-level prevalence.

Important: the global prevalence-weighted metric uses sample-level prevalence. Later low-HOMD filtering uses subject-level prevalence.

## 5. Low-HOMD ASV filtering

Initial low-HOMD definition used the relaxed 97% threshold:

- no HOMD hit, or
- HOMD identity <97%, or
- query coverage <99%.

Filtering steps:

- Initial low-HOMD dereplicated sequences: 1,898,185
- Retain `n_bioprojects >= 3` and `total_reads >= 100`: 26,126
- Recalculate prevalence using `subject_id`, not sample ID.
- If subject ID was missing, sample ID was used as a conservative fallback.
- Retain `subject_prevalence >= 50`: 1,416

## 6. GG2 annotation and bacterial QC

- The 1,416 recurrent subject-prevalent low-HOMD ASVs were mapped to GG2.
- A permissive VSEARCH search floor was used only to recover best hits for annotation/QC, not as a biological capture threshold.
- Because some GG2 FASTA IDs lacked taxonomy entries, two hit concepts were retained:
  - absolute best alignment hit for audit;
  - best taxonomy-resolved hit for QC and taxonomic summaries.
- Removed QC classes: unassigned, chloroplast, mitochondria, archaea, eukaryote/non-bacterial.
- Final retained bacterial low-HOMD ASVs: 736.
- These 736 ASVs represented 6,076 source ASVs, 27,216,002 reads, and subject-prevalence burden 89,464.

## 7. Low-HOMD bacterial ASV visualization

- ASVs were ordered by subject-level prevalence.
- Heatmap tracks included:
  - number of BioProjects;
  - log10(total reads + 1);
  - HOMD best identity bin;
  - GG2 best identity bin;
  - geographic region detection;
  - oral sampling site detection.
- Identity bins: <90%, 90-95%, 95-97%, 97-99%, >=99%.
- Top 10 ASVs were labelled using GG2 species if available, otherwise GG2 genus or phylum.

## 8. Cohort-level accumulation analysis

- BioProjects were randomly ordered 100 times.
- After each added BioProject, cumulative all ASVs, HOMD-captured ASVs, and HOMD-uncaptured ASVs were counted.
- Curves report median and 95% interval across random orders.
- Primary threshold: identity >=99% and query coverage >=99%.
- Supplementary threshold: identity >=97% and query coverage >=99%.
- Caveat: because BioProjects differ in sample size, this should be interpreted as cohort-level accumulation, not strict subject-level rarefaction.

## 9. Reproduction order

Typical reproduction order from scripts in `code/`:

1. `01_collect_asvs_and_abundance.py`
2. `02_align_homd_and_summarize.py` or `03_summarize_existing_homd_coverage.py`
3. `06_plot_homd_coverage.py`
4. `14_filter_low_identity_subject_prevalence_upset.py`
5. `15_annotate_filtered_low_identity_subject_asvs.py`
6. `16_plot_bacterial_low_identity_region_upset.py`
7. `17_plot_bacterial_asv_prevalence_profile.R`
8. `18_homd_representation_saturation.py`

R plotting environment:

```bash
module load miniconda/4.9.2
source /home/software/anaconda3/etc/profile.d/conda.sh
conda activate renv
```
