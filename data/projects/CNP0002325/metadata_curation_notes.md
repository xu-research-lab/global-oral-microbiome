# Metadata Curation Notes - CNP0002325

## Paper
- **Title**: From Mouth to Brain: Distinct Supragingival Plaque Microbiota Composition in Cerebral Palsy Children With Caries
- **Journal**: Frontiers in Cellular and Infection Microbiology, April 2022
- **DOI**: 10.3389/fcimb.2022.814473

## Field Sources
- **bioproject**: From initial_meta column 'project_accession' = CNP0002325
- **sample_id**: From initial_meta column 'run_accession' (e.g., CNR0462543)
- **subject_id**: From initial_meta column 'sample_name' (e.g., B1, S15, C12, D7). Paper Methods states "each subject contributed a single supragingival plaque sample", establishing 1:1 mapping between sample_name and subject.
- **age**: NA — paper reports group means only (CPCS: 10.86+/-2.28, NCPCS: 10.17+/-1.99, CPCM: 11.65+/-1.85, CPCF: 11.44+/-2.12 years), not individual values
- **age_info**: "6-14" — Paper Methods: inclusion criterion 6-14 years old
- **bmi**: NA — not reported in paper or metadata
- **sex**: NA — paper reports group totals only (Table 1), not per-sample assignment
- **sequencing_region**: "V3-V4" — Paper Methods: "16S rRNA V3-V4 region was constructed by the PCR amplification"
- **platform**: "Illumina MiSeq" — Paper Methods (NOT from metadata column 'instrument_model' which says HiSeq 2500; paper is authoritative)
- **country**: "China" — CNP database (China National GeneBank), study conducted in China
- **site**: "supragingival plaque" — Paper Methods: "supragingival plaque on the buccal surface of the right maxillary first molar"
- **site_group**: "oral" — supragingival plaque is an oral site
- **disease_study**: "yes" — study compares cerebral palsy children with caries vs non-CP controls
- **status**: Mapped from sample_name prefix per paper Table 1:
  - B prefix (n=23) -> "NCPCS" (non-CP children with severe caries)
  - S prefix (n=27) -> "CPCS" (CP children with severe caries)
  - C prefix (n=12) -> "CPCM" (CP children with mild-to-moderate caries)
  - D prefix (n=16) -> "CPCF" (CP children without caries)
- **public_time**: NA — not available in CNP metadata
- **pcr_primers**: NA — paper only states "V3-V4 region", does not specify primer names (e.g., 341F-806R)
- **matched_samples**: "no" — Paper Methods: each subject contributed a single supragingival plaque sample
- **is_time_series**: "no" — Paper: cross-sectional design, samples collected November-December 2020

## Control Keywords
- None found in metadata or paper

## Non-human Keywords
- None found; all samples are human

## Notes
- Metadata file is in CNP TSV format (not SRA CSV), columns differ from typical NCBI downloads
- Total samples: 78 (23 NCPCS + 27 CPCS + 12 CPCM + 16 CPCF)
- Platform discrepancy: metadata says "Illumina HiSeq 2500" but paper says "Illumina MiSeq" — using paper as authoritative
- No individual-level age, sex, or BMI data available; only group-level statistics in paper Table 1
