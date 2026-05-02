# Metadata Review: PRJEB30836

**Paper**: "Response of the human gut and saliva microbiome to urbanization in Cameroon"
**File**: `data/projects/PRJEB30836/metadata_curation_v2.Rmd`

## Issues Found

### 1. Missing `base_dir` variable; hardcoded `setwd()` path
- **Location**: Line 20
- **Issue**: Uses `setwd("/home/dongbiao/word_embedding_microbiome/all_data/oral/global-oral-microbiome/data/projects/PRJEB30836")` directly instead of a configurable `base_dir` + `project_id` pattern. Should follow the standard pattern:
  ```r
  base_dir <- "/home/dongbiao/word_embedding_microbiome/all_data/oral/global-oral-microbiome/data/projects"
  setwd(file.path(base_dir, project_id))
  ```

### 2. `PCR_primers` set to `NA`
- **Location**: Line 74
- **Issue**: `metadata$PCR_primers <- NA` but the paper Methods and SRA `PCR_primers` column clearly show 515F/806rcbc0 primers targeting the V4 region. Should be `"515F-806R"`.

### 3. `matched_samples` hardcoded to `"YES"` for all samples
- **Location**: Line 91
- **Issue**: Set to `"YES"` unconditionally. Should be computed dynamically based on whether subjects have multiple samples (which they do — each subject has both feces and saliva samples).

### 4. Value format inconsistency
- **Location**: Lines 76, 78, 91
- **Issue**: Uses uppercase `"NO"`, `"YES"`, `"HC"`. The standard meta-analysis template uses lowercase: `"no"`, `"yes"`.

### 5. Platform uses raw `Instrument` column instead of standardized value
- **Location**: Line 59
- **Issue**: `Platform = initial_meta$Instrument` reads directly from the SRA table. The paper specifies "Illumina MiSeq 300". Should use a consistent standardized value.

### 6. Missing columns from standard template
- **Location**: Throughout
- **Issue**: The v2 Rmd does not include these columns from the standard metadata template:
  - `age_info` — age range/grouping information
  - `is_time_series` — whether the study is longitudinal

### 7. `Site` column uses raw values without standardization
- **Location**: Line 61
- **Issue**: `Site = initial_meta$host_body_site` — the raw values are "feces" and "saliva" which are already clean, but this should be verified and documented rather than assumed.

### 8. BMI threshold too low
- **Location**: Line 67
- **Issue**: `metadata$BMI[metadata$BMI > 50] <- NA` — threshold of 50 may incorrectly remove valid extreme BMI values. The standard template uses 100.

### 9. `Age` threshold at 110
- **Location**: Line 66
- **Issue**: `metadata$Age[metadata$Age > 110] <- NA` — standard template uses 100. The SRA data contains values of 999 (placeholder) which are above both thresholds, so this works but should be consistent.

## Summary

| # | Issue | Severity |
|---|-------|----------|
| 1 | Hardcoded `setwd()` path | Medium |
| 2 | `PCR_primers` incorrectly set to `NA` | High |
| 3 | `matched_samples` hardcoded | Medium |
| 4 | Value format (uppercase vs lowercase) | Low |
| 5 | Platform from raw column | Low |
| 6 | Missing standard columns | Medium |
| 7 | Site not verified | Low |
| 8 | BMI threshold inconsistent | Low |
| 9 | Age threshold inconsistent | Low |
