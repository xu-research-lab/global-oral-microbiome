R Notebook: Metadata Processing Workflow
================

# Step 1: Obtain Initial Metadata from Each Project

This step sets the project directory, loads the raw SRA metadata, and
optionally merges supplement tables.

``` r
# =============================================================
# Curation Notes - PRJNA552608
# bioproject     : From initial_meta column 'BioProject'
# sample_id      : From initial_meta column 'Run' (SRR IDs)
# subject_id     : From 'Sample Name' letter prefix; single letter A-P
# raw_subject_id : original 'Sample Name' (e.g., A_1_1, B_3_4)
# age            : NA in SRA metadata; paper Section 2.1 gives age range but not per-subject
# age_info       : Paper Section 2.1 - overall range 12y9m-33y3m, recorded as "12-33"
# bmi            : NA - not reported in paper or metadata
# sex            : NA - SRA has no sex column; paper gives group-level counts only
# sequencing_region : Paper Section 2.4 - "V1-3 region of 16S rRNA"
# platform       : Paper Section 2.4 - "MiSeq platform using v3-600 cycle reagents"
# country        : Paper affiliation - UCSF, San Francisco, CA, USA
# site           : From 'Sample Name' last number (1-4); mapped to subgingival sites
#                  per paper Figure 1A: 1=34B, 2=46B, 3=34L, 4=41F
# site_group     : "oral" - all subgingival samples
# public_time    : From initial_meta column 'ReleaseDate'
# pcr_primers    : NA - paper does NOT list primer sequences
# title          : Paper title
# disease_study  : "yes" - studies effect of orthodontic treatment
# status         : "HC" for control (J-N); "full-fixed orthodontic treatment" for braces (A-I, O-P)
# matched_samples: "yes" - same subjects at 3 timepoints, 4 sites each
# is_time_series : "no" - intervention study (before/after comparison), not natural time series
# controls       : 5 control subjects (letters J-N); no negative controls in SRA table
# non-human      : none - all Homo sapiens
# =============================================================

# Load libraries
library(dplyr)
library(ggplot2)
library(cowplot)

# (1) Setting working directory
# EDIT HERE: change project_id for each study
project_id <- "PRJNA552608"

base_dir <- "D:/metadata_collect/oral_metadata_new/update"
setwd(file.path(base_dir, project_id))

# (2) Loading raw SRA metadata file
sra_file <- "SraRunTable.csv"
initial_meta <- read.csv(sra_file, header = TRUE)

# (3) OPTIONAL: If extra metadata (e.g., Table.x or supplemental materials) are available from the article, you may need to merge the file with initial_meta manually. Below is example code for merging. Depending on the specific case of each study, you can customize the code as needed.
# supp_table <- read.csv("supplemental_table.csv", header = TRUE)
# initial_meta <- left_join(initial_meta, supp_table, by = merge_col_name)
```

# Step 2: Data Overview

First, you need to get an overview of the data, which includes two
steps: reading the article and reading the initial metadata. This step
is very IMPORTANT!

``` r
# (1) Read the article in the following order: abstract, the beginning of the results section, and the methods. Understand the sample size, experimental design, and disease group overview. Clarify the sampling sites; for oral and gut samples, specify the exact sites like supragingival plaque, tongue, etc.

# (2) IMPORTANT: MUST go through each column of initial_meta
head(initial_meta, 10)
```

    ##           Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample
    ## 1  SRR9654724   AMPLICON        600 23353925 PRJNA552608 SAMN12213982
    ## 2  SRR9654725   AMPLICON        601 15566240 PRJNA552608 SAMN12213983
    ## 3  SRR9654726   AMPLICON        601  8076463 PRJNA552608 SAMN12213984
    ## 4  SRR9654727   AMPLICON        601 24296536 PRJNA552608 SAMN12213985
    ## 5  SRR9654728   AMPLICON        601 16995408 PRJNA552608 SAMN12213986
    ## 6  SRR9654729   AMPLICON        600  3821500 PRJNA552608 SAMN12213987
    ## 7  SRR9654730   AMPLICON        601 14411824 PRJNA552608 SAMN12213988
    ## 8  SRR9654731   AMPLICON        600    31231 PRJNA552608 SAMN12213989
    ## 9  SRR9654732   AMPLICON        600 26696745 PRJNA552608 SAMN12213990
    ## 10 SRR9654733   AMPLICON        600 20556249 PRJNA552608 SAMN12213991
    ##                 BioSampleModel    Bytes               Center.Name
    ## 1  Metagenome or environmental 14585848 UNIVERSITY OF CONNECTICUT
    ## 2  Metagenome or environmental  9220792 UNIVERSITY OF CONNECTICUT
    ## 3  Metagenome or environmental  4776330 UNIVERSITY OF CONNECTICUT
    ## 4  Metagenome or environmental 14336751 UNIVERSITY OF CONNECTICUT
    ## 5  Metagenome or environmental 10402737 UNIVERSITY OF CONNECTICUT
    ## 6  Metagenome or environmental  2491116 UNIVERSITY OF CONNECTICUT
    ## 7  Metagenome or environmental  9056544 UNIVERSITY OF CONNECTICUT
    ## 8  Metagenome or environmental   103366 UNIVERSITY OF CONNECTICUT
    ## 9  Metagenome or environmental 16508461 UNIVERSITY OF CONNECTICUT
    ## 10 Metagenome or environmental 13074726 UNIVERSITY OF CONNECTICUT
    ##    Collection_Date Consent DATASTORE.filetype DATASTORE.provider
    ## 1               NA  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 2               NA  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 3               NA  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 4               NA  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5               NA  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 6               NA  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 7               NA  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 8               NA  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 9               NA  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 10              NA  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                        DATASTORE.region Experiment geo_loc_name_country
    ## 1  gs.us-east1,ncbi.public,s3.us-east-1 SRX6416250         uncalculated
    ## 2  gs.us-east1,ncbi.public,s3.us-east-1 SRX6416249         uncalculated
    ## 3  gs.us-east1,ncbi.public,s3.us-east-1 SRX6416248         uncalculated
    ## 4  gs.us-east1,ncbi.public,s3.us-east-1 SRX6416247         uncalculated
    ## 5  gs.us-east1,ncbi.public,s3.us-east-1 SRX6416246         uncalculated
    ## 6  gs.us-east1,ncbi.public,s3.us-east-1 SRX6416245         uncalculated
    ## 7  gs.us-east1,ncbi.public,s3.us-east-1 SRX6416244         uncalculated
    ## 8  gs.us-east1,ncbi.public,s3.us-east-1 SRX6416243         uncalculated
    ## 9  gs.us-east1,ncbi.public,s3.us-east-1 SRX6416242         uncalculated
    ## 10 gs.us-east1,ncbi.public,s3.us-east-1 SRX6416241         uncalculated
    ##    geo_loc_name_country_continent geo_loc_name         HOST     Instrument
    ## 1                    uncalculated           NA Homo sapiens Illumina MiSeq
    ## 2                    uncalculated           NA Homo sapiens Illumina MiSeq
    ## 3                    uncalculated           NA Homo sapiens Illumina MiSeq
    ## 4                    uncalculated           NA Homo sapiens Illumina MiSeq
    ## 5                    uncalculated           NA Homo sapiens Illumina MiSeq
    ## 6                    uncalculated           NA Homo sapiens Illumina MiSeq
    ## 7                    uncalculated           NA Homo sapiens Illumina MiSeq
    ## 8                    uncalculated           NA Homo sapiens Illumina MiSeq
    ## 9                    uncalculated           NA Homo sapiens Illumina MiSeq
    ## 10                   uncalculated           NA Homo sapiens Illumina MiSeq
    ##    lat_lon Library.Name LibraryLayout LibrarySelection LibrarySource
    ## 1       NA        B_3_1        PAIRED              PCR   METAGENOMIC
    ## 2       NA        B_3_2        PAIRED              PCR   METAGENOMIC
    ## 3       NA        B_3_3        PAIRED              PCR   METAGENOMIC
    ## 4       NA        B_3_4        PAIRED              PCR   METAGENOMIC
    ## 5       NA        C_1_1        PAIRED              PCR   METAGENOMIC
    ## 6       NA        C_1_2        PAIRED              PCR   METAGENOMIC
    ## 7       NA        C_1_3        PAIRED              PCR   METAGENOMIC
    ## 8       NA        C_1_4        PAIRED              PCR   METAGENOMIC
    ## 9       NA        C_2_1        PAIRED              PCR   METAGENOMIC
    ## 10      NA        C_2_2        PAIRED              PCR   METAGENOMIC
    ##                 Organism Platform          ReleaseDate          create_date
    ## 1  human oral metagenome ILLUMINA 2019-12-31T00:00:00Z 2019-07-08T15:06:00Z
    ## 2  human oral metagenome ILLUMINA 2019-12-31T00:00:00Z 2019-07-08T15:06:00Z
    ## 3  human oral metagenome ILLUMINA 2019-12-31T00:00:00Z 2019-07-08T15:06:00Z
    ## 4  human oral metagenome ILLUMINA 2019-12-31T00:00:00Z 2019-07-08T15:06:00Z
    ## 5  human oral metagenome ILLUMINA 2019-12-31T00:00:00Z 2019-07-08T15:06:00Z
    ## 6  human oral metagenome ILLUMINA 2019-12-31T00:00:00Z 2019-07-08T15:06:00Z
    ## 7  human oral metagenome ILLUMINA 2019-12-31T00:00:00Z 2019-07-08T15:06:00Z
    ## 8  human oral metagenome ILLUMINA 2019-12-31T00:00:00Z 2019-07-08T15:06:00Z
    ## 9  human oral metagenome ILLUMINA 2019-12-31T00:00:00Z 2019-07-08T15:06:00Z
    ## 10 human oral metagenome ILLUMINA 2019-12-31T00:00:00Z 2019-07-08T15:06:00Z
    ##    version Sample.Name source_material_id SRA.Study
    ## 1        1       B_3_1                 21 SRP213699
    ## 2        1       B_3_2                 22 SRP213699
    ## 3        1       B_3_3                 23 SRP213699
    ## 4        1       B_3_4                 24 SRP213699
    ## 5        1       C_1_1                 25 SRP213699
    ## 6        1       C_1_2                 26 SRP213699
    ## 7        1       C_1_3                 27 SRP213699
    ## 8        1       C_1_4                 28 SRP213699
    ## 9        1       C_2_1                 29 SRP213699
    ## 10       1       C_2_2                 30 SRP213699

``` r
# IMPORTANT: Check the unique values in the column of interest or query column.
# No body_habitat/body_site columns in this SRA table — site info must be derived from Sample Name
table(initial_meta$Organism)          # all "human oral metagenome"
```

    ## 
    ## human oral metagenome 
    ##                   192

``` r
table(initial_meta$HOST)              # all "Homo sapiens"
```

    ## 
    ## Homo sapiens 
    ##          192

``` r
table(initial_meta$Platform)          # all "ILLUMINA"
```

    ## 
    ## ILLUMINA 
    ##      192

``` r
table(initial_meta$Instrument)        # all "Illumina MiSeq"
```

    ## 
    ## Illumina MiSeq 
    ##            192

``` r
table(initial_meta$geo_loc_name_country)  # all "uncalculated" — use paper for country
```

    ## 
    ## uncalculated 
    ##          192

``` r
unique(initial_meta$Sample.Name)      # pattern: {Letter}_{SubjectNum}_{ReplicateNum}
```

    ##   [1] "B_3_1" "B_3_2" "B_3_3" "B_3_4" "C_1_1" "C_1_2" "C_1_3" "C_1_4" "C_2_1"
    ##  [10] "C_2_2" "H_2_2" "H_1_2" "L_3_1" "L_3_2" "L_3_3" "L_3_4" "M_1_1" "M_1_2"
    ##  [19] "M_1_3" "M_1_4" "M_2_1" "M_2_2" "B_2_4" "B_2_3" "B_1_4" "B_1_3" "B_2_2"
    ##  [28] "B_2_1" "A_3_4" "A_3_3" "B_1_2" "B_1_1" "H_3_4" "H_3_3" "I_1_2" "I_1_1"
    ##  [37] "H_2_4" "H_2_3" "N_1_4" "N_1_3" "N_1_2" "N_1_1" "M_3_4" "M_3_3" "M_3_2"
    ##  [46] "M_3_1" "M_2_4" "M_2_3" "I_1_4" "I_1_3" "A_1_4" "H_1_3" "H_1_4" "O_1_1"
    ##  [55] "O_1_2" "G_3_3" "G_3_4" "G_3_1" "G_3_2" "N_2_3" "N_2_4" "N_2_1" "N_2_2"
    ##  [64] "N_3_3" "N_3_4" "N_3_1" "N_3_2" "G_2_4" "G_2_3" "G_1_2" "G_1_1" "F_3_4"
    ##  [73] "F_3_3" "G_2_2" "G_2_1" "G_1_4" "G_1_3" "H_3_2" "O_3_4" "O_3_3" "O_1_4"
    ##  [82] "O_1_3" "O_2_2" "O_2_1" "O_2_4" "O_2_3" "O_3_2" "O_3_1" "H_3_1" "F_2_1"
    ##  [91] "F_2_2" "F_2_3" "F_2_4" "F_1_1" "F_1_2" "F_1_3" "F_1_4" "F_3_1" "F_3_2"
    ## [100] "P_2_1" "P_2_2" "P_2_3" "P_2_4" "P_1_1" "P_1_2" "P_1_3" "P_1_4" "P_3_1"
    ## [109] "P_3_2" "A_2_4" "A_2_3" "A_2_2" "A_2_1" "J_1_1" "J_1_2" "A_1_2" "A_1_1"
    ## [118] "I_3_1" "I_3_2" "I_3_3" "I_3_4" "I_2_1" "I_2_2" "I_2_3" "I_2_4" "E_1_4"
    ## [127] "E_1_3" "A_1_3" "E_2_1" "E_2_4" "E_2_3" "E_3_2" "E_3_1" "E_3_4" "E_3_3"
    ## [136] "H_2_1" "P_3_4" "P_3_3" "J_3_4" "J_3_3" "J_2_2" "J_2_1" "J_1_4" "J_1_3"
    ## [145] "J_3_2" "J_3_1" "J_2_4" "J_2_3" "D_2_3" "D_2_4" "D_2_1" "D_2_2" "D_3_3"
    ## [154] "D_3_4" "D_3_1" "D_3_2" "E_1_1" "E_1_2" "E_2_2" "K_2_3" "K_2_4" "K_2_1"
    ## [163] "K_2_2" "K_1_3" "K_1_4" "K_1_1" "K_1_2" "K_3_1" "K_3_2" "A_3_2" "D_1_2"
    ## [172] "D_1_1" "C_3_4" "C_3_3" "C_3_2" "C_3_1" "C_2_4" "C_2_3" "D_1_4" "D_1_3"
    ## [181] "A_3_1" "L_1_4" "L_1_3" "L_2_2" "L_2_1" "K_3_4" "K_3_3" "L_1_2" "L_1_1"
    ## [190] "H_1_1" "L_2_4" "L_2_3"

``` r
table(initial_meta$LibraryLayout)     # all "PAIRED"
```

    ## 
    ## PAIRED 
    ##    192

``` r
table(initial_meta$LibrarySelection)  # all "PCR"
```

    ## 
    ## PCR 
    ## 192

``` r
table(initial_meta$Assay.Type)        # all "AMPLICON"
```

    ## 
    ## AMPLICON 
    ##      192

# Step 3: Removing Control and Non-human Samples

First, we need to remove control samples(including both negative
controls and positive controls) and non-human samples.

``` r
# (1) EDIT HERE: Parameter setup
# No explicit control column in metadata. Controls identified by paper: 5 control subjects
# (no braces). No negative/positive control rows visible in SRA table.
# If negative controls exist in SRA they would be in Sample Name; none found.
control_column <- ""   # no control column in initial_meta
control_keywords <- c()     # no keywords needed

# All samples are Homo sapiens per HOST column; no non-human samples
species_column <- "Organism"
unique(initial_meta[[species_column]])
```

    ## [1] "human oral metagenome"

``` r
non_human_keywords <- c()            # no non-human samples — all Homo sapiens


# (2) Remove controls and non-human samples. DO NOT need to modify code in this part, just run it.
n_samples <- nrow(initial_meta)
# Control evaluation
if (length(control_keywords) > 0) {
  if (is.null(control_column) || control_column == "") {
    stop("Error: control_keywords provided but control_column is missing.")
  }
  
  if (!(control_column %in% colnames(initial_meta))) {
    stop(sprintf("Error: Column '%s' not found in initial_meta.", control_column))
  }
  
  control_pattern <- paste(control_keywords, collapse = "|")
  control_text <- initial_meta[[control_column]]
  is_control <- grepl(control_pattern, control_text)
  
} else {
  is_control <- rep(FALSE, n_samples)
}

# Non-human evaluation
if (length(non_human_keywords) > 0) {
  if (is.null(species_column) || species_column == "") {
    stop("Error: non_human_keywords provided but species_column is missing.")
  }
  
  if (!(species_column %in% colnames(initial_meta))) {
    stop(sprintf("Error: Column '%s' not found in initial_meta.", species_column))
  }
  
  non_human_pattern <- paste(non_human_keywords, collapse = "|")
  species_text <- initial_meta[[species_column]]
  is_non_human <- grepl(non_human_pattern, species_text)
} else {
  is_non_human <- rep(FALSE, n_samples)
}

# Execute removal
samples_to_remove <- is_control | is_non_human
initial_meta <- initial_meta[!samples_to_remove, ]

cat("Initial samples counts:",n_samples,"\n")
```

    ## Initial samples counts: 192

``` r
cat("Removed control samples:", sum(is_control), "\n")
```

    ## Removed control samples: 0

``` r
cat("Removed non-human samples:", sum(is_non_human), "\n")
```

    ## Removed non-human samples: 0

``` r
cat("Total samples removed ", sum(samples_to_remove), "\n")
```

    ## Total samples removed  0

``` r
cat("Remaining samples after Step 2:", nrow(initial_meta), "\n")
```

    ## Remaining samples after Step 2: 192

# Step 4: Creating a Uniform Metadata Template

This step converts the raw SRA table into the standard metadata format
used in our meta-analysis. The uniform metadata should include the
following columns: 1. **bioproject**: Project ID from SRA or other
database e.g. `PRJNA1234`, `CNP1234`, ‘OEP1234’, ‘qiita_1234’ 2.
**sample_id**: SRA sample/run ID. Unique id for each sample!
e.g. `DRR1234567`, `ERR1234567` and `SRR1234567`. If no such column
exists, can use column like sample_name. 3. **subject_id**: subject id
for each person if available. e.g. `patient_01`, `Host_03_before`,
`subject12` 4. **age**: Age of the subject if available. e.g. `35`, `42`
\| not `adult`, `chilren` 5. **age_info**: If the age is not explicitly
stated but approximate age range information is available, mark it in
this column. e.g., `18-55`, `children`, `infant`, `elder` 6. **bmi**:
Body Mass Index if available. e.g. `22.5`, `18.9` \| not `normal`,
`overweight` 7. **sex**: Sex of the subject. Mandatory to write as M/F
\| not `female`, `woman` 8. **sequencing_region**: Targeted 16S rRNA
region. e.g. `V3-V4`, `V4` 9. **platform**: Sequencing platform.
e.g. `Illumina MiSeq`, `454 GS FLX Titanium`, `Ion Torrent` 10.
**country**: Country of sample origin. e.g. `China`, `Japan` \| not city
name like `Beijing` 11. **site**: Sampling site. For oral and gut sites,
specify the exact sites as detailed as possible. For example, if it is
known to be supragingival plaque, do not record it as ‘dental plaque’,
unless the original text provides only a vague description.
e.g. `saliva`, `tongue coating` \| not broad term like `oral` 12.
**site_group**: only accepts one of following three values: ‘oral’,
‘gut’, or ‘other’. 13. **public_time**: Data release date.
e.g. `2021-06-22` \| not `June 2021` 14. **pcr_primers**: PCR primers
used. e.g. `515F-806R`, `27F-338R` \| not primer sequences 15.
**title**: Title of the associated publication. 16. **disease_study**:
Whether the study involves disease samples. If it’s not clear, can
record as NA. Only three values accepted: ‘yes’, ‘no’ or NA 17.
**status**: Disease status for each sample if available, else record as
NA. Only accept following abbreviations, Write full names for all other
diseases. HC: healthy control CD: Crohn’s disease CRC: colorectal cancer
COPD: chronic obstructive pulmonary disease HIV: human immunodeficiency
virus ICU: intensive care unit OSCC: oral squamous cell carcinoma RA:
rheumatoid arthritis UC: ulcerative colitis 18. **matched_samples**:
Whether the subject has matched samples. Only accept ‘yes’ or ‘no’. 19.
**is_time_series**: Whether the study is a time series study. Only
accept ‘yes’ or ‘no’.

``` r
# Extract columns from initial_meta as follows: (1) No modification needed → extract directly. (2) Modification needed → create a new column with the prefix raw_ and mark the original column as “pending”. (3) Need to find or determine based on the article.
metadata <- data.frame(
  bioproject = initial_meta$BioProject,
  # From initial_meta column 'BioProject'
  sample_id = initial_meta$Run,
  # From initial_meta column 'Run' (SRR IDs)
  subject_id = initial_meta$Sample.Name,
  # From initial_meta 'Sample Name'; will extract letter portion below
  raw_subject_id = initial_meta$Sample.Name,
  # From initial_meta 'Sample Name'; original sample name as raw
  age = NA,
  # NA - SRA has no age column; paper gives group-level age only
  age_info = "12-33",
  # Paper Section 2.1 - overall age range 12y9m to 33y3m; no per-subject age available
  bmi = NA,
  # NA - not reported
  sex = NA_character_,
  # NA - SRA has no sex column; paper gives group-level sex counts but no per-sample mapping
  sequencing_region = "V1-V3",
  # Paper Section 2.4 - "V1-3 region of 16S rRNA"
  platform = "Illumina MiSeq",
  # Paper Section 2.4 - "MiSeq platform using v3-600 cycle reagents (Illumina)"
  country = "USA",
  # Paper: UCSF, San Francisco, CA, USA
  raw_site = initial_meta$Sample.Name,
  # From 'Sample Name'; letter prefix maps to tooth site per paper Figure 1A
  site = "pending",
  site_group = "pending",
  public_time = initial_meta$ReleaseDate,
  # From initial_meta column 'ReleaseDate'
  pcr_primers = NA_character_,
  # Paper Section 2.4 does NOT list primer sequences, only "V1-3 region of 16S rRNA"
  title = "Alterations in subgingival microbiota during full-fixed appliance orthodontic treatment- A prospective study",
  # Paper title
  disease_study = "yes",
  # Paper studies effect of orthodontic treatment on subgingival microbiota
  status = "pending",
  # Control subjects = HC; braces subjects = full-fixed orthodontic treatment
  matched_samples = "pending",
  is_time_series = "no",
  # Paper is an intervention study (before vs after treatment comparison), not a natural time series
  # Subjects sampled at 3 timepoints to assess intervention effect, not longitudinal dynamics
  stringsAsFactors = FALSE
)
metadata$public_time <- as.Date(metadata$public_time)
```

# Step 5: Manually solving pending issues

This step involves two parts: (1) handle the columns with the `raw_`
prefix, and next (2) handle the `site_group` and `matched_samples`
columns.

``` r
# (1) Handle the columns with the `raw_` prefix.
# Sample Name pattern: {Letter}_{Timepoint}_{Site}
# Letter (A-P) = subject identifier (16 letters; paper says 14 subjects + 2 extra/excluded)
# Middle number (1-3) = timepoint (0wk, 6wk, 12wk)
# Last number (1-4) = site (34B, 46B, 34L, 41F)
# Each letter = 1 subject with 3 timepoints x 4 sites = 12 samples
# Paper Section 2.1: 9 braces + 5 control = 14 subjects; 1 excluded (Patient F)
# 16 letters in SRA = 14 real subjects + 2 extra (excluded subject / additional)
# Assignment: A-I (9) = braces group; J-N (5) = control group; O-P (2) = excluded/extra

# Extract subject_id from the letter prefix of Sample Name
metadata$subject_id <- sapply(strsplit(as.character(metadata$raw_site), "_"), function(x) {
  x[1]
}) # From Sample Name letter prefix; single letter A-P as subject identifier

# Map site from Sample Name: last number = site index
# Paper Figure 1A: 4 mandibular sites - 34B, 46B, 34L, 41F
metadata$site <- sapply(strsplit(as.character(metadata$raw_site), "_"), function(x) {
  site_num <- as.integer(x[3])
  if (site_num == 1) {
    "subgingival_34B"   # site 34 buccal
  } else if (site_num == 2) {
    "subgingival_46B"   # site 46 buccal
  } else if (site_num == 3) {
    "subgingival_34L"   # site 34 lingual
  } else if (site_num == 4) {
    "subgingival_41F"   # site 41 facial
  } else {
    NA_character_
  }
}) # From Sample Name last number; mapped to subgingival tooth sites per paper Figure 1A

metadata$raw_site <- NULL
```

``` r
# (2) Filling the 'site_group' and 'matched_samples' columns
# editing site_group
unique(metadata$site)
```

    ## [1] "subgingival_34B" "subgingival_46B" "subgingival_34L" "subgingival_41F"

``` r
# EDIT HERE: add oral sites — all subgingival sites are oral
Oral_site = c("subgingival_34B", "subgingival_46B", "subgingival_34L", "subgingival_41F")
# EDIT HERE: add gut sites
Gut_site = c()
# just run
metadata <- metadata %>%
  mutate(site_group = case_when(
    site %in% Oral_site ~ "oral",
    site %in% Gut_site ~ "gut",
    TRUE ~ "other"
  ))

# editing matched_samples: subjects have 4 sites x 3 timepoints = multiple samples per subject
metadata$matched_samples <- ifelse(
  ave(metadata$subject_id, metadata$subject_id, FUN = length) > 1,
  "yes",
  "no"
)

# Assign disease status based on group assignment
# Paper Section 2.1: 9 braces subjects + 5 control subjects = 14 total
# Letters A-I (9) = braces group; J-N (5) = control group; O-P (2) = excluded/extra
# Control subjects are healthy (HC); braces subjects undergoing full-fixed orthodontic treatment
# One control subject showed signs of compromised periodontal health at 6w/12w but cannot be identified per sample
metadata$status <- ifelse(
  metadata$subject_id %in% toupper(letters[10:14]),
  "HC",
  "full-fixed orthodontic treatment"
) # Control subjects (J to N) = HC; braces subjects = full-fixed orthodontic treatment; O-P = same
```

``` r
# (3) check if there is pending value in metadata. If there is pending value, go on modify the metadata.
# Function to check if a column contains "pending"
has_pending <- function(col) {
  if (is.character(col) || is.factor(col)) {
    any(col == "pending", na.rm = TRUE)
  } else {
    FALSE
  }
}

# Check for any pending value
any_pending_value <- any(sapply(metadata, has_pending))

if (!any_pending_value) {
  message("Good! No cell contains the value 'pending'.")
} else {
  pending_cols <- names(metadata)[sapply(metadata, has_pending)]
  message("Columns containing 'pending': ", paste(pending_cols, collapse = ", "))
  message("Please continue to resolve these pending columns!")
}
```

    ## Good! No cell contains the value 'pending'.

# Step 6: Visualize the Metadata

This step makes a quick set of plots for checking the distribution of
key metadata fields.

``` r
# Select columns for visualization, just run
plot_data <- metadata %>%
  select(-bioproject, -sample_id, -subject_id, -title)

create_simple_plots <- function(data) {
  plots <- list()
  continuous_vars <- c("age", "bmi")

  for (col in colnames(data)) {
    if (all(is.na(data[[col]]))) {
      p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "this column is NA", size = 6) +
        labs(title = col) +
        theme_void()
    } else if (col %in% continuous_vars && is.numeric(data[[col]])) {
      p <- ggplot(data, aes(x = .data[[col]])) +
        geom_histogram(fill = "blue", alpha = 0.6, bins = 20, na.rm = TRUE) +
        labs(title = col, x = col, y = "Count") +
        theme_minimal()
    } else {
      count_data <- as.data.frame(table(data[[col]], useNA = "ifany"))
      p <- ggplot(count_data, aes(x = reorder(Var1, -Freq), y = Freq)) +
        geom_col(fill = "red", alpha = 0.6) +
        labs(title = col, x = NULL, y = "Count") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 15, hjust = 1))
    }
    plots[[col]] <- p
  }

  return(plots)
}

# Generate and arrange plots
plots <- create_simple_plots(plot_data)
p <- plot_grid(plots$age,
               plots$age_info,
               plots$bmi,
               plots$sex,
               plots$site,
               plots$site_group,
               plots$sequencing_region,
               plots$platform,
               plots$country,
               plots$pcr_primers,
               plots$status,
               plots$matched_samples,
               ncol = 3, align = "hv", axis = "lb")
ggsave("Metadata_Summary_v2.png", plot = p, width = 8, height = 6, dpi = 150)
```

``` r
# Insert Metadata_Summary_v2.png into the current Rmd output
knitr::include_graphics("Metadata_Summary_v2.png")
```

![](Metadata_Summary_v2.png)<!-- -->

# Step 7: Sort and export the metadata

This step saves the cleaned metadata table as the final `metadata.csv`
file.

``` r
# Reorder the data frame
fixed_cols <- c("bioproject", "sample_id", "subject_id", "age", "age_info",
                "bmi", "sex", "platform", "sequencing_region", "pcr_primers", 
                "country", "site", "site_group", "disease_study", "status", 
                "title", "public_time", "matched_samples", "is_time_series")  
metadata <- metadata %>%
  select(all_of(fixed_cols), everything())

# save metadata 
# head(metadata)
write.csv(metadata, file = "metadata_v2.csv", quote = TRUE)
```
