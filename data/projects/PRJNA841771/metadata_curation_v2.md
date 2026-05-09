R Notebook: Metadata Processing Workflow
================

# Step 1: Obtain Initial Metadata from Each Project

This step sets the project directory, loads the raw SRA metadata, and
optionally merges supplement tables.

``` r
# Load libraries
library(dplyr)
library(ggplot2)
library(cowplot)

# (1) Setting working directory
# EDIT HERE: change project_id for each study
project_id <- "PRJNA841771" 

base_dir <- "E:/oral/26/metadata-2604/Oral-Metadata-YL-02"
setwd(file.path(base_dir, project_id))

# (2) Loading raw SRA metadata file
sra_file <- "SraRunTable.csv"
initial_meta <- read.csv(sra_file, header = TRUE)

# (3) OPTIONAL: If extra metadata (e.g., Table.x or supplemental materials) are available from the article, you may need to merge the file with initial_meta manually. Below is example code for merging. Depending on the specific case of each study, you can customize the code as needed.
#supp_table <- read.csv("supplemental_table_v2.csv", header = TRUE)
#initial_meta <- left_join(initial_meta, supp_table, by = "Sample.Name") 
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

    ##            Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample
    ## 1  SRR19370359   AMPLICON        426 42826040 PRJNA841771 SAMN28614607
    ## 2  SRR19370361   AMPLICON        425 39071076 PRJNA841771 SAMN28614605
    ## 3  SRR19370362   AMPLICON        425 33200097 PRJNA841771 SAMN28614604
    ## 4  SRR19370363   AMPLICON        423 37374803 PRJNA841771 SAMN28614603
    ## 5  SRR19370365   AMPLICON        427 35248253 PRJNA841771 SAMN28614602
    ## 6  SRR19370366   AMPLICON        426 34796910 PRJNA841771 SAMN28614601
    ## 7  SRR19370367   AMPLICON        425 30914847 PRJNA841771 SAMN28614600
    ## 8  SRR19370368   AMPLICON        421 35651565 PRJNA841771 SAMN28614599
    ## 9  SRR19370369   AMPLICON        421 34193139 PRJNA841771 SAMN28614598
    ## 10 SRR19370370   AMPLICON        424 31160711 PRJNA841771 SAMN28614597
    ##                                 BioSampleModel    Bytes
    ## 1  MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 13795394
    ## 2  MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 12646566
    ## 3  MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 10674463
    ## 4  MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 12395752
    ## 5  MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 11412076
    ## 6  MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 11821928
    ## 7  MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 10123695
    ## 8  MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 11783946
    ## 9  MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 11376260
    ## 10 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 10321536
    ##                                             Center.Name Collection_Date Consent
    ## 1  CENTER FOR LIFE SCIENCES, NATIONAL LABORATORY ASTANA      2020-08-27  public
    ## 2  CENTER FOR LIFE SCIENCES, NATIONAL LABORATORY ASTANA      2020-08-26  public
    ## 3  CENTER FOR LIFE SCIENCES, NATIONAL LABORATORY ASTANA      2020-08-26  public
    ## 4  CENTER FOR LIFE SCIENCES, NATIONAL LABORATORY ASTANA      2020-08-26  public
    ## 5  CENTER FOR LIFE SCIENCES, NATIONAL LABORATORY ASTANA      2020-08-26  public
    ## 6  CENTER FOR LIFE SCIENCES, NATIONAL LABORATORY ASTANA      2020-08-25  public
    ## 7  CENTER FOR LIFE SCIENCES, NATIONAL LABORATORY ASTANA      2020-08-25  public
    ## 8  CENTER FOR LIFE SCIENCES, NATIONAL LABORATORY ASTANA      2020-08-25  public
    ## 9  CENTER FOR LIFE SCIENCES, NATIONAL LABORATORY ASTANA      2020-08-25  public
    ## 10 CENTER FOR LIFE SCIENCES, NATIONAL LABORATORY ASTANA      2020-08-25  public
    ##    DATASTORE.filetype DATASTORE.provider                     DATASTORE.region
    ## 1    fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 2    fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 3    fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 4    fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 5    fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 6    fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 7    fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 8    fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 9    fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 10   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ##    env_broad_scale env_local_scale    env_medium ETHNICITY  Experiment
    ## 1    ENVO:00000428   ENVO:00000364 ENVO:00010483    kazakh SRX15428860
    ## 2    ENVO:00000428   ENVO:00000362 ENVO:00010483    kazakh SRX15428858
    ## 3    ENVO:00000428   ENVO:00000361 ENVO:00010483    kazakh SRX15428857
    ## 4    ENVO:00000428   ENVO:00000360 ENVO:00010483    kazakh SRX15428856
    ## 5    ENVO:00000428   ENVO:00000359 ENVO:00010483    kazakh SRX15428854
    ## 6    ENVO:00000428   ENVO:00000358 ENVO:00010483    kazakh SRX15428853
    ## 7    ENVO:00000428   ENVO:00000357 ENVO:00010483    kazakh SRX15428852
    ## 8    ENVO:00000428   ENVO:00000356 ENVO:00010483    kazakh SRX15428851
    ## 9    ENVO:00000428   ENVO:00000355 ENVO:00010483    kazakh SRX15428850
    ## 10   ENVO:00000428   ENVO:00000354 ENVO:00010483    kazakh SRX15428849
    ##    geo_loc_name_country geo_loc_name_country_continent           geo_loc_name
    ## 1            Kazakhstan                           Asia Kazakhstan: Nur-Sultan
    ## 2            Kazakhstan                           Asia Kazakhstan: Nur-Sultan
    ## 3            Kazakhstan                           Asia Kazakhstan: Nur-Sultan
    ## 4            Kazakhstan                           Asia Kazakhstan: Nur-Sultan
    ## 5            Kazakhstan                           Asia Kazakhstan: Nur-Sultan
    ## 6            Kazakhstan                           Asia Kazakhstan: Nur-Sultan
    ## 7            Kazakhstan                           Asia Kazakhstan: Nur-Sultan
    ## 8            Kazakhstan                           Asia Kazakhstan: Nur-Sultan
    ## 9            Kazakhstan                           Asia Kazakhstan: Nur-Sultan
    ## 10           Kazakhstan                           Asia Kazakhstan: Nur-Sultan
    ##    Host_age host_body_mass_index host_body_product         Host
    ## 1        33             20.44914       oral cavity Homo sapiens
    ## 2        39             21.64412       oral cavity Homo sapiens
    ## 3        54             21.35991       oral cavity Homo sapiens
    ## 4        51             23.03005       oral cavity Homo sapiens
    ## 5        55             29.51594       oral cavity Homo sapiens
    ## 6        48             27.00513       oral cavity Homo sapiens
    ## 7        34             25.45236       oral cavity Homo sapiens
    ## 8        32             28.12500       oral cavity Homo sapiens
    ## 9        46             21.87500       oral cavity Homo sapiens
    ## 10       41             23.03005       oral cavity Homo sapiens
    ##               Instrument         lat_lon Library.Name LibraryLayout
    ## 1  Illumina NovaSeq 6000 51.10 N 71.26 E       RA116C        SINGLE
    ## 2  Illumina NovaSeq 6000 51.10 N 71.26 E       RA114C        SINGLE
    ## 3  Illumina NovaSeq 6000 51.10 N 71.26 E       RA113C        SINGLE
    ## 4  Illumina NovaSeq 6000 51.10 N 71.26 E       RA112C        SINGLE
    ## 5  Illumina NovaSeq 6000 51.10 N 71.26 E       RA111C        SINGLE
    ## 6  Illumina NovaSeq 6000 51.10 N 71.26 E       RA110C        SINGLE
    ## 7  Illumina NovaSeq 6000 51.10 N 71.26 E       RA109C        SINGLE
    ## 8  Illumina NovaSeq 6000 51.10 N 71.26 E       RA108C        SINGLE
    ## 9  Illumina NovaSeq 6000 51.10 N 71.26 E       RA107C        SINGLE
    ## 10 Illumina NovaSeq 6000 51.10 N 71.26 E       RA106C        SINGLE
    ##    LibrarySelection LibrarySource   Organism Platform          ReleaseDate
    ## 1               PCR   METAGENOMIC metagenome ILLUMINA 2022-05-24T00:00:00Z
    ## 2               PCR   METAGENOMIC metagenome ILLUMINA 2022-05-24T00:00:00Z
    ## 3               PCR   METAGENOMIC metagenome ILLUMINA 2022-05-24T00:00:00Z
    ## 4               PCR   METAGENOMIC metagenome ILLUMINA 2022-05-24T00:00:00Z
    ## 5               PCR   METAGENOMIC metagenome ILLUMINA 2022-05-24T00:00:00Z
    ## 6               PCR   METAGENOMIC metagenome ILLUMINA 2022-05-24T00:00:00Z
    ## 7               PCR   METAGENOMIC metagenome ILLUMINA 2022-05-24T00:00:00Z
    ## 8               PCR   METAGENOMIC metagenome ILLUMINA 2022-05-24T00:00:00Z
    ## 9               PCR   METAGENOMIC metagenome ILLUMINA 2022-05-24T00:00:00Z
    ## 10              PCR   METAGENOMIC metagenome ILLUMINA 2022-05-24T00:00:00Z
    ##             create_date version Sample.Name SRA.Study
    ## 1  2022-05-24T03:18:00Z       1      RA116C SRP376757
    ## 2  2022-05-24T03:19:00Z       1      RA114C SRP376757
    ## 3  2022-05-24T03:18:00Z       1      RA113C SRP376757
    ## 4  2022-05-24T03:18:00Z       1      RA112C SRP376757
    ## 5  2022-05-24T03:18:00Z       1      RA111C SRP376757
    ## 6  2022-05-24T03:18:00Z       1      RA110C SRP376757
    ## 7  2022-05-24T03:18:00Z       1      RA109C SRP376757
    ## 8  2022-05-24T03:18:00Z       1      RA108C SRP376757
    ## 9  2022-05-24T03:18:00Z       1      RA107C SRP376757
    ## 10 2022-05-24T03:18:00Z       1      RA106C SRP376757

``` r
# IMPORTANT: Check the unique values in the column of interest or query column.
lapply(initial_meta[c('host_body_product')], table) # site
```

    ## $host_body_product
    ## 
    ## oral cavity 
    ##         189

``` r
lapply(initial_meta[c('Organism', 'Host')], table) # organism
```

    ## $Organism
    ## 
    ## metagenome 
    ##        189 
    ## 
    ## $Host
    ## 
    ## Homo sapiens 
    ##          189

``` r
table(initial_meta$Sample.Name) # subject id
```

    ## 
    ##  RA001 RA001C  RA002 RA002C  RA003 RA003C  RA004 RA004C  RA005 RA005C  RA006 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## RA006C  RA007 RA007C  RA008 RA008C  RA009 RA009C  RA010  RA011 RA011C  RA012 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  RA013 RA013C  RA014 RA014C  RA015 RA015C  RA016 RA016C  RA017 RA017C  RA018 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## RA018C  RA019 RA019C  RA020 RA020C  RA021 RA021C  RA022 RA022C  RA023 RA023C 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  RA024 RA024C  RA025 RA025C  RA026 RA026C RA027C  RA028 RA028C  RA029 RA029C 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  RA030 RA030C  RA031 RA031C  RA032 RA032C  RA033 RA033C  RA034 RA034C  RA035 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## RA035C  RA036 RA036C  RA037 RA037C RA038C RA039C  RA040 RA040C  RA041 RA041C 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  RA042 RA042C  RA043 RA043C  RA044 RA044C  RA045 RA045C  RA046 RA046C  RA047 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## RA047C  RA048 RA048C  RA049 RA049C RA050C  RA051 RA051C  RA052 RA052C  RA053 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## RA053C  RA054 RA054C  RA055 RA055C  RA056 RA056C  RA057 RA057C  RA058 RA058C 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  RA059 RA059C RA060C  RA061 RA061C  RA062 RA062C  RA063  RA064 RA064C  RA065 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## RA065C  RA066 RA066C  RA067 RA067C  RA068 RA068C  RA069 RA069C  RA070 RA070C 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  RA071 RA071C  RA072 RA072C  RA073 RA073C  RA074 RA074C  RA075 RA075C RA076C 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  RA077 RA077C  RA078 RA078C  RA079 RA079C  RA080 RA080C RA081C  RA082 RA083C 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## RA084C RA086C RA087C RA088C RA090C RA091C RA092C RA093C RA094C RA095C RA096C 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## RA097C RA098C RA099C RA100C RA101C RA102C RA103C RA104C RA105C RA106C RA107C 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## RA108C RA109C RA110C RA111C RA112C RA113C RA114C RA115C RA116C RA117C RA118C 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## RA119C RA120C 
    ##      1      1

``` r
#table(initial_meta$sequencing_meth..exp.) # output is pyrosequencing
table(initial_meta$Platform) # output is Illumina. according to the methods, Illumina is correct.
```

    ## 
    ## ILLUMINA 
    ##      189

``` r
# table(initial_meta$target_subfragment..exp.)
```

# Step 3: Removing Control and Non-human Samples

First, we need to remove control samples(including both negative
controls and positive controls) and non-human samples.

``` r
# (1) EDIT HERE: Parameter setup
control_column <- "" # column name in initial_meta
# unique(initial_meta[[control_column]])
control_keywords <- c() # Enter exact keywords

species_column <- "" # column name in initial_meta
unique(initial_meta[[species_column]])
```

    ## NULL

``` r
non_human_keywords <- c() # Enter exact keywords

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

    ## Initial samples counts: 189

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

    ## Remaining samples after Step 2: 189

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
  sample_id = initial_meta$Run,
  subject_id = initial_meta$Sample.Name,
  age = initial_meta$Host_age,
  age_info = NA,
  bmi = initial_meta$host_body_mass_index,
  sex = NA,
  sequencing_region = NA, # extract from methods section
  platform = "Illumina NovaSeq 6000", # extract from methods section, NOT initial_meta
  country = initial_meta$geo_loc_name_country,
  raw_site = initial_meta$host_body_product,
  site = "pending", # modify later
  site_group = "pending", # modify later
  public_time = initial_meta$ReleaseDate,
  pcr_primers = NA, # Format: xxF-xxR
  title = "Oral microbial signature of rheumatoid arthritis in female patients", # extract from article
  disease_study = "yes", # determine according to the article
  raw_status = initial_meta$Sample.Name,  
  status = "pending", # modify later
  matched_samples = "pending", # modify later
  is_time_series = "no", # determine according to the article
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
# In this part, you need to customize the code based on the situation of the raw_ columns to obtain the content of the target columns. IMPORTANT: Add comments if the modification is not straightforward!!
metadata <- metadata %>%
  mutate(site = case_when(
    grepl("oral cavity", raw_site) ~ "oral scraping", # During the visit, scrapings were collected from the surface of the tongue, gums, tonsils, and plaque.
    TRUE ~ NA_character_
  )) # according to the 'host_body_product' column of initial_meta

metadata <- metadata %>%
  mutate(status = case_when(
    !grepl("C$", raw_status) ~ "RA", # rheumatoid arthritis (RA)
    grepl("C$", raw_status) ~ "HC", 
    TRUE ~ NA_character_
  )) # according to the 'Sample.Name' column of initial_meta
```

``` r
# (2) Filling the 'site_group' and 'matched_samples' columns
# editing site_group
unique(metadata$site)
```

    ## [1] "oral scraping"

``` r
# EDIT HERE: add oral sites
Oral_site = unique(metadata$site)
# EDIT HERE: add gut sites
Gut_site = c()
# just run
metadata <- metadata %>%
  mutate(site_group = case_when(
    site %in% Oral_site ~ "oral",
    site %in% Gut_site ~ "gut",
    TRUE ~ "other" 
  ))

# editing matched_samples, just run.
metadata$matched_samples <- ifelse(
  ave(metadata$subject_id, metadata$subject_id, FUN = length) > 1, 
  "yes", 
  "no"
)
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
# Insert Metadata_Summary.png into the current Rmd output
knitr::include_graphics("Metadata_Summary_v2.png")
```

<img src="Metadata_Summary_v2.png" alt="" width="1200" />

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
