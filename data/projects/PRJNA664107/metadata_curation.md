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
project_id <- "PRJNA664107" 

base_dir <- "E:/oral/26/metadata-2604/Oral-Metadata-YL-02"
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

    ##            Run Assay.Type AvgSpotLen      Bases  BioProject    BioSample
    ## 1  SRR12660203   AMPLICON        500   32815500 PRJNA664107 SAMN16198071
    ## 2  SRR12660204   AMPLICON        500   16811500 PRJNA664107 SAMN16198070
    ## 3  SRR12660205   AMPLICON        500   25061000 PRJNA664107 SAMN16198069
    ## 4  SRR12660206   AMPLICON        500   29310500 PRJNA664107 SAMN16198068
    ## 5  SRR12660207   AMPLICON        500   12234000 PRJNA664107 SAMN16198067
    ## 6  SRR12660208        WGS        500 5483142500 PRJNA664107 SAMN16198110
    ## 7  SRR12660209        WGS        500 5881458500 PRJNA664107 SAMN16198109
    ## 8  SRR12660210        WGS        500 1894430000 PRJNA664107 SAMN16198108
    ## 9  SRR12660211   AMPLICON        500   16388000 PRJNA664107 SAMN16198107
    ## 10 SRR12660212   AMPLICON        500   22247000 PRJNA664107 SAMN16198106
    ##                          BioSampleModel      Bytes
    ## 1  MIMS.me,MIGS/MIMS/MIMARKS.human-oral   19348214
    ## 2  MIMS.me,MIGS/MIMS/MIMARKS.human-oral    9953828
    ## 3  MIMS.me,MIGS/MIMS/MIMARKS.human-oral   15716106
    ## 4  MIMS.me,MIGS/MIMS/MIMARKS.human-oral   17636552
    ## 5  MIMS.me,MIGS/MIMS/MIMARKS.human-oral    7135309
    ## 6  MIMS.me,MIGS/MIMS/MIMARKS.human-oral 3208546902
    ## 7  MIMS.me,MIGS/MIMS/MIMARKS.human-oral 3335032999
    ## 8  MIMS.me,MIGS/MIMS/MIMARKS.human-oral 1089474618
    ## 9  MIMS.me,MIGS/MIMS/MIMARKS.human-oral    9715961
    ## 10 MIMS.me,MIGS/MIMS/MIMARKS.human-oral   13796139
    ##                                                           Center.Name
    ## 1  FEDERAL RESEARCH AND CLINICAL CENTRE OF PHYSICAL-CHEMICAL MEDICINE
    ## 2  FEDERAL RESEARCH AND CLINICAL CENTRE OF PHYSICAL-CHEMICAL MEDICINE
    ## 3  FEDERAL RESEARCH AND CLINICAL CENTRE OF PHYSICAL-CHEMICAL MEDICINE
    ## 4  FEDERAL RESEARCH AND CLINICAL CENTRE OF PHYSICAL-CHEMICAL MEDICINE
    ## 5  FEDERAL RESEARCH AND CLINICAL CENTRE OF PHYSICAL-CHEMICAL MEDICINE
    ## 6  FEDERAL RESEARCH AND CLINICAL CENTRE OF PHYSICAL-CHEMICAL MEDICINE
    ## 7  FEDERAL RESEARCH AND CLINICAL CENTRE OF PHYSICAL-CHEMICAL MEDICINE
    ## 8  FEDERAL RESEARCH AND CLINICAL CENTRE OF PHYSICAL-CHEMICAL MEDICINE
    ## 9  FEDERAL RESEARCH AND CLINICAL CENTRE OF PHYSICAL-CHEMICAL MEDICINE
    ## 10 FEDERAL RESEARCH AND CLINICAL CENTRE OF PHYSICAL-CHEMICAL MEDICINE
    ##    Collection_Date Consent DATASTORE.filetype DATASTORE.provider
    ## 1       2018-01-10  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 2       2018-01-09  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 3       2018-01-08  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 4       2018-01-07  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5       2018-01-06  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 6       2018-02-19  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 7       2018-02-18  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 8       2018-02-17  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 9       2018-02-16  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 10      2018-02-15  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                        DATASTORE.region env_broad_scale env_local_scale
    ## 1  gs.us-east1,ncbi.public,s3.us-east-1 oral microbiome gingival sulcus
    ## 2  gs.us-east1,ncbi.public,s3.us-east-1 oral microbiome gingival sulcus
    ## 3  gs.us-east1,ncbi.public,s3.us-east-1 oral microbiome gingival sulcus
    ## 4  gs.us-east1,ncbi.public,s3.us-east-1 oral microbiome gingival sulcus
    ## 5  gs.us-east1,ncbi.public,s3.us-east-1 oral microbiome gingival sulcus
    ## 6  gs.us-east1,ncbi.public,s3.us-east-1 oral microbiome gingival sulcus
    ## 7  gs.us-east1,ncbi.public,s3.us-east-1 oral microbiome gingival sulcus
    ## 8  gs.us-east1,ncbi.public,s3.us-east-1 oral microbiome gingival sulcus
    ## 9  gs.us-east1,ncbi.public,s3.us-east-1 oral microbiome gingival sulcus
    ## 10 gs.us-east1,ncbi.public,s3.us-east-1 oral microbiome gingival sulcus
    ##       env_medium Experiment geo_loc_name_country geo_loc_name_country_continent
    ## 1  ENVO:00010483 SRX9141303               Russia                         Europe
    ## 2  ENVO:00010483 SRX9141302               Russia                         Europe
    ## 3  ENVO:00010483 SRX9141301               Russia                         Europe
    ## 4  ENVO:00010483 SRX9141300               Russia                         Europe
    ## 5  ENVO:00010483 SRX9141299               Russia                         Europe
    ## 6  ENVO:00010483 SRX9141298               Russia                         Europe
    ## 7  ENVO:00010483 SRX9141297               Russia                         Europe
    ## 8  ENVO:00010483 SRX9141296               Russia                         Europe
    ## 9  ENVO:00010483 SRX9141295               Russia                         Europe
    ## 10 ENVO:00010483 SRX9141294               Russia                         Europe
    ##      geo_loc_name         HOST          Instrument isolation_source
    ## 1  Russia: Moscow Homo sapiens Illumina HiSeq 2500  gingival sulcus
    ## 2  Russia: Moscow Homo sapiens Illumina HiSeq 2500  gingival sulcus
    ## 3  Russia: Moscow Homo sapiens Illumina HiSeq 2500  gingival sulcus
    ## 4  Russia: Moscow Homo sapiens Illumina HiSeq 2500  gingival sulcus
    ## 5  Russia: Moscow Homo sapiens Illumina HiSeq 2500  gingival sulcus
    ## 6  Russia: Moscow Homo sapiens Illumina HiSeq 2500  gingival sulcus
    ## 7  Russia: Moscow Homo sapiens Illumina HiSeq 2500  gingival sulcus
    ## 8  Russia: Moscow Homo sapiens Illumina HiSeq 2500  gingival sulcus
    ## 9  Russia: Moscow Homo sapiens Illumina HiSeq 2500  gingival sulcus
    ## 10 Russia: Moscow Homo sapiens Illumina HiSeq 2500  gingival sulcus
    ##                    lat_lon Library.Name LibraryLayout LibrarySelection
    ## 1  55.751244 N 37.618423 E          13C        PAIRED              PCR
    ## 2  55.751244 N 37.618423 E        12SDP        PAIRED              PCR
    ## 3  55.751244 N 37.618423 E          12P        PAIRED              PCR
    ## 4  55.751244 N 37.618423 E          12C        PAIRED              PCR
    ## 5  55.751244 N 37.618423 E        11SDP        PAIRED              PCR
    ## 6  55.751244 N 37.618423 E       CPT2DM        PAIRED              PCR
    ## 7  55.751244 N 37.618423 E           CP        PAIRED              PCR
    ## 8  55.751244 N 37.618423 E      Control        PAIRED              PCR
    ## 9  55.751244 N 37.618423 E         9SDP        PAIRED              PCR
    ## 10 55.751244 N 37.618423 E           9P        PAIRED              PCR
    ##    LibrarySource              Organism Platform          ReleaseDate
    ## 1    METAGENOMIC human oral metagenome ILLUMINA 2020-09-17T00:00:00Z
    ## 2    METAGENOMIC human oral metagenome ILLUMINA 2020-09-17T00:00:00Z
    ## 3    METAGENOMIC human oral metagenome ILLUMINA 2020-09-17T00:00:00Z
    ## 4    METAGENOMIC human oral metagenome ILLUMINA 2020-09-17T00:00:00Z
    ## 5    METAGENOMIC human oral metagenome ILLUMINA 2020-09-17T00:00:00Z
    ## 6    METAGENOMIC human oral metagenome ILLUMINA 2020-09-17T00:00:00Z
    ## 7    METAGENOMIC human oral metagenome ILLUMINA 2020-09-17T00:00:00Z
    ## 8    METAGENOMIC human oral metagenome ILLUMINA 2020-09-17T00:00:00Z
    ## 9    METAGENOMIC human oral metagenome ILLUMINA 2020-09-17T00:00:00Z
    ## 10   METAGENOMIC human oral metagenome ILLUMINA 2020-09-17T00:00:00Z
    ##             create_date version Sample.Name SRA.Study
    ## 1  2020-09-17T12:10:00Z       1     13C_S43 SRP282880
    ## 2  2020-09-17T12:09:00Z       1   12SDP_S27 SRP282880
    ## 3  2020-09-17T12:10:00Z       1     12P_S12 SRP282880
    ## 4  2020-09-17T12:10:00Z       1     12C_S42 SRP282880
    ## 5  2020-09-17T12:10:00Z       1   11SDP_S26 SRP282880
    ## 6  2020-09-17T12:18:00Z       1  CPT2DM_S21 SRP282880
    ## 7  2020-09-17T12:16:00Z       1      CP_S14 SRP282880
    ## 8  2020-09-17T12:12:00Z       1 Control_S12 SRP282880
    ## 9  2020-09-17T12:10:00Z       1    9SDP_S24 SRP282880
    ## 10 2020-09-17T12:10:00Z       1       9P_S9 SRP282880

``` r
# IMPORTANT: Check the unique values in the column of interest or query column.
lapply(initial_meta[c('isolation_source')], table) # site
```

    ## $isolation_source
    ## 
    ## gingival sulcus 
    ##              49

``` r
lapply(initial_meta[c('env_broad_scale', 'HOST', 'env_local_scale', 'Organism')], table) # organism
```

    ## $env_broad_scale
    ## 
    ## oral microbiome 
    ##              49 
    ## 
    ## $HOST
    ## 
    ## Homo sapiens 
    ##           49 
    ## 
    ## $env_local_scale
    ## 
    ## gingival sulcus 
    ##              49 
    ## 
    ## $Organism
    ## 
    ## human oral metagenome 
    ##                    49

``` r
table(initial_meta$Library.Name) # subject id
```

    ## 
    ##     10C     10P   10SDP     11C     11P   11SDP     12C     12P   12SDP     13C 
    ##       1       1       1       1       1       1       1       1       1       1 
    ##     13P   13SDP     14C     14P   14SDP     15C     15P   15SDP     16C      1C 
    ##       1       1       1       1       1       1       1       1       1       1 
    ##      1P    1SDP      2C      2P    2SDP      3C      3P    3SDP      4C      4P 
    ##       1       1       1       1       1       1       1       1       1       1 
    ##    4SDP      5C      5P    5SDP      6C      6P    6SDP      7C      7P    7SDP 
    ##       1       1       1       1       1       1       1       1       1       1 
    ##      8C      8P    8SDP      9C      9P    9SDP Control      CP  CPT2DM 
    ##       1       1       1       1       1       1       1       1       1

``` r
#table(initial_meta$sequencing_meth..exp.) # output is Sequencing by synthesis
table(initial_meta$Platform) # output is Illumina. according to the methods,  Illumina is correct.
```

    ## 
    ## ILLUMINA 
    ##       49

``` r
#table(initial_meta$target_subfragment..exp.)
```

# Step 3: Removing Control and Non-human Samples

First, we need to remove control samples(including both negative
controls and positive controls) and non-human samples.

``` r
# (1) EDIT HERE: Parameter setup
control_column <- "Assay.Type" # column name in initial_meta
# unique(initial_meta[[control_column]])
control_keywords <- c("WGS") # Enter exact keywords

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

    ## Initial samples counts: 49

``` r
cat("Removed control samples:", sum(is_control), "\n")
```

    ## Removed control samples: 3

``` r
cat("Removed non-human samples:", sum(is_non_human), "\n")
```

    ## Removed non-human samples: 0

``` r
cat("Total samples removed ", sum(samples_to_remove), "\n")
```

    ## Total samples removed  3

``` r
cat("Remaining samples after Step 2:", nrow(initial_meta), "\n")
```

    ## Remaining samples after Step 2: 46

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
  subject_id = initial_meta$Library.Name,
  age = NA,
  age_info = "53.9±7.3", # extract from abstract section
  bmi = NA,
  sex = NA,
  sequencing_region = "V3-V4", # extract from methods section
  platform = "Illumina HiSeq 2500", # extract from methods section, NOT initial_meta
  country = initial_meta$geo_loc_name_country,
  raw_site = initial_meta$Library.Name,
  site = "pending", # modify later
  site_group = "pending", # modify later
  public_time = initial_meta$ReleaseDate,
  pcr_primers = NA, # Format: xxF-xxR
  title = "Drift of the subgingival periodontal microbiome during chronic periodontitis in type 2 diabetes mellitus patients", # extract from article
  disease_study = "yes", # determine according to the article
  raw_status = initial_meta$Library.Name,
  status = "pending", # modify later 
  matched_samples = "pending", # modify later
  is_time_series = "no", 
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
    grepl("C", raw_site) ~ "gingival sulcus", # healthy subjects (Control)
    grepl("P|SDP", raw_site) ~ "periodontal pocket", # chronic periodontitis (CP), chronic periodontitis associated with type 2 diabetes mellitus (CPT2DM)
    TRUE ~ NA_character_
  )) # according to the 'Library.Name' column of initial_meta, and the article: The contents of the periodontal pocket in patients with chronic periodontitis (CP), chronic periodontitis associated with type 2 diabetes mellitus (CPT2DM) and the contents of the gingival sulcus in control subjects were the study material.

metadata <- metadata %>%
  mutate(status = case_when(
    grepl("SDP$", raw_status) ~ "chronic periodontitis associated with type 2 diabetes mellitus",
    grepl("C$", raw_status) ~ "HC",
    grepl("P$", raw_status) ~ "chronic periodontitis",
    TRUE ~ NA_character_
  )) # according to the 'Library.Name' column of initial_meta, and the article: We examined 15 patients with chronic periodontitis (CP), 15 patients with chronic periodontitis associated with type 2 diabetes mellitus (CPT2DM), and 16 healthy subjects (Control).
```

``` r
# (2) Filling the 'site_group' and 'matched_samples' columns
# editing site_group
unique(metadata$site)
```

    ## [1] "gingival sulcus"    "periodontal pocket"

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
