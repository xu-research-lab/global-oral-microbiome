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

# =============================================================
# Curation Notes - PRJNA753008
# age          : From initial_meta column 'host_age_(days)', converted from days to years
# age_info     : "infant" — paper Table 1; subjects are newborns 0-180 days old
# sex          : From initial_meta 'host_sex'; "male"/"female" -> "M"/"F"
# subject_id   : From initial_meta 'Sample_code'; each code maps to one subject with paired oral+gut samples
# bmi          : NA — not available in CSV or paper
# platform     : Paper Methods section 2.2 — Illumina MiSeq
# region       : Paper Methods section 2.2 — V3-V4
# pcr_primers  : Paper Methods section 2.2 — 341F/805R (Illumina 16S Metagenomic Library prep)
# site         : From initial_meta 'env_medium'; "Stool" -> "stool", "saliva" -> "saliva"
# site_group   : saliva -> "oral"; Stool -> "gut"
# disease_study: "no" — paper title describes healthy infants cohort
# is_time_series: "no" — cross-sectional design; different babies enrolled at each time point
#                  (paper section 2.1: "we enrolled different babies at each time point")
# status       : "HC" — all subjects are healthy newborns per paper title and Table 1
# matched_samples: "yes" — 80 subjects have both oral and gut samples per Sample_code
# title        : Paper full title
# country      : From initial_meta 'geo_loc_name_country' — Italy
# public_time  : From initial_meta 'ReleaseDate'
# sample_id    : From initial_meta 'Run' column
# controls     : none found in metadata or paper for SRA samples
# non-human    : none — HOST column all "Homo sapiens"
# =============================================================

# (1) Setting working directory
# EDIT HERE: change project_id for each study
project_id <- "PRJNA753008"

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

# Paper: "A Parallel Tracking of Salivary and Gut Microbiota Profiles Can Reveal
# Maturation and Interplay of Early Life Microbial Communities in Healthy Infants"
# Microorganisms 2022, 10, 468. 116 healthy newborns, cross-sectional.
# Time points: 0, 7, 15, 30, 90, 180 days. Different babies at each time point.
# Sampling sites: saliva (collected by sterile cotton swab), stool/feces (from infant diapers).
# All subjects: healthy, breastfed, term delivery (37-42 weeks), birth weight 2.5-4.5 kg.
# Sequencing: Illumina MiSeq, V3-V4 region, primers per Illumina 16S Metagenomic protocol (341F/805R).

# (2) IMPORTANT: MUST go through each column of initial_meta
# head(initial_meta, 10)

# IMPORTANT: Check the unique values in the column of interest or query column.
# lapply(initial_meta[c('env_broad_scale', 'env_local_scale', 'env_medium')], table) # site: "saliva biome"/"fecal biome" / "saliva"/"Stool"
# table(initial_meta$HOST) # all "Homo sapiens" — no non-human filtering needed
# table(initial_meta$host_sex) # "male" / "female"
# table(initial_meta$host_age_(days)) # 7, 15, 30, 90, 180
# table(initial_meta$geo_loc_name_country) # "Italy"
# table(initial_meta$LibraryLayout) # all "PAIRED"
# table(initial_meta$platform..exp.) # all "ILLUMINA"
# table(initial_meta$sequencing_meth..exp.) # all "PCR"
```

# Step 3: Removing Control and Non-human Samples

First, we need to remove control samples(including both negative
controls and positive controls) and non-human samples.

``` r
# (1) EDIT HERE: Parameter setup
# No control samples in SRA metadata for this study; all are biological samples
control_column <- ""   # column name in initial_meta
# unique(initial_meta[[control_column]])
control_keywords <- c()     # Enter exact keywords

# All samples are Homo sapiens; no non-human filtering needed
species_column <- "HOST"          # column name in initial_meta
unique(initial_meta[[species_column]])
```

    ## [1] "Homo sapiens"

``` r
non_human_keywords <- c()            # Enter exact keywords


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

    ## Initial samples counts: 162

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

    ## Remaining samples after Step 2: 162

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
# Extract columns from initial_meta as follows: (1) No modification needed -> extract directly. (2) Modification needed -> create a new column with the prefix raw_ and mark the original column as "pending". (3) Need to find or determine based on the article.
metadata <- data.frame(
  bioproject = initial_meta$BioProject, # From initial_meta column 'BioProject'
  sample_id = initial_meta$Run, # From initial_meta column 'Run' (SRR IDs)
  subject_id = initial_meta[, "Sample_code"], # From initial_meta 'Sample_code'; same subject has both oral and gut samples
  raw_age = initial_meta[, "host_age_.days."], # From initial_meta column 'host_age_(days)'; days -> convert to years
  age = "pending",
  age_info = "infant", # From paper Table 1; all subjects are newborns 0-180 days old
  bmi = NA, # NA — not available in CSV or paper
  raw_sex = initial_meta$host_sex, # From initial_meta 'host_sex'; "male"/"female"
  sex = "pending",
  sequencing_region = "V3-V4", # From paper Methods section 2.2: "variable region V3-V4 of the 16S rRNA locus"
  platform = "Illumina MiSeq", # From paper Methods section 2.2: "sequenced on an Illumina MiSeqTM platform"
  country = initial_meta$geo_loc_name_country, # From initial_meta column 'geo_loc_name_country' — "Italy"
  raw_site = initial_meta$env_medium, # From initial_meta 'env_medium'; "Stool"/"saliva"
  site = "pending",
  site_group = "pending", # modify later
  public_time = initial_meta$ReleaseDate, # From initial_meta column 'ReleaseDate'
  pcr_primers = "341F-805R", # From paper Methods 2.2: Illumina 16S Metagenomic protocol, V3-V4 primers (341F/805R)
  title = "A Parallel Tracking of Salivary and Gut Microbiota Profiles Can Reveal Maturation and Interplay of Early Life Microbial Communities in Healthy Infants", # From paper title
  disease_study = "no", # From paper title and abstract — healthy infant cohort
  status = "HC", # From paper title and Table 1 — all subjects are healthy newborns, so healthy control
  matched_samples = "pending", # modify later
  is_time_series = "no", # From paper section 2.1: "we enrolled different babies at each time point" — cross-sectional, not longitudinal
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

# Convert age from days to years
metadata$age <- round(metadata$raw_age / 365, 4) # From initial_meta 'host_age_(days)'; days / 365 -> years (e.g., 7 days = 0.0192, 180 days = 0.4932)

metadata <- metadata %>%
  mutate(sex = case_when(
    raw_sex == "male" ~ "M",
    raw_sex == "female" ~ "F",
    TRUE ~ NA_character_
  )) # From initial_meta 'host_sex' column: "male" -> "M", "female" -> "F"

metadata$site <- tolower(metadata$raw_site) # From initial_meta 'env_medium': "Stool" -> "stool", "saliva" -> "saliva"
```

``` r
# (2) Filling the 'site_group' and 'matched_samples' columns
# editing site_group
unique(metadata$site)
```

    ## [1] "stool"  "saliva"

``` r
# EDIT HERE: add oral sites
# From paper Methods 2.1: saliva collected by swabbing infant cheeks; recorded as "saliva" in env_medium
Oral_site = c("saliva")
# EDIT HERE: add gut sites
# From paper Methods 2.1: stool/feces collected from infant diapers; recorded as "Stool" in env_medium
Gut_site = c("stool")
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
