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
project_id <- "PRJNA731066" 

base_dir <- "E:/oral/26/metadata-2604/Oral-Metadata-YL-03"
setwd(file.path(base_dir, project_id))

# (2) Loading raw SRA metadata file
sra_file <- "SraRunTable.csv"
initial_meta <- read.csv(sra_file, header = TRUE)

# (3) OPTIONAL: If extra metadata (e.g., Table.x or supplemental materials) are available from the article, you may need to merge the file with initial_meta manually. Below is example code for merging. Depending on the specific case of each study, you can customize the code as needed.
supp_table <- read.csv("supplemental_table_v2.csv", header = TRUE)
initial_meta <- left_join(initial_meta, supp_table, by = c("Sample.Name" = "Patient"))
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
    ## 1  SRR14585963   AMPLICON        501 72704748 PRJNA731066 SAMN19267559
    ## 2  SRR14585964   AMPLICON        501 43002635 PRJNA731066 SAMN19267558
    ## 3  SRR14585965   AMPLICON        501 19478039 PRJNA731066 SAMN19267557
    ## 4  SRR14585966   AMPLICON        501  5614168 PRJNA731066 SAMN19267556
    ## 5  SRR14585967   AMPLICON        501 21256063 PRJNA731066 SAMN19267555
    ## 6  SRR14585968   AMPLICON        499  2215535 PRJNA731066 SAMN19267554
    ## 7  SRR14585969   AMPLICON        501  9027672 PRJNA731066 SAMN19267553
    ## 8  SRR14585970   AMPLICON        501 21513199 PRJNA731066 SAMN19267552
    ## 9  SRR14585971   AMPLICON        501 25573924 PRJNA731066 SAMN19267569
    ## 10 SRR14585972   AMPLICON        501 29959580 PRJNA731066 SAMN19267568
    ##                 BioSampleModel    Bytes                      Center.Name
    ## 1  Metagenome or environmental 51448407 JUSTUS LIEBIG UNIVERSITY GIESSEN
    ## 2  Metagenome or environmental 30190798 JUSTUS LIEBIG UNIVERSITY GIESSEN
    ## 3  Metagenome or environmental 13530303 JUSTUS LIEBIG UNIVERSITY GIESSEN
    ## 4  Metagenome or environmental  3973501 JUSTUS LIEBIG UNIVERSITY GIESSEN
    ## 5  Metagenome or environmental 15016443 JUSTUS LIEBIG UNIVERSITY GIESSEN
    ## 6  Metagenome or environmental  1643653 JUSTUS LIEBIG UNIVERSITY GIESSEN
    ## 7  Metagenome or environmental  6429327 JUSTUS LIEBIG UNIVERSITY GIESSEN
    ## 8  Metagenome or environmental 15202395 JUSTUS LIEBIG UNIVERSITY GIESSEN
    ## 9  Metagenome or environmental 17805237 JUSTUS LIEBIG UNIVERSITY GIESSEN
    ## 10 Metagenome or environmental 20957774 JUSTUS LIEBIG UNIVERSITY GIESSEN
    ##    Collection_Date Consent DATASTORE.filetype DATASTORE.provider
    ## 1             2019  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 2             2019  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 3             2019  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 4             2019  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5             2019  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 6             2019  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 7             2019  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 8             2019  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 9             2019  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 10            2019  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                        DATASTORE.region  Experiment geo_loc_name_country
    ## 1  gs.us-east1,ncbi.public,s3.us-east-1 SRX10933206              Germany
    ## 2  gs.us-east1,ncbi.public,s3.us-east-1 SRX10933205              Germany
    ## 3  gs.us-east1,ncbi.public,s3.us-east-1 SRX10933204              Germany
    ## 4  gs.us-east1,ncbi.public,s3.us-east-1 SRX10933203              Germany
    ## 5  gs.us-east1,ncbi.public,s3.us-east-1 SRX10933202              Germany
    ## 6  gs.us-east1,ncbi.public,s3.us-east-1 SRX10933201              Germany
    ## 7  gs.us-east1,ncbi.public,s3.us-east-1 SRX10933200              Germany
    ## 8  gs.us-east1,ncbi.public,s3.us-east-1 SRX10933199              Germany
    ## 9  gs.us-east1,ncbi.public,s3.us-east-1 SRX10933198              Germany
    ## 10 gs.us-east1,ncbi.public,s3.us-east-1 SRX10933197              Germany
    ##    geo_loc_name_country_continent geo_loc_name         HOST host_subject
    ## 1                          Europe      Germany Homo sapiens        CHD10
    ## 2                          Europe      Germany Homo sapiens        CHD09
    ## 3                          Europe      Germany Homo sapiens        CHD08
    ## 4                          Europe      Germany Homo sapiens        CHD07
    ## 5                          Europe      Germany Homo sapiens        CHD06
    ## 6                          Europe      Germany Homo sapiens        CHD05
    ## 7                          Europe      Germany Homo sapiens        CHD04
    ## 8                          Europe      Germany Homo sapiens        CHD03
    ## 9                          Europe      Germany Homo sapiens        HCG10
    ## 10                         Europe      Germany Homo sapiens        HCG09
    ##        Instrument isolation_source        lat_lon Library.Name LibraryLayout
    ## 1  Illumina MiSeq carious material not applicable        CHD10        PAIRED
    ## 2  Illumina MiSeq carious material not applicable        CHD09        PAIRED
    ## 3  Illumina MiSeq carious material not applicable        CHD08        PAIRED
    ## 4  Illumina MiSeq carious material not applicable        CHD07        PAIRED
    ## 5  Illumina MiSeq carious material not applicable        CHD06        PAIRED
    ## 6  Illumina MiSeq carious material not applicable        CHD05        PAIRED
    ## 7  Illumina MiSeq carious material not applicable        CHD04        PAIRED
    ## 8  Illumina MiSeq carious material not applicable        CHD03        PAIRED
    ## 9  Illumina MiSeq carious material not applicable        HCG10        PAIRED
    ## 10 Illumina MiSeq carious material not applicable        HCG09        PAIRED
    ##    LibrarySelection LibrarySource              Organism Platform
    ## 1               PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 2               PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 3               PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 4               PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 5               PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 6               PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 7               PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 8               PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 9               PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 10              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ##             ReleaseDate          create_date version Sample.Name SRA.Study Sex
    ## 1  2021-09-06T00:00:00Z 2021-05-19T13:28:00Z       1       CHD10 SRP320499   m
    ## 2  2021-09-06T00:00:00Z 2021-05-19T13:28:00Z       1       CHD09 SRP320499   m
    ## 3  2021-09-06T00:00:00Z 2021-05-19T13:28:00Z       1       CHD08 SRP320499   f
    ## 4  2021-09-06T00:00:00Z 2021-05-19T13:28:00Z       1       CHD07 SRP320499   m
    ## 5  2021-09-06T00:00:00Z 2021-05-19T13:28:00Z       1       CHD06 SRP320499   f
    ## 6  2021-09-06T00:00:00Z 2021-05-19T13:28:00Z       1       CHD05 SRP320499   f
    ## 7  2021-09-06T00:00:00Z 2021-05-19T13:28:00Z       1       CHD04 SRP320499   m
    ## 8  2021-09-06T00:00:00Z 2021-05-19T13:28:00Z       1       CHD03 SRP320499   m
    ## 9  2021-09-06T00:00:00Z 2021-05-19T13:28:00Z       1       HCG10 SRP320499   f
    ## 10 2021-09-06T00:00:00Z 2021-05-19T13:28:00Z       1       HCG09 SRP320499   f
    ##    Age..in.months. dmf.t fistula.abscess
    ## 1               61     5               -
    ## 2               62    10               -
    ## 3               62    11               -
    ## 4               36     6               -
    ## 5               36     9               +
    ## 6               60     5               -
    ## 7               24    18               -
    ## 8               67    12               -
    ## 9               62     6               +
    ## 10              62    11               +
    ##                                  CHD.Warnes.Classification.
    ## 1                                                   ASD (1)
    ## 2                                                   ASD (1)
    ## 3                    VSD, Truncus arteriosus \ncommunis (3)
    ## 4                                   Tetralogy of Fallot (3)
    ## 5                           Mital valve \ninsufficiency (3)
    ## 6                  Hypoplastic left heart syndrome/ VSD (3)
    ## 7                                     Pulmonary atresia (3)
    ## 8  Pulmonary \ninsufficiency/ \ninterrupted aortic arch (3)
    ## 9                                                      n.a.
    ## 10                                                     n.a.
    ##                                           Heart.surgery
    ## 1                                                     -
    ## 2                                                     -
    ## 3                    VSD-Closure, Aortic reconstruction
    ## 4                     Conduit replacement, VSD -Closure
    ## 5                                       LIMA-LAD Bypass
    ## 6              Glenn \nanastomoses, Norwood \noperation
    ## 7                                                     -
    ## 8  Fontan-Operation,\nComprehensive stage 2 \nprocedure
    ## 9                                                  n.a.
    ## 10                                                 n.a.
    ##                                             medication
    ## 1                                                    -
    ## 2                                 Acetylsalicylic acid
    ## 3                                                    -
    ## 4                              Beta-Blocker, Diuretics
    ## 5                                                    -
    ## 6                                 Acetylsalicylic acid
    ## 7                                            Diuretics
    ## 8  Acetylsalicylic acid, ACE \ninhibitor, Beta-Blocker
    ## 9                                                 n.a.
    ## 10                                                n.a.

``` r
# IMPORTANT: Check the unique values in the column of interest or query column.
lapply(initial_meta[c('isolation_source')], table) # site
```

    ## $isolation_source
    ## 
    ## carious material 
    ##               20

``` r
lapply(initial_meta[c('HOST', 'Organism')], table) # organism
```

    ## $HOST
    ## 
    ## Homo sapiens 
    ##           20 
    ## 
    ## $Organism
    ## 
    ## human oral metagenome 
    ##                    20

``` r
table(initial_meta$Sample.Name) # subject id
```

    ## 
    ## CHD01 CHD02 CHD03 CHD04 CHD05 CHD06 CHD07 CHD08 CHD09 CHD10 CHD11 HCG01 HCG02 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## HCG04 HCG05 HCG06 HCG07 HCG08 HCG09 HCG10 
    ##     1     1     1     1     1     1     1

``` r
#table(initial_meta$sequencing_meth..exp.) # output is Sequencing by synthesis
table(initial_meta$Platform) # output is Illumina. according to the methods,  Illumina is correct.
```

    ## 
    ## ILLUMINA 
    ##       20

``` r
#table(initial_meta$target_subfragment) # output is V4
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

    ## Initial samples counts: 20

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

    ## Remaining samples after Step 2: 20

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
  age = round(initial_meta$Age..in.months./12, 2),
  age_info = "2-6", # extract from methods section 
  bmi = NA,
  raw_sex = initial_meta$Sex,
  sex = "pending", # modify later
  sequencing_region = "V4", # extract from methods section
  platform = "Illumina MiSeq", # extract from methods section, NOT initial_meta
  country = initial_meta$geo_loc_name_country,
  site = initial_meta$isolation_source, # Article: A complete dental examination was performed after a supragingival polish, in which carious lesions, and missing and filled teeth (dmf-t) indexes were analyzed. During therapy, in the course of removal of caries and the treatment need for filling, the carious material was collected by a sterile excavator.
  site_group = "pending", # modify later
  public_time = initial_meta$ReleaseDate,
  pcr_primers = "515F-806R", # Format: xxF-xxR
  title = "Microbiome analysis of carious lesions in pre-school children with early childhood caries and congenital heart disease", # extract from article
  disease_study = "yes", # determine according to the article
  raw_status = initial_meta$Sample.Name, 
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
  mutate(sex = case_when(
    grepl("m", raw_sex) ~ "M",
    grepl("f", raw_sex) ~ "F",
    TRUE ~ NA_character_
  )) # according to the 'Sex' column of initial_meta

metadata <- metadata %>%
  mutate(status = case_when(
    grepl("HCG", raw_status) ~ "early childhood caries",
    grepl("CHD", raw_status) ~ "early childhood caries and congenital heart disease",
    TRUE ~ NA_character_
  )) # according to the 'Sample.Name' column of initial_meta
```

``` r
# (2) Filling the 'site_group' and 'matched_samples' columns
# editing site_group
unique(metadata$site)
```

    ## [1] "carious material"

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
