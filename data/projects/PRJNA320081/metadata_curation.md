R Notebook: Metadata Processing Workflow
================

# Step 1: Obtain Initial Metadata from Each Project

To obtain the initial metadata: 1. Search for the BioProject ID (e.g.,
“PRJEBXXXX”) on the [NCBI SRA
website](https://www.ncbi.nlm.nih.gov/Traces/study/). 2. In the “Select”
table, click on “Metadata” and download the file. 3. The downloaded file
is typically named “SraRunTable.csv” and contains metadata for the
project.

``` r
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)

# Set working directory (update this path to your local directory)
setwd("E:/bioproject_YL/PRJNA320081")

# Read the initial metadata files
initial_meta <- read.csv("SraRunTable.csv", header = T) 

# Display the first few rows of the metadata
head(initial_meta)
```

    ##          Run Assay.Type AvgSpotLen     Bases  BioProject    BioSample
    ## 1 SRR3467146        WXS        502  23621610 PRJNA320081 SAMN04915549
    ## 2 SRR3467147        WXS        502 159826760 PRJNA320081 SAMN04915550
    ## 3 SRR3467148        WXS        502  14751772 PRJNA320081 SAMN04915559
    ## 4 SRR3467149        WXS        502  59857476 PRJNA320081 SAMN04915560
    ## 5 SRR3467150        WXS        502  21340522 PRJNA320081 SAMN04915561
    ## 6 SRR3467151        WXS        502 157291158 PRJNA320081 SAMN04915562
    ##                BioSampleModel    Bytes                    Center.Name
    ## 1 Metagenome or environmental 14134612 LOS ALAMOS NATIONAL LABORATORY
    ## 2 Metagenome or environmental 94362692 LOS ALAMOS NATIONAL LABORATORY
    ## 3 Metagenome or environmental  8615645 LOS ALAMOS NATIONAL LABORATORY
    ## 4 Metagenome or environmental 35422290 LOS ALAMOS NATIONAL LABORATORY
    ## 5 Metagenome or environmental 12574210 LOS ALAMOS NATIONAL LABORATORY
    ## 6 Metagenome or environmental 91432625 LOS ALAMOS NATIONAL LABORATORY
    ##   Collection_Date Consent DATASTORE.filetype DATASTORE.provider
    ## 1            2013  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 2            2013  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 3            2013  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 4            2013  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5            2013  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 6            2013  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                       DATASTORE.region Experiment geo_loc_name_country
    ## 1 gs.us-east1,ncbi.public,s3.us-east-1 SRX1736927              Bolivia
    ## 2 gs.us-east1,ncbi.public,s3.us-east-1 SRX1736928              Bolivia
    ## 3 gs.us-east1,ncbi.public,s3.us-east-1 SRX1736929              Bolivia
    ## 4 gs.us-east1,ncbi.public,s3.us-east-1 SRX1736930              Bolivia
    ## 5 gs.us-east1,ncbi.public,s3.us-east-1 SRX1736931              Bolivia
    ## 6 gs.us-east1,ncbi.public,s3.us-east-1 SRX1736932              Bolivia
    ##   geo_loc_name_country_continent geo_loc_name         HOST     Instrument
    ## 1                  South America      Bolivia Homo sapiens Illumina MiSeq
    ## 2                  South America      Bolivia Homo sapiens Illumina MiSeq
    ## 3                  South America      Bolivia Homo sapiens Illumina MiSeq
    ## 4                  South America      Bolivia Homo sapiens Illumina MiSeq
    ## 5                  South America      Bolivia Homo sapiens Illumina MiSeq
    ## 6                  South America      Bolivia Homo sapiens Illumina MiSeq
    ##      isolation_source Library.Name LibraryLayout LibrarySelection LibrarySource
    ## 1 Human infant saliva        BSA15        PAIRED           RANDOM   METAGENOMIC
    ## 2  Human adult saliva        MSA15        PAIRED           RANDOM   METAGENOMIC
    ## 3 Human infant saliva        BSA10        PAIRED           RANDOM   METAGENOMIC
    ## 4  Human adult saliva        MSA10        PAIRED           RANDOM   METAGENOMIC
    ## 5 Human infant saliva         BSA9        PAIRED           RANDOM   METAGENOMIC
    ## 6  Human adult saliva         MSA9        PAIRED           RANDOM   METAGENOMIC
    ##                  Organism Platform          ReleaseDate          create_date
    ## 1 human saliva metagenome ILLUMINA 2016-05-04T00:00:00Z 2016-04-29T11:57:00Z
    ## 2 human saliva metagenome ILLUMINA 2016-05-04T00:00:00Z 2016-04-29T11:57:00Z
    ## 3 human saliva metagenome ILLUMINA 2016-05-04T00:00:00Z 2016-04-29T11:57:00Z
    ## 4 human saliva metagenome ILLUMINA 2016-05-04T00:00:00Z 2016-04-29T11:57:00Z
    ## 5 human saliva metagenome ILLUMINA 2016-05-04T00:00:00Z 2016-04-29T11:58:00Z
    ## 6 human saliva metagenome ILLUMINA 2016-05-04T00:00:00Z 2016-04-29T11:58:00Z
    ##   version Sample.Name SRA.Study ref_biomaterial
    ## 1       1       BSA15 SRP074186              15
    ## 2       1       MSA15 SRP074186              15
    ## 3       1       BSA10 SRP074186              10
    ## 4       1       MSA10 SRP074186              10
    ## 5       1        BSA9 SRP074186               9
    ## 6       1        MSA9 SRP074186               9

# Step 2: Create a Uniform Metadata Template

The uniform metadata should include the following columns:

1.  **BioProject**: Project ID (e.g., “PRJNA339212”)
2.  **Run**: Sample ID (e.g., “ERRXXXXXX” or “SRRXXXXXX”)
3.  **sample_name**: Original sample name submitted by the authors (not
    the database-assigned ID)
4.  **Age**: Age of the subject (if available)
5.  **BMI**: Body Mass Index (if available)
6.  **Sex**: Sex of the subject (“F” or “M”)
7.  **Region**: Targeted 16S rRNA region (e.g., “V3-V4”)
8.  **Platform**: Sequencing platform (e.g., “Illumina MiSeq”)
9.  **Country**: Country of sample origin (e.g., “China”)
10. **Site**: Sample site (e.g., “saliva”, “dental plaque”)
11. **public_time**: Data release date (format: “YYYY-MM-DD”)
12. **PCR_primers**: PCR primers used (e.g., “515F-806R”)
13. **title**: Title of the associated publication
14. **Disease_study**: Whether the study involves a disease (“YES” or
    “NO”)
15. **Status**: Disease status (e.g., “T2DM”, “HC”). Note: If the study
    includes both disease and healthy samples, this must be specified
    per sample.
16. **matched_samples**: matched samples (“YES” or “NO”)

``` r
# Create a uniform metadata dataframe
metadata <- data.frame(
  BioProject = initial_meta$BioProject,
  Run = initial_meta$Run,
  sample_name = initial_meta$Sample.Name,
  BMI = NA,
  Sex = initial_meta$Sample.Name,
  Region = "V4",
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  Site = initial_meta$Sample.Name,
  public_time = initial_meta$ReleaseDate,
  stringsAsFactors = FALSE
)

# Assign Site based on Sample.Name patterns
metadata <- metadata %>%
  mutate(Site = case_when(
    # Match BSA series and MSA1-MSA11
    str_detect(Site, "^BSA\\d+$") ~ "buccal swab",
    str_detect(Site, "^MSA(1[0-1][1-9])$") ~ "buccal swab",
    # Match MSA12-MSA15 and FD series
    str_detect(Site, "^MSA1[2-5]$") ~ "expectoration",
    str_detect(Site, "^FD\\d+$") ~ "expectoration",
    TRUE ~ Site
  ))

# Assign Sex based on Sample.Name patterns
metadata <- metadata %>%
  mutate(Sex = ifelse(grepl("BSA1[0-5]|BSA[48]", Sex), "M", "F"))

# Set the corresponding Age value for each sample_name value
age_map <- c(
  "BSA15" = 0.9, "MSA15" = 25.2, "BSA10" = 1.7, "MSA10" = 24.6,
  "BSA9" = 1.7, "MSA9" = 22.7, "BSA8" = 2.0, "MSA8" = 29.6,
  "BSA4" = 1.5, "MSA4" = 21.3, "BSA3" = 0.8, "MSA3" = 36.3,
  "BSA14" = 1.1, "BSA2" = 1.1, "MSA2" = 18.8, "BSA1" = 1.3,
  "MSA1" = 20.1, "RFD" = 24.9, "MFD" = 31.6, "MSA14" = 31.6,
  "BSA13" = 1.2, "MSA13" = 27.7, "BSA12" = 0.8, "MSA12" = 24.9,
  "BSA11" = 1.1, "MSA11" = 26.8
)
metadata$Age <- age_map[metadata$sample_name]
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$PCR_primers <- "515F-806R"
metadata$title <- "Salivary microbiomes of indigenous tsimane mothers and infants are distinct despite frequent premastication"
metadata$Disease_study <- "NO"
# Assign disease status based on sample_name patterns
metadata$Status <- "HC"

# Check initial data distribution
cat("Initial sample distribution:\n")
```

    ## Initial sample distribution:

``` r
table(metadata$sample_name)
```

    ## 
    ##  BSA1 BSA10 BSA11 BSA12 BSA13 BSA14 BSA15  BSA2  BSA3  BSA4  BSA8  BSA9   MFD 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##  MSA1 MSA10 MSA11 MSA12 MSA13 MSA14 MSA15  MSA2  MSA3  MSA4  MSA8  MSA9   RFD 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1

``` r
table(metadata$Site)
```

    ## 
    ##   buccal swab expectoration           MFD          MSA1         MSA10 
    ##            12             4             1             1             1 
    ##         MSA11          MSA2          MSA3          MSA4          MSA8 
    ##             1             1             1             1             1 
    ##          MSA9           RFD 
    ##             1             1

``` r
# Calculate sampling frequency per subject
sample_counts <- table(metadata$sample_name)
cat("\nNumber of samples per subject:\n")
```

    ## 
    ## Number of samples per subject:

``` r
print(sample_counts)
```

    ## 
    ##  BSA1 BSA10 BSA11 BSA12 BSA13 BSA14 BSA15  BSA2  BSA3  BSA4  BSA8  BSA9   MFD 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##  MSA1 MSA10 MSA11 MSA12 MSA13 MSA14 MSA15  MSA2  MSA3  MSA4  MSA8  MSA9   RFD 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1

``` r
# Create matched_samples column
metadata$matched_samples <- ifelse(metadata$sample_name %in% names(sample_counts[sample_counts >= 2]), "YES", "NO")

cat("\nMatched samples distribution:\n")
```

    ## 
    ## Matched samples distribution:

``` r
table(metadata$matched_samples)
```

    ## 
    ## NO 
    ## 26

``` r
# Identify samples for removal: single samples excluding buccal swab and expectoration
single_samples <- names(sample_counts[sample_counts == 1])
samples_to_remove <- metadata$sample_name %in% single_samples & 
                     !(metadata$Site %in% c("buccal swab", "expectoration"))

cat("\nNumber of samples to be removed:", sum(samples_to_remove), "\n")
```

    ## 
    ## Number of samples to be removed: 10

``` r
# Remove identified samples
metadata_clean <- metadata[!samples_to_remove, ]

cat("\nFinal sample distribution after cleaning:\n")
```

    ## 
    ## Final sample distribution after cleaning:

``` r
table(metadata_clean$sample_name)
```

    ## 
    ##  BSA1 BSA10 BSA11 BSA12 BSA13 BSA14 BSA15  BSA2  BSA3  BSA4  BSA8  BSA9 MSA12 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## MSA13 MSA14 MSA15 
    ##     1     1     1

``` r
table(metadata_clean$Site)
```

    ## 
    ##   buccal swab expectoration 
    ##            12             4

``` r
table(metadata_clean$matched_samples)
```

    ## 
    ## NO 
    ## 16

``` r
# Use cleaned dataset for subsequent analysis
metadata <- metadata_clean
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 16"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 16"

# Step 5: Visualize the Metadata

``` r
knitr::include_graphics("Metadata_Summary.png")
```

<img src="Metadata_Summary.png" width="1200" />

# Step 6: Finalize and Export the Metadata

``` r
# Reorder columns to match the desired format
metadata <- metadata[, c(
  "BioProject", "Run", "sample_name", "Age", "BMI", "Sex", 
  "Region", "Platform", "Country", "Site", "public_time", 
  "PCR_primers", "title", "Disease_study", "Status", "matched_samples"
)]

# Export the uniform metadata to a CSV file
write.csv(metadata, file = "metadata.csv", quote = FALSE)
head(metadata)
```

    ##    BioProject        Run sample_name  Age BMI Sex Region       Platform Country
    ## 1 PRJNA320081 SRR3467146       BSA15  0.9  NA   M     V4 Illumina MiSeq Bolivia
    ## 2 PRJNA320081 SRR3467147       MSA15 25.2  NA   F     V4 Illumina MiSeq Bolivia
    ## 3 PRJNA320081 SRR3467148       BSA10  1.7  NA   M     V4 Illumina MiSeq Bolivia
    ## 5 PRJNA320081 SRR3467150        BSA9  1.7  NA   F     V4 Illumina MiSeq Bolivia
    ## 7 PRJNA320081 SRR3467152        BSA8  2.0  NA   M     V4 Illumina MiSeq Bolivia
    ## 9 PRJNA320081 SRR3467154        BSA4  1.5  NA   M     V4 Illumina MiSeq Bolivia
    ##            Site          public_time PCR_primers
    ## 1   buccal swab 2016-05-04T00:00:00Z   515F-806R
    ## 2 expectoration 2016-05-04T00:00:00Z   515F-806R
    ## 3   buccal swab 2016-05-04T00:00:00Z   515F-806R
    ## 5   buccal swab 2016-05-04T00:00:00Z   515F-806R
    ## 7   buccal swab 2016-05-04T00:00:00Z   515F-806R
    ## 9   buccal swab 2016-05-04T00:00:00Z   515F-806R
    ##                                                                                                         title
    ## 1 Salivary microbiomes of indigenous tsimane mothers and infants are distinct despite frequent premastication
    ## 2 Salivary microbiomes of indigenous tsimane mothers and infants are distinct despite frequent premastication
    ## 3 Salivary microbiomes of indigenous tsimane mothers and infants are distinct despite frequent premastication
    ## 5 Salivary microbiomes of indigenous tsimane mothers and infants are distinct despite frequent premastication
    ## 7 Salivary microbiomes of indigenous tsimane mothers and infants are distinct despite frequent premastication
    ## 9 Salivary microbiomes of indigenous tsimane mothers and infants are distinct despite frequent premastication
    ##   Disease_study Status matched_samples
    ## 1            NO     HC              NO
    ## 2            NO     HC              NO
    ## 3            NO     HC              NO
    ## 5            NO     HC              NO
    ## 7            NO     HC              NO
    ## 9            NO     HC              NO

## Notes for Students:

1.  **Metadata Accuracy**: Always cross-check the metadata with the
    original publication to ensure accuracy.
2.  **Manual Entries**: Columns like `PCR_primers`, and `Status` often
    require manual entry based on the publication.
3.  **Disease Status**: Be careful when assigning `Status`—some studies
    include both disease and healthy samples.
4.  **Visualization**: The plots help quickly identify the distribution
    of values in each column, which is useful for quality control.
5.  **Reproducibility**: Update the working directory and file paths to
    match your local environment.
