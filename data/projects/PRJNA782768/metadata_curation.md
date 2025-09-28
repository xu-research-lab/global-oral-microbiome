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

# Set working directory (update this path to your local directory)
setwd("D:/oral_project/PRJNA782768")

# Read the initial metadata files
initial_meta <- read.csv("SraRunTable.csv", header = TRUE)
initial_meta1 <- read.csv("SraRunTable1.csv", header = TRUE)

# Combine the two metadata files
meta <- rbind(initial_meta, initial_meta1)

# Remove entries where Organism is "gut metagenome" (non-oral samples)
#meta <- subset(meta, Organism != "gut metagenome")

# Check the number of samples containing "HT" in Sample.Name (e.g., healthy controls)
#nrow(subset(meta, grepl("HT", Sample.Name)))

# Display the first few rows of the metadata
head(meta)
```

    ##           Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample
    ## 1 SRR17012663   AMPLICON        495 39118365 PRJNA782768 SAMN23394823
    ## 2 SRR17012664   AMPLICON        493 48684243 PRJNA782768 SAMN23394825
    ## 3 SRR17012665   AMPLICON        492 42691824 PRJNA782768 SAMN23394822
    ## 4 SRR17012666   AMPLICON        493 40515726 PRJNA782768 SAMN23394821
    ## 5 SRR17012667   AMPLICON        494 49553634 PRJNA782768 SAMN23394820
    ## 6 SRR17012668   AMPLICON        495 52647210 PRJNA782768 SAMN23394819
    ##                BioSampleModel    Bytes              Center.Name Collection_Date
    ## 1 Metagenome or environmental 13124801 SCHOOL OF BASIC MEDICINE      2021-01-06
    ## 2 Metagenome or environmental 16341806 SCHOOL OF BASIC MEDICINE      2021-01-06
    ## 3 Metagenome or environmental 14354884 SCHOOL OF BASIC MEDICINE      2021-01-06
    ## 4 Metagenome or environmental 13825252 SCHOOL OF BASIC MEDICINE      2021-01-06
    ## 5 Metagenome or environmental 16530367 SCHOOL OF BASIC MEDICINE      2021-01-06
    ## 6 Metagenome or environmental 17690855 SCHOOL OF BASIC MEDICINE      2021-01-06
    ##   Consent DATASTORE.filetype DATASTORE.provider
    ## 1  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 2  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 3  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 4  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 6  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                       DATASTORE.region  Experiment geo_loc_name_country
    ## 1 gs.us-east1,ncbi.public,s3.us-east-1 SRX13202890                China
    ## 2 gs.us-east1,ncbi.public,s3.us-east-1 SRX13202883                China
    ## 3 gs.us-east1,ncbi.public,s3.us-east-1 SRX13202889                China
    ## 4 gs.us-east1,ncbi.public,s3.us-east-1 SRX13202888                China
    ## 5 gs.us-east1,ncbi.public,s3.us-east-1 SRX13202887                China
    ## 6 gs.us-east1,ncbi.public,s3.us-east-1 SRX13202886                China
    ##   geo_loc_name_country_continent   geo_loc_name    HOST     Instrument
    ## 1                           Asia China:Shanghai missing Illumina MiSeq
    ## 2                           Asia China:Shanghai missing Illumina MiSeq
    ## 3                           Asia China:Shanghai missing Illumina MiSeq
    ## 4                           Asia China:Shanghai missing Illumina MiSeq
    ## 5                           Asia China:Shanghai missing Illumina MiSeq
    ## 6                           Asia China:Shanghai missing Illumina MiSeq
    ##   isolation_source          lat_lon Library.Name LibraryLayout LibrarySelection
    ## 1          missing 31.19 N 121.60 E          F20        PAIRED              PCR
    ## 2          missing 31.19 N 121.60 E          F22        PAIRED              PCR
    ## 3          missing 31.19 N 121.60 E          F19        PAIRED              PCR
    ## 4          missing 31.19 N 121.60 E          F18        PAIRED              PCR
    ## 5          missing 31.19 N 121.60 E          F17        PAIRED              PCR
    ## 6          missing 31.19 N 121.60 E          F16        PAIRED              PCR
    ##   LibrarySource       Organism Platform          ReleaseDate
    ## 1   METAGENOMIC gut metagenome ILLUMINA 2022-03-19T00:00:00Z
    ## 2   METAGENOMIC gut metagenome ILLUMINA 2022-03-19T00:00:00Z
    ## 3   METAGENOMIC gut metagenome ILLUMINA 2022-03-19T00:00:00Z
    ## 4   METAGENOMIC gut metagenome ILLUMINA 2022-03-19T00:00:00Z
    ## 5   METAGENOMIC gut metagenome ILLUMINA 2022-03-19T00:00:00Z
    ## 6   METAGENOMIC gut metagenome ILLUMINA 2022-03-19T00:00:00Z
    ##            create_date version Sample.Name SRA.Study
    ## 1 2021-11-23T02:11:00Z       1         F20 SRP347352
    ## 2 2021-11-23T02:11:00Z       1         F22 SRP347352
    ## 3 2021-11-23T02:11:00Z       1         F19 SRP347352
    ## 4 2021-11-23T02:11:00Z       1         F18 SRP347352
    ## 5 2021-11-23T02:11:00Z       1         F17 SRP347352
    ## 6 2021-11-23T02:11:00Z       1         F16 SRP347352

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

``` r
# Create a uniform metadata dataframe
metadata <- data.frame(
  BioProject = meta$BioProject,
  Run = meta$Run,
  sample_name = meta$Sample.Name,
  Age = NA,  # To be filled Update based on you metadata
  BMI = NA,  # To be filled based on you metadata
  Sex = NA,  # To be filled based on you metadata
  Region = "V3-V4",  # Update if different for your study
  Platform = meta$Instrument,
  Country = meta$geo_loc_name_country,
  Site = meta$isolation_source, # other names of this columns such as env_median, env_material
  public_time = meta$ReleaseDate,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$PCR_primers <- "338F-806R"  # Update based on the publication
metadata$title <- "Distribution characteristics of oral microbiota and its relationship with intestinal microbiota in patients with type 2 diabetes mellitus"
metadata$Disease_study <- "YES"

# Assign disease status based on sample_name patterns
metadata <- metadata %>%
  mutate(Status = case_when(
    grepl("DT", sample_name) ~ "T2DM",
    grepl("HT|Tb", sample_name) ~ "HC"
  ))
#get subjects id according to sample id
metadata<-metadata%>%
  mutate(Site=case_when(grepl("F",sample_name)~"Stool",grepl("T",sample_name)~"Tongue coating"))

metadata$sample_name<-gsub("DT|DF|T|F","",metadata$sample_name)
#get subjects id with Saliva available
keep_subjects<-subset(metadata,Site %in% "Tongue coating")
#remove subjects without Saliva samples
metadata<-subset(metadata,sample_name %in% keep_subjects$sample_name)
#get subjects sampling counts
sample_counts <- table(metadata$sample_name)
# Create matched_samples column, any subjects have >=2 samples will be label "YES", otherwise "NO"
metadata$matched_samples <- ifelse(metadata$sample_name %in% names(sample_counts[sample_counts >= 2]), "YES", "NO")
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 426"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 298"

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
  "PCR_primers", "title", "Disease_study", "Status"
)]

# Export the uniform metadata to a CSV file
write.csv(metadata, file = "metadata.csv", quote = FALSE)
head(metadata)
```

    ##     BioProject         Run sample_name Age BMI Sex Region       Platform
    ## 11 PRJNA782768 SRR17012673         11a  NA  NA  NA  V3-V4 Illumina MiSeq
    ## 12 PRJNA782768 SRR17012674         02a  NA  NA  NA  V3-V4 Illumina MiSeq
    ## 13 PRJNA782768 SRR17012675         01a  NA  NA  NA  V3-V4 Illumina MiSeq
    ## 14 PRJNA782768 SRR17012676         b20  NA  NA  NA  V3-V4 Illumina MiSeq
    ## 15 PRJNA782768 SRR17012677         b19  NA  NA  NA  V3-V4 Illumina MiSeq
    ## 16 PRJNA782768 SRR17012678         b18  NA  NA  NA  V3-V4 Illumina MiSeq
    ##    Country           Site          public_time PCR_primers
    ## 11   China Tongue coating 2022-03-19T00:00:00Z   338F-806R
    ## 12   China Tongue coating 2022-03-19T00:00:00Z   338F-806R
    ## 13   China Tongue coating 2022-03-19T00:00:00Z   338F-806R
    ## 14   China Tongue coating 2022-03-19T00:00:00Z   338F-806R
    ## 15   China Tongue coating 2022-03-19T00:00:00Z   338F-806R
    ## 16   China Tongue coating 2022-03-19T00:00:00Z   338F-806R
    ##                                                                                                                                        title
    ## 11 Distribution characteristics of oral microbiota and its relationship with intestinal microbiota in patients with type 2 diabetes mellitus
    ## 12 Distribution characteristics of oral microbiota and its relationship with intestinal microbiota in patients with type 2 diabetes mellitus
    ## 13 Distribution characteristics of oral microbiota and its relationship with intestinal microbiota in patients with type 2 diabetes mellitus
    ## 14 Distribution characteristics of oral microbiota and its relationship with intestinal microbiota in patients with type 2 diabetes mellitus
    ## 15 Distribution characteristics of oral microbiota and its relationship with intestinal microbiota in patients with type 2 diabetes mellitus
    ## 16 Distribution characteristics of oral microbiota and its relationship with intestinal microbiota in patients with type 2 diabetes mellitus
    ##    Disease_study Status
    ## 11           YES   <NA>
    ## 12           YES   <NA>
    ## 13           YES   <NA>
    ## 14           YES     HC
    ## 15           YES     HC
    ## 16           YES     HC

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
