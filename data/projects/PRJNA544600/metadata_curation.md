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
setwd("/home/hby/oral_project/PRJNA544600")

# Read the initial metadata files
initial_meta <- read.csv("SraRunTable.csv", header = T)

# Combine the two metadata files
# meta <- rbind(initial_meta, initial_meta1)

# Remove entries non-oral samples
# initial_meta <- subset(initial_meta,isolate %in% c("saliva","tongue coat"))

# Check the number of samples containing "HT" in Sample.Name (e.g., healthy controls)
# nrow(subset(meta, grepl("HT", Sample.Name)))

# Display the first few rows of the metadata
head(initial_meta)
```

    ##          Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample              BioSampleModel   Bytes                                        Center.Name Collection_Date Consent
    ## 1 SRR9165946   AMPLICON        418 14345622 PRJNA544600 SAMN11843730 Metagenome or environmental 8007310 PEKING UNIVERSITY SCHOOL AND HOSPITAL OF STOMATOLO      2009/10/10  public
    ## 2 SRR9165947   AMPLICON        412 13885131 PRJNA544600 SAMN11843730 Metagenome or environmental 8106886 PEKING UNIVERSITY SCHOOL AND HOSPITAL OF STOMATOLO      2009/10/10  public
    ## 3 SRR9165948   AMPLICON        421 14530605 PRJNA544600 SAMN11843730 Metagenome or environmental 8184295 PEKING UNIVERSITY SCHOOL AND HOSPITAL OF STOMATOLO      2009/10/10  public
    ## 4 SRR9165949   AMPLICON        423 14351305 PRJNA544600 SAMN11843730 Metagenome or environmental 7860433 PEKING UNIVERSITY SCHOOL AND HOSPITAL OF STOMATOLO      2009/10/10  public
    ## 5 SRR9165950   AMPLICON        423 15955295 PRJNA544600 SAMN11843730 Metagenome or environmental 8822233 PEKING UNIVERSITY SCHOOL AND HOSPITAL OF STOMATOLO      2009/10/10  public
    ## 6 SRR9165951   AMPLICON        421 14817774 PRJNA544600 SAMN11843730 Metagenome or environmental 8173689 PEKING UNIVERSITY SCHOOL AND HOSPITAL OF STOMATOLO      2009/10/10  public
    ##   DATASTORE.filetype DATASTORE.provider                     DATASTORE.region Experiment geo_loc_name_country geo_loc_name_country_continent  geo_loc_name         HOST     Instrument
    ## 1   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX5938738                China                           Asia China:Beijing Homo sapiens Illumina MiSeq
    ## 2   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX5938737                China                           Asia China:Beijing Homo sapiens Illumina MiSeq
    ## 3   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX5938736                China                           Asia China:Beijing Homo sapiens Illumina MiSeq
    ## 4   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX5938735                China                           Asia China:Beijing Homo sapiens Illumina MiSeq
    ## 5   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX5938734                China                           Asia China:Beijing Homo sapiens Illumina MiSeq
    ## 6   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX5938733                China                           Asia China:Beijing Homo sapiens Illumina MiSeq
    ##   isolation_source          lat_lon state                 Library_Name LibraryLayout LibrarySelection LibrarySource              Organism Platform          ReleaseDate
    ## 1           saliva 39.26 N 115.27 E  OSAS salivary microbiome in osas8        SINGLE              PCR       GENOMIC human oral metagenome ILLUMINA 2019-06-03T00:00:00Z
    ## 2           saliva 39.26 N 115.27 E  OSAS salivary microbiome in osas9        SINGLE              PCR       GENOMIC human oral metagenome ILLUMINA 2019-06-03T00:00:00Z
    ## 3           saliva 39.26 N 115.27 E  OSAS salivary microbiome in osas6        SINGLE              PCR       GENOMIC human oral metagenome ILLUMINA 2019-06-03T00:00:00Z
    ## 4           saliva 39.26 N 115.27 E  OSAS salivary microbiome in osas7        SINGLE              PCR       GENOMIC human oral metagenome ILLUMINA 2019-06-03T00:00:00Z
    ## 5           saliva 39.26 N 115.27 E  OSAS salivary microbiome in osas4        SINGLE              PCR       GENOMIC human oral metagenome ILLUMINA 2019-06-03T00:00:00Z
    ## 6           saliva 39.26 N 115.27 E  OSAS salivary microbiome in osas5        SINGLE              PCR       GENOMIC human oral metagenome ILLUMINA 2019-06-03T00:00:00Z
    ##            create_date version                                                      Sample.Name SRA.Study
    ## 1 2019-05-31T12:52:00Z       1 salivary microbiome in Obstructive Sleep Apnea syndrome patients SRP200013
    ## 2 2019-05-31T12:52:00Z       1 salivary microbiome in Obstructive Sleep Apnea syndrome patients SRP200013
    ## 3 2019-05-31T12:52:00Z       1 salivary microbiome in Obstructive Sleep Apnea syndrome patients SRP200013
    ## 4 2019-05-31T12:52:00Z       1 salivary microbiome in Obstructive Sleep Apnea syndrome patients SRP200013
    ## 5 2019-05-31T12:52:00Z       1 salivary microbiome in Obstructive Sleep Apnea syndrome patients SRP200013
    ## 6 2019-05-31T12:52:00Z       1 salivary microbiome in Obstructive Sleep Apnea syndrome patients SRP200013

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
  sample_name = initial_meta$Library_Name,
  Platform = initial_meta$Instrument,
  Status = initial_meta$state,
  Site = initial_meta$isolation_source,
  public_time = initial_meta$ReleaseDate,
  Country = initial_meta$geo_loc_name_country,
  Sex = NA,
  Age = NA,
  BMI = NA
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$Region <- "V3-V4"
metadata$PCR_primers <- "341F-806R"
metadata$title <- "Analysis of the Salivary Microbiome in Obstructive Sleep Apnea Syndrome Patients"
metadata$Disease_study <- "YES"
metadata$matched_samples <- "NO"
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 24"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 24"

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

    ##    BioProject        Run                  sample_name Age BMI Sex Region       Platform Country   Site          public_time PCR_primers
    ## 1 PRJNA544600 SRR9165946 salivary microbiome in osas8  NA  NA  NA  V3-V4 Illumina MiSeq   China saliva 2019-06-03T00:00:00Z   341F-806R
    ## 2 PRJNA544600 SRR9165947 salivary microbiome in osas9  NA  NA  NA  V3-V4 Illumina MiSeq   China saliva 2019-06-03T00:00:00Z   341F-806R
    ## 3 PRJNA544600 SRR9165948 salivary microbiome in osas6  NA  NA  NA  V3-V4 Illumina MiSeq   China saliva 2019-06-03T00:00:00Z   341F-806R
    ## 4 PRJNA544600 SRR9165949 salivary microbiome in osas7  NA  NA  NA  V3-V4 Illumina MiSeq   China saliva 2019-06-03T00:00:00Z   341F-806R
    ## 5 PRJNA544600 SRR9165950 salivary microbiome in osas4  NA  NA  NA  V3-V4 Illumina MiSeq   China saliva 2019-06-03T00:00:00Z   341F-806R
    ## 6 PRJNA544600 SRR9165951 salivary microbiome in osas5  NA  NA  NA  V3-V4 Illumina MiSeq   China saliva 2019-06-03T00:00:00Z   341F-806R
    ##                                                                              title Disease_study Status matched_samples
    ## 1 Analysis of the Salivary Microbiome in Obstructive Sleep Apnea Syndrome Patients           YES   OSAS              NO
    ## 2 Analysis of the Salivary Microbiome in Obstructive Sleep Apnea Syndrome Patients           YES   OSAS              NO
    ## 3 Analysis of the Salivary Microbiome in Obstructive Sleep Apnea Syndrome Patients           YES   OSAS              NO
    ## 4 Analysis of the Salivary Microbiome in Obstructive Sleep Apnea Syndrome Patients           YES   OSAS              NO
    ## 5 Analysis of the Salivary Microbiome in Obstructive Sleep Apnea Syndrome Patients           YES   OSAS              NO
    ## 6 Analysis of the Salivary Microbiome in Obstructive Sleep Apnea Syndrome Patients           YES   OSAS              NO

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
