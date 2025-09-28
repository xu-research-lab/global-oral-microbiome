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
setwd("/home/hby/oral_project/PRJNA356414")

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

    ##          Run Assay.Type AvgSpotLen     Bases  BioProject    BioSample                                      BioSampleModel     Bytes            Center.Name Collection_Date Consent
    ## 1 SRR5080276   AMPLICON        600 313997611 PRJNA356414 SAMN06113483 MIMARKS.specimen,MIGS/MIMS/MIMARKS.human-associated 186722602 UNIVERSITY OF FLORENCE   not collected  public
    ## 2 SRR5080277   AMPLICON        600  24553001 PRJNA356414 SAMN06113492 MIMARKS.specimen,MIGS/MIMS/MIMARKS.human-associated  17101414 UNIVERSITY OF FLORENCE   not collected  public
    ## 3 SRR5080246   AMPLICON        599  47879676 PRJNA356414 SAMN06113491 MIMARKS.specimen,MIGS/MIMS/MIMARKS.human-associated  33164741 UNIVERSITY OF FLORENCE   not collected  public
    ## 4 SRR5080272   AMPLICON        601  75381173 PRJNA356414 SAMN06113501 MIMARKS.specimen,MIGS/MIMS/MIMARKS.human-associated  48889137 UNIVERSITY OF FLORENCE   not collected  public
    ## 5 SRR5080244   AMPLICON        600 658337866 PRJNA356414 SAMN06113482 MIMARKS.specimen,MIGS/MIMS/MIMARKS.human-associated 415062525 UNIVERSITY OF FLORENCE   not collected  public
    ## 6 SRR5080271   AMPLICON        600  67727468 PRJNA356414 SAMN06113493 MIMARKS.specimen,MIGS/MIMS/MIMARKS.human-associated  46684715 UNIVERSITY OF FLORENCE   not collected  public
    ##   DATASTORE.filetype DATASTORE.provider                     DATASTORE.region env_biome      env_feature env_material Experiment geo_loc_name_country geo_loc_name_country_continent
    ## 1   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1     human digestive system        stool SRX2398690                Italy                         Europe
    ## 2   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1     human      oral cavity       saliva SRX2398691                Italy                         Europe
    ## 3   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1     human digestive system        stool SRX2398660                Italy                         Europe
    ## 4   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1     human      oral cavity       saliva SRX2398686                Italy                         Europe
    ## 5   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1     human      oral cavity       saliva SRX2398658                Italy                         Europe
    ## 6   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1     human digestive system        stool SRX2398685                Italy                         Europe
    ##            geo_loc_name         HOST status host_subject_id   BMI     Instrument isol_growth_condt isolate                 lat_lon   Library.Name LibraryLayout LibrarySelection
    ## 1 Italy: Florence\\, FI Homo sapiens      H            CFP1 21.41 Illumina MiSeq    not applicable missing 43.803903 N 11.246223 E FAS-1 PILOT 02        PAIRED              PCR
    ## 2 Italy: Florence\\, FI Homo sapiens      H            CFP1 21.41 Illumina MiSeq    not applicable missing 43.803903 N 11.246223 E FAS-1 PILOT 11        PAIRED              PCR
    ## 3 Italy: Florence\\, FI Homo sapiens      H           CFP10 21.50 Illumina MiSeq    not applicable missing 43.803903 N 11.246223 E FAS-1 PILOT 10        PAIRED              PCR
    ## 4 Italy: Florence\\, FI Homo sapiens      H           CFP10 21.50 Illumina MiSeq    not applicable missing 43.803903 N 11.246223 E FAS-1 PILOT 20        PAIRED              PCR
    ## 5 Italy: Florence\\, FI Homo sapiens      H           CFP11 29.41 Illumina MiSeq    not applicable missing 43.803903 N 11.246223 E FAS-1 PILOT 01        PAIRED              PCR
    ## 6 Italy: Florence\\, FI Homo sapiens      H           CFP11 29.41 Illumina MiSeq    not applicable missing 43.803903 N 11.246223 E FAS-1 PILOT 12        PAIRED              PCR
    ##   LibrarySource         Organism Platform          ReleaseDate          create_date version     Sample_Name SRA.Study
    ## 1   METAGENOMIC human metagenome ILLUMINA 2016-12-11T00:00:00Z 2016-12-06T12:00:00Z       1 FAS_RFa11_pilot SRP094636
    ## 2   METAGENOMIC human metagenome ILLUMINA 2016-12-11T00:00:00Z 2016-12-06T11:59:00Z       1  FAS_RFa1_pilot SRP094636
    ## 3   METAGENOMIC human metagenome ILLUMINA 2016-12-11T00:00:00Z 2016-12-06T11:57:00Z       1 FAS_RFa19_pilot SRP094636
    ## 4   METAGENOMIC human metagenome ILLUMINA 2016-12-11T00:00:00Z 2016-12-06T11:59:00Z       1  FAS_RFa9_pilot SRP094636
    ## 5   METAGENOMIC human metagenome ILLUMINA 2016-12-11T00:00:00Z 2016-12-06T11:59:00Z       1 FAS_RFa10_pilot SRP094636
    ## 6   METAGENOMIC human metagenome ILLUMINA 2016-12-11T00:00:00Z 2016-12-06T11:59:00Z       1 FAS_RFa20_pilot SRP094636

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
  sample_name = initial_meta$host_subject_id,
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  Site = initial_meta$env_material,
  public_time = initial_meta$ReleaseDate,
  Status = initial_meta$status,
  BMI = initial_meta$BMI,
  stringsAsFactors = FALSE,
  Sex = NA,
  Age = NA
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$Region <- "V3-V4"
metadata$PCR_primers <- "515F-806R"
metadata$title <- "Prospective correlation between the patient microbiome with response to and development of immune-mediated adverse effects to immunotherapy in lung cancer"
metadata$Disease_study <- "YES"


# extract subject id by remove prefix of sample_name
metadata$subject_id <- metadata$sample_name

all_oral_sites <- c("saliva", "tongue coat", "buccal mucosa", "lower jaw mucosa", "lower labial mucosa", "upper jaw mucosa", "tongue", "buccal")
# Check initial data distribution
cat("Initial sample distribution:\n")
```

    ## Initial sample distribution:

``` r
table(metadata$subject_id)
```

    ## 
    ##  CFP1 CFP10 CFP11  CFP2  CFP3  CFP4  CFP6  CFP7  CFP8  CFP9  CM10  CM11  CM18  CM19  CM20  CM22  CM23  CM24   CM7   CM8 
    ##     2     2     2     2     2     2     2     2     2     2     3     3     3     3     3     3     3     3     3     3

``` r
table(metadata$Site)
```

    ## 
    ## biopsy saliva  stool 
    ##     10     20     20

``` r
subjects_with_oral_samples <- metadata %>%
  filter(Site %in% all_oral_sites) %>%
  distinct(subject_id) %>%
  pull(subject_id)

metadata_clean <- metadata %>%
  filter(subject_id %in% subjects_with_oral_samples)


sample_counts <- metadata_clean %>%
  group_by(subject_id) %>%
  summarise(n_samples = n(), .groups = "drop")

metadata_clean <- merge(metadata_clean, sample_counts, by = "subject_id")

# Create matched_samples column
metadata_clean$matched_samples <- ifelse(metadata_clean$n_samples >= 2, "YES", "NO")

cat("\nFinal sample distribution after cleaning:\n")
```

    ## 
    ## Final sample distribution after cleaning:

``` r
print(table(metadata_clean$Site))
```

    ## 
    ## biopsy saliva  stool 
    ##     10     20     20

``` r
cat("\nFinal sample distribution after cleaning:\n")
```

    ## 
    ## Final sample distribution after cleaning:

``` r
print(table(metadata_clean$matched_samples))
```

    ## 
    ## YES 
    ##  50

``` r
metadata_clean$subject_id <- NULL
metadata_clean$n_samples <- NULL

# Use cleaned dataset for subsequent analysis
metadata <- metadata_clean
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 50"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 20"

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

    ##    BioProject        Run sample_name Age   BMI Sex Region       Platform Country   Site          public_time PCR_primers
    ## 1 PRJNA356414 SRR5080276        CFP1  NA 21.41  NA  V3-V4 Illumina MiSeq   Italy  stool 2016-12-11T00:00:00Z   515F-806R
    ## 2 PRJNA356414 SRR5080277        CFP1  NA 21.41  NA  V3-V4 Illumina MiSeq   Italy saliva 2016-12-11T00:00:00Z   515F-806R
    ## 3 PRJNA356414 SRR5080246       CFP10  NA 21.50  NA  V3-V4 Illumina MiSeq   Italy  stool 2016-12-11T00:00:00Z   515F-806R
    ## 4 PRJNA356414 SRR5080272       CFP10  NA 21.50  NA  V3-V4 Illumina MiSeq   Italy saliva 2016-12-11T00:00:00Z   515F-806R
    ## 5 PRJNA356414 SRR5080244       CFP11  NA 29.41  NA  V3-V4 Illumina MiSeq   Italy saliva 2016-12-11T00:00:00Z   515F-806R
    ## 6 PRJNA356414 SRR5080271       CFP11  NA 29.41  NA  V3-V4 Illumina MiSeq   Italy  stool 2016-12-11T00:00:00Z   515F-806R
    ##                                                                                                                                                        title Disease_study Status
    ## 1 Prospective correlation between the patient microbiome with response to and development of immune-mediated adverse effects to immunotherapy in lung cancer           YES      H
    ## 2 Prospective correlation between the patient microbiome with response to and development of immune-mediated adverse effects to immunotherapy in lung cancer           YES      H
    ## 3 Prospective correlation between the patient microbiome with response to and development of immune-mediated adverse effects to immunotherapy in lung cancer           YES      H
    ## 4 Prospective correlation between the patient microbiome with response to and development of immune-mediated adverse effects to immunotherapy in lung cancer           YES      H
    ## 5 Prospective correlation between the patient microbiome with response to and development of immune-mediated adverse effects to immunotherapy in lung cancer           YES      H
    ## 6 Prospective correlation between the patient microbiome with response to and development of immune-mediated adverse effects to immunotherapy in lung cancer           YES      H
    ##   matched_samples
    ## 1             YES
    ## 2             YES
    ## 3             YES
    ## 4             YES
    ## 5             YES
    ## 6             YES

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
