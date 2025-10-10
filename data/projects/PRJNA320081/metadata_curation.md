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
setwd("/home/hby/oral_project/PRJNA320081")

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

    ##          Run Assay.Type AvgSpotLen     Bases  BioProject    BioSample              BioSampleModel    Bytes                    Center.Name Collection_Date Consent DATASTORE.filetype
    ## 1 SRR3467161        WXS        502  12745780 PRJNA320081 SAMN04915571 Metagenome or environmental  7314834 LOS ALAMOS NATIONAL LABORATORY            2013  public   fastq,run.zq,sra
    ## 2 SRR3467148        WXS        502  14751772 PRJNA320081 SAMN04915559 Metagenome or environmental  8615645 LOS ALAMOS NATIONAL LABORATORY            2013  public   fastq,run.zq,sra
    ## 3 SRR3467170        WXS        502  11121810 PRJNA320081 SAMN04915557 Metagenome or environmental  6486082 LOS ALAMOS NATIONAL LABORATORY            2013  public   fastq,run.zq,sra
    ## 4 SRR3467168        WXS        502  77706086 PRJNA320081 SAMN04915555 Metagenome or environmental 47020728 LOS ALAMOS NATIONAL LABORATORY            2013  public   fastq,run.zq,sra
    ## 5 SRR3467166        WXS        502 158260018 PRJNA320081 SAMN04915553 Metagenome or environmental 92623701 LOS ALAMOS NATIONAL LABORATORY            2013  public   fastq,run.zq,sra
    ## 6 SRR3467158        WXS        502  32249986 PRJNA320081 SAMN04915551 Metagenome or environmental 18606768 LOS ALAMOS NATIONAL LABORATORY            2013  public   fastq,run.zq,sra
    ##   DATASTORE.provider                     DATASTORE.region Experiment geo_loc_name_country geo_loc_name_country_continent geo_loc_name         HOST     Instrument    isolation_source
    ## 1         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX1736942              Bolivia                  South America      Bolivia Homo sapiens Illumina MiSeq Human infant saliva
    ## 2         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX1736929              Bolivia                  South America      Bolivia Homo sapiens Illumina MiSeq Human infant saliva
    ## 3         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX1736951              Bolivia                  South America      Bolivia Homo sapiens Illumina MiSeq Human infant saliva
    ## 4         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX1736949              Bolivia                  South America      Bolivia Homo sapiens Illumina MiSeq Human infant saliva
    ## 5         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX1736947              Bolivia                  South America      Bolivia Homo sapiens Illumina MiSeq Human infant saliva
    ## 6         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX1736939              Bolivia                  South America      Bolivia Homo sapiens Illumina MiSeq Human infant saliva
    ##   Library.Name LibraryLayout LibrarySelection LibrarySource                Organism Platform          ReleaseDate          create_date version Host_Age host_sex Sample_Name SRA.Study
    ## 1         BSA1        PAIRED           RANDOM   METAGENOMIC human saliva metagenome ILLUMINA 2016-05-04T00:00:00Z 2016-04-29T11:59:00Z       1      1.3        F        BSA1 SRP074186
    ## 2        BSA10        PAIRED           RANDOM   METAGENOMIC human saliva metagenome ILLUMINA 2016-05-04T00:00:00Z 2016-04-29T11:57:00Z       1      1.7        M       BSA10 SRP074186
    ## 3        BSA11        PAIRED           RANDOM   METAGENOMIC human saliva metagenome ILLUMINA 2016-05-04T00:00:00Z 2016-04-29T11:58:00Z       1      1.1        M       BSA11 SRP074186
    ## 4        BSA12        PAIRED           RANDOM   METAGENOMIC human saliva metagenome ILLUMINA 2016-05-04T00:00:00Z 2016-04-29T11:59:00Z       1      0.8        M       BSA12 SRP074186
    ## 5        BSA13        PAIRED           RANDOM   METAGENOMIC human saliva metagenome ILLUMINA 2016-05-04T00:00:00Z 2016-04-29T12:00:00Z       1      1.2        M       BSA13 SRP074186
    ## 6        BSA14        PAIRED           RANDOM   METAGENOMIC human saliva metagenome ILLUMINA 2016-05-04T00:00:00Z 2016-04-29T11:58:00Z       1      1.1        M       BSA14 SRP074186
    ##   ref_biomaterial
    ## 1               1
    ## 2              10
    ## 3              11
    ## 4              12
    ## 5              13
    ## 6              14

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
  sample_name = initial_meta$Sample_Name,
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  public_time = initial_meta$ReleaseDate,
  Sex = initial_meta$host_sex,
  Age = initial_meta$Host_Age,
  BMI = NA
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$Site <- "saliva"
metadata$Status <- "H"
metadata$Region <- "V4"
metadata$PCR_primers <- "515F-806R"
metadata$title <- "Salivary microbiomes of indigenous Tsimane mothers and infants are distinct despite frequent premastication"
metadata$Disease_study <- "NO"
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

    ##    BioProject        Run sample_name Age BMI Sex Region       Platform Country   Site          public_time PCR_primers
    ## 1 PRJNA320081 SRR3467161        BSA1 1.3  NA   F     V4 Illumina MiSeq Bolivia saliva 2016-05-04T00:00:00Z   515F-806R
    ## 2 PRJNA320081 SRR3467148       BSA10 1.7  NA   M     V4 Illumina MiSeq Bolivia saliva 2016-05-04T00:00:00Z   515F-806R
    ## 3 PRJNA320081 SRR3467170       BSA11 1.1  NA   M     V4 Illumina MiSeq Bolivia saliva 2016-05-04T00:00:00Z   515F-806R
    ## 4 PRJNA320081 SRR3467168       BSA12 0.8  NA   M     V4 Illumina MiSeq Bolivia saliva 2016-05-04T00:00:00Z   515F-806R
    ## 5 PRJNA320081 SRR3467166       BSA13 1.2  NA   M     V4 Illumina MiSeq Bolivia saliva 2016-05-04T00:00:00Z   515F-806R
    ## 6 PRJNA320081 SRR3467158       BSA14 1.1  NA   M     V4 Illumina MiSeq Bolivia saliva 2016-05-04T00:00:00Z   515F-806R
    ##                                                                                                         title Disease_study Status matched_samples
    ## 1 Salivary microbiomes of indigenous Tsimane mothers and infants are distinct despite frequent premastication            NO      H              NO
    ## 2 Salivary microbiomes of indigenous Tsimane mothers and infants are distinct despite frequent premastication            NO      H              NO
    ## 3 Salivary microbiomes of indigenous Tsimane mothers and infants are distinct despite frequent premastication            NO      H              NO
    ## 4 Salivary microbiomes of indigenous Tsimane mothers and infants are distinct despite frequent premastication            NO      H              NO
    ## 5 Salivary microbiomes of indigenous Tsimane mothers and infants are distinct despite frequent premastication            NO      H              NO
    ## 6 Salivary microbiomes of indigenous Tsimane mothers and infants are distinct despite frequent premastication            NO      H              NO

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
