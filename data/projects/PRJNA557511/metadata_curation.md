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
setwd("D:/oral_project/PRJNA557511")

# Read the initial metadata files
initial_meta<-read.csv("SraRunTable.csv",header = T)

# Combine the two metadata files
#meta <- rbind(initial_meta, initial_meta1)

#Remove entries where Organism is "gut metagenome" (non-oral samples)
#initial_meta <- subset(initial_meta,isolation_source!= "feces")

# Check the number of samples containing "HT" in Sample.Name (e.g., healthy controls)
#nrow(subset(meta, grepl("HT", Sample.Name)))

# Display the first few rows of the metadata
head(initial_meta)
```

    ##           Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample
    ## 1 SRR10188199   AMPLICON        602 37266810 PRJNA557511 SAMN12826803
    ## 2 SRR10188200   AMPLICON        602 28186844 PRJNA557511 SAMN12826802
    ## 3 SRR10188201   AMPLICON        602 34899746 PRJNA557511 SAMN12826801
    ## 4 SRR10188202   AMPLICON        602 55844530 PRJNA557511 SAMN12826800
    ## 5 SRR10188203   AMPLICON        602 31148082 PRJNA557511 SAMN12826799
    ## 6 SRR10188204   AMPLICON        602 35811776 PRJNA557511 SAMN12826798
    ##                BioSampleModel    Bytes
    ## 1 Metagenome or environmental 20236127
    ## 2 Metagenome or environmental 15594036
    ## 3 Metagenome or environmental 19169390
    ## 4 Metagenome or environmental 30563896
    ## 5 Metagenome or environmental 17239379
    ## 6 Metagenome or environmental 19661010
    ##                                             Center.Name Collection_Date Consent
    ## 1 THE FIRST AFFILIATED HOSPITAL OF ZHENGZHOU UNIVERSITY         2019-03  public
    ## 2 THE FIRST AFFILIATED HOSPITAL OF ZHENGZHOU UNIVERSITY         2019-03  public
    ## 3 THE FIRST AFFILIATED HOSPITAL OF ZHENGZHOU UNIVERSITY         2019-03  public
    ## 4 THE FIRST AFFILIATED HOSPITAL OF ZHENGZHOU UNIVERSITY         2019-03  public
    ## 5 THE FIRST AFFILIATED HOSPITAL OF ZHENGZHOU UNIVERSITY         2019-03  public
    ## 6 THE FIRST AFFILIATED HOSPITAL OF ZHENGZHOU UNIVERSITY         2019-03  public
    ##   DATASTORE.filetype DATASTORE.provider                     DATASTORE.region
    ## 1   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 2   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 3   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 4   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 5   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 6   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ##   Experiment geo_loc_name_country geo_loc_name_country_continent
    ## 1 SRX6908792                China                           Asia
    ## 2 SRX6908791                China                           Asia
    ## 3 SRX6908790                China                           Asia
    ## 4 SRX6908789                China                           Asia
    ## 5 SRX6908788                China                           Asia
    ## 6 SRX6908787                China                           Asia
    ##      geo_loc_name      HOST     Instrument isolation_source          lat_lon
    ## 1 China:Zhengzhou Sample100 Illumina MiSeq           saliva 34.44 N 113.42 E
    ## 2 China:Zhengzhou Sample099 Illumina MiSeq           saliva 34.44 N 113.42 E
    ## 3 China:Zhengzhou Sample098 Illumina MiSeq           saliva 34.44 N 113.42 E
    ## 4 China:Zhengzhou Sample097 Illumina MiSeq           saliva 34.44 N 113.42 E
    ## 5 China:Zhengzhou Sample096 Illumina MiSeq           saliva 34.44 N 113.42 E
    ## 6 China:Zhengzhou Sample095 Illumina MiSeq           saliva 34.44 N 113.42 E
    ##   Library.Name LibraryLayout LibrarySelection LibrarySource
    ## 1       HC__31        PAIRED              PCR         OTHER
    ## 2       HC__30        PAIRED              PCR         OTHER
    ## 3       HC__29        PAIRED              PCR         OTHER
    ## 4       HC__28        PAIRED              PCR         OTHER
    ## 5       HC__27        PAIRED              PCR         OTHER
    ## 6       HC__26        PAIRED              PCR         OTHER
    ##                Organism Platform          ReleaseDate          create_date
    ## 1 human oral metagenome ILLUMINA 2020-10-23T00:00:00Z 2019-09-26T00:28:00Z
    ## 2 human oral metagenome ILLUMINA 2020-10-23T00:00:00Z 2019-09-26T00:28:00Z
    ## 3 human oral metagenome ILLUMINA 2020-10-23T00:00:00Z 2019-09-26T00:28:00Z
    ## 4 human oral metagenome ILLUMINA 2020-10-23T00:00:00Z 2019-09-26T00:28:00Z
    ## 5 human oral metagenome ILLUMINA 2020-10-23T00:00:00Z 2019-09-26T00:28:00Z
    ## 6 human oral metagenome ILLUMINA 2020-10-23T00:00:00Z 2019-09-26T00:28:00Z
    ##   version Sample.Name SRA.Study
    ## 1       1      HC__31 SRP223308
    ## 2       1      HC__30 SRP223308
    ## 3       1      HC__29 SRP223308
    ## 4       1      HC__28 SRP223308
    ## 5       1      HC__27 SRP223308
    ## 6       1      HC__26 SRP223308

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
  BioProject = initial_meta$BioProject,
  Run = initial_meta$Run,
  sample_name = initial_meta$Library.Name,
  Age = NA,
  BMI = NA,
  Sex=NA,
  Region = "V3-V4",
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  Site =initial_meta$isolation_source,
  public_time = initial_meta$ReleaseDate,
  Status="Healthy Control",
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$PCR_primers <- "341F-805R"
metadata$title <- "Elevated Dietary Carbohydrate and Glycemic Intake Associate with an Altered Oral Microbial Ecosystem in Two Large U.S. Cohorts"
#metadata<-subset(uniform_meta,Site!="feces")
metadata$Disease_study <- "YES"
# Assign disease status based on sample_name patterns
#extract subject id by remove prefix of sample_name
#metadata$sample_name <- gsub("[FSP].*", "", metadata$sample_name)
#get subjects id with Saliva available
#keep_subjects<-subset(metadata,Site %in% c("Saliva"))
#remove subjects without Saliva samples
#metadata<-subset(metadata,sample_name %in% keep_subjects$sample_name)
#get subjects sampling counts
#sample_counts <- table(metadata$sample_name)
# Create matched_samples column, any subjects have >=2 samples will be label "YES", otherwise "NO"
#metadata$matched_samples <- ifelse(metadata$sample_name %in% names(sample_counts[sample_counts >= 2]), "YES", "NO")
#metadata$BMI<-gsub("not provided",NA,metadata$BMI)
#metadata$BMI<-as.numeric(metadata$BMI)
metadata$matched_samples<-"NO"
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 384"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 384"

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
  "PCR_primers", "title", "Disease_study", "Status","matched_samples"
)]

# Export the uniform metadata to a CSV file
write.csv(metadata, file = "metadata.csv", quote = FALSE)
head(metadata)
```

    ##    BioProject         Run sample_name Age BMI Sex Region       Platform Country
    ## 1 PRJNA557511 SRR10188199      HC__31  NA  NA  NA  V3-V4 Illumina MiSeq   China
    ## 2 PRJNA557511 SRR10188200      HC__30  NA  NA  NA  V3-V4 Illumina MiSeq   China
    ## 3 PRJNA557511 SRR10188201      HC__29  NA  NA  NA  V3-V4 Illumina MiSeq   China
    ## 4 PRJNA557511 SRR10188202      HC__28  NA  NA  NA  V3-V4 Illumina MiSeq   China
    ## 5 PRJNA557511 SRR10188203      HC__27  NA  NA  NA  V3-V4 Illumina MiSeq   China
    ## 6 PRJNA557511 SRR10188204      HC__26  NA  NA  NA  V3-V4 Illumina MiSeq   China
    ##     Site          public_time PCR_primers
    ## 1 saliva 2020-10-23T00:00:00Z   341F-805R
    ## 2 saliva 2020-10-23T00:00:00Z   341F-805R
    ## 3 saliva 2020-10-23T00:00:00Z   341F-805R
    ## 4 saliva 2020-10-23T00:00:00Z   341F-805R
    ## 5 saliva 2020-10-23T00:00:00Z   341F-805R
    ## 6 saliva 2020-10-23T00:00:00Z   341F-805R
    ##                                                                                                                            title
    ## 1 Elevated Dietary Carbohydrate and Glycemic Intake Associate with an Altered Oral Microbial Ecosystem in Two Large U.S. Cohorts
    ## 2 Elevated Dietary Carbohydrate and Glycemic Intake Associate with an Altered Oral Microbial Ecosystem in Two Large U.S. Cohorts
    ## 3 Elevated Dietary Carbohydrate and Glycemic Intake Associate with an Altered Oral Microbial Ecosystem in Two Large U.S. Cohorts
    ## 4 Elevated Dietary Carbohydrate and Glycemic Intake Associate with an Altered Oral Microbial Ecosystem in Two Large U.S. Cohorts
    ## 5 Elevated Dietary Carbohydrate and Glycemic Intake Associate with an Altered Oral Microbial Ecosystem in Two Large U.S. Cohorts
    ## 6 Elevated Dietary Carbohydrate and Glycemic Intake Associate with an Altered Oral Microbial Ecosystem in Two Large U.S. Cohorts
    ##   Disease_study          Status matched_samples
    ## 1           YES Healthy Control              NO
    ## 2           YES Healthy Control              NO
    ## 3           YES Healthy Control              NO
    ## 4           YES Healthy Control              NO
    ## 5           YES Healthy Control              NO
    ## 6           YES Healthy Control              NO

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
