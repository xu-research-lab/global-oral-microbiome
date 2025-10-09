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
setwd("D:/oral_project/PRJNA893253")

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
    ## 1 SRR22052817   AMPLICON        500 40816500 PRJNA893253 SAMN31416460
    ## 2 SRR22052818   AMPLICON        500 40078000 PRJNA893253 SAMN31416459
    ## 3 SRR22052819   AMPLICON        500 40733500 PRJNA893253 SAMN31416458
    ## 4 SRR22052820   AMPLICON        500 40407000 PRJNA893253 SAMN31416457
    ## 5 SRR22052821   AMPLICON        500 40434500 PRJNA893253 SAMN31416456
    ## 6 SRR22052825   AMPLICON        500 40860500 PRJNA893253 SAMN31416467
    ##                BioSampleModel    Bytes                 Center.Name
    ## 1 Metagenome or environmental 14781482 THE SECOND XIANGYA HOSPITAL
    ## 2 Metagenome or environmental 14021600 THE SECOND XIANGYA HOSPITAL
    ## 3 Metagenome or environmental 14905985 THE SECOND XIANGYA HOSPITAL
    ## 4 Metagenome or environmental 14025573 THE SECOND XIANGYA HOSPITAL
    ## 5 Metagenome or environmental 13406938 THE SECOND XIANGYA HOSPITAL
    ## 6 Metagenome or environmental 14577835 THE SECOND XIANGYA HOSPITAL
    ##   Collection_Date Consent DATASTORE.filetype DATASTORE.provider
    ## 1      2021-03-01  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 2      2021-03-01  public   sra,fastq,run.zq         gs,ncbi,s3
    ## 3      2021-03-01  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 4      2021-03-01  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5      2021-03-01  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 6      2021-03-01  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                       DATASTORE.region  Experiment geo_loc_name_country
    ## 1 gs.us-east1,ncbi.public,s3.us-east-1 SRX18033894                China
    ## 2 gs.us-east1,ncbi.public,s3.us-east-1 SRX18033893                China
    ## 3 gs.us-east1,ncbi.public,s3.us-east-1 SRX18033892                China
    ## 4 gs.us-east1,ncbi.public,s3.us-east-1 SRX18033891                China
    ## 5 gs.us-east1,ncbi.public,s3.us-east-1 SRX18033890                China
    ## 6 gs.us-east1,ncbi.public,s3.us-east-1 SRX18033886                China
    ##   geo_loc_name_country_continent   geo_loc_name          HOST
    ## 1                           Asia China:Changsha not collected
    ## 2                           Asia China:Changsha not collected
    ## 3                           Asia China:Changsha not collected
    ## 4                           Asia China:Changsha not collected
    ## 5                           Asia China:Changsha not collected
    ## 6                           Asia China:Changsha not collected
    ##              Instrument   isolation_source           lat_lon Library.Name
    ## 1 Illumina NovaSeq 6000 subgingival plaque 28.19 N 112.104 E            7
    ## 2 Illumina NovaSeq 6000 subgingival plaque 28.19 N 112.103 E            6
    ## 3 Illumina NovaSeq 6000 subgingival plaque 28.19 N 112.102 E            5
    ## 4 Illumina NovaSeq 6000 subgingival plaque 28.19 N 112.101 E            4
    ## 5 Illumina NovaSeq 6000 subgingival plaque 28.19 N 112.100 E            3
    ## 6 Illumina NovaSeq 6000 subgingival plaque 28.19 N 112.111 E           14
    ##   LibraryLayout LibrarySelection LibrarySource        Organism Platform
    ## 1        PAIRED              PCR   METAGENOMIC oral metagenome ILLUMINA
    ## 2        PAIRED              PCR   METAGENOMIC oral metagenome ILLUMINA
    ## 3        PAIRED              PCR   METAGENOMIC oral metagenome ILLUMINA
    ## 4        PAIRED              PCR   METAGENOMIC oral metagenome ILLUMINA
    ## 5        PAIRED              PCR   METAGENOMIC oral metagenome ILLUMINA
    ## 6        PAIRED              PCR   METAGENOMIC oral metagenome ILLUMINA
    ##            ReleaseDate          create_date version Sample.Name SRA.Study
    ## 1 2022-10-26T00:00:00Z 2022-10-26T15:52:00Z       1       G1.09 SRP404688
    ## 2 2022-10-26T00:00:00Z 2022-10-26T15:52:00Z       1       G1.08 SRP404688
    ## 3 2022-10-26T00:00:00Z 2022-10-26T15:52:00Z       1       G1.07 SRP404688
    ## 4 2022-10-26T00:00:00Z 2022-10-26T15:52:00Z       1       G1.05 SRP404688
    ## 5 2022-10-26T00:00:00Z 2022-10-26T15:53:00Z       1       G1.04 SRP404688
    ## 6 2022-10-26T00:00:00Z 2022-10-26T15:53:00Z       1       G1.17 SRP404688

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
  sample_name = initial_meta$Sample.Name,
  Age = NA,
  BMI = NA,
  Sex=NA,
  Region = "V3-V4",
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  Site = initial_meta$isolation_source,
  public_time = initial_meta$ReleaseDate,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$PCR_primers <- "343F-798R"
metadata$title <- "Prevotella genus and its related NOD-like receptor signaling pathway in young males with stage III periodontitis"
#metadata<-subset(uniform_meta,Site!="feces")
metadata$Disease_study <- "YES"
# Assign disease status based on sample_name patterns
metadata$Status <- "periodontitis"
#extract subject id by remove prefix of sample_name
#metadata$sample_name <- gsub("S-Saliva-.*|F", "", metadata$sample_name)
#get subjects id with Saliva available
#keep_subjects<-subset(metadata,Site %in% c("Saliva"))
#remove subjects without Saliva samples
#metadata<-subset(metadata,sample_name %in% keep_subjects$sample_name)
#get subjects sampling counts
#sample_counts <- table(metadata$sample_name)
# Create matched_samples column, any subjects have >=2 samples will be label "YES", otherwise "NO"
#metadata$matched_samples <- ifelse(metadata$sample_name %in% names(sample_counts[sample_counts >= 2]), "YES", "NO")
metadata$matched_samples<-"NO"
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 17"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 17"

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

    ##    BioProject         Run sample_name Age BMI Sex Region              Platform
    ## 1 PRJNA893253 SRR22052817       G1.09  NA  NA  NA  V3-V4 Illumina NovaSeq 6000
    ## 2 PRJNA893253 SRR22052818       G1.08  NA  NA  NA  V3-V4 Illumina NovaSeq 6000
    ## 3 PRJNA893253 SRR22052819       G1.07  NA  NA  NA  V3-V4 Illumina NovaSeq 6000
    ## 4 PRJNA893253 SRR22052820       G1.05  NA  NA  NA  V3-V4 Illumina NovaSeq 6000
    ## 5 PRJNA893253 SRR22052821       G1.04  NA  NA  NA  V3-V4 Illumina NovaSeq 6000
    ## 6 PRJNA893253 SRR22052825       G1.17  NA  NA  NA  V3-V4 Illumina NovaSeq 6000
    ##   Country               Site          public_time PCR_primers
    ## 1   China subgingival plaque 2022-10-26T00:00:00Z   343F-798R
    ## 2   China subgingival plaque 2022-10-26T00:00:00Z   343F-798R
    ## 3   China subgingival plaque 2022-10-26T00:00:00Z   343F-798R
    ## 4   China subgingival plaque 2022-10-26T00:00:00Z   343F-798R
    ## 5   China subgingival plaque 2022-10-26T00:00:00Z   343F-798R
    ## 6   China subgingival plaque 2022-10-26T00:00:00Z   343F-798R
    ##                                                                                                              title
    ## 1 Prevotella genus and its related NOD-like receptor signaling pathway in young males with stage III periodontitis
    ## 2 Prevotella genus and its related NOD-like receptor signaling pathway in young males with stage III periodontitis
    ## 3 Prevotella genus and its related NOD-like receptor signaling pathway in young males with stage III periodontitis
    ## 4 Prevotella genus and its related NOD-like receptor signaling pathway in young males with stage III periodontitis
    ## 5 Prevotella genus and its related NOD-like receptor signaling pathway in young males with stage III periodontitis
    ## 6 Prevotella genus and its related NOD-like receptor signaling pathway in young males with stage III periodontitis
    ##   Disease_study        Status matched_samples
    ## 1           YES periodontitis              NO
    ## 2           YES periodontitis              NO
    ## 3           YES periodontitis              NO
    ## 4           YES periodontitis              NO
    ## 5           YES periodontitis              NO
    ## 6           YES periodontitis              NO

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
