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
setwd("D:/oral_project/PRJNA545410")

# Read the initial metadata files
initial_meta<-read.csv("SraRunTable.csv",header = T)

# Combine the two metadata files
#meta <- rbind(initial_meta, initial_meta1)

#Remove entries non-oral samples
#initial_meta <- subset(initial_meta,isolate %in% c("saliva","tongue coat"))

# Check the number of samples containing "HT" in Sample.Name (e.g., healthy controls)
#nrow(subset(meta, grepl("HT", Sample.Name)))

# Display the first few rows of the metadata
head(initial_meta)
```

    ##          Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample
    ## 1 SRR9203169   AMPLICON        598 72829384 PRJNA545410 SAMN11891204
    ## 2 SRR9203170   AMPLICON        598 56954680 PRJNA545410 SAMN11891202
    ## 3 SRR9203171   AMPLICON        598 61831219 PRJNA545410 SAMN11891202
    ## 4 SRR9203172   AMPLICON        598 64675222 PRJNA545410 SAMN11891205
    ## 5 SRR9203173   AMPLICON        598 60968801 PRJNA545410 SAMN11891204
    ## 6 SRR9203174   AMPLICON        598 67037067 PRJNA545410 SAMN11891202
    ##                BioSampleModel    Bytes Center.Name Collection_Date Consent
    ## 1 Metagenome or environmental 43298990     FISABIO      2018-01-01  public
    ## 2 Metagenome or environmental 33568295     FISABIO      2018-01-01  public
    ## 3 Metagenome or environmental 37368563     FISABIO      2018-01-01  public
    ## 4 Metagenome or environmental 39054846     FISABIO      2018-01-01  public
    ## 5 Metagenome or environmental 36155173     FISABIO      2018-01-01  public
    ## 6 Metagenome or environmental 40955645     FISABIO      2018-01-01  public
    ##   DATASTORE.filetype DATASTORE.provider                     DATASTORE.region
    ## 1   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 2   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 3   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 4   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 5   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 6   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ##   Experiment geo_loc_name_country geo_loc_name_country_continent geo_loc_name
    ## 1 SRX5974447                Spain                         Europe        Spain
    ## 2 SRX5974446                Spain                         Europe        Spain
    ## 3 SRX5974445                Spain                         Europe        Spain
    ## 4 SRX5974444                Spain                         Europe        Spain
    ## 5 SRX5974443                Spain                         Europe        Spain
    ## 6 SRX5974442                Spain                         Europe        Spain
    ##           HOST     Instrument             lat_lon    Library.Name LibraryLayout
    ## 1 Homo sapiens Illumina MiSeq 40.4165 N 3.70256 W 115_R1.fastq.gz        PAIRED
    ## 2 Homo sapiens Illumina MiSeq 40.4165 N 3.70256 W 116_R1.fastq.gz        PAIRED
    ## 3 Homo sapiens Illumina MiSeq 40.4165 N 3.70256 W 117_R1.fastq.gz        PAIRED
    ## 4 Homo sapiens Illumina MiSeq 40.4165 N 3.70256 W 118_R1.fastq.gz        PAIRED
    ## 5 Homo sapiens Illumina MiSeq 40.4165 N 3.70256 W 119_R1.fastq.gz        PAIRED
    ## 6 Homo sapiens Illumina MiSeq 40.4165 N 3.70256 W 120_R1.fastq.gz        PAIRED
    ##   LibrarySelection LibrarySource              Organism Platform
    ## 1              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 2              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 3              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 4              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 5              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 6              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ##            ReleaseDate          create_date version
    ## 1 2019-06-05T00:00:00Z 2019-06-05T10:06:00Z       1
    ## 2 2019-06-05T00:00:00Z 2019-06-05T10:06:00Z       1
    ## 3 2019-06-05T00:00:00Z 2019-06-05T10:06:00Z       1
    ## 4 2019-06-05T00:00:00Z 2019-06-05T10:07:00Z       1
    ## 5 2019-06-05T00:00:00Z 2019-06-05T10:07:00Z       1
    ## 6 2019-06-05T00:00:00Z 2019-06-05T10:06:00Z       1
    ##                              Sample.Name SRA.Study isolation_source
    ## 1    Dental plaque vestibular microbiota SRP200474                 
    ## 2 Dental plaque interproximal microbiota SRP200474    dental plaque
    ## 3 Dental plaque interproximal microbiota SRP200474    dental plaque
    ## 4       Dental plaque lingual microbiota SRP200474                 
    ## 5    Dental plaque vestibular microbiota SRP200474                 
    ## 6 Dental plaque interproximal microbiota SRP200474    dental plaque

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
  sample_name = initial_meta$Library.Name,
  Age = NA,
  BMI = NA,
  Sex = NA,
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
metadata$PCR_primers <- NA
metadata$title <- "High-throughput DNA sequencing of microbiota at interproximal sites"
metadata$Disease_study <- "NO"
# Assign disease status based on sample_name patterns
metadata$Status<-"YES"

#get subjects id according to sample id
metadata$sample_name<-gsub("_R1.fastq.gz","",metadata$sample_name)
#get subjects id with Saliva available
keep_subjects<-subset(metadata,Site %in% "dental plaque")
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

    ## [1] "Number of unique runs: 20"

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
  "PCR_primers", "title", "Disease_study", "Status","matched_samples"
)]

# Export the uniform metadata to a CSV file
write.csv(metadata, file = "metadata.csv", quote = TRUE)
head(metadata)
```

    ##     BioProject        Run sample_name Age BMI Sex Region       Platform Country
    ## 2  PRJNA545410 SRR9203170         116  NA  NA  NA  V3-V4 Illumina MiSeq   Spain
    ## 3  PRJNA545410 SRR9203171         117  NA  NA  NA  V3-V4 Illumina MiSeq   Spain
    ## 6  PRJNA545410 SRR9203174         120  NA  NA  NA  V3-V4 Illumina MiSeq   Spain
    ## 7  PRJNA545410 SRR9203175         121  NA  NA  NA  V3-V4 Illumina MiSeq   Spain
    ## 10 PRJNA545410 SRR9203178         124  NA  NA  NA  V3-V4 Illumina MiSeq   Spain
    ## 24 PRJNA545410 SRR9203192         112  NA  NA  NA  V3-V4 Illumina MiSeq   Spain
    ##             Site          public_time PCR_primers
    ## 2  dental plaque 2019-06-05T00:00:00Z          NA
    ## 3  dental plaque 2019-06-05T00:00:00Z          NA
    ## 6  dental plaque 2019-06-05T00:00:00Z          NA
    ## 7  dental plaque 2019-06-05T00:00:00Z          NA
    ## 10 dental plaque 2019-06-05T00:00:00Z          NA
    ## 24 dental plaque 2019-06-05T00:00:00Z          NA
    ##                                                                  title
    ## 2  High-throughput DNA sequencing of microbiota at interproximal sites
    ## 3  High-throughput DNA sequencing of microbiota at interproximal sites
    ## 6  High-throughput DNA sequencing of microbiota at interproximal sites
    ## 7  High-throughput DNA sequencing of microbiota at interproximal sites
    ## 10 High-throughput DNA sequencing of microbiota at interproximal sites
    ## 24 High-throughput DNA sequencing of microbiota at interproximal sites
    ##    Disease_study Status matched_samples
    ## 2             NO    YES              NO
    ## 3             NO    YES              NO
    ## 6             NO    YES              NO
    ## 7             NO    YES              NO
    ## 10            NO    YES              NO
    ## 24            NO    YES              NO

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
