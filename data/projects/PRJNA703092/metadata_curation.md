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
setwd("D:/oral_project/PRJNA703092")

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
    ## 1 SRR13754630   AMPLICON        394 26493889 PRJNA703092 SAMN18011321
    ## 2 SRR13754636   AMPLICON        394 14267763 PRJNA703092 SAMN18011345
    ## 3 SRR13754644   AMPLICON        395 29226730 PRJNA703092 SAMN18011343
    ## 4 SRR13754650   AMPLICON        394 26371121 PRJNA703092 SAMN18011313
    ## 5 SRR13754656   AMPLICON        392 17872345 PRJNA703092 SAMN18011366
    ## 6 SRR13754657   AMPLICON        393 19176083 PRJNA703092 SAMN18011332
    ##                BioSampleModel   Bytes
    ## 1 Metagenome or environmental 8110688
    ## 2 Metagenome or environmental 3778578
    ## 3 Metagenome or environmental 9075324
    ## 4 Metagenome or environmental 8119341
    ## 5 Metagenome or environmental 5353540
    ## 6 Metagenome or environmental 5928996
    ##                                                                               Center.Name
    ## 1 SHANGHAI NINTH PEOPLE'S HOSPITAL, COLLEGE OF STOMATOLOGY, SHANGHAI JIAO TONG UNIVERSITY
    ## 2 SHANGHAI NINTH PEOPLE'S HOSPITAL, COLLEGE OF STOMATOLOGY, SHANGHAI JIAO TONG UNIVERSITY
    ## 3 SHANGHAI NINTH PEOPLE'S HOSPITAL, COLLEGE OF STOMATOLOGY, SHANGHAI JIAO TONG UNIVERSITY
    ## 4 SHANGHAI NINTH PEOPLE'S HOSPITAL, COLLEGE OF STOMATOLOGY, SHANGHAI JIAO TONG UNIVERSITY
    ## 5 SHANGHAI NINTH PEOPLE'S HOSPITAL, COLLEGE OF STOMATOLOGY, SHANGHAI JIAO TONG UNIVERSITY
    ## 6 SHANGHAI NINTH PEOPLE'S HOSPITAL, COLLEGE OF STOMATOLOGY, SHANGHAI JIAO TONG UNIVERSITY
    ##   Collection_Date Consent DATASTORE.filetype DATASTORE.provider
    ## 1      2013-04-12  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 2      2013-04-14  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 3      2013-04-14  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 4      2013-04-12  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5      2013-10-17  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 6      2013-04-15  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                       DATASTORE.region  Experiment geo_loc_name_country
    ## 1 gs.us-east1,ncbi.public,s3.us-east-1 SRX10141464                China
    ## 2 gs.us-east1,ncbi.public,s3.us-east-1 SRX10141458                China
    ## 3 gs.us-east1,ncbi.public,s3.us-east-1 SRX10141450                China
    ## 4 gs.us-east1,ncbi.public,s3.us-east-1 SRX10141444                China
    ## 5 gs.us-east1,ncbi.public,s3.us-east-1 SRX10141438                China
    ## 6 gs.us-east1,ncbi.public,s3.us-east-1 SRX10141437                China
    ##   geo_loc_name_country_continent    geo_loc_name         HOST     Instrument
    ## 1                           Asia China: Shanghai Homo sapiens Illumina MiSeq
    ## 2                           Asia China: Shanghai Homo sapiens Illumina MiSeq
    ## 3                           Asia China: Shanghai Homo sapiens Illumina MiSeq
    ## 4                           Asia China: Shanghai Homo sapiens Illumina MiSeq
    ## 5                           Asia China: Shanghai Homo sapiens Illumina MiSeq
    ## 6                           Asia China: Shanghai Homo sapiens Illumina MiSeq
    ##   isolation_source        lat_lon Library.Name LibraryLayout LibrarySelection
    ## 1    dental plaque not applicable         1052        SINGLE              PCR
    ## 2    dental plaque not applicable         2172        SINGLE              PCR
    ## 3    dental plaque not applicable         2132        SINGLE              PCR
    ## 4    dental plaque not applicable         1042        SINGLE              PCR
    ## 5    dental plaque not applicable         2063        SINGLE              PCR
    ## 6    dental plaque not applicable         2062        SINGLE              PCR
    ##   LibrarySource              Organism Platform          ReleaseDate
    ## 1   METAGENOMIC human oral metagenome ILLUMINA 2021-02-24T00:00:00Z
    ## 2   METAGENOMIC human oral metagenome ILLUMINA 2021-02-24T00:00:00Z
    ## 3   METAGENOMIC human oral metagenome ILLUMINA 2021-02-24T00:00:00Z
    ## 4   METAGENOMIC human oral metagenome ILLUMINA 2021-02-24T00:00:00Z
    ## 5   METAGENOMIC human oral metagenome ILLUMINA 2021-02-24T00:00:00Z
    ## 6   METAGENOMIC human oral metagenome ILLUMINA 2021-02-24T00:00:00Z
    ##            create_date version Sample.Name SRA.Study
    ## 1 2021-02-24T02:10:00Z       1        1052 SRP307227
    ## 2 2021-02-24T02:12:00Z       1        2172 SRP307227
    ## 3 2021-02-24T02:13:00Z       1        2132 SRP307227
    ## 4 2021-02-24T02:13:00Z       1        1042 SRP307227
    ## 5 2021-02-24T02:14:00Z       1        2063 SRP307227
    ## 6 2021-02-24T02:14:00Z       1        2062 SRP307227

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
  Region = "V4-V5",
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  Site =initial_meta$isolation_source,
  public_time = initial_meta$ReleaseDate,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$PCR_primers <- "515F-907R"
metadata$title <- "Dynamic Observation of the Effect of Maternal Caries on the Oral Microbiota of Infants Aged 12–24 Months"
#metadata<-subset(uniform_meta,Site!="feces")
metadata$Disease_study <- "YES"
# Assign disease status based on sample_name patterns
metadata$Status <- NA
#extract subject id by remove prefix of sample_name
#metadata$sample_name <- gsub("[FSP].*", "", metadata$sample_name)
#get subjects id with Saliva available
#keep_subjects<-subset(metadata,Site %in% c("oral swab"))
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

    ## [1] "Number of unique runs: 108"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 108"

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
    ## 1 PRJNA703092 SRR13754630        1052  NA  NA  NA  V4-V5 Illumina MiSeq   China
    ## 2 PRJNA703092 SRR13754636        2172  NA  NA  NA  V4-V5 Illumina MiSeq   China
    ## 3 PRJNA703092 SRR13754644        2132  NA  NA  NA  V4-V5 Illumina MiSeq   China
    ## 4 PRJNA703092 SRR13754650        1042  NA  NA  NA  V4-V5 Illumina MiSeq   China
    ## 5 PRJNA703092 SRR13754656        2063  NA  NA  NA  V4-V5 Illumina MiSeq   China
    ## 6 PRJNA703092 SRR13754657        2062  NA  NA  NA  V4-V5 Illumina MiSeq   China
    ##            Site          public_time PCR_primers
    ## 1 dental plaque 2021-02-24T00:00:00Z   515F-907R
    ## 2 dental plaque 2021-02-24T00:00:00Z   515F-907R
    ## 3 dental plaque 2021-02-24T00:00:00Z   515F-907R
    ## 4 dental plaque 2021-02-24T00:00:00Z   515F-907R
    ## 5 dental plaque 2021-02-24T00:00:00Z   515F-907R
    ## 6 dental plaque 2021-02-24T00:00:00Z   515F-907R
    ##                                                                                                      title
    ## 1 Dynamic Observation of the Effect of Maternal Caries on the Oral Microbiota of Infants Aged 12–24 Months
    ## 2 Dynamic Observation of the Effect of Maternal Caries on the Oral Microbiota of Infants Aged 12–24 Months
    ## 3 Dynamic Observation of the Effect of Maternal Caries on the Oral Microbiota of Infants Aged 12–24 Months
    ## 4 Dynamic Observation of the Effect of Maternal Caries on the Oral Microbiota of Infants Aged 12–24 Months
    ## 5 Dynamic Observation of the Effect of Maternal Caries on the Oral Microbiota of Infants Aged 12–24 Months
    ## 6 Dynamic Observation of the Effect of Maternal Caries on the Oral Microbiota of Infants Aged 12–24 Months
    ##   Disease_study Status matched_samples
    ## 1           YES     NA              NO
    ## 2           YES     NA              NO
    ## 3           YES     NA              NO
    ## 4           YES     NA              NO
    ## 5           YES     NA              NO
    ## 6           YES     NA              NO

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
