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
library(readxl)

# Set working directory (update this path to your local directory)
setwd("/home/dongbiao/word_embedding_microbiome/all_data/oral/global-oral-microbiome/data/projects/PRJNA684772")

# Read the initial metadata files
initial_meta<-read.csv("SraRunTable.csv") 

head(initial_meta)
```

    ##           Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample
    ## 1 SRR13261442   AMPLICON        602  9216620 PRJNA684772 SAMN17069167
    ## 2 SRR13261443   AMPLICON        602 11063556 PRJNA684772 SAMN17069166
    ## 3 SRR13261444   AMPLICON        602 18284546 PRJNA684772 SAMN17069165
    ## 4 SRR13261445   AMPLICON        602 13328882 PRJNA684772 SAMN17069164
    ## 5 SRR13261446   AMPLICON        602  9217824 PRJNA684772 SAMN17069163
    ## 6 SRR13261447   AMPLICON        602 16902354 PRJNA684772 SAMN17069162
    ##                                BioSampleModel    Bytes
    ## 1 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral  5661958
    ## 2 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral  6835002
    ## 3 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 11229051
    ## 4 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral  8175732
    ## 5 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral  5677628
    ## 6 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 10339748
    ##                                Center.Name Collection_Date Consent
    ## 1 SCHOOL OF FOOD SCIENCE AND BIOTECHNOLOGY      2018-05-25  public
    ## 2 SCHOOL OF FOOD SCIENCE AND BIOTECHNOLOGY      2018-05-25  public
    ## 3 SCHOOL OF FOOD SCIENCE AND BIOTECHNOLOGY      2018-05-25  public
    ## 4 SCHOOL OF FOOD SCIENCE AND BIOTECHNOLOGY      2018-05-25  public
    ## 5 SCHOOL OF FOOD SCIENCE AND BIOTECHNOLOGY      2018-05-25  public
    ## 6 SCHOOL OF FOOD SCIENCE AND BIOTECHNOLOGY      2018-05-25  public
    ##   DATASTORE.filetype DATASTORE.provider                     DATASTORE.region
    ## 1   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 2   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 3   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 4   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 5   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 6   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ##   env_broad_scale env_local_scale env_medium Experiment geo_loc_name_country
    ## 1     oral cavity     oral helath     saliva SRX9692032                China
    ## 2     oral cavity     oral helath     saliva SRX9692031                China
    ## 3     oral cavity     oral helath     saliva SRX9692030                China
    ## 4     oral cavity     oral helath     saliva SRX9692029                China
    ## 5     oral cavity     oral helath     saliva SRX9692028                China
    ## 6     oral cavity     oral helath     saliva SRX9692027                China
    ##   geo_loc_name_country_continent   geo_loc_name         HOST     Instrument
    ## 1                           Asia China:Hangzhou Homo sapiens Illumina MiSeq
    ## 2                           Asia China:Hangzhou Homo sapiens Illumina MiSeq
    ## 3                           Asia China:Hangzhou Homo sapiens Illumina MiSeq
    ## 4                           Asia China:Hangzhou Homo sapiens Illumina MiSeq
    ## 5                           Asia China:Hangzhou Homo sapiens Illumina MiSeq
    ## 6                           Asia China:Hangzhou Homo sapiens Illumina MiSeq
    ##      lat_lon Library.Name LibraryLayout LibrarySelection LibrarySource
    ## 1 30 N 120 E         SC06        PAIRED              PCR   METAGENOMIC
    ## 2 30 N 120 E         SC05        PAIRED              PCR   METAGENOMIC
    ## 3 30 N 120 E         SC04        PAIRED              PCR   METAGENOMIC
    ## 4 30 N 120 E         SC03        PAIRED              PCR   METAGENOMIC
    ## 5 30 N 120 E         SC02        PAIRED              PCR   METAGENOMIC
    ## 6 30 N 120 E         SC01        PAIRED              PCR   METAGENOMIC
    ##                Organism Platform          ReleaseDate          create_date
    ## 1 human oral metagenome ILLUMINA 2020-12-19T00:00:00Z 2020-12-19T20:37:00Z
    ## 2 human oral metagenome ILLUMINA 2020-12-19T00:00:00Z 2020-12-19T20:37:00Z
    ## 3 human oral metagenome ILLUMINA 2020-12-19T00:00:00Z 2020-12-19T20:37:00Z
    ## 4 human oral metagenome ILLUMINA 2020-12-19T00:00:00Z 2020-12-19T20:37:00Z
    ## 5 human oral metagenome ILLUMINA 2020-12-19T00:00:00Z 2020-12-19T20:37:00Z
    ## 6 human oral metagenome ILLUMINA 2020-12-19T00:00:00Z 2020-12-19T20:37:00Z
    ##   version Sample.Name SRA.Study
    ## 1       1        SC06 SRP298210
    ## 2       1        SC05 SRP298210
    ## 3       1        SC04 SRP298210
    ## 4       1        SC03 SRP298210
    ## 5       1        SC02 SRP298210
    ## 6       1        SC01 SRP298210

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
  Country = initial_meta$geo_loc_name,
  Site = initial_meta$env_medium,
  public_time = initial_meta$ReleaseDate,
  stringsAsFactors = FALSE
)
metadata<-unique(metadata)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$PCR_primers <- "338F-806R"
metadata$title <- "Tooth brushing using toothpaste containing theaflavins reduces the oral pathogenic bacteria in healthy adults"
metadata$Disease_study <- "NO"

# Assign disease status based on sample_name patterns
metadata$Status <- "HC"
metadata$matched_samples <- "NO"
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 80"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 80"

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
metadata<-unique(metadata)
write.csv(metadata, file = "metadata.csv", quote = FALSE)
head(metadata)
```

    ##    BioProject         Run sample_name Age BMI Sex Region       Platform
    ## 1 PRJNA684772 SRR13261442        SC06  NA  NA  NA  V3-V4 Illumina MiSeq
    ## 2 PRJNA684772 SRR13261443        SC05  NA  NA  NA  V3-V4 Illumina MiSeq
    ## 3 PRJNA684772 SRR13261444        SC04  NA  NA  NA  V3-V4 Illumina MiSeq
    ## 4 PRJNA684772 SRR13261445        SC03  NA  NA  NA  V3-V4 Illumina MiSeq
    ## 5 PRJNA684772 SRR13261446        SC02  NA  NA  NA  V3-V4 Illumina MiSeq
    ## 6 PRJNA684772 SRR13261447        SC01  NA  NA  NA  V3-V4 Illumina MiSeq
    ##          Country   Site          public_time PCR_primers
    ## 1 China:Hangzhou saliva 2020-12-19T00:00:00Z   338F-806R
    ## 2 China:Hangzhou saliva 2020-12-19T00:00:00Z   338F-806R
    ## 3 China:Hangzhou saliva 2020-12-19T00:00:00Z   338F-806R
    ## 4 China:Hangzhou saliva 2020-12-19T00:00:00Z   338F-806R
    ## 5 China:Hangzhou saliva 2020-12-19T00:00:00Z   338F-806R
    ## 6 China:Hangzhou saliva 2020-12-19T00:00:00Z   338F-806R
    ##                                                                                                           title
    ## 1 Tooth brushing using toothpaste containing theaflavins reduces the oral pathogenic bacteria in healthy adults
    ## 2 Tooth brushing using toothpaste containing theaflavins reduces the oral pathogenic bacteria in healthy adults
    ## 3 Tooth brushing using toothpaste containing theaflavins reduces the oral pathogenic bacteria in healthy adults
    ## 4 Tooth brushing using toothpaste containing theaflavins reduces the oral pathogenic bacteria in healthy adults
    ## 5 Tooth brushing using toothpaste containing theaflavins reduces the oral pathogenic bacteria in healthy adults
    ## 6 Tooth brushing using toothpaste containing theaflavins reduces the oral pathogenic bacteria in healthy adults
    ##   Disease_study Status matched_samples
    ## 1            NO     HC              NO
    ## 2            NO     HC              NO
    ## 3            NO     HC              NO
    ## 4            NO     HC              NO
    ## 5            NO     HC              NO
    ## 6            NO     HC              NO

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
