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
setwd("/home/dongbiao/word_embedding_microbiome/all_data/oral/global-oral-microbiome/data/projects/PRJDB8324")

# Read the initial metadata files
initial_meta<-read.csv("SraRunTable.csv")

head(initial_meta)
```

    ##         Run Assay.Type AvgSpotLen   Bases BioProject    BioSample
    ## 1 DRR180822   AMPLICON        355 5757076  PRJDB8324 SAMD00173638
    ## 2 DRR180823   AMPLICON        356 7680406  PRJDB8324 SAMD00173639
    ## 3 DRR180824   AMPLICON        354 6356123  PRJDB8324 SAMD00173640
    ## 4 DRR180825   AMPLICON        354 7286306  PRJDB8324 SAMD00173641
    ## 5 DRR180826   AMPLICON        354 9636865  PRJDB8324 SAMD00173642
    ## 6 DRR180827   AMPLICON        352 6925790  PRJDB8324 SAMD00173643
    ##              BioSampleModel   Bytes Center.Name Collection_Date Consent
    ## 1 MIMARKS.survey.human-oral 4325872      KYUSHU            2015  public
    ## 2 MIMARKS.survey.human-oral 5574910      KYUSHU            2015  public
    ## 3 MIMARKS.survey.human-oral 4618314      KYUSHU            2015  public
    ## 4 MIMARKS.survey.human-oral 5200276      KYUSHU            2015  public
    ## 5 MIMARKS.survey.human-oral 6978951      KYUSHU            2015  public
    ## 6 MIMARKS.survey.human-oral 4940888      KYUSHU            2015  public
    ##   DATASTORE.filetype DATASTORE.provider                     DATASTORE.region
    ## 1         run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 2         run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 3         run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 4         run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 5         run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 6         run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ##     env_broad_scale env_local_scale env_medium Experiment geo_loc_name_country
    ## 1 human oral cavity          Saliva     Saliva  DRX171328                Japan
    ## 2 human oral cavity          Saliva     Saliva  DRX171329                Japan
    ## 3 human oral cavity          Saliva     Saliva  DRX171330                Japan
    ## 4 human oral cavity          Saliva     Saliva  DRX171331                Japan
    ## 5 human oral cavity          Saliva     Saliva  DRX171332                Japan
    ## 6 human oral cavity          Saliva     Saliva  DRX171333                Japan
    ##   geo_loc_name_country_continent geo_loc_name         HOST      Instrument
    ## 1                           Asia        Japan Homo sapiens Ion Torrent PGM
    ## 2                           Asia        Japan Homo sapiens Ion Torrent PGM
    ## 3                           Asia        Japan Homo sapiens Ion Torrent PGM
    ## 4                           Asia        Japan Homo sapiens Ion Torrent PGM
    ## 5                           Asia        Japan Homo sapiens Ion Torrent PGM
    ## 6                           Asia        Japan Homo sapiens Ion Torrent PGM
    ##      lat_lon Library.Name LibraryLayout LibrarySelection LibrarySource
    ## 1 33 N 130 E           NA        SINGLE              PCR       GENOMIC
    ## 2 33 N 130 E           NA        SINGLE              PCR       GENOMIC
    ## 3 33 N 130 E           NA        SINGLE              PCR       GENOMIC
    ## 4 33 N 130 E           NA        SINGLE              PCR       GENOMIC
    ## 5 33 N 130 E           NA        SINGLE              PCR       GENOMIC
    ## 6 33 N 130 E           NA        SINGLE              PCR       GENOMIC
    ##                Organism    Platform
    ## 1 human oral metagenome ION_TORRENT
    ## 2 human oral metagenome ION_TORRENT
    ## 3 human oral metagenome ION_TORRENT
    ## 4 human oral metagenome ION_TORRENT
    ## 5 human oral metagenome ION_TORRENT
    ## 6 human oral metagenome ION_TORRENT
    ##                                                project_name
    ## 1 Characteristics of salivary microbiota in cancer patients
    ## 2 Characteristics of salivary microbiota in cancer patients
    ## 3 Characteristics of salivary microbiota in cancer patients
    ## 4 Characteristics of salivary microbiota in cancer patients
    ## 5 Characteristics of salivary microbiota in cancer patients
    ## 6 Characteristics of salivary microbiota in cancer patients
    ##            ReleaseDate          create_date version  Sample.Name
    ## 1 2019-06-12T00:00:00Z 2019-06-14T14:18:00Z       1 SAMD00173638
    ## 2 2019-06-12T00:00:00Z 2019-06-14T14:20:00Z       1 SAMD00173639
    ## 3 2019-06-12T00:00:00Z 2019-06-14T14:19:00Z       1 SAMD00173640
    ## 4 2019-06-12T00:00:00Z 2019-06-14T14:18:00Z       1 SAMD00173641
    ## 5 2019-06-12T00:00:00Z 2019-06-14T14:20:00Z       1 SAMD00173642
    ## 6 2019-06-12T00:00:00Z 2019-06-14T14:19:00Z       1 SAMD00173643
    ##                                        Sample_name SRA.Study
    ## 1 16S amplicon analysis of oral microbiome CAN_001 DRP005151
    ## 2 16S amplicon analysis of oral microbiome CAN_002 DRP005151
    ## 3 16S amplicon analysis of oral microbiome CAN_003 DRP005151
    ## 4 16S amplicon analysis of oral microbiome CAN_004 DRP005151
    ## 5 16S amplicon analysis of oral microbiome CAN_005 DRP005151
    ## 6 16S amplicon analysis of oral microbiome CAN_006 DRP005151

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
  sample_name = initial_meta$Sample_name,
  Age = NA,
  BMI = NA,
  Sex=NA,
  Region = "V1-V2",
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name,
  Site = "Saliva",
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
metadata$PCR_primers <- "8F-338R"
metadata$title <- "Characteristics of the Salivary Microbiota in Patients With Various Digestive Tract Cancers"
metadata$Disease_study <- "YES"

# Assign disease status based on sample_name patterns
metadata$Status <- NA
metadata$matched_samples <- "NO"
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 177"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 177"

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

    ##   BioProject       Run                                      sample_name Age BMI
    ## 1  PRJDB8324 DRR180822 16S amplicon analysis of oral microbiome CAN_001  NA  NA
    ## 2  PRJDB8324 DRR180823 16S amplicon analysis of oral microbiome CAN_002  NA  NA
    ## 3  PRJDB8324 DRR180824 16S amplicon analysis of oral microbiome CAN_003  NA  NA
    ## 4  PRJDB8324 DRR180825 16S amplicon analysis of oral microbiome CAN_004  NA  NA
    ## 5  PRJDB8324 DRR180826 16S amplicon analysis of oral microbiome CAN_005  NA  NA
    ## 6  PRJDB8324 DRR180827 16S amplicon analysis of oral microbiome CAN_006  NA  NA
    ##   Sex Region        Platform Country   Site          public_time PCR_primers
    ## 1  NA  V1-V2 Ion Torrent PGM   Japan Saliva 2019-06-12T00:00:00Z     8F-338R
    ## 2  NA  V1-V2 Ion Torrent PGM   Japan Saliva 2019-06-12T00:00:00Z     8F-338R
    ## 3  NA  V1-V2 Ion Torrent PGM   Japan Saliva 2019-06-12T00:00:00Z     8F-338R
    ## 4  NA  V1-V2 Ion Torrent PGM   Japan Saliva 2019-06-12T00:00:00Z     8F-338R
    ## 5  NA  V1-V2 Ion Torrent PGM   Japan Saliva 2019-06-12T00:00:00Z     8F-338R
    ## 6  NA  V1-V2 Ion Torrent PGM   Japan Saliva 2019-06-12T00:00:00Z     8F-338R
    ##                                                                                         title
    ## 1 Characteristics of the Salivary Microbiota in Patients With Various Digestive Tract Cancers
    ## 2 Characteristics of the Salivary Microbiota in Patients With Various Digestive Tract Cancers
    ## 3 Characteristics of the Salivary Microbiota in Patients With Various Digestive Tract Cancers
    ## 4 Characteristics of the Salivary Microbiota in Patients With Various Digestive Tract Cancers
    ## 5 Characteristics of the Salivary Microbiota in Patients With Various Digestive Tract Cancers
    ## 6 Characteristics of the Salivary Microbiota in Patients With Various Digestive Tract Cancers
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
