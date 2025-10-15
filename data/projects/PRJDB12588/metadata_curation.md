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
setwd("/home/hby/oral_project/PRJDB12588")

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

    ##         Run Assay.Type AvgSpotLen    Bases BioProject    BioSample            BioSampleModel   Bytes Center.Name Collection_Date Consent DATASTORE.filetype DATASTORE.provider
    ## 1 DRR329029   AMPLICON       1471 13143526 PRJDB12588 SAMD00423397 MIMARKS.survey.human-oral 5574232      KYUSHU            2020  public         run.zq,sra         gs,ncbi,s3
    ## 2 DRR329030   AMPLICON       1473 12192077 PRJDB12588 SAMD00423398 MIMARKS.survey.human-oral 5145900      KYUSHU            2020  public         run.zq,sra         gs,ncbi,s3
    ## 3 DRR329031   AMPLICON       1472 13325748 PRJDB12588 SAMD00423399 MIMARKS.survey.human-oral 5643648      KYUSHU            2020  public         run.zq,sra         gs,ncbi,s3
    ## 4 DRR329032   AMPLICON       1475 16292893 PRJDB12588 SAMD00423400 MIMARKS.survey.human-oral 6712729      KYUSHU            2020  public         run.zq,sra         gs,ncbi,s3
    ## 5 DRR329033   AMPLICON       1472  9795597 PRJDB12588 SAMD00423401 MIMARKS.survey.human-oral 4175019      KYUSHU            2020  public         run.zq,sra         gs,ncbi,s3
    ## 6 DRR329034   AMPLICON       1468  6477262 PRJDB12588 SAMD00423402 MIMARKS.survey.human-oral 2328577      KYUSHU            2020  public         run.zq,sra         gs,ncbi,s3
    ##                       DATASTORE.region   env_broad_scale env_local_scale  env_medium Experiment geo_loc_name_country geo_loc_name_country_continent geo_loc_name         HOST Instrument
    ## 1 gs.us-east1,ncbi.public,s3.us-east-1 human oral cavity     Tongue swab Tongue swab  DRX318033                Japan                           Asia        Japan Homo sapiens  Sequel II
    ## 2 gs.us-east1,ncbi.public,s3.us-east-1 human oral cavity     Tongue swab Tongue swab  DRX318034                Japan                           Asia        Japan Homo sapiens  Sequel II
    ## 3 gs.us-east1,ncbi.public,s3.us-east-1 human oral cavity     Tongue swab Tongue swab  DRX318035                Japan                           Asia        Japan Homo sapiens  Sequel II
    ## 4 gs.us-east1,ncbi.public,s3.us-east-1 human oral cavity     Tongue swab Tongue swab  DRX318036                Japan                           Asia        Japan Homo sapiens  Sequel II
    ## 5 gs.us-east1,ncbi.public,s3.us-east-1 human oral cavity     Tongue swab Tongue swab  DRX318037                Japan                           Asia        Japan Homo sapiens  Sequel II
    ## 6 gs.us-east1,ncbi.public,s3.us-east-1 human oral cavity     Tongue swab Tongue swab  DRX318038                Japan                           Asia        Japan Homo sapiens  Sequel II
    ##      lat_lon Library.Name LibraryLayout LibrarySelection LibrarySource              Organism    Platform                                              project_name          ReleaseDate
    ## 1 33 N 130 E           NA        SINGLE              PCR       GENOMIC human oral metagenome PACBIO_SMRT Characteristics of oral microbiota in mothers and infants 2021-12-10T00:00:00Z
    ## 2 33 N 130 E           NA        SINGLE              PCR       GENOMIC human oral metagenome PACBIO_SMRT Characteristics of oral microbiota in mothers and infants 2021-12-10T00:00:00Z
    ## 3 33 N 130 E           NA        SINGLE              PCR       GENOMIC human oral metagenome PACBIO_SMRT Characteristics of oral microbiota in mothers and infants 2021-12-10T00:00:00Z
    ## 4 33 N 130 E           NA        SINGLE              PCR       GENOMIC human oral metagenome PACBIO_SMRT Characteristics of oral microbiota in mothers and infants 2021-12-10T00:00:00Z
    ## 5 33 N 130 E           NA        SINGLE              PCR       GENOMIC human oral metagenome PACBIO_SMRT Characteristics of oral microbiota in mothers and infants 2021-12-10T00:00:00Z
    ## 6 33 N 130 E           NA        SINGLE              PCR       GENOMIC human oral metagenome PACBIO_SMRT Characteristics of oral microbiota in mothers and infants 2021-12-10T00:00:00Z
    ##            create_date version  Sample.Name                                    Sample_name SRA.Study
    ## 1 2021-12-10T14:29:00Z       1 SAMD00423397 16S amplicon analysis of oral microbiome FH001 DRP007969
    ## 2 2021-12-10T14:26:00Z       1 SAMD00423398 16S amplicon analysis of oral microbiome FH002 DRP007969
    ## 3 2021-12-10T14:30:00Z       1 SAMD00423399 16S amplicon analysis of oral microbiome FH003 DRP007969
    ## 4 2021-12-10T14:56:00Z       1 SAMD00423400 16S amplicon analysis of oral microbiome FH004 DRP007969
    ## 5 2021-12-10T14:24:00Z       1 SAMD00423401 16S amplicon analysis of oral microbiome FH005 DRP007969
    ## 6 2021-12-10T14:22:00Z       1 SAMD00423402 16S amplicon analysis of oral microbiome FH006 DRP007969

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
  sample_name = initial_meta$Sample_name,
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  public_time = initial_meta$ReleaseDate,
  Site = initial_meta$env_medium,
  Sex = NA,
  Age = NA,
  BMI = NA,
  Region = NA,
  PCR_primers = NA,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$Status <- "H"
metadata$title <- "Oral Microbiome Transmission and Infant Feeding Habits"
metadata$Disease_study <- "NO"
metadata$matched_samples <- "NO"
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 892"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 892"

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

    ##   BioProject       Run                                    sample_name Age BMI Sex Region  Platform Country        Site          public_time PCR_primers
    ## 1 PRJDB12588 DRR329029 16S amplicon analysis of oral microbiome FH001  NA  NA  NA     NA Sequel II   Japan Tongue swab 2021-12-10T00:00:00Z          NA
    ## 2 PRJDB12588 DRR329030 16S amplicon analysis of oral microbiome FH002  NA  NA  NA     NA Sequel II   Japan Tongue swab 2021-12-10T00:00:00Z          NA
    ## 3 PRJDB12588 DRR329031 16S amplicon analysis of oral microbiome FH003  NA  NA  NA     NA Sequel II   Japan Tongue swab 2021-12-10T00:00:00Z          NA
    ## 4 PRJDB12588 DRR329032 16S amplicon analysis of oral microbiome FH004  NA  NA  NA     NA Sequel II   Japan Tongue swab 2021-12-10T00:00:00Z          NA
    ## 5 PRJDB12588 DRR329033 16S amplicon analysis of oral microbiome FH005  NA  NA  NA     NA Sequel II   Japan Tongue swab 2021-12-10T00:00:00Z          NA
    ## 6 PRJDB12588 DRR329034 16S amplicon analysis of oral microbiome FH006  NA  NA  NA     NA Sequel II   Japan Tongue swab 2021-12-10T00:00:00Z          NA
    ##                                                    title Disease_study Status matched_samples
    ## 1 Oral Microbiome Transmission and Infant Feeding Habits            NO      H              NO
    ## 2 Oral Microbiome Transmission and Infant Feeding Habits            NO      H              NO
    ## 3 Oral Microbiome Transmission and Infant Feeding Habits            NO      H              NO
    ## 4 Oral Microbiome Transmission and Infant Feeding Habits            NO      H              NO
    ## 5 Oral Microbiome Transmission and Infant Feeding Habits            NO      H              NO
    ## 6 Oral Microbiome Transmission and Infant Feeding Habits            NO      H              NO

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
