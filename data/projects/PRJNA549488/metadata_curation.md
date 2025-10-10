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
setwd("/home/hby/oral_project/PRJNA549488")

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

    ##          Run       Assay.Type AvgSpotLen    Bases  BioProject    BioSample              BioSampleModel    Bytes CaseControl Center.Name Collection_Date Consent DATASTORE.filetype
    ## 1 SRR9321177 Targeted-Capture        501 53504436 PRJNA549488 SAMN12087176 Metagenome or environmental 35128648     Control         NCI              NA  public   fastq,run.zq,sra
    ## 2 SRR9321178 Targeted-Capture        500 50982893 PRJNA549488 SAMN12087175 Metagenome or environmental 32953658     Control         NCI              NA  public   fastq,run.zq,sra
    ## 3 SRR9321179 Targeted-Capture        500 61418685 PRJNA549488 SAMN12087170 Metagenome or environmental 39970393     Control         NCI              NA  public   fastq,run.zq,sra
    ## 4 SRR9321180 Targeted-Capture        500 72829273 PRJNA549488 SAMN12087169 Metagenome or environmental 47159847        Case         NCI              NA  public   fastq,run.zq,sra
    ## 5 SRR9321181 Targeted-Capture        500 60707673 PRJNA549488 SAMN12087172 Metagenome or environmental 38103524        Case         NCI              NA  public   fastq,run.zq,sra
    ## 6 SRR9321182 Targeted-Capture        500 59190089 PRJNA549488 SAMN12087171 Metagenome or environmental 38209738     Control         NCI              NA  public   fastq,run.zq,sra
    ##   DATASTORE.provider                     DATASTORE.region Experiment geo_loc_name_country geo_loc_name_country_continent geo_loc_name         HOST     Instrument lat_lon
    ## 1         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX6088606         uncalculated                   uncalculated           NA Homo sapiens Illumina MiSeq      NA
    ## 2         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX6088605         uncalculated                   uncalculated           NA Homo sapiens Illumina MiSeq      NA
    ## 3         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX6088604         uncalculated                   uncalculated           NA Homo sapiens Illumina MiSeq      NA
    ## 4         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX6088603         uncalculated                   uncalculated           NA Homo sapiens Illumina MiSeq      NA
    ## 5         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX6088602         uncalculated                   uncalculated           NA Homo sapiens Illumina MiSeq      NA
    ## 6         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 SRX6088601         uncalculated                   uncalculated           NA Homo sapiens Illumina MiSeq      NA
    ##   lib_const_method Library_Name LibraryLayout LibrarySelection LibrarySource              Organism Platform          ReleaseDate          create_date version Sample.Name    seq_methods
    ## 1       paired-end     SC086788        PAIRED              PCR       GENOMIC human oral metagenome ILLUMINA 2019-10-24T00:00:00Z 2019-06-18T16:25:00Z       1    SC086788 Illumina MiSeq
    ## 2       paired-end     SC086787        PAIRED              PCR       GENOMIC human oral metagenome ILLUMINA 2019-10-24T00:00:00Z 2019-06-18T16:25:00Z       1    SC086787 Illumina MiSeq
    ## 3       paired-end     SC086778        PAIRED              PCR       GENOMIC human oral metagenome ILLUMINA 2019-10-24T00:00:00Z 2019-06-18T16:25:00Z       1    SC086778 Illumina MiSeq
    ## 4       paired-end     SC086777        PAIRED              PCR       GENOMIC human oral metagenome ILLUMINA 2019-10-24T00:00:00Z 2019-06-18T16:25:00Z       1    SC086777 Illumina MiSeq
    ## 5       paired-end     SC086784        PAIRED              PCR       GENOMIC human oral metagenome ILLUMINA 2019-10-24T00:00:00Z 2019-06-18T16:25:00Z       1    SC086784 Illumina MiSeq
    ## 6       paired-end     SC086783        PAIRED              PCR       GENOMIC human oral metagenome ILLUMINA 2019-10-24T00:00:00Z 2019-06-18T16:25:00Z       1    SC086783 Illumina MiSeq
    ##   SRA.Study Target_Gene
    ## 1 SRP201785    16S rRNA
    ## 2 SRP201785    16S rRNA
    ## 3 SRP201785    16S rRNA
    ## 4 SRP201785    16S rRNA
    ## 5 SRP201785    16S rRNA
    ## 6 SRP201785    16S rRNA

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
  Status = initial_meta$CaseControl,
  public_time = initial_meta$ReleaseDate,
  Sex = NA,
  Age = NA,
  BMI = NA
) %>%
  mutate(
    Status = case_when(
      Status == "Control" ~ "H",
      Status == "Case" ~ "pancreatic cancer",
      TRUE ~ as.character(Status)
    )
  )
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$Country <- "Iran"
metadata$Site <- "saliva"
metadata$Region <- "V4"
metadata$PCR_primers <- "515F-806R"
metadata$title <- "Cigarette Smoking and Opium Use in Relation to the Oral Microbiota in Iran"
metadata$Disease_study <- "YES"
metadata$matched_samples <- "NO"
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 558"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 558"

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
    ## 1 PRJNA549488 SRR9321177    SC086788  NA  NA  NA     V4 Illumina MiSeq    Iran saliva 2019-10-24T00:00:00Z   515F-806R
    ## 2 PRJNA549488 SRR9321178    SC086787  NA  NA  NA     V4 Illumina MiSeq    Iran saliva 2019-10-24T00:00:00Z   515F-806R
    ## 3 PRJNA549488 SRR9321179    SC086778  NA  NA  NA     V4 Illumina MiSeq    Iran saliva 2019-10-24T00:00:00Z   515F-806R
    ## 4 PRJNA549488 SRR9321180    SC086777  NA  NA  NA     V4 Illumina MiSeq    Iran saliva 2019-10-24T00:00:00Z   515F-806R
    ## 5 PRJNA549488 SRR9321181    SC086784  NA  NA  NA     V4 Illumina MiSeq    Iran saliva 2019-10-24T00:00:00Z   515F-806R
    ## 6 PRJNA549488 SRR9321182    SC086783  NA  NA  NA     V4 Illumina MiSeq    Iran saliva 2019-10-24T00:00:00Z   515F-806R
    ##                                                                        title Disease_study            Status matched_samples
    ## 1 Cigarette Smoking and Opium Use in Relation to the Oral Microbiota in Iran           YES                 H              NO
    ## 2 Cigarette Smoking and Opium Use in Relation to the Oral Microbiota in Iran           YES                 H              NO
    ## 3 Cigarette Smoking and Opium Use in Relation to the Oral Microbiota in Iran           YES                 H              NO
    ## 4 Cigarette Smoking and Opium Use in Relation to the Oral Microbiota in Iran           YES pancreatic cancer              NO
    ## 5 Cigarette Smoking and Opium Use in Relation to the Oral Microbiota in Iran           YES pancreatic cancer              NO
    ## 6 Cigarette Smoking and Opium Use in Relation to the Oral Microbiota in Iran           YES                 H              NO

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
