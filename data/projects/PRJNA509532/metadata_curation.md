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
setwd("D:/oral_project/PRJNA509532")

# Read the initial metadata files
initial_meta<-read.csv("SraRunTable.csv",header = T)

# Combine the two metadata files
#meta <- rbind(initial_meta, initial_meta1)

#Remove entries where Organism is "gut metagenome" (non-oral samples)
#initial_meta <- subset(initial_meta,host_body_product=="oral cavity")

# Check the number of samples containing "HT" in Sample.Name (e.g., healthy controls)
#nrow(subset(meta, grepl("HT", Sample.Name)))

# Display the first few rows of the metadata
head(initial_meta)
```

    ##          Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample
    ## 1 SRR8303833   AMPLICON        602 26771542 PRJNA509532 SAMN10584284
    ## 2 SRR8303834   AMPLICON        602 24641064 PRJNA509532 SAMN10584285
    ## 3 SRR8303835   AMPLICON        602 18165350 PRJNA509532 SAMN10584278
    ## 4 SRR8303836   AMPLICON        602 26794418 PRJNA509532 SAMN10584279
    ## 5 SRR8303837   AMPLICON        602 26053958 PRJNA509532 SAMN10584276
    ## 6 SRR8303838   AMPLICON        602 24148026 PRJNA509532 SAMN10584277
    ##                BioSampleModel    Bytes                             Center.Name
    ## 1 Metagenome or environmental 13895091 PEKING UNIVERSITY HEALTH SCIENCE CENTER
    ## 2 Metagenome or environmental 12071479 PEKING UNIVERSITY HEALTH SCIENCE CENTER
    ## 3 Metagenome or environmental  9554432 PEKING UNIVERSITY HEALTH SCIENCE CENTER
    ## 4 Metagenome or environmental 13756910 PEKING UNIVERSITY HEALTH SCIENCE CENTER
    ## 5 Metagenome or environmental 13619276 PEKING UNIVERSITY HEALTH SCIENCE CENTER
    ## 6 Metagenome or environmental 12514228 PEKING UNIVERSITY HEALTH SCIENCE CENTER
    ##   Collection_Date Consent DATASTORE.filetype DATASTORE.provider
    ## 1            2017  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 2            2017  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 3            2017  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 4            2017  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5            2017  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 6            2017  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                       DATASTORE.region Experiment geo_loc_name_country
    ## 1 gs.us-east1,ncbi.public,s3.us-east-1 SRX5118479                China
    ## 2 gs.us-east1,ncbi.public,s3.us-east-1 SRX5118478                China
    ## 3 gs.us-east1,ncbi.public,s3.us-east-1 SRX5118477                China
    ## 4 gs.us-east1,ncbi.public,s3.us-east-1 SRX5118476                China
    ## 5 gs.us-east1,ncbi.public,s3.us-east-1 SRX5118475                China
    ## 6 gs.us-east1,ncbi.public,s3.us-east-1 SRX5118474                China
    ##   geo_loc_name_country_continent geo_loc_name                             HOST
    ## 1                           Asia        China chronic periodontitis patient3.5
    ## 2                           Asia        China chronic periodontitis patient4.1
    ## 3                           Asia        China chronic periodontitis patient2.5
    ## 4                           Asia        China chronic periodontitis patient2.6
    ## 5                           Asia        China chronic periodontitis patient2.3
    ## 6                           Asia        China chronic periodontitis patient2.4
    ##       Instrument isolation_source       lat_lon Library.Name LibraryLayout
    ## 1 Illumina MiSeq        ChP_Sub21 not collected        C3_c2        PAIRED
    ## 2 Illumina MiSeq           ChP_B6 not collected         C4_B        PAIRED
    ## 3 Illumina MiSeq        ChP_Sub16 not collected        C2_b4        PAIRED
    ## 4 Illumina MiSeq        ChP_Sub17 not collected        C2_c2        PAIRED
    ## 5 Illumina MiSeq        ChP_Sub14 not collected        C2_b2        PAIRED
    ## 6 Illumina MiSeq        ChP_Sub15 not collected        C2_b3        PAIRED
    ##   LibrarySelection LibrarySource              Organism Platform
    ## 1              PCR         OTHER human oral metagenome ILLUMINA
    ## 2              PCR         OTHER human oral metagenome ILLUMINA
    ## 3              PCR         OTHER human oral metagenome ILLUMINA
    ## 4              PCR         OTHER human oral metagenome ILLUMINA
    ## 5              PCR         OTHER human oral metagenome ILLUMINA
    ## 6              PCR         OTHER human oral metagenome ILLUMINA
    ##            ReleaseDate          create_date version Sample.Name SRA.Study
    ## 1 2020-01-12T00:00:00Z 2018-12-12T12:06:00Z       1       C3_c2 SRP173111
    ## 2 2020-01-12T00:00:00Z 2018-12-12T12:06:00Z       1        C4_B SRP173111
    ## 3 2020-01-12T00:00:00Z 2018-12-12T12:06:00Z       1       C2_b4 SRP173111
    ## 4 2020-01-12T00:00:00Z 2018-12-12T12:06:00Z       1       C2_c2 SRP173111
    ## 5 2020-01-12T00:00:00Z 2018-12-12T12:06:00Z       1       C2_b2 SRP173111
    ## 6 2020-01-12T00:00:00Z 2018-12-12T12:07:00Z       1       C2_b3 SRP173111

# Step 2: Create a Uniform Metadata Template

The uniform metadata should include the following columns:

1.  **BioProject**: Project ID (e.g., “PRJNA339212”)
2.  **Run**: Sample ID (e.g., “ERRXXXXXX” or “SRRXXXXXX”)
3.  **sample_name**: Original sample name (if subject id is available,
    use subject id) submitted by the authors (not the database-assigned
    ID)
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
16. **matched_samples**: whether the study have matched samples from
    oral and other sites (“YES”,“NO”)

``` r
# Create a uniform metadata dataframe
metadata <- data.frame(
  BioProject = initial_meta$BioProject,
  Run = initial_meta$Run,
  sample_name = initial_meta$isolation_source,
  Age = NA,
  BMI = NA,
  Sex=NA,
  Region = "V4",
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

metadata$title <- "Alterations and Correlations in Microbial Community and Metabolome Characteristics in Generalized Aggressive Periodontitis"
metadata$Disease_study <- "YES"
# Assign disease status based on sample_name patterns
metadata <- metadata %>%
  mutate(
    Site = case_when(
      grepl("_B", Site) ~ "Gingival crevicular fluid",  
      grepl("_Sub", Site) ~ "Subgingival plaque"        
    ),
    Status = case_when(
      grepl("AgP", sample_name) ~ "aggressive periodontitis",
      grepl("ChP", sample_name) ~ "chronic periodontitis", 
      grepl("HP", sample_name) ~ "health subject"
    )
  )

metadata$matched_samples="NO"
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 145"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 145"

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

    ##    BioProject        Run sample_name Age BMI Sex Region       Platform Country
    ## 1 PRJNA509532 SRR8303833   ChP_Sub21  NA  NA  NA     V4 Illumina MiSeq   China
    ## 2 PRJNA509532 SRR8303834      ChP_B6  NA  NA  NA     V4 Illumina MiSeq   China
    ## 3 PRJNA509532 SRR8303835   ChP_Sub16  NA  NA  NA     V4 Illumina MiSeq   China
    ## 4 PRJNA509532 SRR8303836   ChP_Sub17  NA  NA  NA     V4 Illumina MiSeq   China
    ## 5 PRJNA509532 SRR8303837   ChP_Sub14  NA  NA  NA     V4 Illumina MiSeq   China
    ## 6 PRJNA509532 SRR8303838   ChP_Sub15  NA  NA  NA     V4 Illumina MiSeq   China
    ##                        Site          public_time PCR_primers
    ## 1        Subgingival plaque 2020-01-12T00:00:00Z          NA
    ## 2 Gingival crevicular fluid 2020-01-12T00:00:00Z          NA
    ## 3        Subgingival plaque 2020-01-12T00:00:00Z          NA
    ## 4        Subgingival plaque 2020-01-12T00:00:00Z          NA
    ## 5        Subgingival plaque 2020-01-12T00:00:00Z          NA
    ## 6        Subgingival plaque 2020-01-12T00:00:00Z          NA
    ##                                                                                                                        title
    ## 1 Alterations and Correlations in Microbial Community and Metabolome Characteristics in Generalized Aggressive Periodontitis
    ## 2 Alterations and Correlations in Microbial Community and Metabolome Characteristics in Generalized Aggressive Periodontitis
    ## 3 Alterations and Correlations in Microbial Community and Metabolome Characteristics in Generalized Aggressive Periodontitis
    ## 4 Alterations and Correlations in Microbial Community and Metabolome Characteristics in Generalized Aggressive Periodontitis
    ## 5 Alterations and Correlations in Microbial Community and Metabolome Characteristics in Generalized Aggressive Periodontitis
    ## 6 Alterations and Correlations in Microbial Community and Metabolome Characteristics in Generalized Aggressive Periodontitis
    ##   Disease_study                Status matched_samples
    ## 1           YES chronic periodontitis              NO
    ## 2           YES chronic periodontitis              NO
    ## 3           YES chronic periodontitis              NO
    ## 4           YES chronic periodontitis              NO
    ## 5           YES chronic periodontitis              NO
    ## 6           YES chronic periodontitis              NO

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
