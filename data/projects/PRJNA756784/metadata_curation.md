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
setwd("D:/studylearning/董彪师兄项目/任务/张辰飞/PRJNA756784")

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

    ##           Run Assay.Type AvgSpotLen     Bases  BioProject    BioSample
    ## 1 SRR15560887   AMPLICON        300 195696000 PRJNA756784 SAMN20926076
    ## 2 SRR15560888   AMPLICON        300 202776600 PRJNA756784 SAMN20926075
    ## 3 SRR15560889   AMPLICON        300 171236400 PRJNA756784 SAMN20926074
    ## 4 SRR15560890   AMPLICON        300 186323100 PRJNA756784 SAMN20926073
    ## 5 SRR15560891   AMPLICON        300 223891200 PRJNA756784 SAMN20926126
    ## 6 SRR15560892   AMPLICON        300 227489100 PRJNA756784 SAMN20926125
    ##                         BioSampleModel    Bytes
    ## 1 MIMS.me,MIGS/MIMS/MIMARKS.human-oral 81167346
    ## 2 MIMS.me,MIGS/MIMS/MIMARKS.human-oral 80524214
    ## 3 MIMS.me,MIGS/MIMS/MIMARKS.human-oral 70755584
    ## 4 MIMS.me,MIGS/MIMS/MIMARKS.human-oral 77023989
    ## 5 MIMS.me,MIGS/MIMS/MIMARKS.human-oral 89719844
    ## 6 MIMS.me,MIGS/MIMS/MIMARKS.human-oral 90913520
    ##                             Center.Name Collection_Date Consent
    ## 1 BABASAHEB BHIMRAO AMBEDAKR UNIVERSITY         2016-04  public
    ## 2 BABASAHEB BHIMRAO AMBEDAKR UNIVERSITY         2016-04  public
    ## 3 BABASAHEB BHIMRAO AMBEDAKR UNIVERSITY         2016-04  public
    ## 4 BABASAHEB BHIMRAO AMBEDAKR UNIVERSITY         2016-04  public
    ## 5 BABASAHEB BHIMRAO AMBEDAKR UNIVERSITY         2016-03  public
    ## 6 BABASAHEB BHIMRAO AMBEDAKR UNIVERSITY         2016-02  public
    ##   DATASTORE.filetype DATASTORE.provider                     DATASTORE.region
    ## 1   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 2   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 3   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 4   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 5   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 6   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ##   env_broad_scale env_local_scale           env_medium  Experiment
    ## 1 Oral Microbiome      Mouth Wash Sterile Saline Water SRX11859069
    ## 2 Oral Microbiome      Mouth Wash Sterile Saline Water SRX11859068
    ## 3 Oral Microbiome      Mouth Wash Sterile Saline Water SRX11859067
    ## 4 Oral Microbiome      Mouth Wash Sterile Saline Water SRX11859066
    ## 5 Oral Microbiome      Mouth Wash Sterile Saline Water SRX11859065
    ## 6 Oral Microbiome      Mouth Wash Sterile Saline Water SRX11859064
    ##   geo_loc_name_country geo_loc_name_country_continent geo_loc_name         HOST
    ## 1                India                           Asia        India Homo sapiens
    ## 2                India                           Asia        India Homo sapiens
    ## 3                India                           Asia        India Homo sapiens
    ## 4                India                           Asia        India Homo sapiens
    ## 5                India                           Asia        India Homo sapiens
    ## 6                India                           Asia        India Homo sapiens
    ##            Instrument             lat_lon Library.Name LibraryLayout
    ## 1 Illumina HiSeq 2500     28.61 N 77.03 E     Lib_HI11        PAIRED
    ## 2 Illumina HiSeq 2500     28.61 N 77.03 E     Lib_HI10        PAIRED
    ## 3 Illumina HiSeq 2500     28.61 N 77.03 E      Lib_HI9        PAIRED
    ## 4 Illumina HiSeq 2500     28.61 N 77.03 E      Lib_HI8        PAIRED
    ## 5 Illumina HiSeq 2500 28.5659 N 77.2111 E     Lib_CP20        PAIRED
    ## 6 Illumina HiSeq 2500 28.5659 N 77.2111 E     Lib_CP19        PAIRED
    ##   LibrarySelection LibrarySource              Organism Platform
    ## 1              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 2              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 3              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 4              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 5              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 6              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ##            ReleaseDate          create_date version Sample.Name SRA.Study
    ## 1 2021-08-22T00:00:00Z 2021-08-22T19:38:00Z       1        HI11 SRP333603
    ## 2 2021-08-22T00:00:00Z 2021-08-22T19:38:00Z       1        HI10 SRP333603
    ## 3 2021-08-22T00:00:00Z 2021-08-22T19:38:00Z       1         HI9 SRP333603
    ## 4 2021-08-22T00:00:00Z 2021-08-22T19:38:00Z       1         HI8 SRP333603
    ## 5 2021-08-22T00:00:00Z 2021-08-22T19:38:00Z       1        CP21 SRP333603
    ## 6 2021-08-22T00:00:00Z 2021-08-22T19:38:00Z       1        CP19 SRP333603
    ##   hypervariable_region sample_id        sample_site sample_type sequence_length
    ## 1            V3 region      HI11 New Delhi\\, India   Oral Wash           150bp
    ## 2            V3 region      HI10 New Delhi\\, India   Oral Wash           150bp
    ## 3            V3 region       HI9 New Delhi\\, India   Oral Wash           150bp
    ## 4            V3 region       HI8 New Delhi\\, India   Oral Wash           150bp
    ## 5            V3 region      CP20 New Delhi\\, India   Oral Wash           150bp
    ## 6            V3 region      CP19 New Delhi\\, India   Oral Wash           150bp
    ##   sequencing_platform
    ## 1      Illumina Miseq
    ## 2      Illumina Miseq
    ## 3      Illumina Miseq
    ## 4      Illumina Miseq
    ## 5      Illumina Miseq
    ## 6      Illumina Miseq

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
  sample_name = initial_meta$Sample.Name,
  Age = NA,
  BMI = NA,
  Sex= NA,
  Region = "V3",
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  Site = "oral wash",
  public_time = initial_meta$ReleaseDate,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$PCR_primers <- NA_character_
metadata$title <- "Comparative and analytical characterization of the oral bacteriome of smokeless tobacco users with oral squamous cell carcinoma"
metadata$Disease_study <- "YES"
# Assign disease status based on sample_name patterns
metadata <- metadata %>% 
  mutate(Status = case_when(
    grepl("H|T",sample_name) ~ "HC",
    grepl("CP",sample_name) ~ "OSCC"
  ))


#extract subject id by remove prefix of sample_name
#metadata$sample_name <- gsub("T|CR|CS|C|N|S", "", metadata$sample_name)

# Check initial data distribution
cat("Initial sample distribution:\n")
```

    ## Initial sample distribution:

``` r
table(metadata$sample_name)
```

    ## 
    ##   CP1  CP10  CP11  CP12  CP13  CP14  CP15  CP16  CP17  CP18  CP19   CP2  CP21 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##   CP3   CP4   CP5   CP6   CP7   CP8   CP9   HI1  HI10  HI11  HI12  HI13  HI14 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##  HI16  HI17  HI18  HI19   HI2  HI20  HI21  HI22   HI3   HI4   HI5   HI7   HI8 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##   HI9  HTC1 HTC10 HTC11 HTC12 HTC13 HTC14 HTC15 HTC16 HTC17 HTC19  HTC2 HTC20 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##  HTC3  HTC4  HTC5  HTC6  HTC7  HTC8  HTC9  TC18 
    ##     1     1     1     1     1     1     1     1

``` r
table(metadata$Site)
```

    ## 
    ## oral wash 
    ##        60

``` r
# Calculate sampling frequency per subject
sample_counts <- table(metadata$sample_name)
cat("\nNumber of samples per subject:\n")
```

    ## 
    ## Number of samples per subject:

``` r
print(sample_counts)
```

    ## 
    ##   CP1  CP10  CP11  CP12  CP13  CP14  CP15  CP16  CP17  CP18  CP19   CP2  CP21 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##   CP3   CP4   CP5   CP6   CP7   CP8   CP9   HI1  HI10  HI11  HI12  HI13  HI14 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##  HI16  HI17  HI18  HI19   HI2  HI20  HI21  HI22   HI3   HI4   HI5   HI7   HI8 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##   HI9  HTC1 HTC10 HTC11 HTC12 HTC13 HTC14 HTC15 HTC16 HTC17 HTC19  HTC2 HTC20 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##  HTC3  HTC4  HTC5  HTC6  HTC7  HTC8  HTC9  TC18 
    ##     1     1     1     1     1     1     1     1

``` r
# Create matched_samples column
metadata$matched_samples <- ifelse(metadata$sample_name %in% names(sample_counts[sample_counts >= 2]), "YES", "NO")

cat("\nMatched samples distribution:\n")
```

    ## 
    ## Matched samples distribution:

``` r
table(metadata$matched_samples)
```

    ## 
    ## NO 
    ## 60

``` r
# Identify samples for removal: single samples excluding oral wash
single_samples <- names(sample_counts[sample_counts == 1])
samples_to_remove <- metadata$sample_name %in% single_samples & !(metadata$Site %in% c("oral wash"))

cat("\nNumber of samples to be removed:", sum(samples_to_remove), "\n")
```

    ## 
    ## Number of samples to be removed: 0

``` r
# Remove identified samples
metadata_clean <- metadata[!samples_to_remove, ]

cat("\nFinal sample distribution after cleaning:\n")
```

    ## 
    ## Final sample distribution after cleaning:

``` r
table(metadata_clean$sample_name)
```

    ## 
    ##   CP1  CP10  CP11  CP12  CP13  CP14  CP15  CP16  CP17  CP18  CP19   CP2  CP21 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##   CP3   CP4   CP5   CP6   CP7   CP8   CP9   HI1  HI10  HI11  HI12  HI13  HI14 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##  HI16  HI17  HI18  HI19   HI2  HI20  HI21  HI22   HI3   HI4   HI5   HI7   HI8 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##   HI9  HTC1 HTC10 HTC11 HTC12 HTC13 HTC14 HTC15 HTC16 HTC17 HTC19  HTC2 HTC20 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##  HTC3  HTC4  HTC5  HTC6  HTC7  HTC8  HTC9  TC18 
    ##     1     1     1     1     1     1     1     1

``` r
table(metadata_clean$Site)
```

    ## 
    ## oral wash 
    ##        60

``` r
table(metadata_clean$matched_samples)
```

    ## 
    ## NO 
    ## 60

``` r
# Use cleaned dataset for subsequent analysis
metadata <- metadata_clean
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 60"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 60"

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

    ##    BioProject         Run sample_name Age BMI Sex Region            Platform
    ## 1 PRJNA756784 SRR15560887        HI11  NA  NA  NA     V3 Illumina HiSeq 2500
    ## 2 PRJNA756784 SRR15560888        HI10  NA  NA  NA     V3 Illumina HiSeq 2500
    ## 3 PRJNA756784 SRR15560889         HI9  NA  NA  NA     V3 Illumina HiSeq 2500
    ## 4 PRJNA756784 SRR15560890         HI8  NA  NA  NA     V3 Illumina HiSeq 2500
    ## 5 PRJNA756784 SRR15560891        CP21  NA  NA  NA     V3 Illumina HiSeq 2500
    ## 6 PRJNA756784 SRR15560892        CP19  NA  NA  NA     V3 Illumina HiSeq 2500
    ##   Country      Site          public_time PCR_primers
    ## 1   India oral wash 2021-08-22T00:00:00Z        <NA>
    ## 2   India oral wash 2021-08-22T00:00:00Z        <NA>
    ## 3   India oral wash 2021-08-22T00:00:00Z        <NA>
    ## 4   India oral wash 2021-08-22T00:00:00Z        <NA>
    ## 5   India oral wash 2021-08-22T00:00:00Z        <NA>
    ## 6   India oral wash 2021-08-22T00:00:00Z        <NA>
    ##                                                                                                                             title
    ## 1 Comparative and analytical characterization of the oral bacteriome of smokeless tobacco users with oral squamous cell carcinoma
    ## 2 Comparative and analytical characterization of the oral bacteriome of smokeless tobacco users with oral squamous cell carcinoma
    ## 3 Comparative and analytical characterization of the oral bacteriome of smokeless tobacco users with oral squamous cell carcinoma
    ## 4 Comparative and analytical characterization of the oral bacteriome of smokeless tobacco users with oral squamous cell carcinoma
    ## 5 Comparative and analytical characterization of the oral bacteriome of smokeless tobacco users with oral squamous cell carcinoma
    ## 6 Comparative and analytical characterization of the oral bacteriome of smokeless tobacco users with oral squamous cell carcinoma
    ##   Disease_study Status matched_samples
    ## 1           YES     HC              NO
    ## 2           YES     HC              NO
    ## 3           YES     HC              NO
    ## 4           YES     HC              NO
    ## 5           YES   OSCC              NO
    ## 6           YES   OSCC              NO

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
