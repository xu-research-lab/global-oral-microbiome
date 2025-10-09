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
setwd("D:/studylearning/董彪师兄项目/任务/张辰飞/PRJNA681486")

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
    ## 1 SRR13194555   AMPLICON        598  35138546 PRJNA681486 SAMN16954691
    ## 2 SRR13194556   AMPLICON        598  52194781 PRJNA681486 SAMN16954740
    ## 3 SRR13194557   AMPLICON        597  57942994 PRJNA681486 SAMN16954740
    ## 4 SRR13194558   AMPLICON        598  33690852 PRJNA681486 SAMN16954691
    ## 5 SRR13194559   AMPLICON        598 125882358 PRJNA681486 SAMN16954691
    ## 6 SRR13194560   AMPLICON        597  27333249 PRJNA681486 SAMN16954740
    ##                BioSampleModel    Bytes Center.Name Collection_Date Consent
    ## 1 Metagenome or environmental 21916539     FISABIO      2020-01-01  public
    ## 2 Metagenome or environmental 32973757     FISABIO      2020-01-01  public
    ## 3 Metagenome or environmental 36060264     FISABIO      2020-01-01  public
    ## 4 Metagenome or environmental 21366333     FISABIO      2020-01-01  public
    ## 5 Metagenome or environmental 78648832     FISABIO      2020-01-01  public
    ## 6 Metagenome or environmental 17031694     FISABIO      2020-01-01  public
    ##   DATASTORE.filetype DATASTORE.provider                     DATASTORE.region
    ## 1   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 2   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 3   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 4   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 5   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 6   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ##   Experiment geo_loc_name_country geo_loc_name_country_continent geo_loc_name
    ## 1 SRX9628763               Sweden                         Europe       Sweden
    ## 2 SRX9628762               Sweden                         Europe       Sweden
    ## 3 SRX9628761               Sweden                         Europe       Sweden
    ## 4 SRX9628760               Sweden                         Europe       Sweden
    ## 5 SRX9628759               Sweden                         Europe       Sweden
    ## 6 SRX9628758               Sweden                         Europe       Sweden
    ##       Instrument isolation_source                 lat_lon   Library.Name
    ## 1 Illumina MiSeq    dental plaque 58.298584 N 12.961619 E 10_R1.fastq.gz
    ## 2 Illumina MiSeq    dental plaque 58.298584 N 12.961619 E  9_R1.fastq.gz
    ## 3 Illumina MiSeq    dental plaque 58.298584 N 12.961619 E  8_R1.fastq.gz
    ## 4 Illumina MiSeq    dental plaque 58.298584 N 12.961619 E  7_R1.fastq.gz
    ## 5 Illumina MiSeq    dental plaque 58.298584 N 12.961619 E  6_R1.fastq.gz
    ## 6 Illumina MiSeq    dental plaque 58.298584 N 12.961619 E  5_R1.fastq.gz
    ##   LibraryLayout LibrarySelection LibrarySource              Organism Platform
    ## 1        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 2        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 3        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 4        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 5        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 6        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ##            ReleaseDate          create_date version
    ## 1 2020-12-03T00:00:00Z 2020-12-03T16:56:00Z       1
    ## 2 2020-12-03T00:00:00Z 2020-12-03T16:56:00Z       1
    ## 3 2020-12-03T00:00:00Z 2020-12-03T16:56:00Z       1
    ## 4 2020-12-03T00:00:00Z 2020-12-03T16:56:00Z       1
    ## 5 2020-12-03T00:00:00Z 2020-12-03T16:57:00Z       1
    ## 6 2020-12-03T00:00:00Z 2020-12-03T16:56:00Z       1
    ##                                  Sample.Name SRA.Study
    ## 1      Dental plaque microbiota caries group SRP295713
    ## 2 Dental plaque microbiota caries-free group SRP295713
    ## 3 Dental plaque microbiota caries-free group SRP295713
    ## 4      Dental plaque microbiota caries group SRP295713
    ## 5      Dental plaque microbiota caries group SRP295713
    ## 6 Dental plaque microbiota caries-free group SRP295713

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
  Sex= NA,
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
metadata$PCR_primers <- "Forward (CCTACGGGNGGCWGCAG) and Reverse (GACTACHVGGGTATCTAATCC) "
metadata$title <- "Bacterial Composition and Metabolomics of Dental Plaque From Adolescents"
metadata$Disease_study <- "YES"
# Assign disease status based on sample_name patterns
metadata$Status <- initial_meta$Sample.Name
metadata <- metadata %>% 
  mutate(Status = case_when(
    grepl("free", Status) ~ "HC",
    TRUE ~ "Dental_Caries"
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
    ##  1_R1.fastq.gz 10_R1.fastq.gz 11_R1.fastq.gz 12_R1.fastq.gz 13_R1.fastq.gz 
    ##              1              1              1              1              1 
    ## 14_R1.fastq.gz 15_R1.fastq.gz 16_R1.fastq.gz 17_R1.fastq.gz 18_R1.fastq.gz 
    ##              1              1              1              1              1 
    ## 19_R1.fastq.gz  2_R1.fastq.gz 20_R1.fastq.gz 21_R1.fastq.gz 22_R1.fastq.gz 
    ##              1              1              1              1              1 
    ## 23_R1.fastq.gz 24_R1.fastq.gz 25_R1.fastq.gz 26_R1.fastq.gz 27_R1.fastq.gz 
    ##              1              1              1              1              1 
    ## 28_R1.fastq.gz 29_R1.fastq.gz  3_R1.fastq.gz 30_R1.fastq.gz 31_R1.fastq.gz 
    ##              1              1              1              1              1 
    ## 32_R1.fastq.gz 33_R1.fastq.gz 34_R1.fastq.gz 35_R1.fastq.gz 36_R1.fastq.gz 
    ##              1              1              1              1              1 
    ## 37_R1.fastq.gz 38_R1.fastq.gz 39_R1.fastq.gz  4_R1.fastq.gz 40_R1.fastq.gz 
    ##              1              1              1              1              1 
    ##  5_R1.fastq.gz  6_R1.fastq.gz  7_R1.fastq.gz  8_R1.fastq.gz  9_R1.fastq.gz 
    ##              1              1              1              1              1

``` r
table(metadata$Site)
```

    ## 
    ## dental plaque 
    ##            40

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
    ##  1_R1.fastq.gz 10_R1.fastq.gz 11_R1.fastq.gz 12_R1.fastq.gz 13_R1.fastq.gz 
    ##              1              1              1              1              1 
    ## 14_R1.fastq.gz 15_R1.fastq.gz 16_R1.fastq.gz 17_R1.fastq.gz 18_R1.fastq.gz 
    ##              1              1              1              1              1 
    ## 19_R1.fastq.gz  2_R1.fastq.gz 20_R1.fastq.gz 21_R1.fastq.gz 22_R1.fastq.gz 
    ##              1              1              1              1              1 
    ## 23_R1.fastq.gz 24_R1.fastq.gz 25_R1.fastq.gz 26_R1.fastq.gz 27_R1.fastq.gz 
    ##              1              1              1              1              1 
    ## 28_R1.fastq.gz 29_R1.fastq.gz  3_R1.fastq.gz 30_R1.fastq.gz 31_R1.fastq.gz 
    ##              1              1              1              1              1 
    ## 32_R1.fastq.gz 33_R1.fastq.gz 34_R1.fastq.gz 35_R1.fastq.gz 36_R1.fastq.gz 
    ##              1              1              1              1              1 
    ## 37_R1.fastq.gz 38_R1.fastq.gz 39_R1.fastq.gz  4_R1.fastq.gz 40_R1.fastq.gz 
    ##              1              1              1              1              1 
    ##  5_R1.fastq.gz  6_R1.fastq.gz  7_R1.fastq.gz  8_R1.fastq.gz  9_R1.fastq.gz 
    ##              1              1              1              1              1

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
    ## 40

``` r
# Identify samples for removal: single samples excluding dental plaque
single_samples <- names(sample_counts[sample_counts == 1])
samples_to_remove <- metadata$sample_name %in% single_samples & !(metadata$Site %in% c("dental plaque"))

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
    ##  1_R1.fastq.gz 10_R1.fastq.gz 11_R1.fastq.gz 12_R1.fastq.gz 13_R1.fastq.gz 
    ##              1              1              1              1              1 
    ## 14_R1.fastq.gz 15_R1.fastq.gz 16_R1.fastq.gz 17_R1.fastq.gz 18_R1.fastq.gz 
    ##              1              1              1              1              1 
    ## 19_R1.fastq.gz  2_R1.fastq.gz 20_R1.fastq.gz 21_R1.fastq.gz 22_R1.fastq.gz 
    ##              1              1              1              1              1 
    ## 23_R1.fastq.gz 24_R1.fastq.gz 25_R1.fastq.gz 26_R1.fastq.gz 27_R1.fastq.gz 
    ##              1              1              1              1              1 
    ## 28_R1.fastq.gz 29_R1.fastq.gz  3_R1.fastq.gz 30_R1.fastq.gz 31_R1.fastq.gz 
    ##              1              1              1              1              1 
    ## 32_R1.fastq.gz 33_R1.fastq.gz 34_R1.fastq.gz 35_R1.fastq.gz 36_R1.fastq.gz 
    ##              1              1              1              1              1 
    ## 37_R1.fastq.gz 38_R1.fastq.gz 39_R1.fastq.gz  4_R1.fastq.gz 40_R1.fastq.gz 
    ##              1              1              1              1              1 
    ##  5_R1.fastq.gz  6_R1.fastq.gz  7_R1.fastq.gz  8_R1.fastq.gz  9_R1.fastq.gz 
    ##              1              1              1              1              1

``` r
table(metadata_clean$Site)
```

    ## 
    ## dental plaque 
    ##            40

``` r
table(metadata_clean$matched_samples)
```

    ## 
    ## NO 
    ## 40

``` r
# Use cleaned dataset for subsequent analysis
metadata <- metadata_clean
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 40"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 40"

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

    ##    BioProject         Run    sample_name Age BMI Sex Region       Platform
    ## 1 PRJNA681486 SRR13194555 10_R1.fastq.gz  NA  NA  NA  V3-V4 Illumina MiSeq
    ## 2 PRJNA681486 SRR13194556  9_R1.fastq.gz  NA  NA  NA  V3-V4 Illumina MiSeq
    ## 3 PRJNA681486 SRR13194557  8_R1.fastq.gz  NA  NA  NA  V3-V4 Illumina MiSeq
    ## 4 PRJNA681486 SRR13194558  7_R1.fastq.gz  NA  NA  NA  V3-V4 Illumina MiSeq
    ## 5 PRJNA681486 SRR13194559  6_R1.fastq.gz  NA  NA  NA  V3-V4 Illumina MiSeq
    ## 6 PRJNA681486 SRR13194560  5_R1.fastq.gz  NA  NA  NA  V3-V4 Illumina MiSeq
    ##   Country          Site          public_time
    ## 1  Sweden dental plaque 2020-12-03T00:00:00Z
    ## 2  Sweden dental plaque 2020-12-03T00:00:00Z
    ## 3  Sweden dental plaque 2020-12-03T00:00:00Z
    ## 4  Sweden dental plaque 2020-12-03T00:00:00Z
    ## 5  Sweden dental plaque 2020-12-03T00:00:00Z
    ## 6  Sweden dental plaque 2020-12-03T00:00:00Z
    ##                                                        PCR_primers
    ## 1 Forward (CCTACGGGNGGCWGCAG) and Reverse (GACTACHVGGGTATCTAATCC) 
    ## 2 Forward (CCTACGGGNGGCWGCAG) and Reverse (GACTACHVGGGTATCTAATCC) 
    ## 3 Forward (CCTACGGGNGGCWGCAG) and Reverse (GACTACHVGGGTATCTAATCC) 
    ## 4 Forward (CCTACGGGNGGCWGCAG) and Reverse (GACTACHVGGGTATCTAATCC) 
    ## 5 Forward (CCTACGGGNGGCWGCAG) and Reverse (GACTACHVGGGTATCTAATCC) 
    ## 6 Forward (CCTACGGGNGGCWGCAG) and Reverse (GACTACHVGGGTATCTAATCC) 
    ##                                                                      title
    ## 1 Bacterial Composition and Metabolomics of Dental Plaque From Adolescents
    ## 2 Bacterial Composition and Metabolomics of Dental Plaque From Adolescents
    ## 3 Bacterial Composition and Metabolomics of Dental Plaque From Adolescents
    ## 4 Bacterial Composition and Metabolomics of Dental Plaque From Adolescents
    ## 5 Bacterial Composition and Metabolomics of Dental Plaque From Adolescents
    ## 6 Bacterial Composition and Metabolomics of Dental Plaque From Adolescents
    ##   Disease_study        Status matched_samples
    ## 1           YES Dental_Caries              NO
    ## 2           YES            HC              NO
    ## 3           YES            HC              NO
    ## 4           YES Dental_Caries              NO
    ## 5           YES Dental_Caries              NO
    ## 6           YES            HC              NO

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
