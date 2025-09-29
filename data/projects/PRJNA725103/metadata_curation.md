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
setwd("D:/studylearning/董彪师兄项目/任务/张辰飞/PRJNA725103")

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

    ##           Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample
    ## 1 SRR14325659   AMPLICON        598 48967753 PRJNA725103 SAMN18869930
    ## 2 SRR14325660   AMPLICON        598 30280792 PRJNA725103 SAMN18869930
    ## 3 SRR14325661   AMPLICON        597 24285221 PRJNA725103 SAMN18869930
    ## 4 SRR14325662   AMPLICON        598 52630644 PRJNA725103 SAMN18869930
    ## 5 SRR14325663   AMPLICON        598 52176368 PRJNA725103 SAMN18869930
    ## 6 SRR14325664   AMPLICON        598 67502544 PRJNA725103 SAMN18869929
    ##                BioSampleModel    Bytes Center.Name Collection_Date Consent
    ## 1 Metagenome or environmental 29579767     FISABIO      2017-11-01  public
    ## 2 Metagenome or environmental 17352318     FISABIO      2017-11-01  public
    ## 3 Metagenome or environmental 14426617     FISABIO      2017-11-01  public
    ## 4 Metagenome or environmental 31210045     FISABIO      2017-11-01  public
    ## 5 Metagenome or environmental 30542696     FISABIO      2017-11-01  public
    ## 6 Metagenome or environmental 40037028     FISABIO      2017-08-01  public
    ##   DATASTORE.filetype DATASTORE.provider                     DATASTORE.region
    ## 1   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 2   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 3   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 4   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 5   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 6   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ##    Experiment geo_loc_name_country geo_loc_name_country_continent
    ## 1 SRX10680749       United Kingdom                         Europe
    ## 2 SRX10680748       United Kingdom                         Europe
    ## 3 SRX10680747       United Kingdom                         Europe
    ## 4 SRX10680746       United Kingdom                         Europe
    ## 5 SRX10680745       United Kingdom                         Europe
    ## 6 SRX10680744       United Kingdom                         Europe
    ##                                                                geo_loc_name
    ## 1 United Kingdom: Scotland\\, Glasgow Dental Hospital and School\\, Glasgow
    ## 2 United Kingdom: Scotland\\, Glasgow Dental Hospital and School\\, Glasgow
    ## 3 United Kingdom: Scotland\\, Glasgow Dental Hospital and School\\, Glasgow
    ## 4 United Kingdom: Scotland\\, Glasgow Dental Hospital and School\\, Glasgow
    ## 5 United Kingdom: Scotland\\, Glasgow Dental Hospital and School\\, Glasgow
    ## 6 United Kingdom: Scotland\\, Glasgow Dental Hospital and School\\, Glasgow
    ##           HOST     Instrument   isolation_source                  lat_lon
    ## 1 Homo sapiens Illumina MiSeq Subgingival plaque 55.8656865 N 4.2658059 W
    ## 2 Homo sapiens Illumina MiSeq Subgingival plaque 55.8656865 N 4.2658059 W
    ## 3 Homo sapiens Illumina MiSeq Subgingival plaque 55.8656865 N 4.2658059 W
    ## 4 Homo sapiens Illumina MiSeq Subgingival plaque 55.8656865 N 4.2658059 W
    ## 5 Homo sapiens Illumina MiSeq Subgingival plaque 55.8656865 N 4.2658059 W
    ## 6 Homo sapiens Illumina MiSeq Subgingival plaque 55.8656865 N 4.2658059 W
    ##            Library.Name LibraryLayout LibrarySelection LibrarySource
    ## 1 54_R1_and_R2.fastq.gz        PAIRED              PCR   METAGENOMIC
    ## 2 53_R1_and_R2.fastq.gz        PAIRED              PCR   METAGENOMIC
    ## 3 52_R1_and_R2.fastq.gz        PAIRED              PCR   METAGENOMIC
    ## 4 51_R1_and_R2.fastq.gz        PAIRED              PCR   METAGENOMIC
    ## 5 50_R1_and_R2.fastq.gz        PAIRED              PCR   METAGENOMIC
    ## 6  8_R1_and_R2.fastq.gz        PAIRED              PCR   METAGENOMIC
    ##           Organism Platform          ReleaseDate          create_date version
    ## 1 human metagenome ILLUMINA 2021-04-27T00:00:00Z 2021-04-27T04:28:00Z       1
    ## 2 human metagenome ILLUMINA 2021-04-27T00:00:00Z 2021-04-27T04:28:00Z       1
    ## 3 human metagenome ILLUMINA 2021-04-27T00:00:00Z 2021-04-27T04:28:00Z       1
    ## 4 human metagenome ILLUMINA 2021-04-27T00:00:00Z 2021-04-27T04:28:00Z       1
    ## 5 human metagenome ILLUMINA 2021-04-27T00:00:00Z 2021-04-27T04:29:00Z       1
    ## 6 human metagenome ILLUMINA 2021-04-27T00:00:00Z 2021-04-27T04:28:00Z       1
    ##                                                           Sample.Name SRA.Study
    ## 1     Subgingival_plaque_periodontitis_perio2.0_day90_after_treatment SRP316367
    ## 2     Subgingival_plaque_periodontitis_perio2.0_day90_after_treatment SRP316367
    ## 3     Subgingival_plaque_periodontitis_perio2.0_day90_after_treatment SRP316367
    ## 4     Subgingival_plaque_periodontitis_perio2.0_day90_after_treatment SRP316367
    ## 5     Subgingival_plaque_periodontitis_perio2.0_day90_after_treatment SRP316367
    ## 6 Subgingival_plaque_periodontitis_perio2.0_baseline_before_treatment SRP316367

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
metadata$PCR_primers <- "341F-805R"
metadata$title <- "Mechanical biofilm disruption causes microbial and immunological shifts in periodontitis patients"
metadata$Disease_study <- "YES"
# Assign disease status based on sample_name patterns
metadata$Status <- initial_meta$Sample.Name
metadata <- metadata %>% 
  mutate(Status = case_when(
    grepl("before",Status) ~ "PD_before_treatment",
    grepl("after",Status) ~ "PD_after_treatment"
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
    ##    1_R1_and_R2.fastq.gz   10_R1_and_R2.fastq.gz   11_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   12_R1_and_R2.fastq.gz   13_R1_and_R2.fastq.gz   14_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   15_R1_and_R2.fastq.gz   16_R1_and_R2.fastq.gz   17_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   18_R1_and_R2.fastq.gz   19_R1_and_R2.fastq.gz    2_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   20_R1_and_R2.fastq.gz   21_R1_and_R2.fastq.gz   22_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   23_R1_and_R2.fastq.gz   24_R1_and_R2.fastq.gz   25_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   26_R1_and_R2.fastq.gz   27_R1_and_R2.fastq.gz   28_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   29_R1_and_R2.fastq.gz   30_R1_and_R2.fastq.gz   31_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   32_R1_and_R2.fastq.gz   33_R1_and_R2.fastq.gz   34_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   35_R1_and_R2.fastq.gz   36_R1_and_R2.fastq.gz   37_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   38_R1_and_R2.fastq.gz   39_R1_and_R2.fastq.gz    4_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   40_R1_and_R2.fastq.gz   41_R1_and_R2.fastq.gz   43_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   45_R1_and_R2.fastq.gz   46_R1_and_R2.fastq.gz   47_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   48_R1_and_R2.fastq.gz   49_R1_and_R2.fastq.gz   50_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   51_R1_and_R2.fastq.gz   52_R1_and_R2.fastq.gz   53_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   54_R1_and_R2.fastq.gz   55_R1_and_R2.fastq.gz   56_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   57_R1_and_R2.fastq.gz   58_R1_and_R2.fastq.gz   59_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##    6_R1_and_R2.fastq.gz   60_R1_and_R2.fastq.gz   61_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   62_R1_and_R2.fastq.gz   63_R1_and_R2.fastq.gz   64_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   65_R1_and_R2.fastq.gz   66_R1_and_R2.fastq.gz   67_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   68_R1_and_R2.fastq.gz   69_R1_and_R2.fastq.gz    7_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   70_R1_and_R2.fastq.gz   71_R1_and_R2.fastq.gz   72_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   73_R1_and_R2.fastq.gz   74_R1_and_R2.fastq.gz   75_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   76_R1_and_R2.fastq.gz   77_R1_and_R2.fastq.gz   78_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##    8_R1_and_R2.fastq.gz    9_R1_and_R2.fastq.gz    P141_R1_and_R2.fastq 
    ##                       1                       1                       1 
    ##    P142_R1_and_R2.fastq    P143_R1_and_R2.fastq    P144_R1_and_R2.fastq 
    ##                       1                       1                       1 
    ##    P145_R1_and_R2.fastq P241_R1_and_R2.fastq.gz P242_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ## P243_R1_and_R2.fastq.gz P244_R1_and_R2.fastq.gz P245_R1_and_R2.fastq.gz 
    ##                       1                       1                       1

``` r
table(metadata$Site)
```

    ## 
    ## Subgingival plaque 
    ##                 84

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
    ##    1_R1_and_R2.fastq.gz   10_R1_and_R2.fastq.gz   11_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   12_R1_and_R2.fastq.gz   13_R1_and_R2.fastq.gz   14_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   15_R1_and_R2.fastq.gz   16_R1_and_R2.fastq.gz   17_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   18_R1_and_R2.fastq.gz   19_R1_and_R2.fastq.gz    2_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   20_R1_and_R2.fastq.gz   21_R1_and_R2.fastq.gz   22_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   23_R1_and_R2.fastq.gz   24_R1_and_R2.fastq.gz   25_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   26_R1_and_R2.fastq.gz   27_R1_and_R2.fastq.gz   28_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   29_R1_and_R2.fastq.gz   30_R1_and_R2.fastq.gz   31_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   32_R1_and_R2.fastq.gz   33_R1_and_R2.fastq.gz   34_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   35_R1_and_R2.fastq.gz   36_R1_and_R2.fastq.gz   37_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   38_R1_and_R2.fastq.gz   39_R1_and_R2.fastq.gz    4_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   40_R1_and_R2.fastq.gz   41_R1_and_R2.fastq.gz   43_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   45_R1_and_R2.fastq.gz   46_R1_and_R2.fastq.gz   47_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   48_R1_and_R2.fastq.gz   49_R1_and_R2.fastq.gz   50_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   51_R1_and_R2.fastq.gz   52_R1_and_R2.fastq.gz   53_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   54_R1_and_R2.fastq.gz   55_R1_and_R2.fastq.gz   56_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   57_R1_and_R2.fastq.gz   58_R1_and_R2.fastq.gz   59_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##    6_R1_and_R2.fastq.gz   60_R1_and_R2.fastq.gz   61_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   62_R1_and_R2.fastq.gz   63_R1_and_R2.fastq.gz   64_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   65_R1_and_R2.fastq.gz   66_R1_and_R2.fastq.gz   67_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   68_R1_and_R2.fastq.gz   69_R1_and_R2.fastq.gz    7_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   70_R1_and_R2.fastq.gz   71_R1_and_R2.fastq.gz   72_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   73_R1_and_R2.fastq.gz   74_R1_and_R2.fastq.gz   75_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   76_R1_and_R2.fastq.gz   77_R1_and_R2.fastq.gz   78_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##    8_R1_and_R2.fastq.gz    9_R1_and_R2.fastq.gz    P141_R1_and_R2.fastq 
    ##                       1                       1                       1 
    ##    P142_R1_and_R2.fastq    P143_R1_and_R2.fastq    P144_R1_and_R2.fastq 
    ##                       1                       1                       1 
    ##    P145_R1_and_R2.fastq P241_R1_and_R2.fastq.gz P242_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ## P243_R1_and_R2.fastq.gz P244_R1_and_R2.fastq.gz P245_R1_and_R2.fastq.gz 
    ##                       1                       1                       1

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
    ## 84

``` r
# Identify samples for removal: single samples excluding Subgingival plaque
single_samples <- names(sample_counts[sample_counts == 1])
samples_to_remove <- metadata$sample_name %in% single_samples & !(metadata$Site %in% c("Subgingival plaque"))

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
    ##    1_R1_and_R2.fastq.gz   10_R1_and_R2.fastq.gz   11_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   12_R1_and_R2.fastq.gz   13_R1_and_R2.fastq.gz   14_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   15_R1_and_R2.fastq.gz   16_R1_and_R2.fastq.gz   17_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   18_R1_and_R2.fastq.gz   19_R1_and_R2.fastq.gz    2_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   20_R1_and_R2.fastq.gz   21_R1_and_R2.fastq.gz   22_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   23_R1_and_R2.fastq.gz   24_R1_and_R2.fastq.gz   25_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   26_R1_and_R2.fastq.gz   27_R1_and_R2.fastq.gz   28_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   29_R1_and_R2.fastq.gz   30_R1_and_R2.fastq.gz   31_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   32_R1_and_R2.fastq.gz   33_R1_and_R2.fastq.gz   34_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   35_R1_and_R2.fastq.gz   36_R1_and_R2.fastq.gz   37_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   38_R1_and_R2.fastq.gz   39_R1_and_R2.fastq.gz    4_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   40_R1_and_R2.fastq.gz   41_R1_and_R2.fastq.gz   43_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   45_R1_and_R2.fastq.gz   46_R1_and_R2.fastq.gz   47_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   48_R1_and_R2.fastq.gz   49_R1_and_R2.fastq.gz   50_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   51_R1_and_R2.fastq.gz   52_R1_and_R2.fastq.gz   53_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   54_R1_and_R2.fastq.gz   55_R1_and_R2.fastq.gz   56_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   57_R1_and_R2.fastq.gz   58_R1_and_R2.fastq.gz   59_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##    6_R1_and_R2.fastq.gz   60_R1_and_R2.fastq.gz   61_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   62_R1_and_R2.fastq.gz   63_R1_and_R2.fastq.gz   64_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   65_R1_and_R2.fastq.gz   66_R1_and_R2.fastq.gz   67_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   68_R1_and_R2.fastq.gz   69_R1_and_R2.fastq.gz    7_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   70_R1_and_R2.fastq.gz   71_R1_and_R2.fastq.gz   72_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   73_R1_and_R2.fastq.gz   74_R1_and_R2.fastq.gz   75_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##   76_R1_and_R2.fastq.gz   77_R1_and_R2.fastq.gz   78_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ##    8_R1_and_R2.fastq.gz    9_R1_and_R2.fastq.gz    P141_R1_and_R2.fastq 
    ##                       1                       1                       1 
    ##    P142_R1_and_R2.fastq    P143_R1_and_R2.fastq    P144_R1_and_R2.fastq 
    ##                       1                       1                       1 
    ##    P145_R1_and_R2.fastq P241_R1_and_R2.fastq.gz P242_R1_and_R2.fastq.gz 
    ##                       1                       1                       1 
    ## P243_R1_and_R2.fastq.gz P244_R1_and_R2.fastq.gz P245_R1_and_R2.fastq.gz 
    ##                       1                       1                       1

``` r
table(metadata_clean$Site)
```

    ## 
    ## Subgingival plaque 
    ##                 84

``` r
table(metadata_clean$matched_samples)
```

    ## 
    ## NO 
    ## 84

``` r
# Use cleaned dataset for subsequent analysis
metadata <- metadata_clean
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 84"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 84"

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

    ##    BioProject         Run           sample_name Age BMI Sex Region
    ## 1 PRJNA725103 SRR14325659 54_R1_and_R2.fastq.gz  NA  NA  NA  V3-V4
    ## 2 PRJNA725103 SRR14325660 53_R1_and_R2.fastq.gz  NA  NA  NA  V3-V4
    ## 3 PRJNA725103 SRR14325661 52_R1_and_R2.fastq.gz  NA  NA  NA  V3-V4
    ## 4 PRJNA725103 SRR14325662 51_R1_and_R2.fastq.gz  NA  NA  NA  V3-V4
    ## 5 PRJNA725103 SRR14325663 50_R1_and_R2.fastq.gz  NA  NA  NA  V3-V4
    ## 6 PRJNA725103 SRR14325664  8_R1_and_R2.fastq.gz  NA  NA  NA  V3-V4
    ##         Platform        Country               Site          public_time
    ## 1 Illumina MiSeq United Kingdom Subgingival plaque 2021-04-27T00:00:00Z
    ## 2 Illumina MiSeq United Kingdom Subgingival plaque 2021-04-27T00:00:00Z
    ## 3 Illumina MiSeq United Kingdom Subgingival plaque 2021-04-27T00:00:00Z
    ## 4 Illumina MiSeq United Kingdom Subgingival plaque 2021-04-27T00:00:00Z
    ## 5 Illumina MiSeq United Kingdom Subgingival plaque 2021-04-27T00:00:00Z
    ## 6 Illumina MiSeq United Kingdom Subgingival plaque 2021-04-27T00:00:00Z
    ##   PCR_primers
    ## 1   341F-805R
    ## 2   341F-805R
    ## 3   341F-805R
    ## 4   341F-805R
    ## 5   341F-805R
    ## 6   341F-805R
    ##                                                                                               title
    ## 1 Mechanical biofilm disruption causes microbial and immunological shifts in periodontitis patients
    ## 2 Mechanical biofilm disruption causes microbial and immunological shifts in periodontitis patients
    ## 3 Mechanical biofilm disruption causes microbial and immunological shifts in periodontitis patients
    ## 4 Mechanical biofilm disruption causes microbial and immunological shifts in periodontitis patients
    ## 5 Mechanical biofilm disruption causes microbial and immunological shifts in periodontitis patients
    ## 6 Mechanical biofilm disruption causes microbial and immunological shifts in periodontitis patients
    ##   Disease_study              Status matched_samples
    ## 1           YES  PD_after_treatment              NO
    ## 2           YES  PD_after_treatment              NO
    ## 3           YES  PD_after_treatment              NO
    ## 4           YES  PD_after_treatment              NO
    ## 5           YES  PD_after_treatment              NO
    ## 6           YES PD_before_treatment              NO

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
