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
setwd("D:/studylearning/董彪师兄项目/任务/张辰飞/PRJNA666972")

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
    ## 1 SRR12771660   AMPLICON        498 13614448 PRJNA666972 SAMN16361678
    ## 2 SRR12771661   AMPLICON        495 27770355 PRJNA666972 SAMN16361684
    ## 3 SRR12771662   AMPLICON        496 21871683 PRJNA666972 SAMN16361686
    ## 4 SRR12771663   AMPLICON        495 30704065 PRJNA666972 SAMN16361688
    ## 5 SRR12771664   AMPLICON        498 16962632 PRJNA666972 SAMN16361694
    ## 6 SRR12771665   AMPLICON        495 27844091 PRJNA666972 SAMN16361696
    ##                BioSampleModel    Bytes        Center.Name Collection_Date
    ## 1 Metagenome or environmental  7560185 SICHUAN UNIVERSITY      2020-05-01
    ## 2 Metagenome or environmental 17180056 SICHUAN UNIVERSITY      2020-05-01
    ## 3 Metagenome or environmental 12963266 SICHUAN UNIVERSITY      2020-05-01
    ## 4 Metagenome or environmental 17548026 SICHUAN UNIVERSITY      2020-05-01
    ## 5 Metagenome or environmental  9402500 SICHUAN UNIVERSITY      2020-05-01
    ## 6 Metagenome or environmental 15947570 SICHUAN UNIVERSITY      2020-05-01
    ##   Consent DATASTORE.filetype DATASTORE.provider
    ## 1  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 2  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 3  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 4  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 6  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                       DATASTORE.region Experiment geo_loc_name_country
    ## 1 gs.us-east1,ncbi.public,s3.us-east-1 SRX9242330                China
    ## 2 gs.us-east1,ncbi.public,s3.us-east-1 SRX9242329                China
    ## 3 gs.us-east1,ncbi.public,s3.us-east-1 SRX9242328                China
    ## 4 gs.us-east1,ncbi.public,s3.us-east-1 SRX9242327                China
    ## 5 gs.us-east1,ncbi.public,s3.us-east-1 SRX9242326                China
    ## 6 gs.us-east1,ncbi.public,s3.us-east-1 SRX9242325                China
    ##   geo_loc_name_country_continent   geo_loc_name     Instrument isolation_source
    ## 1                           Asia China: Chengdu Illumina MiSeq    HNC18969067_1
    ## 2                           Asia China: Chengdu Illumina MiSeq    HNC19003579_1
    ## 3                           Asia China: Chengdu Illumina MiSeq    HNC18751267_1
    ## 4                           Asia China: Chengdu Illumina MiSeq    HNC02554362_1
    ## 5                           Asia China: Chengdu Illumina MiSeq    HNC18676892_1
    ## 6                           Asia China: Chengdu Illumina MiSeq    HNC17760127_1
    ##            lat_lon Library.Name LibraryLayout LibrarySelection LibrarySource
    ## 1 30.67 N 104.08 E         L161        PAIRED              PCR   METAGENOMIC
    ## 2 30.67 N 104.08 E         L167        PAIRED              PCR   METAGENOMIC
    ## 3 30.67 N 104.08 E         L169        PAIRED              PCR   METAGENOMIC
    ## 4 30.67 N 104.08 E         L171        PAIRED              PCR   METAGENOMIC
    ## 5 30.67 N 104.08 E         L177        PAIRED              PCR   METAGENOMIC
    ## 6 30.67 N 104.08 E         L179        PAIRED              PCR   METAGENOMIC
    ##                  Organism Platform          ReleaseDate          create_date
    ## 1 human saliva metagenome ILLUMINA 2024-10-04T00:00:00Z 2020-10-05T08:24:00Z
    ## 2 human saliva metagenome ILLUMINA 2024-10-04T00:00:00Z 2020-10-05T08:24:00Z
    ## 3 human saliva metagenome ILLUMINA 2024-10-04T00:00:00Z 2020-10-05T08:24:00Z
    ## 4 human saliva metagenome ILLUMINA 2024-10-04T00:00:00Z 2020-10-05T08:24:00Z
    ## 5 human saliva metagenome ILLUMINA 2024-10-04T00:00:00Z 2020-10-05T08:24:00Z
    ## 6 human saliva metagenome ILLUMINA 2024-10-04T00:00:00Z 2020-10-05T08:24:00Z
    ##   version  Sample.Name SRA.Study
    ## 1       1 S62_R1.fastq SRP286350
    ## 2       1 S66_R1.fastq SRP286350
    ## 3       1 S69_R1.fastq SRP286350
    ## 4       1 S73_R1.fastq SRP286350
    ## 5       1 S79_R1.fastq SRP286350
    ## 6       1 S80_R1.fastq SRP286350

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
  Sex= "M",
  Region = "V3-V4",
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  Site = "saliva",
  public_time = initial_meta$ReleaseDate,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$PCR_primers <- "341F-806R"
metadata$title <- "Study on the Salivary Microbial Alteration of Men With Head and Neck Cancer and Its Relationship With Symptoms in Southwest China"
metadata$Disease_study <- "YES"
# Assign disease status based on sample_name patterns
metadata <- metadata %>% 
  mutate(Status = case_when(
    grepl("N",sample_name) ~ "HC",
    grepl("S",sample_name) ~ "HNC"
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
    ## N101_R1.fastq N102_R1.fastq N103_R1.fastq N104_R1.fastq N105_R1.fastq 
    ##             1             1             1             1             1 
    ## N107_R1.fastq N109_R1.fastq N111_R1.fastq  N12_R1.fastq  N13_R1.fastq 
    ##             1             1             1             1             1 
    ##  N17_R1.fastq  N19_R1.fastq  N20_R1.fastq  N21_R1.fastq  N23_R1.fastq 
    ##             1             1             1             1             1 
    ##  N26_R1.fastq  N27_R1.fastq   N3_R1.fastq  N31_R1.fastq  N35_R1.fastq 
    ##             1             1             1             1             1 
    ##  N36_R1.fastq  N37_R1.fastq  N38_R1.fastq  N41_R1.fastq  N42_R1.fastq 
    ##             1             1             1             1             1 
    ##  N43_R1.fastq  N44_R1.fastq  N45_R1.fastq  N46_R1.fastq  N47_R1.fastq 
    ##             1             1             1             1             1 
    ##  N48_R1.fastq   N5_R1.fastq  N50_R1.fastq  N52_R1.fastq  N53_R1.fastq 
    ##             1             1             1             1             1 
    ##  N54_R1.fastq  N55_R1.fastq  N56_R1.fastq  N57_R1.fastq  N62_R1.fastq 
    ##             1             1             1             1             1 
    ##  N63_R1.fastq  N64_R1.fastq  N65_R1.fastq  N66_R1.fastq  N67_R1.fastq 
    ##             1             1             1             1             1 
    ##  N68_R1.fastq  N69_R1.fastq  N70_R1.fastq  N71_R1.fastq  N72_R1.fastq 
    ##             1             1             1             1             1 
    ##  N73_R1.fastq  N74_R1.fastq  N75_R1.fastq  N76_R1.fastq  N80_R1.fastq 
    ##             1             1             1             1             1 
    ##  N81_R1.fastq  N84_R1.fastq  N85_R1.fastq  N86_R1.fastq  N87_R1.fastq 
    ##             1             1             1             1             1 
    ##  N88_R1.fastq  N89_R1.fastq  N90_R1.fastq  N92_R1.fastq S100_R1.fastq 
    ##             1             1             1             1             1 
    ## S104_R1.fastq S106_R1.fastq S111_R1.fastq S112_R1.fastq S114_R1.fastq 
    ##             1             1             1             1             1 
    ##  S12_R1.fastq S126_R1.fastq S127_R1.fastq S128_R1.fastq S150_R1.fastq 
    ##             1             1             1             1             1 
    ## S157_R1.fastq S168_R1.fastq S169_R1.fastq  S17_R1.fastq S170_R1.fastq 
    ##             1             1             1             1             1 
    ## S172_R1.fastq S174_R1.fastq  S18_R1.fastq S187_R1.fastq  S19_R1.fastq 
    ##             1             1             1             1             1 
    ##  S20_R1.fastq S205_R1.fastq S206_R1.fastq S208_R1.fastq  S21_R1.fastq 
    ##             1             1             1             1             1 
    ## S213_R1.fastq S214_R1.fastq S219_R1.fastq  S26_R1.fastq  S27_R1.fastq 
    ##             1             1             1             1             1 
    ##  S29_R1.fastq  S30_R1.fastq  S31_R1.fastq  S49_R1.fastq  S50_R1.fastq 
    ##             1             1             1             1             1 
    ##  S51_R1.fastq  S59_R1.fastq  S62_R1.fastq  S63_R1.fastq  S64_R1.fastq 
    ##             1             1             1             1             1 
    ##  S66_R1.fastq  S69_R1.fastq  S73_R1.fastq  S75_R1.fastq  S77_R1.fastq 
    ##             1             1             1             1             1 
    ##  S79_R1.fastq  S80_R1.fastq  S82_R1.fastq  S83_R1.fastq  S86_R1.fastq 
    ##             1             1             1             1             1 
    ##  S87_R1.fastq  S89_R1.fastq   S9_R1.fastq  S98_R1.fastq  S99_R1.fastq 
    ##             1             1             1             1             1

``` r
table(metadata$Site)
```

    ## 
    ## saliva 
    ##    120

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
    ## N101_R1.fastq N102_R1.fastq N103_R1.fastq N104_R1.fastq N105_R1.fastq 
    ##             1             1             1             1             1 
    ## N107_R1.fastq N109_R1.fastq N111_R1.fastq  N12_R1.fastq  N13_R1.fastq 
    ##             1             1             1             1             1 
    ##  N17_R1.fastq  N19_R1.fastq  N20_R1.fastq  N21_R1.fastq  N23_R1.fastq 
    ##             1             1             1             1             1 
    ##  N26_R1.fastq  N27_R1.fastq   N3_R1.fastq  N31_R1.fastq  N35_R1.fastq 
    ##             1             1             1             1             1 
    ##  N36_R1.fastq  N37_R1.fastq  N38_R1.fastq  N41_R1.fastq  N42_R1.fastq 
    ##             1             1             1             1             1 
    ##  N43_R1.fastq  N44_R1.fastq  N45_R1.fastq  N46_R1.fastq  N47_R1.fastq 
    ##             1             1             1             1             1 
    ##  N48_R1.fastq   N5_R1.fastq  N50_R1.fastq  N52_R1.fastq  N53_R1.fastq 
    ##             1             1             1             1             1 
    ##  N54_R1.fastq  N55_R1.fastq  N56_R1.fastq  N57_R1.fastq  N62_R1.fastq 
    ##             1             1             1             1             1 
    ##  N63_R1.fastq  N64_R1.fastq  N65_R1.fastq  N66_R1.fastq  N67_R1.fastq 
    ##             1             1             1             1             1 
    ##  N68_R1.fastq  N69_R1.fastq  N70_R1.fastq  N71_R1.fastq  N72_R1.fastq 
    ##             1             1             1             1             1 
    ##  N73_R1.fastq  N74_R1.fastq  N75_R1.fastq  N76_R1.fastq  N80_R1.fastq 
    ##             1             1             1             1             1 
    ##  N81_R1.fastq  N84_R1.fastq  N85_R1.fastq  N86_R1.fastq  N87_R1.fastq 
    ##             1             1             1             1             1 
    ##  N88_R1.fastq  N89_R1.fastq  N90_R1.fastq  N92_R1.fastq S100_R1.fastq 
    ##             1             1             1             1             1 
    ## S104_R1.fastq S106_R1.fastq S111_R1.fastq S112_R1.fastq S114_R1.fastq 
    ##             1             1             1             1             1 
    ##  S12_R1.fastq S126_R1.fastq S127_R1.fastq S128_R1.fastq S150_R1.fastq 
    ##             1             1             1             1             1 
    ## S157_R1.fastq S168_R1.fastq S169_R1.fastq  S17_R1.fastq S170_R1.fastq 
    ##             1             1             1             1             1 
    ## S172_R1.fastq S174_R1.fastq  S18_R1.fastq S187_R1.fastq  S19_R1.fastq 
    ##             1             1             1             1             1 
    ##  S20_R1.fastq S205_R1.fastq S206_R1.fastq S208_R1.fastq  S21_R1.fastq 
    ##             1             1             1             1             1 
    ## S213_R1.fastq S214_R1.fastq S219_R1.fastq  S26_R1.fastq  S27_R1.fastq 
    ##             1             1             1             1             1 
    ##  S29_R1.fastq  S30_R1.fastq  S31_R1.fastq  S49_R1.fastq  S50_R1.fastq 
    ##             1             1             1             1             1 
    ##  S51_R1.fastq  S59_R1.fastq  S62_R1.fastq  S63_R1.fastq  S64_R1.fastq 
    ##             1             1             1             1             1 
    ##  S66_R1.fastq  S69_R1.fastq  S73_R1.fastq  S75_R1.fastq  S77_R1.fastq 
    ##             1             1             1             1             1 
    ##  S79_R1.fastq  S80_R1.fastq  S82_R1.fastq  S83_R1.fastq  S86_R1.fastq 
    ##             1             1             1             1             1 
    ##  S87_R1.fastq  S89_R1.fastq   S9_R1.fastq  S98_R1.fastq  S99_R1.fastq 
    ##             1             1             1             1             1

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
    ##  NO 
    ## 120

``` r
# Identify samples for removal: single samples excluding saliva 
single_samples <- names(sample_counts[sample_counts == 1])
samples_to_remove <- metadata$sample_name %in% single_samples & !(metadata$Site %in% c("saliva"))

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
    ## N101_R1.fastq N102_R1.fastq N103_R1.fastq N104_R1.fastq N105_R1.fastq 
    ##             1             1             1             1             1 
    ## N107_R1.fastq N109_R1.fastq N111_R1.fastq  N12_R1.fastq  N13_R1.fastq 
    ##             1             1             1             1             1 
    ##  N17_R1.fastq  N19_R1.fastq  N20_R1.fastq  N21_R1.fastq  N23_R1.fastq 
    ##             1             1             1             1             1 
    ##  N26_R1.fastq  N27_R1.fastq   N3_R1.fastq  N31_R1.fastq  N35_R1.fastq 
    ##             1             1             1             1             1 
    ##  N36_R1.fastq  N37_R1.fastq  N38_R1.fastq  N41_R1.fastq  N42_R1.fastq 
    ##             1             1             1             1             1 
    ##  N43_R1.fastq  N44_R1.fastq  N45_R1.fastq  N46_R1.fastq  N47_R1.fastq 
    ##             1             1             1             1             1 
    ##  N48_R1.fastq   N5_R1.fastq  N50_R1.fastq  N52_R1.fastq  N53_R1.fastq 
    ##             1             1             1             1             1 
    ##  N54_R1.fastq  N55_R1.fastq  N56_R1.fastq  N57_R1.fastq  N62_R1.fastq 
    ##             1             1             1             1             1 
    ##  N63_R1.fastq  N64_R1.fastq  N65_R1.fastq  N66_R1.fastq  N67_R1.fastq 
    ##             1             1             1             1             1 
    ##  N68_R1.fastq  N69_R1.fastq  N70_R1.fastq  N71_R1.fastq  N72_R1.fastq 
    ##             1             1             1             1             1 
    ##  N73_R1.fastq  N74_R1.fastq  N75_R1.fastq  N76_R1.fastq  N80_R1.fastq 
    ##             1             1             1             1             1 
    ##  N81_R1.fastq  N84_R1.fastq  N85_R1.fastq  N86_R1.fastq  N87_R1.fastq 
    ##             1             1             1             1             1 
    ##  N88_R1.fastq  N89_R1.fastq  N90_R1.fastq  N92_R1.fastq S100_R1.fastq 
    ##             1             1             1             1             1 
    ## S104_R1.fastq S106_R1.fastq S111_R1.fastq S112_R1.fastq S114_R1.fastq 
    ##             1             1             1             1             1 
    ##  S12_R1.fastq S126_R1.fastq S127_R1.fastq S128_R1.fastq S150_R1.fastq 
    ##             1             1             1             1             1 
    ## S157_R1.fastq S168_R1.fastq S169_R1.fastq  S17_R1.fastq S170_R1.fastq 
    ##             1             1             1             1             1 
    ## S172_R1.fastq S174_R1.fastq  S18_R1.fastq S187_R1.fastq  S19_R1.fastq 
    ##             1             1             1             1             1 
    ##  S20_R1.fastq S205_R1.fastq S206_R1.fastq S208_R1.fastq  S21_R1.fastq 
    ##             1             1             1             1             1 
    ## S213_R1.fastq S214_R1.fastq S219_R1.fastq  S26_R1.fastq  S27_R1.fastq 
    ##             1             1             1             1             1 
    ##  S29_R1.fastq  S30_R1.fastq  S31_R1.fastq  S49_R1.fastq  S50_R1.fastq 
    ##             1             1             1             1             1 
    ##  S51_R1.fastq  S59_R1.fastq  S62_R1.fastq  S63_R1.fastq  S64_R1.fastq 
    ##             1             1             1             1             1 
    ##  S66_R1.fastq  S69_R1.fastq  S73_R1.fastq  S75_R1.fastq  S77_R1.fastq 
    ##             1             1             1             1             1 
    ##  S79_R1.fastq  S80_R1.fastq  S82_R1.fastq  S83_R1.fastq  S86_R1.fastq 
    ##             1             1             1             1             1 
    ##  S87_R1.fastq  S89_R1.fastq   S9_R1.fastq  S98_R1.fastq  S99_R1.fastq 
    ##             1             1             1             1             1

``` r
table(metadata_clean$Site)
```

    ## 
    ## saliva 
    ##    120

``` r
table(metadata_clean$matched_samples)
```

    ## 
    ##  NO 
    ## 120

``` r
# Use cleaned dataset for subsequent analysis
metadata <- metadata_clean
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 120"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 120"

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

    ##    BioProject         Run  sample_name Age BMI Sex Region       Platform
    ## 1 PRJNA666972 SRR12771660 S62_R1.fastq  NA  NA   M  V3-V4 Illumina MiSeq
    ## 2 PRJNA666972 SRR12771661 S66_R1.fastq  NA  NA   M  V3-V4 Illumina MiSeq
    ## 3 PRJNA666972 SRR12771662 S69_R1.fastq  NA  NA   M  V3-V4 Illumina MiSeq
    ## 4 PRJNA666972 SRR12771663 S73_R1.fastq  NA  NA   M  V3-V4 Illumina MiSeq
    ## 5 PRJNA666972 SRR12771664 S79_R1.fastq  NA  NA   M  V3-V4 Illumina MiSeq
    ## 6 PRJNA666972 SRR12771665 S80_R1.fastq  NA  NA   M  V3-V4 Illumina MiSeq
    ##   Country   Site          public_time PCR_primers
    ## 1   China saliva 2024-10-04T00:00:00Z   341F-806R
    ## 2   China saliva 2024-10-04T00:00:00Z   341F-806R
    ## 3   China saliva 2024-10-04T00:00:00Z   341F-806R
    ## 4   China saliva 2024-10-04T00:00:00Z   341F-806R
    ## 5   China saliva 2024-10-04T00:00:00Z   341F-806R
    ## 6   China saliva 2024-10-04T00:00:00Z   341F-806R
    ##                                                                                                                               title
    ## 1 Study on the Salivary Microbial Alteration of Men With Head and Neck Cancer and Its Relationship With Symptoms in Southwest China
    ## 2 Study on the Salivary Microbial Alteration of Men With Head and Neck Cancer and Its Relationship With Symptoms in Southwest China
    ## 3 Study on the Salivary Microbial Alteration of Men With Head and Neck Cancer and Its Relationship With Symptoms in Southwest China
    ## 4 Study on the Salivary Microbial Alteration of Men With Head and Neck Cancer and Its Relationship With Symptoms in Southwest China
    ## 5 Study on the Salivary Microbial Alteration of Men With Head and Neck Cancer and Its Relationship With Symptoms in Southwest China
    ## 6 Study on the Salivary Microbial Alteration of Men With Head and Neck Cancer and Its Relationship With Symptoms in Southwest China
    ##   Disease_study Status matched_samples
    ## 1           YES    HNC              NO
    ## 2           YES    HNC              NO
    ## 3           YES    HNC              NO
    ## 4           YES    HNC              NO
    ## 5           YES    HNC              NO
    ## 6           YES    HNC              NO

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
