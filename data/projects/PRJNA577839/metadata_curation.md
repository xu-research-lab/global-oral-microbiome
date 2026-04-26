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
setwd("D:/studylearning/董彪师兄项目/任务/修改20250922/PRJNA577839")

# Read the initial metadata files
initial_meta<-read.csv("SraRunTable.csv",header = T)

# Display the first few rows of the metadata
head(initial_meta)
```

    ##           Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample
    ## 1 SRR10294619   AMPLICON        602 21746648 PRJNA577839 SAMN13038908
    ## 2 SRR10294620   AMPLICON        602 21419762 PRJNA577839 SAMN13038907
    ## 3 SRR10294621   AMPLICON        602 36098930 PRJNA577839 SAMN13038973
    ## 4 SRR10294622   AMPLICON        602 40451992 PRJNA577839 SAMN13038972
    ## 5 SRR10294623   AMPLICON        602 51254882 PRJNA577839 SAMN13038971
    ## 6 SRR10294624   AMPLICON        602 48892032 PRJNA577839 SAMN13038970
    ##                                BioSampleModel    Bytes
    ## 1 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 15139708
    ## 2 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 14615134
    ## 3 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 25056655
    ## 4 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 27285783
    ## 5 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 35431822
    ## 6 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 33769472
    ##                                         Center.Name Collection_Date Consent
    ## 1 UNIVERSITY HOSPITAL FREIBURG, FACULTY OF MEDICINE         2016-10  public
    ## 2 UNIVERSITY HOSPITAL FREIBURG, FACULTY OF MEDICINE         2016-10  public
    ## 3 UNIVERSITY HOSPITAL FREIBURG, FACULTY OF MEDICINE         2013-12  public
    ## 4 UNIVERSITY HOSPITAL FREIBURG, FACULTY OF MEDICINE         2013-12  public
    ## 5 UNIVERSITY HOSPITAL FREIBURG, FACULTY OF MEDICINE         2013-11  public
    ## 6 UNIVERSITY HOSPITAL FREIBURG, FACULTY OF MEDICINE         2013-11  public
    ##   DATASTORE.filetype DATASTORE.provider                     DATASTORE.region
    ## 1   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 2   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 3   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 4   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 5   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 6   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ##   env_broad_scale env_local_scale env_medium Experiment geo_loc_name_country
    ## 1  human organism     oral cavity     saliva SRX7007521              Germany
    ## 2  human organism     oral cavity     saliva SRX7007520              Germany
    ## 3  human organism     oral cavity     saliva SRX7007519              Germany
    ## 4  human organism     oral cavity     saliva SRX7007518              Germany
    ## 5  human organism     oral cavity     saliva SRX7007517              Germany
    ## 6  human organism     oral cavity     saliva SRX7007516              Germany
    ##   geo_loc_name_country_continent        geo_loc_name         Host
    ## 1                         Europe Germany: Heidelberg Homo sapiens
    ## 2                         Europe Germany: Heidelberg Homo sapiens
    ## 3                         Europe Germany: Heidelberg Homo sapiens
    ## 4                         Europe Germany: Heidelberg Homo sapiens
    ## 5                         Europe Germany: Heidelberg Homo sapiens
    ## 6                         Europe Germany: Heidelberg Homo sapiens
    ##   host_tissue_sampled     Instrument        lat_lon        Library.Name
    ## 1              saliva Illumina MiSeq 49.40 N 8.69 E 17003713-F_16.fastq
    ## 2              saliva Illumina MiSeq 49.40 N 8.69 E 17003717-F_15.fastq
    ## 3              saliva Illumina MiSeq 49.40 N 8.69 E 17003752-S_56.fastq
    ## 4              saliva Illumina MiSeq 49.40 N 8.69 E 17003755-S_55.fastq
    ## 5              saliva Illumina MiSeq 49.40 N 8.69 E 17003780-S_54.fastq
    ## 6              saliva Illumina MiSeq 49.40 N 8.69 E 17003742-S_53.fastq
    ##   LibraryLayout LibrarySelection LibrarySource              Organism Platform
    ## 1        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 2        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 3        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 4        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 5        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 6        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ##            ReleaseDate          create_date version      Sample.Name SRA.Study
    ## 1 2020-02-03T00:00:00Z 2019-10-16T18:17:00Z       1  10_F_16_TG_Sn_3 SRP225957
    ## 2 2020-02-03T00:00:00Z 2019-10-16T18:17:00Z       1   9_F_15_TG_Sn_3 SRP225957
    ## 3 2020-02-03T00:00:00Z 2019-10-16T18:16:00Z       1 75_S_56_CG_nSn_0 SRP225957
    ## 4 2020-02-03T00:00:00Z 2019-10-16T18:17:00Z       1 74_S_55_CG_nSn_0 SRP225957
    ## 5 2020-02-03T00:00:00Z 2019-10-16T18:17:00Z       1 73_S_54_CG_nSn_0 SRP225957
    ## 6 2020-02-03T00:00:00Z 2019-10-16T18:17:00Z       1 72_S_53_CG_nSn_0 SRP225957

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
  Sex= NA,
  Region = "V3-V4",
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  Site = initial_meta$env_medium,
  public_time = initial_meta$ReleaseDate,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$PCR_primers <- "341F-806R"  # Update based on the publication
metadata$title <- "Influence of the long-term use of oral hygiene products containing stannous ions on the salivary microbiome – a randomized controlled trial"
metadata$Disease_study <- "NO"
metadata <- metadata %>% mutate(Status=case_when(grepl("TG",sample_name)~"TG",grepl("CG",sample_name)~"CG"))

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
    ##   1_F_01_TG_Sn_3  10_F_16_TG_Sn_3  11_F_17_TG_Sn_3  12_F_19_TG_Sn_3 
    ##                1                1                1                1 
    ##  13_F_23_TG_Sn_3  14_F_24_TG_Sn_3  15_F_25_TG_Sn_3  16_F_26_TG_Sn_3 
    ##                1                1                1                1 
    ## 17_F_32_CG_nSn_3 18_F_33_CG_nSn_3 19_F_34_CG_nSn_3   2_F_02_TG_Sn_3 
    ##                1                1                1                1 
    ## 20_F_35_CG_nSn_3 21_F_37_CG_nSn_3 22_F_38_CG_nSn_3 23_F_39_CG_nSn_3 
    ##                1                1                1                1 
    ## 24_F_40_CG_nSn_3 25_F_41_CG_nSn_3 26_F_42_CG_nSn_3 27_F_43_CG_nSn_3 
    ##                1                1                1                1 
    ## 28_F_44_CG_nSn_3 29_F_46_CG_nSn_3   3_F_04_TG_Sn_3 30_F_47_CG_nSn_3 
    ##                1                1                1                1 
    ## 31_F_48_CG_nSn_3 32_F_49_CG_nSn_3 33_F_51_CG_nSn_3 34_F_52_CG_nSn_3 
    ##                1                1                1                1 
    ## 35_F_53_CG_nSn_3 36_F_54_CG_nSn_3 37_F_55_CG_nSn_3 38_F_56_CG_nSn_3 
    ##                1                1                1                1 
    ##  39_S_01_TG_Sn_0   4_F_05_TG_Sn_3  40_S_02_TG_Sn_0  41_S_04_TG_Sn_0 
    ##                1                1                1                1 
    ##  42_S_05_TG_Sn_0  43_S_06_TG_Sn_0  44_S_08_TG_Sn_0  45_S_09_TG_Sn_0 
    ##                1                1                1                1 
    ##  46_S_12_TG_Sn_0  47_S_16_TG_Sn_0  48_S_17_TG_Sn_0  49_S_19_TG_Sn_0 
    ##                1                1                1                1 
    ##   5_F_06_TG_Sn_3  50_S_23_TG_Sn_0  51_S_24_TG_Sn_0  52_S_25_TG_Sn_0 
    ##                1                1                1                1 
    ##  53_S_26_TG_Sn_0 54_S_32_CG_nSn_0 55_S_33_CG_nSn_0 56_S_34_CG_nSn_0 
    ##                1                1                1                1 
    ## 57_S_35_CG_nSn_0 58_S_37_CG_nSn_0 59_S_38_CG_nSn_0   6_F_08_TG_Sn_3 
    ##                1                1                1                1 
    ## 60_S_39_CG_nSn_0 61_S_40_CG_nSn_0 62_S_41_CG_nSn_0 63_S_42_CG_nSn_0 
    ##                1                1                1                1 
    ## 64_S_43_CG_nSn_0 65_S_44_CG_nSn_0 66_S_46_CG_nSn_0 67_S_47_CG_nSn_0 
    ##                1                1                1                1 
    ## 68_S_48_CG_nSn_0 69_S_49_CG_nSn_0   7_F_09_TG_Sn_3 70_S_51_CG_nSn_0 
    ##                1                1                1                1 
    ## 71_S_52_CG_nSn_0 72_S_53_CG_nSn_0 73_S_54_CG_nSn_0 74_S_55_CG_nSn_0 
    ##                1                1                1                1 
    ## 75_S_56_CG_nSn_0   8_F_12_TG_Sn_3   9_F_15_TG_Sn_3 
    ##                1                1                1

``` r
table(metadata$Site)
```

    ## 
    ## saliva 
    ##     75

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
    ##   1_F_01_TG_Sn_3  10_F_16_TG_Sn_3  11_F_17_TG_Sn_3  12_F_19_TG_Sn_3 
    ##                1                1                1                1 
    ##  13_F_23_TG_Sn_3  14_F_24_TG_Sn_3  15_F_25_TG_Sn_3  16_F_26_TG_Sn_3 
    ##                1                1                1                1 
    ## 17_F_32_CG_nSn_3 18_F_33_CG_nSn_3 19_F_34_CG_nSn_3   2_F_02_TG_Sn_3 
    ##                1                1                1                1 
    ## 20_F_35_CG_nSn_3 21_F_37_CG_nSn_3 22_F_38_CG_nSn_3 23_F_39_CG_nSn_3 
    ##                1                1                1                1 
    ## 24_F_40_CG_nSn_3 25_F_41_CG_nSn_3 26_F_42_CG_nSn_3 27_F_43_CG_nSn_3 
    ##                1                1                1                1 
    ## 28_F_44_CG_nSn_3 29_F_46_CG_nSn_3   3_F_04_TG_Sn_3 30_F_47_CG_nSn_3 
    ##                1                1                1                1 
    ## 31_F_48_CG_nSn_3 32_F_49_CG_nSn_3 33_F_51_CG_nSn_3 34_F_52_CG_nSn_3 
    ##                1                1                1                1 
    ## 35_F_53_CG_nSn_3 36_F_54_CG_nSn_3 37_F_55_CG_nSn_3 38_F_56_CG_nSn_3 
    ##                1                1                1                1 
    ##  39_S_01_TG_Sn_0   4_F_05_TG_Sn_3  40_S_02_TG_Sn_0  41_S_04_TG_Sn_0 
    ##                1                1                1                1 
    ##  42_S_05_TG_Sn_0  43_S_06_TG_Sn_0  44_S_08_TG_Sn_0  45_S_09_TG_Sn_0 
    ##                1                1                1                1 
    ##  46_S_12_TG_Sn_0  47_S_16_TG_Sn_0  48_S_17_TG_Sn_0  49_S_19_TG_Sn_0 
    ##                1                1                1                1 
    ##   5_F_06_TG_Sn_3  50_S_23_TG_Sn_0  51_S_24_TG_Sn_0  52_S_25_TG_Sn_0 
    ##                1                1                1                1 
    ##  53_S_26_TG_Sn_0 54_S_32_CG_nSn_0 55_S_33_CG_nSn_0 56_S_34_CG_nSn_0 
    ##                1                1                1                1 
    ## 57_S_35_CG_nSn_0 58_S_37_CG_nSn_0 59_S_38_CG_nSn_0   6_F_08_TG_Sn_3 
    ##                1                1                1                1 
    ## 60_S_39_CG_nSn_0 61_S_40_CG_nSn_0 62_S_41_CG_nSn_0 63_S_42_CG_nSn_0 
    ##                1                1                1                1 
    ## 64_S_43_CG_nSn_0 65_S_44_CG_nSn_0 66_S_46_CG_nSn_0 67_S_47_CG_nSn_0 
    ##                1                1                1                1 
    ## 68_S_48_CG_nSn_0 69_S_49_CG_nSn_0   7_F_09_TG_Sn_3 70_S_51_CG_nSn_0 
    ##                1                1                1                1 
    ## 71_S_52_CG_nSn_0 72_S_53_CG_nSn_0 73_S_54_CG_nSn_0 74_S_55_CG_nSn_0 
    ##                1                1                1                1 
    ## 75_S_56_CG_nSn_0   8_F_12_TG_Sn_3   9_F_15_TG_Sn_3 
    ##                1                1                1

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
    ## 75

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
    ##   1_F_01_TG_Sn_3  10_F_16_TG_Sn_3  11_F_17_TG_Sn_3  12_F_19_TG_Sn_3 
    ##                1                1                1                1 
    ##  13_F_23_TG_Sn_3  14_F_24_TG_Sn_3  15_F_25_TG_Sn_3  16_F_26_TG_Sn_3 
    ##                1                1                1                1 
    ## 17_F_32_CG_nSn_3 18_F_33_CG_nSn_3 19_F_34_CG_nSn_3   2_F_02_TG_Sn_3 
    ##                1                1                1                1 
    ## 20_F_35_CG_nSn_3 21_F_37_CG_nSn_3 22_F_38_CG_nSn_3 23_F_39_CG_nSn_3 
    ##                1                1                1                1 
    ## 24_F_40_CG_nSn_3 25_F_41_CG_nSn_3 26_F_42_CG_nSn_3 27_F_43_CG_nSn_3 
    ##                1                1                1                1 
    ## 28_F_44_CG_nSn_3 29_F_46_CG_nSn_3   3_F_04_TG_Sn_3 30_F_47_CG_nSn_3 
    ##                1                1                1                1 
    ## 31_F_48_CG_nSn_3 32_F_49_CG_nSn_3 33_F_51_CG_nSn_3 34_F_52_CG_nSn_3 
    ##                1                1                1                1 
    ## 35_F_53_CG_nSn_3 36_F_54_CG_nSn_3 37_F_55_CG_nSn_3 38_F_56_CG_nSn_3 
    ##                1                1                1                1 
    ##  39_S_01_TG_Sn_0   4_F_05_TG_Sn_3  40_S_02_TG_Sn_0  41_S_04_TG_Sn_0 
    ##                1                1                1                1 
    ##  42_S_05_TG_Sn_0  43_S_06_TG_Sn_0  44_S_08_TG_Sn_0  45_S_09_TG_Sn_0 
    ##                1                1                1                1 
    ##  46_S_12_TG_Sn_0  47_S_16_TG_Sn_0  48_S_17_TG_Sn_0  49_S_19_TG_Sn_0 
    ##                1                1                1                1 
    ##   5_F_06_TG_Sn_3  50_S_23_TG_Sn_0  51_S_24_TG_Sn_0  52_S_25_TG_Sn_0 
    ##                1                1                1                1 
    ##  53_S_26_TG_Sn_0 54_S_32_CG_nSn_0 55_S_33_CG_nSn_0 56_S_34_CG_nSn_0 
    ##                1                1                1                1 
    ## 57_S_35_CG_nSn_0 58_S_37_CG_nSn_0 59_S_38_CG_nSn_0   6_F_08_TG_Sn_3 
    ##                1                1                1                1 
    ## 60_S_39_CG_nSn_0 61_S_40_CG_nSn_0 62_S_41_CG_nSn_0 63_S_42_CG_nSn_0 
    ##                1                1                1                1 
    ## 64_S_43_CG_nSn_0 65_S_44_CG_nSn_0 66_S_46_CG_nSn_0 67_S_47_CG_nSn_0 
    ##                1                1                1                1 
    ## 68_S_48_CG_nSn_0 69_S_49_CG_nSn_0   7_F_09_TG_Sn_3 70_S_51_CG_nSn_0 
    ##                1                1                1                1 
    ## 71_S_52_CG_nSn_0 72_S_53_CG_nSn_0 73_S_54_CG_nSn_0 74_S_55_CG_nSn_0 
    ##                1                1                1                1 
    ## 75_S_56_CG_nSn_0   8_F_12_TG_Sn_3   9_F_15_TG_Sn_3 
    ##                1                1                1

``` r
table(metadata_clean$Site)
```

    ## 
    ## saliva 
    ##     75

``` r
table(metadata_clean$matched_samples)
```

    ## 
    ## NO 
    ## 75

``` r
# Use cleaned dataset for subsequent analysis
metadata <- metadata_clean
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 75"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 75"

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

    ##    BioProject         Run      sample_name Age BMI Sex Region       Platform
    ## 1 PRJNA577839 SRR10294619  10_F_16_TG_Sn_3  NA  NA  NA  V3-V4 Illumina MiSeq
    ## 2 PRJNA577839 SRR10294620   9_F_15_TG_Sn_3  NA  NA  NA  V3-V4 Illumina MiSeq
    ## 3 PRJNA577839 SRR10294621 75_S_56_CG_nSn_0  NA  NA  NA  V3-V4 Illumina MiSeq
    ## 4 PRJNA577839 SRR10294622 74_S_55_CG_nSn_0  NA  NA  NA  V3-V4 Illumina MiSeq
    ## 5 PRJNA577839 SRR10294623 73_S_54_CG_nSn_0  NA  NA  NA  V3-V4 Illumina MiSeq
    ## 6 PRJNA577839 SRR10294624 72_S_53_CG_nSn_0  NA  NA  NA  V3-V4 Illumina MiSeq
    ##   Country   Site          public_time PCR_primers
    ## 1 Germany saliva 2020-02-03T00:00:00Z   341F-806R
    ## 2 Germany saliva 2020-02-03T00:00:00Z   341F-806R
    ## 3 Germany saliva 2020-02-03T00:00:00Z   341F-806R
    ## 4 Germany saliva 2020-02-03T00:00:00Z   341F-806R
    ## 5 Germany saliva 2020-02-03T00:00:00Z   341F-806R
    ## 6 Germany saliva 2020-02-03T00:00:00Z   341F-806R
    ##                                                                                                                                         title
    ## 1 Influence of the long-term use of oral hygiene products containing stannous ions on the salivary microbiome – a randomized controlled trial
    ## 2 Influence of the long-term use of oral hygiene products containing stannous ions on the salivary microbiome – a randomized controlled trial
    ## 3 Influence of the long-term use of oral hygiene products containing stannous ions on the salivary microbiome – a randomized controlled trial
    ## 4 Influence of the long-term use of oral hygiene products containing stannous ions on the salivary microbiome – a randomized controlled trial
    ## 5 Influence of the long-term use of oral hygiene products containing stannous ions on the salivary microbiome – a randomized controlled trial
    ## 6 Influence of the long-term use of oral hygiene products containing stannous ions on the salivary microbiome – a randomized controlled trial
    ##   Disease_study Status matched_samples
    ## 1            NO     TG              NO
    ## 2            NO     TG              NO
    ## 3            NO     CG              NO
    ## 4            NO     CG              NO
    ## 5            NO     CG              NO
    ## 6            NO     CG              NO

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
