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
```

    ## Warning: package 'dplyr' was built under R version 4.5.1

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(ggplot2)
```

    ## Warning: package 'ggplot2' was built under R version 4.5.1

``` r
library(cowplot)
```

    ## Warning: package 'cowplot' was built under R version 4.5.1

``` r
setwd("D:/PRJNA670573")
initial_meta<-read.csv("SraRunTable .csv",header = T)
head(initial_meta)
```

    ##           Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample
    ## 1 SRR12875167   AMPLICON        501 10891807 PRJNA670573 SAMN16513238
    ## 2 SRR12875168   AMPLICON        600  2483805 PRJNA670573 SAMN16513237
    ## 3 SRR12875169   AMPLICON        601   772765 PRJNA670573 SAMN16513236
    ## 4 SRR12875170   AMPLICON        497  1592126 PRJNA670573 SAMN16513235
    ## 5 SRR12875171   AMPLICON        601 21880635 PRJNA670573 SAMN16513234
    ## 6 SRR12875172   AMPLICON        500   697275 PRJNA670573 SAMN16513233
    ##                                BioSampleModel    Bytes
    ## 1 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral  6355154
    ## 2 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral  1909762
    ## 3 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral   646343
    ## 4 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral   973147
    ## 5 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 15844641
    ## 6 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral   485947
    ##                                      Center.Name Collection_Date Consent
    ## 1 INSTITUTE OF MICROBIOLOGY OF THE CAS, V. V. I.         missing  public
    ## 2 INSTITUTE OF MICROBIOLOGY OF THE CAS, V. V. I.         missing  public
    ## 3 INSTITUTE OF MICROBIOLOGY OF THE CAS, V. V. I.         missing  public
    ## 4 INSTITUTE OF MICROBIOLOGY OF THE CAS, V. V. I.         missing  public
    ## 5 INSTITUTE OF MICROBIOLOGY OF THE CAS, V. V. I.         missing  public
    ## 6 INSTITUTE OF MICROBIOLOGY OF THE CAS, V. V. I.         missing  public
    ##   Cultivar DATASTORE.filetype DATASTORE.provider
    ## 1  missing   fastq,run.zq,sra         gs,ncbi,s3
    ## 2  missing   fastq,run.zq,sra         gs,ncbi,s3
    ## 3  missing   fastq,run.zq,sra         gs,ncbi,s3
    ## 4  missing   fastq,run.zq,sra         gs,ncbi,s3
    ## 5  missing   fastq,run.zq,sra         gs,ncbi,s3
    ## 6  missing   fastq,run.zq,sra         gs,ncbi,s3
    ##                       DATASTORE.region Ecotype env_broad_scale
    ## 1 gs.us-east1,ncbi.public,s3.us-east-1 missing oral microbiome
    ## 2 gs.us-east1,ncbi.public,s3.us-east-1 missing oral microbiome
    ## 3 gs.us-east1,ncbi.public,s3.us-east-1 missing oral microbiome
    ## 4 gs.us-east1,ncbi.public,s3.us-east-1 missing oral microbiome
    ## 5 gs.us-east1,ncbi.public,s3.us-east-1 missing oral microbiome
    ## 6 gs.us-east1,ncbi.public,s3.us-east-1 missing oral microbiome
    ##      env_local_scale         env_medium Experiment geo_loc_name_country
    ## 1 periodontal pocket subgingival plaque SRX9341343       Czech Republic
    ## 2 periodontal pocket subgingival plaque SRX9341342       Czech Republic
    ## 3 periodontal pocket subgingival plaque SRX9341341       Czech Republic
    ## 4 periodontal pocket subgingival plaque SRX9341340       Czech Republic
    ## 5 periodontal pocket subgingival plaque SRX9341339       Czech Republic
    ## 6 periodontal pocket subgingival plaque SRX9341338       Czech Republic
    ##   geo_loc_name_country_continent           geo_loc_name         HOST
    ## 1                         Europe Czech Republic: Prague Homo sapiens
    ## 2                         Europe Czech Republic: Prague Homo sapiens
    ## 3                         Europe Czech Republic: Prague Homo sapiens
    ## 4                         Europe Czech Republic: Prague Homo sapiens
    ## 5                         Europe Czech Republic: Prague Homo sapiens
    ## 6                         Europe Czech Republic: Prague Homo sapiens
    ##       Instrument isol_growth_condt isolate                   lat_lon
    ## 1 Illumina MiSeq           missing missing 50.0746114 N 14.4265931 E
    ## 2 Illumina MiSeq           missing missing 50.0746114 N 14.4265931 E
    ## 3 Illumina MiSeq           missing missing 50.0746114 N 14.4265931 E
    ## 4 Illumina MiSeq           missing missing 50.0746114 N 14.4265931 E
    ## 5 Illumina MiSeq           missing missing 50.0746114 N 14.4265931 E
    ## 6 Illumina MiSeq           missing missing 50.0746114 N 14.4265931 E
    ##   Library.Name LibraryLayout LibrarySelection LibrarySource
    ## 1   library 26        PAIRED              PCR   METAGENOMIC
    ## 2   library 25        PAIRED              PCR   METAGENOMIC
    ## 3   library 24        PAIRED              PCR   METAGENOMIC
    ## 4   library 23        PAIRED              PCR   METAGENOMIC
    ## 5   library 22        PAIRED              PCR   METAGENOMIC
    ## 6   library 21        PAIRED              PCR   METAGENOMIC
    ##                Organism Platform          ReleaseDate          create_date
    ## 1 human oral metagenome ILLUMINA 2020-11-16T00:00:00Z 2020-10-22T02:25:00Z
    ## 2 human oral metagenome ILLUMINA 2020-11-16T00:00:00Z 2020-10-22T02:25:00Z
    ## 3 human oral metagenome ILLUMINA 2020-11-16T00:00:00Z 2020-10-22T02:25:00Z
    ## 4 human oral metagenome ILLUMINA 2020-11-16T00:00:00Z 2020-10-22T02:25:00Z
    ## 5 human oral metagenome ILLUMINA 2020-11-16T00:00:00Z 2020-10-22T02:26:00Z
    ## 6 human oral metagenome ILLUMINA 2020-11-16T00:00:00Z 2020-10-22T02:25:00Z
    ##   version Sample.Name source_material_id SRA.Study  strain
    ## 1       1       AP_26   sequence file 26 SRP288098 missing
    ## 2       1       AP_25   sequence file 25 SRP288098 missing
    ## 3       1       AP_24   sequence file 24 SRP288098 missing
    ## 4       1       AP_23   sequence file 23 SRP288098 missing
    ## 5       1       AP_22   sequence file 22 SRP288098 missing
    ## 6       1       AP_21   sequence file 21 SRP288098 missing

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
metadata <- data.frame(
  BioProject = initial_meta$BioProject,
  Run = initial_meta$Run,
  sample_name = initial_meta$Sample.Name,
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  public_time = initial_meta$ReleaseDate,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
  metadata$Age <- NA
  metadata$BMI <- NA
  metadata$Sex <- NA
  metadata$PCR_primers <-"530F-907R"
  metadata$Region <-"V4-V5"
  metadata$Site <-"oral"
  metadata$title <-"The Oral Microbiome in Periodontal Health"
  metadata$Disease_study <- "YES"
  metadata<-metadata%>%
    mutate(Status=case_when(grepl("AP",sample_name)~"AP",grepl("HY|HO",sample_name)~"HC"))
  metadata$matched_samples <-"NO"
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check initial data distribution
cat("Initial sample distribution:\n")
```

    ## Initial sample distribution:

``` r
table(metadata$sample_name)
```

    ## 
    ##  AP_1 AP_10 AP_11 AP_12 AP_13 AP_14 AP_15 AP_16 AP_17 AP_18 AP_19  AP_2 AP_20 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## AP_21 AP_22 AP_23 AP_24 AP_25 AP_26 AP_27 AP_28 AP_29  AP_3 AP_30 AP_31 AP_32 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## AP_33 AP_34 AP_35 AP_36 AP_37 AP_38 AP_39  AP_4 AP_40 AP_41 AP_42 AP_43 AP_44 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## AP_45  AP_5  AP_6  AP_7  AP_8  AP_9  HO_1 HO_10 HO_11 HO_12 HO_13 HO_14 HO_15 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## HO_16 HO_17  HO_2  HO_3  HO_4  HO_5  HO_6  HO_7  HO_8  HO_9  HY_1 HY_10 HY_11 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## HY_12 HY_13 HY_14 HY_15 HY_16 HY_17 HY_18 HY_19  HY_2 HY_20 HY_21 HY_22 HY_23 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## HY_24 HY_25 HY_26 HY_27 HY_28 HY_29  HY_3 HY_30 HY_31 HY_32 HY_33 HY_34 HY_35 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## HY_36 HY_37 HY_38 HY_39  HY_4 HY_40 HY_41 HY_42 HY_43 HY_44 HY_45 HY_46 HY_47 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## HY_48 HY_49  HY_5 HY_50 HY_51 HY_52 HY_53 HY_54 HY_55 HY_56 HY_57 HY_58 HY_59 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##  HY_6 HY_60 HY_61 HY_62 HY_63 HY_64 HY_65 HY_66 HY_67 HY_68 HY_69  HY_7 HY_70 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## HY_71 HY_72 HY_73 HY_74 HY_75 HY_76 HY_77 HY_78 HY_79  HY_8 HY_80 HY_81 HY_82 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## HY_83 HY_84 HY_85 HY_86 HY_87 HY_88 HY_89  HY_9 HY_90 HY_91 
    ##     1     1     1     1     1     1     1     1     1     1

``` r
table(metadata$Site)
```

    ## 
    ## oral 
    ##  153

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
    ##  AP_1 AP_10 AP_11 AP_12 AP_13 AP_14 AP_15 AP_16 AP_17 AP_18 AP_19  AP_2 AP_20 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## AP_21 AP_22 AP_23 AP_24 AP_25 AP_26 AP_27 AP_28 AP_29  AP_3 AP_30 AP_31 AP_32 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## AP_33 AP_34 AP_35 AP_36 AP_37 AP_38 AP_39  AP_4 AP_40 AP_41 AP_42 AP_43 AP_44 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## AP_45  AP_5  AP_6  AP_7  AP_8  AP_9  HO_1 HO_10 HO_11 HO_12 HO_13 HO_14 HO_15 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## HO_16 HO_17  HO_2  HO_3  HO_4  HO_5  HO_6  HO_7  HO_8  HO_9  HY_1 HY_10 HY_11 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## HY_12 HY_13 HY_14 HY_15 HY_16 HY_17 HY_18 HY_19  HY_2 HY_20 HY_21 HY_22 HY_23 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## HY_24 HY_25 HY_26 HY_27 HY_28 HY_29  HY_3 HY_30 HY_31 HY_32 HY_33 HY_34 HY_35 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## HY_36 HY_37 HY_38 HY_39  HY_4 HY_40 HY_41 HY_42 HY_43 HY_44 HY_45 HY_46 HY_47 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## HY_48 HY_49  HY_5 HY_50 HY_51 HY_52 HY_53 HY_54 HY_55 HY_56 HY_57 HY_58 HY_59 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ##  HY_6 HY_60 HY_61 HY_62 HY_63 HY_64 HY_65 HY_66 HY_67 HY_68 HY_69  HY_7 HY_70 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## HY_71 HY_72 HY_73 HY_74 HY_75 HY_76 HY_77 HY_78 HY_79  HY_8 HY_80 HY_81 HY_82 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## HY_83 HY_84 HY_85 HY_86 HY_87 HY_88 HY_89  HY_9 HY_90 HY_91 
    ##     1     1     1     1     1     1     1     1     1     1

# Step 5: Visualize the Metadata

``` r
knitr::include_graphics("Metadata_Summary.png")
```

![](Metadata_Summary.png)<!-- -->

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

    ##    BioProject         Run sample_name Age BMI Sex Region       Platform
    ## 1 PRJNA670573 SRR12875167       AP_26  NA  NA  NA  V4-V5 Illumina MiSeq
    ## 2 PRJNA670573 SRR12875168       AP_25  NA  NA  NA  V4-V5 Illumina MiSeq
    ## 3 PRJNA670573 SRR12875169       AP_24  NA  NA  NA  V4-V5 Illumina MiSeq
    ## 4 PRJNA670573 SRR12875170       AP_23  NA  NA  NA  V4-V5 Illumina MiSeq
    ## 5 PRJNA670573 SRR12875171       AP_22  NA  NA  NA  V4-V5 Illumina MiSeq
    ## 6 PRJNA670573 SRR12875172       AP_21  NA  NA  NA  V4-V5 Illumina MiSeq
    ##          Country Site          public_time PCR_primers
    ## 1 Czech Republic oral 2020-11-16T00:00:00Z   530F-907R
    ## 2 Czech Republic oral 2020-11-16T00:00:00Z   530F-907R
    ## 3 Czech Republic oral 2020-11-16T00:00:00Z   530F-907R
    ## 4 Czech Republic oral 2020-11-16T00:00:00Z   530F-907R
    ## 5 Czech Republic oral 2020-11-16T00:00:00Z   530F-907R
    ## 6 Czech Republic oral 2020-11-16T00:00:00Z   530F-907R
    ##                                       title Disease_study Status
    ## 1 The Oral Microbiome in Periodontal Health           YES     AP
    ## 2 The Oral Microbiome in Periodontal Health           YES     AP
    ## 3 The Oral Microbiome in Periodontal Health           YES     AP
    ## 4 The Oral Microbiome in Periodontal Health           YES     AP
    ## 5 The Oral Microbiome in Periodontal Health           YES     AP
    ## 6 The Oral Microbiome in Periodontal Health           YES     AP
    ##   matched_samples
    ## 1              NO
    ## 2              NO
    ## 3              NO
    ## 4              NO
    ## 5              NO
    ## 6              NO
