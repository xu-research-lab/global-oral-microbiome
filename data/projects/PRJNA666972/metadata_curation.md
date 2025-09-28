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
setwd("D:/PRJNA666972")
initial_meta<-read.csv("SraRunTable.csv",header = T)
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
metadata <- data.frame(
  BioProject = initial_meta$BioProject,
  Run = initial_meta$Run,
  sample_name = initial_meta$Library.Name,
  Platform = initial_meta$Platform,
  Country = initial_meta$geo_loc_name,
 
  public_time = initial_meta$Collection_Date,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
metadata$Age <-NA
  metadata$BMI <-NA
  metadata$Sex <-"M"
  metadata$PCR_primers <-"341F-806R"
  metadata$Region <- "V3-V4"
  metadata$Site <- "saliva"
  metadata$title <- "Study on the Salivary Microbial Alteration of Men With Head and Neck Cancer and Its Relationship With Symptoms in Southwest China"
  metadata$Disease_study <- "YES"
  metadata$Status <- "HNC"
  metadata$matched_samples <- "NO"
  
  
#extract subject id by remove prefix of sample_name
metadata$sample_name <- gsub("L", "", metadata$sample_name)

# Check initial data distribution
cat("Initial sample distribution:\n")
```

    ## Initial sample distribution:

``` r
table(metadata$sample_name)
```

    ## 
    ##   1 101 103 105 107 109  11 111 113 115 117 119 121 123 125 127 129  13 131 133 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 135 137 139 141 143 145 147 149  15 151 153 155 157 159 161 163 165 167 169  17 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 171 173 175 177 179 181 183 185 187 189  19 191 193 195 197 199 201 203 205 207 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 209  21 211 213 215 217 219 221 223 225 227 229  23 231 233 235 237 239  25  27 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ##  29   3  31  33  35  37  39  41  43  45  47  49   5  51  53  55  57  59  61  63 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ##  65  67  69   7  71  73  75  77  79  81  83  85  87  89   9  91  93  95  97  99 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1

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
    ##   1 101 103 105 107 109  11 111 113 115 117 119 121 123 125 127 129  13 131 133 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 135 137 139 141 143 145 147 149  15 151 153 155 157 159 161 163 165 167 169  17 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 171 173 175 177 179 181 183 185 187 189  19 191 193 195 197 199 201 203 205 207 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ## 209  21 211 213 215 217 219 221 223 225 227 229  23 231 233 235 237 239  25  27 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ##  29   3  31  33  35  37  39  41  43  45  47  49   5  51  53  55  57  59  61  63 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
    ##  65  67  69   7  71  73  75  77  79  81  83  85  87  89   9  91  93  95  97  99 
    ##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1

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

    ##    BioProject         Run sample_name Age BMI Sex Region Platform
    ## 1 PRJNA666972 SRR12771660         161  NA  NA   M  V3-V4 ILLUMINA
    ## 2 PRJNA666972 SRR12771661         167  NA  NA   M  V3-V4 ILLUMINA
    ## 3 PRJNA666972 SRR12771662         169  NA  NA   M  V3-V4 ILLUMINA
    ## 4 PRJNA666972 SRR12771663         171  NA  NA   M  V3-V4 ILLUMINA
    ## 5 PRJNA666972 SRR12771664         177  NA  NA   M  V3-V4 ILLUMINA
    ## 6 PRJNA666972 SRR12771665         179  NA  NA   M  V3-V4 ILLUMINA
    ##          Country   Site public_time PCR_primers
    ## 1 China: Chengdu saliva  2020-05-01   341F-806R
    ## 2 China: Chengdu saliva  2020-05-01   341F-806R
    ## 3 China: Chengdu saliva  2020-05-01   341F-806R
    ## 4 China: Chengdu saliva  2020-05-01   341F-806R
    ## 5 China: Chengdu saliva  2020-05-01   341F-806R
    ## 6 China: Chengdu saliva  2020-05-01   341F-806R
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
