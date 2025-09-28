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
setwd("D:/PRJNA751046")
initial_meta<-read.csv("SraRunTable .csv",header = T)
head(initial_meta)
```

    ##           Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample
    ## 1 SRR21474442   AMPLICON        602 17958262 PRJNA751046 SAMN30722499
    ## 2 SRR21474443   AMPLICON        602 30789290 PRJNA751046 SAMN30722498
    ## 3 SRR21474450   AMPLICON        602 36920058 PRJNA751046 SAMN30722509
    ## 4 SRR21474451   AMPLICON        602 28743694 PRJNA751046 SAMN30722508
    ## 5 SRR21484724   AMPLICON        602 56516362 PRJNA751046 SAMN30725542
    ## 6 SRR21484725   AMPLICON        602 26359774 PRJNA751046 SAMN30725541
    ##                                BioSampleModel    Bytes
    ## 1 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 10385393
    ## 2 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 17344813
    ## 3 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 20145434
    ## 4 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 16293854
    ## 5 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 32135564
    ## 6 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 15129962
    ##                                             Center.Name Collection_Date Consent
    ## 1 NARSEE MONJEE INSTITUTE OF MANAGEMENT STUDIES (NMIMS)      2018-02-10  public
    ## 2 NARSEE MONJEE INSTITUTE OF MANAGEMENT STUDIES (NMIMS)      2018-02-05  public
    ## 3 NARSEE MONJEE INSTITUTE OF MANAGEMENT STUDIES (NMIMS)      2018-04-10  public
    ## 4 NARSEE MONJEE INSTITUTE OF MANAGEMENT STUDIES (NMIMS)      2018-03-27  public
    ## 5 NARSEE MONJEE INSTITUTE OF MANAGEMENT STUDIES (NMIMS)      2018-08-15  public
    ## 6 NARSEE MONJEE INSTITUTE OF MANAGEMENT STUDIES (NMIMS)      2018-08-08  public
    ##   DATASTORE.filetype DATASTORE.provider                     DATASTORE.region
    ## 1   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 2   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 3   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 4   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 5   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 6   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ##         env_broad_scale env_local_scale env_medium  Experiment
    ## 1           oral cavity  not applicable       <NA> SRX17478091
    ## 2           oral cavity  not applicable       <NA> SRX17478090
    ## 3           oral cavity  not applicable       <NA> SRX17478083
    ## 4           oral cavity  not applicable       <NA> SRX17478082
    ## 5 Human Oral microbiome            Oral    missing SRX17488295
    ## 6 Human Oral microbiome            Oral    missing SRX17488294
    ##   geo_loc_name_country geo_loc_name_country_continent geo_loc_name         HOST
    ## 1                India                           Asia        India Homo sapiens
    ## 2                India                           Asia        India Homo sapiens
    ## 3                India                           Asia        India Homo sapiens
    ## 4                India                           Asia        India Homo sapiens
    ## 5                India                           Asia        India Homo sapiens
    ## 6                India                           Asia        India Homo sapiens
    ##       Instrument lat_lon Library.Name LibraryLayout LibrarySelection
    ## 1 Illumina MiSeq    <NA>         S207        PAIRED              PCR
    ## 2 Illumina MiSeq    <NA>         S195        PAIRED              PCR
    ## 3 Illumina MiSeq    <NA>         S244        PAIRED              PCR
    ## 4 Illumina MiSeq    <NA>         S232        PAIRED              PCR
    ## 5 Illumina MiSeq missing         S214        PAIRED              PCR
    ## 6 Illumina MiSeq missing         S202        PAIRED              PCR
    ##   LibrarySource              Organism Platform          ReleaseDate
    ## 1   METAGENOMIC human oral metagenome ILLUMINA 2022-09-08T00:00:00Z
    ## 2   METAGENOMIC human oral metagenome ILLUMINA 2022-09-08T00:00:00Z
    ## 3   METAGENOMIC human oral metagenome ILLUMINA 2022-09-08T00:00:00Z
    ## 4   METAGENOMIC human oral metagenome ILLUMINA 2022-09-08T00:00:00Z
    ## 5   METAGENOMIC human oral metagenome ILLUMINA 2022-09-09T00:00:00Z
    ## 6   METAGENOMIC human oral metagenome ILLUMINA 2022-09-09T00:00:00Z
    ##            create_date version Sample.Name SRA.Study Host_Age
    ## 1 2022-09-08T17:38:00Z       1         C21 SRP330565       NA
    ## 2 2022-09-08T18:18:00Z       1         C20 SRP330565       NA
    ## 3 2022-09-08T17:22:00Z       1         C36 SRP330565       NA
    ## 4 2022-09-08T17:09:00Z       1         C34 SRP330565       NA
    ## 5 2022-09-09T02:43:00Z       1          O7 SRP330565       NA
    ## 6 2022-09-09T02:15:00Z       1          O1 SRP330565       NA

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
  Platform=initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  public_time = initial_meta$ReleaseDate,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
metadata$Age <- NA_character_
metadata$Site <- "oral"
metadata$BMI <- NA_character_
metadata$Sex <- NA_character_
metadata$PCR_primers <- "B969F-BA1406"
metadata$Region<- "V6-V8"
metadata$title <-"Identification & correlation of bacterial diversity in oral cancer and long- term tobacco chewers- A case- control pilot study."
metadata$Disease_study <- "YES"
metadata<-metadata%>%
  mutate(Status=case_when(grepl("O",sample_name)~"OC",grepl("T",sample_name)~"Tobacco",grepl("C",sample_name)~"HC"))
metadata$matched_samples<- "NO"
```

``` r
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 112"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 112"

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

    ##    BioProject         Run sample_name  Age  BMI  Sex Region       Platform
    ## 1 PRJNA751046 SRR21474442         C21 <NA> <NA> <NA>  V6-V8 Illumina MiSeq
    ## 2 PRJNA751046 SRR21474443         C20 <NA> <NA> <NA>  V6-V8 Illumina MiSeq
    ## 3 PRJNA751046 SRR21474450         C36 <NA> <NA> <NA>  V6-V8 Illumina MiSeq
    ## 4 PRJNA751046 SRR21474451         C34 <NA> <NA> <NA>  V6-V8 Illumina MiSeq
    ## 5 PRJNA751046 SRR21484724          O7 <NA> <NA> <NA>  V6-V8 Illumina MiSeq
    ## 6 PRJNA751046 SRR21484725          O1 <NA> <NA> <NA>  V6-V8 Illumina MiSeq
    ##   Country Site          public_time  PCR_primers
    ## 1   India oral 2022-09-08T00:00:00Z B969F-BA1406
    ## 2   India oral 2022-09-08T00:00:00Z B969F-BA1406
    ## 3   India oral 2022-09-08T00:00:00Z B969F-BA1406
    ## 4   India oral 2022-09-08T00:00:00Z B969F-BA1406
    ## 5   India oral 2022-09-09T00:00:00Z B969F-BA1406
    ## 6   India oral 2022-09-09T00:00:00Z B969F-BA1406
    ##                                                                                                                             title
    ## 1 Identification & correlation of bacterial diversity in oral cancer and long- term tobacco chewers- A case- control pilot study.
    ## 2 Identification & correlation of bacterial diversity in oral cancer and long- term tobacco chewers- A case- control pilot study.
    ## 3 Identification & correlation of bacterial diversity in oral cancer and long- term tobacco chewers- A case- control pilot study.
    ## 4 Identification & correlation of bacterial diversity in oral cancer and long- term tobacco chewers- A case- control pilot study.
    ## 5 Identification & correlation of bacterial diversity in oral cancer and long- term tobacco chewers- A case- control pilot study.
    ## 6 Identification & correlation of bacterial diversity in oral cancer and long- term tobacco chewers- A case- control pilot study.
    ##   Disease_study Status matched_samples
    ## 1           YES     HC              NO
    ## 2           YES     HC              NO
    ## 3           YES     HC              NO
    ## 4           YES     HC              NO
    ## 5           YES     OC              NO
    ## 6           YES     OC              NO
