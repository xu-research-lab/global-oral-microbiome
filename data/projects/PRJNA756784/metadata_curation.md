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
setwd("D:/PRJNA756784")
initial_meta<-read.csv("SraRunTable .csv",header = T)
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
    ##   hypervariable_region Sample_ID        sample_site sample_type sequence_length
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
metadata$BMI <- NA_character_
metadata$Sex <- NA_character_

metadata$Site<- "saliva"
metadata$PCR_primers <- "B969F-BA1406R"
metadata$Region<- "V3"
metadata$title <- "Comparative and analytical characterization of the oral bacteriome of smokeless tobacco users with oral squamous cell carcinoma"
metadata$Disease_study <- "YES"
metadata<-metadata%>%
  mutate(Status=case_when(grepl("CP",sample_name)~"OSCC",grepl("HTC",sample_name)~"HC",grepl("HI",sample_name)~"HC"))
metadata$matched_samples<- "NO"
```

# Step 4: Summarize and Validate the Metadata

``` r
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

    ##    BioProject         Run sample_name  Age  BMI  Sex Region            Platform
    ## 1 PRJNA756784 SRR15560887        HI11 <NA> <NA> <NA>     V3 Illumina HiSeq 2500
    ## 2 PRJNA756784 SRR15560888        HI10 <NA> <NA> <NA>     V3 Illumina HiSeq 2500
    ## 3 PRJNA756784 SRR15560889         HI9 <NA> <NA> <NA>     V3 Illumina HiSeq 2500
    ## 4 PRJNA756784 SRR15560890         HI8 <NA> <NA> <NA>     V3 Illumina HiSeq 2500
    ## 5 PRJNA756784 SRR15560891        CP21 <NA> <NA> <NA>     V3 Illumina HiSeq 2500
    ## 6 PRJNA756784 SRR15560892        CP19 <NA> <NA> <NA>     V3 Illumina HiSeq 2500
    ##   Country   Site          public_time   PCR_primers
    ## 1   India saliva 2021-08-22T00:00:00Z B969F-BA1406R
    ## 2   India saliva 2021-08-22T00:00:00Z B969F-BA1406R
    ## 3   India saliva 2021-08-22T00:00:00Z B969F-BA1406R
    ## 4   India saliva 2021-08-22T00:00:00Z B969F-BA1406R
    ## 5   India saliva 2021-08-22T00:00:00Z B969F-BA1406R
    ## 6   India saliva 2021-08-22T00:00:00Z B969F-BA1406R
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
