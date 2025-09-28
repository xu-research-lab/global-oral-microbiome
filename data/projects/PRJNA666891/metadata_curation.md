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
setwd("D:/PRJNA666891")
initial_meta<-read.csv("SraRunTable .csv",header = T)
head(initial_meta)
```

    ##           Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample
    ## 1 SRR12759013        WGS        302 25284044 PRJNA666891 SAMN16330774
    ## 2 SRR12759014        WGS        302 26027568 PRJNA666891 SAMN16330773
    ## 3 SRR12759015        WGS        302 22138110 PRJNA666891 SAMN16330772
    ## 4 SRR12759016   AMPLICON        602    25886 PRJNA666891 SAMN16330834
    ## 5 SRR12759017   AMPLICON        602    21070 PRJNA666891 SAMN16330833
    ## 6 SRR12759018   AMPLICON        602   136654 PRJNA666891 SAMN16330832
    ##                                BioSampleModel    Bytes
    ## 1 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 12645000
    ## 2 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 12851698
    ## 3 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 10916313
    ## 4 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral    93813
    ## 5 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral    90917
    ## 6 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral   149040
    ##                 Center.Name Collection_Date Consent DATASTORE.filetype
    ## 1 J. CRAIG VENTER INSTITUTE      2013-08-09  public   fastq,run.zq,sra
    ## 2 J. CRAIG VENTER INSTITUTE      2013-08-09  public   fastq,run.zq,sra
    ## 3 J. CRAIG VENTER INSTITUTE      2013-08-09  public   fastq,run.zq,sra
    ## 4 J. CRAIG VENTER INSTITUTE      2013-10-30  public   fastq,run.zq,sra
    ## 5 J. CRAIG VENTER INSTITUTE      2013-10-30  public   run.zq,sra,fastq
    ## 6 J. CRAIG VENTER INSTITUTE      2013-10-30  public   fastq,run.zq,sra
    ##   DATASTORE.provider                     DATASTORE.region
    ## 1         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 2         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 3         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 4         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 5         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 6         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ##                    env_broad_scale       env_local_scale        env_medium
    ## 1 Human microbiome\\, WGS Illumina Human oral metagenome Human oral cavity
    ## 2 Human microbiome\\, WGS Illumina Human oral metagenome Human oral cavity
    ## 3 Human microbiome\\, WGS Illumina Human oral metagenome Human oral cavity
    ## 4 Human microbiome\\, 16S Illumina Human oral metagenome Human oral cavity
    ## 5 Human microbiome\\, 16S Illumina Human oral metagenome Human oral cavity
    ## 6 Human microbiome\\, 16S Illumina Human oral metagenome Human oral cavity
    ##   Experiment geo_loc_name_country geo_loc_name_country_continent
    ## 1 SRX9229871               Poland                         Europe
    ## 2 SRX9229870               Poland                         Europe
    ## 3 SRX9229869               Poland                         Europe
    ## 4 SRX9229868               Poland                         Europe
    ## 5 SRX9229867               Poland                         Europe
    ## 6 SRX9229866               Poland                         Europe
    ##                                     geo_loc_name         Host host_subject_id
    ## 1 Poland: Poznan\\, Greater Poland Cancer Center Homo sapiens           P07-Z
    ## 2 Poland: Poznan\\, Greater Poland Cancer Center Homo sapiens           P07-G
    ## 3 Poland: Poznan\\, Greater Poland Cancer Center Homo sapiens           P06-Z
    ## 4 Poland: Poznan\\, Greater Poland Cancer Center Homo sapiens        P0028_TT
    ## 5 Poland: Poznan\\, Greater Poland Cancer Center Homo sapiens       P0028_SAL
    ## 6 Poland: Poznan\\, Greater Poland Cancer Center Homo sapiens        P0028_HT
    ##       Instrument             lat_lon Library.Name LibraryLayout
    ## 1    NextSeq 500 52.4020 N 16.9359 E        P07-Z        PAIRED
    ## 2    NextSeq 500 52.4020 N 16.9359 E        P07-G        PAIRED
    ## 3    NextSeq 500 52.4020 N 16.9359 E        P06-Z        PAIRED
    ## 4 Illumina MiSeq 52.4020 N 16.9359 E     P0028_TT        PAIRED
    ## 5 Illumina MiSeq 52.4020 N 16.9359 E    P0028_SAL        PAIRED
    ## 6 Illumina MiSeq 52.4020 N 16.9359 E     P0028_HT        PAIRED
    ##   LibrarySelection LibrarySource   Organism Platform          ReleaseDate
    ## 1       RANDOM PCR   METAGENOMIC metagenome ILLUMINA 2020-10-01T00:00:00Z
    ## 2       RANDOM PCR   METAGENOMIC metagenome ILLUMINA 2020-10-01T00:00:00Z
    ## 3       RANDOM PCR   METAGENOMIC metagenome ILLUMINA 2020-10-01T00:00:00Z
    ## 4              PCR   METAGENOMIC metagenome ILLUMINA 2020-10-01T00:00:00Z
    ## 5              PCR   METAGENOMIC metagenome ILLUMINA 2020-10-01T00:00:00Z
    ## 6              PCR   METAGENOMIC metagenome ILLUMINA 2020-10-01T00:00:00Z
    ##            create_date version Sample.Name SRA.Study   host_tissue_sampled
    ## 1 2020-10-01T11:31:00Z       1       P07-Z SRP286018                throat
    ## 2 2020-10-01T11:32:00Z       1       P07-G SRP286018                throat
    ## 3 2020-10-01T11:32:00Z       1       P06-Z SRP286018 bottom of oral cavity
    ## 4 2020-10-01T11:32:00Z       1    P0028_TT SRP286018                      
    ## 5 2020-10-01T11:32:00Z       1   P0028_SAL SRP286018                      
    ## 6 2020-10-01T11:32:00Z       1    P0028_HT SRP286018

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
  sample_name =initial_meta$host_subject_id,
  PCR_primers = initial_meta$LibrarySelection,
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  Site = initial_meta$host_tissue_sampled,
  public_time = initial_meta$Collection_Date,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
 metadata$Age <- NA
  metadata$BMI <-NA
  metadata$Sex <-NA
  metadata$Region <- "V4"
  metadata$title <- "Oral Microbial Species and Virulence Factors Associated with Oral Squamous Cell Carcinoma"
  metadata$Disease_study <- "YES"
  
  # Assign disease status based on sample_name patterns
metadata<-metadata%>%
  mutate(Status=case_when(grepl("SAL",sample_name)~"OSCC_Saliva",grepl("TT|G|Z",sample_name)~"OSCC_Tumor",grepl("HT",sample_name)~"OSCC_Non-Tumor"))

#get subjects id according to sample id
metadata<-metadata%>%
  mutate(Site=case_when(grepl("SAL",sample_name)~"Saliva",grepl("TT|G|Z",sample_name)~"Tumor",grepl("HT",sample_name)~"OSCC_Non-Tumor"))


metadata$sample_name<-gsub("-Z|_SAL|_TT|_HT|-G","",metadata$sample_name)
#get subjects id with Saliva available
keep_subjects<-subset(metadata,Site %in% "Saliva")
#remove subjects without Saliva samples
metadata<-subset(metadata,sample_name %in% keep_subjects$sample_name)
#get subjects sampling counts
sample_counts <- table(metadata$sample_name)
# Create matched_samples column, any subjects have >=2 samples will be label "YES", otherwise "NO"
metadata$matched_samples <- ifelse(metadata$sample_name %in% names(sample_counts[sample_counts >= 2]), "YES", "NO")
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 54"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 18"

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

    ##    BioProject         Run sample_name Age BMI Sex Region       Platform Country
    ## 4 PRJNA666891 SRR12759016       P0028  NA  NA  NA     V4 Illumina MiSeq  Poland
    ## 5 PRJNA666891 SRR12759017       P0028  NA  NA  NA     V4 Illumina MiSeq  Poland
    ## 6 PRJNA666891 SRR12759018       P0028  NA  NA  NA     V4 Illumina MiSeq  Poland
    ## 7 PRJNA666891 SRR12759019       P0027  NA  NA  NA     V4 Illumina MiSeq  Poland
    ## 8 PRJNA666891 SRR12759020       P0027  NA  NA  NA     V4 Illumina MiSeq  Poland
    ## 9 PRJNA666891 SRR12759021       P0027  NA  NA  NA     V4 Illumina MiSeq  Poland
    ##             Site public_time PCR_primers
    ## 4          Tumor  2013-10-30         PCR
    ## 5         Saliva  2013-10-30         PCR
    ## 6 OSCC_Non-Tumor  2013-10-30         PCR
    ## 7          Tumor  2013-10-30         PCR
    ## 8         Saliva  2013-10-30         PCR
    ## 9 OSCC_Non-Tumor  2013-10-30         PCR
    ##                                                                                       title
    ## 4 Oral Microbial Species and Virulence Factors Associated with Oral Squamous Cell Carcinoma
    ## 5 Oral Microbial Species and Virulence Factors Associated with Oral Squamous Cell Carcinoma
    ## 6 Oral Microbial Species and Virulence Factors Associated with Oral Squamous Cell Carcinoma
    ## 7 Oral Microbial Species and Virulence Factors Associated with Oral Squamous Cell Carcinoma
    ## 8 Oral Microbial Species and Virulence Factors Associated with Oral Squamous Cell Carcinoma
    ## 9 Oral Microbial Species and Virulence Factors Associated with Oral Squamous Cell Carcinoma
    ##   Disease_study         Status matched_samples
    ## 4           YES     OSCC_Tumor             YES
    ## 5           YES    OSCC_Saliva             YES
    ## 6           YES OSCC_Non-Tumor             YES
    ## 7           YES     OSCC_Tumor             YES
    ## 8           YES    OSCC_Saliva             YES
    ## 9           YES OSCC_Non-Tumor             YES
