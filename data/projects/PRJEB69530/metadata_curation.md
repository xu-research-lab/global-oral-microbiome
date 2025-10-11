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
setwd("D:/oral_project/PRJEB69530")

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

    ##           Run anonymized_name Assay.Type AvgSpotLen barcode..exp.    Bases
    ## 1 ERR12266234             M03   AMPLICON         99  TGCGCTGAATGT  6077187
    ## 2 ERR12266235             M03   AMPLICON         99  TGCCGAGTAATC  8736062
    ## 3 ERR12266238             M03   AMPLICON         99  TGGAGTAGGTGG  6389316
    ## 4 ERR12266240             M03   AMPLICON         99  TAGTAGCACCTG  5202013
    ## 5 ERR12266242             M03   AMPLICON         99  ACCGGTATGTAC  3681628
    ## 6 ERR12266243             M03   AMPLICON         99  TGTGGTGATGTA 16477955
    ##   BioProject      BioSample    Bytes center_name..exp.
    ## 1 PRJEB69530 SAMEA114622785  4386365              CCME
    ## 2 PRJEB69530 SAMEA114622786  6622860              CCME
    ## 3 PRJEB69530 SAMEA114622789  4506689              CCME
    ## 4 PRJEB69530 SAMEA114622791  3822628              CCME
    ## 5 PRJEB69530 SAMEA114622793  2605179              CCME
    ## 6 PRJEB69530 SAMEA114622794 12351935              CCME
    ##                                          Center.Name center_project_name..exp.
    ## 1 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT         Bangladesh_Travel
    ## 2 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT         Bangladesh_Travel
    ## 3 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT         Bangladesh_Travel
    ## 4 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT         Bangladesh_Travel
    ## 5 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT         Bangladesh_Travel
    ## 6 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT         Bangladesh_Travel
    ##   Collection_Date collection_timestamp Consent COUNTRY DATASTORE.filetype
    ## 1      2012-01-13           2012-04-09  public     USA   fastq,run.zq,sra
    ## 2      2012-01-13           2012-04-01  public     USA   fastq,run.zq,sra
    ## 3      2012-01-14           2012-04-11  public     USA   fastq,run.zq,sra
    ## 4      2012-01-14           2012-02-23  public     USA   fastq,run.zq,sra
    ## 5      2012-01-16           2012-02-16  public     USA   fastq,run.zq,sra
    ## 6      2012-01-16           2012-02-17  public     USA   fastq,run.zq,sra
    ##   DATASTORE.provider                     DATASTORE.region
    ## 1         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 2         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 3         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 4         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 5         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 6         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ##   ENA.FIRST.PUBLIC..run.     ena_first_public ENA.LAST.UPDATE..run.
    ## 1             2023-11-16 2023-11-16T01:29:03Z            2023-11-16
    ## 2             2023-11-16 2023-11-16T01:29:03Z            2023-11-16
    ## 3             2023-11-16 2023-11-16T01:29:03Z            2023-11-16
    ## 4             2023-11-16 2023-11-16T01:29:03Z            2023-11-16
    ## 5             2023-11-16 2023-11-16T01:29:03Z            2023-11-16
    ## 6             2023-11-16 2023-11-16T01:29:03Z            2023-11-16
    ##        ena_last_update epoch_time experiment_center..exp.  Experiment
    ## 1 2023-11-16T01:29:03Z 1200096000                    CCME ERX11676478
    ## 2 2023-11-16T01:29:03Z 1200096000                    CCME ERX11676479
    ## 3 2023-11-16T01:29:03Z 1200182400                    CCME ERX11676482
    ## 4 2023-11-16T01:29:03Z 1200182400                    CCME ERX11676484
    ## 5 2023-11-16T01:29:03Z 1200355200                    CCME ERX11676486
    ## 6 2023-11-16T01:29:03Z 1200355200                    CCME ERX11676487
    ##   geo_loc_name_country geo_loc_name_country_continent
    ## 1                  USA                  North America
    ## 2                  USA                  North America
    ## 3                  USA                  North America
    ## 4                  USA                  North America
    ## 5                  USA                  North America
    ## 6                  USA                  North America
    ##   geographic_location_.country_and.or_sea. host_subject_id          Instrument
    ## 1                                      USA        1015:M03 Illumina HiSeq 2000
    ## 2                                      USA        1015:M03 Illumina HiSeq 2000
    ## 3                                      USA        1015:M03 Illumina HiSeq 2000
    ## 4                                      USA        1015:M03 Illumina HiSeq 2000
    ## 5                                      USA        1015:M03 Illumina HiSeq 2000
    ## 6                                      USA        1015:M03 Illumina HiSeq 2000
    ##       Library.Name LibraryLayout LibrarySelection LibrarySource linker..exp.
    ## 1  1015.1.13.M03.M        SINGLE              PCR   METAGENOMIC           GT
    ## 2 1015.1.13.M03.RH        SINGLE              PCR   METAGENOMIC           GT
    ## 3  1015.1.14.M03.M        SINGLE              PCR   METAGENOMIC           GT
    ## 4  1015.1.14.M03.S        SINGLE              PCR   METAGENOMIC           GT
    ## 5  1015.1.16.M03.M        SINGLE              PCR   METAGENOMIC           GT
    ## 6 1015.1.16.M03.RH        SINGLE              PCR   METAGENOMIC           GT
    ##                Organism                                 pcr_primers..exp.
    ## 1 human oral metagenome FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 2 human skin metagenome FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 3 human oral metagenome FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 4  human gut metagenome FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 5 human oral metagenome FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 6 human skin metagenome FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ##   Platform          primer..exp.          ReleaseDate run_center..exp.
    ## 1 ILLUMINA GTGTGCCAGCMGCCGCGGTAA 2024-07-21T00:00:00Z             CCME
    ## 2 ILLUMINA GTGTGCCAGCMGCCGCGGTAA 2024-07-21T00:00:00Z             CCME
    ## 3 ILLUMINA GTGTGCCAGCMGCCGCGGTAA 2024-07-22T00:00:00Z             CCME
    ## 4 ILLUMINA GTGTGCCAGCMGCCGCGGTAA 2024-07-21T00:00:00Z             CCME
    ## 5 ILLUMINA GTGTGCCAGCMGCCGCGGTAA 2024-07-22T00:00:00Z             CCME
    ## 6 ILLUMINA GTGTGCCAGCMGCCGCGGTAA 2024-07-22T00:00:00Z             CCME
    ##   run_date..exp.          create_date version
    ## 1        2/25/13 2024-07-30T21:57:00Z       1
    ## 2        2/25/13 2024-07-30T21:53:00Z       1
    ## 3        2/25/13 2024-07-31T01:53:00Z       1
    ## 4        2/25/13 2024-07-30T21:59:00Z       1
    ## 5        2/25/13 2024-07-30T22:07:00Z       1
    ## 6        2/25/13 2024-07-30T22:05:00Z       1
    ##                     run_prefix..exp. samp_collect_device..exp.
    ## 1     Bangladesh_FO_16s_NoIndex_L003                      swab
    ## 2 Bangladesh_Travel_16s_NoIndex_L004                      swab
    ## 3     Bangladesh_FO_16s_NoIndex_L003                      swab
    ## 4     Bangladesh_FO_16s_NoIndex_L003                      swab
    ## 5     Bangladesh_FO_16s_NoIndex_L003                      swab
    ## 6 Bangladesh_Travel_16s_NoIndex_L004                      swab
    ##   sample_center..exp.    Sample.Name       Scientific_Name
    ## 1                CCME SAMEA114622785 human oral metagenome
    ## 2                CCME SAMEA114622786 human skin metagenome
    ## 3                CCME SAMEA114622789 human oral metagenome
    ## 4                CCME SAMEA114622791  human gut metagenome
    ## 5                CCME SAMEA114622793 human oral metagenome
    ## 6                CCME SAMEA114622794 human skin metagenome
    ##     sequencing_meth..exp. SRA.Study target_gene..exp. target_subfragment..exp.
    ## 1 Sequencing by synthesis ERP154464          16S rRNA                       V4
    ## 2 Sequencing by synthesis ERP154464          16S rRNA                       V4
    ## 3 Sequencing by synthesis ERP154464          16S rRNA                       V4
    ## 4 Sequencing by synthesis ERP154464          16S rRNA                       V4
    ## 5 Sequencing by synthesis ERP154464          16S rRNA                       V4
    ## 6 Sequencing by synthesis ERP154464          16S rRNA                       V4
    ##      external_id INSDC_center_name   INSDC_first_public    INSDC_last_update
    ## 1 SAMEA114622785            UCSDMI 2023-11-16T01:29:03Z 2023-11-16T01:29:03Z
    ## 2 SAMEA114622786            UCSDMI 2023-11-16T01:29:03Z 2023-11-16T01:29:03Z
    ## 3 SAMEA114622789            UCSDMI 2023-11-16T01:29:03Z 2023-11-16T01:29:03Z
    ## 4 SAMEA114622791            UCSDMI 2023-11-16T01:29:03Z 2023-11-16T01:29:03Z
    ## 5 SAMEA114622793            UCSDMI 2023-11-16T01:29:03Z 2023-11-16T01:29:03Z
    ## 6 SAMEA114622794            UCSDMI 2023-11-16T01:29:03Z 2023-11-16T01:29:03Z
    ##   INSDC_status                     Sample_name                    Submitter_Id
    ## 1       public  qiita_sid_1015:1015.1.13.M03.M  qiita_sid_1015:1015.1.13.M03.M
    ## 2       public qiita_sid_1015:1015.1.13.M03.RH qiita_sid_1015:1015.1.13.M03.RH
    ## 3       public  qiita_sid_1015:1015.1.14.M03.M  qiita_sid_1015:1015.1.14.M03.M
    ## 4       public  qiita_sid_1015:1015.1.14.M03.S  qiita_sid_1015:1015.1.14.M03.S
    ## 5       public  qiita_sid_1015:1015.1.16.M03.M  qiita_sid_1015:1015.1.16.M03.M
    ## 6       public qiita_sid_1015:1015.1.16.M03.RH qiita_sid_1015:1015.1.16.M03.RH
    ##   organism
    ## 1         
    ## 2         
    ## 3         
    ## 4         
    ## 5         
    ## 6

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
  sample_name = initial_meta$anonymized_name,
  Age = NA,
  BMI = NA,
  Sex = NA,
  Region = initial_meta$target_subfragment..exp.,
  Platform = initial_meta$Instrument,
  Country = initial_meta$COUNTRY,
  Site = initial_meta$Organism,
  public_time = initial_meta$ReleaseDate,
  PCR_primers = initial_meta$pcr_primers..exp.,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$title <- "Ecological relationships between human gut bacteria predicted from analysis of dense microbiome time series data from US travelers in Bangladesh"
metadata$Disease_study <- "NO"
# Assign disease status based on sample_name patterns
metadata$Status <- "HC"

#get subjects id with Saliva available
keep_subjects<-subset(metadata,Site %in% "human oral metagenome")
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

    ## [1] "Number of unique runs: 1472"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 4"

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

    ##   BioProject         Run sample_name Age BMI Sex Region            Platform
    ## 1 PRJEB69530 ERR12266234         M03  NA  NA  NA     V4 Illumina HiSeq 2000
    ## 2 PRJEB69530 ERR12266235         M03  NA  NA  NA     V4 Illumina HiSeq 2000
    ## 3 PRJEB69530 ERR12266238         M03  NA  NA  NA     V4 Illumina HiSeq 2000
    ## 4 PRJEB69530 ERR12266240         M03  NA  NA  NA     V4 Illumina HiSeq 2000
    ## 5 PRJEB69530 ERR12266242         M03  NA  NA  NA     V4 Illumina HiSeq 2000
    ## 6 PRJEB69530 ERR12266243         M03  NA  NA  NA     V4 Illumina HiSeq 2000
    ##   Country                  Site          public_time
    ## 1     USA human oral metagenome 2024-07-21T00:00:00Z
    ## 2     USA human skin metagenome 2024-07-21T00:00:00Z
    ## 3     USA human oral metagenome 2024-07-22T00:00:00Z
    ## 4     USA  human gut metagenome 2024-07-21T00:00:00Z
    ## 5     USA human oral metagenome 2024-07-22T00:00:00Z
    ## 6     USA human skin metagenome 2024-07-22T00:00:00Z
    ##                                         PCR_primers
    ## 1 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 2 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 3 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 4 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 5 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 6 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ##                                                                                                                                              title
    ## 1 Ecological relationships between human gut bacteria predicted from analysis of dense microbiome time series data from US travelers in Bangladesh
    ## 2 Ecological relationships between human gut bacteria predicted from analysis of dense microbiome time series data from US travelers in Bangladesh
    ## 3 Ecological relationships between human gut bacteria predicted from analysis of dense microbiome time series data from US travelers in Bangladesh
    ## 4 Ecological relationships between human gut bacteria predicted from analysis of dense microbiome time series data from US travelers in Bangladesh
    ## 5 Ecological relationships between human gut bacteria predicted from analysis of dense microbiome time series data from US travelers in Bangladesh
    ## 6 Ecological relationships between human gut bacteria predicted from analysis of dense microbiome time series data from US travelers in Bangladesh
    ##   Disease_study Status matched_samples
    ## 1            NO     HC             YES
    ## 2            NO     HC             YES
    ## 3            NO     HC             YES
    ## 4            NO     HC             YES
    ## 5            NO     HC             YES
    ## 6            NO     HC             YES

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
