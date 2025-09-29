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
setwd("D:/studylearning/董彪师兄项目/任务/张辰飞/PRJNA666891")

# Read the initial metadata files
initial_meta<-read.csv("SraRunTable.csv",header = T)

# Combine the two metadata files
#meta <- rbind(initial_meta, initial_meta1)

#Remove WGS samples
initial_meta <- subset(initial_meta,Assay.Type == "AMPLICON")

# Check the number of samples containing "HT" in Sample.Name (e.g., healthy controls)
#nrow(subset(meta, grepl("HT", Sample.Name)))

# Display the first few rows of the metadata
head(initial_meta)
```

    ##           Run Assay.Type AvgSpotLen  Bases  BioProject    BioSample
    ## 4 SRR12759016   AMPLICON        602  25886 PRJNA666891 SAMN16330834
    ## 5 SRR12759017   AMPLICON        602  21070 PRJNA666891 SAMN16330833
    ## 6 SRR12759018   AMPLICON        602 136654 PRJNA666891 SAMN16330832
    ## 7 SRR12759019   AMPLICON        602  35518 PRJNA666891 SAMN16330831
    ## 8 SRR12759020   AMPLICON        602  25886 PRJNA666891 SAMN16330830
    ## 9 SRR12759021   AMPLICON        602  97524 PRJNA666891 SAMN16330829
    ##                                BioSampleModel  Bytes               Center.Name
    ## 4 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral  93813 J. CRAIG VENTER INSTITUTE
    ## 5 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral  90917 J. CRAIG VENTER INSTITUTE
    ## 6 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 149040 J. CRAIG VENTER INSTITUTE
    ## 7 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral  98731 J. CRAIG VENTER INSTITUTE
    ## 8 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral  94036 J. CRAIG VENTER INSTITUTE
    ## 9 MIMARKS.survey,MIGS/MIMS/MIMARKS.human-oral 132133 J. CRAIG VENTER INSTITUTE
    ##   Collection_Date Consent DATASTORE.filetype DATASTORE.provider
    ## 4      2013-10-30  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5      2013-10-30  public   run.zq,sra,fastq         gs,ncbi,s3
    ## 6      2013-10-30  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 7      2013-10-30  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 8      2013-10-30  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 9      2013-10-30  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                       DATASTORE.region                  env_broad_scale
    ## 4 gs.us-east1,ncbi.public,s3.us-east-1 Human microbiome\\, 16S Illumina
    ## 5 gs.us-east1,ncbi.public,s3.us-east-1 Human microbiome\\, 16S Illumina
    ## 6 gs.us-east1,ncbi.public,s3.us-east-1 Human microbiome\\, 16S Illumina
    ## 7 gs.us-east1,ncbi.public,s3.us-east-1 Human microbiome\\, 16S Illumina
    ## 8 gs.us-east1,ncbi.public,s3.us-east-1 Human microbiome\\, 16S Illumina
    ## 9 gs.us-east1,ncbi.public,s3.us-east-1 Human microbiome\\, 16S Illumina
    ##         env_local_scale        env_medium Experiment geo_loc_name_country
    ## 4 Human oral metagenome Human oral cavity SRX9229868               Poland
    ## 5 Human oral metagenome Human oral cavity SRX9229867               Poland
    ## 6 Human oral metagenome Human oral cavity SRX9229866               Poland
    ## 7 Human oral metagenome Human oral cavity SRX9229865               Poland
    ## 8 Human oral metagenome Human oral cavity SRX9229864               Poland
    ## 9 Human oral metagenome Human oral cavity SRX9229863               Poland
    ##   geo_loc_name_country_continent                                   geo_loc_name
    ## 4                         Europe Poland: Poznan\\, Greater Poland Cancer Center
    ## 5                         Europe Poland: Poznan\\, Greater Poland Cancer Center
    ## 6                         Europe Poland: Poznan\\, Greater Poland Cancer Center
    ## 7                         Europe Poland: Poznan\\, Greater Poland Cancer Center
    ## 8                         Europe Poland: Poznan\\, Greater Poland Cancer Center
    ## 9                         Europe Poland: Poznan\\, Greater Poland Cancer Center
    ##           HOST host_subject_id     Instrument             lat_lon Library.Name
    ## 4 Homo sapiens        P0028_TT Illumina MiSeq 52.4020 N 16.9359 E     P0028_TT
    ## 5 Homo sapiens       P0028_SAL Illumina MiSeq 52.4020 N 16.9359 E    P0028_SAL
    ## 6 Homo sapiens        P0028_HT Illumina MiSeq 52.4020 N 16.9359 E     P0028_HT
    ## 7 Homo sapiens        P0027_TT Illumina MiSeq 52.4020 N 16.9359 E     P0027_TT
    ## 8 Homo sapiens       P0027_SAL Illumina MiSeq 52.4020 N 16.9359 E    P0027_SAL
    ## 9 Homo sapiens        P0027_HT Illumina MiSeq 52.4020 N 16.9359 E     P0027_HT
    ##   LibraryLayout LibrarySelection LibrarySource   Organism Platform
    ## 4        PAIRED              PCR   METAGENOMIC metagenome ILLUMINA
    ## 5        PAIRED              PCR   METAGENOMIC metagenome ILLUMINA
    ## 6        PAIRED              PCR   METAGENOMIC metagenome ILLUMINA
    ## 7        PAIRED              PCR   METAGENOMIC metagenome ILLUMINA
    ## 8        PAIRED              PCR   METAGENOMIC metagenome ILLUMINA
    ## 9        PAIRED              PCR   METAGENOMIC metagenome ILLUMINA
    ##            ReleaseDate          create_date version Sample.Name SRA.Study
    ## 4 2020-10-01T00:00:00Z 2020-10-01T11:32:00Z       1    P0028_TT SRP286018
    ## 5 2020-10-01T00:00:00Z 2020-10-01T11:32:00Z       1   P0028_SAL SRP286018
    ## 6 2020-10-01T00:00:00Z 2020-10-01T11:32:00Z       1    P0028_HT SRP286018
    ## 7 2020-10-01T00:00:00Z 2020-10-01T11:32:00Z       1    P0027_TT SRP286018
    ## 8 2020-10-01T00:00:00Z 2020-10-01T11:32:00Z       1   P0027_SAL SRP286018
    ## 9 2020-10-01T00:00:00Z 2020-10-01T11:32:00Z       1    P0027_HT SRP286018
    ##   host_tissue_sampled
    ## 4                    
    ## 5                    
    ## 6                    
    ## 7                    
    ## 8                    
    ## 9

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
  Sex= NA,
  Region = "V4",
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  Site = NA,
  public_time = initial_meta$ReleaseDate,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$PCR_primers <- NA_character_
metadata$title <- "Oral Microbial Species and Virulence Factors Associated with Oral Squamous Cell Carcinoma"
metadata$Disease_study <- "YES"
metadata$Status <- "OSCC"
# Assign disease status based on sample_name patterns
metadata <- metadata %>% 
  mutate(Site = case_when(
    grepl("SAL",sample_name) ~ "saliva",
    grepl("TT|G|Z",sample_name) ~ "tumor tissue",
    grepl("HT",sample_name) ~ "non-tumor tissue"
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
    ##  P0001_HT P0001_SAL  P0001_TT  P0002_HT P0002_SAL  P0002_TT  P0003_HT P0003_SAL 
    ##         1         1         1         1         1         1         1         1 
    ##  P0003_TT  P0004_HT P0004_SAL  P0004_TT  P0005_HT P0005_SAL  P0005_TT  P0006_HT 
    ##         1         1         1         1         1         1         1         1 
    ## P0006_SAL  P0006_TT  P0007_HT P0007_SAL  P0007_TT  P0008_HT P0008_SAL  P0008_TT 
    ##         1         1         1         1         1         1         1         1 
    ##  P0013_HT P0013_SAL  P0013_TT  P0016_HT P0016_SAL  P0016_TT  P0017_HT P0017_SAL 
    ##         1         1         1         1         1         1         1         1 
    ##  P0017_TT  P0018_HT P0018_SAL  P0018_TT  P0019_HT P0019_SAL  P0019_TT  P0020_HT 
    ##         1         1         1         1         1         1         1         1 
    ## P0020_SAL  P0020_TT  P0024_HT P0024_SAL  P0024_TT  P0025_HT P0025_SAL  P0025_TT 
    ##         1         1         1         1         1         1         1         1 
    ##  P0027_HT P0027_SAL  P0027_TT  P0028_HT P0028_SAL  P0028_TT 
    ##         1         1         1         1         1         1

``` r
table(metadata$Site)
```

    ## 
    ## non-tumor tissue           saliva     tumor tissue 
    ##               18               18               18

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
    ##  P0001_HT P0001_SAL  P0001_TT  P0002_HT P0002_SAL  P0002_TT  P0003_HT P0003_SAL 
    ##         1         1         1         1         1         1         1         1 
    ##  P0003_TT  P0004_HT P0004_SAL  P0004_TT  P0005_HT P0005_SAL  P0005_TT  P0006_HT 
    ##         1         1         1         1         1         1         1         1 
    ## P0006_SAL  P0006_TT  P0007_HT P0007_SAL  P0007_TT  P0008_HT P0008_SAL  P0008_TT 
    ##         1         1         1         1         1         1         1         1 
    ##  P0013_HT P0013_SAL  P0013_TT  P0016_HT P0016_SAL  P0016_TT  P0017_HT P0017_SAL 
    ##         1         1         1         1         1         1         1         1 
    ##  P0017_TT  P0018_HT P0018_SAL  P0018_TT  P0019_HT P0019_SAL  P0019_TT  P0020_HT 
    ##         1         1         1         1         1         1         1         1 
    ## P0020_SAL  P0020_TT  P0024_HT P0024_SAL  P0024_TT  P0025_HT P0025_SAL  P0025_TT 
    ##         1         1         1         1         1         1         1         1 
    ##  P0027_HT P0027_SAL  P0027_TT  P0028_HT P0028_SAL  P0028_TT 
    ##         1         1         1         1         1         1

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
    ## 54

``` r
# Identify samples for removal: single samples excluding saliva,tumor tissue and non-tumor tissue
single_samples <- names(sample_counts[sample_counts == 1])
samples_to_remove <- metadata$sample_name %in% single_samples & !(metadata$Site %in% c("saliva", "tumor tissue","non-tumor tissue"))

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
    ##  P0001_HT P0001_SAL  P0001_TT  P0002_HT P0002_SAL  P0002_TT  P0003_HT P0003_SAL 
    ##         1         1         1         1         1         1         1         1 
    ##  P0003_TT  P0004_HT P0004_SAL  P0004_TT  P0005_HT P0005_SAL  P0005_TT  P0006_HT 
    ##         1         1         1         1         1         1         1         1 
    ## P0006_SAL  P0006_TT  P0007_HT P0007_SAL  P0007_TT  P0008_HT P0008_SAL  P0008_TT 
    ##         1         1         1         1         1         1         1         1 
    ##  P0013_HT P0013_SAL  P0013_TT  P0016_HT P0016_SAL  P0016_TT  P0017_HT P0017_SAL 
    ##         1         1         1         1         1         1         1         1 
    ##  P0017_TT  P0018_HT P0018_SAL  P0018_TT  P0019_HT P0019_SAL  P0019_TT  P0020_HT 
    ##         1         1         1         1         1         1         1         1 
    ## P0020_SAL  P0020_TT  P0024_HT P0024_SAL  P0024_TT  P0025_HT P0025_SAL  P0025_TT 
    ##         1         1         1         1         1         1         1         1 
    ##  P0027_HT P0027_SAL  P0027_TT  P0028_HT P0028_SAL  P0028_TT 
    ##         1         1         1         1         1         1

``` r
table(metadata_clean$Site)
```

    ## 
    ## non-tumor tissue           saliva     tumor tissue 
    ##               18               18               18

``` r
table(metadata_clean$matched_samples)
```

    ## 
    ## NO 
    ## 54

``` r
# Use cleaned dataset for subsequent analysis
metadata <- metadata_clean
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

    ## [1] "Number of unique sample names: 54"

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

    ##    BioProject         Run sample_name Age BMI Sex Region       Platform Country
    ## 1 PRJNA666891 SRR12759016    P0028_TT  NA  NA  NA     V4 Illumina MiSeq  Poland
    ## 2 PRJNA666891 SRR12759017   P0028_SAL  NA  NA  NA     V4 Illumina MiSeq  Poland
    ## 3 PRJNA666891 SRR12759018    P0028_HT  NA  NA  NA     V4 Illumina MiSeq  Poland
    ## 4 PRJNA666891 SRR12759019    P0027_TT  NA  NA  NA     V4 Illumina MiSeq  Poland
    ## 5 PRJNA666891 SRR12759020   P0027_SAL  NA  NA  NA     V4 Illumina MiSeq  Poland
    ## 6 PRJNA666891 SRR12759021    P0027_HT  NA  NA  NA     V4 Illumina MiSeq  Poland
    ##               Site          public_time PCR_primers
    ## 1     tumor tissue 2020-10-01T00:00:00Z        <NA>
    ## 2           saliva 2020-10-01T00:00:00Z        <NA>
    ## 3 non-tumor tissue 2020-10-01T00:00:00Z        <NA>
    ## 4     tumor tissue 2020-10-01T00:00:00Z        <NA>
    ## 5           saliva 2020-10-01T00:00:00Z        <NA>
    ## 6 non-tumor tissue 2020-10-01T00:00:00Z        <NA>
    ##                                                                                       title
    ## 1 Oral Microbial Species and Virulence Factors Associated with Oral Squamous Cell Carcinoma
    ## 2 Oral Microbial Species and Virulence Factors Associated with Oral Squamous Cell Carcinoma
    ## 3 Oral Microbial Species and Virulence Factors Associated with Oral Squamous Cell Carcinoma
    ## 4 Oral Microbial Species and Virulence Factors Associated with Oral Squamous Cell Carcinoma
    ## 5 Oral Microbial Species and Virulence Factors Associated with Oral Squamous Cell Carcinoma
    ## 6 Oral Microbial Species and Virulence Factors Associated with Oral Squamous Cell Carcinoma
    ##   Disease_study Status matched_samples
    ## 1           YES   OSCC              NO
    ## 2           YES   OSCC              NO
    ## 3           YES   OSCC              NO
    ## 4           YES   OSCC              NO
    ## 5           YES   OSCC              NO
    ## 6           YES   OSCC              NO

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
