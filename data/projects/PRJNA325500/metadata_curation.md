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
setwd("D:/oral_project/PRJNA325500")

# Read the initial metadata files
initial_meta<-read.csv("SraRunTable.csv",header = T)

# Combine the two metadata files
#meta <- rbind(initial_meta, initial_meta1)

#Remove entries where Organism is "gut metagenome" (non-oral samples)
#initial_meta <- subset(initial_meta,isolation_source!= "feces")

# Check the number of samples containing "HT" in Sample.Name (e.g., healthy controls)
#nrow(subset(meta, grepl("HT", Sample.Name)))

# Display the first few rows of the metadata
head(initial_meta)
```

    ##          Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample
    ## 1 SRR3667148   AMPLICON        300 27532200 PRJNA325500 SAMN05236073
    ## 2 SRR3667149   AMPLICON        300 23912100 PRJNA325500 SAMN05236174
    ## 3 SRR3667150   AMPLICON        300 24828900 PRJNA325500 SAMN05236057
    ## 4 SRR3667151   AMPLICON        300 23225100 PRJNA325500 SAMN05236166
    ## 5 SRR3667152   AMPLICON        300  8148300 PRJNA325500 SAMN05236046
    ## 6 SRR3667153   AMPLICON        300 21676200 PRJNA325500 SAMN05236111
    ##                         BioSampleModel    Bytes
    ## 1 MIMS.me,MIGS/MIMS/MIMARKS.human-oral 15846043
    ## 2 MIMS.me,MIGS/MIMS/MIMARKS.human-oral 13802706
    ## 3 MIMS.me,MIGS/MIMS/MIMARKS.human-oral 14053138
    ## 4 MIMS.me,MIGS/MIMS/MIMARKS.human-oral 13470490
    ## 5 MIMS.me,MIGS/MIMS/MIMARKS.human-oral  4440886
    ## 6 MIMS.me,MIGS/MIMS/MIMARKS.human-oral 12320991
    ##                                 Center.Name Collection_Date Consent
    ## 1 FEDERAL UNIVERSITY OF MINAS GERAIS - UFMG      2013-02-01  public
    ## 2 FEDERAL UNIVERSITY OF MINAS GERAIS - UFMG      2013-05-13  public
    ## 3 FEDERAL UNIVERSITY OF MINAS GERAIS - UFMG      2013-01-16  public
    ## 4 FEDERAL UNIVERSITY OF MINAS GERAIS - UFMG      2013-05-05  public
    ## 5 FEDERAL UNIVERSITY OF MINAS GERAIS - UFMG      2013-01-05  public
    ## 6 FEDERAL UNIVERSITY OF MINAS GERAIS - UFMG      2013-03-11  public
    ##   DATASTORE.filetype DATASTORE.provider                     DATASTORE.region
    ## 1   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 2   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 3   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 4   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 5   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 6   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ##   Experiment geo_loc_name_country geo_loc_name_country_continent
    ## 1 SRX1845503               Brazil                  South America
    ## 2 SRX1845504               Brazil                  South America
    ## 3 SRX1845505               Brazil                  South America
    ## 4 SRX1845506               Brazil                  South America
    ## 5 SRX1845507               Brazil                  South America
    ## 6 SRX1845508               Brazil                  South America
    ##           geo_loc_name         HOST     Instrument         lat_lon Library.Name
    ## 1 Brazil: Minas Gerais Homo sapiens Illumina MiSeq 19.48 S 43.57 W    A38_cbc39
    ## 2 Brazil: Minas Gerais Homo sapiens Illumina MiSeq 19.48 S 43.57 W    L55_cbc43
    ## 3 Brazil: Minas Gerais Homo sapiens Illumina MiSeq 19.48 S 43.57 W    A22_cbc81
    ## 4 Brazil: Minas Gerais Homo sapiens Illumina MiSeq 19.48 S 43.57 W    L47_cbc10
    ## 5 Brazil: Minas Gerais Homo sapiens Illumina MiSeq 19.48 S 43.57 W   A12_cbc121
    ## 6 Brazil: Minas Gerais Homo sapiens Illumina MiSeq 19.48 S 43.57 W    C44_cbc08
    ##   LibraryLayout LibrarySelection LibrarySource              Organism Platform
    ## 1        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 2        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 3        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 4        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 5        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ## 6        PAIRED              PCR   METAGENOMIC human oral metagenome ILLUMINA
    ##            ReleaseDate Sample.Name SRA.Study                env_biome
    ## 1 2016-06-19T00:00:00Z   A38_cbc39 SRP076567 human-associated habitat
    ## 2 2016-06-19T00:00:00Z   L55_cbc43 SRP076567 human-associated habitat
    ## 3 2016-06-19T00:00:00Z   A22_cbc81 SRP076567 human-associated habitat
    ## 4 2016-06-19T00:00:00Z   L47_cbc10 SRP076567 human-associated habitat
    ## 5 2016-06-19T00:00:00Z  A12_cbc121 SRP076567 human-associated habitat
    ## 6 2016-06-19T00:00:00Z   C44_cbc08 SRP076567 human-associated habitat
    ##   env_feature env_material          create_date version isolate sample_type
    ## 1        oral      biofilm 2016-06-14T17:36:00Z       1                    
    ## 2        oral      biofilm 2016-06-14T17:36:00Z       1                    
    ## 3        oral      biofilm 2016-06-14T17:36:00Z       1                    
    ## 4        oral      biofilm 2016-06-14T17:36:00Z       1                    
    ## 5        oral      biofilm 2016-06-14T17:36:00Z       1                    
    ## 6        oral      biofilm 2016-06-14T17:36:00Z       1

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
  Sex=NA,
  Region = "V4",
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  Site = "subgingival plaque",
  public_time = initial_meta$ReleaseDate,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$PCR_primers <- "515F-806R"
metadata$title <- "Subgingival microbiota dysbiosis in systemic lupus erythematosus association with periodontal status"
#metadata<-subset(uniform_meta,Site!="feces")
metadata$Disease_study <- "YES"
# Assign disease status based on sample_name patterns
metadata$Status <- NA
#extract subject id by remove prefix of sample_name
#metadata$sample_name <- gsub("[FSP].*", "", metadata$sample_name)
#get subjects id with Saliva available
#keep_subjects<-subset(metadata,Site %in% c("human hard palate swab","human Keratinized gingiva swab","human oral wash","human saliva","human Supragingival plaque scrape","human tongue dorsum swab"))
#remove subjects without Saliva samples
#metadata<-subset(metadata,sample_name %in% keep_subjects$sample_name)
#get subjects sampling counts
#sample_counts <- table(metadata$sample_name)
# Create matched_samples column, any subjects have >=2 samples will be label "YES", otherwise "NO"
#metadata$matched_samples <- ifelse(metadata$sample_name %in% names(sample_counts[sample_counts >= 2]), "YES", "NO")

#metadata$Site<-gsub("human ","",metadata$Site)
metadata$matched_samples<-"NO"
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 141"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 141"

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

    ##    BioProject        Run sample_name Age BMI Sex Region       Platform Country
    ## 1 PRJNA325500 SRR3667148   A38_cbc39  NA  NA  NA     V4 Illumina MiSeq  Brazil
    ## 2 PRJNA325500 SRR3667149   L55_cbc43  NA  NA  NA     V4 Illumina MiSeq  Brazil
    ## 3 PRJNA325500 SRR3667150   A22_cbc81  NA  NA  NA     V4 Illumina MiSeq  Brazil
    ## 4 PRJNA325500 SRR3667151   L47_cbc10  NA  NA  NA     V4 Illumina MiSeq  Brazil
    ## 5 PRJNA325500 SRR3667152  A12_cbc121  NA  NA  NA     V4 Illumina MiSeq  Brazil
    ## 6 PRJNA325500 SRR3667153   C44_cbc08  NA  NA  NA     V4 Illumina MiSeq  Brazil
    ##                 Site          public_time PCR_primers
    ## 1 subgingival plaque 2016-06-19T00:00:00Z   515F-806R
    ## 2 subgingival plaque 2016-06-19T00:00:00Z   515F-806R
    ## 3 subgingival plaque 2016-06-19T00:00:00Z   515F-806R
    ## 4 subgingival plaque 2016-06-19T00:00:00Z   515F-806R
    ## 5 subgingival plaque 2016-06-19T00:00:00Z   515F-806R
    ## 6 subgingival plaque 2016-06-19T00:00:00Z   515F-806R
    ##                                                                                                  title
    ## 1 Subgingival microbiota dysbiosis in systemic lupus erythematosus association with periodontal status
    ## 2 Subgingival microbiota dysbiosis in systemic lupus erythematosus association with periodontal status
    ## 3 Subgingival microbiota dysbiosis in systemic lupus erythematosus association with periodontal status
    ## 4 Subgingival microbiota dysbiosis in systemic lupus erythematosus association with periodontal status
    ## 5 Subgingival microbiota dysbiosis in systemic lupus erythematosus association with periodontal status
    ## 6 Subgingival microbiota dysbiosis in systemic lupus erythematosus association with periodontal status
    ##   Disease_study Status matched_samples
    ## 1           YES     NA              NO
    ## 2           YES     NA              NO
    ## 3           YES     NA              NO
    ## 4           YES     NA              NO
    ## 5           YES     NA              NO
    ## 6           YES     NA              NO

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
