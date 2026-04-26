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
library(readxl)

# Set working directory (update this path to your local directory)
setwd("/home/dongbiao/word_embedding_microbiome/all_data/oral/global-oral-microbiome/data/projects/OEP000837")

# Read the initial metadata files
initial_meta<-read_xlsx("OEP00000837_Data_1758007425482.xlsx",sheet = "sheet1")

# Combine the two metadata files
#meta <- rbind(initial_meta, initial_meta1)

# Remove entries where Organism is "gut metagenome" (non-oral samples)
#meta <- subset(meta, Organism != "gut metagenome")

# Check the number of samples containing "HT" in Sample.Name (e.g., healthy controls)
#nrow(subset(meta, grepl("HT", Sample.Name)))

# Display the first few rows of the metadata
head(initial_meta)
```

    ## # A tibble: 6 × 14
    ##   experiment_id experiment_name            experiment_type sample_id sample_name
    ##   <chr>         <chr>                      <chr>           <chr>     <chr>      
    ## 1 OEX00003132   the characteristics of or… Metagenomic     OES00028… S001-T     
    ## 2 OEX00003132   the characteristics of or… Metagenomic     OES00028… S001-C     
    ## 3 OEX00003132   the characteristics of or… Metagenomic     OES00028… S002-TY    
    ## 4 OEX00003132   the characteristics of or… Metagenomic     OES00028… S001-T     
    ## 5 OEX00003132   the characteristics of or… Metagenomic     OES00028… S001-C     
    ## 6 OEX00003132   the characteristics of or… Metagenomic     OES00028… S002-TY    
    ## # ℹ 9 more variables: sample_type <chr>, run_id <chr>, run_name <chr>,
    ## #   data_id <chr>, data_name <chr>, data_type <chr>, security <chr>, md5 <chr>,
    ## #   data_size <chr>

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
  BioProject = "OEP000837",
  Run = initial_meta$sample_name,
  sample_name = initial_meta$sample_name,
  Age = NA,
  BMI = NA,
  Sex=NA,
  Region = "V3-V4",
  Platform = "Illumina MiSeq",
  Country = "China",
  #Site = "Saliva",
  public_time = "2020-04-14",
  stringsAsFactors = FALSE
)
metadata<-unique(metadata)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$PCR_primers <- "338F-806R"
metadata$title <- "Oral Microbiota Analysis of Tissue Pairs and Saliva Samples From Patients With Oral Squamous Cell Carcinoma A Pilot Study"
metadata$Disease_study <- "YES"

# Assign disease status based on sample_name patterns
metadata<-metadata%>%
  mutate(Status=case_when(grepl("TY",sample_name)~"OSCC_Saliva",grepl("T|t",sample_name)~"OSCC_Tumor",grepl("C",sample_name)~"OSCC_Paracancerous Tissue"))

#get subjects id according to sample id
metadata<-metadata%>%
  mutate(Site=case_when(grepl("TY",sample_name)~"Saliva",grepl("T|t",sample_name)~"Tumor",grepl("C",sample_name)~"Paracancerous Tissue"))

metadata$sample_name<-gsub("-C|-T|-TY|-t","",metadata$sample_name)
#get subjects id with Saliva available
keep_subjects<-subset(metadata,Site %in% "Saliva")
#remove subjects without Saliva samples
metadata<-subset(metadata,sample_name %in% keep_subjects$sample_name)
#get subjects sampling counts
sample_counts <- table(metadata$sample_name)
# Create matched_samples column, any subjects have >=2 samples will be label "YES", otherwise "NO"
metadata$matched_samples <- ifelse(metadata$sample_name %in% names(sample_counts[sample_counts >= 2]), "YES", "NO")
rownames(metadata)<-1:nrow(metadata)
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 86"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 44"

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
metadata<-unique(metadata)
write.csv(metadata, file = "metadata.csv", quote = FALSE)
head(metadata)
```

    ##   BioProject     Run sample_name Age BMI Sex Region       Platform Country
    ## 1  OEP000837 S002-TY        S002  NA  NA  NA  V3-V4 Illumina MiSeq   China
    ## 2  OEP000837 S003-TY        S003  NA  NA  NA  V3-V4 Illumina MiSeq   China
    ## 3  OEP000837 S004-TY        S004  NA  NA  NA  V3-V4 Illumina MiSeq   China
    ## 4  OEP000837  S006-C        S006  NA  NA  NA  V3-V4 Illumina MiSeq   China
    ## 5  OEP000837 S005-TY        S005  NA  NA  NA  V3-V4 Illumina MiSeq   China
    ## 6  OEP000837 S006-TY        S006  NA  NA  NA  V3-V4 Illumina MiSeq   China
    ##                   Site public_time PCR_primers
    ## 1               Saliva  2020-04-14   338F-806R
    ## 2               Saliva  2020-04-14   338F-806R
    ## 3               Saliva  2020-04-14   338F-806R
    ## 4 Paracancerous Tissue  2020-04-14   338F-806R
    ## 5               Saliva  2020-04-14   338F-806R
    ## 6               Saliva  2020-04-14   338F-806R
    ##                                                                                                                       title
    ## 1 Oral Microbiota Analysis of Tissue Pairs and Saliva Samples From Patients With Oral Squamous Cell Carcinoma A Pilot Study
    ## 2 Oral Microbiota Analysis of Tissue Pairs and Saliva Samples From Patients With Oral Squamous Cell Carcinoma A Pilot Study
    ## 3 Oral Microbiota Analysis of Tissue Pairs and Saliva Samples From Patients With Oral Squamous Cell Carcinoma A Pilot Study
    ## 4 Oral Microbiota Analysis of Tissue Pairs and Saliva Samples From Patients With Oral Squamous Cell Carcinoma A Pilot Study
    ## 5 Oral Microbiota Analysis of Tissue Pairs and Saliva Samples From Patients With Oral Squamous Cell Carcinoma A Pilot Study
    ## 6 Oral Microbiota Analysis of Tissue Pairs and Saliva Samples From Patients With Oral Squamous Cell Carcinoma A Pilot Study
    ##   Disease_study                    Status matched_samples
    ## 1           YES               OSCC_Saliva              NO
    ## 2           YES               OSCC_Saliva              NO
    ## 3           YES               OSCC_Saliva              NO
    ## 4           YES OSCC_Paracancerous Tissue             YES
    ## 5           YES               OSCC_Saliva              NO
    ## 6           YES               OSCC_Saliva             YES

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
