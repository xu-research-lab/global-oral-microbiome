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
setwd("D:/oral_project/PRJNA750937")

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

    ##           Run      AGE Assay.Type AvgSpotLen    Bases
    ## 1 SRR15307426 48 years   AMPLICON        602 19548144
    ## 2 SRR15307427 35 years   AMPLICON        602 34841954
    ## 3 SRR15307428 40 years   AMPLICON        602 31089688
    ## 4 SRR15307429 23 years   AMPLICON        602 30522002
    ## 5 SRR15307430 50 years   AMPLICON        602 40086578
    ## 6 SRR15307431 49 years   AMPLICON        602 22103032
    ##                                                                                                                                                                                      BIOMATERIAL_PROVIDER
    ## 1 Neurosurgery Clinical Medical Research Center of Sichuan Province\\, Prof. Dr. Ligang Chen\\, Department of Neurosurgery\\, the Affiliated Hospital of Southwest Medical University\\, Lu Zhou\\, China
    ## 2 Neurosurgery Clinical Medical Research Center of Sichuan Province\\, Prof. Dr. Ligang Chen\\, Department of Neurosurgery\\, the Affiliated Hospital of Southwest Medical University\\, Lu Zhou\\, China
    ## 3 Neurosurgery Clinical Medical Research Center of Sichuan Province\\, Prof. Dr. Ligang Chen\\, Department of Neurosurgery\\, the Affiliated Hospital of Southwest Medical University\\, Lu Zhou\\, China
    ## 4 Neurosurgery Clinical Medical Research Center of Sichuan Province\\, Prof. Dr. Ligang Chen\\, Department of Neurosurgery\\, the Affiliated Hospital of Southwest Medical University\\, Lu Zhou\\, China
    ## 5 Neurosurgery Clinical Medical Research Center of Sichuan Province\\, Prof. Dr. Ligang Chen\\, Department of Neurosurgery\\, the Affiliated Hospital of Southwest Medical University\\, Lu Zhou\\, China
    ## 6 Neurosurgery Clinical Medical Research Center of Sichuan Province\\, Prof. Dr. Ligang Chen\\, Department of Neurosurgery\\, the Affiliated Hospital of Southwest Medical University\\, Lu Zhou\\, China
    ##    BioProject    BioSample BioSampleModel    Bytes
    ## 1 PRJNA750937 SAMN20502004          Human 11152744
    ## 2 PRJNA750937 SAMN20502022          Human 20083556
    ## 3 PRJNA750937 SAMN20502026          Human 17284858
    ## 4 PRJNA750937 SAMN20502027          Human 16946039
    ## 5 PRJNA750937 SAMN20502031          Human 22154301
    ## 6 PRJNA750937 SAMN20502005          Human 12923609
    ##                                               Center.Name Consent
    ## 1 THE AFFILIATED HOSPITAL OF SOUTHWEST MEDICAL UNIVERSITY  public
    ## 2 THE AFFILIATED HOSPITAL OF SOUTHWEST MEDICAL UNIVERSITY  public
    ## 3 THE AFFILIATED HOSPITAL OF SOUTHWEST MEDICAL UNIVERSITY  public
    ## 4 THE AFFILIATED HOSPITAL OF SOUTHWEST MEDICAL UNIVERSITY  public
    ## 5 THE AFFILIATED HOSPITAL OF SOUTHWEST MEDICAL UNIVERSITY  public
    ## 6 THE AFFILIATED HOSPITAL OF SOUTHWEST MEDICAL UNIVERSITY  public
    ##   DATASTORE.filetype DATASTORE.provider                     DATASTORE.region
    ## 1   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 2   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 3   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 4   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 5   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 6   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ##    Experiment health_state     Instrument                   isolate
    ## 1 SRX11612107       health Illumina MiSeq   health control 3 saliva
    ## 2 SRX11612106       health Illumina MiSeq  health control 21 saliva
    ## 3 SRX11612105      disease Illumina MiSeq low glioma group 1 saliva
    ## 4 SRX11612104      disease Illumina MiSeq low glioma group 2 saliva
    ## 5 SRX11612103      disease Illumina MiSeq low glioma group 6 saliva
    ## 6 SRX11612102       health Illumina MiSeq   health control 4 saliva
    ##      Library.Name LibraryLayout LibrarySelection LibrarySource     Organism
    ## 1  ="sample  "&A4        PAIRED              PCR   METAGENOMIC Homo sapiens
    ## 2 ="sample  "&A22        PAIRED              PCR   METAGENOMIC Homo sapiens
    ## 3 ="sample  "&A26        PAIRED              PCR   METAGENOMIC Homo sapiens
    ## 4 ="sample  "&A27        PAIRED              PCR   METAGENOMIC Homo sapiens
    ## 5 ="sample  "&A31        PAIRED              PCR   METAGENOMIC Homo sapiens
    ## 6  ="sample  "&A5        PAIRED              PCR   METAGENOMIC Homo sapiens
    ##   Platform population          ReleaseDate          create_date version
    ## 1 ILLUMINA  for human 2021-10-01T00:00:00Z 2021-08-01T05:35:00Z       1
    ## 2 ILLUMINA  for human 2021-10-01T00:00:00Z 2021-08-01T05:35:00Z       1
    ## 3 ILLUMINA  for human 2021-10-01T00:00:00Z 2021-08-01T05:35:00Z       1
    ## 4 ILLUMINA  for human 2021-10-01T00:00:00Z 2021-08-08T21:45:00Z       1
    ## 5 ILLUMINA  for human 2021-10-01T00:00:00Z 2021-08-08T21:46:00Z       1
    ## 6 ILLUMINA  for human 2021-10-01T00:00:00Z 2021-08-08T21:46:00Z       1
    ##   Sample.Name    sex SRA.Study tissue sample_type
    ## 1         T71   male SRP330605 saliva            
    ## 2         T91   male SRP330605 saliva            
    ## 3          T2   male SRP330605 saliva            
    ## 4          T7   male SRP330605 saliva            
    ## 5         T25   male SRP330605 saliva            
    ## 6         T72 female SRP330605 saliva

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
  Age = initial_meta$AGE,
  BMI = NA,
  Sex=initial_meta$sex,
  Region = "V3-V4",
  Platform = initial_meta$Instrument,
  Country = "China",
  Site = initial_meta$tissue,
  public_time = initial_meta$ReleaseDate,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$PCR_primers <- "338F-806R"
metadata$title <- "Association Between Oral Microbiota and Human Brain Glioma Grade A Case-Control Study"
#metadata<-subset(uniform_meta,Site!="feces")
metadata$Disease_study <- "NO"
# Assign disease status based on sample_name patterns
metadata$Status <- initial_meta$health_state
#extract subject id by remove prefix of sample_name
#metadata$sample_name <- gsub("[FSP].*", "", metadata$sample_name)
#get subjects id with Saliva available
keep_subjects<-subset(metadata,Site %in% c("saliva"))
#remove subjects without Saliva samples
metadata<-subset(metadata,sample_name %in% keep_subjects$sample_name)
#get subjects sampling counts
sample_counts <- table(metadata$sample_name)
# Create matched_samples column, any subjects have >=2 samples will be label "YES", otherwise "NO"
metadata$matched_samples <- ifelse(metadata$sample_name %in% names(sample_counts[sample_counts >= 2]), "YES", "NO")
metadata$Age<-gsub(" years","",metadata$Age)
metadata$Age<-as.numeric(metadata$Age)
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 59"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 59"

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

    ##    BioProject         Run sample_name Age BMI    Sex Region       Platform
    ## 1 PRJNA750937 SRR15307426         T71  48  NA   male  V3-V4 Illumina MiSeq
    ## 2 PRJNA750937 SRR15307427         T91  35  NA   male  V3-V4 Illumina MiSeq
    ## 3 PRJNA750937 SRR15307428          T2  40  NA   male  V3-V4 Illumina MiSeq
    ## 4 PRJNA750937 SRR15307429          T7  23  NA   male  V3-V4 Illumina MiSeq
    ## 5 PRJNA750937 SRR15307430         T25  50  NA   male  V3-V4 Illumina MiSeq
    ## 6 PRJNA750937 SRR15307431         T72  49  NA female  V3-V4 Illumina MiSeq
    ##   Country   Site          public_time PCR_primers
    ## 1   China saliva 2021-10-01T00:00:00Z   338F-806R
    ## 2   China saliva 2021-10-01T00:00:00Z   338F-806R
    ## 3   China saliva 2021-10-01T00:00:00Z   338F-806R
    ## 4   China saliva 2021-10-01T00:00:00Z   338F-806R
    ## 5   China saliva 2021-10-01T00:00:00Z   338F-806R
    ## 6   China saliva 2021-10-01T00:00:00Z   338F-806R
    ##                                                                                   title
    ## 1 Association Between Oral Microbiota and Human Brain Glioma Grade A Case-Control Study
    ## 2 Association Between Oral Microbiota and Human Brain Glioma Grade A Case-Control Study
    ## 3 Association Between Oral Microbiota and Human Brain Glioma Grade A Case-Control Study
    ## 4 Association Between Oral Microbiota and Human Brain Glioma Grade A Case-Control Study
    ## 5 Association Between Oral Microbiota and Human Brain Glioma Grade A Case-Control Study
    ## 6 Association Between Oral Microbiota and Human Brain Glioma Grade A Case-Control Study
    ##   Disease_study  Status matched_samples
    ## 1            NO  health              NO
    ## 2            NO  health              NO
    ## 3            NO disease              NO
    ## 4            NO disease              NO
    ## 5            NO disease              NO
    ## 6            NO  health              NO

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
