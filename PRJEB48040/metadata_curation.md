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
setwd("E:/bioproject/PRJEB48040")

# Read the initial metadata files
initial_meta <- read.csv("SraRunTable.csv", header = TRUE)
initial_meta <- initial_meta[initial_meta$env_material != "", ]

# Display the first few rows of the metadata
head(initial_meta)
```

    ##           Run Assay.Type AvgSpotLen barcode..exp.   Bases BioProject
    ## 5  ERR7013439   AMPLICON        150  GTTGACCATCGC 3371250 PRJEB48040
    ## 6  ERR7013440   AMPLICON        150  AGACACCAATGT 5207100 PRJEB48040
    ## 7  ERR7013441   AMPLICON        150  TCCTAGGTCCGA 3781050 PRJEB48040
    ## 8  ERR7013442   AMPLICON        150  CGCAGATTAGTA 6031950 PRJEB48040
    ## 9  ERR7013443   AMPLICON        150  GGCCTATAAGTC 3344700 PRJEB48040
    ## 10 ERR7013444   AMPLICON        150  CTCTATTCCACC 2933550 PRJEB48040
    ##        BioSample   Bytes  center_name..exp.
    ## 5  SAMEA10332186 2161160 Rutgers University
    ## 6  SAMEA10332187 3444051 Rutgers University
    ## 7  SAMEA10332188 2427374 Rutgers University
    ## 8  SAMEA10332189 3700523 Rutgers University
    ## 9  SAMEA10332190 2232579 Rutgers University
    ## 10 SAMEA10332191 1942843 Rutgers University
    ##                                           Center.Name center_project_name..exp.
    ## 5  UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT      Saliva methods pilot
    ## 6  UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT      Saliva methods pilot
    ## 7  UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT      Saliva methods pilot
    ## 8  UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT      Saliva methods pilot
    ## 9  UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT      Saliva methods pilot
    ## 10 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT      Saliva methods pilot
    ##    collection_timestamp colletion_method Consent DATASTORE.filetype
    ## 5               9/16/20  No preservative  public   fastq,run.zq,sra
    ## 6               9/16/20     Preservative  public   fastq,run.zq,sra
    ## 7               9/16/20     Preservative  public   fastq,run.zq,sra
    ## 8               9/16/20  No preservative  public   fastq,run.zq,sra
    ## 9               9/16/20  No preservative  public   fastq,run.zq,sra
    ## 10              9/16/20     Preservative  public   fastq,run.zq,sra
    ##    DATASTORE.provider                     DATASTORE.region dna_extracted
    ## 5          gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1          TRUE
    ## 6          gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1          TRUE
    ## 7          gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1          TRUE
    ## 8          gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1          TRUE
    ## 9          gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1          TRUE
    ## 10         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1          TRUE
    ##    elevation    emp_status          empo_1            empo_2           empo_3
    ## 5         21 EMP_Processed Host-associated Animal-associated Animal secretion
    ## 6         21 EMP_Processed Host-associated Animal-associated Animal secretion
    ## 7         21 EMP_Processed Host-associated Animal-associated Animal secretion
    ## 8         21 EMP_Processed Host-associated Animal-associated Animal secretion
    ## 9         21 EMP_Processed Host-associated Animal-associated Animal secretion
    ## 10        21 EMP_Processed Host-associated Animal-associated Animal secretion
    ##    ENA.FIRST.PUBLIC..run. ENA_first_public ENA.LAST.UPDATE..run.
    ## 5              2021-11-17       2021-11-17            2021-11-17
    ## 6              2021-11-17       2021-11-17            2021-11-17
    ## 7              2021-11-17       2021-11-17            2021-11-17
    ## 8              2021-11-17       2021-11-17            2021-11-17
    ## 9              2021-11-17       2021-11-17            2021-11-17
    ## 10             2021-11-17       2021-11-17            2021-11-17
    ##    ena_last_update env_package Experiment   external_id host_subject_id
    ## 5       2021-11-17  human-oral ERX6583672 SAMEA10332186            10FN
    ## 6       2021-11-17  human-oral ERX6583673 SAMEA10332187            10SP
    ## 7       2021-11-17  human-oral ERX6583674 SAMEA10332188            11FP
    ## 8       2021-11-17  human-oral ERX6583675 SAMEA10332189            11SN
    ## 9       2021-11-17  human-oral ERX6583676 SAMEA10332190            12FN
    ## 10      2021-11-17  human-oral ERX6583677 SAMEA10332191            12SP
    ##    INSDC_center_alias                                        INSDC_center_name
    ## 5              UCSDMI University of California San Diego Microbiome Initiative
    ## 6              UCSDMI University of California San Diego Microbiome Initiative
    ## 7              UCSDMI University of California San Diego Microbiome Initiative
    ## 8              UCSDMI University of California San Diego Microbiome Initiative
    ## 9              UCSDMI University of California San Diego Microbiome Initiative
    ## 10             UCSDMI University of California San Diego Microbiome Initiative
    ##      INSDC_first_public    INSDC_last_update INSDC_status     Instrument
    ## 5  2021-11-17T16:28:35Z 2021-11-17T16:28:35Z       public Illumina MiSeq
    ## 6  2021-11-17T16:28:35Z 2021-11-17T16:28:35Z       public Illumina MiSeq
    ## 7  2021-11-17T16:28:35Z 2021-11-17T16:28:35Z       public Illumina MiSeq
    ## 8  2021-11-17T16:28:35Z 2021-11-17T16:28:35Z       public Illumina MiSeq
    ## 9  2021-11-17T16:28:35Z 2021-11-17T16:28:35Z       public Illumina MiSeq
    ## 10 2021-11-17T16:28:35Z 2021-11-17T16:28:35Z       public Illumina MiSeq
    ##    Lattitude  Library.Name LibraryLayout LibrarySelection LibrarySource
    ## 5   40.524 N 13807.SSP10FN        SINGLE              PCR   METAGENOMIC
    ## 6   40.524 N 13807.SSP10SP        SINGLE              PCR   METAGENOMIC
    ## 7   40.524 N 13807.SSP11FP        SINGLE              PCR   METAGENOMIC
    ## 8   40.524 N 13807.SSP11SN        SINGLE              PCR   METAGENOMIC
    ## 9   40.524 N 13807.SSP12FN        SINGLE              PCR   METAGENOMIC
    ## 10  40.524 N 13807.SSP12SP        SINGLE              PCR   METAGENOMIC
    ##    lognitude                Organism
    ## 5   74.469 W human saliva metagenome
    ## 6   74.469 W human saliva metagenome
    ## 7   74.469 W human saliva metagenome
    ## 8   74.469 W human saliva metagenome
    ## 9   74.469 W human saliva metagenome
    ## 10  74.469 W human saliva metagenome
    ##                                    pcr_primers..exp. physical_specimen_location
    ## 5  FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT         Rutgers University
    ## 6  FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT         Rutgers University
    ## 7  FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT         Rutgers University
    ## 8  FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT         Rutgers University
    ## 9  FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT         Rutgers University
    ## 10 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT         Rutgers University
    ##    physical_specimen_remaining Platform           primer..exp.
    ## 5                         TRUE ILLUMINA CCGGACTACHVGGGTWTCTAAT
    ## 6                         TRUE ILLUMINA CCGGACTACHVGGGTWTCTAAT
    ## 7                         TRUE ILLUMINA CCGGACTACHVGGGTWTCTAAT
    ## 8                         TRUE ILLUMINA CCGGACTACHVGGGTWTCTAAT
    ## 9                         TRUE ILLUMINA CCGGACTACHVGGGTWTCTAAT
    ## 10                        TRUE ILLUMINA CCGGACTACHVGGGTWTCTAAT
    ##             ReleaseDate   run_center..exp. run_date..exp.          create_date
    ## 5  2021-11-28T00:00:00Z Rutgers University        11/1/20 2021-12-23T04:03:00Z
    ## 6  2021-11-28T00:00:00Z Rutgers University        11/1/20 2021-12-23T04:04:00Z
    ## 7  2021-11-28T00:00:00Z Rutgers University        11/1/20 2021-12-23T04:03:00Z
    ## 8  2021-11-28T00:00:00Z Rutgers University        11/1/20 2021-12-23T04:04:00Z
    ## 9  2021-11-28T00:00:00Z Rutgers University        11/1/20 2021-12-23T04:03:00Z
    ## 10 2021-11-28T00:00:00Z Rutgers University        11/1/20 2021-12-25T18:35:00Z
    ##    version run_prefix..exp.   Sample.Name                   Sample_name
    ## 5        1           CK6254 SAMEA10332186 qiita_sid_13807:13807.SSP10FN
    ## 6        1           CK6254 SAMEA10332187 qiita_sid_13807:13807.SSP10SP
    ## 7        1           CK6254 SAMEA10332188 qiita_sid_13807:13807.SSP11FP
    ## 8        1           CK6254 SAMEA10332189 qiita_sid_13807:13807.SSP11SN
    ## 9        1           CK6254 SAMEA10332190 qiita_sid_13807:13807.SSP12FN
    ## 10       1           CK6254 SAMEA10332191 qiita_sid_13807:13807.SSP12SP
    ##    sample_type sequencing_method..exp. SRA.Study                  Submitter_Id
    ## 5       saliva Sequencing by Synthesis ERP132364 qiita_sid_13807:13807.SSP10FN
    ## 6       saliva Sequencing by Synthesis ERP132364 qiita_sid_13807:13807.SSP10SP
    ## 7       saliva Sequencing by Synthesis ERP132364 qiita_sid_13807:13807.SSP11FP
    ## 8       saliva Sequencing by Synthesis ERP132364 qiita_sid_13807:13807.SSP11SN
    ## 9       saliva Sequencing by Synthesis ERP132364 qiita_sid_13807:13807.SSP12FN
    ## 10      saliva Sequencing by Synthesis ERP132364 qiita_sid_13807:13807.SSP12SP
    ##    target_gene..exp. target_subfragment..exp.   env_biome
    ## 5           16S rRNA                       V4 Urban Biome
    ## 6           16S rRNA                       V4 Urban Biome
    ## 7           16S rRNA                       V4 Urban Biome
    ## 8           16S rRNA                       V4 Urban Biome
    ## 9           16S rRNA                       V4 Urban Biome
    ## 10          16S rRNA                       V4 Urban Biome
    ##                 env_feature env_material Host_age host_age_units
    ## 5  Human-associated habitat       saliva       70          years
    ## 6  Human-associated habitat       saliva       70          years
    ## 7  Human-associated habitat       saliva       25          years
    ## 8  Human-associated habitat       saliva       25          years
    ## 9  Human-associated habitat       saliva       27          years
    ## 10 Human-associated habitat       saliva       27          years
    ##    host_common_name host_scientific_name host_sex host_taxid Sample_Day
    ## 5             human         Homo sapiens     male       9606          1
    ## 6             human         Homo sapiens     male       9606          1
    ## 7             human         Homo sapiens   female       9606          1
    ## 8             human         Homo sapiens   female       9606          1
    ## 9             human         Homo sapiens   female       9606          1
    ## 10            human         Homo sapiens   female       9606          1
    ##    sample_order
    ## 5             1
    ## 6             2
    ## 7             1
    ## 8             2
    ## 9             1
    ## 10            2

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
16. **matched_samples** : Matched samples across sites

``` r
# Create a uniform metadata dataframe
metadata <- data.frame(
  BioProject = initial_meta$BioProject,
  Run = initial_meta$Run,
  sample_name = initial_meta$Sample.Name,
  Age = initial_meta$Host_age,
  BMI = NA,  
  Sex = initial_meta$host_sex,
  Region = initial_meta$target_subfragment..exp.,
  Platform = initial_meta$Instrument,
  Site = initial_meta$env_material,
  public_time = initial_meta$ReleaseDate,
  stringsAsFactors=FALSE
)

metadata <- metadata %>%
  mutate(Sex = recode(Sex, "female" = "F", "male" = "M"))
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$Country <- "USA"
metadata$PCR_primers <- "338F-518R"
metadata$title <- "Assessing saliva microbiome collection and processing methods"
metadata$Disease_study <- "NO"

# Assign disease status based on sample_name patterns
metadata$Status <- "HC"

metadata$matched_samples <- "NO"
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 56"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 56"

# Step 5: Visualize the Metadata

``` r
# Select columns for visualization
plot_data <- metadata %>%
  select(-BioProject, -Run, -sample_name, -title)

# Function to create simple plots for each column
create_simple_plots <- function(data) {
  plots <- list()
  
  for (col in colnames(data)) {
    # Check if the column is all NA
    if (all(is.na(data[[col]]))) {
      # Create an empty plot with text indicating all NA
      p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "this column is NA", size = 6) +
        labs(title = col) +
        theme_void()
    } else if (is.numeric(data[[col]])) {
      # Numeric variable: histogram
      p <- ggplot(data, aes(x = .data[[col]])) +
        geom_histogram(fill = "blue", alpha = 0.6, bins = 20, na.rm = TRUE) +
        labs(title = col) +
        theme_minimal()
    } else {
      # Categorical variable: bar plot (remove NA values)
      count_data <- as.data.frame(table(data[[col]], useNA = "ifany"))
      p <- ggplot(count_data, aes(x = reorder(Var1, -Freq), y = Freq)) +
        geom_col(fill = "red", alpha = 0.6) +
        labs(title = col, x = NULL, y = "Count") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 15, hjust = 1))
    }
    plots[[col]] <- p
  }
  
  return(plots)
}

# Generate and arrange plots
plots <- create_simple_plots(plot_data)
# grid.arrange(grobs = plots, ncol = 3, top = "Metadata Summary", align = "hv")
plot_grid(plots$Age, plots$BMI, plots$Sex, plots$Region, plots$Platform, 
          plots$Country, plots$Site, plots$public_time, plots$PCR_primers, 
          ncol = 3, align = "hv", axis = "lb")
```

<img src="metadata_curation_files/figure-gfm/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />
\# Step 6: Finalize and Export the Metadata

``` r
# Reorder columns to match the desired format
metadata <- metadata[, c(
  "BioProject", "Run", "sample_name", "Age", "BMI", "Sex", 
  "Region", "Platform", "Country", "Site", "public_time", 
  "PCR_primers", "title", "Disease_study", "Status", "matched_samples"
)]

# Export the uniform metadata to a CSV file
write.csv(metadata, file = "metadata.csv", quote = FALSE)
```

# Step 6: Finalize and Export the Metadata

``` r
# Reorder columns to match the desired format
metadata <- metadata[, c(
  "BioProject", "Run", "sample_name", "Age", "BMI", "Sex", 
  "Region", "Platform", "Country", "Site", "public_time", 
  "PCR_primers", "title", "Disease_study", "Status", "matched_samples"
)]

# Export the uniform metadata to a CSV file
write.csv(metadata, file = "metadata.csv", quote = FALSE)
```

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
