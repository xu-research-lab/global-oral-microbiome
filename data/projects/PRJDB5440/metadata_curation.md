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
library(gridExtra)

# Set working directory (update this path to your local directory)
setwd("D:/oral_project/PRJDB5440")

# Read the initial metadata files
initial_meta<-read.csv("SraRunTable.csv",header = T)

# Combine the two metadata files
#meta <- rbind(initial_meta, initial_meta1)

# Remove entries where Organism is "gut metagenome" (non-oral samples)
#meta <- subset(meta, Organism != "gut metagenome")

# Check the number of samples containing "HT" in Sample.Name (e.g., healthy controls)
#nrow(subset(meta, grepl("HT", Sample.Name)))

# Display the first few rows of the metadata
head(initial_meta)
```

    ##         Run Assay.Type AvgSpotLen    Bases BioProject    BioSample
    ## 1 DRR082910   AMPLICON        602 13702122  PRJDB5440 SAMD00070851
    ## 2 DRR082911   AMPLICON        602 13623260  PRJDB5440 SAMD00070852
    ## 3 DRR082912   AMPLICON        602 13467944  PRJDB5440 SAMD00070853
    ## 4 DRR082913   AMPLICON        602 13281324  PRJDB5440 SAMD00070854
    ## 5 DRR082914   AMPLICON        602 14331212  PRJDB5440 SAMD00070855
    ## 6 DRR082915   AMPLICON        602 14205394  PRJDB5440 SAMD00070856
    ##       BioSampleModel   Bytes Center.Name Collection_Date Consent
    ## 1 MIMS.me.human-oral 7489539          UB      2015-03-23  public
    ## 2 MIMS.me.human-oral 7341031          UB      2015-03-23  public
    ## 3 MIMS.me.human-oral 7369055          UB      2015-03-23  public
    ## 4 MIMS.me.human-oral 7161874          UB      2015-03-23  public
    ## 5 MIMS.me.human-oral 7657480          UB      2015-03-23  public
    ## 6 MIMS.me.human-oral 7619319          UB      2015-03-23  public
    ##   DATASTORE.filetype DATASTORE.provider                     DATASTORE.region
    ## 1         run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 2         run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 3         run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 4         run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 5         run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 6         run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ##   env_broad_scale env_local_scale   env_medium Experiment geo_loc_name_country
    ## 1  not applicable  not applicable human saliva  DRX076741             Thailand
    ## 2  not applicable  not applicable human saliva  DRX076742             Thailand
    ## 3  not applicable  not applicable human saliva  DRX076743             Thailand
    ## 4  not applicable  not applicable human saliva  DRX076744             Thailand
    ## 5  not applicable  not applicable human saliva  DRX076745             Thailand
    ## 6  not applicable  not applicable human saliva  DRX076746             Thailand
    ##   geo_loc_name_country_continent                           geo_loc_name
    ## 1                           Asia Thailand:Bangkok\\, Mahidol university
    ## 2                           Asia Thailand:Bangkok\\, Mahidol university
    ## 3                           Asia Thailand:Bangkok\\, Mahidol university
    ## 4                           Asia Thailand:Bangkok\\, Mahidol university
    ## 5                           Asia Thailand:Bangkok\\, Mahidol university
    ## 6                           Asia Thailand:Bangkok\\, Mahidol university
    ##           Host     Instrument        lat_lon Library.Name LibraryLayout
    ## 1 Homo sapiens Illumina MiSeq not applicable           NA        PAIRED
    ## 2 Homo sapiens Illumina MiSeq not applicable           NA        PAIRED
    ## 3 Homo sapiens Illumina MiSeq not applicable           NA        PAIRED
    ## 4 Homo sapiens Illumina MiSeq not applicable           NA        PAIRED
    ## 5 Homo sapiens Illumina MiSeq not applicable           NA        PAIRED
    ## 6 Homo sapiens Illumina MiSeq not applicable           NA        PAIRED
    ##   LibrarySelection LibrarySource                Organism Platform
    ## 1              PCR       GENOMIC human saliva metagenome ILLUMINA
    ## 2              PCR       GENOMIC human saliva metagenome ILLUMINA
    ## 3              PCR       GENOMIC human saliva metagenome ILLUMINA
    ## 4              PCR       GENOMIC human saliva metagenome ILLUMINA
    ## 5              PCR       GENOMIC human saliva metagenome ILLUMINA
    ## 6              PCR       GENOMIC human saliva metagenome ILLUMINA
    ##                      project_name          ReleaseDate  Sample.Name Sample_name
    ## 1 Exploring Human Oral Microbiome 2017-09-22T00:00:00Z SAMD00070851       L3010
    ## 2 Exploring Human Oral Microbiome 2017-09-22T00:00:00Z SAMD00070852       L3011
    ## 3 Exploring Human Oral Microbiome 2017-09-22T00:00:00Z SAMD00070853       L3012
    ## 4 Exploring Human Oral Microbiome 2017-09-22T00:00:00Z SAMD00070854       L3013
    ## 5 Exploring Human Oral Microbiome 2017-09-22T00:00:00Z SAMD00070855       L3014
    ## 6 Exploring Human Oral Microbiome 2017-09-22T00:00:00Z SAMD00070856       L3015
    ##   SRA.Study          create_date version
    ## 1 DRP003829 2017-09-22T16:06:00Z       1
    ## 2 DRP003829 2017-09-22T16:08:00Z       1
    ## 3 DRP003829 2017-09-22T16:08:00Z       1
    ## 4 DRP003829 2017-09-22T16:05:00Z       1
    ## 5 DRP003829 2017-09-22T16:07:00Z       1
    ## 6 DRP003829 2017-09-22T16:06:00Z       1

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
  sample_name = initial_meta$Sample_name,
  Age = "children",
  BMI = NA,
  Sex=NA,
  Region = "V3-V4",
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  Site = initial_meta$env_medium,
  public_time = initial_meta$ReleaseDate,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$PCR_primers <- NA
metadata$title <- "Exploring the salivary microbiome of children stratified by the oral hygiene index"
metadata$Disease_study <- "YES"

# Assign disease status based on sample_name patterns
metadata$Status <- "OSCC"
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 140"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 140"

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
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
    plots[[col]] <- p
  }
  
  return(plots)
}

# Generate and arrange plots
plots <- create_simple_plots(plot_data)
grid.arrange(grobs = plots, ncol = 3, top = "Metadata Summary")
```

![](Metadata_Summary.png)<!-- -->

# Step 6: Finalize and Export the Metadata

``` r
# Reorder columns to match the desired format
metadata <- metadata[, c(
  "BioProject", "Run", "sample_name", "Age", "BMI", "Sex", 
  "Region", "Platform", "Country", "Site", "public_time", 
  "PCR_primers", "title", "Disease_study", "Status"
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
