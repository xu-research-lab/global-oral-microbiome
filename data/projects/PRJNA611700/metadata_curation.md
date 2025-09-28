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
setwd("/Users/dongbiao/Documents/oral/global-oral-microbiome/data/projects/PRJNA611700/")

# Read the initial metadata files
initial_meta<-read.csv("SraRunTable.csv",header = T)

head(initial_meta)
```

    ##           Run     AGE Assay.Type AvgSpotLen    Bases
    ## 1 SRR11279612 missing   AMPLICON        500 20896500
    ## 2 SRR11279613 missing   AMPLICON        499 28452916
    ## 3 SRR11279614 missing   AMPLICON        499 25701954
    ## 4 SRR11279615 missing   AMPLICON        499 39914736
    ## 5 SRR11279616 missing   AMPLICON        499 26555166
    ## 6 SRR11279617 missing   AMPLICON        499 39873212
    ##                                              BIOMATERIAL_PROVIDER  BioProject
    ## 1 School of Biological Science and Medical Engineering\\, BUAA-45 PRJNA611700
    ## 2 School of Biological Science and Medical Engineering\\, BUAA-44 PRJNA611700
    ## 3 School of Biological Science and Medical Engineering\\, BUAA-43 PRJNA611700
    ## 4 School of Biological Science and Medical Engineering\\, BUAA-42 PRJNA611700
    ## 5 School of Biological Science and Medical Engineering\\, BUAA-41 PRJNA611700
    ## 6  School of Biological Science and Medical Engineering\\, BUAA-5 PRJNA611700
    ##      BioSample BioSampleModel    Bytes
    ## 1 SAMN14342524          Human 11653816
    ## 2 SAMN14342523          Human 15773332
    ## 3 SAMN14342522          Human 14208667
    ## 4 SAMN14342521          Human 23017116
    ## 5 SAMN14342520          Human 14513541
    ## 6 SAMN14342484          Human 21964054
    ##                                                  Center.Name Consent
    ## 1 SCHOOL OF BIOLOGICAL SCIENCE AND MEDICAL ENGINEERING, BUAA  public
    ## 2 SCHOOL OF BIOLOGICAL SCIENCE AND MEDICAL ENGINEERING, BUAA  public
    ## 3 SCHOOL OF BIOLOGICAL SCIENCE AND MEDICAL ENGINEERING, BUAA  public
    ## 4 SCHOOL OF BIOLOGICAL SCIENCE AND MEDICAL ENGINEERING, BUAA  public
    ## 5 SCHOOL OF BIOLOGICAL SCIENCE AND MEDICAL ENGINEERING, BUAA  public
    ## 6 SCHOOL OF BIOLOGICAL SCIENCE AND MEDICAL ENGINEERING, BUAA  public
    ##   DATASTORE.filetype DATASTORE.provider                     DATASTORE.region
    ## 1   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 2   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 3   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 4   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 5   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 6   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ##   Experiment          Instrument isolate Library.Name LibraryLayout
    ## 1 SRX7885923 Illumina HiSeq 2500  saliva         SC12        PAIRED
    ## 2 SRX7885922 Illumina HiSeq 2500  saliva         SC11        PAIRED
    ## 3 SRX7885921 Illumina HiSeq 2500  saliva         SC09        PAIRED
    ## 4 SRX7885920 Illumina HiSeq 2500  saliva         SC08        PAIRED
    ## 5 SRX7885919 Illumina HiSeq 2500  saliva         SC07        PAIRED
    ## 6 SRX7885918 Illumina HiSeq 2500  saliva         SA06        PAIRED
    ##   LibrarySelection LibrarySource     Organism Platform          ReleaseDate
    ## 1              PCR   METAGENOMIC Homo sapiens ILLUMINA 2020-03-10T00:00:00Z
    ## 2              PCR   METAGENOMIC Homo sapiens ILLUMINA 2020-03-10T00:00:00Z
    ## 3              PCR   METAGENOMIC Homo sapiens ILLUMINA 2020-03-10T00:00:00Z
    ## 4              PCR   METAGENOMIC Homo sapiens ILLUMINA 2020-03-10T00:00:00Z
    ## 5              PCR   METAGENOMIC Homo sapiens ILLUMINA 2020-03-10T00:00:00Z
    ## 6              PCR   METAGENOMIC Homo sapiens ILLUMINA 2020-03-10T00:00:00Z
    ##            create_date version Sample.Name      sex SRA.Study   tissue
    ## 1 2020-03-10T08:51:00Z       1        SC12 intersex SRP252174 saliva45
    ## 2 2020-03-10T08:51:00Z       1        SC11 intersex SRP252174 saliva44
    ## 3 2020-03-10T08:51:00Z       1        SC09 intersex SRP252174 saliva43
    ## 4 2020-03-10T08:52:00Z       1        SC08 intersex SRP252174 saliva42
    ## 5 2020-03-10T08:52:00Z       1        SC07 intersex SRP252174 saliva41
    ## 6 2020-03-10T08:52:00Z       1        SA06 intersex SRP252174  saliva5

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
16. **matched_samples** : matched samples (“YES” or “NO”)

``` r
# Create a uniform metadata dataframe
metadata <- data.frame(
  BioProject = initial_meta$BioProject,
  Run = initial_meta$Run,
  sample_name = initial_meta$`Sample.Name`,
  Age = NA,
  BMI = NA,
  Sex=NA,
  Region = "V3-V4",
  Platform = initial_meta$Instrument,
  Country = "China",
  Site = initial_meta$isolate,
  public_time = initial_meta$ReleaseDate,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$PCR_primers <- "341F-806R"
metadata$title <- "Effects of Long-Term Enclosed Environment on Human Health Based on the Analysis of Salivary Microbiota and Cytokines"
metadata$Disease_study <- "NO"

# Assign disease status based on sample_name patterns
metadata$Status <- "H"

# Create matched_samples column
metadata$matched_samples <- "NO"

# Check initial data distribution
cat("Initial sample distribution:\n")
```

    ## Initial sample distribution:

``` r
table(metadata$sample_name)
```

    ## 
    ## SA01 SA03 SA04 SA05 SA06 SA07 SA09 SA11 SA12 SA13 SA14 SA15 SA16 SA17 SA18 SA19 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## SB01 SB02 SB03 SB04 SB05 SB06 SB07 SB08 SB09 SB10 SB11 SB12 SB13 SB15 SB16 SB17 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## SB18 SB19 SC01 SC02 SC03 SC04 SC05 SC06 SC07 SC08 SC09 SC11 SC12 SC13 SC14 SC15 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## SC17 SC18 SC19 SD01 SD02 SD03 SD04 SD05 SD06 SD07 SD08 SD09 SD10 SD11 SD12 SD13 
    ##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
    ## SD14 SD15 SD16 SD17 SD18 SD19 
    ##    1    1    1    1    1    1

``` r
table(metadata$Site)
```

    ## 
    ## saliva 
    ##     70

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 70"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 70"

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
p<-plot_grid(plots$Age, 
          plots$BMI, 
          plots$Sex, 
          plots$Region, 
          plots$Platform, 
          plots$Country, 
          plots$Site,
          plots$public_time, 
          plots$PCR_primers,
          plots$Disease_study,
          plots$Status,
          plots$matched_samples,# add 
          ncol = 3, align = "hv", axis = "lb")
ggsave("Metadata_Summary.png", plot = p, width = 8, height = 6, dpi = 150)
```

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

    ##    BioProject         Run sample_name Age BMI Sex Region            Platform
    ## 1 PRJNA611700 SRR11279612        SC12  NA  NA  NA  V3-V4 Illumina HiSeq 2500
    ## 2 PRJNA611700 SRR11279613        SC11  NA  NA  NA  V3-V4 Illumina HiSeq 2500
    ## 3 PRJNA611700 SRR11279614        SC09  NA  NA  NA  V3-V4 Illumina HiSeq 2500
    ## 4 PRJNA611700 SRR11279615        SC08  NA  NA  NA  V3-V4 Illumina HiSeq 2500
    ## 5 PRJNA611700 SRR11279616        SC07  NA  NA  NA  V3-V4 Illumina HiSeq 2500
    ## 6 PRJNA611700 SRR11279617        SA06  NA  NA  NA  V3-V4 Illumina HiSeq 2500
    ##   Country   Site          public_time PCR_primers
    ## 1   China saliva 2020-03-10T00:00:00Z   341F-806R
    ## 2   China saliva 2020-03-10T00:00:00Z   341F-806R
    ## 3   China saliva 2020-03-10T00:00:00Z   341F-806R
    ## 4   China saliva 2020-03-10T00:00:00Z   341F-806R
    ## 5   China saliva 2020-03-10T00:00:00Z   341F-806R
    ## 6   China saliva 2020-03-10T00:00:00Z   341F-806R
    ##                                                                                                                  title
    ## 1 Effects of Long-Term Enclosed Environment on Human Health Based on the Analysis of Salivary Microbiota and Cytokines
    ## 2 Effects of Long-Term Enclosed Environment on Human Health Based on the Analysis of Salivary Microbiota and Cytokines
    ## 3 Effects of Long-Term Enclosed Environment on Human Health Based on the Analysis of Salivary Microbiota and Cytokines
    ## 4 Effects of Long-Term Enclosed Environment on Human Health Based on the Analysis of Salivary Microbiota and Cytokines
    ## 5 Effects of Long-Term Enclosed Environment on Human Health Based on the Analysis of Salivary Microbiota and Cytokines
    ## 6 Effects of Long-Term Enclosed Environment on Human Health Based on the Analysis of Salivary Microbiota and Cytokines
    ##   Disease_study Status matched_samples
    ## 1            NO      H              NO
    ## 2            NO      H              NO
    ## 3            NO      H              NO
    ## 4            NO      H              NO
    ## 5            NO      H              NO
    ## 6            NO      H              NO

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
