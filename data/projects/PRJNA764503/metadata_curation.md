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
setwd("D:/oral_project/PRJNA764503")

# Read the initial metadata files
initial_meta<-read.csv("SraRunTable.csv",header = T)

# Combine the two metadata files
#meta <- rbind(initial_meta, initial_meta1)

#Remove entries where Organism is "gut metagenome" (non-oral samples)
initial_meta <- subset(initial_meta,isolation_source!= "feces")

# Check the number of samples containing "HT" in Sample.Name (e.g., healthy controls)
#nrow(subset(meta, grepl("HT", Sample.Name)))

# Display the first few rows of the metadata
head(initial_meta)
```

    ##            Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample
    ## 10 SRR15974357   AMPLICON        495 57234375 PRJNA764503 SAMN21504824
    ## 21 SRR15974368   AMPLICON        492 60240972 PRJNA764503 SAMN21504823
    ## 32 SRR15974379   AMPLICON        496 17064662 PRJNA764503 SAMN21504822
    ## 43 SRR15974390   AMPLICON        496 24734571 PRJNA764503 SAMN21504821
    ## 44 SRR15974391   AMPLICON        494 18393982 PRJNA764503 SAMN21504785
    ## 45 SRR15974392   AMPLICON        495 52319520 PRJNA764503 SAMN21505180
    ##                 BioSampleModel    Bytes
    ## 10 Metagenome or environmental 18245781
    ## 21 Metagenome or environmental 19221909
    ## 32 Metagenome or environmental  8743975
    ## 43 Metagenome or environmental  8757843
    ## 44 Metagenome or environmental  9581297
    ## 45 Metagenome or environmental 16595212
    ##                                         Center.Name      Collection_Date
    ## 10 SHANGHAI JIAO TONG UNIVERSITY SCHOOL OF MEDICINE 2020-10-11T17:53:46Z
    ## 21 SHANGHAI JIAO TONG UNIVERSITY SCHOOL OF MEDICINE 2020-10-11T17:53:45Z
    ## 32 SHANGHAI JIAO TONG UNIVERSITY SCHOOL OF MEDICINE 2020-10-11T17:53:44Z
    ## 43 SHANGHAI JIAO TONG UNIVERSITY SCHOOL OF MEDICINE 2020-10-11T17:53:43Z
    ## 44 SHANGHAI JIAO TONG UNIVERSITY SCHOOL OF MEDICINE 2020-10-11T17:53:07Z
    ## 45 SHANGHAI JIAO TONG UNIVERSITY SCHOOL OF MEDICINE 2020-10-11T13:46:47Z
    ##    Consent DATASTORE.filetype DATASTORE.provider
    ## 10  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 21  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 32  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 43  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 44  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 45  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                        DATASTORE.region  Experiment geo_loc_name_country
    ## 10 gs.us-east1,ncbi.public,s3.us-east-1 SRX12263460                China
    ## 21 gs.us-east1,ncbi.public,s3.us-east-1 SRX12263449                China
    ## 32 gs.us-east1,ncbi.public,s3.us-east-1 SRX12263438                China
    ## 43 gs.us-east1,ncbi.public,s3.us-east-1 SRX12263427                China
    ## 44 gs.us-east1,ncbi.public,s3.us-east-1 SRX12263426                China
    ## 45 gs.us-east1,ncbi.public,s3.us-east-1 SRX12263425                China
    ##    geo_loc_name_country_continent   geo_loc_name         HOST
    ## 10                           Asia China:Shanghai Homo sapiens
    ## 21                           Asia China:Shanghai Homo sapiens
    ## 32                           Asia China:Shanghai Homo sapiens
    ## 43                           Asia China:Shanghai Homo sapiens
    ## 44                           Asia China:Shanghai Homo sapiens
    ## 45                           Asia China:Shanghai Homo sapiens
    ##               Instrument   isolation_source       lat_lon Library.Name
    ## 10 Illumina NovaSeq 6000             saliva not collected    1907127S3
    ## 21 Illumina NovaSeq 6000             saliva not collected    1907127S1
    ## 32 Illumina NovaSeq 6000             saliva not collected     190713S2
    ## 43 Illumina NovaSeq 6000             saliva not collected     190706S1
    ## 44 Illumina NovaSeq 6000             saliva not collected     181110S1
    ## 45 Illumina NovaSeq 6000 subgingival palque not collected     X00516P1
    ##    LibraryLayout LibrarySelection LibrarySource              Organism Platform
    ## 10        PAIRED              PCR       GENOMIC human oral metagenome ILLUMINA
    ## 21        PAIRED              PCR       GENOMIC human oral metagenome ILLUMINA
    ## 32        PAIRED              PCR       GENOMIC human oral metagenome ILLUMINA
    ## 43        PAIRED              PCR       GENOMIC human oral metagenome ILLUMINA
    ## 44        PAIRED              PCR       GENOMIC human oral metagenome ILLUMINA
    ## 45        PAIRED              PCR       GENOMIC human oral metagenome ILLUMINA
    ##             ReleaseDate          create_date version Sample.Name SRA.Study
    ## 10 2022-03-20T00:00:00Z 2021-09-20T14:29:00Z       1   1907127S3 SRP337928
    ## 21 2022-03-20T00:00:00Z 2021-09-20T14:29:00Z       1   1907127S1 SRP337928
    ## 32 2022-03-20T00:00:00Z 2021-09-20T14:29:00Z       1    190713S2 SRP337928
    ## 43 2022-03-20T00:00:00Z 2021-09-20T14:30:00Z       1    190706S1 SRP337928
    ## 44 2022-03-20T00:00:00Z 2021-09-20T14:30:00Z       1    181110S1 SRP337928
    ## 45 2022-03-20T00:00:00Z 2021-09-20T14:29:00Z       1    X00516P1 SRP337928

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
metadata <- uniform_meta <- data.frame(
  BioProject = initial_meta$BioProject,
  Run = initial_meta$Run,
  sample_name = initial_meta$Sample.Name,
  Age = NA,
  BMI = NA,
  Sex=NA,
  Region = "V3-V4",
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  Site = initial_meta$isolation_source,
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
metadata$title <- "Integrated Omic Analysis of Human Plasma Metabolites and Microbiota in a Hypertension Cohort"
#metadata<-subset(uniform_meta,Site!="feces")
metadata$Disease_study <- "YES"
# Assign disease status based on sample_name patterns
metadata$Status <- NA
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 400"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 400"

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
