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
setwd("D:/oral_project/PRJNA866676")

# Read the initial metadata files
initial_meta<-read.csv("SraRunTable.csv",header = T)

# Combine the two metadata files
#meta <- rbind(initial_meta, initial_meta1)

#Remove entries non-oral samples
initial_meta <- subset(initial_meta,isolate %in% c("saliva","tongue coat"))

# Check the number of samples containing "HT" in Sample.Name (e.g., healthy controls)
#nrow(subset(meta, grepl("HT", Sample.Name)))

# Display the first few rows of the metadata
head(initial_meta)
```

    ##           Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample
    ## 1 SRR20888871   AMPLICON        411 42243572 PRJNA866676 SAMN30170762
    ## 2 SRR20888872   AMPLICON        416 45077663 PRJNA866676 SAMN30170761
    ## 3 SRR20888873   AMPLICON        412 42087917 PRJNA866676 SAMN30170760
    ## 4 SRR20888874   AMPLICON        413 37886669 PRJNA866676 SAMN30170759
    ## 5 SRR20888875   AMPLICON        417 46423989 PRJNA866676 SAMN30170758
    ## 6 SRR20888876   AMPLICON        417 46180323 PRJNA866676 SAMN30170757
    ##                                  BioSampleModel    Bytes
    ## 1 MIMARKS.specimen,MIGS/MIMS/MIMARKS.human-oral 34114909
    ## 2 MIMARKS.specimen,MIGS/MIMS/MIMARKS.human-oral 36652557
    ## 3 MIMARKS.specimen,MIGS/MIMS/MIMARKS.human-oral 33828766
    ## 4 MIMARKS.specimen,MIGS/MIMS/MIMARKS.human-oral 30801457
    ## 5 MIMARKS.specimen,MIGS/MIMS/MIMARKS.human-oral 37815434
    ## 6 MIMARKS.specimen,MIGS/MIMS/MIMARKS.human-oral 36865956
    ##                                                                             Center.Name
    ## 1 SCHOOL AND HOSPITAL OF STOMATOLOGY, CHEELOO COLLEGE OF  MEDICINE, SHANDONG UNIVERSITY
    ## 2 SCHOOL AND HOSPITAL OF STOMATOLOGY, CHEELOO COLLEGE OF  MEDICINE, SHANDONG UNIVERSITY
    ## 3 SCHOOL AND HOSPITAL OF STOMATOLOGY, CHEELOO COLLEGE OF  MEDICINE, SHANDONG UNIVERSITY
    ## 4 SCHOOL AND HOSPITAL OF STOMATOLOGY, CHEELOO COLLEGE OF  MEDICINE, SHANDONG UNIVERSITY
    ## 5 SCHOOL AND HOSPITAL OF STOMATOLOGY, CHEELOO COLLEGE OF  MEDICINE, SHANDONG UNIVERSITY
    ## 6 SCHOOL AND HOSPITAL OF STOMATOLOGY, CHEELOO COLLEGE OF  MEDICINE, SHANDONG UNIVERSITY
    ##        Collection_Date Consent DATASTORE.filetype DATASTORE.provider
    ## 1 2018-10-10T08:37:00Z  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 2 2018-10-10T08:24:00Z  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 3 2018-10-10T08:13:00Z  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 4 2018-10-10T05:20:00Z  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5 2018-09-20T07:29:00Z  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 6 2018-09-09T03:23:00Z  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                       DATASTORE.region                        env_broad_scale
    ## 1 gs.us-east1,ncbi.public,s3.us-east-1 mucosa of oral region [UBERON:0003446]
    ## 2 gs.us-east1,ncbi.public,s3.us-east-1 mucosa of oral region [UBERON:0003445]
    ## 3 gs.us-east1,ncbi.public,s3.us-east-1 mucosa of oral region [UBERON:0003444]
    ## 4 gs.us-east1,ncbi.public,s3.us-east-1 mucosa of oral region [UBERON:0003443]
    ## 5 gs.us-east1,ncbi.public,s3.us-east-1 mucosa of oral region [UBERON:0003442]
    ## 6 gs.us-east1,ncbi.public,s3.us-east-1 mucosa of oral region [UBERON:0003441]
    ##                          env_local_scale                             env_medium
    ## 1 mucosa of oral region [UBERON:0003446] mucosa of oral region [UBERON:0003446]
    ## 2 mucosa of oral region [UBERON:0003445] mucosa of oral region [UBERON:0003445]
    ## 3 mucosa of oral region [UBERON:0003444] mucosa of oral region [UBERON:0003444]
    ## 4 mucosa of oral region [UBERON:0003443] mucosa of oral region [UBERON:0003443]
    ## 5 mucosa of oral region [UBERON:0003442] mucosa of oral region [UBERON:0003442]
    ## 6 mucosa of oral region [UBERON:0003441] mucosa of oral region [UBERON:0003441]
    ##    Experiment geo_loc_name_country geo_loc_name_country_continent geo_loc_name
    ## 1 SRX16907766                China                           Asia China: Jinan
    ## 2 SRX16907765                China                           Asia China: Jinan
    ## 3 SRX16907764                China                           Asia China: Jinan
    ## 4 SRX16907763                China                           Asia China: Jinan
    ## 5 SRX16907762                China                           Asia China: Jinan
    ## 6 SRX16907761                China                           Asia China: Jinan
    ##           Host        Instrument                          isol_growth_condt
    ## 1 Homo sapiens Ion Torrent S5 XL https://doi.org/10.1016/j.sjbs.2021.08.272
    ## 2 Homo sapiens Ion Torrent S5 XL https://doi.org/10.1016/j.sjbs.2021.08.271
    ## 3 Homo sapiens Ion Torrent S5 XL https://doi.org/10.1016/j.sjbs.2021.08.270
    ## 4 Homo sapiens Ion Torrent S5 XL https://doi.org/10.1016/j.sjbs.2021.08.269
    ## 5 Homo sapiens Ion Torrent S5 XL https://doi.org/10.1016/j.sjbs.2021.08.268
    ## 6 Homo sapiens Ion Torrent S5 XL https://doi.org/10.1016/j.sjbs.2021.08.267
    ##       isolate      lat_lon Library.Name LibraryLayout LibrarySelection
    ## 1 tongue coat 36.4 N 117 E          177        SINGLE              PCR
    ## 2 tongue coat 36.4 N 117 E          176        SINGLE              PCR
    ## 3 tongue coat 36.4 N 117 E          175        SINGLE              PCR
    ## 4 tongue coat 36.4 N 117 E          174        SINGLE              PCR
    ## 5 tongue coat 36.4 N 117 E          173        SINGLE              PCR
    ## 6 tongue coat 36.4 N 117 E          172        SINGLE              PCR
    ##   LibrarySource              Organism    Platform          ReleaseDate
    ## 1       GENOMIC human oral metagenome ION_TORRENT 2022-08-06T00:00:00Z
    ## 2       GENOMIC human oral metagenome ION_TORRENT 2022-08-06T00:00:00Z
    ## 3       GENOMIC human oral metagenome ION_TORRENT 2022-08-06T00:00:00Z
    ## 4       GENOMIC human oral metagenome ION_TORRENT 2022-08-06T00:00:00Z
    ## 5       GENOMIC human oral metagenome ION_TORRENT 2022-08-06T00:00:00Z
    ## 6       GENOMIC human oral metagenome ION_TORRENT 2022-08-06T00:00:00Z
    ##            create_date version Sample.Name SRA.Study
    ## 1 2022-08-06T06:06:00Z       1         T64 SRP390217
    ## 2 2022-08-06T06:06:00Z       1         T63 SRP390217
    ## 3 2022-08-06T06:06:00Z       1         T62 SRP390217
    ## 4 2022-08-06T06:06:00Z       1         T61 SRP390217
    ## 5 2022-08-06T06:06:00Z       1         T60 SRP390217
    ## 6 2022-08-06T06:06:00Z       1         T59 SRP390217

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
  Region = "V3-V4",
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
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
metadata$title <- "Characteristics of Microbial Distribution in Different Oral Niches of Oral Squamous Cell Carcinoma"
metadata$Disease_study <- "YES"
# Assign disease status based on sample_name patterns
metadata$Status<-"OSCC"
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 133"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 133"

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
