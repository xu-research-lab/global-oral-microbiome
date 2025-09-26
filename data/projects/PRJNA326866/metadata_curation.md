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
library(stringr)
# Set working directory (update this path to your local directory)
setwd("/Users/dongbiao/Documents/oral/PRJNA326866/")

# Read the initial metadata files
initial_meta <- read.csv("SraRunTable.csv",header = T) %>% filter(`Sample.Name` != "Mock.DNA.even")
 
head(initial_meta)
```

    ##          Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample
    ## 1 SRR4457864   AMPLICON        501 56459952 PRJNA326866 SAMN05948754
    ## 2 SRR4457865   AMPLICON        501 61386260 PRJNA326866 SAMN05948763
    ## 3 SRR4457866   AMPLICON        501 64607576 PRJNA326866 SAMN05948764
    ## 4 SRR4457867   AMPLICON        501 56925045 PRJNA326866 SAMN05948765
    ## 5 SRR4457868   AMPLICON        501 58006225 PRJNA326866 SAMN05948766
    ## 6 SRR4457869   AMPLICON        501 51925132 PRJNA326866 SAMN05948759
    ##                        BioSampleModel    Bytes             Center.Name Consent
    ## 1 MIMS.me,MIGS/MIMS/MIMARKS.human-gut 35080616 MEDICAL UNIVERSITY GRAZ  public
    ## 2 MIMS.me,MIGS/MIMS/MIMARKS.human-gut 36722309 MEDICAL UNIVERSITY GRAZ  public
    ## 3 MIMS.me,MIGS/MIMS/MIMARKS.human-gut 38811864 MEDICAL UNIVERSITY GRAZ  public
    ## 4 MIMS.me,MIGS/MIMS/MIMARKS.human-gut 33634236 MEDICAL UNIVERSITY GRAZ  public
    ## 5 MIMS.me,MIGS/MIMS/MIMARKS.human-gut 34403093 MEDICAL UNIVERSITY GRAZ  public
    ## 6 MIMS.me,MIGS/MIMS/MIMARKS.human-gut 30660594 MEDICAL UNIVERSITY GRAZ  public
    ##   DATASTORE.filetype DATASTORE.provider                     DATASTORE.region
    ## 1   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 2   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 3   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 4   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 5   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 6   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ##   Experiment     Instrument Library.Name LibraryLayout LibrarySelection
    ## 1 SRX2273786 Illumina MiSeq         C021        PAIRED              PCR
    ## 2 SRX2273787 Illumina MiSeq         P006        PAIRED              PCR
    ## 3 SRX2273788 Illumina MiSeq         P007        PAIRED              PCR
    ## 4 SRX2273789 Illumina MiSeq         P009        PAIRED              PCR
    ## 5 SRX2273790 Illumina MiSeq         P010        PAIRED              PCR
    ## 6 SRX2273791 Illumina MiSeq         P002        PAIRED              PCR
    ##   LibrarySource         Organism Platform          ReleaseDate
    ## 1   METAGENOMIC human metagenome ILLUMINA 2016-12-12T00:00:00Z
    ## 2   METAGENOMIC human metagenome ILLUMINA 2016-12-12T00:00:00Z
    ## 3   METAGENOMIC human metagenome ILLUMINA 2016-12-12T00:00:00Z
    ## 4   METAGENOMIC human metagenome ILLUMINA 2016-12-12T00:00:00Z
    ## 5   METAGENOMIC human metagenome ILLUMINA 2016-12-12T00:00:00Z
    ## 6   METAGENOMIC human metagenome ILLUMINA 2016-12-12T00:00:00Z
    ##            create_date version    Sample.Name SRA.Study AGE
    ## 1 2016-10-27T12:44:00Z       1 C021.T0.Stool1 SRP077213    
    ## 2 2016-10-27T12:44:00Z       1 P006.T0.Stool1 SRP077213    
    ## 3 2016-10-27T12:44:00Z       1 P007.T0.Stool1 SRP077213    
    ## 4 2016-10-27T12:44:00Z       1 P009.T0.Stool1 SRP077213    
    ## 5 2016-10-27T12:45:00Z       1 P010.T0.Stool1 SRP077213    
    ## 6 2016-10-27T12:45:00Z       1 P002.T0.Stool1 SRP077213    
    ##   BIOMATERIAL_PROVIDER isolate sex tissue Collection_Date env_biome env_feature
    ## 1                                                 2014-10     Human   Intestine
    ## 2                                                 2014-08     Human   Intestine
    ## 3                                                 2014-08     Human   Intestine
    ## 4                                                 2014-09     Human   Intestine
    ## 5                                                 2014-10     Human   Intestine
    ## 6                                                 2014-06     Human   Intestine
    ##   env_material geo_loc_name_country geo_loc_name_country_continent
    ## 1        Stool              Austria                         Europe
    ## 2        Stool              Austria                         Europe
    ## 3        Stool              Austria                         Europe
    ## 4        Stool              Austria                         Europe
    ## 5        Stool              Austria                         Europe
    ## 6        Stool              Austria                         Europe
    ##    geo_loc_name Host_disease         Host host_subject_id             lat_lon
    ## 1 Austria: Graz      Control Homo sapiens            C021 47.0707 N 15.4395 E
    ## 2 Austria: Graz         PCOS Homo sapiens            P006 47.0707 N 15.4395 E
    ## 3 Austria: Graz         PCOS Homo sapiens            P007 47.0707 N 15.4395 E
    ## 4 Austria: Graz         PCOS Homo sapiens            P009 47.0707 N 15.4395 E
    ## 5 Austria: Graz         PCOS Homo sapiens            P010 47.0707 N 15.4395 E
    ## 6 Austria: Graz         PCOS Homo sapiens            P002 47.0707 N 15.4395 E

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
subject_id <- gsub(".Stool1|.Saliva1", "", initial_meta$Sample.Name)
sites <- sapply(strsplit(initial_meta$Sample.Name, "\\."), function(x) x[3])
sites <- gsub("[0-9]", "", sites)
metadata <- data.frame(
  BioProject = initial_meta$BioProject,
  Run = initial_meta$Run,
  sample_name = subject_id,
  Age = as.integer(gsub("\\\\,[0-9]|\\\\[0-9]", "", initial_meta$AGE)),
  BMI = NA,
  Sex= "F",
  Region = "V1-V2",
  Platform = initial_meta$Instrument,
  Country = "Austria: Graz",
  Site = str_to_lower(sites),
  public_time = initial_meta$ReleaseDate,
  stringsAsFactors = FALSE
)

meta_age <- metadata[complete.cases(metadata$Age), ]
my_dict <- setNames(meta_age$Age, meta_age$sample_name)
metadata$Age <- my_dict[metadata$sample_name]
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$PCR_primers <- NA
metadata$title <- "Alterations in Gut Microbiome Composition and Barrier Function Are Associated with Reproductive and Metabolic Defects in Women with Polycystic Ovary Syndrome (PCOS): A Pilot Study"
metadata$Disease_study <- "PCOS"

# Assign disease status based on sample_name patterns
status <- substr(subject_id, 1, 1)
metadata <- mutate(metadata, Status = if_else(status == "C", "H", "PCOS"))

# Create matched_samples column
metadata$matched_samples <- "YES"

# Check initial data distribution
cat("Initial sample distribution:\n")
```

    ## Initial sample distribution:

``` r
table(metadata$sample_name)
```

    ## 
    ## C001.T0 C003.T0 C004.T0 C005.T0 C007.T0 C008.T0 C009.T0 C011.T0 C012.T0 C014.T0 
    ##       2       2       2       2       2       2       2       2       2       2 
    ## C015.T0 C016.T0 C017.T0 C018.T0 C019.T0 C020.T0 C021.T0 C022.T0 C024.T0 C025.T0 
    ##       2       2       2       2       1       2       2       2       2       2 
    ## P001.T0 P002.T0 P003.T0 P004.T0 P005.T0 P006.T0 P007.T0 P009.T0 P010.T0 P011.T0 
    ##       2       2       2       2       2       2       2       2       2       2 
    ## P012.T0 P013.T0 P014.T0 P015.T0 P016.T0 P017.T0 P018.T0 P019.T0 P020.T0 P021.T0 
    ##       2       2       2       2       2       2       2       2       2       2 
    ## P022.T0 P023.T0 P024.T0 P025.T0 
    ##       2       2       2       2

``` r
table(metadata$Site)
```

    ## 
    ## saliva  stool 
    ##     44     43

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 87"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 44"

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

    ##    BioProject        Run sample_name Age BMI Sex Region       Platform
    ## 1 PRJNA326866 SRR4457864     C021.T0  26  NA   F  V1-V2 Illumina MiSeq
    ## 2 PRJNA326866 SRR4457865     P006.T0  42  NA   F  V1-V2 Illumina MiSeq
    ## 3 PRJNA326866 SRR4457866     P007.T0  35  NA   F  V1-V2 Illumina MiSeq
    ## 4 PRJNA326866 SRR4457867     P009.T0  24  NA   F  V1-V2 Illumina MiSeq
    ## 5 PRJNA326866 SRR4457868     P010.T0  36  NA   F  V1-V2 Illumina MiSeq
    ## 6 PRJNA326866 SRR4457869     P002.T0  25  NA   F  V1-V2 Illumina MiSeq
    ##         Country  Site          public_time PCR_primers
    ## 1 Austria: Graz stool 2016-12-12T00:00:00Z          NA
    ## 2 Austria: Graz stool 2016-12-12T00:00:00Z          NA
    ## 3 Austria: Graz stool 2016-12-12T00:00:00Z          NA
    ## 4 Austria: Graz stool 2016-12-12T00:00:00Z          NA
    ## 5 Austria: Graz stool 2016-12-12T00:00:00Z          NA
    ## 6 Austria: Graz stool 2016-12-12T00:00:00Z          NA
    ##                                                                                                                                                                                 title
    ## 1 Alterations in Gut Microbiome Composition and Barrier Function Are Associated with Reproductive and Metabolic Defects in Women with Polycystic Ovary Syndrome (PCOS): A Pilot Study
    ## 2 Alterations in Gut Microbiome Composition and Barrier Function Are Associated with Reproductive and Metabolic Defects in Women with Polycystic Ovary Syndrome (PCOS): A Pilot Study
    ## 3 Alterations in Gut Microbiome Composition and Barrier Function Are Associated with Reproductive and Metabolic Defects in Women with Polycystic Ovary Syndrome (PCOS): A Pilot Study
    ## 4 Alterations in Gut Microbiome Composition and Barrier Function Are Associated with Reproductive and Metabolic Defects in Women with Polycystic Ovary Syndrome (PCOS): A Pilot Study
    ## 5 Alterations in Gut Microbiome Composition and Barrier Function Are Associated with Reproductive and Metabolic Defects in Women with Polycystic Ovary Syndrome (PCOS): A Pilot Study
    ## 6 Alterations in Gut Microbiome Composition and Barrier Function Are Associated with Reproductive and Metabolic Defects in Women with Polycystic Ovary Syndrome (PCOS): A Pilot Study
    ##   Disease_study Status matched_samples
    ## 1          PCOS      H             YES
    ## 2          PCOS   PCOS             YES
    ## 3          PCOS   PCOS             YES
    ## 4          PCOS   PCOS             YES
    ## 5          PCOS   PCOS             YES
    ## 6          PCOS   PCOS             YES

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
