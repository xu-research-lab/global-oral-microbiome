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
setwd("E:/bioproject/PRJEB48623")

# Read the initial metadata files
initial_meta <- read.csv("SraRunTable.csv", header = TRUE) 

# Display the first few rows of the metadata
head(initial_meta)
```

    ##          Run Assay.Type AvgSpotLen barcode..exp.   Bases BioProject
    ## 1 ERR7247444   AMPLICON        151  ATCGTAGTGGTC 3425435 PRJEB48623
    ## 2 ERR7247445   AMPLICON        151  GCTGCCCACCTA 1274591 PRJEB48623
    ## 3 ERR7247446   AMPLICON        151  TCGAGAGTTTGC 2579533 PRJEB48623
    ## 4 ERR7247447   AMPLICON        151  CGTTCTGGTGGT 3366696 PRJEB48623
    ## 5 ERR7247448   AMPLICON        151  ACGGTACCCTAC 1859263 PRJEB48623
    ## 6 ERR7247449   AMPLICON        151  GTCACTCCGAAC 3027852 PRJEB48623
    ##       BioSample   Bytes center_name..exp.
    ## 1 SAMEA10793111 2026879              CCME
    ## 2 SAMEA10793112  759262              CCME
    ## 3 SAMEA10793113 1526630              CCME
    ## 4 SAMEA10793114 2006675              CCME
    ## 5 SAMEA10793115 1137693              CCME
    ## 6 SAMEA10793116 1824799              CCME
    ##                                          Center.Name center_project_name..exp.
    ## 1 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT                      CCME
    ## 2 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT                      CCME
    ## 3 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT                      CCME
    ## 4 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT                      CCME
    ## 5 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT                      CCME
    ## 6 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT                      CCME
    ##   collection_timestamp Consent DATASTORE.filetype DATASTORE.provider
    ## 1           2014-09-15  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 2           2014-09-15  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 3           2014-09-15  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 4           2014-09-15  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5           2014-09-15  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 6           2014-09-15  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                       DATASTORE.region Day elevation          empo_1 empo_2
    ## 1 gs.us-east1,ncbi.public,s3.us-east-1   1    1624.1 Host-associated Animal
    ## 2 gs.us-east1,ncbi.public,s3.us-east-1   1    1624.1 Host-associated Animal
    ## 3 gs.us-east1,ncbi.public,s3.us-east-1   1    1624.1 Host-associated Animal
    ## 4 gs.us-east1,ncbi.public,s3.us-east-1   1    1624.1 Host-associated Animal
    ## 5 gs.us-east1,ncbi.public,s3.us-east-1   1    1624.1 Host-associated Animal
    ## 6 gs.us-east1,ncbi.public,s3.us-east-1   1    1624.1 Host-associated Animal
    ##             empo_3 ENA.FIRST.PUBLIC..run. ENA.LAST.UPDATE..run.   env_biome
    ## 1 Animal secretion             2021-11-10            2021-11-10 urban biome
    ## 2 Animal secretion             2021-11-10            2021-11-10 urban biome
    ## 3 Animal secretion             2021-11-10            2021-11-10 urban biome
    ## 4 Animal secretion             2021-11-10            2021-11-10 urban biome
    ## 5 Animal secretion             2021-11-10            2021-11-10 urban biome
    ## 6 Animal secretion             2021-11-10            2021-11-10 urban biome
    ##                env_feature env_material env_package Experiment   external_id
    ## 1 human-associated habitat       saliva  human-oral ERX6816903 SAMEA10793111
    ## 2 human-associated habitat       saliva  human-oral ERX6816904 SAMEA10793112
    ## 3 human-associated habitat       saliva  human-oral ERX6816905 SAMEA10793113
    ## 4 human-associated habitat       saliva  human-oral ERX6816906 SAMEA10793114
    ## 5 human-associated habitat       saliva  human-oral ERX6816907 SAMEA10793115
    ## 6 human-associated habitat       saliva  human-oral ERX6816908 SAMEA10793116
    ##        Host_age host_age_units host_body_habitat host_body_mass_index
    ## 1 not collected          years     UBERON:saliva        not collected
    ## 2 not collected          years     UBERON:saliva        not collected
    ## 3 not collected          years     UBERON:saliva        not collected
    ## 4 not collected          years     UBERON:saliva        not collected
    ## 5 not collected          years     UBERON:saliva        not collected
    ## 6 not collected          years     UBERON:saliva        not collected
    ##   host_body_product host_body_site host_common_name   host_height
    ## 1     UBERON:saliva         saliva            human not collected
    ## 2     UBERON:saliva         saliva            human not collected
    ## 3     UBERON:saliva         saliva            human not collected
    ## 4     UBERON:saliva         saliva            human not collected
    ## 5     UBERON:saliva         saliva            human not collected
    ## 6     UBERON:saliva         saliva            human not collected
    ##   host_height_units host_scientific_name host_subject_id host_taxid
    ## 1                 m         Homo sapiens              MA       9606
    ## 2                 m         Homo sapiens              MA       9606
    ## 3                 m         Homo sapiens              MA       9606
    ## 4                 m         Homo sapiens              MA       9606
    ## 5                 m         Homo sapiens              MA       9606
    ## 6                 m         Homo sapiens              MA       9606
    ##     host_weight host_weight_units Hour
    ## 1 not collected                kg    7
    ## 2 not collected                kg    7
    ## 3 not collected                kg    8
    ## 4 not collected                kg    9
    ## 5 not collected                kg    9
    ## 6 not collected                kg   10
    ##                                          INSDC_center_name   INSDC_first_public
    ## 1 University of California San Diego Microbiome Initiative 2021-11-10T16:23:23Z
    ## 2 University of California San Diego Microbiome Initiative 2021-11-10T16:23:23Z
    ## 3 University of California San Diego Microbiome Initiative 2021-11-10T16:23:23Z
    ## 4 University of California San Diego Microbiome Initiative 2021-11-10T16:23:23Z
    ## 5 University of California San Diego Microbiome Initiative 2021-11-10T16:23:23Z
    ## 6 University of California San Diego Microbiome Initiative 2021-11-10T16:23:23Z
    ##      INSDC_last_update INSDC_status          Instrument LATITUDE Library.Name
    ## 1 2021-11-10T16:23:23Z       public Illumina HiSeq 4000    40.01     14062.19
    ## 2 2021-11-10T16:23:23Z       public Illumina HiSeq 4000    40.01     14062.19
    ## 3 2021-11-10T16:23:23Z       public Illumina HiSeq 4000    40.01     14062.19
    ## 4 2021-11-10T16:23:23Z       public Illumina HiSeq 4000    40.01     14062.19
    ## 5 2021-11-10T16:23:23Z       public Illumina HiSeq 4000    40.01     14062.19
    ## 6 2021-11-10T16:23:23Z       public Illumina HiSeq 4000    40.01     14062.19
    ##   LibraryLayout LibrarySelection LibrarySource longitude min
    ## 1        SINGLE              PCR   METAGENOMIC   -105.27   0
    ## 2        SINGLE              PCR   METAGENOMIC   -105.27  30
    ## 3        SINGLE              PCR   METAGENOMIC   -105.27  40
    ## 4        SINGLE              PCR   METAGENOMIC   -105.27   5
    ## 5        SINGLE              PCR   METAGENOMIC   -105.27  35
    ## 6        SINGLE              PCR   METAGENOMIC   -105.27   0
    ##                Organism                                 pcr_primers..exp.
    ## 1 human oral metagenome FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 2 human oral metagenome FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 3 human oral metagenome FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 4 human oral metagenome FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 5 human oral metagenome FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 6 human oral metagenome FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ##   physical_specimen_location Platform        primer..exp.          ReleaseDate
    ## 1                     UCSDMI ILLUMINA GTGCCAGCMGCCGCGGTAA 2021-11-22T00:00:00Z
    ## 2                     UCSDMI ILLUMINA GTGCCAGCMGCCGCGGTAA 2021-11-22T00:00:00Z
    ## 3                     UCSDMI ILLUMINA GTGCCAGCMGCCGCGGTAA 2021-11-22T00:00:00Z
    ## 4                     UCSDMI ILLUMINA GTGCCAGCMGCCGCGGTAA 2021-11-22T00:00:00Z
    ## 5                     UCSDMI ILLUMINA GTGCCAGCMGCCGCGGTAA 2021-11-22T00:00:00Z
    ## 6                     UCSDMI ILLUMINA GTGCCAGCMGCCGCGGTAA 2021-11-22T00:00:00Z
    ##   run_center..exp.          create_date version   Sample.Name
    ## 1             CCME 2021-11-16T17:20:00Z       1 SAMEA10793111
    ## 2             CCME 2021-11-16T17:20:00Z       1 SAMEA10793112
    ## 3             CCME 2021-11-16T17:20:00Z       1 SAMEA10793113
    ## 4             CCME 2021-11-16T17:20:00Z       1 SAMEA10793114
    ## 5             CCME 2021-11-16T17:20:00Z       1 SAMEA10793115
    ## 6             CCME 2021-11-16T17:21:00Z       1 SAMEA10793116
    ##                   Sample_name sample_type sequencing_method..exp. SRA.Study
    ## 1 qiita_sid_14062:14062.19122      saliva Sequencing by synthesis ERP133012
    ## 2 qiita_sid_14062:14062.19123      saliva Sequencing by synthesis ERP133012
    ## 3 qiita_sid_14062:14062.19124      saliva Sequencing by synthesis ERP133012
    ## 4 qiita_sid_14062:14062.19125      saliva Sequencing by synthesis ERP133012
    ## 5 qiita_sid_14062:14062.19126      saliva Sequencing by synthesis ERP133012
    ## 6 qiita_sid_14062:14062.19127      saliva Sequencing by synthesis ERP133012
    ##                  Submitter_Id target_gene..exp. target_subfragment..exp.
    ## 1 qiita_sid_14062:14062.19122          16S rRNA                       V4
    ## 2 qiita_sid_14062:14062.19123          16S rRNA                       V4
    ## 3 qiita_sid_14062:14062.19124          16S rRNA                       V4
    ## 4 qiita_sid_14062:14062.19125          16S rRNA                       V4
    ## 5 qiita_sid_14062:14062.19126          16S rRNA                       V4
    ## 6 qiita_sid_14062:14062.19127          16S rRNA                       V4
    ##   ENA_first_public ena_last_update INSDC_center_alias events  Other sleep
    ## 1       2021-11-10      2021-11-10             UCSDMI  sleep          yes
    ## 2       2021-11-10      2021-11-10             UCSDMI                    
    ## 3       2021-11-10      2021-11-10             UCSDMI  sleep          yes
    ## 4       2021-11-10      2021-11-10             UCSDMI                    
    ## 5       2021-11-10      2021-11-10             UCSDMI        shower      
    ## 6       2021-11-10      2021-11-10             UCSDMI          talk      
    ##   ena_first_public ena_last_update.1 Scientific_Name
    ## 1                                                   
    ## 2                                                   
    ## 3                                                   
    ## 4                                                   
    ## 5                                                   
    ## 6

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
  Age = "20-45",
  BMI = NA, 
  Sex = NA,
  Region = initial_meta$target_subfragment..exp.,
  Platform = initial_meta$Instrument,
  Site = initial_meta$env_material,
  public_time = initial_meta$ReleaseDate,
  matched_samples = NA,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$Country <- "USA"
metadata$PCR_primers <- "515F-806R"
metadata$title <- "Diurnal and eating-associated microbial patterns revealed via high-frequency saliva sampling"
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

    ## [1] "Number of unique runs: 866"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 866"

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
