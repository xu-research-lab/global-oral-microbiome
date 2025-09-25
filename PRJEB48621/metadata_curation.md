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
setwd("E:/bioproject/PRJEB48621")

# Read the initial metadata files
initial_meta <- read.csv("SraRunTable.csv" , header = TRUE)

# Display the first few rows of the metadata
head(initial_meta)
```

    ##          Run Assay.Type AvgSpotLen barcode..exp.   Bases BioProject
    ## 1 ERR7247129   AMPLICON        151  TCGTAAGCCGTC 1058812 PRJEB48621
    ## 2 ERR7247130   AMPLICON        151  AGCACTTTGAGA 2240689 PRJEB48621
    ## 3 ERR7247131   AMPLICON        151  TTCTCGGTTCTC 3351898 PRJEB48621
    ## 4 ERR7247132   AMPLICON        151  TCTGGGCATTGA 1951222 PRJEB48621
    ## 5 ERR7247133   AMPLICON        151  TGGAGAGGAGAT 2843179 PRJEB48621
    ## 6 ERR7247134   AMPLICON        151  GTTATGACGGAT 2322380 PRJEB48621
    ##       BioSample   Bytes center_name..exp.
    ## 1 SAMEA10792800  782299              CCME
    ## 2 SAMEA10792801 1252096              CCME
    ## 3 SAMEA10792802 2171434              CCME
    ## 4 SAMEA10792803 1060791              CCME
    ## 5 SAMEA10792804 1566934              CCME
    ## 6 SAMEA10792805 1292410              CCME
    ##                                          Center.Name center_project_name..exp.
    ## 1 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT                      CCME
    ## 2 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT                      CCME
    ## 3 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT                      CCME
    ## 4 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT                      CCME
    ## 5 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT                      CCME
    ## 6 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT                      CCME
    ##   collection_timestamp Consent DATASTORE.filetype DATASTORE.provider
    ## 1           2014-09-21  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 2           2014-09-21  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 3           2014-09-21  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 4           2014-09-21  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5           2014-09-21  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 6           2014-09-21  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                       DATASTORE.region elapsed_time elevation          empo_1
    ## 1 gs.us-east1,ncbi.public,s3.us-east-1          -46    1624.1 Host-associated
    ## 2 gs.us-east1,ncbi.public,s3.us-east-1          -36    1624.1 Host-associated
    ## 3 gs.us-east1,ncbi.public,s3.us-east-1          -24    1624.1 Host-associated
    ## 4 gs.us-east1,ncbi.public,s3.us-east-1           25    1624.1 Host-associated
    ## 5 gs.us-east1,ncbi.public,s3.us-east-1           37    1624.1 Host-associated
    ## 6 gs.us-east1,ncbi.public,s3.us-east-1           47    1624.1 Host-associated
    ##   empo_2           empo_3 ENA.FIRST.PUBLIC..run. ENA_first_public
    ## 1 Animal Animal secretion             2021-11-10       2021-11-10
    ## 2 Animal Animal secretion             2021-11-10       2021-11-10
    ## 3 Animal Animal secretion             2021-11-10       2021-11-10
    ## 4 Animal Animal secretion             2021-11-10       2021-11-10
    ## 5 Animal Animal secretion             2021-11-10       2021-11-10
    ## 6 Animal Animal secretion             2021-11-10       2021-11-10
    ##   ENA.LAST.UPDATE..run. ena_last_update   env_biome              env_feature
    ## 1            2021-11-10      2021-11-10 urban biome human-associated habitat
    ## 2            2021-11-10      2021-11-10 urban biome human-associated habitat
    ## 3            2021-11-10      2021-11-10 urban biome human-associated habitat
    ## 4            2021-11-10      2021-11-10 urban biome human-associated habitat
    ## 5            2021-11-10      2021-11-10 urban biome human-associated habitat
    ## 6            2021-11-10      2021-11-10 urban biome human-associated habitat
    ##   env_material env_package Experiment   external_id gender Host_age
    ## 1       saliva  human-oral ERX6816588 SAMEA10792800   male       27
    ## 2       saliva  human-oral ERX6816589 SAMEA10792801   male       27
    ## 3       saliva  human-oral ERX6816590 SAMEA10792802   male       27
    ## 4       saliva  human-oral ERX6816591 SAMEA10792803   male       27
    ## 5       saliva  human-oral ERX6816592 SAMEA10792804   male       27
    ## 6       saliva  human-oral ERX6816593 SAMEA10792805   male       27
    ##   host_age_units host_body_habitat host_body_mass_index host_body_product
    ## 1          years     UBERON:saliva        not collected     UBERON:saliva
    ## 2          years     UBERON:saliva        not collected     UBERON:saliva
    ## 3          years     UBERON:saliva        not collected     UBERON:saliva
    ## 4          years     UBERON:saliva        not collected     UBERON:saliva
    ## 5          years     UBERON:saliva        not collected     UBERON:saliva
    ## 6          years     UBERON:saliva        not collected     UBERON:saliva
    ##   host_body_site host_common_name   host_height host_height_units
    ## 1         saliva            human not collected                 m
    ## 2         saliva            human not collected                 m
    ## 3         saliva            human not collected                 m
    ## 4         saliva            human not collected                 m
    ## 5         saliva            human not collected                 m
    ## 6         saliva            human not collected                 m
    ##   host_scientific_name host_subject_id host_taxid host_weight host_weight_units
    ## 1         Homo sapiens              10       9606         115             catty
    ## 2         Homo sapiens              10       9606         115             catty
    ## 3         Homo sapiens              10       9606         115             catty
    ## 4         Homo sapiens              10       9606         115             catty
    ## 5         Homo sapiens              10       9606         115             catty
    ## 6         Homo sapiens              10       9606         115             catty
    ##   INSDC_center_alias                                        INSDC_center_name
    ## 1             UCSDMI University of California San Diego Microbiome Initiative
    ## 2             UCSDMI University of California San Diego Microbiome Initiative
    ## 3             UCSDMI University of California San Diego Microbiome Initiative
    ## 4             UCSDMI University of California San Diego Microbiome Initiative
    ## 5             UCSDMI University of California San Diego Microbiome Initiative
    ## 6             UCSDMI University of California San Diego Microbiome Initiative
    ##     INSDC_first_public    INSDC_last_update INSDC_status          Instrument
    ## 1 2021-11-10T16:23:23Z 2021-11-10T16:23:23Z       public Illumina HiSeq 4000
    ## 2 2021-11-10T16:23:23Z 2021-11-10T16:23:23Z       public Illumina HiSeq 4000
    ## 3 2021-11-10T16:23:23Z 2021-11-10T16:23:23Z       public Illumina HiSeq 4000
    ## 4 2021-11-10T16:23:23Z 2021-11-10T16:23:23Z       public Illumina HiSeq 4000
    ## 5 2021-11-10T16:23:23Z 2021-11-10T16:23:23Z       public Illumina HiSeq 4000
    ## 6 2021-11-10T16:23:23Z 2021-11-10T16:23:23Z       public Illumina HiSeq 4000
    ##   LATITUDE Library.Name LibraryLayout LibrarySelection LibrarySource longitude
    ## 1    40.01    14057.10A        SINGLE              PCR   METAGENOMIC   -105.27
    ## 2    40.01    14057.10B        SINGLE              PCR   METAGENOMIC   -105.27
    ## 3    40.01    14057.10C        SINGLE              PCR   METAGENOMIC   -105.27
    ## 4    40.01    14057.10D        SINGLE              PCR   METAGENOMIC   -105.27
    ## 5    40.01    14057.10E        SINGLE              PCR   METAGENOMIC   -105.27
    ## 6    40.01    14057.10F        SINGLE              PCR   METAGENOMIC   -105.27
    ##                Organism                                 pcr_primers..exp.
    ## 1 human oral metagenome FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 2 human oral metagenome FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 3 human oral metagenome FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 4 human oral metagenome FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 5 human oral metagenome FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 6 human oral metagenome FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ##   physical_specimen_location Platform        primer..exp.          ReleaseDate
    ## 1                     UCSDMI ILLUMINA GTGCCAGCMGCCGCGGTAA 2021-12-02T00:00:00Z
    ## 2                     UCSDMI ILLUMINA GTGCCAGCMGCCGCGGTAA 2021-12-02T00:00:00Z
    ## 3                     UCSDMI ILLUMINA GTGCCAGCMGCCGCGGTAA 2021-12-02T00:00:00Z
    ## 4                     UCSDMI ILLUMINA GTGCCAGCMGCCGCGGTAA 2021-12-02T00:00:00Z
    ## 5                     UCSDMI ILLUMINA GTGCCAGCMGCCGCGGTAA 2021-12-02T00:00:00Z
    ## 6                     UCSDMI ILLUMINA GTGCCAGCMGCCGCGGTAA 2021-12-02T00:00:00Z
    ##   run_center..exp.          create_date version   Sample.Name
    ## 1             CCME 2021-11-16T17:03:00Z       1 SAMEA10792800
    ## 2             CCME 2021-11-16T17:03:00Z       1 SAMEA10792801
    ## 3             CCME 2021-11-16T17:03:00Z       1 SAMEA10792802
    ## 4             CCME 2021-11-16T17:03:00Z       1 SAMEA10792803
    ## 5             CCME 2021-11-16T17:04:00Z       1 SAMEA10792804
    ## 6             CCME 2021-11-16T17:05:00Z       1 SAMEA10792805
    ##                 Sample_name sample_type sequencing_method..exp. SRA.Study
    ## 1 qiita_sid_14057:14057.10A      saliva Sequencing by synthesis ERP133010
    ## 2 qiita_sid_14057:14057.10B      saliva Sequencing by synthesis ERP133010
    ## 3 qiita_sid_14057:14057.10C      saliva Sequencing by synthesis ERP133010
    ## 4 qiita_sid_14057:14057.10D      saliva Sequencing by synthesis ERP133010
    ## 5 qiita_sid_14057:14057.10E      saliva Sequencing by synthesis ERP133010
    ## 6 qiita_sid_14057:14057.10F      saliva Sequencing by synthesis ERP133010
    ##                Submitter_Id target_gene..exp. target_subfragment..exp.
    ## 1 qiita_sid_14057:14057.10A          16S rRNA                       V4
    ## 2 qiita_sid_14057:14057.10B          16S rRNA                       V4
    ## 3 qiita_sid_14057:14057.10C          16S rRNA                       V4
    ## 4 qiita_sid_14057:14057.10D          16S rRNA                       V4
    ## 5 qiita_sid_14057:14057.10E          16S rRNA                       V4
    ## 6 qiita_sid_14057:14057.10F          16S rRNA                       V4
    ##   timepoint
    ## 1         A
    ## 2         B
    ## 3         C
    ## 4         D
    ## 5         E
    ## 6         F

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
metadata<-data.frame(
  BioProject = initial_meta$BioProject,
  Run = initial_meta$Run,
  sample_name = initial_meta$Sample.Name,
  Age = initial_meta$Host_age,
  BMI = NA,  
  Sex = initial_meta$gender,
  Region = initial_meta$target_subfragment..exp.,
  Platform = initial_meta$Instrument,
  Site = initial_meta$env_material,
  public_time = initial_meta$ReleaseDate,
  matched_samples = NA,
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

    ## [1] "Number of unique runs: 229"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 229"

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
