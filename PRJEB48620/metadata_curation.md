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
setwd("E:/bioproject/PRJEB48620")

# Read the initial metadata files
initial_meta <- read.csv("SraRunTable.csv", header = TRUE)

# Display the first few rows of the metadata
head(initial_meta)
```

    ##          Run Assay.Type AvgSpotLen barcode..exp.   Bases BioProject
    ## 1 ERR7247035   AMPLICON        151  TTCTCGGTTCTC 2997048 PRJEB48620
    ## 2 ERR7247036   AMPLICON        151  TTATCCAGTCCT 3211166 PRJEB48620
    ## 3 ERR7247037   AMPLICON        151  CACCGTGACACT 2622719 PRJEB48620
    ## 4 ERR7247038   AMPLICON        151  CAGAGCTAATTG 2195389 PRJEB48620
    ## 5 ERR7247039   AMPLICON        151  AGGTCATCTTGG 1902751 PRJEB48620
    ## 6 ERR7247040   AMPLICON        151  TGACGCCTCCAA 2332648 PRJEB48620
    ##       BioSample   Bytes center_name..exp.
    ## 1 SAMEA10792706 2000939              CCME
    ## 2 SAMEA10792707 2101166              CCME
    ## 3 SAMEA10792708 1694555              CCME
    ## 4 SAMEA10792709 1481510              CCME
    ## 5 SAMEA10792710 1277859              CCME
    ## 6 SAMEA10792711 1525907              CCME
    ##                                          Center.Name center_project_name..exp.
    ## 1 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT                      CCME
    ## 2 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT                      CCME
    ## 3 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT                      CCME
    ## 4 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT                      CCME
    ## 5 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT                      CCME
    ## 6 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT                      CCME
    ##   collection_timestamp Consent DATASTORE.filetype DATASTORE.provider
    ## 1           2013/10/16  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 2           2013/10/16  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 3           2013/10/16  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 4           2013/10/16  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5           2013/10/16  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 6           2013/10/16  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                       DATASTORE.region Day     Diet elevation          empo_1
    ## 1 gs.us-east1,ncbi.public,s3.us-east-1   1 omnivore    1624.1 Host-associated
    ## 2 gs.us-east1,ncbi.public,s3.us-east-1   1 omnivore    1624.1 Host-associated
    ## 3 gs.us-east1,ncbi.public,s3.us-east-1   1 omnivore    1624.1 Host-associated
    ## 4 gs.us-east1,ncbi.public,s3.us-east-1   1 omnivore    1624.1 Host-associated
    ## 5 gs.us-east1,ncbi.public,s3.us-east-1   1 omnivore    1624.1 Host-associated
    ## 6 gs.us-east1,ncbi.public,s3.us-east-1   1 omnivore    1624.1 Host-associated
    ##   empo_2           empo_3 ENA.FIRST.PUBLIC..run. ENA_first_public
    ## 1 Animal Animal secretion             2021/11/10       2021/11/10
    ## 2 Animal Animal secretion             2021/11/10       2021/11/10
    ## 3 Animal Animal secretion             2021/11/10       2021/11/10
    ## 4 Animal Animal secretion             2021/11/10       2021/11/10
    ## 5 Animal Animal secretion             2021/11/10       2021/11/10
    ## 6 Animal Animal secretion             2021/11/10       2021/11/10
    ##   ENA.LAST.UPDATE..run. ena_last_update        env_biome
    ## 1            2021/11/10      2021/11/10 ENVO:urban biome
    ## 2            2021/11/10      2021/11/10 ENVO:urban biome
    ## 3            2021/11/10      2021/11/10 ENVO:urban biome
    ## 4            2021/11/10      2021/11/10 ENVO:urban biome
    ## 5            2021/11/10      2021/11/10 ENVO:urban biome
    ## 6            2021/11/10      2021/11/10 ENVO:urban biome
    ##                     env_feature env_material env_package Experiment
    ## 1 ENVO:human-associated habitat       saliva  human-oral ERX6816494
    ## 2 ENVO:human-associated habitat       saliva  human-oral ERX6816495
    ## 3 ENVO:human-associated habitat       saliva  human-oral ERX6816496
    ## 4 ENVO:human-associated habitat       saliva  human-oral ERX6816497
    ## 5 ENVO:human-associated habitat       saliva  human-oral ERX6816498
    ## 6 ENVO:human-associated habitat       saliva  human-oral ERX6816499
    ##     external_id Host_age host_age_units host_body_habitat host_body_mass_index
    ## 1 SAMEA10792706       40          years     UBERON:saliva                 23.1
    ## 2 SAMEA10792707       40          years     UBERON:saliva                 23.1
    ## 3 SAMEA10792708       40          years     UBERON:saliva                 23.1
    ## 4 SAMEA10792709       40          years     UBERON:saliva                 23.1
    ## 5 SAMEA10792710       40          years     UBERON:saliva                 23.1
    ## 6 SAMEA10792711       40          years     UBERON:saliva                 23.1
    ##   host_body_product host_body_site host_common_name host_height
    ## 1     UBERON:saliva   UBERON:mouth            human         1.8
    ## 2     UBERON:saliva   UBERON:mouth            human         1.8
    ## 3     UBERON:saliva   UBERON:mouth            human         1.8
    ## 4     UBERON:saliva   UBERON:mouth            human         1.8
    ## 5     UBERON:saliva   UBERON:mouth            human         1.8
    ## 6     UBERON:saliva   UBERON:mouth            human         1.8
    ##   host_height_units host_scientific_name host_subject_id host_taxid
    ## 1                 m         Homo sapiens        2052:M11       9606
    ## 2                 m         Homo sapiens        2052:M11       9606
    ## 3                 m         Homo sapiens        2052:M11       9606
    ## 4                 m         Homo sapiens        2052:M11       9606
    ## 5                 m         Homo sapiens        2052:M11       9606
    ## 6                 m         Homo sapiens        2052:M11       9606
    ##     host_weight host_weight_units Hour INSDC_center_alias
    ## 1 not collected                kg    8             UCSDMI
    ## 2 not collected                kg   10             UCSDMI
    ## 3 not collected                kg   10             UCSDMI
    ## 4 not collected                kg   10             UCSDMI
    ## 5 not collected                kg   10             UCSDMI
    ## 6 not collected                kg   11             UCSDMI
    ##                                          INSDC_center_name   INSDC_first_public
    ## 1 University of California San Diego Microbiome Initiative 2021-11-10T16:23:23Z
    ## 2 University of California San Diego Microbiome Initiative 2021-11-10T16:23:23Z
    ## 3 University of California San Diego Microbiome Initiative 2021-11-10T16:23:23Z
    ## 4 University of California San Diego Microbiome Initiative 2021-11-10T16:23:23Z
    ## 5 University of California San Diego Microbiome Initiative 2021-11-10T16:23:23Z
    ## 6 University of California San Diego Microbiome Initiative 2021-11-10T16:23:23Z
    ##      INSDC_last_update INSDC_status          Instrument LATITUDE
    ## 1 2021-11-10T16:23:23Z       public Illumina HiSeq 4000    40.01
    ## 2 2021-11-10T16:23:23Z       public Illumina HiSeq 4000    40.01
    ## 3 2021-11-10T16:23:23Z       public Illumina HiSeq 4000    40.01
    ## 4 2021-11-10T16:23:23Z       public Illumina HiSeq 4000    40.01
    ## 5 2021-11-10T16:23:23Z       public Illumina HiSeq 4000    40.01
    ## 6 2021-11-10T16:23:23Z       public Illumina HiSeq 4000    40.01
    ##                  Library.Name LibraryLayout LibrarySelection LibrarySource
    ## 1  14058.1.M11.saliva.1168567        SINGLE              PCR   METAGENOMIC
    ## 2 14058.10.M11.saliva.1168573        SINGLE              PCR   METAGENOMIC
    ## 3 14058.11.M11.saliva.1168556        SINGLE              PCR   METAGENOMIC
    ## 4 14058.12.M11.saliva.1168590        SINGLE              PCR   METAGENOMIC
    ## 5 14058.13.M11.saliva.1168512        SINGLE              PCR   METAGENOMIC
    ## 6 14058.14.M11.saliva.1168536        SINGLE              PCR   METAGENOMIC
    ##   life_stage longitude mf_sample_number min              Organism
    ## 1      adult   -105.27                1  25 human oral metagenome
    ## 2      adult   -105.27               10   4 human oral metagenome
    ## 3      adult   -105.27               11  19 human oral metagenome
    ## 4      adult   -105.27               12  33 human oral metagenome
    ## 5      adult   -105.27               13  48 human oral metagenome
    ## 6      adult   -105.27               14   3 human oral metagenome
    ##                                   pcr_primers..exp. physical_specimen_location
    ## 1 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT                     UCSDMI
    ## 2 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT                     UCSDMI
    ## 3 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT                     UCSDMI
    ## 4 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT                     UCSDMI
    ## 5 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT                     UCSDMI
    ## 6 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT                     UCSDMI
    ##   Platform        primer..exp.          ReleaseDate run_center..exp.
    ## 1 ILLUMINA GTGCCAGCMGCCGCGGTAA 2021-12-02T00:00:00Z             CCME
    ## 2 ILLUMINA GTGCCAGCMGCCGCGGTAA 2021-12-02T00:00:00Z             CCME
    ## 3 ILLUMINA GTGCCAGCMGCCGCGGTAA 2021-12-02T00:00:00Z             CCME
    ## 4 ILLUMINA GTGCCAGCMGCCGCGGTAA 2021-12-02T00:00:00Z             CCME
    ## 5 ILLUMINA GTGCCAGCMGCCGCGGTAA 2021-12-02T00:00:00Z             CCME
    ## 6 ILLUMINA GTGCCAGCMGCCGCGGTAA 2021-12-02T00:00:00Z             CCME
    ##            create_date version   Sample.Name
    ## 1 2021-11-16T17:00:00Z       1 SAMEA10792706
    ## 2 2021-11-16T17:00:00Z       1 SAMEA10792707
    ## 3 2021-11-16T17:02:00Z       1 SAMEA10792708
    ## 4 2021-11-16T17:00:00Z       1 SAMEA10792709
    ## 5 2021-11-16T17:00:00Z       1 SAMEA10792710
    ## 6 2021-11-16T17:03:00Z       1 SAMEA10792711
    ##                                   Sample_name sample_type
    ## 1  qiita_sid_14058:14058.1.M11.saliva.1168567      saliva
    ## 2 qiita_sid_14058:14058.10.M11.saliva.1168573      saliva
    ## 3 qiita_sid_14058:14058.11.M11.saliva.1168556      saliva
    ## 4 qiita_sid_14058:14058.12.M11.saliva.1168590      saliva
    ## 5 qiita_sid_14058:14058.13.M11.saliva.1168512      saliva
    ## 6 qiita_sid_14058:14058.14.M11.saliva.1168536      saliva
    ##   sequencing_method..exp.  sex smoker SRA.Study
    ## 1 Sequencing by synthesis male      n ERP133009
    ## 2 Sequencing by synthesis male      n ERP133009
    ## 3 Sequencing by synthesis male      n ERP133009
    ## 4 Sequencing by synthesis male      n ERP133009
    ## 5 Sequencing by synthesis male      n ERP133009
    ## 6 Sequencing by synthesis male      n ERP133009
    ##                                  Submitter_Id target_gene..exp.
    ## 1  qiita_sid_14058:14058.1.M11.saliva.1168567          16S rRNA
    ## 2 qiita_sid_14058:14058.10.M11.saliva.1168573          16S rRNA
    ## 3 qiita_sid_14058:14058.11.M11.saliva.1168556          16S rRNA
    ## 4 qiita_sid_14058:14058.12.M11.saliva.1168590          16S rRNA
    ## 5 qiita_sid_14058:14058.13.M11.saliva.1168512          16S rRNA
    ## 6 qiita_sid_14058:14058.14.M11.saliva.1168536          16S rRNA
    ##   target_subfragment..exp. activity
    ## 1                       V4    sleep
    ## 2                       V4         
    ## 3                       V4      muc
    ## 4                       V4      eat
    ## 5                       V4         
    ## 6                       V4

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
  BMI = initial_meta$host_body_mass_index,  
  Sex = initial_meta$sex,
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

    ## [1] "Number of unique runs: 94"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 94"

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
