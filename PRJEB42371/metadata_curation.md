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
setwd("E:/bioproject/PRJEB42371")

# Read the initial metadata files
initial_meta <- read.csv("SraRunTable.csv", header = TRUE) 
initial_meta <- initial_meta[initial_meta$env_material != "sterile water", ]

# Display the first few rows of the metadata
head(initial_meta)
```

    ##          Run Assay.Type AvgSpotLen barcode..exp.    Bases BioProject
    ## 3 ERR5050934   AMPLICON        250         AAAAA  8649606 PRJEB42371
    ## 4 ERR5050935   AMPLICON        250         AAAAA 12412209 PRJEB42371
    ## 5 ERR5050936   AMPLICON        250         AAAAA  5886944 PRJEB42371
    ## 6 ERR5050937   AMPLICON        250         AAAAA  7570303 PRJEB42371
    ## 7 ERR5050938   AMPLICON        250         AAAAA  5406358 PRJEB42371
    ## 8 ERR5050939   AMPLICON        250         AAAAA  7028778 PRJEB42371
    ##      BioSample   Bytes center_name..exp.
    ## 3 SAMEA7767094 6297508     U Puerto Rico
    ## 4 SAMEA7767095 9072165     U Puerto Rico
    ## 5 SAMEA7767096 3212096     U Puerto Rico
    ## 6 SAMEA7767097 4162448     U Puerto Rico
    ## 7 SAMEA7767098 3142792     U Puerto Rico
    ## 8 SAMEA7767099 3949993     U Puerto Rico
    ##                                          Center.Name
    ## 3 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT
    ## 4 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT
    ## 5 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT
    ## 6 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT
    ## 7 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT
    ## 8 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT
    ##     center_project_name..exp. collection_timestamp Consent     COUNTRY
    ## 3 ORAL MICROBIOTA PUERTO RICO                 2019  public Puerto Rico
    ## 4 ORAL MICROBIOTA PUERTO RICO                 2019  public Puerto Rico
    ## 5 ORAL MICROBIOTA PUERTO RICO                 2019  public Puerto Rico
    ## 6 ORAL MICROBIOTA PUERTO RICO                 2019  public Puerto Rico
    ## 7 ORAL MICROBIOTA PUERTO RICO                 2019  public Puerto Rico
    ## 8 ORAL MICROBIOTA PUERTO RICO                 2019  public Puerto Rico
    ##   DATASTORE.filetype DATASTORE.provider                     DATASTORE.region
    ## 3   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 4   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 5   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 6   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 7   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ## 8   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1
    ##   dna_extracted elevation          empo_1            empo_2           empo_3
    ## 3          TRUE         8 Host-associated Animal-associated Animal secretion
    ## 4          TRUE         8 Host-associated Animal-associated Animal secretion
    ## 5          TRUE         8 Host-associated Animal-associated Animal secretion
    ## 6          TRUE         8 Host-associated Animal-associated Animal secretion
    ## 7          TRUE         8 Host-associated Animal-associated Animal secretion
    ## 8          TRUE         8 Host-associated Animal-associated Animal secretion
    ##   ENA.FIRST.PUBLIC..run. ENA.LAST.UPDATE..run.   env_biome
    ## 3              2021/1/13              2021/1/6 urban biome
    ## 4              2021/1/13              2021/1/6 urban biome
    ## 5              2021/1/13              2021/1/6 urban biome
    ## 6              2021/1/13              2021/1/6 urban biome
    ## 7              2021/1/13              2021/1/6 urban biome
    ## 8              2021/1/13              2021/1/6 urban biome
    ##                env_feature env_material env_package Experiment  external_id
    ## 3 human-associated habitat       saliva  human-oral ERX4857266 SAMEA7767094
    ## 4 human-associated habitat       saliva  human-oral ERX4857267 SAMEA7767095
    ## 5 human-associated habitat       saliva  human-oral ERX4857268 SAMEA7767096
    ## 6 human-associated habitat       saliva  human-oral ERX4857269 SAMEA7767097
    ## 7 human-associated habitat       saliva  human-oral ERX4857270 SAMEA7767098
    ## 8 human-associated habitat       saliva  human-oral ERX4857271 SAMEA7767099
    ##   geo_loc_name_country geo_loc_name_country_continent         geo_loc_name
    ## 3          Puerto Rico                  North America Puerto Rico:San Juan
    ## 4          Puerto Rico                  North America Puerto Rico:San Juan
    ## 5          Puerto Rico                  North America Puerto Rico:San Juan
    ## 6          Puerto Rico                  North America Puerto Rico:San Juan
    ## 7          Puerto Rico                  North America Puerto Rico:San Juan
    ## 8          Puerto Rico                  North America Puerto Rico:San Juan
    ##   host_age_units host_height_units host_subject_id host_weight_units
    ## 3          years                cm              O1                kg
    ## 4          years                cm             O10                kg
    ## 5          years                cm            O100                kg
    ## 6          years                cm            O101                kg
    ## 7          years                cm            O102                kg
    ## 8          years                cm            O103                kg
    ##                                          INSDC_center_name   INSDC_first_public
    ## 3 University of California San Diego Microbiome Initiative 2021-01-13T04:11:46Z
    ## 4 University of California San Diego Microbiome Initiative 2021-01-13T04:11:46Z
    ## 5 University of California San Diego Microbiome Initiative 2021-01-13T04:11:46Z
    ## 6 University of California San Diego Microbiome Initiative 2021-01-13T04:11:46Z
    ## 7 University of California San Diego Microbiome Initiative 2021-01-13T04:11:46Z
    ## 8 University of California San Diego Microbiome Initiative 2021-01-13T04:11:46Z
    ##      INSDC_last_update INSDC_status     Instrument latitude Library.Name
    ## 3 2021-01-06T13:48:14Z       public Illumina MiSeq 18.39837     13193.O1
    ## 4 2021-01-06T13:48:14Z       public Illumina MiSeq 18.39837    13193.O10
    ## 5 2021-01-06T13:48:14Z       public Illumina MiSeq 18.39837   13193.O100
    ## 6 2021-01-06T13:48:14Z       public Illumina MiSeq 18.39837   13193.O101
    ## 7 2021-01-06T13:48:14Z       public Illumina MiSeq 18.39837   13193.O102
    ## 8 2021-01-06T13:48:14Z       public Illumina MiSeq 18.39837   13193.O103
    ##   LibraryLayout LibrarySelection LibrarySource longitude
    ## 3        SINGLE              PCR   METAGENOMIC -66.07453
    ## 4        SINGLE              PCR   METAGENOMIC -66.07453
    ## 5        SINGLE              PCR   METAGENOMIC -66.07453
    ## 6        SINGLE              PCR   METAGENOMIC -66.07453
    ## 7        SINGLE              PCR   METAGENOMIC -66.07453
    ## 8        SINGLE              PCR   METAGENOMIC -66.07453
    ##                  Organism                         physical_specimen_location
    ## 3 human saliva metagenome Microbial Ecology and Genomics Lab RCM San Juan PR
    ## 4 human saliva metagenome Microbial Ecology and Genomics Lab RCM San Juan PR
    ## 5 human saliva metagenome Microbial Ecology and Genomics Lab RCM San Juan PR
    ## 6 human saliva metagenome Microbial Ecology and Genomics Lab RCM San Juan PR
    ## 7 human saliva metagenome Microbial Ecology and Genomics Lab RCM San Juan PR
    ## 8 human saliva metagenome Microbial Ecology and Genomics Lab RCM San Juan PR
    ##   physical_specimen_remaining Platform          ReleaseDate
    ## 3                        TRUE ILLUMINA 2021-01-18T00:00:00Z
    ## 4                        TRUE ILLUMINA 2021-01-18T00:00:00Z
    ## 5                        TRUE ILLUMINA 2021-01-18T00:00:00Z
    ## 6                        TRUE ILLUMINA 2021-01-18T00:00:00Z
    ## 7                        TRUE ILLUMINA 2021-01-18T00:00:00Z
    ## 8                        TRUE ILLUMINA 2021-01-18T00:00:00Z
    ##             run_center..exp.          create_date version
    ## 3 Louisiana State University 2021-01-28T10:32:00Z       1
    ## 4 Louisiana State University 2021-01-28T10:32:00Z       1
    ## 5 Louisiana State University 2021-01-28T10:31:00Z       1
    ## 6 Louisiana State University 2021-01-28T10:31:00Z       1
    ## 7 Louisiana State University 2021-01-28T10:31:00Z       1
    ## 8 Louisiana State University 2021-01-28T10:32:00Z       1
    ##           run_prefix..exp.  Sample.Name                Sample_name sample_type
    ## 3      GodO01_S175_L001_R1 SAMEA7767094   qiita_sid_13193:13193.O1      saliva
    ## 4      GodO10_S184_L001_R1 SAMEA7767095  qiita_sid_13193:13193.O10      saliva
    ## 5 VitO100_S221_L001_R1_001 SAMEA7767096 qiita_sid_13193:13193.O100      saliva
    ## 6 VitO101_S233_L001_R1_001 SAMEA7767097 qiita_sid_13193:13193.O101      saliva
    ## 7 VitO102_S245_L001_R1_001 SAMEA7767098 qiita_sid_13193:13193.O102      saliva
    ## 8 VitO103_S257_L001_R1_001 SAMEA7767099 qiita_sid_13193:13193.O103      saliva
    ##     sequencing_meth..exp. SRA.Study               Submitter_Id
    ## 3 Sequencing by synthesis ERP126217   qiita_sid_13193:13193.O1
    ## 4 Sequencing by synthesis ERP126217  qiita_sid_13193:13193.O10
    ## 5 Sequencing by synthesis ERP126217 qiita_sid_13193:13193.O100
    ## 6 Sequencing by synthesis ERP126217 qiita_sid_13193:13193.O101
    ## 7 Sequencing by synthesis ERP126217 qiita_sid_13193:13193.O102
    ## 8 Sequencing by synthesis ERP126217 qiita_sid_13193:13193.O103
    ##   target_gene..exp. target_subfragment..exp. age_category
    ## 3          16S rRNA                       V4        31_40
    ## 4          16S rRNA                       V4        41_49
    ## 5          16S rRNA                       V4        21_30
    ## 6          16S rRNA                       V4        31_40
    ## 7          16S rRNA                       V4        31_40
    ## 8          16S rRNA                       V4        21_30
    ##   alcohol_consumption_last_month     Host_age host_body_habitat
    ## 3                            yes not provided     UBERON:saliva
    ## 4                             no not provided     UBERON:saliva
    ## 5                            yes not provided     UBERON:saliva
    ## 6                            yes not provided     UBERON:saliva
    ## 7                            yes not provided     UBERON:saliva
    ## 8                            yes not provided     UBERON:saliva
    ##   host_body_mass_index host_body_product host_body_site host_common_name
    ## 3         not provided     UBERON:saliva  UBERON:saliva            human
    ## 4         not provided     UBERON:saliva  UBERON:saliva            human
    ## 5         not provided     UBERON:saliva  UBERON:saliva            human
    ## 6         not provided     UBERON:saliva  UBERON:saliva            human
    ## 7         not provided     UBERON:saliva  UBERON:saliva            human
    ## 8         not provided     UBERON:saliva  UBERON:saliva            human
    ##    host_height host_life_stage host_scientific_name host_taxid  host_weight
    ## 3 not provided    not provided         Homo sapiens       9606 not provided
    ## 4 not provided    not provided         Homo sapiens       9606 not provided
    ## 5 not provided    not provided         Homo sapiens       9606 not provided
    ## 6 not provided    not provided         Homo sapiens       9606 not provided
    ## 7 not provided    not provided         Homo sapiens       9606 not provided
    ## 8 not provided    not provided         Homo sapiens       9606 not provided
    ##   perio_3_categories    sex Smoking perio_4_categories ENA_first_public
    ## 3                 no female      no                 no        2021/1/13
    ## 4                 no female      no                 no        2021/1/13
    ## 5                 no   male      no                 no        2021/1/13
    ## 6                 no   male      no                 no        2021/1/13
    ## 7                 no   male      no                 no        2021/1/13
    ## 8                 no   male     yes                 no        2021/1/13
    ##   ena_last_update INSDC_center_alias oral_sex_practices        primer..exp.
    ## 3        2021/1/6             UCSDMI                yes GTGCCAGCMGCCGCGGTAA
    ## 4        2021/1/6             UCSDMI                yes GTGCCAGCMGCCGCGGTAA
    ## 5        2021/1/6             UCSDMI                yes GTGCCAGCMGCCGCGGTAA
    ## 6        2021/1/6             UCSDMI                yes GTGCCAGCMGCCGCGGTAA
    ## 7        2021/1/6             UCSDMI                yes GTGCCAGCMGCCGCGGTAA
    ## 8        2021/1/6             UCSDMI                yes GTGCCAGCMGCCGCGGTAA
    ##   linker..exp.                                  PCR_Primer..exp.
    ## 3           GT FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT
    ## 4           GT FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT
    ## 5           GT FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT
    ## 6           GT FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT
    ## 7           GT FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT
    ## 8           GT FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT
    ##   barcodesequence..exp. ena_first_public ena_last_update.1 Scientific_Name
    ## 3                                                                         
    ## 4                                                                         
    ## 5                                                                         
    ## 6                                                                         
    ## 7                                                                         
    ## 8

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
  Age = initial_meta$age_category, 
  BMI = NA,  
  Sex = initial_meta$sex,
  Region = initial_meta$target_subfragment..exp.,
  Platform = initial_meta$Instrument,
  Country = initial_meta$COUNTRY,
  Site = initial_meta$env_material,
  public_time = initial_meta$ReleaseDate,
  stringsAsFactors = FALSE
)

metadata <- metadata %>%
  mutate(Sex = recode(Sex, "female" = "F", "male" = "M"))
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$PCR_primers <- "515F-806R"
metadata$title <- "Ecological competition in the oral mycobiome of Hispanic adults living in Puerto Rico associates with periodontitis"
metadata$Disease_study <- "YES"

# Assign disease status based on sample_name patterns
metadata$Status <- "NA"

metadata$matched_samples <- "NO"
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 356"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 136"

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
