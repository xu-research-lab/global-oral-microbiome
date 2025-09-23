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
setwd("D:/oral_project/qiita_id1928")

# Read the initial metadata files
initial_meta<-read.csv("SraRunTable.csv",header = T)

# Combine the two metadata files
#meta <- rbind(initial_meta, initial_meta1)

#Remove entries non-oral samples
#initial_meta <- subset(initial_meta,isolate %in% c("saliva","tongue coat"))

# Check the number of samples containing "HT" in Sample.Name (e.g., healthy controls)
#nrow(subset(meta, grepl("HT", Sample.Name)))

# Display the first few rows of the metadata
head(initial_meta)
```

    ##         Run Assay.Type AvgSpotLen     Bases BioProject    BioSample     Bytes
    ## 1 SRR052699      OTHER        561 261276500 PRJNA48333 SAMN00010557 568521673
    ## 2 SRR057609      OTHER        565 267471073 PRJNA48333 SAMN00010557 605430345
    ## 3 SRR057663      OTHER        552 246151846 PRJNA48333 SAMN00010557 553350153
    ## 4 SRR058084      OTHER        569 283226778 PRJNA48333 SAMN00010557 623404038
    ## 5 SRR058085      OTHER        561 267408754 PRJNA48333 SAMN00010557 600603916
    ## 6 SRR058086      OTHER        560 257324330 PRJNA48333 SAMN00010557 556681911
    ##   Center.Name Consent DATASTORE.filetype DATASTORE.provider
    ## 1        JCVI  public         run.zq,sra         gs,ncbi,s3
    ## 2        JCVI  public         run.zq,sra         gs,ncbi,s3
    ## 3        JCVI  public         run.zq,sra         gs,ncbi,s3
    ## 4        JCVI  public         run.zq,sra         gs,ncbi,s3
    ## 5        JCVI  public         run.zq,sra         gs,ncbi,s3
    ## 6        JCVI  public         run.zq,sra         gs,ncbi,s3
    ##                                  DATASTORE.region Experiment
    ## 1 gs.us-east1,ncbi.dbgap,ncbi.public,s3.us-east-1  SRX021470
    ## 2 gs.us-east1,ncbi.dbgap,ncbi.public,s3.us-east-1  SRX022036
    ## 3 gs.us-east1,ncbi.dbgap,ncbi.public,s3.us-east-1  SRX022097
    ## 4 gs.us-east1,ncbi.dbgap,ncbi.public,s3.us-east-1  SRX022221
    ## 5 gs.us-east1,ncbi.dbgap,ncbi.public,s3.us-east-1  SRX022222
    ## 6 gs.us-east1,ncbi.dbgap,ncbi.public,s3.us-east-1  SRX022223
    ##              gene..exp.          Instrument        Library.Name
    ## 1 16S rRNA V3-V5 region 454 GS FLX Titanium 2ZAL8V5V3-BA-01-652
    ## 2 16S rRNA V3-V5 region 454 GS FLX Titanium 2ZARYV5V3-BA-01-659
    ## 3 16S rRNA V1-V3 region 454 GS FLX Titanium  2ZAGGHMP-BA-01-630
    ## 4 16S rRNA V3-V5 region 454 GS FLX Titanium 2ZB85V5V3-BA-01-675
    ## 5 16S rRNA V3-V5 region 454 GS FLX Titanium 2ZARXV5V3-BA-01-655
    ## 6 16S rRNA V3-V5 region 454 GS FLX Titanium 2ZAL7V5V3-BA-01-658
    ##   library_strategy..exp. LibraryLayout LibrarySelection LibrarySource
    ## 1         targeted-locus        SINGLE              PCR       GENOMIC
    ## 2         targeted-locus        SINGLE              PCR       GENOMIC
    ## 3         targeted-locus        SINGLE              PCR       GENOMIC
    ## 4         targeted-locus        SINGLE              PCR       GENOMIC
    ## 5         targeted-locus        SINGLE              PCR       GENOMIC
    ## 6         targeted-locus        SINGLE              PCR       GENOMIC
    ##       Organism Platform          ReleaseDate          create_date version
    ## 1 unidentified    LS454 2010-05-27T00:00:00Z 2018-09-20T17:51:00Z       3
    ## 2 unidentified    LS454 2010-06-15T00:00:00Z 2018-09-19T12:50:00Z       3
    ## 3 unidentified    LS454 2010-06-16T00:00:00Z 2018-09-19T12:52:00Z       3
    ## 4 unidentified    LS454 2010-06-17T00:00:00Z 2018-09-19T12:49:00Z       3
    ## 5 unidentified    LS454 2010-06-17T00:00:00Z 2018-09-19T12:52:00Z       3
    ## 6 unidentified    LS454 2010-06-17T00:00:00Z 2018-09-19T12:49:00Z       3
    ##              Sample.Name SRA.Study analyte_type biospecimen_repository
    ## 1 unidentified-protected SRP002395                                    
    ## 2 unidentified-protected SRP002395                                    
    ## 3 unidentified-protected SRP002395                                    
    ## 4 unidentified-protected SRP002395                                    
    ## 5 unidentified-protected SRP002395                                    
    ## 6 unidentified-protected SRP002395                                    
    ##   biospecimen_repository_sample_id Consent.code dbGaP.accession host_sex
    ## 1                               NA           NA                         
    ## 2                               NA           NA                         
    ## 3                               NA           NA                         
    ## 4                               NA           NA                         
    ## 5                               NA           NA                         
    ## 6                               NA           NA                         
    ##   isolation_source study_design study_name submitted_subject_id
    ## 1                                                            NA
    ## 2                                                            NA
    ## 3                                                            NA
    ## 4                                                            NA
    ## 5                                                            NA
    ## 6                                                            NA
    ##   instrument_name..run. project..run. region..run. run_barcode..run.
    ## 1                                               NA                NA
    ## 2                                               NA                NA
    ## 3                                               NA                NA
    ## 4                                               NA                NA
    ## 5                                               NA                NA
    ## 6                                               NA                NA
    ##   run_name..run. gssr_id..run. lsid..run. molecular_indexing_scheme..run.
    ## 1                           NA                                           
    ## 2                           NA                                           
    ## 3                           NA                                           
    ## 4                           NA                                           
    ## 5                           NA                                           
    ## 6                           NA                                           
    ##   work_request..run. BI_project_name..exp. BI_work_request_ID..exp.
    ## 1                 NA                                             NA
    ## 2                 NA                                             NA
    ## 3                 NA                                             NA
    ## 4                 NA                                             NA
    ## 5                 NA                                             NA
    ## 6                 NA                                             NA
    ##   BI_GSSR_sample_IDs..exp. BI_GSSR_sample_LSIDs..exp. material_type..exp.
    ## 1                                                                        
    ## 2                                                                        
    ## 3                                                                        
    ## 4                                                                        
    ## 5                                                                        
    ## 6                                                                        
    ##   project..exp. gssr_id..exp. lsid..exp. work_request..exp. work_requests..exp.
    ## 1                          NA                            NA                    
    ## 2                          NA                            NA                    
    ## 3                          NA                            NA                    
    ## 4                          NA                            NA                    
    ## 5                          NA                            NA                    
    ## 6                          NA                            NA                    
    ##   work_requests..run.
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

``` r
# Create a uniform metadata dataframe
metadata <- data.frame(
  BioProject = initial_meta$BioProject,
  Run = initial_meta$Run,
  sample_name = initial_meta$Sample.Name,
  Age = NA,
  BMI = NA,
  Sex=initial_meta$host_sex,
  Region = initial_meta$gene..exp.,
  Platform = initial_meta$Instrument,
  Country = "USA",
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
metadata$PCR_primers <- NA
metadata$title <- "Structure function and diversity of the
healthy human microbiome"
metadata$Disease_study <- "NO"
# Assign disease status based on sample_name patterns
metadata$Status<- NA

table(metadata$Site)
```

    ## 
    ##                                                  G_DNA_Anterior nares 
    ##                                 80                                473 
    ## G_DNA_Attached/Keratinized gingiva                G_DNA_Buccal mucosa 
    ##                                452                                453 
    ##                  G_DNA_Hard palate          G_DNA_L_Antecubital fossa 
    ##                                460                                467 
    ##      G_DNA_L_Retroauricular crease                   G_DNA_Mid vagina 
    ##                                490                                197 
    ##             G_DNA_Palatine Tonsils             G_DNA_Posterior fornix 
    ##                                461                                193 
    ##          G_DNA_R_Antecubital fossa      G_DNA_R_Retroauricular crease 
    ##                                468                                502 
    ##                       G_DNA_Saliva                        G_DNA_Stool 
    ##                                393                                455 
    ##           G_DNA_Subgingival plaque         G_DNA_Supragingival plaque 
    ##                                445                                453 
    ##                       G_DNA_Throat                G_DNA_Tongue dorsum 
    ##                                449                                451 
    ##            G_DNA_Vaginal introitus 
    ##                                185

``` r
### subset oral samples

metadata<-subset(metadata,Site %in% c("G_DNA_Saliva","G_DNA_Supragingival plaque","G_DNA_Hard palate","G_DNA_Attached/Keratinized gingiva","G_DNA_Subgingival plaque","G_DNA_Tongue dorsum"))
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 2654"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 1465"

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
