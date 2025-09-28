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
16. **matched_samples**: matched samples (“YES” or “NO”)

``` r
# Create a uniform metadata dataframe
metadata <- data.frame(
  BioProject = initial_meta$BioProject,
  Run = initial_meta$Run,
  sample_name = initial_meta$submitted_subject_id,
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

# Check initial data distribution
cat("Initial sample distribution:\n")
```

    ## Initial sample distribution:

``` r
table(metadata$sample_name)
```

    ## 
    ## 103092734 132902142 147406386 158013734 158114885 158155345 158216035 158236265 
    ##         5        21        28        39        42        29        30        39 
    ## 158256496 158276726 158337416 158357646 158398106 158418336 158438567 158458797 
    ##        30        38        26        18        37        34        32        39 
    ## 158479027 158499257 158721788 158742018 158802708 158822939 158883629 158924089 
    ##        17        29        20        36        83        55        38        86 
    ## 158944319 158964549 158984779 159005010 159085930 159146620 159166850 159207311 
    ##        36       101        21       123       106       100        93        57 
    ## 159227541 159247771 159268001 159288231 159308461 159328691 159369152 159389382 
    ##        73        99        37        35        53        21        41        57 
    ## 159429842 159450072 159470302 159490532 159510762 159551223 159571453 159591683 
    ##        23        37        75        36        57        66        50        33 
    ## 159611913 159632143 159672603 159713063 159733294 159753524 159814214 159915365 
    ##        39        35        51       113        77        57        37        41 
    ## 160016515 160036745 160056975 160097436 160137896 160158126 160178356 160218816 
    ##        30        23        33        38        25        46        31        49 
    ## 160239046 160259276 160319967 160340197 160380657 160400887 160421117 160441347 
    ##        16        24        41        16        40        47        38        19 
    ## 160461578 160481808 160502038 160542498 160582958 160603188 160643649 160663879 
    ##        28        31        55        19        43        37        57        30 
    ## 160684109 160704339 160765029 160845950 178713055 184349034 208027353 250593475 
    ##        19        44        37        48        26        14        33         2 
    ## 256789458 257905678 275382046 289996019 338793263 355657046 370027359 375450439 
    ##        24        18        18        45         7        27        30        57 
    ## 404239096 451588811 465578759 508703490 516889361 517810313 550534656 553359145 
    ##        28         3        31        33        14        29         8        31 
    ## 561079058 612472597 643185023 658594300 663835652 668248235 675950834 682102541 
    ##         7        34        24        17        27        29         5        10 
    ## 682449369 706846339 737052003 739574095 763395383 763435843 763456073 763476303 
    ##        33         8        29        12        24        15        30        34 
    ## 763516763 763536994 763557224 763597684 763638144 763678604 763698834 763719065 
    ##        58       120        26       123        58        60       140        63 
    ## 763721271 763759525 763820215 763840445 763860675 763880905 763901136 763921366 
    ##        18        63        56        70        60        73        70        37 
    ## 763982056 764002286 764042746 764062976 764083206 764143897 764184357 764224817 
    ##        50        65        57        53       120        52        71        44 
    ## 764245047 764285508 764305738 764325968 764346198 764366428 764447348 764467579 
    ##        34        20        60        43        64        32        43        51 
    ## 764487809 764508039 764588959 764649650 764669880 764710340 764750800 764811490 
    ##        44        31        31        37        61        18        27        26 
    ## 764831721 764872181 764892411 764953101 765013792 765034022 765074482 765094712 
    ##        26        27        30        24        34        35        45        54 
    ## 765135172 765155402 765195863 765216093 765256553 765276783 765317243 765337473 
    ##        29        19        22        21        21        24        18        25 
    ## 765377934 765398164 765519544 765539774 765560005 765620695 765640925 765661155 
    ##        14         5         8         6        30         6         4        10 
    ## 765701615 823052294 857980665 861967750 863126187 866653606 872848689 875002022 
    ##        28         4         6        28        22        18        18         8 
    ## 901775393 937495960 938202701 953045535 
    ##         4         4         2        28

``` r
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
# Define oral sampling sites
oral_sites <- c("G_DNA_Saliva","G_DNA_Supragingival plaque","G_DNA_Hard palate",
                "G_DNA_Attached/Keratinized gingiva","G_DNA_Subgingival plaque",
                "G_DNA_Tongue dorsum")

# Identify subjects with at least one oral sample
subjects_with_oral <- unique(metadata$sample_name[metadata$Site %in% oral_sites])
cat("\nSubjects with oral sampling:", length(subjects_with_oral), "\n")
```

    ## 
    ## Subjects with oral sampling: 174

``` r
# Keep only subjects that have at least one oral sample
metadata_clean <- metadata[metadata$sample_name %in% subjects_with_oral, ]

# Recalculate sample counts after filtering
sample_counts_clean <- table(metadata_clean$sample_name)

# For each subject, check if they have both oral and non-oral samples
metadata_clean$matched_samples <- "NO"

for (subject in unique(metadata_clean$sample_name)) {
    subject_data <- metadata_clean[metadata_clean$sample_name == subject, ]
    
    has_oral <- any(subject_data$Site %in% oral_sites)
    has_non_oral <- any(!subject_data$Site %in% oral_sites)
    
    if (has_oral && has_non_oral) {
        metadata_clean$matched_samples[metadata_clean$sample_name == subject] <- "YES"
    }
}

cat("\nFinal sample distribution after cleaning:\n")
```

    ## 
    ## Final sample distribution after cleaning:

``` r
cat("Number of subjects:", length(unique(metadata_clean$sample_name)), "\n")
```

    ## Number of subjects: 174

``` r
cat("Total samples:", nrow(metadata_clean), "\n")
```

    ## Total samples: 7298

``` r
cat("\nSample distribution by subject:\n")
```

    ## 
    ## Sample distribution by subject:

``` r
print(table(metadata_clean$sample_name))
```

    ## 
    ## 132902142 147406386 158013734 158114885 158155345 158216035 158236265 158256496 
    ##        21        28        39        42        29        30        39        30 
    ## 158276726 158337416 158357646 158398106 158418336 158438567 158458797 158479027 
    ##        38        26        18        37        34        32        39        17 
    ## 158499257 158721788 158742018 158802708 158822939 158883629 158924089 158944319 
    ##        29        20        36        83        55        38        86        36 
    ## 158964549 158984779 159005010 159085930 159146620 159166850 159207311 159227541 
    ##       101        21       123       106       100        93        57        73 
    ## 159247771 159268001 159288231 159308461 159328691 159369152 159389382 159429842 
    ##        99        37        35        53        21        41        57        23 
    ## 159450072 159470302 159490532 159510762 159551223 159571453 159591683 159611913 
    ##        37        75        36        57        66        50        33        39 
    ## 159632143 159672603 159713063 159733294 159753524 159814214 159915365 160016515 
    ##        35        51       113        77        57        37        41        30 
    ## 160036745 160056975 160097436 160137896 160158126 160178356 160218816 160239046 
    ##        23        33        38        25        46        31        49        16 
    ## 160259276 160319967 160340197 160380657 160400887 160421117 160441347 160461578 
    ##        24        41        16        40        47        38        19        28 
    ## 160481808 160502038 160542498 160582958 160603188 160643649 160663879 160684109 
    ##        31        55        19        43        37        57        30        19 
    ## 160704339 160765029 160845950 178713055 184349034 208027353 256789458 257905678 
    ##        44        37        48        26        14        33        24        18 
    ## 275382046 289996019 355657046 370027359 375450439 404239096 451588811 465578759 
    ##        18        45        27        30        57        28         3        31 
    ## 508703490 516889361 517810313 553359145 612472597 643185023 658594300 663835652 
    ##        33        14        29        31        34        24        17        27 
    ## 668248235 682449369 737052003 763395383 763435843 763456073 763476303 763516763 
    ##        29        33        29        24        15        30        34        58 
    ## 763536994 763557224 763597684 763638144 763678604 763698834 763719065 763721271 
    ##       120        26       123        58        60       140        63        18 
    ## 763759525 763820215 763840445 763860675 763880905 763901136 763921366 763982056 
    ##        63        56        70        60        73        70        37        50 
    ## 764002286 764042746 764062976 764083206 764143897 764184357 764224817 764245047 
    ##        65        57        53       120        52        71        44        34 
    ## 764285508 764305738 764325968 764346198 764366428 764447348 764467579 764487809 
    ##        20        60        43        64        32        43        51        44 
    ## 764508039 764588959 764649650 764669880 764750800 764811490 764831721 764872181 
    ##        31        31        37        61        27        26        26        27 
    ## 764892411 764953101 765013792 765034022 765074482 765094712 765135172 765155402 
    ##        30        24        34        35        45        54        29        19 
    ## 765195863 765216093 765256553 765276783 765317243 765337473 765377934 765560005 
    ##        22        21        21        24        18        25        14        30 
    ## 765701615 861967750 863126187 866653606 872848689 953045535 
    ##        28        28        22        18        18        28

``` r
cat("\nSite distribution:\n")
```

    ## 
    ## Site distribution:

``` r
print(table(metadata_clean$Site))
```

    ## 
    ##               G_DNA_Anterior nares G_DNA_Attached/Keratinized gingiva 
    ##                                456                                452 
    ##                G_DNA_Buccal mucosa                  G_DNA_Hard palate 
    ##                                453                                460 
    ##          G_DNA_L_Antecubital fossa      G_DNA_L_Retroauricular crease 
    ##                                448                                467 
    ##                   G_DNA_Mid vagina             G_DNA_Palatine Tonsils 
    ##                                187                                461 
    ##             G_DNA_Posterior fornix          G_DNA_R_Antecubital fossa 
    ##                                186                                448 
    ##      G_DNA_R_Retroauricular crease                       G_DNA_Saliva 
    ##                                479                                393 
    ##                        G_DNA_Stool           G_DNA_Subgingival plaque 
    ##                                435                                445 
    ##         G_DNA_Supragingival plaque                       G_DNA_Throat 
    ##                                453                                448 
    ##                G_DNA_Tongue dorsum            G_DNA_Vaginal introitus 
    ##                                451                                176

``` r
cat("\nMatched samples distribution:\n")
```

    ## 
    ## Matched samples distribution:

``` r
print(table(metadata_clean$matched_samples))
```

    ## 
    ##  YES 
    ## 7298

``` r
# Use cleaned dataset for subsequent analysis
metadata <- metadata_clean
# renames sequencing regions
metadata$Region<-gsub("16s rRNA | region|16S rRNA |16S ","",metadata$Region)
metadata$Region<-gsub(",","/",metadata$Region)
metadata$Site<-gsub("G_DNA_","",metadata$Site)
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 7298"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 174"

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
        theme(axis.text.x = element_text(angle = 35, hjust = 1))
    }
    plots[[col]] <- p
  }
  
  return(plots)
}

# Generate and arrange plots
plots <- create_simple_plots(plot_data)
row1 <- plot_grid(plots$Age, plots$BMI, plots$Sex, nrow = 1, rel_widths = c(1, 1, 1))
row2 <- plot_grid(plots$Region, plots$Platform, plots$Country, nrow = 1, rel_widths = c(3, 1, 1))
row3 <- plot_grid(plots$Site, plots$public_time, plots$PCR_primers, nrow = 1, rel_widths = c(2, 2, 1))
row4<- plot_grid(plots$Disease_study, plots$Status, plots$matched_samples,nrow = 1, rel_widths = c(1,1,1))

p <- plot_grid(row1, row2, row3, row4,ncol = 1)
ggsave("Metadata_Summary.png", plot = p, width = 10, height = 8, dpi = 150)
```

``` r
knitr::include_graphics("Metadata_Summary.png")
```

<img src="Metadata_Summary.png" width="1500" />

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

    ##    BioProject       Run sample_name Age BMI    Sex Region            Platform
    ## 57 PRJNA48333 SRR040950   159247771  NA  NA female  V3-V5 454 GS FLX Titanium
    ## 58 PRJNA48333 SRR040951   160158126  NA  NA   male  V3-V5 454 GS FLX Titanium
    ## 59 PRJNA48333 SRR040952   160016515  NA  NA female  V3-V5 454 GS FLX Titanium
    ## 60 PRJNA48333 SRR040953   159450072  NA  NA   male  V3-V5 454 GS FLX Titanium
    ## 61 PRJNA48333 SRR040954   159369152  NA  NA   male  V3-V5 454 GS FLX Titanium
    ## 62 PRJNA48333 SRR040955   159369152  NA  NA   male  V3-V5 454 GS FLX Titanium
    ##    Country                 Site          public_time PCR_primers
    ## 57     USA           Mid vagina 2010-05-11T00:00:00Z          NA
    ## 58     USA        Buccal mucosa 2010-05-11T00:00:00Z          NA
    ## 59     USA               Throat 2010-05-11T00:00:00Z          NA
    ## 60     USA                Stool 2010-05-11T00:00:00Z          NA
    ## 61     USA Supragingival plaque 2010-05-11T00:00:00Z          NA
    ## 62     USA   Subgingival plaque 2010-05-11T00:00:00Z          NA
    ##                                                                title
    ## 57 Structure function and diversity of the\nhealthy human microbiome
    ## 58 Structure function and diversity of the\nhealthy human microbiome
    ## 59 Structure function and diversity of the\nhealthy human microbiome
    ## 60 Structure function and diversity of the\nhealthy human microbiome
    ## 61 Structure function and diversity of the\nhealthy human microbiome
    ## 62 Structure function and diversity of the\nhealthy human microbiome
    ##    Disease_study Status matched_samples
    ## 57            NO     NA             YES
    ## 58            NO     NA             YES
    ## 59            NO     NA             YES
    ## 60            NO     NA             YES
    ## 61            NO     NA             YES
    ## 62            NO     NA             YES

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
