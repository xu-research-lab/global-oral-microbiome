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
setwd("D:/metadata/PRJEB39420")

# Read the initial metadata files
initial_meta<-read.csv("SraRunTable.csv",header = T)

# Combine the two metadata files
#meta <- rbind(initial_meta, initial_meta1)

#Remove entries non-oral samples
#initial_meta <- subset(initial_meta, Organism == "human oral metagenome")

# Check the number of samples containing "HT" in Sample.Name (e.g., healthy controls)
#nrow(subset(meta, grepl("HT", Sample.Name)))

# Display the first few rows of the metadata
head(initial_meta)
```

    ##          Run Assay.Type Barcode..exp. BioProject    BioSample center_name..exp.
    ## 1 ERR4345071      OTHER  AGCCTTCGTCGC PRJEB39420 SAMEA7062799            UCSDMI
    ## 2 ERR4345072      OTHER  CGTATAAATGCG PRJEB39420 SAMEA7062800            UCSDMI
    ## 3 ERR4345073      OTHER  TGACTAATGGCC PRJEB39420 SAMEA7062801            UCSDMI
    ## 4 ERR4345074      OTHER  GTGGAGTCTCAT PRJEB39420 SAMEA7062802            UCSDMI
    ## 5 ERR4345075      OTHER  TGATGTGCTAAG PRJEB39420 SAMEA7062803            UCSDMI
    ## 6 ERR4345076      OTHER  TGTGCACGCCAT PRJEB39420 SAMEA7062804            UCSDMI
    ##                                          Center.Name center_project_name..exp.
    ## 1 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT                Megan_Doty
    ## 2 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT                Megan_Doty
    ## 3 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT                Megan_Doty
    ## 4 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT                Megan_Doty
    ## 5 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT                Megan_Doty
    ## 6 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT                Megan_Doty
    ##   Consent ENA_FIRST_PUBLIC     ENA_first_public ENA.LAST.UPDATE..run.
    ## 1  public        2020/7/16 2020-07-16T17:05:32Z             2020/7/16
    ## 2  public        2020/7/16 2020-07-16T17:05:32Z             2020/7/16
    ## 3  public        2020/7/16 2020-07-16T17:05:32Z             2020/7/16
    ## 4  public        2020/7/16 2020-07-16T17:05:32Z             2020/7/16
    ## 5  public        2020/7/16 2020-07-16T17:05:32Z             2020/7/16
    ## 6  public        2020/7/16 2020-07-16T17:05:32Z             2020/7/16
    ##        ENA_last_update Experiment  external_id extraction_robot..exp.
    ## 1 2020-07-16T09:42:09Z ERX4289735 SAMEA7062799               HOWE_KF3
    ## 2 2020-07-16T09:42:09Z ERX4289736 SAMEA7062800               HOWE_KF3
    ## 3 2020-07-16T09:42:09Z ERX4289737 SAMEA7062801               HOWE_KF3
    ## 4 2020-07-16T09:42:09Z ERX4289738 SAMEA7062802               HOWE_KF3
    ## 5 2020-07-16T09:42:09Z ERX4289739 SAMEA7062803               HOWE_KF3
    ## 6 2020-07-16T09:42:09Z ERX4289740 SAMEA7062804               HOWE_KF3
    ##   extractionkit_lot..exp. host_subject_id INSDC_center_name
    ## 1               154050095       BLANK1.1A            UCSDMI
    ## 2               154050095       BLANK1.1B            UCSDMI
    ## 3               154050095       BLANK1.1C            UCSDMI
    ## 4               154050095       BLANK1.1D            UCSDMI
    ## 5               154050095       BLANK1.1E            UCSDMI
    ## 6               154050095       BLANK1.1F            UCSDMI
    ##     INSDC_first_public    INSDC_last_update INSDC_status          Instrument
    ## 1 2020-07-16T17:05:32Z 2020-07-16T09:42:09Z       public Illumina HiSeq 2500
    ## 2 2020-07-16T17:05:32Z 2020-07-16T09:42:09Z       public Illumina HiSeq 2500
    ## 3 2020-07-16T17:05:32Z 2020-07-16T09:42:09Z       public Illumina HiSeq 2500
    ## 4 2020-07-16T17:05:32Z 2020-07-16T09:42:09Z       public Illumina HiSeq 2500
    ## 5 2020-07-16T17:05:32Z 2020-07-16T09:42:09Z       public Illumina HiSeq 2500
    ## 6 2020-07-16T17:05:32Z 2020-07-16T09:42:09Z       public Illumina HiSeq 2500
    ##   lane..exp.    Library.Name LibraryLayout LibrarySelection LibrarySource
    ## 1          1 11712.BLANK1.1A        SINGLE              PCR   METAGENOMIC
    ## 2          1 11712.BLANK1.1B        SINGLE              PCR   METAGENOMIC
    ## 3          1 11712.BLANK1.1C        SINGLE              PCR   METAGENOMIC
    ## 4          1 11712.BLANK1.1D        SINGLE              PCR   METAGENOMIC
    ## 5          1 11712.BLANK1.1E        SINGLE              PCR   METAGENOMIC
    ## 6          1 11712.BLANK1.1F        SINGLE              PCR   METAGENOMIC
    ##   mastermix_lot..exp.   Organism orig_name..exp.
    ## 1              443912 metagenome       BLANK1.1A
    ## 2              443912 metagenome       BLANK1.1B
    ## 3              443912 metagenome       BLANK1.1C
    ## 4              443912 metagenome       BLANK1.1D
    ## 5              443912 metagenome       BLANK1.1E
    ## 6              443912 metagenome       BLANK1.1F
    ##                                   pcr_primers..exp. Platform plating..exp.
    ## 1 FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT ILLUMINA           LDG
    ## 2 FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT ILLUMINA           LDG
    ## 3 FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT ILLUMINA           LDG
    ## 4 FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT ILLUMINA           LDG
    ## 5 FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT ILLUMINA           LDG
    ## 6 FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT ILLUMINA           LDG
    ##   primer_date..exp.          primer..exp. primer_plate..exp.
    ## 1           1/19/17 GTGTGYCAGCMGCCGCGGTAA                  1
    ## 2           1/19/17 GTGTGYCAGCMGCCGCGGTAA                  1
    ## 3           1/19/17 GTGTGYCAGCMGCCGCGGTAA                  1
    ## 4           1/19/17 GTGTGYCAGCMGCCGCGGTAA                  1
    ## 5           1/19/17 GTGTGYCAGCMGCCGCGGTAA                  1
    ## 6           1/19/17 GTGTGYCAGCMGCCGCGGTAA                  1
    ##   processing_robot..exp. project_name..exp.          ReleaseDate run_date..exp.
    ## 1                   RIKE         Megan_Doty 2020-07-17T00:00:00Z      2003/7/17
    ## 2                   RIKE         Megan_Doty 2020-07-17T00:00:00Z      2003/7/17
    ## 3                   RIKE         Megan_Doty 2020-07-17T00:00:00Z      2003/7/17
    ## 4                   RIKE         Megan_Doty 2020-07-17T00:00:00Z      2003/7/17
    ## 5                   RIKE         Megan_Doty 2020-07-17T00:00:00Z      2003/7/17
    ## 6                   RIKE         Megan_Doty 2020-07-17T00:00:00Z      2003/7/17
    ##       run_prefix..exp.                          runid..exp.  Sample.Name
    ## 1 Megan_Doty_1_S1_L001 170307_D00611_0447_BHFFTFBCXY_Knight SAMEA7062799
    ## 2 Megan_Doty_1_S1_L001 170307_D00611_0447_BHFFTFBCXY_Knight SAMEA7062800
    ## 3 Megan_Doty_1_S1_L001 170307_D00611_0447_BHFFTFBCXY_Knight SAMEA7062801
    ## 4 Megan_Doty_1_S1_L001 170307_D00611_0447_BHFFTFBCXY_Knight SAMEA7062802
    ## 5 Megan_Doty_1_S1_L001 170307_D00611_0447_BHFFTFBCXY_Knight SAMEA7062803
    ## 6 Megan_Doty_1_S1_L001 170307_D00611_0447_BHFFTFBCXY_Knight SAMEA7062804
    ##                       Sample_name sample_plate..exp. Scientific_Name
    ## 1 qiita_sid_11712:11712.BLANK1.1A       Megan_Doty_1      metagenome
    ## 2 qiita_sid_11712:11712.BLANK1.1B       Megan_Doty_1      metagenome
    ## 3 qiita_sid_11712:11712.BLANK1.1C       Megan_Doty_1      metagenome
    ## 4 qiita_sid_11712:11712.BLANK1.1D       Megan_Doty_1      metagenome
    ## 5 qiita_sid_11712:11712.BLANK1.1E       Megan_Doty_1      metagenome
    ## 6 qiita_sid_11712:11712.BLANK1.1F       Megan_Doty_1      metagenome
    ##     sequencing_meth..exp. SRA.Study                    Submitter_Id
    ## 1 Sequencing by synthesis ERP122936 qiita_sid_11712:11712.BLANK1.1A
    ## 2 Sequencing by synthesis ERP122936 qiita_sid_11712:11712.BLANK1.1B
    ## 3 Sequencing by synthesis ERP122936 qiita_sid_11712:11712.BLANK1.1C
    ## 4 Sequencing by synthesis ERP122936 qiita_sid_11712:11712.BLANK1.1D
    ## 5 Sequencing by synthesis ERP122936 qiita_sid_11712:11712.BLANK1.1E
    ## 6 Sequencing by synthesis ERP122936 qiita_sid_11712:11712.BLANK1.1F
    ##   target_gene..exp. target_subfragment tm1000_8_tool..exp. tm300_8_tool..exp.
    ## 1          16S rRNA                 V4             108379Z            311318B
    ## 2          16S rRNA                 V4             108379Z            311318B
    ## 3          16S rRNA                 V4             108379Z            311318B
    ## 4          16S rRNA                 V4             108379Z            311318B
    ## 5          16S rRNA                 V4             108379Z            311318B
    ## 6          16S rRNA                 V4             108379Z            311318B
    ##   tm50_8_tool..exp. water_lot..exp. AvgSpotLen  Bases  Bytes DATASTORE.filetype
    ## 1           311426B         1798402        150   2250  29663         run.zq,sra
    ## 2           311426B         1798402        150    750  28938         run.zq,sra
    ## 3           311426B         1798402        150   2700  29820         run.zq,sra
    ## 4           311426B         1798402        150   1500  29368         run.zq,sra
    ## 5           311426B         1798402        150   2100  29892         run.zq,sra
    ## 6           311426B         1798402        150 607950 331978         run.zq,sra
    ##   DATASTORE.provider                     DATASTORE.region          create_date
    ## 1         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 2020-07-31T15:45:00Z
    ## 2         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 2020-07-31T15:44:00Z
    ## 3         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 2020-07-19T18:21:00Z
    ## 4         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 2020-07-19T18:23:00Z
    ## 5         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 2020-07-19T18:24:00Z
    ## 6         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 2020-07-19T18:24:00Z
    ##   version birth_length birth_weight Collection_time dc_wt dhm_wk1 formula_wk1
    ## 1       1           NA           NA                    NA      NA          NA
    ## 2       1           NA           NA                    NA      NA          NA
    ## 3       1           NA           NA                    NA      NA          NA
    ## 4       1           NA           NA                    NA      NA          NA
    ## 5       1           NA           NA                    NA      NA          NA
    ## 6       1           NA           NA                    NA      NA          NA
    ##   Host_age host_body_mass_index hr_1stfeed maternal_age mbm_wk1 ml_dhm_wk1
    ## 1       NA                   NA         NA           NA                 NA
    ## 2       NA                   NA         NA           NA                 NA
    ## 3       NA                   NA         NA           NA                 NA
    ## 4       NA                   NA         NA           NA                 NA
    ## 5       NA                   NA         NA           NA                 NA
    ## 6       NA                   NA         NA           NA                 NA
    ##   ml_formula_wk1 ml_mbm_wk1 nadir_wt pct_bw_lost pct_dc_wt_abovebw time_only_bm
    ## 1             NA                  NA          NA                NA           NA
    ## 2             NA                  NA          NA                NA           NA
    ## 3             NA                  NA          NA                NA           NA
    ## 4             NA                  NA          NA                NA           NA
    ## 5             NA                  NA          NA                NA           NA
    ## 6             NA                  NA          NA                NA           NA
    ##   hm_only_wk1 pct_formula_wk1 LOS ms_plate_num dol_1stformula avg_hero_score
    ## 1          NA              NA  NA           NA             NA             NA
    ## 2          NA              NA  NA           NA             NA             NA
    ## 3          NA              NA  NA           NA             NA             NA
    ## 4          NA              NA  NA           NA             NA             NA
    ## 5          NA              NA  NA           NA             NA             NA
    ## 6          NA              NA  NA           NA             NA             NA
    ##   high_hero_score
    ## 1              NA
    ## 2              NA
    ## 3              NA
    ## 4              NA
    ## 5              NA
    ## 6              NA

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
  sample_name = initial_meta$Sample.Name,
  Age = initial_meta$Host_age / 365,
  BMI = initial_meta$host_body_mass_index,
  Sex=NA,
  Region = "V4",
  Platform = initial_meta$Instrument,
  Country = NA,
  Site = "oral",
  public_time = initial_meta$ReleaseDate,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$PCR_primers <- "515F-806R"
metadata$title <- "Multiomic Analyses of Nascent Preterm Infant Microbiomes Differentiation Suggest Opportunities for Targeted Intervention"
metadata$Disease_study <- "YES"
# Assign disease status based on sample_name patterns
metadata$Status<-"preterm infant"
keep_subjects <- unique(subset(metadata, Site == "oral")$sample_name)
metadata <- subset(metadata, sample_name %in% keep_subjects)

sample_counts <- table(metadata$sample_name)

metadata$matched_samples <-
  ifelse(metadata$sample_name %in% names(sample_counts[sample_counts >= 2]),
         "YES", "NO")

cat("Subjects retained:", length(unique(metadata$sample_name)), "\n")
```

    ## Subjects retained: 1385

``` r
cat("Matched samples distribution:\n")
```

    ## Matched samples distribution:

``` r
table(metadata$matched_samples)
```

    ## 
    ##   NO 
    ## 1385

``` r
cat("Sample distribution by Site:\n")
```

    ## Sample distribution by Site:

``` r
table(metadata$Site)
```

    ## 
    ## oral 
    ## 1385

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 1385"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 1385"

# Step 5: Visualize the Metadata

``` r
knitr::include_graphics("metadata_summary.png")
```

![](metadata_summary.png)<!-- -->

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

    ##   BioProject        Run  sample_name Age BMI Sex Region            Platform
    ## 1 PRJEB39420 ERR4345071 SAMEA7062799  NA  NA  NA     V4 Illumina HiSeq 2500
    ## 2 PRJEB39420 ERR4345072 SAMEA7062800  NA  NA  NA     V4 Illumina HiSeq 2500
    ## 3 PRJEB39420 ERR4345073 SAMEA7062801  NA  NA  NA     V4 Illumina HiSeq 2500
    ## 4 PRJEB39420 ERR4345074 SAMEA7062802  NA  NA  NA     V4 Illumina HiSeq 2500
    ## 5 PRJEB39420 ERR4345075 SAMEA7062803  NA  NA  NA     V4 Illumina HiSeq 2500
    ## 6 PRJEB39420 ERR4345076 SAMEA7062804  NA  NA  NA     V4 Illumina HiSeq 2500
    ##   Country Site          public_time PCR_primers
    ## 1      NA oral 2020-07-17T00:00:00Z   515F-806R
    ## 2      NA oral 2020-07-17T00:00:00Z   515F-806R
    ## 3      NA oral 2020-07-17T00:00:00Z   515F-806R
    ## 4      NA oral 2020-07-17T00:00:00Z   515F-806R
    ## 5      NA oral 2020-07-17T00:00:00Z   515F-806R
    ## 6      NA oral 2020-07-17T00:00:00Z   515F-806R
    ##                                                                                                                      title
    ## 1 Multiomic Analyses of Nascent Preterm Infant Microbiomes Differentiation Suggest Opportunities for Targeted Intervention
    ## 2 Multiomic Analyses of Nascent Preterm Infant Microbiomes Differentiation Suggest Opportunities for Targeted Intervention
    ## 3 Multiomic Analyses of Nascent Preterm Infant Microbiomes Differentiation Suggest Opportunities for Targeted Intervention
    ## 4 Multiomic Analyses of Nascent Preterm Infant Microbiomes Differentiation Suggest Opportunities for Targeted Intervention
    ## 5 Multiomic Analyses of Nascent Preterm Infant Microbiomes Differentiation Suggest Opportunities for Targeted Intervention
    ## 6 Multiomic Analyses of Nascent Preterm Infant Microbiomes Differentiation Suggest Opportunities for Targeted Intervention
    ##   Disease_study         Status matched_samples
    ## 1           YES preterm infant              NO
    ## 2           YES preterm infant              NO
    ## 3           YES preterm infant              NO
    ## 4           YES preterm infant              NO
    ## 5           YES preterm infant              NO
    ## 6           YES preterm infant              NO

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
