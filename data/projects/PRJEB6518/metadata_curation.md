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
setwd("D:/oral_project/PRJEB6518")

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

    ##         Run AGE age_unit altitude anonymized_name Assay.Type assigned_from_geo
    ## 1 ERR531441  26       78        0        Stool314   AMPLICON                 n
    ## 2 ERR531442  36       78        0        Stool499   AMPLICON                 n
    ## 3 ERR531443  26       78        0        Stool129   AMPLICON                 n
    ## 4 ERR531444  26       78        0        Stool348   AMPLICON                 n
    ## 5 ERR531445  36       78        0        Stool437   AMPLICON                 n
    ## 6 ERR531446  36       78        0        Stool417   AMPLICON                 n
    ##   AvgSpotLen    Bases BioProject    BioSample body_habitat body_mass_index
    ## 1         98 28485870  PRJEB6518 SAMEA2592369 UBERON:feces            19.4
    ## 2         98  9782737  PRJEB6518 SAMEA2592370 UBERON:feces            21.5
    ## 3         97 24212559  PRJEB6518 SAMEA2592371 UBERON:feces            19.4
    ## 4         98  8517721  PRJEB6518 SAMEA2592372 UBERON:feces            19.4
    ## 5         99  7849894  PRJEB6518 SAMEA2592373 UBERON:feces            21.5
    ## 6         99  7643232  PRJEB6518 SAMEA2592374 UBERON:feces            21.5
    ##   body_product    body_site    Bytes   Center.Name Collection_Date
    ## 1 UBERON:feces UBERON:feces 22725729 CCME-COLORADO            2011
    ## 2 UBERON:feces UBERON:feces  7714083 CCME-COLORADO            2011
    ## 3 UBERON:feces UBERON:feces 20031948 CCME-COLORADO            2011
    ## 4 UBERON:feces UBERON:feces  6761782 CCME-COLORADO            2011
    ## 5 UBERON:feces UBERON:feces  6131123 CCME-COLORADO            2011
    ## 6 UBERON:feces UBERON:feces  5967749 CCME-COLORADO            2011
    ##            common_name Consent DATASTORE.filetype DATASTORE.provider
    ## 1 human gut metagenome  public   fastq,run.zq,sra     ena,gs,ncbi,s3
    ## 2 human gut metagenome  public   fastq,run.zq,sra     ena,gs,ncbi,s3
    ## 3 human gut metagenome  public   fastq,run.zq,sra     ena,gs,ncbi,s3
    ## 4 human gut metagenome  public   fastq,run.zq,sra     ena,gs,ncbi,s3
    ## 5 human gut metagenome  public   fastq,run.zq,sra     ena,gs,ncbi,s3
    ## 6 human gut metagenome  public   fastq,run.zq,sra     ena,gs,ncbi,s3
    ##                           DATASTORE.region Depth     Diet elevation
    ## 1 ena,gs.us-east1,ncbi.public,s3.us-east-1     0 omnivore    256.05
    ## 2 ena,gs.us-east1,ncbi.public,s3.us-east-1     0 omnivore    256.05
    ## 3 ena,gs.us-east1,ncbi.public,s3.us-east-1     0 omnivore    256.05
    ## 4 ena,gs.us-east1,ncbi.public,s3.us-east-1     0 omnivore    256.05
    ## 5 ena,gs.us-east1,ncbi.public,s3.us-east-1     0 omnivore    256.05
    ## 6 ena,gs.us-east1,ncbi.public,s3.us-east-1     0 omnivore    256.05
    ##   ENA.FIRST.PUBLIC..run.     ENA_first_public ENA.LAST.UPDATE..run.
    ## 1             2014-06-14 2014-06-14T17:03:01Z            2018-11-16
    ## 2             2014-06-14 2014-06-14T17:03:01Z            2018-11-16
    ## 3             2014-06-14 2014-06-14T17:03:01Z            2018-11-16
    ## 4             2014-06-14 2014-06-14T17:03:01Z            2018-11-16
    ## 5             2014-06-14 2014-06-14T17:03:01Z            2018-11-16
    ## 6             2014-06-14 2014-06-14T17:03:01Z            2018-11-16
    ##        ENA_last_update        env_biome                   env_feature
    ## 1 2018-03-09T11:14:32Z ENVO:urban biome ENVO:human-associated habitat
    ## 2 2018-03-09T11:14:32Z ENVO:urban biome ENVO:human-associated habitat
    ## 3 2018-03-09T11:14:32Z ENVO:urban biome ENVO:human-associated habitat
    ## 4 2018-03-09T11:14:32Z ENVO:urban biome ENVO:human-associated habitat
    ## 5 2018-03-09T11:14:32Z ENVO:urban biome ENVO:human-associated habitat
    ## 6 2018-03-09T11:14:32Z ENVO:urban biome ENVO:human-associated habitat
    ##   env_matter experiment_center..exp. experiment_design_description..exp.
    ## 1 ENVO:feces                    CCME        daily microbiota time-series
    ## 2 ENVO:feces                    CCME        daily microbiota time-series
    ## 3 ENVO:feces                    CCME        daily microbiota time-series
    ## 4 ENVO:feces                    CCME        daily microbiota time-series
    ## 5 ENVO:feces                    CCME        daily microbiota time-series
    ## 6 ENVO:feces                    CCME        daily microbiota time-series
    ##   Experiment                 experiment_title..exp.  external_id
    ## 1  ERX496644 MIT daily human microbiota time-series SAMEA2592369
    ## 2  ERX496645 MIT daily human microbiota time-series SAMEA2592370
    ## 3  ERX496646 MIT daily human microbiota time-series SAMEA2592371
    ## 4  ERX496647 MIT daily human microbiota time-series SAMEA2592372
    ## 5  ERX496648 MIT daily human microbiota time-series SAMEA2592373
    ## 6  ERX496649 MIT daily human microbiota time-series SAMEA2592374
    ##   geo_loc_name_country geo_loc_name_country_continent
    ## 1                  USA                  North America
    ## 2                  USA                  North America
    ## 3                  USA                  North America
    ## 4                  USA                  North America
    ## 5                  USA                  North America
    ## 6                  USA                  North America
    ##   geographic_location_.country_and.or_sea. host_subject_id host_taxid
    ## 1                                 USA: GAZ     2202:DonorA       9606
    ## 2                                 USA: GAZ     2202:DonorB       9606
    ## 3                                 USA: GAZ     2202:DonorA       9606
    ## 4                                 USA: GAZ     2202:DonorA       9606
    ## 5                                 USA: GAZ     2202:DonorB       9606
    ## 6                                 USA: GAZ     2202:DonorB       9606
    ##   illumina_technology INSDC_center_name   INSDC_first_public
    ## 1               GAIIx     CCME-COLORADO 2014-06-14T17:03:01Z
    ## 2               GAIIx     CCME-COLORADO 2014-06-14T17:03:01Z
    ## 3               GAIIx     CCME-COLORADO 2014-06-14T17:03:01Z
    ## 4               GAIIx     CCME-COLORADO 2014-06-14T17:03:01Z
    ## 5               GAIIx     CCME-COLORADO 2014-06-14T17:03:01Z
    ## 6               GAIIx     CCME-COLORADO 2014-06-14T17:03:01Z
    ##      INSDC_last_update INSDC_status  Instrument LATITUDE Library.Name
    ## 1 2018-03-09T11:14:32Z       public unspecified    37.09 Stool314:511
    ## 2 2018-03-09T11:14:32Z       public unspecified    37.09 Stool499:715
    ## 3 2018-03-09T11:14:32Z       public unspecified    37.09 Stool129:319
    ## 4 2018-03-09T11:14:32Z       public unspecified    37.09 Stool348:548
    ## 5 2018-03-09T11:14:32Z       public unspecified    37.09 Stool437:647
    ## 6 2018-03-09T11:14:32Z       public unspecified    37.09 Stool417:625
    ##   LibraryLayout LibrarySelection LibrarySource LONGITUDE num_sequences..exp.
    ## 1        SINGLE              PCR   METAGENOMIC    -95.71              289890
    ## 2        SINGLE              PCR   METAGENOMIC    -95.71               98906
    ## 3        SINGLE              PCR   METAGENOMIC    -95.71              248353
    ## 4        SINGLE              PCR   METAGENOMIC    -95.71               86568
    ## 5        SINGLE              PCR   METAGENOMIC    -95.71               78520
    ## 6        SINGLE              PCR   METAGENOMIC    -95.71               76721
    ##               Organism                                 pcr_primers..exp.
    ## 1 human gut metagenome FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 2 human gut metagenome FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 3 human gut metagenome FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 4 human gut metagenome FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 5 human gut metagenome FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 6 human gut metagenome FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ##   platform..exp. Platform public          ReleaseDate row_number..exp.
    ## 1       Illumina ILLUMINA      n 2014-06-14T00:00:00Z              511
    ## 2       Illumina ILLUMINA      n 2014-06-14T00:00:00Z              715
    ## 3       Illumina ILLUMINA      n 2014-06-14T00:00:00Z              319
    ## 4       Illumina ILLUMINA      n 2014-06-14T00:00:00Z              548
    ## 5       Illumina ILLUMINA      n 2014-06-14T00:00:00Z              647
    ## 6       Illumina ILLUMINA      n 2014-06-14T00:00:00Z              625
    ##   run_center..exp. run_date..exp. samp_size sample_center..exp. Sample_ID
    ## 1           CGS-GL        5/14/11    .1\\,g                 MIT   1247321
    ## 2           CGS-GL        5/14/11    .1\\,g                 MIT   1247322
    ## 3           CGS-GL        5/14/11    .1\\,g                 MIT   1247323
    ## 4           CGS-GL        5/14/11    .1\\,g                 MIT   1247324
    ## 5           CGS-GL        5/14/11    .1\\,g                 MIT   1247325
    ## 6           CGS-GL        5/14/11    .1\\,g                 MIT   1247326
    ##    Sample.Name Sample_name      scientific_name   sequencing_meth..exp. sex
    ## 1 SAMEA2592369    Stool314 human gut metagenome sequencing by synthesis  47
    ## 2 SAMEA2592370    Stool499 human gut metagenome sequencing by synthesis  47
    ## 3 SAMEA2592371    Stool129 human gut metagenome sequencing by synthesis  47
    ## 4 SAMEA2592372    Stool348 human gut metagenome sequencing by synthesis  47
    ## 5 SAMEA2592373    Stool437 human gut metagenome sequencing by synthesis  47
    ## 6 SAMEA2592374    Stool417 human gut metagenome sequencing by synthesis  47
    ##   SRA.Study study_center..exp.              Submitter_Id target_gene..exp.
    ## 1 ERP006059                MIT qiime_study_2202:Stool314          16S rRNA
    ## 2 ERP006059                MIT qiime_study_2202:Stool499          16S rRNA
    ## 3 ERP006059                MIT qiime_study_2202:Stool129          16S rRNA
    ## 4 ERP006059                MIT qiime_study_2202:Stool348          16S rRNA
    ## 5 ERP006059                MIT qiime_study_2202:Stool437          16S rRNA
    ## 6 ERP006059                MIT qiime_study_2202:Stool417          16S rRNA
    ##   target_subfragment..exp. taxon_id collection_day          create_date version
    ## 1                       V4   408170            233 2014-06-15T00:01:00Z       1
    ## 2                       V4   408170            169 2014-06-15T00:00:00Z       1
    ## 3                       V4   408170            138 2014-06-15T00:20:00Z       1
    ## 4                       V4   408170            267 2014-06-15T00:20:00Z       1
    ## 5                       V4   408170             96 2014-06-15T00:00:00Z       1
    ## 6                       V4   408170             63 2014-06-15T00:00:00Z       1
    ##   nutrition_calcium nutrition_calorie nutrition_carb nutrition_cholesterol
    ## 1                NA                NA             NA                    NA
    ## 2                NA                NA             NA                    NA
    ## 3           1820.42          1124.263       158.2823                    82
    ## 4                NA                NA             NA                    NA
    ## 5                NA                NA             NA                    NA
    ## 6                NA                NA             NA                    NA
    ##   nutrition_fat nutrition_fiber nutrition_protein nutrition_satfat
    ## 1            NA              NA                NA               NA
    ## 2            NA              NA                NA               NA
    ## 3       30.8004          21.481          40.98985         12.10026
    ## 4            NA              NA                NA               NA
    ## 5            NA              NA                NA               NA
    ## 6            NA              NA                NA               NA
    ##   nutrition_sodium nutrition_sugar
    ## 1               NA              NA
    ## 2               NA              NA
    ## 3         1609.289        25.36069
    ## 4               NA              NA
    ## 5               NA              NA
    ## 6               NA              NA

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
  sample_name = initial_meta$Sample_name,
  Age = "26",
  BMI = initial_meta$body_mass_index,
  Sex=initial_meta$sex,
  Region = initial_meta$target_subfragment..exp.,
  Platform = initial_meta$platform..exp.,
  Country = initial_meta$geo_loc_name_country,
  Site = gsub("UBERON:", "", initial_meta$body_site),
  public_time = initial_meta$ReleaseDate,
  PCR_primers = initial_meta$pcr_primers..exp.,
  stringsAsFactors = FALSE
)
```

``` r
metadata <- metadata %>%
  mutate(Sex = recode(Sex, "48" = "F", "47" = "M"))
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$title <- "Host lifestyle affects human microbiota on daily timescales"
metadata$Disease_study <- "NO"
# Assign disease status based on sample_name patterns
metadata$Status <- "HC"

#get subjects id with Saliva available
keep_subjects<-subset(metadata,Site %in% "saliva")
#remove subjects without Saliva samples
metadata<-subset(metadata,sample_name %in% keep_subjects$sample_name)
#get subjects sampling counts
sample_counts <- table(metadata$sample_name)
# Create matched_samples column, any subjects have >=2 samples will be label "YES", otherwise "NO"
metadata$matched_samples <- ifelse(metadata$sample_name %in% names(sample_counts[sample_counts >= 2]), "YES", "NO")
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 286"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 286"

# Step 5: Visualize the Metadata

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

    ##    BioProject       Run sample_name Age  BMI Sex Region Platform Country   Site
    ## 14  PRJEB6518 ERR531454   Saliva184  26 19.4   M     V4 Illumina     USA saliva
    ## 15  PRJEB6518 ERR531455    Saliva84  26 19.4   M     V4 Illumina     USA saliva
    ## 16  PRJEB6518 ERR531456   Saliva207  26 19.4   M     V4 Illumina     USA saliva
    ## 23  PRJEB6518 ERR531463    Saliva64  26 19.4   M     V4 Illumina     USA saliva
    ## 26  PRJEB6518 ERR531466   Saliva209  26 19.4   M     V4 Illumina     USA saliva
    ## 27  PRJEB6518 ERR531467   Saliva144  26 19.4   M     V4 Illumina     USA saliva
    ##             public_time                                       PCR_primers
    ## 14 2014-06-14T00:00:00Z FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 15 2014-06-14T00:00:00Z FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 16 2014-06-14T00:00:00Z FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 23 2014-06-14T00:00:00Z FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 26 2014-06-14T00:00:00Z FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 27 2014-06-14T00:00:00Z FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ##                                                          title Disease_study
    ## 14 Host lifestyle affects human microbiota on daily timescales            NO
    ## 15 Host lifestyle affects human microbiota on daily timescales            NO
    ## 16 Host lifestyle affects human microbiota on daily timescales            NO
    ## 23 Host lifestyle affects human microbiota on daily timescales            NO
    ## 26 Host lifestyle affects human microbiota on daily timescales            NO
    ## 27 Host lifestyle affects human microbiota on daily timescales            NO
    ##    Status matched_samples
    ## 14     HC              NO
    ## 15     HC              NO
    ## 16     HC              NO
    ## 23     HC              NO
    ## 26     HC              NO
    ## 27     HC              NO

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
