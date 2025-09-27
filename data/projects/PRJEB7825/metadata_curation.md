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
setwd("D:/oral_project/PRJEB7825")

# Read the initial metadata files
initial_meta<-read.csv("SraRunTable.csv",header = T)

# Combine the two metadata files
#meta <- rbind(initial_meta, initial_meta1)

#Remove entries non-oral samples
#initial_meta <- subset(initial_meta,isolate %in% c("saliva","tongue coat"))
initial_meta <- initial_meta[initial_meta$Organism != "human gut metagenome", ]
initial_meta <- initial_meta[initial_meta$Organism != "human skin metagenome", ]
# Check the number of samples containing "HT" in Sample.Name (e.g., healthy controls)
#nrow(subset(meta, grepl("HT", Sample.Name)))

# Display the first few rows of the metadata
head(initial_meta)
```

    ##          Run age_unit altitude anonymized_name Assay.Type assigned_from_geo
    ## 40 ERR687975    years        0          MG.O10      OTHER             False
    ## 41 ERR687976    years        0          MG.O11      OTHER             False
    ## 42 ERR687977    years        0          MG.O12      OTHER             False
    ## 43 ERR687978    years        0          MG.O13      OTHER             False
    ## 44 ERR687979    years        0          MG.O14      OTHER             False
    ## 45 ERR687980    years        0          MG.O16      OTHER             False
    ##    AvgSpotLen barcodesequence..exp.    Bases BioProject    BioSample
    ## 40        127          AGTTCAGACGCT  1825000  PRJEB7825 SAMEA3139776
    ## 41        126          ATCACGTAGCGG  1872743  PRJEB7825 SAMEA3139777
    ## 42        128          ATCTGGTGCTAT  4107688  PRJEB7825 SAMEA3139778
    ## 43        128          ATGGCGTGCACA  4795309  PRJEB7825 SAMEA3139779
    ## 44        128          CAAGTGAGAGAG  7549485  PRJEB7825 SAMEA3139780
    ## 45        130          CACTCTGATTAG 14038609  PRJEB7825 SAMEA3139781
    ##          body_habitat  body_product    body_site    Bytes center_name..exp.
    ## 40 UBERON:oral cavity UBERON:saliva UBERON:mouth  1510817              CCME
    ## 41 UBERON:oral cavity UBERON:saliva UBERON:mouth  1548981              CCME
    ## 42 UBERON:oral cavity UBERON:saliva UBERON:mouth  3308666              CCME
    ## 43 UBERON:oral cavity UBERON:saliva UBERON:mouth  3844150              CCME
    ## 44 UBERON:oral cavity UBERON:saliva UBERON:mouth  5854868              CCME
    ## 45 UBERON:oral cavity UBERON:saliva UBERON:mouth 11064456              CCME
    ##                                           Center.Name
    ## 40 COLORADO CENTER FOR MICROBIAL ECOLOGY, UNIVERSITY 
    ## 41 COLORADO CENTER FOR MICROBIAL ECOLOGY, UNIVERSITY 
    ## 42 COLORADO CENTER FOR MICROBIAL ECOLOGY, UNIVERSITY 
    ## 43 COLORADO CENTER FOR MICROBIAL ECOLOGY, UNIVERSITY 
    ## 44 COLORADO CENTER FOR MICROBIAL ECOLOGY, UNIVERSITY 
    ## 45 COLORADO CENTER FOR MICROBIAL ECOLOGY, UNIVERSITY 
    ##                center_project_name..exp. Collection_Date collection_timestamp
    ## 40 Dominguez_Perisipiwei_skin_oral_fecal         4/25/09  2009-04-25 00:00:00
    ## 41 Dominguez_Perisipiwei_skin_oral_fecal         4/25/09  2009-04-25 00:00:00
    ## 42 Dominguez_Perisipiwei_skin_oral_fecal         4/25/09  2009-04-25 00:00:00
    ## 43 Dominguez_Perisipiwei_skin_oral_fecal         4/25/09  2009-04-25 00:00:00
    ## 44 Dominguez_Perisipiwei_skin_oral_fecal         4/25/09  2009-04-25 00:00:00
    ## 45 Dominguez_Perisipiwei_skin_oral_fecal         4/25/09  2009-04-25 00:00:00
    ##              common_name Consent        COUNTRY DATASTORE.filetype
    ## 40 human oral metagenome  public Venezuela: GAZ   fastq,run.zq,sra
    ## 41 human oral metagenome  public Venezuela: GAZ   fastq,run.zq,sra
    ## 42 human oral metagenome  public Venezuela: GAZ   fastq,run.zq,sra
    ## 43 human oral metagenome  public Venezuela: GAZ   fastq,run.zq,sra
    ## 44 human oral metagenome  public Venezuela: GAZ   fastq,run.zq,sra
    ## 45 human oral metagenome  public Venezuela: GAZ   fastq,run.zq,sra
    ##    DATASTORE.provider                         DATASTORE.region Depth
    ## 40     ena,gs,ncbi,s3 ena,gs.us-east1,ncbi.public,s3.us-east-1     0
    ## 41     ena,gs,ncbi,s3 ena,gs.us-east1,ncbi.public,s3.us-east-1     0
    ## 42     ena,gs,ncbi,s3 ena,gs.us-east1,ncbi.public,s3.us-east-1     0
    ## 43     ena,gs,ncbi,s3 ena,gs.us-east1,ncbi.public,s3.us-east-1     0
    ## 44     ena,gs,ncbi,s3 ena,gs.us-east1,ncbi.public,s3.us-east-1     0
    ## 45     ena,gs,ncbi,s3 ena,gs.us-east1,ncbi.public,s3.us-east-1     0
    ##    dna_extracted elevation emp_status..exp. ENA.FIRST.PUBLIC..run.
    ## 40          True  89.84744          NOT_EMP             2015-04-16
    ## 41          True  89.84744          NOT_EMP             2015-04-16
    ## 42          True  89.84744          NOT_EMP             2015-04-16
    ## 43          True  89.84744          NOT_EMP             2015-04-16
    ## 44          True  89.84744          NOT_EMP             2015-04-16
    ## 45          True  89.84744          NOT_EMP             2015-04-16
    ##        ENA_first_public ENA.LAST.UPDATE..run.      ENA_last_update
    ## 40 2015-04-16T17:02:02Z            2018-11-16 2018-03-09T12:10:32Z
    ## 41 2015-04-16T17:02:02Z            2018-11-16 2018-03-09T12:10:32Z
    ## 42 2015-04-16T17:02:02Z            2018-11-16 2018-03-09T12:10:32Z
    ## 43 2015-04-16T17:02:02Z            2018-11-16 2018-03-09T12:10:32Z
    ## 44 2015-04-16T17:02:02Z            2018-11-16 2018-03-09T12:10:32Z
    ## 45 2015-04-16T17:02:02Z            2018-11-16 2018-03-09T12:10:32Z
    ##             env_biome                   env_feature  env_matter
    ## 40 ENVO:village biome ENVO:human-associated habitat ENVO:saliva
    ## 41 ENVO:village biome ENVO:human-associated habitat ENVO:saliva
    ## 42 ENVO:village biome ENVO:human-associated habitat ENVO:saliva
    ## 43 ENVO:village biome ENVO:human-associated habitat ENVO:saliva
    ## 44 ENVO:village biome ENVO:human-associated habitat ENVO:saliva
    ## 45 ENVO:village biome ENVO:human-associated habitat ENVO:saliva
    ##                                                           experiment_design_description..exp.
    ## 40 samples of, skin, saliva and feces from Perisipiwei as part of western acculturation study
    ## 41 samples of, skin, saliva and feces from Perisipiwei as part of western acculturation study
    ## 42 samples of, skin, saliva and feces from Perisipiwei as part of western acculturation study
    ## 43 samples of, skin, saliva and feces from Perisipiwei as part of western acculturation study
    ## 44 samples of, skin, saliva and feces from Perisipiwei as part of western acculturation study
    ## 45 samples of, skin, saliva and feces from Perisipiwei as part of western acculturation study
    ##    Experiment  external_id
    ## 40  ERX632969 SAMEA3139776
    ## 41  ERX632970 SAMEA3139777
    ## 42  ERX632971 SAMEA3139778
    ## 43  ERX632972 SAMEA3139779
    ## 44  ERX632973 SAMEA3139780
    ## 45  ERX632974 SAMEA3139781
    ##                                                  file_path..exp.
    ## 40 /projects/qiita_data/working_dir/6jcIl1/10052.MG.O10.fastq.gz
    ## 41 /projects/qiita_data/working_dir/6jcIl1/10052.MG.O11.fastq.gz
    ## 42 /projects/qiita_data/working_dir/6jcIl1/10052.MG.O12.fastq.gz
    ## 43 /projects/qiita_data/working_dir/6jcIl1/10052.MG.O13.fastq.gz
    ## 44 /projects/qiita_data/working_dir/6jcIl1/10052.MG.O14.fastq.gz
    ## 45 /projects/qiita_data/working_dir/6jcIl1/10052.MG.O16.fastq.gz
    ##    file_type..exp. host_common_name host_subject_id host_taxid
    ## 40           fastq            human               9       9606
    ## 41           fastq            human              10       9606
    ## 42           fastq            human              11       9606
    ## 43           fastq            human              12       9606
    ## 44           fastq            human              13       9606
    ## 45           fastq            human              14       9606
    ##    illumina_technology..exp. INSDC_center_name   INSDC_first_public
    ## 40                     GAIIx     CCME-COLORADO 2015-04-16T17:02:02Z
    ## 41                     GAIIx     CCME-COLORADO 2015-04-16T17:02:02Z
    ## 42                     GAIIx     CCME-COLORADO 2015-04-16T17:02:02Z
    ## 43                     GAIIx     CCME-COLORADO 2015-04-16T17:02:02Z
    ## 44                     GAIIx     CCME-COLORADO 2015-04-16T17:02:02Z
    ## 45                     GAIIx     CCME-COLORADO 2015-04-16T17:02:02Z
    ##       INSDC_last_update INSDC_status  Instrument key_seq..exp. LATITUDE
    ## 40 2018-03-09T12:10:32Z       public unspecified           NaN  6.42375
    ## 41 2018-03-09T12:10:32Z       public unspecified           NaN  6.42375
    ## 42 2018-03-09T12:10:32Z       public unspecified           NaN  6.42375
    ## 43 2018-03-09T12:10:32Z       public unspecified           NaN  6.42375
    ## 44 2018-03-09T12:10:32Z       public unspecified           NaN  6.42375
    ## 45 2018-03-09T12:10:32Z       public unspecified           NaN  6.42375
    ##    library_construction_protocol..exp. Library.Name LibraryLayout
    ## 40                      Illumina GAIIx 10052.MG.O10        SINGLE
    ## 41                      Illumina GAIIx 10052.MG.O11        SINGLE
    ## 42                      Illumina GAIIx 10052.MG.O12        SINGLE
    ## 43                      Illumina GAIIx 10052.MG.O13        SINGLE
    ## 44                      Illumina GAIIx 10052.MG.O14        SINGLE
    ## 45                      Illumina GAIIx 10052.MG.O16        SINGLE
    ##    LibrarySelection LibrarySource linker..exp. linkerprimersequence..exp.
    ## 40      unspecified   METAGENOMIC           TC       AGAGTTTGATCCTGGCTCAG
    ## 41      unspecified   METAGENOMIC           TC       AGAGTTTGATCCTGGCTCAG
    ## 42      unspecified   METAGENOMIC           TC       AGAGTTTGATCCTGGCTCAG
    ## 43      unspecified   METAGENOMIC           TC       AGAGTTTGATCCTGGCTCAG
    ## 44      unspecified   METAGENOMIC           TC       AGAGTTTGATCCTGGCTCAG
    ## 45      unspecified   METAGENOMIC           TC       AGAGTTTGATCCTGGCTCAG
    ##    LONGITUDE              Organism
    ## 40 -66.58973 human oral metagenome
    ## 41 -66.58973 human oral metagenome
    ## 42 -66.58973 human oral metagenome
    ## 43 -66.58973 human oral metagenome
    ## 44 -66.58973 human oral metagenome
    ## 45 -66.58973 human oral metagenome
    ##                                    pcr_primers..exp. physical_specimen_location
    ## 40 FWD:AGAGTTTGATCCTGGCTCAG; REV:TGCTGCCTCCCGTAGGAGT                        UPR
    ## 41 FWD:AGAGTTTGATCCTGGCTCAG; REV:TGCTGCCTCCCGTAGGAGT                        UPR
    ## 42 FWD:AGAGTTTGATCCTGGCTCAG; REV:TGCTGCCTCCCGTAGGAGT                        UPR
    ## 43 FWD:AGAGTTTGATCCTGGCTCAG; REV:TGCTGCCTCCCGTAGGAGT                        UPR
    ## 44 FWD:AGAGTTTGATCCTGGCTCAG; REV:TGCTGCCTCCCGTAGGAGT                        UPR
    ## 45 FWD:AGAGTTTGATCCTGGCTCAG; REV:TGCTGCCTCCCGTAGGAGT                        UPR
    ##    physical_specimen_remaining platform..exp. Platform public region..exp.
    ## 40                        True       Illumina ILLUMINA   True          NaN
    ## 41                        True       Illumina ILLUMINA   True          NaN
    ## 42                        True       Illumina ILLUMINA   True          NaN
    ## 43                        True       Illumina ILLUMINA   True          NaN
    ## 44                        True       Illumina ILLUMINA   True          NaN
    ## 45                        True       Illumina ILLUMINA   True          NaN
    ##             ReleaseDate required_sample_info_status run_center..exp.
    ## 40 2015-04-17T00:00:00Z                   completed             CCME
    ## 41 2015-04-17T00:00:00Z                   completed             CCME
    ## 42 2015-04-17T00:00:00Z                   completed             CCME
    ## 43 2015-04-17T00:00:00Z                   completed             CCME
    ## 44 2015-04-17T00:00:00Z                   completed             CCME
    ## 45 2015-04-17T00:00:00Z                   completed             CCME
    ##    run_date..exp.          create_date version run_prefix..exp. samp_size..exp.
    ## 40        5/14/11 2017-08-03T12:06:00Z       1     s_4_sequence            .1,g
    ## 41        5/14/11 2017-08-03T12:12:00Z       1     s_4_sequence            .1,g
    ## 42        5/14/11 2017-08-03T12:13:00Z       1     s_4_sequence            .1,g
    ## 43        5/14/11 2017-08-03T12:14:00Z       1     s_4_sequence            .1,g
    ## 44        5/14/11 2016-02-10T22:40:00Z       1     s_4_sequence            .1,g
    ## 45        5/14/11 2017-08-03T12:15:00Z       1     s_4_sequence            .1,g
    ##          sample_center..exp.  Sample.Name                Sample_name
    ## 40 University of Puerto Rico SAMEA3139776 qiime_ppdid_4:10052.MG.O10
    ## 41 University of Puerto Rico SAMEA3139777 qiime_ppdid_4:10052.MG.O11
    ## 42 University of Puerto Rico SAMEA3139778 qiime_ppdid_4:10052.MG.O12
    ## 43 University of Puerto Rico SAMEA3139779 qiime_ppdid_4:10052.MG.O13
    ## 44 University of Puerto Rico SAMEA3139780 qiime_ppdid_4:10052.MG.O14
    ## 45 University of Puerto Rico SAMEA3139781 qiime_ppdid_4:10052.MG.O16
    ##    sample_type       scientific_name   sequencing_meth..exp. SRA.Study
    ## 40      Saliva human oral metagenome sequencing by synthesis ERP008799
    ## 41      Saliva human oral metagenome sequencing by synthesis ERP008799
    ## 42      Saliva human oral metagenome sequencing by synthesis ERP008799
    ## 43      Saliva human oral metagenome sequencing by synthesis ERP008799
    ## 44      Saliva human oral metagenome sequencing by synthesis ERP008799
    ## 45      Saliva human oral metagenome sequencing by synthesis ERP008799
    ##           study_center..exp.               Submitter_Id target_gene..exp.
    ## 40 University of Puerto Rico qiime_ppdid_4:10052.MG.O10          16S rRNA
    ## 41 University of Puerto Rico qiime_ppdid_4:10052.MG.O11          16S rRNA
    ## 42 University of Puerto Rico qiime_ppdid_4:10052.MG.O12          16S rRNA
    ## 43 University of Puerto Rico qiime_ppdid_4:10052.MG.O13          16S rRNA
    ## 44 University of Puerto Rico qiime_ppdid_4:10052.MG.O14          16S rRNA
    ## 45 University of Puerto Rico qiime_ppdid_4:10052.MG.O16          16S rRNA
    ##    target_subfragment..exp.                                         diet_type
    ## 40                       V2 hunter gatherer diet, cassava crops, local fruits
    ## 41                       V2 hunter gatherer diet, cassava crops, local fruits
    ## 42                       V2 hunter gatherer diet, cassava crops, local fruits
    ## 43                       V2 hunter gatherer diet, cassava crops, local fruits
    ## 44                       V2 hunter gatherer diet, cassava crops, local fruits
    ## 45                       V2 hunter gatherer diet, cassava crops, local fruits
    ##    ethnic_group village_shabono_name AGE    sex
    ## 40     Yanomami        Yanomani 2008  40   male
    ## 41     Yanomami        Yanomani 2008  11   male
    ## 42     Yanomami        Yanomani 2008  45 female
    ## 43     Yanomami        Yanomani 2008  48 female
    ## 44     Yanomami        Yanomani 2008  47 female
    ## 45     Yanomami        Yanomani 2008  45 female

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
  sample_name = initial_meta$Library.Name,
  Age = initial_meta$AGE,
  BMI = NA,
  Sex=initial_meta$sex,
  Region = initial_meta$target_subfragment..exp.,
  Platform = initial_meta$platform..exp.,
  Country = initial_meta$COUNTRY,
  Site = initial_meta$Organism,
  public_time = "2015-04-17",
  PCR_primers = "27F-338R",
  stringsAsFactors = FALSE
)
```

``` r
metadata <- metadata %>%
  mutate(Sex = recode(Sex, "female" = "F", "male" = "M"))
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$title <- "The microbiome of uncontacted Amerindians"
metadata$Disease_study <- "NO"
# Assign disease status based on sample_name patterns
metadata$Status <- "HC"

#extract subject id by remove prefix of sample_name
metadata$sample_name <- gsub("SAMEA", "", metadata$sample_name)

# Check initial data distribution
cat("Initial sample distribution:\n")
```

    ## Initial sample distribution:

``` r
table(metadata$sample_name)
```

    ## 
    ## 10052.MG.O10 10052.MG.O11 10052.MG.O12 10052.MG.O13 10052.MG.O14 10052.MG.O16 
    ##            1            1            1            1            1            1 
    ## 10052.MG.O17 10052.MG.O18 10052.MG.O19  10052.MG.O2 10052.MG.O20 10052.MG.O21 
    ##            1            1            1            1            1            1 
    ## 10052.MG.O22 10052.MG.O23 10052.MG.O24 10052.MG.O25  10052.MG.O3  10052.MG.O4 
    ##            1            1            1            1            1            1 
    ## 10052.MG.O49  10052.MG.O5 10052.MG.O50 10052.MG.O51 10052.MG.O54  10052.MG.O6 
    ##            1            1            1            1            1            1 
    ##  10052.MG.O7  10052.MG.O8  10052.MG.O9 
    ##            1            1            1

``` r
table(metadata$Site)
```

    ## 
    ## human oral metagenome 
    ##                    27

``` r
# Calculate sampling frequency per subject
sample_counts <- table(metadata$sample_name)
cat("\nNumber of samples per subject:\n")
```

    ## 
    ## Number of samples per subject:

``` r
print(sample_counts)
```

    ## 
    ## 10052.MG.O10 10052.MG.O11 10052.MG.O12 10052.MG.O13 10052.MG.O14 10052.MG.O16 
    ##            1            1            1            1            1            1 
    ## 10052.MG.O17 10052.MG.O18 10052.MG.O19  10052.MG.O2 10052.MG.O20 10052.MG.O21 
    ##            1            1            1            1            1            1 
    ## 10052.MG.O22 10052.MG.O23 10052.MG.O24 10052.MG.O25  10052.MG.O3  10052.MG.O4 
    ##            1            1            1            1            1            1 
    ## 10052.MG.O49  10052.MG.O5 10052.MG.O50 10052.MG.O51 10052.MG.O54  10052.MG.O6 
    ##            1            1            1            1            1            1 
    ##  10052.MG.O7  10052.MG.O8  10052.MG.O9 
    ##            1            1            1

``` r
# Create matched_samples column
metadata$matched_samples <- ifelse(metadata$sample_name %in% names(sample_counts[sample_counts >= 2]), "YES", "NO")
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 27"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 27"

# Step 5: Visualize the Metadata

``` r
knitr::include_graphics("Metadata_Summary.png")
```

![](Metadata_Summary.png)<!-- -->

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

    ##   BioProject       Run  sample_name Age BMI Sex Region Platform        Country
    ## 1  PRJEB7825 ERR687975 10052.MG.O10  40  NA   M     V2 Illumina Venezuela: GAZ
    ## 2  PRJEB7825 ERR687976 10052.MG.O11  11  NA   M     V2 Illumina Venezuela: GAZ
    ## 3  PRJEB7825 ERR687977 10052.MG.O12  45  NA   F     V2 Illumina Venezuela: GAZ
    ## 4  PRJEB7825 ERR687978 10052.MG.O13  48  NA   F     V2 Illumina Venezuela: GAZ
    ## 5  PRJEB7825 ERR687979 10052.MG.O14  47  NA   F     V2 Illumina Venezuela: GAZ
    ## 6  PRJEB7825 ERR687980 10052.MG.O16  45  NA   F     V2 Illumina Venezuela: GAZ
    ##                    Site public_time PCR_primers
    ## 1 human oral metagenome  2015-04-17    27F-338R
    ## 2 human oral metagenome  2015-04-17    27F-338R
    ## 3 human oral metagenome  2015-04-17    27F-338R
    ## 4 human oral metagenome  2015-04-17    27F-338R
    ## 5 human oral metagenome  2015-04-17    27F-338R
    ## 6 human oral metagenome  2015-04-17    27F-338R
    ##                                       title Disease_study Status
    ## 1 The microbiome of uncontacted Amerindians            NO     HC
    ## 2 The microbiome of uncontacted Amerindians            NO     HC
    ## 3 The microbiome of uncontacted Amerindians            NO     HC
    ## 4 The microbiome of uncontacted Amerindians            NO     HC
    ## 5 The microbiome of uncontacted Amerindians            NO     HC
    ## 6 The microbiome of uncontacted Amerindians            NO     HC
    ##   matched_samples
    ## 1              NO
    ## 2              NO
    ## 3              NO
    ## 4              NO
    ## 5              NO
    ## 6              NO

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
