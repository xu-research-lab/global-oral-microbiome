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
setwd("D:/oral_project/PRJEB11419")

# Read the initial metadata files
initial_meta<-read.csv("SraRunTable.csv",header = T)

# Combine the two metadata files
#meta <- rbind(initial_meta, initial_meta1)

#Remove entries non-oral samples
#initial_meta <- subset(initial_meta,isolate %in% c("saliva","tongue coat"))
initial_meta <- initial_meta[initial_meta$Organism != "feces metagenome
", ]
initial_meta <- initial_meta[initial_meta$Organism != "food metagenome
", ]
initial_meta <- initial_meta[initial_meta$Organism != "gut metagenome
", ]
initial_meta <- initial_meta[initial_meta$Organism != "human blood metagenome", ]
initial_meta <- initial_meta[initial_meta$Organism != "human eye metagenome", ]
initial_meta <- initial_meta[initial_meta$Organism != "human gut metagenome", ]
initial_meta <- initial_meta[initial_meta$Organism != "human metagenome", ]
initial_meta <- initial_meta[initial_meta$Organism != "human nasopharyngeal metagenome", ]
initial_meta <- initial_meta[initial_meta$Organism != "human skin metagenome", ]
initial_meta <- initial_meta[initial_meta$Organism != "human vaginal metagenome", ]
initial_meta <- initial_meta[initial_meta$Organism != "metagenome", ]
# Check the number of samples containing "HT" in Sample.Name (e.g., healthy controls)
#nrow(subset(meta, grepl("HT", Sample.Name)))

# Display the first few rows of the metadata
head(initial_meta)
```

    ##           Run Assay.Type BioProject    BioSample center_name..exp.
    ## 2  ERR1077495      OTHER PRJEB11419 SAMEA3612924            UCSDMI
    ## 4  ERR1077497      OTHER PRJEB11419 SAMEA3612926            UCSDMI
    ## 8  ERR1077501      OTHER PRJEB11419 SAMEA3612930            UCSDMI
    ## 9  ERR1077502      OTHER PRJEB11419 SAMEA3612931            UCSDMI
    ## 11 ERR1077504      OTHER PRJEB11419 SAMEA3612933            UCSDMI
    ## 16 ERR1077509      OTHER PRJEB11419 SAMEA3612938            UCSDMI
    ##                                           Center.Name Consent
    ## 2  UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT  public
    ## 4  UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT  public
    ## 8  UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT  public
    ## 9  UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT  public
    ## 11 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT  public
    ## 16 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT  public
    ##    ENA.FIRST.PUBLIC..run. ENA.LAST.UPDATE..run. Experiment
    ## 2              2015-10-15            2018-11-16 ERX1157510
    ## 4              2015-10-15            2018-11-16 ERX1157512
    ## 8              2015-10-15            2018-11-16 ERX1157516
    ## 9              2015-10-15            2018-11-16 ERX1157517
    ## 11             2015-10-15            2018-11-16 ERX1157519
    ## 16             2015-10-15            2018-11-16 ERX1157524
    ##                                           INSDC_center_name
    ## 2  University of California San Diego Microbiome Initiative
    ## 4  University of California San Diego Microbiome Initiative
    ## 8  University of California San Diego Microbiome Initiative
    ## 9  University of California San Diego Microbiome Initiative
    ## 11 University of California San Diego Microbiome Initiative
    ## 16 University of California San Diego Microbiome Initiative
    ##      INSDC_first_public    INSDC_last_update INSDC_status     Instrument
    ## 2  2015-10-15T17:01:25Z 2017-08-10T13:39:52Z       public Illumina MiSeq
    ## 4  2015-10-15T17:01:20Z 2017-08-10T13:39:52Z       public Illumina MiSeq
    ## 8  2015-10-15T17:01:22Z 2017-08-10T13:39:52Z       public Illumina MiSeq
    ## 9  2015-10-15T17:01:24Z 2017-08-10T13:39:52Z       public Illumina MiSeq
    ## 11 2015-10-15T17:01:20Z 2017-08-10T13:39:52Z       public Illumina MiSeq
    ## 16 2015-10-15T17:01:25Z 2017-08-10T13:39:53Z       public Illumina MiSeq
    ##       Library.Name LibraryLayout LibrarySelection LibrarySource
    ## 2  10317.000007751        SINGLE              PCR   METAGENOMIC
    ## 4  10317.000007753        SINGLE              PCR   METAGENOMIC
    ## 8  10317.000007758        SINGLE              PCR   METAGENOMIC
    ## 9  10317.000007759        SINGLE              PCR   METAGENOMIC
    ## 11 10317.000007761        SINGLE              PCR   METAGENOMIC
    ## 16 10317.000007767        SINGLE              PCR   METAGENOMIC
    ##                 Organism Platform          ReleaseDate  Sample.Name
    ## 2  human oral metagenome ILLUMINA 2015-10-28T00:00:00Z SAMEA3612924
    ## 4  human oral metagenome ILLUMINA 2015-10-28T00:00:00Z SAMEA3612926
    ## 8  human oral metagenome ILLUMINA 2015-10-28T00:00:00Z SAMEA3612930
    ## 9  human oral metagenome ILLUMINA 2015-10-28T00:00:00Z SAMEA3612931
    ## 11 human oral metagenome ILLUMINA 2015-10-28T00:00:00Z SAMEA3612933
    ## 16 human oral metagenome ILLUMINA 2015-10-28T00:00:00Z SAMEA3612938
    ##                        Sample_name   sequencing_meth..exp. SRA.Study   Bytes
    ## 2  qiita_sid_10317:10317.000007751 Sequencing by synthesis ERP012803  370128
    ## 4  qiita_sid_10317:10317.000007753 Sequencing by synthesis ERP012803  590452
    ## 8  qiita_sid_10317:10317.000007758 Sequencing by synthesis ERP012803 1002326
    ## 9  qiita_sid_10317:10317.000007759 Sequencing by synthesis ERP012803  598346
    ## 11 qiita_sid_10317:10317.000007761 Sequencing by synthesis ERP012803  563072
    ## 16 qiita_sid_10317:10317.000007767 Sequencing by synthesis ERP012803 1269547
    ##    DATASTORE.filetype DATASTORE.provider
    ## 2    fastq,run.zq,sra     ena,gs,ncbi,s3
    ## 4    fastq,run.zq,sra     ena,gs,ncbi,s3
    ## 8    fastq,run.zq,sra     ena,gs,ncbi,s3
    ## 9    fastq,run.zq,sra     ena,gs,ncbi,s3
    ## 11   fastq,run.zq,sra     ena,gs,ncbi,s3
    ## 16   fastq,run.zq,sra     ena,gs,ncbi,s3
    ##                            DATASTORE.region well_description..exp.
    ## 2  ena,gs.us-east1,ncbi.public,s3.us-east-1         AG_41_7751_C12
    ## 4  ena,gs.us-east1,ncbi.public,s3.us-east-1         AG_41_7753_E12
    ## 8  ena,gs.us-east1,ncbi.public,s3.us-east-1         AG_41_7758_F11
    ## 9  ena,gs.us-east1,ncbi.public,s3.us-east-1         AG_41_7759_B10
    ## 11 ena,gs.us-east1,ncbi.public,s3.us-east-1         AG_41_7761_E11
    ## 16 ena,gs.us-east1,ncbi.public,s3.us-east-1          AG_41_7767_E9
    ##             create_date version center_project_name..exp.     run_prefix..exp.
    ## 2  2015-12-08T15:40:00Z       1                      AG16 Undetermined_S0_L001
    ## 4  2015-12-08T15:40:00Z       1                      AG16 Undetermined_S0_L001
    ## 8  2015-12-08T15:39:00Z       1                      AG16 Undetermined_S0_L001
    ## 9  2015-12-08T15:39:00Z       1                      AG16 Undetermined_S0_L001
    ## 11 2015-12-08T15:39:00Z       1                      AG16 Undetermined_S0_L001
    ## 16 2015-12-08T15:39:00Z       1                      AG16 Undetermined_S0_L001
    ##    AvgSpotLen   Bases run_date..exp. run_center..exp. sample_plate..exp.
    ## 2         151  625442        8/21/14             CCME              AG_41
    ## 4         151 1017589        8/21/14             CCME              AG_41
    ## 8         151 1791162        8/21/14             CCME              AG_41
    ## 9         151 1014418        8/21/14             CCME              AG_41
    ## 11        151  959152        8/21/14             CCME              AG_41
    ## 16        151 2213811        8/21/14             CCME              AG_41
    ##     external_id                    Submitter_Id prepared_meals_frequency
    ## 2  SAMEA3612924 qiita_sid_10317:10317.000007751             Not provided
    ## 4  SAMEA3612926 qiita_sid_10317:10317.000007753             Not provided
    ## 8  SAMEA3612930 qiita_sid_10317:10317.000007758             Not provided
    ## 9  SAMEA3612931 qiita_sid_10317:10317.000007759             Not provided
    ## 11 SAMEA3612933 qiita_sid_10317:10317.000007761             Not provided
    ## 16 SAMEA3612938 qiita_sid_10317:10317.000007767             Not provided
    ##    country_of_birth Birth_Year red_meat_frequency diet_type      Race
    ## 2           Belgium     2006.0       Not provided  Omnivore Caucasian
    ## 4           Belgium     1978.0       Not provided  Omnivore Caucasian
    ## 8           Belgium     2009.0       Not provided  Omnivore Caucasian
    ## 9           Belgium     1978.0       Not provided  Omnivore Caucasian
    ## 11          Belgium     1977.0       Not provided  Omnivore Caucasian
    ## 16          Belgium     1978.0       Not provided  Omnivore Caucasian
    ##                                antibiotic_history            exercise_frequency
    ## 2                                            Year Occasionally (1-2 times/week)
    ## 4  I have not taken antibiotics in the past year.                         Never
    ## 8  I have not taken antibiotics in the past year. Occasionally (1-2 times/week)
    ## 9  I have not taken antibiotics in the past year.                         Never
    ## 11 I have not taken antibiotics in the past year. Occasionally (1-2 times/week)
    ## 16 I have not taken antibiotics in the past year.                         Never
    ##    vitamin_b_supplement_frequency vitamin_d_supplement_frequency
    ## 2                    Not provided                   Not provided
    ## 4                    Not provided                   Not provided
    ## 8                    Not provided                   Not provided
    ## 9                    Not provided                   Not provided
    ## 11                   Not provided                   Not provided
    ## 16                   Not provided                   Not provided
    ##    bowel_movement_frequency             contraceptive whole_grain_frequency
    ## 2              Not provided              Not provided          Not provided
    ## 4              Not provided Yes, I am taking the pill          Not provided
    ## 8              Not provided                        No          Not provided
    ## 9              Not provided Yes, I am taking the pill          Not provided
    ## 11             Not provided              Not provided          Not provided
    ## 16             Not provided Yes, I am taking the pill          Not provided
    ##    age_cat probiotic_frequency sugary_sweets_frequency cosmetics_frequency
    ## 2    child        Not provided            Not provided               Never
    ## 4      30s        Not provided            Not provided               Daily
    ## 8    child        Not provided            Not provided               Never
    ## 9      30s        Not provided            Not provided               Daily
    ## 11     30s        Not provided            Not provided               Daily
    ## 16     30s        Not provided            Not provided               Daily
    ##                   pool_frequency   whole_eggs exercise_location types_of_plants
    ## 2  Occasionally (1-2 times/week) Not provided           Indoors         6 to 10
    ## 4                          Never Not provided None of the above     Less than 5
    ## 8     Rarely (a few times/month) Not provided           Indoors         6 to 10
    ## 9                          Never Not provided None of the above     Less than 5
    ## 11                         Never Not provided           Indoors         6 to 10
    ## 16                         Never Not provided None of the above     Less than 5
    ##    vivid_dreams    olive_oil  last_travel level_of_education
    ## 2  Not provided Not provided Not provided       Not provided
    ## 4  Not provided Not provided        Month       Not provided
    ## 8  Not provided Not provided Not provided       Not provided
    ## 9  Not provided Not provided        Month       Not provided
    ## 11 Not provided Not provided        Month       Not provided
    ## 16 Not provided Not provided        Month       Not provided
    ##    fermented_plant_frequency          alcohol_frequency drinks_per_session
    ## 2               Not provided                      Never       Not provided
    ## 4               Not provided Regularly (3-5 times/week)       Not provided
    ## 8               Not provided                      Never       Not provided
    ## 9               Not provided Regularly (3-5 times/week)       Not provided
    ## 11              Not provided Regularly (3-5 times/week)       Not provided
    ## 16              Not provided Regularly (3-5 times/week)       Not provided
    ##    salted_snacks_frequency sugar_sweetened_drink_frequency
    ## 2             Not provided                    Not provided
    ## 4             Not provided                    Not provided
    ## 8             Not provided                    Not provided
    ## 9             Not provided                    Not provided
    ## 11            Not provided                    Not provided
    ## 16            Not provided                    Not provided
    ##    teethbrushing_frequency one_liter_of_water_a_day_frequency    last_move
    ## 2                    Daily                       Not provided Not provided
    ## 4                    Daily                       Not provided Not provided
    ## 8                    Daily                       Not provided Not provided
    ## 9                    Daily                       Not provided Not provided
    ## 11                   Daily                       Not provided Not provided
    ## 16                   Daily                       Not provided Not provided
    ##    poultry_frequency ready_to_eat_meals_frequency seafood_frequency
    ## 2       Not provided                 Not provided      Not provided
    ## 4       Not provided                 Not provided      Not provided
    ## 8       Not provided                 Not provided      Not provided
    ## 9       Not provided                 Not provided      Not provided
    ## 11      Not provided                 Not provided      Not provided
    ## 16      Not provided                 Not provided      Not provided
    ##       sleep_duration smoking_frequency      thyroid vegetable_frequency
    ## 2    8 or more hours             Never Not provided        Not provided
    ## 4          7-8 hours             Never Not provided        Not provided
    ## 8    8 or more hours             Never Not provided        Not provided
    ## 9          7-8 hours             Never Not provided        Not provided
    ## 11 Less than 5 hours             Never Not provided        Not provided
    ## 16         7-8 hours             Never Not provided        Not provided
    ##                                                deodorant_use
    ## 2                I do not use deodorant or an antiperspirant
    ## 4  Not sure, but I use some form of deodorant/antiperspirant
    ## 8                I do not use deodorant or an antiperspirant
    ## 9  Not sure, but I use some form of deodorant/antiperspirant
    ## 11 Not sure, but I use some form of deodorant/antiperspirant
    ## 16 Not sure, but I use some form of deodorant/antiperspirant
    ##    drinking_water_source fed_as_infant flossing_frequency
    ## 2               Filtered  Not provided              Never
    ## 4               Filtered  Not provided              Never
    ## 8               Filtered  Not provided              Never
    ## 9               Filtered  Not provided              Never
    ## 11              Filtered  Not provided              Never
    ## 16              Filtered  Not provided              Never
    ##                                       flu_vaccine_date fruit_frequency
    ## 2  I have not gotten the flu vaccine in the past year.    Not provided
    ## 4  I have not gotten the flu vaccine in the past year.    Not provided
    ## 8  I have not gotten the flu vaccine in the past year.    Not provided
    ## 9  I have not gotten the flu vaccine in the past year.    Not provided
    ## 11 I have not gotten the flu vaccine in the past year.    Not provided
    ## 16 I have not gotten the flu vaccine in the past year.    Not provided
    ##    fungal_overgrowth homecooked_meals_frequency meat_eggs_frequency
    ## 2       Not provided               Not provided        Not provided
    ## 4       Not provided               Not provided        Not provided
    ## 8       Not provided               Not provided        Not provided
    ## 9       Not provided               Not provided        Not provided
    ## 11      Not provided               Not provided        Not provided
    ## 16      Not provided               Not provided        Not provided
    ##        add_adhd   autoimmune frozen_dessert_frequency
    ## 2  Not provided Not provided             Not provided
    ## 4  Not provided Not provided             Not provided
    ## 8  Not provided Not provided             Not provided
    ## 9  Not provided Not provided             Not provided
    ## 11 Not provided Not provided             Not provided
    ## 16 Not provided Not provided             Not provided
    ##    high_fat_red_meat_frequency milk_cheese_frequency artificial_sweeteners
    ## 2                 Not provided          Not provided          Not provided
    ## 4                 Not provided          Not provided          Not provided
    ## 8                 Not provided          Not provided          Not provided
    ## 9                 Not provided          Not provided          Not provided
    ## 11                Not provided          Not provided          Not provided
    ## 16                Not provided          Not provided          Not provided
    ##    milk_substitute_frequency bowel_movement_quality gluten
    ## 2               Not provided           Not provided     No
    ## 4               Not provided           Not provided     No
    ## 8               Not provided           Not provided     No
    ## 9               Not provided           Not provided     No
    ## 11              Not provided           Not provided     No
    ## 16              Not provided           Not provided     No
    ##                        migraine          ibs roommates_in_study   weight_change
    ## 2  I do not have this condition Not provided       Not provided Remained stable
    ## 4  I do not have this condition Not provided       Not provided Remained stable
    ## 8  I do not have this condition Not provided       Not provided Remained stable
    ## 9  I do not have this condition Not provided       Not provided Remained stable
    ## 11 I do not have this condition Not provided       Not provided Remained stable
    ## 16 I do not have this condition Not provided       Not provided Remained stable
    ##    tonsils_removed               skin_condition appendix_removed
    ## 2            false I do not have this condition            false
    ## 4            false I do not have this condition            false
    ## 8            false I do not have this condition            false
    ## 9            false I do not have this condition            false
    ## 11            true I do not have this condition             true
    ## 16           false I do not have this condition            false
    ##    consume_animal_products_abx     csection     dominant_hand livingwith
    ## 2                 Not provided         true I am right handed        Yes
    ## 4                 Not provided        false I am right handed        Yes
    ## 8                 Not provided         true I am right handed        Yes
    ## 9                 Not provided        false I am right handed        Yes
    ## 11                Not provided Not provided I am right handed        Yes
    ## 16                Not provided        false I am right handed        Yes
    ##             asd chickenpox  acid_reflux breastmilk_formula_ensure
    ## 2  Not provided         No Not provided              Not provided
    ## 4  Not provided        Yes Not provided              Not provided
    ## 8  Not provided         No Not provided              Not provided
    ## 9  Not provided        Yes Not provided              Not provided
    ## 11 Not provided        Yes Not provided              Not provided
    ## 16 Not provided        Yes Not provided              Not provided
    ##    clinical_condition     pregnant         sibo                          ibd
    ## 2        Not provided Not provided Not provided I do not have this condition
    ## 4        Not provided        false Not provided I do not have this condition
    ## 8        Not provided        false Not provided I do not have this condition
    ## 9        Not provided        false Not provided I do not have this condition
    ## 11       Not provided Not provided Not provided I do not have this condition
    ## 16       Not provided        false Not provided I do not have this condition
    ##                        Diabetes cancer_treatment height_cm weight_kg
    ## 2  I do not have this condition     Not provided     120.0      18.0
    ## 4  I do not have this condition     Not provided     169.0      51.0
    ## 8  I do not have this condition     Not provided      96.0      15.0
    ## 9  I do not have this condition     Not provided     169.0      51.0
    ## 11 I do not have this condition     Not provided     177.0      72.0
    ## 16 I do not have this condition     Not provided     169.0      51.0
    ##    INSDC_center_alias          sex country_residence ENA_first_public
    ## 2              UCSDMI Not provided      Not provided       2015-10-15
    ## 4              UCSDMI       female      Not provided       2015-10-15
    ## 8              UCSDMI       female      Not provided       2015-10-15
    ## 9              UCSDMI       female      Not provided       2015-10-15
    ## 11             UCSDMI         male      Not provided       2015-10-15
    ## 16             UCSDMI       female      Not provided       2015-10-15
    ##    ENA_last_update cardiovascular_disease                 lung_disease
    ## 2       2017-08-10           Not provided I do not have this condition
    ## 4       2017-08-10           Not provided I do not have this condition
    ## 8       2017-08-10           Not provided I do not have this condition
    ## 9       2017-08-10           Not provided I do not have this condition
    ## 11      2017-08-10           Not provided I do not have this condition
    ## 16      2017-08-10           Not provided I do not have this condition
    ##    diabetes_type ibd_diagnosis_refined Barcode..exp. linker..exp.
    ## 2   Not provided          Not provided  TACGAGCCCTAA           GT
    ## 4   Not provided          Not provided  ACAATAGACACC           GT
    ## 8   Not provided          Not provided  GATTCCGGCTCA           GT
    ## 9   Not provided          Not provided  CAGCTCATCAGC           GT
    ## 11  Not provided          Not provided  AGTTGAGGCATT           GT
    ## 16  Not provided          Not provided  CATCCCTCTACT           GT
    ##                                    pcr_primers..exp.          primer..exp.
    ## 2  FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT GTGTGCCAGCMGCCGCGGTAA
    ## 4  FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT GTGTGCCAGCMGCCGCGGTAA
    ## 8  FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT GTGTGCCAGCMGCCGCGGTAA
    ## 9  FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT GTGTGCCAGCMGCCGCGGTAA
    ## 11 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT GTGTGCCAGCMGCCGCGGTAA
    ## 16 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT GTGTGCCAGCMGCCGCGGTAA
    ##    target_gene..exp. target_subfragment..exp. liver_disease kidney_disease
    ## 2           16S rRNA                       V4  Not provided   Not provided
    ## 4           16S rRNA                       V4  Not provided   Not provided
    ## 8           16S rRNA                       V4  Not provided   Not provided
    ## 9           16S rRNA                       V4  Not provided   Not provided
    ## 11          16S rRNA                       V4  Not provided   Not provided
    ## 16          16S rRNA                       V4  Not provided   Not provided
    ##           cdiff acne_medication other_supplement_frequency softener
    ## 2  Not provided           false                      false    false
    ## 4  Not provided           false                      false    false
    ## 8  Not provided           false                      false    false
    ## 9  Not provided           false                      false    false
    ## 11 Not provided           false                       true    false
    ## 16 Not provided           false                      false    false
    ##    seasonal_allergies acne_medication_otc nail_biter
    ## 2               false               false      false
    ## 4               false               false      false
    ## 8               false               false      false
    ## 9               false               false      false
    ## 11              false                true      false
    ## 16              false               false      false

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
  Age = NA,
  BMI = "18.5-30",
  Sex = initial_meta$sex,
  Region = initial_meta$target_subfragment..exp.,
  Platform = initial_meta$Instrument,
  Country = initial_meta$country_of_birth,
  Site = initial_meta$Organism,
  public_time = initial_meta$ReleaseDate,
  PCR_primers = "515F-806R",
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
metadata$title <- "American Gut: an Open Platform for Citizen Science Microbiome Research"
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
    ##  10317.000001164  10317.000001179  10317.000001454  10317.000001505 
    ##                1                1                1                1 
    ##  10317.000001697  10317.000001933  10317.000001970  10317.000002037 
    ##                1                1                1                1 
    ##  10317.000002049  10317.000002102  10317.000002178  10317.000002182 
    ##                1                3                1                1 
    ##  10317.000002207  10317.000002221  10317.000002248  10317.000002251 
    ##                1                1                3                3 
    ##  10317.000002259  10317.000002303  10317.000002316  10317.000002320 
    ##                1                1                1                1 
    ##  10317.000002335  10317.000002360  10317.000002381  10317.000002383 
    ##                1                1                1                1 
    ##  10317.000002388  10317.000002390  10317.000002392  10317.000002395 
    ##                1                1                1                1 
    ##  10317.000002413  10317.000002438  10317.000002445  10317.000002476 
    ##                1                1                1                1 
    ##  10317.000002491  10317.000002493  10317.000002504  10317.000002505 
    ##                1                1                1                1 
    ##  10317.000002530  10317.000002532  10317.000002533  10317.000002608 
    ##                1                1                1                1 
    ##  10317.000002620  10317.000002626  10317.000002641  10317.000002642 
    ##                1                1                1                1 
    ##  10317.000002646  10317.000002684  10317.000002688  10317.000002704 
    ##                1                1                2                1 
    ##  10317.000002771  10317.000002794  10317.000002808  10317.000002813 
    ##                1                1                1                1 
    ##  10317.000002820  10317.000002821 10317.000002824A 10317.000002824B 
    ##                1                1                1                1 
    ##  10317.000002841  10317.000002842  10317.000002854  10317.000002864 
    ##                1                1                1                1 
    ##  10317.000002868  10317.000002872  10317.000002881  10317.000002940 
    ##                1                1                1                1 
    ##  10317.000003007  10317.000003027  10317.000003030  10317.000003033 
    ##                1                1                1                1 
    ##  10317.000003080  10317.000003087  10317.000003101  10317.000003104 
    ##                1                1                1                1 
    ##  10317.000003118  10317.000003131  10317.000003166  10317.000003174 
    ##                1                1                1                1 
    ##  10317.000003176  10317.000003205  10317.000003224  10317.000003228 
    ##                1                1                2                1 
    ##  10317.000003236  10317.000003267  10317.000003298  10317.000003321 
    ##                3                1                1                1 
    ##  10317.000003332  10317.000003334  10317.000003338  10317.000003351 
    ##                1                1                1                1 
    ##  10317.000003359  10317.000003367  10317.000003409  10317.000003422 
    ##                1                1                1                1 
    ##  10317.000003437  10317.000003440  10317.000003443  10317.000003511 
    ##                1                1                1                1 
    ##  10317.000003551  10317.000003557  10317.000003558  10317.000003601 
    ##                1                2                3                1 
    ##  10317.000003604  10317.000003639  10317.000003652  10317.000003702 
    ##                1                1                1                1 
    ##  10317.000003720  10317.000003818  10317.000003819  10317.000003949 
    ##                1                1                1                1 
    ##  10317.000003950  10317.000003995  10317.000004010  10317.000004012 
    ##                1                1                1                1 
    ##  10317.000004054  10317.000004060  10317.000004062  10317.000004208 
    ##                1                1                1                1 
    ##  10317.000004239  10317.000004293  10317.000004689  10317.000004808 
    ##                1                1                1                1 
    ##  10317.000004849  10317.000004865  10317.000004873  10317.000004879 
    ##                1                1                1                1 
    ##  10317.000004896  10317.000004903  10317.000004904  10317.000004921 
    ##                1                1                1                1 
    ##  10317.000004925  10317.000004931  10317.000004938  10317.000004953 
    ##                1                1                1                1 
    ##  10317.000004986  10317.000005062  10317.000005064  10317.000005072 
    ##                1                1                1                1 
    ##  10317.000005082  10317.000005087  10317.000005162  10317.000005164 
    ##                1                1                1                1 
    ##  10317.000005220  10317.000005226  10317.000005237  10317.000005259 
    ##                1                1                1                1 
    ##  10317.000005292  10317.000005362  10317.000005422  10317.000005429 
    ##                1                1                1                1 
    ##  10317.000005455  10317.000005552  10317.000005558  10317.000005567 
    ##                1                1                1                1 
    ##  10317.000005572  10317.000005575  10317.000005595  10317.000005635 
    ##                1                1                1                1 
    ##  10317.000005647  10317.000005649  10317.000005658  10317.000005717 
    ##                1                1                1                1 
    ##  10317.000005726  10317.000005789  10317.000005806  10317.000005812 
    ##                1                1                1                1 
    ##  10317.000005813  10317.000005842  10317.000005843  10317.000005891 
    ##                1                1                1                1 
    ##  10317.000005899  10317.000005901  10317.000005911  10317.000005921 
    ##                2                1                1                1 
    ##  10317.000005946  10317.000005972  10317.000006679  10317.000006788 
    ##                1                1                1                1 
    ##  10317.000006825  10317.000006852  10317.000006880  10317.000006914 
    ##                1                1                1                1 
    ##  10317.000006939  10317.000006971  10317.000006974  10317.000006975 
    ##                1                1                1                1 
    ##  10317.000007010  10317.000007021  10317.000007040  10317.000007043 
    ##                1                1                1                1 
    ##  10317.000007052  10317.000007076  10317.000007081  10317.000007082 
    ##                1                1                1                1 
    ##  10317.000007086  10317.000007088  10317.000007089  10317.000007090 
    ##                1                1                1                1 
    ##  10317.000007091  10317.000007092  10317.000007093  10317.000007094 
    ##                1                1                1                1 
    ##  10317.000007095  10317.000007096  10317.000007097  10317.000007098 
    ##                1                1                1                1 
    ##  10317.000007099  10317.000007100  10317.000007101  10317.000007102 
    ##                1                1                1                1 
    ##  10317.000007107  10317.000007109  10317.000007116  10317.000007124 
    ##                1                1                1                1 
    ##  10317.000007125  10317.000007126  10317.000007130  10317.000007132 
    ##                1                1                1                1 
    ##  10317.000007133  10317.000007134  10317.000007135  10317.000007152 
    ##                1                1                1                1 
    ##  10317.000007187  10317.000007196  10317.000007198  10317.000007200 
    ##                1                1                1                1 
    ##  10317.000007216  10317.000007221  10317.000007226  10317.000007237 
    ##                1                1                1                1 
    ##  10317.000007238  10317.000007239  10317.000007264  10317.000007272 
    ##                1                1                1                1 
    ##  10317.000007275  10317.000007280  10317.000007289  10317.000007309 
    ##                1                1                1                1 
    ##  10317.000007312  10317.000007324  10317.000007327  10317.000007332 
    ##                1                1                1                1 
    ##  10317.000007336  10317.000007343  10317.000007348  10317.000007364 
    ##                1                1                1                1 
    ##  10317.000007369  10317.000007379  10317.000007381  10317.000007386 
    ##                1                1                1                1 
    ##  10317.000007391  10317.000007397  10317.000007406  10317.000007409 
    ##                1                1                1                1 
    ##  10317.000007413  10317.000007428  10317.000007433  10317.000007438 
    ##                1                1                1                1 
    ##  10317.000007441  10317.000007447  10317.000007453  10317.000007460 
    ##                1                1                1                1 
    ##  10317.000007465  10317.000007467  10317.000007478  10317.000007479 
    ##                1                1                1                1 
    ##  10317.000007480  10317.000007485  10317.000007490  10317.000007494 
    ##                1                1                1                1 
    ##  10317.000007501  10317.000007508  10317.000007514  10317.000007519 
    ##                1                1                1                1 
    ##  10317.000007523  10317.000007545  10317.000007546  10317.000007550 
    ##                1                1                1                1 
    ##  10317.000007567  10317.000007568  10317.000007571  10317.000007575 
    ##                1                1                1                1 
    ##  10317.000007596  10317.000007600  10317.000007602  10317.000007606 
    ##                1                1                1                1 
    ##  10317.000007612  10317.000007617  10317.000007635  10317.000007672 
    ##                1                1                1                1 
    ##  10317.000007674  10317.000007676  10317.000007684  10317.000007689 
    ##                1                1                1                1 
    ##  10317.000007694  10317.000007727  10317.000007751  10317.000007753 
    ##                1                1                1                1 
    ##  10317.000007758  10317.000007759  10317.000007761  10317.000007767 
    ##                1                1                1                1 
    ##  10317.000007769  10317.000007771  10317.000007772  10317.000007774 
    ##                1                1                1                1 
    ##  10317.000007781  10317.000007784  10317.000007788  10317.000007789 
    ##                1                1                1                1 
    ##  10317.000007791  10317.000007810  10317.000007818  10317.000007822 
    ##                1                1                1                1 
    ##  10317.000007830  10317.000007833  10317.000007834  10317.000007835 
    ##                1                1                1                1 
    ##  10317.000007837  10317.000007843  10317.000007844  10317.000007846 
    ##                1                1                1                1 
    ##  10317.000007850  10317.000007851  10317.000007853  10317.000008195 
    ##                1                1                1                1 
    ##  10317.000008197  10317.000008199  10317.000008203  10317.000009037 
    ##                1                1                1                1 
    ##  10317.000009041  10317.000009051  10317.000009058  10317.000009072 
    ##                1                1                1                1 
    ##  10317.000009076  10317.000009081  10317.000009092  10317.000009106 
    ##                1                1                1                1 
    ##  10317.000009159  10317.000009161  10317.000009181  10317.000009367 
    ##                1                1                1                1 
    ##  10317.000009458  10317.000009471  10317.000009474  10317.000009494 
    ##                1                1                1                1 
    ##  10317.000009539  10317.000009634  10317.000009641  10317.000009649 
    ##                1                1                1                1 
    ##  10317.000009674  10317.000009700  10317.000009714  10317.000009715 
    ##                1                1                1                1 
    ##  10317.000009723  10317.000009756  10317.000010122  10317.000010171 
    ##                1                1                1                1 
    ##  10317.000010253  10317.000010255  10317.000010257  10317.000010313 
    ##                1                3                1                1 
    ##  10317.000010325  10317.000010333  10317.000010335  10317.000010345 
    ##                1                1                1                1 
    ##  10317.000010351  10317.000010357  10317.000010392  10317.000010395 
    ##                1                1                1                1 
    ##  10317.000010401  10317.000010405  10317.000010410  10317.000010415 
    ##                1                1                1                1 
    ##  10317.000010418  10317.000010424  10317.000010433  10317.000010438 
    ##                1                1                1                1 
    ##  10317.000010445  10317.000010464  10317.000010468  10317.000010470 
    ##                1                1                1                1 
    ##  10317.000010739  10317.000010740  10317.000010745  10317.000010810 
    ##                1                1                1                1 
    ##  10317.000010837  10317.000010878  10317.000010915  10317.000010977 
    ##                1                1                1                1 
    ##  10317.000011052  10317.000011233  10317.000011360  10317.000011377 
    ##                1                2                1                1 
    ##  10317.000011938  10317.000011978  10317.000012117  10317.000012312 
    ##                1                1                1                2 
    ##  10317.000012330  10317.000012333  10317.000012340  10317.000012341 
    ##                1                1                1                1 
    ##  10317.000012419  10317.000012435  10317.000012442  10317.000012445 
    ##                1                1                1                1 
    ##  10317.000012450  10317.000012451  10317.000012463  10317.000012473 
    ##                1                1                1                1 
    ##  10317.000012478  10317.000012495  10317.000012505  10317.000012511 
    ##                1                1                1                1 
    ##  10317.000012513  10317.000012522  10317.000012528  10317.000012535 
    ##                1                1                1                1 
    ##  10317.000012542  10317.000012551  10317.000012588  10317.000012604 
    ##                1                1                1                1 
    ##  10317.000012605  10317.000012611  10317.000012619  10317.000012623 
    ##                1                1                1                1 
    ##  10317.000012628  10317.000012662  10317.000012690  10317.000012693 
    ##                1                1                1                1 
    ##  10317.000012698  10317.000012705  10317.000012710  10317.000012742 
    ##                1                1                1                1 
    ##  10317.000012745  10317.000012795  10317.000012806  10317.000012812 
    ##                1                1                1                1 
    ##  10317.000012817  10317.000012825  10317.000012840  10317.000012845 
    ##                1                1                1                1 
    ##  10317.000012850  10317.000012864  10317.000012953  10317.000012997 
    ##                1                1                1                1 
    ##  10317.000013011  10317.000013024  10317.000013050  10317.000013052 
    ##                1                1                1                1 
    ##  10317.000013183  10317.000013280  10317.000013300  10317.000013311 
    ##                3                1                1                1 
    ##  10317.000013341  10317.000013358  10317.000013360  10317.000013376 
    ##                1                1                1                3 
    ##  10317.000013396  10317.000013418  10317.000013512  10317.000013561 
    ##                1                1                1                1 
    ##  10317.000013643  10317.000013687  10317.000013688  10317.000013731 
    ##                1                1                1                1 
    ##  10317.000013833  10317.000013863  10317.000013865  10317.000013866 
    ##                1                1                1                1 
    ##  10317.000013867  10317.000013868  10317.000013869  10317.000013870 
    ##                1                1                1                1 
    ##  10317.000013871  10317.000013872  10317.000013873  10317.000013874 
    ##                1                1                1                1 
    ##  10317.000013875  10317.000013876  10317.000013877  10317.000013878 
    ##                1                1                1                1 
    ##  10317.000013879  10317.000013880  10317.000013881  10317.000013882 
    ##                1                1                1                1 
    ##  10317.000013884  10317.000013885  10317.000013886  10317.000013887 
    ##                1                1                1                1 
    ##  10317.000013888  10317.000013889  10317.000013967  10317.000013998 
    ##                1                1                3                1 
    ##  10317.000014002  10317.000014055  10317.000014147  10317.000014159 
    ##                1                1                2                1 
    ##  10317.000014193  10317.000014205  10317.000014213  10317.000014221 
    ##                1                1                1                1 
    ##  10317.000014225  10317.000014647  10317.000014659  10317.000014683 
    ##                1                1                1                1 
    ##  10317.000014889  10317.000014938  10317.000014971  10317.000014974 
    ##                2                3                1                1 
    ##  10317.000014975  10317.000015011  10317.000015017  10317.000015056 
    ##                1                1                2                2 
    ##  10317.000015063  10317.000015119  10317.000015137  10317.000015153 
    ##                1                1                1                1 
    ##  10317.000015162  10317.000015175  10317.000015177  10317.000015181 
    ##                1                1                1                3 
    ##  10317.000015189  10317.000015206  10317.000015213  10317.000015218 
    ##                2                1                1                1 
    ##  10317.000015225  10317.000015251  10317.000015277  10317.000015281 
    ##                1                1                1                3 
    ##  10317.000015282  10317.000015294  10317.000015304  10317.000015327 
    ##                1                1                3                1 
    ##  10317.000015344  10317.000015348  10317.000015353  10317.000015362 
    ##                2                1                1                3 
    ##  10317.000015384  10317.000015395  10317.000015402  10317.000015446 
    ##                1                1                1                1 
    ##  10317.000015454  10317.000015462  10317.000015493  10317.000015504 
    ##                3                1                3                3 
    ##  10317.000015509  10317.000015521  10317.000015557  10317.000015561 
    ##                3                1                1                3 
    ##  10317.000015648  10317.000015653  10317.000015661  10317.000015679 
    ##                1                1                1                1 
    ##  10317.000015693  10317.000015697  10317.000015698  10317.000015708 
    ##                1                1                1                1 
    ##  10317.000015737  10317.000015747  10317.000015770  10317.000015775 
    ##                1                1                1                1 
    ##  10317.000015929  10317.000015930  10317.000015960  10317.000015975 
    ##                1                1                1                1 
    ##  10317.000016349  10317.000016447  10317.000016459  10317.000016471 
    ##                1                1                3                1 
    ##  10317.000017059  10317.000017141  10317.000017142  10317.000017302 
    ##                1                2                1                1 
    ##  10317.000017303  10317.000017317  10317.000017318  10317.000017689 
    ##                1                3                3                2 
    ##  10317.000017704  10317.000017757  10317.000017767  10317.000017769 
    ##                1                2                3                1 
    ##  10317.000017778  10317.000018053  10317.000018077  10317.000018079 
    ##                1                2                2                3 
    ##  10317.000018313  10317.000020343  10317.000020447  10317.000020485 
    ##                1                3                3                1 
    ##  10317.000020490  10317.000020502  10317.000020604  10317.000020611 
    ##                3                3                3                1 
    ##  10317.000020695  10317.000020696  10317.000020707  10317.000020746 
    ##                3                3                3                1 
    ##  10317.000020748  10317.000020766  10317.000021209  10317.000021246 
    ##                1                3                1                3 
    ##  10317.000021307  10317.000021308  10317.000021595  10317.000021700 
    ##                1                1                1                1 
    ##  10317.000021701  10317.000021800  10317.000022080  10317.000022091 
    ##                1                1                3                1 
    ##  10317.000022095  10317.000022130  10317.000022131  10317.000022140 
    ##                1                3                3                3 
    ##  10317.000022188  10317.000022310  10317.000022331  10317.000022373 
    ##                3                1                1                1 
    ##  10317.000022374  10317.000022398  10317.000022446  10317.000022500 
    ##                1                1                3                1 
    ##  10317.000022535  10317.000022565  10317.000022631  10317.000022687 
    ##                3                1                1                1 
    ##  10317.000022695  10317.000022906  10317.000022911  10317.000022948 
    ##                3                3                3                3 
    ##  10317.000022970  10317.000022992  10317.000023008  10317.000023010 
    ##                3                1                3                3 
    ##  10317.000023018  10317.000023039  10317.000023043  10317.000023045 
    ##                3                3                3                2 
    ##  10317.000023058  10317.000023062  10317.000023075  10317.000023127 
    ##                2                3                3                3 
    ##  10317.000023140  10317.000023174  10317.000023175  10317.000023301 
    ##                1                2                3                3 
    ##  10317.000023302  10317.000023330  10317.000023610  10317.000023645 
    ##                3                3                1                1 
    ##  10317.000023646  10317.000023664  10317.000023665  10317.000023666 
    ##                1                1                1                1 
    ##  10317.000023669  10317.000023729  10317.000023730  10317.000023757 
    ##                1                1                1                1 
    ##  10317.000023816  10317.000023839  10317.000023862  10317.000023875 
    ##                1                1                1                2 
    ##  10317.000023912  10317.000023964  10317.000026043  10317.000026058 
    ##                3                1                3                3 
    ##  10317.000026436  10317.000026509  10317.000026616  10317.000027219 
    ##                1                3                2                1 
    ##  10317.000027270  10317.000027314  10317.000027317  10317.000027335 
    ##                1                1                1                1 
    ##  10317.000027336  10317.000027350  10317.000027357  10317.000027365 
    ##                1                3                3                1 
    ##  10317.000027595  10317.000027660  10317.000027728  10317.000027737 
    ##                1                1                3                1 
    ##  10317.000027738  10317.000027750  10317.000027757  10317.000027764 
    ##                1                1                1                1 
    ##  10317.000027802  10317.000028097  10317.000028148  10317.000028161 
    ##                1                1                1                1 
    ##  10317.000028162  10317.000028240  10317.000028263  10317.000028264 
    ##                1                1                1                1 
    ##  10317.000028288  10317.000028313  10317.000028372  10317.000028373 
    ##                1                1                1                1 
    ##  10317.000028374  10317.000028573  10317.000028665  10317.000028691 
    ##                1                1                2                1 
    ##  10317.000028764  10317.000028774  10317.000028822  10317.000028844 
    ##                3                2                1                2 
    ##  10317.000028853  10317.000028876  10317.000028903  10317.000028917 
    ##                2                3                3                2 
    ##  10317.000028921  10317.000028929  10317.000028934  10317.000028939 
    ##                1                2                1                1 
    ##  10317.000028941  10317.000029244  10317.000029263  10317.000029363 
    ##                2                1                2                1 
    ##  10317.000029372  10317.000029406  10317.000029408  10317.000029438 
    ##                1                1                1                1 
    ##  10317.000029506  10317.000029534  10317.000029555  10317.000029562 
    ##                1                1                1                1 
    ##  10317.000029616  10317.000029617  10317.000029988  10317.000030115 
    ##                2                1                1                1 
    ##  10317.000030121  10317.000030137  10317.000030247  10317.000030308 
    ##                1                1                1                1 
    ##  10317.000030340  10317.000030341  10317.000030347  10317.000030367 
    ##                1                1                1                1 
    ##  10317.000030410  10317.000030488  10317.000030521  10317.000030523 
    ##                1                1                1                1 
    ##  10317.000030545  10317.000030569  10317.000030582  10317.000030584 
    ##                1                1                1                1 
    ##  10317.000030625  10317.000030626  10317.000030968  10317.000030981 
    ##                1                1                1                1 
    ##  10317.000031012  10317.000031032  10317.000031067  10317.000031071 
    ##                1                1                1                1 
    ##  10317.000031279  10317.000031366  10317.000031384  10317.000031442 
    ##                1                1                1                1 
    ##  10317.000031497  10317.000031535  10317.000031816  10317.000031822 
    ##                2                1                1                1 
    ##  10317.000031838  10317.000031867  10317.000031868  10317.000031877 
    ##                1                1                1                2 
    ##  10317.000031880  10317.000031899  10317.000031911  10317.000031923 
    ##                1                1                1                1 
    ##  10317.000032624  10317.000032637  10317.000032820  10317.000032882 
    ##                2                2                1                1 
    ##  10317.000032930  10317.000032934  10317.000032975  10317.000032994 
    ##                1                1                1                2 
    ##  10317.000033023  10317.000033106  10317.000033127  10317.000033217 
    ##                1                2                1                1 
    ##  10317.000033237  10317.000033321  10317.000033392  10317.000033558 
    ##                1                1                1                1 
    ##  10317.000033560  10317.000033567  10317.000033568  10317.000033585 
    ##                1                1                1                1 
    ##  10317.000033592  10317.000033593  10317.000033664  10317.000033674 
    ##                1                1                1                1 
    ##  10317.000033696  10317.000033710  10317.000033736  10317.000033737 
    ##                1                1                1                1 
    ##  10317.000033739  10317.000033793  10317.000033802  10317.000034165 
    ##                4                1                1                1 
    ##  10317.000034599  10317.000034637  10317.000035164  10317.000035167 
    ##                1                1                1                1 
    ##  10317.000035233  10317.000035384  10317.000035444  10317.000035471 
    ##                1                1                1                1 
    ##  10317.000036127  10317.000036171  10317.000036208  10317.000036442 
    ##                1                2                2                1 
    ##  10317.000036449  10317.000036450  10317.000036681  10317.000036684 
    ##                1                1                1                1 
    ##  10317.000036723  10317.000036737  10317.000036883  10317.000037537 
    ##                1                1                1                2 
    ##  10317.000037585  10317.000037587  10317.000037600  10317.000037603 
    ##                1                1                1                1 
    ##  10317.000037945  10317.000038188  10317.000038200  10317.000038219 
    ##                1                2                1                1 
    ##  10317.000038223  10317.000038282  10317.000038363  10317.000038394 
    ##                1                1                1                1 
    ##  10317.000038427  10317.000038436  10317.000038445  10317.000039931 
    ##                1                1                1                2 
    ##  10317.000040134  10317.000040146  10317.000040181  10317.000040183 
    ##                1                1                1                1 
    ##  10317.000040508  10317.000040545  10317.000041055  10317.000041071 
    ##                1                1                1                1 
    ##  10317.000041106  10317.000041306  10317.000041344  10317.000041363 
    ##                1                1                1                1 
    ##  10317.000041365  10317.000041370  10317.000041639  10317.000041668 
    ##                1                1                1                1 
    ##  10317.000041742  10317.000041762  10317.000041810  10317.000041834 
    ##                1                1                1                1 
    ##  10317.000042608  10317.000042609  10317.000042637  10317.000042662 
    ##                1                1                1                1 
    ##  10317.000042693  10317.000042757  10317.000042771  10317.000042792 
    ##                1                1                1                1 
    ##  10317.000042822  10317.000042823  10317.000042898  10317.000042900 
    ##                1                1                1                1 
    ##  10317.000042927  10317.000043103  10317.000043178  10317.000043193 
    ##                1                2                1                1 
    ##  10317.000043206  10317.000044325  10317.000044362  10317.000045002 
    ##                1                1                1                1 
    ##  10317.000045010  10317.000046088  10317.000046126  10317.000046202 
    ##                1                1                1                1 
    ##  10317.000046211  10317.000046221  10317.000046268  10317.000046365 
    ##                1                1                2                1 
    ##  10317.000046381  10317.000047145  10317.000047182  10317.000047213 
    ##                2                1                1                1 
    ##  10317.000047227  10317.000047281  10317.000047319  10317.000047352 
    ##                1                1                1                1 
    ##  10317.000047399  10317.000047487  10317.000047488  10317.000047489 
    ##                1                1                1                1 
    ##  10317.000047559  10317.000047565  10317.000047587  10317.000047590 
    ##                1                1                1                1 
    ##  10317.000047644  10317.000047657  10317.000047757  10317.000047780 
    ##                1                1                1                1 
    ##  10317.000047831  10317.000047836  10317.000047839  10317.000047842 
    ##                2                1                1                2 
    ##  10317.000047854  10317.000047857  10317.000047860  10317.000047863 
    ##                1                1                1                1 
    ##  10317.000047866  10317.000047869  10317.000047879  10317.000047884 
    ##                1                1                1                1 
    ##  10317.000047887  10317.000047899  10317.000047902  10317.000047905 
    ##                1                1                1                1 
    ##  10317.000047932  10317.000047950  10317.000048090  10317.000048249 
    ##                1                1                1                1 
    ##  10317.000048269  10317.000048301  10317.000048361  10317.000048363 
    ##                2                1                2                1 
    ##  10317.000048370  10317.000048381  10317.000048383  10317.000048391 
    ##                1                1                1                1 
    ##  10317.000048403  10317.000048935  10317.000048945  10317.000048955 
    ##                1                1                1                1 
    ##  10317.000049006  10317.000049028  10317.000049120  10317.000049122 
    ##                1                1                1                1 
    ##  10317.000049133  10317.000049220  10317.000049242  10317.000049296 
    ##                2                1                1                1 
    ##  10317.000049731  10317.000049740  10317.000049744  10317.000049749 
    ##                1                1                1                1 
    ##  10317.000049750  10317.000049751  10317.000049752  10317.000049753 
    ##                1                1                1                1 
    ##  10317.000049760  10317.000049765  10317.000049776  10317.000049789 
    ##                1                1                1                1 
    ##  10317.000049805  10317.000049820  10317.000049827  10317.000049830 
    ##                1                1                1                1 
    ##  10317.000049842  10317.000049847  10317.000049878  10317.000049890 
    ##                1                1                1                1 
    ##  10317.000049904  10317.000049910  10317.000049942  10317.000049953 
    ##                1                1                1                1 
    ##  10317.000049969  10317.000050171  10317.000050192  10317.000050216 
    ##                2                1                1                1 
    ##  10317.000050222  10317.000050276  10317.000050332  10317.000050339 
    ##                1                1                1                1 
    ##  10317.000050353  10317.000050360  10317.000050362  10317.000050364 
    ##                1                1                1                1 
    ##  10317.000050368  10317.000050369  10317.000050371  10317.000050377 
    ##                1                1                1                1 
    ##  10317.000050380  10317.000050382  10317.000050388  10317.000050391 
    ##                1                1                1                1 
    ##  10317.000050400  10317.000050405  10317.000050412  10317.000050420 
    ##                1                1                1                1 
    ##  10317.000050438  10317.000050440  10317.000050452  10317.000050459 
    ##                1                1                1                1 
    ##  10317.000050466  10317.000050478  10317.000050495  10317.000050500 
    ##                1                1                1                1 
    ##  10317.000050506  10317.000050520  10317.000050527  10317.000050536 
    ##                1                1                1                1 
    ##  10317.000050538  10317.000050549  10317.000050557  10317.000050565 
    ##                1                1                1                1 
    ##  10317.000050566  10317.000050571  10317.000050577  10317.000050581 
    ##                1                1                1                1 
    ##  10317.000050582  10317.000050592  10317.000050606  10317.000050607 
    ##                1                1                1                1 
    ##  10317.000050614  10317.000050619  10317.000050621  10317.000050627 
    ##                1                1                1                1 
    ##  10317.000050638  10317.000050639  10317.000050661  10317.000050663 
    ##                1                1                1                1 
    ##  10317.000050729  10317.000050748  10317.000050774  10317.000050794 
    ##                1                1                1                1 
    ##  10317.000050804  10317.000050856  10317.000050877  10317.000050890 
    ##                1                1                1                1 
    ##  10317.000050903  10317.000050935  10317.000050975  10317.000050981 
    ##                1                1                1                1 
    ##  10317.000050991  10317.000051010  10317.000051038  10317.000051047 
    ##                1                1                1                1 
    ##  10317.000051058  10317.000051061  10317.000051067  10317.000051142 
    ##                1                1                1                1 
    ##  10317.000051187  10317.000051221  10317.000051226  10317.000051254 
    ##                1                1                1                1 
    ##  10317.000051272  10317.000051533  10317.000051538  10317.000051988 
    ##                1                1                1                1 
    ##  10317.000052035  10317.000052103  10317.000052118  10317.000052179 
    ##                1                1                1                1 
    ##  10317.000052290  10317.000052291  10317.000052294  10317.000052298 
    ##                1                1                1                1 
    ##  10317.000052299  10317.000052300  10317.000052301  10317.000052302 
    ##                1                1                1                1 
    ##  10317.000052304  10317.000052305  10317.000052307  10317.000052308 
    ##                1                1                1                1 
    ##  10317.000052309  10317.000052311  10317.000052314  10317.000052315 
    ##                1                1                1                1 
    ##  10317.000052316  10317.000052317  10317.000052318  10317.000052371 
    ##                1                1                1                2 
    ##  10317.000052379  10317.000052473  10317.000052477  10317.000052515 
    ##                1                1                1                1 
    ##  10317.000052862  10317.000053438  10317.000053445  10317.000053456 
    ##                1                2                1                1 
    ##  10317.000053467  10317.000054020  10317.000054048  10317.000054110 
    ##                1                1                1                1 
    ##  10317.000054111  10317.000054269  10317.000054278  10317.000058339 
    ##                1                1                1                1 
    ##  10317.000058356  10317.000058403  10317.000058482  10317.000058490 
    ##                1                1                1                1 
    ##  10317.000058544  10317.000058546  10317.000058584  10317.000058784 
    ##                1                1                1                1 
    ##  10317.000058867  10317.000058952  10317.000058953  10317.000058954 
    ##                1                1                1                1 
    ##  10317.000059074  10317.000059202  10317.000059222  10317.000060025 
    ##                1                1                1                1 
    ##  10317.000060033  10317.000060263  10317.000060268  10317.000060324 
    ##                1                1                1                1 
    ##  10317.000060341  10317.000062111  10317.000062135  10317.000062167 
    ##                1                1                1                1 
    ##  10317.000062168  10317.000062194  10317.000062916  10317.000062945 
    ##                1                1                1                1 
    ##  10317.000062988  10317.000063075  10317.000063102  10317.000063103 
    ##                1                1                1                1 
    ##  10317.000063104  10317.000063105  10317.000063156  10317.000063168 
    ##                1                1                1                1 
    ##  10317.000063180  10317.000063192  10317.000063216  10317.000063228 
    ##                1                1                1                1 
    ##  10317.000063240  10317.000063252  10317.000063264  10317.000063276 
    ##                1                1                1                1 
    ##  10317.000063288  10317.000063300  10317.000063312  10317.000063324 
    ##                1                1                1                1 
    ##  10317.000063336  10317.000063348  10317.000063384  10317.000063419 
    ##                1                1                1                1 
    ##  10317.000063444  10317.000063468  10317.000063480  10317.000063504 
    ##                1                1                1                1 
    ##  10317.000063572  10317.000065414  10317.000065440  10317.000065444 
    ##                1                1                1                1 
    ##  10317.000065469  10317.000065470  10317.000065481  10317.000065493 
    ##                1                1                1                1 
    ##  10317.000065537  10317.000065553  10317.000065584  10317.000065676 
    ##                1                1                1                1 
    ##  10317.000065693  10317.000065773  10317.000066550  10317.000066566 
    ##                1                1                1                1 
    ##  10317.000066650  10317.000066828  10317.000066840  10317.000066852 
    ##                1                1                1                1 
    ##  10317.000066864  10317.000066909  10317.000066917  10317.000066920 
    ##                1                1                1                1 
    ##  10317.000066993  10317.000067069  10317.000067073  10317.000067226 
    ##                1                1                1                1 
    ##  10317.000067343  10317.000067546  10317.000067685  10317.000068020 
    ##                1                1                1                1 
    ##  10317.000068137  10317.000068140  10317.000068254  10317.000068501 
    ##                1                1                1                1 
    ##  10317.000068562  10317.000068647  10317.000068648  10317.000068649 
    ##                1                1                1                1 
    ##  10317.000068650  10317.000068651  10317.000068652  10317.000068653 
    ##                1                1                1                1 
    ##  10317.000068654  10317.000068659  10317.000068660  10317.000068661 
    ##                1                1                1                1 
    ##  10317.000068662  10317.000068704  10317.000068709  10317.000068767 
    ##                1                1                1                1 
    ##  10317.000068779  10317.000068790  10317.000068803  10317.000068815 
    ##                1                1                1                1 
    ##  10317.000068827  10317.000068839  10317.000068851  10317.000068923 
    ##                1                1                1                1 
    ##  10317.000068947  10317.000068983  10317.000068995  10317.000069007 
    ##                1                1                1                1 
    ##  10317.000069019  10317.000069031  10317.000069043  10317.000069055 
    ##                1                1                1                1 
    ##  10317.000069067  10317.000069079  10317.000069091  10317.000069103 
    ##                1                1                1                1 
    ##  10317.000069115  10317.000069139  10317.000069158  10317.000069161 
    ##                1                1                1                1 
    ##  10317.000069181  10317.000069182  10317.000069185  10317.000069190 
    ##                1                1                1                1 
    ##  10317.000069238  10317.000069239  10317.000069266  10317.000069282 
    ##                1                1                1                1 
    ##  10317.000069321  10317.000069322  10317.000069346  10317.000069367 
    ##                1                1                1                1 
    ##  10317.000069368  10317.000069369  10317.000069405  10317.000069419 
    ##                1                1                1                1 
    ##  10317.000069436  10317.000069785  10317.000069786  10317.000069801 
    ##                1                1                1                1 
    ##  10317.000069802  10317.000069820  10317.000069843  10317.000069844 
    ##                1                1                1                1 
    ##  10317.000069857  10317.000069858  10317.000070819  10317.000071492 
    ##                1                1                1                1 
    ##  10317.000071523  10317.000071538  10317.000071892  10317.000071893 
    ##                1                1                1                1 
    ##  10317.000071895  10317.000072002  10317.000072013  10317.000072014 
    ##                1                1                1                1 
    ##  10317.000072039  10317.000072058  10317.000072127  10317.000072187 
    ##                1                1                1                1 
    ##  10317.000072263  10317.000072271  10317.000072297  10317.000072303 
    ##                1                1                1                1 
    ##  10317.000072314  10317.000072329  10317.000072331  10317.000072338 
    ##                1                1                1                1 
    ##  10317.000072359  10317.000072400  10317.000072570  10317.000072625 
    ##                1                1                1                1 
    ##  10317.000072638  10317.000073257  10317.000073317  10317.000073322 
    ##                1                1                1                1 
    ##  10317.000073327  10317.000073347  10317.000073401  10317.000073591 
    ##                1                1                1                1 
    ##  10317.000073626  10317.000073633  10317.000073760  10317.000073761 
    ##                1                1                1                1 
    ##  10317.000073851  10317.000073852  10317.000073853  10317.000074166 
    ##                1                1                1                1 
    ##  10317.000074185  10317.000074317  10317.000074340  10317.000074413 
    ##                1                1                1                1 
    ##  10317.000074414  10317.000074497  10317.000074574  10317.000074595 
    ##                2                1                1                1 
    ##  10317.000074596  10317.000074719  10317.000074759  10317.000074777 
    ##                1                1                1                1 
    ##  10317.000076216  10317.000076219  10317.000076364  10317.000076411 
    ##                1                1                1                1 
    ##  10317.000076463  10317.000076464  10317.000076477  10317.000076546 
    ##                1                1                1                1 
    ##  10317.000076547  10317.000076548  10317.000076577  10317.000076578 
    ##                1                1                1                1 
    ##  10317.000076580  10317.000076596  10317.000076631  10317.000076653 
    ##                1                1                1                1 
    ##  10317.000076699  10317.000076707  10317.000076723  10317.000076724 
    ##                1                1                1                1 
    ##  10317.000076725  10317.000076726  10317.000079063  10317.000079412 
    ##                1                1                1                1 
    ##  10317.000079413  10317.000079425  10317.000079428  10317.000079429 
    ##                1                1                1                1 
    ##  10317.000079765  10317.000079770  10317.000079771  10317.000079774 
    ##                1                1                1                1 
    ##  10317.000079777  10317.000079783  10317.000079784  10317.000079786 
    ##                1                1                1                1 
    ##  10317.000079789  10317.000079790  10317.000079791  10317.000079827 
    ##                1                1                1                1 
    ##  10317.000082576  10317.000082691  10317.000082702  10317.000082709 
    ##                1                1                1                1 
    ##  10317.000082936  10317.000084431  10317.000087098  10317.000087099 
    ##                2                1                1                1 
    ##  10317.000087111  10317.000087606  10317.000089400  10317.000089580 
    ##                1                1                1                1 
    ##  10317.000089587  10317.000089616  10317.000090045  10317.000090074 
    ##                1                1                1                1 
    ##  10317.000090104  10317.000090122  10317.000090153  10317.000090163 
    ##                1                1                1                1 
    ##  10317.000090532  10317.000090543  10317.000090548  10317.000090555 
    ##                1                1                1                1 
    ##  10317.000091028  10317.000092609  10317.000092743  10317.000092751 
    ##                1                1                1                1 
    ##  10317.000092764  10317.000092810  10317.000092820  10317.000092844 
    ##                1                1                1                1 
    ##  10317.000092845  10317.000092878  10317.000092889  10317.000092928 
    ##                1                1                1                1 
    ##  10317.000092964  10317.000092967  10317.000093002  10317.000093003 
    ##                1                1                1                1 
    ##  10317.000093026  10317.000093027  10317.000093117  10317.000093120 
    ##                1                1                1                1 
    ##  10317.000093591  10317.000093598  10317.000093599  10317.000094710 
    ##                1                1                1                1 
    ##  10317.000094717  10317.000094720  10317.000094724  10317.000094745 
    ##                1                1                1                1 
    ##  10317.000094751  10317.000096818  10317.000096826  10317.000096834 
    ##                1                1                1                1 
    ##  10317.000096842  10317.000096850  10317.000096857  10317.000096866 
    ##                1                1                1                1 
    ##  10317.000096874  10317.000096882  10317.000096890  10317.000096898 
    ##                1                1                1                1 
    ##  10317.000096906  10317.000096914  10317.000096922  10317.000096930 
    ##                1                1                1                1 
    ##  10317.000096938  10317.000096946  10317.000096954  10317.000096962 
    ##                1                1                1                1 
    ##  10317.000096970  10317.000096978  10317.000096986  10317.000096994 
    ##                1                1                1                1 
    ##  10317.000097002  10317.000097010  10317.000097018  10317.000097026 
    ##                1                1                1                1 
    ##  10317.000097034  10317.000097042  10317.000097050  10317.000097058 
    ##                1                1                1                1 
    ##  10317.000097066  10317.000097074  10317.000097142  10317.000097150 
    ##                1                1                1                1 
    ##  10317.000097158  10317.000097166  10317.000097174  10317.000097182 
    ##                1                1                1                1 
    ##  10317.000097203  10317.000097266  10317.000097310  10317.000097318 
    ##                1                1                1                1 
    ##  10317.000097365  10317.000097388  10317.000097397  10317.000097412 
    ##                1                1                1                1 
    ##  10317.000097476  10317.000097478  10317.000097482  10317.000097708 
    ##                1                1                1                1 
    ##  10317.000097789  10317.000097800  10317.000097893  10317.000097899 
    ##                1                1                1                1 
    ##  10317.000097913  10317.000098089  10317.000098093  10317.000098095 
    ##                1                1                1                1 
    ##  10317.000098336  10317.000098344  10317.000098352  10317.000098360 
    ##                1                1                1                1 
    ##  10317.000098368  10317.000098376  10317.000098384  10317.000098392 
    ##                1                1                1                1 
    ##  10317.000098400  10317.000098408  10317.000098416  10317.000098424 
    ##                1                1                1                1 
    ##  10317.000098432  10317.000098440  10317.000098448  10317.000098456 
    ##                1                1                1                1 
    ##  10317.000098464  10317.000098472  10317.000098480  10317.000098488 
    ##                1                1                1                1 
    ##  10317.000098496  10317.000098820  10317.000098888  10317.000098979 
    ##                1                1                1                1 
    ##  10317.000099117  10317.000099121  10317.000099227  10317.000099235 
    ##                1                1                1                1 
    ##  10317.000099238  10317.000099254  10317.000100378  10317.000100481 
    ##                1                1                1                1 
    ##  10317.000100488  10317.000101009  10317.000101338  10317.000102652 
    ##                1                1                1                1 
    ##  10317.000102672  10317.000102685  10317.000102691  10317.000102743 
    ##                1                1                1                1 
    ##  10317.000102764  10317.000102767  10317.000102788  10317.000102833 
    ##                1                1                1                1 
    ##  10317.000102971  10317.000103057  10317.000103094  10317.000103732 
    ##                1                1                1                1 
    ##  10317.000103757  10317.000103812  10317.000103816  10317.000103833 
    ##                1                1                1                1 
    ##  10317.000103836  10317.000103840  10317.000103868  10317.000103908 
    ##                1                1                1                1 
    ##  10317.000105341  10317.000105437  10317.000105438  10317.000105525 
    ##                1                1                1                1 
    ##  10317.000105526  10317.000105574  10317.000105575  10317.000105718 
    ##                1                1                1                1 
    ##  10317.000106415  10317.000106460  10317.000106469  10317.000106846 
    ##                1                1                1                1 
    ##  10317.000106847  10317.000107232  10317.000107438  10317.000107492 
    ##                1                1                1                1 
    ##  10317.000107522  10317.000107523  10317.000107526  10317.000107569 
    ##                1                1                1                1 
    ##  10317.000107658  10317.000107949  10317.000107983  10317.000108027 
    ##                1                1                1                1 
    ##  10317.000108034  10317.000108039  10317.000108094  10317.000108096 
    ##                1                1                1                1 
    ##  10317.000108399  10317.000108478  10317.000108555  10317.000108769 
    ##                1                1                1                1 
    ##  10317.000109032  10317.000109047  10317.000109159  10317.000109245 
    ##                1                1                1                1 
    ##  10317.000109452  10317.000109779  10317.000109793  10317.000109843 
    ##                1                1                1                1 
    ##  10317.000109890  10317.000110364  10317.000112116  10317.000112124 
    ##                1                1                1                1 
    ##  10317.000112191  10317.000112222  10317.000112665  10317.000112957 
    ##                1                1                1                1 
    ##  10317.000113027  10317.000113061  10317.000113091  10317.000113595 
    ##                1                1                1                1 
    ##  10317.000113686  10317.000113687  10317.000113702  10317.000113750 
    ##                1                1                1                1 
    ##  10317.000113783  10317.000113805  10317.000114419  10317.000114431 
    ##                1                1                1                1 
    ##  10317.000114469  10317.000114521  10317.000114829  10317.000114862 
    ##                1                1                1                1 
    ##  10317.000114866  10317.000115021  10317.000115069  10317.000145614 
    ##                1                1                1                2 
    ##  10317.000145664  10317.000145667  10317.000145734  10317.000145769 
    ##                2                2                2                2 
    ##  10317.000145841  10317.000145947  10317.000146018  10317.000146250 
    ##                2                2                2                2 
    ##  10317.000146273  10317.000146277  10317.000146411  10317.000146424 
    ##                1                2                2                2 
    ##  10317.000146438  10317.000146445  10317.000146532  10317.000146733 
    ##                2                2                2                2 
    ##  10317.000146735  10317.000146745  10317.000147489  10317.000147491 
    ##                2                1                2                1 
    ##  10317.000147498  10317.000147555  10317.000147720  10317.000147729 
    ##                2                2                2                1 
    ##  10317.000147744  10317.000147842  10317.000147881  10317.000147891 
    ##                2                2                2                2 
    ##  10317.000147894  10317.000147897  10317.000147906  10317.000147920 
    ##                2                2                2                2 
    ##  10317.000147945  10317.000147956  10317.000147960  10317.000147963 
    ##                2                2                2                2 
    ##  10317.000147992  10317.000148025  10317.000148107  10317.000148136 
    ##                2                2                2                2 
    ##  10317.000148200  10317.000148278  10317.000148280  10317.000148316 
    ##                1                2                2                2 
    ##  10317.000148407  10317.000148416  10317.000148418  10317.000148421 
    ##                2                2                2                1 
    ##  10317.000148425  10317.000148500  10317.000148503  10317.000148725 
    ##                2                1                1                2 
    ##  10317.000148793  10317.000148800  10317.000148805  10317.000148808 
    ##                2                2                1                2 
    ##  10317.000148824  10317.000148938  10317.000148974  10317.000148980 
    ##                2                2                2                2 
    ##  10317.000149018  10317.000149036  10317.000149039  10317.000149216 
    ##                2                2                2                2 
    ##  10317.000149220  10317.000149226  10317.000149315  10317.000149343 
    ##                2                2                2                2 
    ##  10317.000149348  10317.000149417  10317.000149424  10317.000149658 
    ##                2                2                2                2 
    ##  10317.000149660  10317.000149664  10317.000149882  10317.000149883 
    ##                2                2                2                2 
    ##  10317.000149975  10317.000149984  10317.000150020  10317.000150042 
    ##                2                2                2                1 
    ##  10317.000150066  10317.000150090  10317.000150146  10317.000150197 
    ##                2                2                1                1 
    ##  10317.000150201  10317.000150222  10317.000150239  10317.000150245 
    ##                2                2                1                1 
    ##  10317.000150249  10317.000150252  10317.000150263  10317.000150270 
    ##                2                1                1                1 
    ##  10317.000150272  10317.000150300  10317.000150486  10317.000150492 
    ##                1                2                1                2 
    ##  10317.000150501  10317.000150510  10317.000150516  10317.000150518 
    ##                2                2                2                2 
    ##  10317.000150521  10317.000150531  10317.000150626  10317.000150633 
    ##                2                2                2                1 
    ##  10317.000150644  10317.000150704  10317.000150990  10317.000150993 
    ##                2                2                2                1 
    ##  10317.000151004  10317.000151037  10317.000151055  10317.000151118 
    ##                1                2                2                2 
    ##  10317.000151124  10317.000151145  10317.000151218  10317.000151220 
    ##                2                2                2                2 
    ##  10317.000151239  10317.000151251  10317.000151256  10317.000151259 
    ##                2                1                2                2 
    ##  10317.000151275  10317.000151287  10317.000151289  10317.000151295 
    ##                2                1                1                1 
    ##  10317.000151301  10317.000151305  10317.000151311  10317.000151317 
    ##                1                1                1                2 
    ##  10317.000151319  10317.000151326  10317.000151329  10317.000151412 
    ##                2                1                1                2 
    ##  10317.000151445  10317.000151745  10317.000151782  10317.000151805 
    ##                2                1                1                1 
    ##  10317.000151812  10317.000151827  10317.000151851  10317.000151865 
    ##                2                2                3                2 
    ##  10317.000151889  10317.000151895  10317.000151905  10317.000151908 
    ##                1                1                2                1 
    ##  10317.000152007  10317.000152232  10317.000152253  10317.000152306 
    ##                1                1                2                1 
    ##  10317.000152340  10317.000152364  10317.000152373  10317.000152418 
    ##                1                2                2                2 
    ##  10317.000152426  10317.000152430  10317.000152474  10317.000152478 
    ##                2                2                2                2 
    ##  10317.000152765  10317.000152951  10317.000153062  10317.000153093 
    ##                2                2                2                2 
    ##  10317.000153168  10317.000154305  10317.000154308  10317.000154338 
    ##                2                1                2                2 
    ##  10317.000154377  10317.000154383  10317.000154386  10317.000154397 
    ##                1                2                2                2 
    ##  10317.000154415  10317.000154427  10317.000154434  10317.000154443 
    ##                2                2                1                2 
    ##  10317.000154457  10317.000154460  10317.000154482  10317.000154509 
    ##                2                2                2                2 
    ##  10317.000154512  10317.000154523  10317.000154526  10317.000154554 
    ##                2                1                1                2 
    ##  10317.000154559  10317.000154571  10317.000154575  10317.000154577 
    ##                2                2                2                2 
    ##  10317.000154596  10317.000154605  10317.000154614  10317.000154616 
    ##                2                2                2                1 
    ##  10317.000154631  10317.000154656  10317.000154662  10317.000154679 
    ##                2                2                2                2 
    ##  10317.000154704  10317.000154724  10317.000154727  10317.000154742 
    ##                1                2                1                2 
    ##  10317.000154785  10317.000154791  10317.000154893  10317.000154896 
    ##                2                2                2                2 
    ##  10317.000154914  10317.000154923  10317.000154946  10317.000155010 
    ##                1                2                1                1 
    ##  10317.000155016  10317.000155039  10317.000155043  10317.000155045 
    ##                2                2                2                2 
    ##  10317.000155058  10317.000155060  10317.000155064  10317.000155102 
    ##                2                2                3                2 
    ##  10317.000155109  10317.000155121  10317.000155132  10317.000155141 
    ##                2                2                2                2 
    ##  10317.000155150  10317.000155184  10317.000155196  10317.000155205 
    ##                2                2                2                2 
    ##  10317.000155208  10317.000155235  10317.000155238  10317.000155252 
    ##                2                2                2                2 
    ##  10317.000155268  10317.000155295  10317.000155298  10317.000155304 
    ##                2                2                2                1 
    ##  10317.000155318  10317.000155372  10317.000155388  10317.000155391 
    ##                2                2                2                1 
    ##  10317.000155396  10317.000155403  10317.000155408  10317.000155411 
    ##                2                2                2                1 
    ##  10317.000155436  10317.000155439  10317.000155580  10317.000155585 
    ##                2                2                1                1 
    ##  10317.000155595  10317.000155600  10317.000155606  10317.000155633 
    ##                2                2                1                2 
    ##  10317.000155651  10317.000155660  10317.000155718  10317.000155852 
    ##                1                2                1                2 
    ##  10317.000155856  10317.000155859  10317.000155862  10317.000155864 
    ##                1                2                2                2 
    ##  10317.000155871  10317.000155907  10317.000155915  10317.000155922 
    ##                2                2                2                2 
    ##  10317.000155930  10317.000155963  10317.000155969  10317.000155976 
    ##                2                2                1                1 
    ##  10317.000155991  10317.000155993  10317.000156009  10317.000156011 
    ##                2                2                2                2 
    ##  10317.000156014  10317.000156035  10317.000156074  10317.000156080 
    ##                1                2                1                2 
    ##  10317.000156083  10317.000156095  10317.000156107  10317.000156123 
    ##                2                2                2                2 
    ##  10317.000156125  10317.000156704  10317.000157472  10317.000157527 
    ##                2                1                2                2 
    ##  10317.000157551  10317.000157566  10317.000157572  10317.000157593 
    ##                2                2                1                2 
    ##  10317.000157617  10317.000157680  10317.000157778  10317.000157812 
    ##                2                2                2                2 
    ##  10317.000157824  10317.000157835  10317.000157844  10317.000157929 
    ##                2                2                2                2 
    ##  10317.000157935  10317.000157983  10317.000157986  10317.000157989 
    ##                2                1                2                1 
    ##  10317.000157994  10317.000158001  10317.000158018  10317.000158031 
    ##                2                1                2                2 
    ##  10317.000158036  10317.000158042  10317.000158057  10317.000158064 
    ##                2                2                1                2 
    ##  10317.000158073  10317.000158090  10317.X00160958  10317.X00160990 
    ##                2                2                1                2 
    ##  10317.X00161003  10317.X00161094  10317.X00161177  10317.X00161226 
    ##                2                1                2                1 
    ##  10317.X00161234  10317.X00161238  10317.X00161282  10317.X00161442 
    ##                1                1                1                1 
    ##  10317.X00161471  10317.X00161504  10317.X00161540  10317.X00161645 
    ##                2                1                1                2 
    ##  10317.X00161684  10317.X00161691  10317.X00161728  10317.X00161746 
    ##                1                1                2                1 
    ##  10317.X00161755  10317.X00161791  10317.X00161827  10317.X00161849 
    ##                1                1                1                1 
    ##  10317.X00161853  10317.X00161854  10317.X00161858  10317.X00161878 
    ##                1                1                1                1 
    ##  10317.X00161897  10317.X00161956  10317.X00162026  10317.X00162061 
    ##                1                1                1                1 
    ##  10317.X00162078  10317.X00162114  10317.X00162143  10317.X00162145 
    ##                1                1                1                2 
    ##  10317.X00162169  10317.X00162295  10317.X00162298  10317.X00162347 
    ##                1                2                2                1 
    ##  10317.X00162352  10317.X00162408  10317.X00162414  10317.X00162439 
    ##                1                1                1                2 
    ##  10317.X00162551  10317.X00162634  10317.X00162738  10317.X00162784 
    ##                1                1                1                1 
    ##  10317.X00162795  10317.X00162841  10317.X00162899  10317.X00162916 
    ##                2                1                1                2 
    ##  10317.X00162972  10317.X00163001  10317.X00163009  10317.X00163024 
    ##                1                1                1                2 
    ##  10317.X00163060  10317.X00163067  10317.X00163077  10317.X00163080 
    ##                1                1                1                1 
    ##  10317.X00163180  10317.X00163195  10317.X00163294  10317.X00163337 
    ##                1                1                1                1 
    ##  10317.X00163361  10317.X00163410  10317.X00163422  10317.X00163570 
    ##                1                2                1                1 
    ##  10317.X00163588  10317.X00163631  10317.X00163642  10317.X00163693 
    ##                1                1                1                1 
    ##  10317.X00163733  10317.X00163778  10317.X00163859  10317.X00163899 
    ##                1                1                2                1 
    ##  10317.X00163921  10317.X00163966  10317.X00163979  10317.X00163983 
    ##                1                1                1                1 
    ##  10317.X00164025  10317.X00164028  10317.X00164032  10317.X00164100 
    ##                1                1                1                2 
    ##  10317.X00164115  10317.X00164163  10317.X00164172  10317.X00164197 
    ##                1                1                1                1 
    ##  10317.X00164237  10317.X00164248  10317.X00164282  10317.X00164368 
    ##                1                1                1                1 
    ##  10317.X00164406  10317.X00164556  10317.X00164607  10317.X00164610 
    ##                1                1                2                1 
    ##  10317.X00164681  10317.X00164685  10317.X00164695  10317.X00164773 
    ##                1                1                1                1 
    ##  10317.X00164793  10317.X00164805  10317.X00164826  10317.X00164843 
    ##                2                1                1                1 
    ##  10317.X00164869  10317.X00164945  10317.X00164981  10317.X00165054 
    ##                2                2                1                1 
    ##  10317.X00165119  10317.X00165142  10317.X00165204  10317.X00165210 
    ##                1                1                2                1 
    ##  10317.X00165213  10317.X00165297  10317.X00165369  10317.X00165391 
    ##                1                1                1                1 
    ##  10317.X00165402  10317.X00165496  10317.X00165505  10317.X00165612 
    ##                1                1                1                1 
    ##  10317.X00165614  10317.X00165628  10317.X00165672  10317.X00165766 
    ##                1                1                1                1 
    ##  10317.X00165788  10317.X00165803  10317.X00165809  10317.X00165822 
    ##                1                1                1                1 
    ##  10317.X00165892  10317.X00166063  10317.X00166131  10317.X00166144 
    ##                1                1                1                1 
    ##  10317.X00166293  10317.X00166298  10317.X00166306  10317.X00166458 
    ##                3                1                1                1 
    ##  10317.X00166463  10317.X00166484  10317.X00166568  10317.X00166575 
    ##                1                1                2                2 
    ##  10317.X00166596  10317.X00166609  10317.X00166656  10317.X00166706 
    ##                2                1                1                1 
    ##  10317.X00166756  10317.X00166759  10317.X00166778  10317.X00166799 
    ##                1                1                3                1 
    ##  10317.X00166817  10317.X00166878  10317.X00167090  10317.X00167097 
    ##                1                1                1                1 
    ##  10317.X00167136  10317.X00167201  10317.X00167222  10317.X00167233 
    ##                1                1                1                1 
    ##  10317.X00167250  10317.X00167287  10317.X00167292  10317.X00167426 
    ##                1                1                1                2 
    ##  10317.X00167446  10317.X00167509  10317.X00167518  10317.X00167638 
    ##                1                1                1                1 
    ##  10317.X00167647  10317.X00167703  10317.X00167706  10317.X00167725 
    ##                1                2                1                1 
    ##  10317.X00167737  10317.X00167794  10317.X00167867  10317.X00167889 
    ##                1                1                2                1 
    ##  10317.X00167897  10317.X00167935  10317.X00168025  10317.X00168035 
    ##                1                2                1                1 
    ##  10317.X00168049  10317.X00168061  10317.X00168068  10317.X00168081 
    ##                1                1                1                1 
    ##  10317.X00168106  10317.X00168165  10317.X00168322  10317.X00168378 
    ##                1                1                1                2 
    ##  10317.X00168400  10317.X00168605  10317.X00168620  10317.X00168784 
    ##                2                1                1                2 
    ##  10317.X00168816  10317.X00168852  10317.X00168905  10317.X00168933 
    ##                2                1                2                2 
    ##  10317.X00168998  10317.X00169002  10317.X00169028  10317.X00169151 
    ##                1                1                1                1 
    ##  10317.X00169174  10317.X00169195  10317.X00169297  10317.X00169317 
    ##                2                1                1                1 
    ##  10317.X00169322  10317.X00169428  10317.X00169446  10317.X00169521 
    ##                1                1                2                1 
    ##  10317.X00169586  10317.X00169631  10317.X00169639  10317.X00169704 
    ##                2                1                1                1 
    ##  10317.X00169792  10317.X00169877  10317.X00169883  10317.X00170015 
    ##                1                1                1                1 
    ##  10317.X00170032  10317.X00170054  10317.X00170066  10317.X00170092 
    ##                1                1                1                2 
    ##  10317.X00170112  10317.X00170114  10317.X00170180  10317.X00170248 
    ##                1                1                1                1 
    ##  10317.X00170257  10317.X00170267  10317.X00170301  10317.X00170417 
    ##                1                1                1                1 
    ##  10317.X00170432  10317.X00170470  10317.X00170474  10317.X00170485 
    ##                2                1                1                1 
    ##  10317.X00170491  10317.X00170500  10317.X00170506  10317.X00170525 
    ##                2                1                1                1 
    ##  10317.X00170551  10317.X00170561  10317.X00170612  10317.X00170622 
    ##                1                1                1                1 
    ##  10317.X00170627  10317.X00170646  10317.X00170663  10317.X00170802 
    ##                1                2                1                1 
    ##  10317.X00170832  10317.X00170869  10317.X00170918  10317.X00170944 
    ##                2                1                2                1 
    ##  10317.X00170968  10317.X00170988  10317.X00171002  10317.X00171025 
    ##                1                1                1                1 
    ##  10317.X00171038  10317.X00171056  10317.X00171088  10317.X00171116 
    ##                1                2                1                1 
    ##  10317.X00171126  10317.X00171251  10317.X00171289  10317.X00171347 
    ##                1                1                1                1 
    ##  10317.X00171376  10317.X00171452  10317.X00171527  10317.X00171534 
    ##                2                2                1                1 
    ##  10317.X00171538  10317.X00171565  10317.X00171709  10317.X00171777 
    ##                1                1                1                1 
    ##  10317.X00171784  10317.X00171803  10317.X00171879  10317.X00171894 
    ##                1                1                1                1 
    ##  10317.X00172002  10317.X00172102  10317.X00172173  10317.X00172177 
    ##                1                1                1                1 
    ##  10317.X00172216  10317.X00172230  10317.X00172254  10317.X00172260 
    ##                2                1                1                1 
    ##  10317.X00172289  10317.X00172577  10317.X00172599  10317.X00172603 
    ##                1                2                1                1 
    ##  10317.X00172632  10317.X00172980  10317.X00172994  10317.X00173001 
    ##                1                1                1                1 
    ##  10317.X00173033  10317.X00173198  10317.X00173218  10317.X00173232 
    ##                1                1                1                1 
    ##  10317.X00173248  10317.X00173299  10317.X00173304  10317.X00173320 
    ##                1                1                1                1 
    ##  10317.X00173359  10317.X00173377  10317.X00173420  10317.X00173451 
    ##                1                1                1                1 
    ##  10317.X00173492  10317.X00173511  10317.X00173550  10317.X00173562 
    ##                2                2                2                1 
    ##  10317.X00173619  10317.X00173629  10317.X00173719  10317.X00173770 
    ##                1                1                1                1 
    ##  10317.X00173879  10317.X00173920  10317.X00173943  10317.X00173960 
    ##                3                1                1                1 
    ##  10317.X00173967  10317.X00173982  10317.X00174046  10317.X00174054 
    ##                1                1                2                1 
    ##  10317.X00174068  10317.X00174071  10317.X00174078  10317.X00174106 
    ##                1                1                1                1 
    ##  10317.X00174145  10317.X00174149  10317.X00174151  10317.X00174248 
    ##                2                1                1                1 
    ##  10317.X00174268  10317.X00174340  10317.X00174347  10317.X00174360 
    ##                1                1                1                1 
    ##  10317.X00174407  10317.X00174427  10317.X00174510  10317.X00174569 
    ##                1                1                1                1 
    ##  10317.X00174573  10317.X00174582  10317.X00174755  10317.X00174885 
    ##                2                1                1                1 
    ##  10317.X00174956  10317.X00174969  10317.X00174995  10317.X00175040 
    ##                1                1                1                1 
    ##  10317.X00175098  10317.X00175241  10317.X00175260  10317.X00175308 
    ##                1                2                1                1 
    ##  10317.X00175353  10317.X00175435  10317.X00175520  10317.X00175593 
    ##                2                1                1                1 
    ##  10317.X00175679  10317.X00175748  10317.X00175868  10317.X00182296 
    ##                1                1                1                1 
    ##  10317.X00191478  10317.X00191497  10317.X00191499  10317.X00191505 
    ##                1                1                1                1 
    ##  10317.X00191523  10317.X00191532  10317.X00191710  10317.X00191736 
    ##                1                1                1                1 
    ##  10317.X00191743  10317.X00191763  10317.X00191767  10317.X00191783 
    ##                1                1                1                1 
    ##  10317.X00191788  10317.X00191814  10317.X00191854  10317.X00191871 
    ##                1                1                1                1 
    ##  10317.X00191878  10317.X00191890  10317.X00191923  10317.X00191928 
    ##                1                1                1                1 
    ##  10317.X00191944  10317.X00191946  10317.X00191953  10317.X00192000 
    ##                1                1                1                1 
    ##  10317.X00192022  10317.X00192027  10317.X00192058  10317.X00192060 
    ##                1                1                1                1 
    ##  10317.X00192084  10317.X00192090  10317.X00192123  10317.X00192144 
    ##                1                1                1                1 
    ##  10317.X00192169  10317.X00192216  10317.X00192231  10317.X00192234 
    ##                1                1                1                1 
    ##  10317.X00192264  10317.X00192309  10317.X00192357  10317.X00192360 
    ##                1                1                1                1 
    ##  10317.X00192381  10317.X00192410  10317.X00192420  10317.X00192427 
    ##                1                1                1                1 
    ##  10317.X00192465  10317.X00192474  10317.X00192490  10317.X00192502 
    ##                1                1                1                1 
    ##  10317.X00192525  10317.X00192550  10317.X00192552  10317.X00192574 
    ##                1                1                1                1 
    ##  10317.X00192585  10317.X00192586  10317.X00192613  10317.X00192616 
    ##                1                1                1                1 
    ##  10317.X00192631  10317.X00192640  10317.X00192642  10317.X00192664 
    ##                1                1                1                1 
    ##  10317.X00192706  10317.X00192774  10317.X00192837  10317.X00192885 
    ##                1                1                1                1 
    ##  10317.X00192900  10317.X00192933 
    ##                1                1

``` r
table(metadata$Site)
```

    ## 
    ##      feces metagenome       food metagenome        gut metagenome 
    ##                     3                    10                    16 
    ## human oral metagenome 
    ##                  2755

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
    ##  10317.000001164  10317.000001179  10317.000001454  10317.000001505 
    ##                1                1                1                1 
    ##  10317.000001697  10317.000001933  10317.000001970  10317.000002037 
    ##                1                1                1                1 
    ##  10317.000002049  10317.000002102  10317.000002178  10317.000002182 
    ##                1                3                1                1 
    ##  10317.000002207  10317.000002221  10317.000002248  10317.000002251 
    ##                1                1                3                3 
    ##  10317.000002259  10317.000002303  10317.000002316  10317.000002320 
    ##                1                1                1                1 
    ##  10317.000002335  10317.000002360  10317.000002381  10317.000002383 
    ##                1                1                1                1 
    ##  10317.000002388  10317.000002390  10317.000002392  10317.000002395 
    ##                1                1                1                1 
    ##  10317.000002413  10317.000002438  10317.000002445  10317.000002476 
    ##                1                1                1                1 
    ##  10317.000002491  10317.000002493  10317.000002504  10317.000002505 
    ##                1                1                1                1 
    ##  10317.000002530  10317.000002532  10317.000002533  10317.000002608 
    ##                1                1                1                1 
    ##  10317.000002620  10317.000002626  10317.000002641  10317.000002642 
    ##                1                1                1                1 
    ##  10317.000002646  10317.000002684  10317.000002688  10317.000002704 
    ##                1                1                2                1 
    ##  10317.000002771  10317.000002794  10317.000002808  10317.000002813 
    ##                1                1                1                1 
    ##  10317.000002820  10317.000002821 10317.000002824A 10317.000002824B 
    ##                1                1                1                1 
    ##  10317.000002841  10317.000002842  10317.000002854  10317.000002864 
    ##                1                1                1                1 
    ##  10317.000002868  10317.000002872  10317.000002881  10317.000002940 
    ##                1                1                1                1 
    ##  10317.000003007  10317.000003027  10317.000003030  10317.000003033 
    ##                1                1                1                1 
    ##  10317.000003080  10317.000003087  10317.000003101  10317.000003104 
    ##                1                1                1                1 
    ##  10317.000003118  10317.000003131  10317.000003166  10317.000003174 
    ##                1                1                1                1 
    ##  10317.000003176  10317.000003205  10317.000003224  10317.000003228 
    ##                1                1                2                1 
    ##  10317.000003236  10317.000003267  10317.000003298  10317.000003321 
    ##                3                1                1                1 
    ##  10317.000003332  10317.000003334  10317.000003338  10317.000003351 
    ##                1                1                1                1 
    ##  10317.000003359  10317.000003367  10317.000003409  10317.000003422 
    ##                1                1                1                1 
    ##  10317.000003437  10317.000003440  10317.000003443  10317.000003511 
    ##                1                1                1                1 
    ##  10317.000003551  10317.000003557  10317.000003558  10317.000003601 
    ##                1                2                3                1 
    ##  10317.000003604  10317.000003639  10317.000003652  10317.000003702 
    ##                1                1                1                1 
    ##  10317.000003720  10317.000003818  10317.000003819  10317.000003949 
    ##                1                1                1                1 
    ##  10317.000003950  10317.000003995  10317.000004010  10317.000004012 
    ##                1                1                1                1 
    ##  10317.000004054  10317.000004060  10317.000004062  10317.000004208 
    ##                1                1                1                1 
    ##  10317.000004239  10317.000004293  10317.000004689  10317.000004808 
    ##                1                1                1                1 
    ##  10317.000004849  10317.000004865  10317.000004873  10317.000004879 
    ##                1                1                1                1 
    ##  10317.000004896  10317.000004903  10317.000004904  10317.000004921 
    ##                1                1                1                1 
    ##  10317.000004925  10317.000004931  10317.000004938  10317.000004953 
    ##                1                1                1                1 
    ##  10317.000004986  10317.000005062  10317.000005064  10317.000005072 
    ##                1                1                1                1 
    ##  10317.000005082  10317.000005087  10317.000005162  10317.000005164 
    ##                1                1                1                1 
    ##  10317.000005220  10317.000005226  10317.000005237  10317.000005259 
    ##                1                1                1                1 
    ##  10317.000005292  10317.000005362  10317.000005422  10317.000005429 
    ##                1                1                1                1 
    ##  10317.000005455  10317.000005552  10317.000005558  10317.000005567 
    ##                1                1                1                1 
    ##  10317.000005572  10317.000005575  10317.000005595  10317.000005635 
    ##                1                1                1                1 
    ##  10317.000005647  10317.000005649  10317.000005658  10317.000005717 
    ##                1                1                1                1 
    ##  10317.000005726  10317.000005789  10317.000005806  10317.000005812 
    ##                1                1                1                1 
    ##  10317.000005813  10317.000005842  10317.000005843  10317.000005891 
    ##                1                1                1                1 
    ##  10317.000005899  10317.000005901  10317.000005911  10317.000005921 
    ##                2                1                1                1 
    ##  10317.000005946  10317.000005972  10317.000006679  10317.000006788 
    ##                1                1                1                1 
    ##  10317.000006825  10317.000006852  10317.000006880  10317.000006914 
    ##                1                1                1                1 
    ##  10317.000006939  10317.000006971  10317.000006974  10317.000006975 
    ##                1                1                1                1 
    ##  10317.000007010  10317.000007021  10317.000007040  10317.000007043 
    ##                1                1                1                1 
    ##  10317.000007052  10317.000007076  10317.000007081  10317.000007082 
    ##                1                1                1                1 
    ##  10317.000007086  10317.000007088  10317.000007089  10317.000007090 
    ##                1                1                1                1 
    ##  10317.000007091  10317.000007092  10317.000007093  10317.000007094 
    ##                1                1                1                1 
    ##  10317.000007095  10317.000007096  10317.000007097  10317.000007098 
    ##                1                1                1                1 
    ##  10317.000007099  10317.000007100  10317.000007101  10317.000007102 
    ##                1                1                1                1 
    ##  10317.000007107  10317.000007109  10317.000007116  10317.000007124 
    ##                1                1                1                1 
    ##  10317.000007125  10317.000007126  10317.000007130  10317.000007132 
    ##                1                1                1                1 
    ##  10317.000007133  10317.000007134  10317.000007135  10317.000007152 
    ##                1                1                1                1 
    ##  10317.000007187  10317.000007196  10317.000007198  10317.000007200 
    ##                1                1                1                1 
    ##  10317.000007216  10317.000007221  10317.000007226  10317.000007237 
    ##                1                1                1                1 
    ##  10317.000007238  10317.000007239  10317.000007264  10317.000007272 
    ##                1                1                1                1 
    ##  10317.000007275  10317.000007280  10317.000007289  10317.000007309 
    ##                1                1                1                1 
    ##  10317.000007312  10317.000007324  10317.000007327  10317.000007332 
    ##                1                1                1                1 
    ##  10317.000007336  10317.000007343  10317.000007348  10317.000007364 
    ##                1                1                1                1 
    ##  10317.000007369  10317.000007379  10317.000007381  10317.000007386 
    ##                1                1                1                1 
    ##  10317.000007391  10317.000007397  10317.000007406  10317.000007409 
    ##                1                1                1                1 
    ##  10317.000007413  10317.000007428  10317.000007433  10317.000007438 
    ##                1                1                1                1 
    ##  10317.000007441  10317.000007447  10317.000007453  10317.000007460 
    ##                1                1                1                1 
    ##  10317.000007465  10317.000007467  10317.000007478  10317.000007479 
    ##                1                1                1                1 
    ##  10317.000007480  10317.000007485  10317.000007490  10317.000007494 
    ##                1                1                1                1 
    ##  10317.000007501  10317.000007508  10317.000007514  10317.000007519 
    ##                1                1                1                1 
    ##  10317.000007523  10317.000007545  10317.000007546  10317.000007550 
    ##                1                1                1                1 
    ##  10317.000007567  10317.000007568  10317.000007571  10317.000007575 
    ##                1                1                1                1 
    ##  10317.000007596  10317.000007600  10317.000007602  10317.000007606 
    ##                1                1                1                1 
    ##  10317.000007612  10317.000007617  10317.000007635  10317.000007672 
    ##                1                1                1                1 
    ##  10317.000007674  10317.000007676  10317.000007684  10317.000007689 
    ##                1                1                1                1 
    ##  10317.000007694  10317.000007727  10317.000007751  10317.000007753 
    ##                1                1                1                1 
    ##  10317.000007758  10317.000007759  10317.000007761  10317.000007767 
    ##                1                1                1                1 
    ##  10317.000007769  10317.000007771  10317.000007772  10317.000007774 
    ##                1                1                1                1 
    ##  10317.000007781  10317.000007784  10317.000007788  10317.000007789 
    ##                1                1                1                1 
    ##  10317.000007791  10317.000007810  10317.000007818  10317.000007822 
    ##                1                1                1                1 
    ##  10317.000007830  10317.000007833  10317.000007834  10317.000007835 
    ##                1                1                1                1 
    ##  10317.000007837  10317.000007843  10317.000007844  10317.000007846 
    ##                1                1                1                1 
    ##  10317.000007850  10317.000007851  10317.000007853  10317.000008195 
    ##                1                1                1                1 
    ##  10317.000008197  10317.000008199  10317.000008203  10317.000009037 
    ##                1                1                1                1 
    ##  10317.000009041  10317.000009051  10317.000009058  10317.000009072 
    ##                1                1                1                1 
    ##  10317.000009076  10317.000009081  10317.000009092  10317.000009106 
    ##                1                1                1                1 
    ##  10317.000009159  10317.000009161  10317.000009181  10317.000009367 
    ##                1                1                1                1 
    ##  10317.000009458  10317.000009471  10317.000009474  10317.000009494 
    ##                1                1                1                1 
    ##  10317.000009539  10317.000009634  10317.000009641  10317.000009649 
    ##                1                1                1                1 
    ##  10317.000009674  10317.000009700  10317.000009714  10317.000009715 
    ##                1                1                1                1 
    ##  10317.000009723  10317.000009756  10317.000010122  10317.000010171 
    ##                1                1                1                1 
    ##  10317.000010253  10317.000010255  10317.000010257  10317.000010313 
    ##                1                3                1                1 
    ##  10317.000010325  10317.000010333  10317.000010335  10317.000010345 
    ##                1                1                1                1 
    ##  10317.000010351  10317.000010357  10317.000010392  10317.000010395 
    ##                1                1                1                1 
    ##  10317.000010401  10317.000010405  10317.000010410  10317.000010415 
    ##                1                1                1                1 
    ##  10317.000010418  10317.000010424  10317.000010433  10317.000010438 
    ##                1                1                1                1 
    ##  10317.000010445  10317.000010464  10317.000010468  10317.000010470 
    ##                1                1                1                1 
    ##  10317.000010739  10317.000010740  10317.000010745  10317.000010810 
    ##                1                1                1                1 
    ##  10317.000010837  10317.000010878  10317.000010915  10317.000010977 
    ##                1                1                1                1 
    ##  10317.000011052  10317.000011233  10317.000011360  10317.000011377 
    ##                1                2                1                1 
    ##  10317.000011938  10317.000011978  10317.000012117  10317.000012312 
    ##                1                1                1                2 
    ##  10317.000012330  10317.000012333  10317.000012340  10317.000012341 
    ##                1                1                1                1 
    ##  10317.000012419  10317.000012435  10317.000012442  10317.000012445 
    ##                1                1                1                1 
    ##  10317.000012450  10317.000012451  10317.000012463  10317.000012473 
    ##                1                1                1                1 
    ##  10317.000012478  10317.000012495  10317.000012505  10317.000012511 
    ##                1                1                1                1 
    ##  10317.000012513  10317.000012522  10317.000012528  10317.000012535 
    ##                1                1                1                1 
    ##  10317.000012542  10317.000012551  10317.000012588  10317.000012604 
    ##                1                1                1                1 
    ##  10317.000012605  10317.000012611  10317.000012619  10317.000012623 
    ##                1                1                1                1 
    ##  10317.000012628  10317.000012662  10317.000012690  10317.000012693 
    ##                1                1                1                1 
    ##  10317.000012698  10317.000012705  10317.000012710  10317.000012742 
    ##                1                1                1                1 
    ##  10317.000012745  10317.000012795  10317.000012806  10317.000012812 
    ##                1                1                1                1 
    ##  10317.000012817  10317.000012825  10317.000012840  10317.000012845 
    ##                1                1                1                1 
    ##  10317.000012850  10317.000012864  10317.000012953  10317.000012997 
    ##                1                1                1                1 
    ##  10317.000013011  10317.000013024  10317.000013050  10317.000013052 
    ##                1                1                1                1 
    ##  10317.000013183  10317.000013280  10317.000013300  10317.000013311 
    ##                3                1                1                1 
    ##  10317.000013341  10317.000013358  10317.000013360  10317.000013376 
    ##                1                1                1                3 
    ##  10317.000013396  10317.000013418  10317.000013512  10317.000013561 
    ##                1                1                1                1 
    ##  10317.000013643  10317.000013687  10317.000013688  10317.000013731 
    ##                1                1                1                1 
    ##  10317.000013833  10317.000013863  10317.000013865  10317.000013866 
    ##                1                1                1                1 
    ##  10317.000013867  10317.000013868  10317.000013869  10317.000013870 
    ##                1                1                1                1 
    ##  10317.000013871  10317.000013872  10317.000013873  10317.000013874 
    ##                1                1                1                1 
    ##  10317.000013875  10317.000013876  10317.000013877  10317.000013878 
    ##                1                1                1                1 
    ##  10317.000013879  10317.000013880  10317.000013881  10317.000013882 
    ##                1                1                1                1 
    ##  10317.000013884  10317.000013885  10317.000013886  10317.000013887 
    ##                1                1                1                1 
    ##  10317.000013888  10317.000013889  10317.000013967  10317.000013998 
    ##                1                1                3                1 
    ##  10317.000014002  10317.000014055  10317.000014147  10317.000014159 
    ##                1                1                2                1 
    ##  10317.000014193  10317.000014205  10317.000014213  10317.000014221 
    ##                1                1                1                1 
    ##  10317.000014225  10317.000014647  10317.000014659  10317.000014683 
    ##                1                1                1                1 
    ##  10317.000014889  10317.000014938  10317.000014971  10317.000014974 
    ##                2                3                1                1 
    ##  10317.000014975  10317.000015011  10317.000015017  10317.000015056 
    ##                1                1                2                2 
    ##  10317.000015063  10317.000015119  10317.000015137  10317.000015153 
    ##                1                1                1                1 
    ##  10317.000015162  10317.000015175  10317.000015177  10317.000015181 
    ##                1                1                1                3 
    ##  10317.000015189  10317.000015206  10317.000015213  10317.000015218 
    ##                2                1                1                1 
    ##  10317.000015225  10317.000015251  10317.000015277  10317.000015281 
    ##                1                1                1                3 
    ##  10317.000015282  10317.000015294  10317.000015304  10317.000015327 
    ##                1                1                3                1 
    ##  10317.000015344  10317.000015348  10317.000015353  10317.000015362 
    ##                2                1                1                3 
    ##  10317.000015384  10317.000015395  10317.000015402  10317.000015446 
    ##                1                1                1                1 
    ##  10317.000015454  10317.000015462  10317.000015493  10317.000015504 
    ##                3                1                3                3 
    ##  10317.000015509  10317.000015521  10317.000015557  10317.000015561 
    ##                3                1                1                3 
    ##  10317.000015648  10317.000015653  10317.000015661  10317.000015679 
    ##                1                1                1                1 
    ##  10317.000015693  10317.000015697  10317.000015698  10317.000015708 
    ##                1                1                1                1 
    ##  10317.000015737  10317.000015747  10317.000015770  10317.000015775 
    ##                1                1                1                1 
    ##  10317.000015929  10317.000015930  10317.000015960  10317.000015975 
    ##                1                1                1                1 
    ##  10317.000016349  10317.000016447  10317.000016459  10317.000016471 
    ##                1                1                3                1 
    ##  10317.000017059  10317.000017141  10317.000017142  10317.000017302 
    ##                1                2                1                1 
    ##  10317.000017303  10317.000017317  10317.000017318  10317.000017689 
    ##                1                3                3                2 
    ##  10317.000017704  10317.000017757  10317.000017767  10317.000017769 
    ##                1                2                3                1 
    ##  10317.000017778  10317.000018053  10317.000018077  10317.000018079 
    ##                1                2                2                3 
    ##  10317.000018313  10317.000020343  10317.000020447  10317.000020485 
    ##                1                3                3                1 
    ##  10317.000020490  10317.000020502  10317.000020604  10317.000020611 
    ##                3                3                3                1 
    ##  10317.000020695  10317.000020696  10317.000020707  10317.000020746 
    ##                3                3                3                1 
    ##  10317.000020748  10317.000020766  10317.000021209  10317.000021246 
    ##                1                3                1                3 
    ##  10317.000021307  10317.000021308  10317.000021595  10317.000021700 
    ##                1                1                1                1 
    ##  10317.000021701  10317.000021800  10317.000022080  10317.000022091 
    ##                1                1                3                1 
    ##  10317.000022095  10317.000022130  10317.000022131  10317.000022140 
    ##                1                3                3                3 
    ##  10317.000022188  10317.000022310  10317.000022331  10317.000022373 
    ##                3                1                1                1 
    ##  10317.000022374  10317.000022398  10317.000022446  10317.000022500 
    ##                1                1                3                1 
    ##  10317.000022535  10317.000022565  10317.000022631  10317.000022687 
    ##                3                1                1                1 
    ##  10317.000022695  10317.000022906  10317.000022911  10317.000022948 
    ##                3                3                3                3 
    ##  10317.000022970  10317.000022992  10317.000023008  10317.000023010 
    ##                3                1                3                3 
    ##  10317.000023018  10317.000023039  10317.000023043  10317.000023045 
    ##                3                3                3                2 
    ##  10317.000023058  10317.000023062  10317.000023075  10317.000023127 
    ##                2                3                3                3 
    ##  10317.000023140  10317.000023174  10317.000023175  10317.000023301 
    ##                1                2                3                3 
    ##  10317.000023302  10317.000023330  10317.000023610  10317.000023645 
    ##                3                3                1                1 
    ##  10317.000023646  10317.000023664  10317.000023665  10317.000023666 
    ##                1                1                1                1 
    ##  10317.000023669  10317.000023729  10317.000023730  10317.000023757 
    ##                1                1                1                1 
    ##  10317.000023816  10317.000023839  10317.000023862  10317.000023875 
    ##                1                1                1                2 
    ##  10317.000023912  10317.000023964  10317.000026043  10317.000026058 
    ##                3                1                3                3 
    ##  10317.000026436  10317.000026509  10317.000026616  10317.000027219 
    ##                1                3                2                1 
    ##  10317.000027270  10317.000027314  10317.000027317  10317.000027335 
    ##                1                1                1                1 
    ##  10317.000027336  10317.000027350  10317.000027357  10317.000027365 
    ##                1                3                3                1 
    ##  10317.000027595  10317.000027660  10317.000027728  10317.000027737 
    ##                1                1                3                1 
    ##  10317.000027738  10317.000027750  10317.000027757  10317.000027764 
    ##                1                1                1                1 
    ##  10317.000027802  10317.000028097  10317.000028148  10317.000028161 
    ##                1                1                1                1 
    ##  10317.000028162  10317.000028240  10317.000028263  10317.000028264 
    ##                1                1                1                1 
    ##  10317.000028288  10317.000028313  10317.000028372  10317.000028373 
    ##                1                1                1                1 
    ##  10317.000028374  10317.000028573  10317.000028665  10317.000028691 
    ##                1                1                2                1 
    ##  10317.000028764  10317.000028774  10317.000028822  10317.000028844 
    ##                3                2                1                2 
    ##  10317.000028853  10317.000028876  10317.000028903  10317.000028917 
    ##                2                3                3                2 
    ##  10317.000028921  10317.000028929  10317.000028934  10317.000028939 
    ##                1                2                1                1 
    ##  10317.000028941  10317.000029244  10317.000029263  10317.000029363 
    ##                2                1                2                1 
    ##  10317.000029372  10317.000029406  10317.000029408  10317.000029438 
    ##                1                1                1                1 
    ##  10317.000029506  10317.000029534  10317.000029555  10317.000029562 
    ##                1                1                1                1 
    ##  10317.000029616  10317.000029617  10317.000029988  10317.000030115 
    ##                2                1                1                1 
    ##  10317.000030121  10317.000030137  10317.000030247  10317.000030308 
    ##                1                1                1                1 
    ##  10317.000030340  10317.000030341  10317.000030347  10317.000030367 
    ##                1                1                1                1 
    ##  10317.000030410  10317.000030488  10317.000030521  10317.000030523 
    ##                1                1                1                1 
    ##  10317.000030545  10317.000030569  10317.000030582  10317.000030584 
    ##                1                1                1                1 
    ##  10317.000030625  10317.000030626  10317.000030968  10317.000030981 
    ##                1                1                1                1 
    ##  10317.000031012  10317.000031032  10317.000031067  10317.000031071 
    ##                1                1                1                1 
    ##  10317.000031279  10317.000031366  10317.000031384  10317.000031442 
    ##                1                1                1                1 
    ##  10317.000031497  10317.000031535  10317.000031816  10317.000031822 
    ##                2                1                1                1 
    ##  10317.000031838  10317.000031867  10317.000031868  10317.000031877 
    ##                1                1                1                2 
    ##  10317.000031880  10317.000031899  10317.000031911  10317.000031923 
    ##                1                1                1                1 
    ##  10317.000032624  10317.000032637  10317.000032820  10317.000032882 
    ##                2                2                1                1 
    ##  10317.000032930  10317.000032934  10317.000032975  10317.000032994 
    ##                1                1                1                2 
    ##  10317.000033023  10317.000033106  10317.000033127  10317.000033217 
    ##                1                2                1                1 
    ##  10317.000033237  10317.000033321  10317.000033392  10317.000033558 
    ##                1                1                1                1 
    ##  10317.000033560  10317.000033567  10317.000033568  10317.000033585 
    ##                1                1                1                1 
    ##  10317.000033592  10317.000033593  10317.000033664  10317.000033674 
    ##                1                1                1                1 
    ##  10317.000033696  10317.000033710  10317.000033736  10317.000033737 
    ##                1                1                1                1 
    ##  10317.000033739  10317.000033793  10317.000033802  10317.000034165 
    ##                4                1                1                1 
    ##  10317.000034599  10317.000034637  10317.000035164  10317.000035167 
    ##                1                1                1                1 
    ##  10317.000035233  10317.000035384  10317.000035444  10317.000035471 
    ##                1                1                1                1 
    ##  10317.000036127  10317.000036171  10317.000036208  10317.000036442 
    ##                1                2                2                1 
    ##  10317.000036449  10317.000036450  10317.000036681  10317.000036684 
    ##                1                1                1                1 
    ##  10317.000036723  10317.000036737  10317.000036883  10317.000037537 
    ##                1                1                1                2 
    ##  10317.000037585  10317.000037587  10317.000037600  10317.000037603 
    ##                1                1                1                1 
    ##  10317.000037945  10317.000038188  10317.000038200  10317.000038219 
    ##                1                2                1                1 
    ##  10317.000038223  10317.000038282  10317.000038363  10317.000038394 
    ##                1                1                1                1 
    ##  10317.000038427  10317.000038436  10317.000038445  10317.000039931 
    ##                1                1                1                2 
    ##  10317.000040134  10317.000040146  10317.000040181  10317.000040183 
    ##                1                1                1                1 
    ##  10317.000040508  10317.000040545  10317.000041055  10317.000041071 
    ##                1                1                1                1 
    ##  10317.000041106  10317.000041306  10317.000041344  10317.000041363 
    ##                1                1                1                1 
    ##  10317.000041365  10317.000041370  10317.000041639  10317.000041668 
    ##                1                1                1                1 
    ##  10317.000041742  10317.000041762  10317.000041810  10317.000041834 
    ##                1                1                1                1 
    ##  10317.000042608  10317.000042609  10317.000042637  10317.000042662 
    ##                1                1                1                1 
    ##  10317.000042693  10317.000042757  10317.000042771  10317.000042792 
    ##                1                1                1                1 
    ##  10317.000042822  10317.000042823  10317.000042898  10317.000042900 
    ##                1                1                1                1 
    ##  10317.000042927  10317.000043103  10317.000043178  10317.000043193 
    ##                1                2                1                1 
    ##  10317.000043206  10317.000044325  10317.000044362  10317.000045002 
    ##                1                1                1                1 
    ##  10317.000045010  10317.000046088  10317.000046126  10317.000046202 
    ##                1                1                1                1 
    ##  10317.000046211  10317.000046221  10317.000046268  10317.000046365 
    ##                1                1                2                1 
    ##  10317.000046381  10317.000047145  10317.000047182  10317.000047213 
    ##                2                1                1                1 
    ##  10317.000047227  10317.000047281  10317.000047319  10317.000047352 
    ##                1                1                1                1 
    ##  10317.000047399  10317.000047487  10317.000047488  10317.000047489 
    ##                1                1                1                1 
    ##  10317.000047559  10317.000047565  10317.000047587  10317.000047590 
    ##                1                1                1                1 
    ##  10317.000047644  10317.000047657  10317.000047757  10317.000047780 
    ##                1                1                1                1 
    ##  10317.000047831  10317.000047836  10317.000047839  10317.000047842 
    ##                2                1                1                2 
    ##  10317.000047854  10317.000047857  10317.000047860  10317.000047863 
    ##                1                1                1                1 
    ##  10317.000047866  10317.000047869  10317.000047879  10317.000047884 
    ##                1                1                1                1 
    ##  10317.000047887  10317.000047899  10317.000047902  10317.000047905 
    ##                1                1                1                1 
    ##  10317.000047932  10317.000047950  10317.000048090  10317.000048249 
    ##                1                1                1                1 
    ##  10317.000048269  10317.000048301  10317.000048361  10317.000048363 
    ##                2                1                2                1 
    ##  10317.000048370  10317.000048381  10317.000048383  10317.000048391 
    ##                1                1                1                1 
    ##  10317.000048403  10317.000048935  10317.000048945  10317.000048955 
    ##                1                1                1                1 
    ##  10317.000049006  10317.000049028  10317.000049120  10317.000049122 
    ##                1                1                1                1 
    ##  10317.000049133  10317.000049220  10317.000049242  10317.000049296 
    ##                2                1                1                1 
    ##  10317.000049731  10317.000049740  10317.000049744  10317.000049749 
    ##                1                1                1                1 
    ##  10317.000049750  10317.000049751  10317.000049752  10317.000049753 
    ##                1                1                1                1 
    ##  10317.000049760  10317.000049765  10317.000049776  10317.000049789 
    ##                1                1                1                1 
    ##  10317.000049805  10317.000049820  10317.000049827  10317.000049830 
    ##                1                1                1                1 
    ##  10317.000049842  10317.000049847  10317.000049878  10317.000049890 
    ##                1                1                1                1 
    ##  10317.000049904  10317.000049910  10317.000049942  10317.000049953 
    ##                1                1                1                1 
    ##  10317.000049969  10317.000050171  10317.000050192  10317.000050216 
    ##                2                1                1                1 
    ##  10317.000050222  10317.000050276  10317.000050332  10317.000050339 
    ##                1                1                1                1 
    ##  10317.000050353  10317.000050360  10317.000050362  10317.000050364 
    ##                1                1                1                1 
    ##  10317.000050368  10317.000050369  10317.000050371  10317.000050377 
    ##                1                1                1                1 
    ##  10317.000050380  10317.000050382  10317.000050388  10317.000050391 
    ##                1                1                1                1 
    ##  10317.000050400  10317.000050405  10317.000050412  10317.000050420 
    ##                1                1                1                1 
    ##  10317.000050438  10317.000050440  10317.000050452  10317.000050459 
    ##                1                1                1                1 
    ##  10317.000050466  10317.000050478  10317.000050495  10317.000050500 
    ##                1                1                1                1 
    ##  10317.000050506  10317.000050520  10317.000050527  10317.000050536 
    ##                1                1                1                1 
    ##  10317.000050538  10317.000050549  10317.000050557  10317.000050565 
    ##                1                1                1                1 
    ##  10317.000050566  10317.000050571  10317.000050577  10317.000050581 
    ##                1                1                1                1 
    ##  10317.000050582  10317.000050592  10317.000050606  10317.000050607 
    ##                1                1                1                1 
    ##  10317.000050614  10317.000050619  10317.000050621  10317.000050627 
    ##                1                1                1                1 
    ##  10317.000050638  10317.000050639  10317.000050661  10317.000050663 
    ##                1                1                1                1 
    ##  10317.000050729  10317.000050748  10317.000050774  10317.000050794 
    ##                1                1                1                1 
    ##  10317.000050804  10317.000050856  10317.000050877  10317.000050890 
    ##                1                1                1                1 
    ##  10317.000050903  10317.000050935  10317.000050975  10317.000050981 
    ##                1                1                1                1 
    ##  10317.000050991  10317.000051010  10317.000051038  10317.000051047 
    ##                1                1                1                1 
    ##  10317.000051058  10317.000051061  10317.000051067  10317.000051142 
    ##                1                1                1                1 
    ##  10317.000051187  10317.000051221  10317.000051226  10317.000051254 
    ##                1                1                1                1 
    ##  10317.000051272  10317.000051533  10317.000051538  10317.000051988 
    ##                1                1                1                1 
    ##  10317.000052035  10317.000052103  10317.000052118  10317.000052179 
    ##                1                1                1                1 
    ##  10317.000052290  10317.000052291  10317.000052294  10317.000052298 
    ##                1                1                1                1 
    ##  10317.000052299  10317.000052300  10317.000052301  10317.000052302 
    ##                1                1                1                1 
    ##  10317.000052304  10317.000052305  10317.000052307  10317.000052308 
    ##                1                1                1                1 
    ##  10317.000052309  10317.000052311  10317.000052314  10317.000052315 
    ##                1                1                1                1 
    ##  10317.000052316  10317.000052317  10317.000052318  10317.000052371 
    ##                1                1                1                2 
    ##  10317.000052379  10317.000052473  10317.000052477  10317.000052515 
    ##                1                1                1                1 
    ##  10317.000052862  10317.000053438  10317.000053445  10317.000053456 
    ##                1                2                1                1 
    ##  10317.000053467  10317.000054020  10317.000054048  10317.000054110 
    ##                1                1                1                1 
    ##  10317.000054111  10317.000054269  10317.000054278  10317.000058339 
    ##                1                1                1                1 
    ##  10317.000058356  10317.000058403  10317.000058482  10317.000058490 
    ##                1                1                1                1 
    ##  10317.000058544  10317.000058546  10317.000058584  10317.000058784 
    ##                1                1                1                1 
    ##  10317.000058867  10317.000058952  10317.000058953  10317.000058954 
    ##                1                1                1                1 
    ##  10317.000059074  10317.000059202  10317.000059222  10317.000060025 
    ##                1                1                1                1 
    ##  10317.000060033  10317.000060263  10317.000060268  10317.000060324 
    ##                1                1                1                1 
    ##  10317.000060341  10317.000062111  10317.000062135  10317.000062167 
    ##                1                1                1                1 
    ##  10317.000062168  10317.000062194  10317.000062916  10317.000062945 
    ##                1                1                1                1 
    ##  10317.000062988  10317.000063075  10317.000063102  10317.000063103 
    ##                1                1                1                1 
    ##  10317.000063104  10317.000063105  10317.000063156  10317.000063168 
    ##                1                1                1                1 
    ##  10317.000063180  10317.000063192  10317.000063216  10317.000063228 
    ##                1                1                1                1 
    ##  10317.000063240  10317.000063252  10317.000063264  10317.000063276 
    ##                1                1                1                1 
    ##  10317.000063288  10317.000063300  10317.000063312  10317.000063324 
    ##                1                1                1                1 
    ##  10317.000063336  10317.000063348  10317.000063384  10317.000063419 
    ##                1                1                1                1 
    ##  10317.000063444  10317.000063468  10317.000063480  10317.000063504 
    ##                1                1                1                1 
    ##  10317.000063572  10317.000065414  10317.000065440  10317.000065444 
    ##                1                1                1                1 
    ##  10317.000065469  10317.000065470  10317.000065481  10317.000065493 
    ##                1                1                1                1 
    ##  10317.000065537  10317.000065553  10317.000065584  10317.000065676 
    ##                1                1                1                1 
    ##  10317.000065693  10317.000065773  10317.000066550  10317.000066566 
    ##                1                1                1                1 
    ##  10317.000066650  10317.000066828  10317.000066840  10317.000066852 
    ##                1                1                1                1 
    ##  10317.000066864  10317.000066909  10317.000066917  10317.000066920 
    ##                1                1                1                1 
    ##  10317.000066993  10317.000067069  10317.000067073  10317.000067226 
    ##                1                1                1                1 
    ##  10317.000067343  10317.000067546  10317.000067685  10317.000068020 
    ##                1                1                1                1 
    ##  10317.000068137  10317.000068140  10317.000068254  10317.000068501 
    ##                1                1                1                1 
    ##  10317.000068562  10317.000068647  10317.000068648  10317.000068649 
    ##                1                1                1                1 
    ##  10317.000068650  10317.000068651  10317.000068652  10317.000068653 
    ##                1                1                1                1 
    ##  10317.000068654  10317.000068659  10317.000068660  10317.000068661 
    ##                1                1                1                1 
    ##  10317.000068662  10317.000068704  10317.000068709  10317.000068767 
    ##                1                1                1                1 
    ##  10317.000068779  10317.000068790  10317.000068803  10317.000068815 
    ##                1                1                1                1 
    ##  10317.000068827  10317.000068839  10317.000068851  10317.000068923 
    ##                1                1                1                1 
    ##  10317.000068947  10317.000068983  10317.000068995  10317.000069007 
    ##                1                1                1                1 
    ##  10317.000069019  10317.000069031  10317.000069043  10317.000069055 
    ##                1                1                1                1 
    ##  10317.000069067  10317.000069079  10317.000069091  10317.000069103 
    ##                1                1                1                1 
    ##  10317.000069115  10317.000069139  10317.000069158  10317.000069161 
    ##                1                1                1                1 
    ##  10317.000069181  10317.000069182  10317.000069185  10317.000069190 
    ##                1                1                1                1 
    ##  10317.000069238  10317.000069239  10317.000069266  10317.000069282 
    ##                1                1                1                1 
    ##  10317.000069321  10317.000069322  10317.000069346  10317.000069367 
    ##                1                1                1                1 
    ##  10317.000069368  10317.000069369  10317.000069405  10317.000069419 
    ##                1                1                1                1 
    ##  10317.000069436  10317.000069785  10317.000069786  10317.000069801 
    ##                1                1                1                1 
    ##  10317.000069802  10317.000069820  10317.000069843  10317.000069844 
    ##                1                1                1                1 
    ##  10317.000069857  10317.000069858  10317.000070819  10317.000071492 
    ##                1                1                1                1 
    ##  10317.000071523  10317.000071538  10317.000071892  10317.000071893 
    ##                1                1                1                1 
    ##  10317.000071895  10317.000072002  10317.000072013  10317.000072014 
    ##                1                1                1                1 
    ##  10317.000072039  10317.000072058  10317.000072127  10317.000072187 
    ##                1                1                1                1 
    ##  10317.000072263  10317.000072271  10317.000072297  10317.000072303 
    ##                1                1                1                1 
    ##  10317.000072314  10317.000072329  10317.000072331  10317.000072338 
    ##                1                1                1                1 
    ##  10317.000072359  10317.000072400  10317.000072570  10317.000072625 
    ##                1                1                1                1 
    ##  10317.000072638  10317.000073257  10317.000073317  10317.000073322 
    ##                1                1                1                1 
    ##  10317.000073327  10317.000073347  10317.000073401  10317.000073591 
    ##                1                1                1                1 
    ##  10317.000073626  10317.000073633  10317.000073760  10317.000073761 
    ##                1                1                1                1 
    ##  10317.000073851  10317.000073852  10317.000073853  10317.000074166 
    ##                1                1                1                1 
    ##  10317.000074185  10317.000074317  10317.000074340  10317.000074413 
    ##                1                1                1                1 
    ##  10317.000074414  10317.000074497  10317.000074574  10317.000074595 
    ##                2                1                1                1 
    ##  10317.000074596  10317.000074719  10317.000074759  10317.000074777 
    ##                1                1                1                1 
    ##  10317.000076216  10317.000076219  10317.000076364  10317.000076411 
    ##                1                1                1                1 
    ##  10317.000076463  10317.000076464  10317.000076477  10317.000076546 
    ##                1                1                1                1 
    ##  10317.000076547  10317.000076548  10317.000076577  10317.000076578 
    ##                1                1                1                1 
    ##  10317.000076580  10317.000076596  10317.000076631  10317.000076653 
    ##                1                1                1                1 
    ##  10317.000076699  10317.000076707  10317.000076723  10317.000076724 
    ##                1                1                1                1 
    ##  10317.000076725  10317.000076726  10317.000079063  10317.000079412 
    ##                1                1                1                1 
    ##  10317.000079413  10317.000079425  10317.000079428  10317.000079429 
    ##                1                1                1                1 
    ##  10317.000079765  10317.000079770  10317.000079771  10317.000079774 
    ##                1                1                1                1 
    ##  10317.000079777  10317.000079783  10317.000079784  10317.000079786 
    ##                1                1                1                1 
    ##  10317.000079789  10317.000079790  10317.000079791  10317.000079827 
    ##                1                1                1                1 
    ##  10317.000082576  10317.000082691  10317.000082702  10317.000082709 
    ##                1                1                1                1 
    ##  10317.000082936  10317.000084431  10317.000087098  10317.000087099 
    ##                2                1                1                1 
    ##  10317.000087111  10317.000087606  10317.000089400  10317.000089580 
    ##                1                1                1                1 
    ##  10317.000089587  10317.000089616  10317.000090045  10317.000090074 
    ##                1                1                1                1 
    ##  10317.000090104  10317.000090122  10317.000090153  10317.000090163 
    ##                1                1                1                1 
    ##  10317.000090532  10317.000090543  10317.000090548  10317.000090555 
    ##                1                1                1                1 
    ##  10317.000091028  10317.000092609  10317.000092743  10317.000092751 
    ##                1                1                1                1 
    ##  10317.000092764  10317.000092810  10317.000092820  10317.000092844 
    ##                1                1                1                1 
    ##  10317.000092845  10317.000092878  10317.000092889  10317.000092928 
    ##                1                1                1                1 
    ##  10317.000092964  10317.000092967  10317.000093002  10317.000093003 
    ##                1                1                1                1 
    ##  10317.000093026  10317.000093027  10317.000093117  10317.000093120 
    ##                1                1                1                1 
    ##  10317.000093591  10317.000093598  10317.000093599  10317.000094710 
    ##                1                1                1                1 
    ##  10317.000094717  10317.000094720  10317.000094724  10317.000094745 
    ##                1                1                1                1 
    ##  10317.000094751  10317.000096818  10317.000096826  10317.000096834 
    ##                1                1                1                1 
    ##  10317.000096842  10317.000096850  10317.000096857  10317.000096866 
    ##                1                1                1                1 
    ##  10317.000096874  10317.000096882  10317.000096890  10317.000096898 
    ##                1                1                1                1 
    ##  10317.000096906  10317.000096914  10317.000096922  10317.000096930 
    ##                1                1                1                1 
    ##  10317.000096938  10317.000096946  10317.000096954  10317.000096962 
    ##                1                1                1                1 
    ##  10317.000096970  10317.000096978  10317.000096986  10317.000096994 
    ##                1                1                1                1 
    ##  10317.000097002  10317.000097010  10317.000097018  10317.000097026 
    ##                1                1                1                1 
    ##  10317.000097034  10317.000097042  10317.000097050  10317.000097058 
    ##                1                1                1                1 
    ##  10317.000097066  10317.000097074  10317.000097142  10317.000097150 
    ##                1                1                1                1 
    ##  10317.000097158  10317.000097166  10317.000097174  10317.000097182 
    ##                1                1                1                1 
    ##  10317.000097203  10317.000097266  10317.000097310  10317.000097318 
    ##                1                1                1                1 
    ##  10317.000097365  10317.000097388  10317.000097397  10317.000097412 
    ##                1                1                1                1 
    ##  10317.000097476  10317.000097478  10317.000097482  10317.000097708 
    ##                1                1                1                1 
    ##  10317.000097789  10317.000097800  10317.000097893  10317.000097899 
    ##                1                1                1                1 
    ##  10317.000097913  10317.000098089  10317.000098093  10317.000098095 
    ##                1                1                1                1 
    ##  10317.000098336  10317.000098344  10317.000098352  10317.000098360 
    ##                1                1                1                1 
    ##  10317.000098368  10317.000098376  10317.000098384  10317.000098392 
    ##                1                1                1                1 
    ##  10317.000098400  10317.000098408  10317.000098416  10317.000098424 
    ##                1                1                1                1 
    ##  10317.000098432  10317.000098440  10317.000098448  10317.000098456 
    ##                1                1                1                1 
    ##  10317.000098464  10317.000098472  10317.000098480  10317.000098488 
    ##                1                1                1                1 
    ##  10317.000098496  10317.000098820  10317.000098888  10317.000098979 
    ##                1                1                1                1 
    ##  10317.000099117  10317.000099121  10317.000099227  10317.000099235 
    ##                1                1                1                1 
    ##  10317.000099238  10317.000099254  10317.000100378  10317.000100481 
    ##                1                1                1                1 
    ##  10317.000100488  10317.000101009  10317.000101338  10317.000102652 
    ##                1                1                1                1 
    ##  10317.000102672  10317.000102685  10317.000102691  10317.000102743 
    ##                1                1                1                1 
    ##  10317.000102764  10317.000102767  10317.000102788  10317.000102833 
    ##                1                1                1                1 
    ##  10317.000102971  10317.000103057  10317.000103094  10317.000103732 
    ##                1                1                1                1 
    ##  10317.000103757  10317.000103812  10317.000103816  10317.000103833 
    ##                1                1                1                1 
    ##  10317.000103836  10317.000103840  10317.000103868  10317.000103908 
    ##                1                1                1                1 
    ##  10317.000105341  10317.000105437  10317.000105438  10317.000105525 
    ##                1                1                1                1 
    ##  10317.000105526  10317.000105574  10317.000105575  10317.000105718 
    ##                1                1                1                1 
    ##  10317.000106415  10317.000106460  10317.000106469  10317.000106846 
    ##                1                1                1                1 
    ##  10317.000106847  10317.000107232  10317.000107438  10317.000107492 
    ##                1                1                1                1 
    ##  10317.000107522  10317.000107523  10317.000107526  10317.000107569 
    ##                1                1                1                1 
    ##  10317.000107658  10317.000107949  10317.000107983  10317.000108027 
    ##                1                1                1                1 
    ##  10317.000108034  10317.000108039  10317.000108094  10317.000108096 
    ##                1                1                1                1 
    ##  10317.000108399  10317.000108478  10317.000108555  10317.000108769 
    ##                1                1                1                1 
    ##  10317.000109032  10317.000109047  10317.000109159  10317.000109245 
    ##                1                1                1                1 
    ##  10317.000109452  10317.000109779  10317.000109793  10317.000109843 
    ##                1                1                1                1 
    ##  10317.000109890  10317.000110364  10317.000112116  10317.000112124 
    ##                1                1                1                1 
    ##  10317.000112191  10317.000112222  10317.000112665  10317.000112957 
    ##                1                1                1                1 
    ##  10317.000113027  10317.000113061  10317.000113091  10317.000113595 
    ##                1                1                1                1 
    ##  10317.000113686  10317.000113687  10317.000113702  10317.000113750 
    ##                1                1                1                1 
    ##  10317.000113783  10317.000113805  10317.000114419  10317.000114431 
    ##                1                1                1                1 
    ##  10317.000114469  10317.000114521  10317.000114829  10317.000114862 
    ##                1                1                1                1 
    ##  10317.000114866  10317.000115021  10317.000115069  10317.000145614 
    ##                1                1                1                2 
    ##  10317.000145664  10317.000145667  10317.000145734  10317.000145769 
    ##                2                2                2                2 
    ##  10317.000145841  10317.000145947  10317.000146018  10317.000146250 
    ##                2                2                2                2 
    ##  10317.000146273  10317.000146277  10317.000146411  10317.000146424 
    ##                1                2                2                2 
    ##  10317.000146438  10317.000146445  10317.000146532  10317.000146733 
    ##                2                2                2                2 
    ##  10317.000146735  10317.000146745  10317.000147489  10317.000147491 
    ##                2                1                2                1 
    ##  10317.000147498  10317.000147555  10317.000147720  10317.000147729 
    ##                2                2                2                1 
    ##  10317.000147744  10317.000147842  10317.000147881  10317.000147891 
    ##                2                2                2                2 
    ##  10317.000147894  10317.000147897  10317.000147906  10317.000147920 
    ##                2                2                2                2 
    ##  10317.000147945  10317.000147956  10317.000147960  10317.000147963 
    ##                2                2                2                2 
    ##  10317.000147992  10317.000148025  10317.000148107  10317.000148136 
    ##                2                2                2                2 
    ##  10317.000148200  10317.000148278  10317.000148280  10317.000148316 
    ##                1                2                2                2 
    ##  10317.000148407  10317.000148416  10317.000148418  10317.000148421 
    ##                2                2                2                1 
    ##  10317.000148425  10317.000148500  10317.000148503  10317.000148725 
    ##                2                1                1                2 
    ##  10317.000148793  10317.000148800  10317.000148805  10317.000148808 
    ##                2                2                1                2 
    ##  10317.000148824  10317.000148938  10317.000148974  10317.000148980 
    ##                2                2                2                2 
    ##  10317.000149018  10317.000149036  10317.000149039  10317.000149216 
    ##                2                2                2                2 
    ##  10317.000149220  10317.000149226  10317.000149315  10317.000149343 
    ##                2                2                2                2 
    ##  10317.000149348  10317.000149417  10317.000149424  10317.000149658 
    ##                2                2                2                2 
    ##  10317.000149660  10317.000149664  10317.000149882  10317.000149883 
    ##                2                2                2                2 
    ##  10317.000149975  10317.000149984  10317.000150020  10317.000150042 
    ##                2                2                2                1 
    ##  10317.000150066  10317.000150090  10317.000150146  10317.000150197 
    ##                2                2                1                1 
    ##  10317.000150201  10317.000150222  10317.000150239  10317.000150245 
    ##                2                2                1                1 
    ##  10317.000150249  10317.000150252  10317.000150263  10317.000150270 
    ##                2                1                1                1 
    ##  10317.000150272  10317.000150300  10317.000150486  10317.000150492 
    ##                1                2                1                2 
    ##  10317.000150501  10317.000150510  10317.000150516  10317.000150518 
    ##                2                2                2                2 
    ##  10317.000150521  10317.000150531  10317.000150626  10317.000150633 
    ##                2                2                2                1 
    ##  10317.000150644  10317.000150704  10317.000150990  10317.000150993 
    ##                2                2                2                1 
    ##  10317.000151004  10317.000151037  10317.000151055  10317.000151118 
    ##                1                2                2                2 
    ##  10317.000151124  10317.000151145  10317.000151218  10317.000151220 
    ##                2                2                2                2 
    ##  10317.000151239  10317.000151251  10317.000151256  10317.000151259 
    ##                2                1                2                2 
    ##  10317.000151275  10317.000151287  10317.000151289  10317.000151295 
    ##                2                1                1                1 
    ##  10317.000151301  10317.000151305  10317.000151311  10317.000151317 
    ##                1                1                1                2 
    ##  10317.000151319  10317.000151326  10317.000151329  10317.000151412 
    ##                2                1                1                2 
    ##  10317.000151445  10317.000151745  10317.000151782  10317.000151805 
    ##                2                1                1                1 
    ##  10317.000151812  10317.000151827  10317.000151851  10317.000151865 
    ##                2                2                3                2 
    ##  10317.000151889  10317.000151895  10317.000151905  10317.000151908 
    ##                1                1                2                1 
    ##  10317.000152007  10317.000152232  10317.000152253  10317.000152306 
    ##                1                1                2                1 
    ##  10317.000152340  10317.000152364  10317.000152373  10317.000152418 
    ##                1                2                2                2 
    ##  10317.000152426  10317.000152430  10317.000152474  10317.000152478 
    ##                2                2                2                2 
    ##  10317.000152765  10317.000152951  10317.000153062  10317.000153093 
    ##                2                2                2                2 
    ##  10317.000153168  10317.000154305  10317.000154308  10317.000154338 
    ##                2                1                2                2 
    ##  10317.000154377  10317.000154383  10317.000154386  10317.000154397 
    ##                1                2                2                2 
    ##  10317.000154415  10317.000154427  10317.000154434  10317.000154443 
    ##                2                2                1                2 
    ##  10317.000154457  10317.000154460  10317.000154482  10317.000154509 
    ##                2                2                2                2 
    ##  10317.000154512  10317.000154523  10317.000154526  10317.000154554 
    ##                2                1                1                2 
    ##  10317.000154559  10317.000154571  10317.000154575  10317.000154577 
    ##                2                2                2                2 
    ##  10317.000154596  10317.000154605  10317.000154614  10317.000154616 
    ##                2                2                2                1 
    ##  10317.000154631  10317.000154656  10317.000154662  10317.000154679 
    ##                2                2                2                2 
    ##  10317.000154704  10317.000154724  10317.000154727  10317.000154742 
    ##                1                2                1                2 
    ##  10317.000154785  10317.000154791  10317.000154893  10317.000154896 
    ##                2                2                2                2 
    ##  10317.000154914  10317.000154923  10317.000154946  10317.000155010 
    ##                1                2                1                1 
    ##  10317.000155016  10317.000155039  10317.000155043  10317.000155045 
    ##                2                2                2                2 
    ##  10317.000155058  10317.000155060  10317.000155064  10317.000155102 
    ##                2                2                3                2 
    ##  10317.000155109  10317.000155121  10317.000155132  10317.000155141 
    ##                2                2                2                2 
    ##  10317.000155150  10317.000155184  10317.000155196  10317.000155205 
    ##                2                2                2                2 
    ##  10317.000155208  10317.000155235  10317.000155238  10317.000155252 
    ##                2                2                2                2 
    ##  10317.000155268  10317.000155295  10317.000155298  10317.000155304 
    ##                2                2                2                1 
    ##  10317.000155318  10317.000155372  10317.000155388  10317.000155391 
    ##                2                2                2                1 
    ##  10317.000155396  10317.000155403  10317.000155408  10317.000155411 
    ##                2                2                2                1 
    ##  10317.000155436  10317.000155439  10317.000155580  10317.000155585 
    ##                2                2                1                1 
    ##  10317.000155595  10317.000155600  10317.000155606  10317.000155633 
    ##                2                2                1                2 
    ##  10317.000155651  10317.000155660  10317.000155718  10317.000155852 
    ##                1                2                1                2 
    ##  10317.000155856  10317.000155859  10317.000155862  10317.000155864 
    ##                1                2                2                2 
    ##  10317.000155871  10317.000155907  10317.000155915  10317.000155922 
    ##                2                2                2                2 
    ##  10317.000155930  10317.000155963  10317.000155969  10317.000155976 
    ##                2                2                1                1 
    ##  10317.000155991  10317.000155993  10317.000156009  10317.000156011 
    ##                2                2                2                2 
    ##  10317.000156014  10317.000156035  10317.000156074  10317.000156080 
    ##                1                2                1                2 
    ##  10317.000156083  10317.000156095  10317.000156107  10317.000156123 
    ##                2                2                2                2 
    ##  10317.000156125  10317.000156704  10317.000157472  10317.000157527 
    ##                2                1                2                2 
    ##  10317.000157551  10317.000157566  10317.000157572  10317.000157593 
    ##                2                2                1                2 
    ##  10317.000157617  10317.000157680  10317.000157778  10317.000157812 
    ##                2                2                2                2 
    ##  10317.000157824  10317.000157835  10317.000157844  10317.000157929 
    ##                2                2                2                2 
    ##  10317.000157935  10317.000157983  10317.000157986  10317.000157989 
    ##                2                1                2                1 
    ##  10317.000157994  10317.000158001  10317.000158018  10317.000158031 
    ##                2                1                2                2 
    ##  10317.000158036  10317.000158042  10317.000158057  10317.000158064 
    ##                2                2                1                2 
    ##  10317.000158073  10317.000158090  10317.X00160958  10317.X00160990 
    ##                2                2                1                2 
    ##  10317.X00161003  10317.X00161094  10317.X00161177  10317.X00161226 
    ##                2                1                2                1 
    ##  10317.X00161234  10317.X00161238  10317.X00161282  10317.X00161442 
    ##                1                1                1                1 
    ##  10317.X00161471  10317.X00161504  10317.X00161540  10317.X00161645 
    ##                2                1                1                2 
    ##  10317.X00161684  10317.X00161691  10317.X00161728  10317.X00161746 
    ##                1                1                2                1 
    ##  10317.X00161755  10317.X00161791  10317.X00161827  10317.X00161849 
    ##                1                1                1                1 
    ##  10317.X00161853  10317.X00161854  10317.X00161858  10317.X00161878 
    ##                1                1                1                1 
    ##  10317.X00161897  10317.X00161956  10317.X00162026  10317.X00162061 
    ##                1                1                1                1 
    ##  10317.X00162078  10317.X00162114  10317.X00162143  10317.X00162145 
    ##                1                1                1                2 
    ##  10317.X00162169  10317.X00162295  10317.X00162298  10317.X00162347 
    ##                1                2                2                1 
    ##  10317.X00162352  10317.X00162408  10317.X00162414  10317.X00162439 
    ##                1                1                1                2 
    ##  10317.X00162551  10317.X00162634  10317.X00162738  10317.X00162784 
    ##                1                1                1                1 
    ##  10317.X00162795  10317.X00162841  10317.X00162899  10317.X00162916 
    ##                2                1                1                2 
    ##  10317.X00162972  10317.X00163001  10317.X00163009  10317.X00163024 
    ##                1                1                1                2 
    ##  10317.X00163060  10317.X00163067  10317.X00163077  10317.X00163080 
    ##                1                1                1                1 
    ##  10317.X00163180  10317.X00163195  10317.X00163294  10317.X00163337 
    ##                1                1                1                1 
    ##  10317.X00163361  10317.X00163410  10317.X00163422  10317.X00163570 
    ##                1                2                1                1 
    ##  10317.X00163588  10317.X00163631  10317.X00163642  10317.X00163693 
    ##                1                1                1                1 
    ##  10317.X00163733  10317.X00163778  10317.X00163859  10317.X00163899 
    ##                1                1                2                1 
    ##  10317.X00163921  10317.X00163966  10317.X00163979  10317.X00163983 
    ##                1                1                1                1 
    ##  10317.X00164025  10317.X00164028  10317.X00164032  10317.X00164100 
    ##                1                1                1                2 
    ##  10317.X00164115  10317.X00164163  10317.X00164172  10317.X00164197 
    ##                1                1                1                1 
    ##  10317.X00164237  10317.X00164248  10317.X00164282  10317.X00164368 
    ##                1                1                1                1 
    ##  10317.X00164406  10317.X00164556  10317.X00164607  10317.X00164610 
    ##                1                1                2                1 
    ##  10317.X00164681  10317.X00164685  10317.X00164695  10317.X00164773 
    ##                1                1                1                1 
    ##  10317.X00164793  10317.X00164805  10317.X00164826  10317.X00164843 
    ##                2                1                1                1 
    ##  10317.X00164869  10317.X00164945  10317.X00164981  10317.X00165054 
    ##                2                2                1                1 
    ##  10317.X00165119  10317.X00165142  10317.X00165204  10317.X00165210 
    ##                1                1                2                1 
    ##  10317.X00165213  10317.X00165297  10317.X00165369  10317.X00165391 
    ##                1                1                1                1 
    ##  10317.X00165402  10317.X00165496  10317.X00165505  10317.X00165612 
    ##                1                1                1                1 
    ##  10317.X00165614  10317.X00165628  10317.X00165672  10317.X00165766 
    ##                1                1                1                1 
    ##  10317.X00165788  10317.X00165803  10317.X00165809  10317.X00165822 
    ##                1                1                1                1 
    ##  10317.X00165892  10317.X00166063  10317.X00166131  10317.X00166144 
    ##                1                1                1                1 
    ##  10317.X00166293  10317.X00166298  10317.X00166306  10317.X00166458 
    ##                3                1                1                1 
    ##  10317.X00166463  10317.X00166484  10317.X00166568  10317.X00166575 
    ##                1                1                2                2 
    ##  10317.X00166596  10317.X00166609  10317.X00166656  10317.X00166706 
    ##                2                1                1                1 
    ##  10317.X00166756  10317.X00166759  10317.X00166778  10317.X00166799 
    ##                1                1                3                1 
    ##  10317.X00166817  10317.X00166878  10317.X00167090  10317.X00167097 
    ##                1                1                1                1 
    ##  10317.X00167136  10317.X00167201  10317.X00167222  10317.X00167233 
    ##                1                1                1                1 
    ##  10317.X00167250  10317.X00167287  10317.X00167292  10317.X00167426 
    ##                1                1                1                2 
    ##  10317.X00167446  10317.X00167509  10317.X00167518  10317.X00167638 
    ##                1                1                1                1 
    ##  10317.X00167647  10317.X00167703  10317.X00167706  10317.X00167725 
    ##                1                2                1                1 
    ##  10317.X00167737  10317.X00167794  10317.X00167867  10317.X00167889 
    ##                1                1                2                1 
    ##  10317.X00167897  10317.X00167935  10317.X00168025  10317.X00168035 
    ##                1                2                1                1 
    ##  10317.X00168049  10317.X00168061  10317.X00168068  10317.X00168081 
    ##                1                1                1                1 
    ##  10317.X00168106  10317.X00168165  10317.X00168322  10317.X00168378 
    ##                1                1                1                2 
    ##  10317.X00168400  10317.X00168605  10317.X00168620  10317.X00168784 
    ##                2                1                1                2 
    ##  10317.X00168816  10317.X00168852  10317.X00168905  10317.X00168933 
    ##                2                1                2                2 
    ##  10317.X00168998  10317.X00169002  10317.X00169028  10317.X00169151 
    ##                1                1                1                1 
    ##  10317.X00169174  10317.X00169195  10317.X00169297  10317.X00169317 
    ##                2                1                1                1 
    ##  10317.X00169322  10317.X00169428  10317.X00169446  10317.X00169521 
    ##                1                1                2                1 
    ##  10317.X00169586  10317.X00169631  10317.X00169639  10317.X00169704 
    ##                2                1                1                1 
    ##  10317.X00169792  10317.X00169877  10317.X00169883  10317.X00170015 
    ##                1                1                1                1 
    ##  10317.X00170032  10317.X00170054  10317.X00170066  10317.X00170092 
    ##                1                1                1                2 
    ##  10317.X00170112  10317.X00170114  10317.X00170180  10317.X00170248 
    ##                1                1                1                1 
    ##  10317.X00170257  10317.X00170267  10317.X00170301  10317.X00170417 
    ##                1                1                1                1 
    ##  10317.X00170432  10317.X00170470  10317.X00170474  10317.X00170485 
    ##                2                1                1                1 
    ##  10317.X00170491  10317.X00170500  10317.X00170506  10317.X00170525 
    ##                2                1                1                1 
    ##  10317.X00170551  10317.X00170561  10317.X00170612  10317.X00170622 
    ##                1                1                1                1 
    ##  10317.X00170627  10317.X00170646  10317.X00170663  10317.X00170802 
    ##                1                2                1                1 
    ##  10317.X00170832  10317.X00170869  10317.X00170918  10317.X00170944 
    ##                2                1                2                1 
    ##  10317.X00170968  10317.X00170988  10317.X00171002  10317.X00171025 
    ##                1                1                1                1 
    ##  10317.X00171038  10317.X00171056  10317.X00171088  10317.X00171116 
    ##                1                2                1                1 
    ##  10317.X00171126  10317.X00171251  10317.X00171289  10317.X00171347 
    ##                1                1                1                1 
    ##  10317.X00171376  10317.X00171452  10317.X00171527  10317.X00171534 
    ##                2                2                1                1 
    ##  10317.X00171538  10317.X00171565  10317.X00171709  10317.X00171777 
    ##                1                1                1                1 
    ##  10317.X00171784  10317.X00171803  10317.X00171879  10317.X00171894 
    ##                1                1                1                1 
    ##  10317.X00172002  10317.X00172102  10317.X00172173  10317.X00172177 
    ##                1                1                1                1 
    ##  10317.X00172216  10317.X00172230  10317.X00172254  10317.X00172260 
    ##                2                1                1                1 
    ##  10317.X00172289  10317.X00172577  10317.X00172599  10317.X00172603 
    ##                1                2                1                1 
    ##  10317.X00172632  10317.X00172980  10317.X00172994  10317.X00173001 
    ##                1                1                1                1 
    ##  10317.X00173033  10317.X00173198  10317.X00173218  10317.X00173232 
    ##                1                1                1                1 
    ##  10317.X00173248  10317.X00173299  10317.X00173304  10317.X00173320 
    ##                1                1                1                1 
    ##  10317.X00173359  10317.X00173377  10317.X00173420  10317.X00173451 
    ##                1                1                1                1 
    ##  10317.X00173492  10317.X00173511  10317.X00173550  10317.X00173562 
    ##                2                2                2                1 
    ##  10317.X00173619  10317.X00173629  10317.X00173719  10317.X00173770 
    ##                1                1                1                1 
    ##  10317.X00173879  10317.X00173920  10317.X00173943  10317.X00173960 
    ##                3                1                1                1 
    ##  10317.X00173967  10317.X00173982  10317.X00174046  10317.X00174054 
    ##                1                1                2                1 
    ##  10317.X00174068  10317.X00174071  10317.X00174078  10317.X00174106 
    ##                1                1                1                1 
    ##  10317.X00174145  10317.X00174149  10317.X00174151  10317.X00174248 
    ##                2                1                1                1 
    ##  10317.X00174268  10317.X00174340  10317.X00174347  10317.X00174360 
    ##                1                1                1                1 
    ##  10317.X00174407  10317.X00174427  10317.X00174510  10317.X00174569 
    ##                1                1                1                1 
    ##  10317.X00174573  10317.X00174582  10317.X00174755  10317.X00174885 
    ##                2                1                1                1 
    ##  10317.X00174956  10317.X00174969  10317.X00174995  10317.X00175040 
    ##                1                1                1                1 
    ##  10317.X00175098  10317.X00175241  10317.X00175260  10317.X00175308 
    ##                1                2                1                1 
    ##  10317.X00175353  10317.X00175435  10317.X00175520  10317.X00175593 
    ##                2                1                1                1 
    ##  10317.X00175679  10317.X00175748  10317.X00175868  10317.X00182296 
    ##                1                1                1                1 
    ##  10317.X00191478  10317.X00191497  10317.X00191499  10317.X00191505 
    ##                1                1                1                1 
    ##  10317.X00191523  10317.X00191532  10317.X00191710  10317.X00191736 
    ##                1                1                1                1 
    ##  10317.X00191743  10317.X00191763  10317.X00191767  10317.X00191783 
    ##                1                1                1                1 
    ##  10317.X00191788  10317.X00191814  10317.X00191854  10317.X00191871 
    ##                1                1                1                1 
    ##  10317.X00191878  10317.X00191890  10317.X00191923  10317.X00191928 
    ##                1                1                1                1 
    ##  10317.X00191944  10317.X00191946  10317.X00191953  10317.X00192000 
    ##                1                1                1                1 
    ##  10317.X00192022  10317.X00192027  10317.X00192058  10317.X00192060 
    ##                1                1                1                1 
    ##  10317.X00192084  10317.X00192090  10317.X00192123  10317.X00192144 
    ##                1                1                1                1 
    ##  10317.X00192169  10317.X00192216  10317.X00192231  10317.X00192234 
    ##                1                1                1                1 
    ##  10317.X00192264  10317.X00192309  10317.X00192357  10317.X00192360 
    ##                1                1                1                1 
    ##  10317.X00192381  10317.X00192410  10317.X00192420  10317.X00192427 
    ##                1                1                1                1 
    ##  10317.X00192465  10317.X00192474  10317.X00192490  10317.X00192502 
    ##                1                1                1                1 
    ##  10317.X00192525  10317.X00192550  10317.X00192552  10317.X00192574 
    ##                1                1                1                1 
    ##  10317.X00192585  10317.X00192586  10317.X00192613  10317.X00192616 
    ##                1                1                1                1 
    ##  10317.X00192631  10317.X00192640  10317.X00192642  10317.X00192664 
    ##                1                1                1                1 
    ##  10317.X00192706  10317.X00192774  10317.X00192837  10317.X00192885 
    ##                1                1                1                1 
    ##  10317.X00192900  10317.X00192933 
    ##                1                1

``` r
# Create matched_samples column
metadata$matched_samples <- ifelse(metadata$sample_name %in% names(sample_counts[sample_counts >= 2]), "YES", "NO")

cat("\nMatched samples distribution:\n")
```

    ## 
    ## Matched samples distribution:

``` r
table(metadata$matched_samples)
```

    ## 
    ##   NO  YES 
    ## 1879  905

``` r
# Identify samples for removal: single samples excluding saliva and tongue coat
single_samples <- names(sample_counts[sample_counts == 1])
samples_to_remove <- metadata$sample_name %in% single_samples & 
                     !(metadata$Site %in% c("saliva", "tongue coat"))

cat("\nNumber of samples to be removed:", sum(samples_to_remove), "\n")
```

    ## 
    ## Number of samples to be removed: 1879

``` r
# Remove identified samples
metadata_clean <- metadata[!samples_to_remove, ]

cat("\nFinal sample distribution after cleaning:\n")
```

    ## 
    ## Final sample distribution after cleaning:

``` r
table(metadata_clean$sample_name)
```

    ## 
    ## 10317.000002102 10317.000002248 10317.000002251 10317.000002688 10317.000003224 
    ##               3               3               3               2               2 
    ## 10317.000003236 10317.000003557 10317.000003558 10317.000005899 10317.000010255 
    ##               3               2               3               2               3 
    ## 10317.000011233 10317.000012312 10317.000013183 10317.000013376 10317.000013967 
    ##               2               2               3               3               3 
    ## 10317.000014147 10317.000014889 10317.000014938 10317.000015017 10317.000015056 
    ##               2               2               3               2               2 
    ## 10317.000015181 10317.000015189 10317.000015281 10317.000015304 10317.000015344 
    ##               3               2               3               3               2 
    ## 10317.000015362 10317.000015454 10317.000015493 10317.000015504 10317.000015509 
    ##               3               3               3               3               3 
    ## 10317.000015561 10317.000016459 10317.000017141 10317.000017317 10317.000017318 
    ##               3               3               2               3               3 
    ## 10317.000017689 10317.000017757 10317.000017767 10317.000018053 10317.000018077 
    ##               2               2               3               2               2 
    ## 10317.000018079 10317.000020343 10317.000020447 10317.000020490 10317.000020502 
    ##               3               3               3               3               3 
    ## 10317.000020604 10317.000020695 10317.000020696 10317.000020707 10317.000020766 
    ##               3               3               3               3               3 
    ## 10317.000021246 10317.000022080 10317.000022130 10317.000022131 10317.000022140 
    ##               3               3               3               3               3 
    ## 10317.000022188 10317.000022446 10317.000022535 10317.000022695 10317.000022906 
    ##               3               3               3               3               3 
    ## 10317.000022911 10317.000022948 10317.000022970 10317.000023008 10317.000023010 
    ##               3               3               3               3               3 
    ## 10317.000023018 10317.000023039 10317.000023043 10317.000023045 10317.000023058 
    ##               3               3               3               2               2 
    ## 10317.000023062 10317.000023075 10317.000023127 10317.000023174 10317.000023175 
    ##               3               3               3               2               3 
    ## 10317.000023301 10317.000023302 10317.000023330 10317.000023875 10317.000023912 
    ##               3               3               3               2               3 
    ## 10317.000026043 10317.000026058 10317.000026509 10317.000026616 10317.000027350 
    ##               3               3               3               2               3 
    ## 10317.000027357 10317.000027728 10317.000028665 10317.000028764 10317.000028774 
    ##               3               3               2               3               2 
    ## 10317.000028844 10317.000028853 10317.000028876 10317.000028903 10317.000028917 
    ##               2               2               3               3               2 
    ## 10317.000028929 10317.000028941 10317.000029263 10317.000029616 10317.000031497 
    ##               2               2               2               2               2 
    ## 10317.000031877 10317.000032624 10317.000032637 10317.000032994 10317.000033106 
    ##               2               2               2               2               2 
    ## 10317.000033739 10317.000036171 10317.000036208 10317.000037537 10317.000038188 
    ##               4               2               2               2               2 
    ## 10317.000039931 10317.000043103 10317.000046268 10317.000046381 10317.000047831 
    ##               2               2               2               2               2 
    ## 10317.000047842 10317.000048269 10317.000048361 10317.000049133 10317.000049969 
    ##               2               2               2               2               2 
    ## 10317.000052371 10317.000053438 10317.000074414 10317.000082936 10317.000145614 
    ##               2               2               2               2               2 
    ## 10317.000145664 10317.000145667 10317.000145734 10317.000145769 10317.000145841 
    ##               2               2               2               2               2 
    ## 10317.000145947 10317.000146018 10317.000146250 10317.000146277 10317.000146411 
    ##               2               2               2               2               2 
    ## 10317.000146424 10317.000146438 10317.000146445 10317.000146532 10317.000146733 
    ##               2               2               2               2               2 
    ## 10317.000146735 10317.000147489 10317.000147498 10317.000147555 10317.000147720 
    ##               2               2               2               2               2 
    ## 10317.000147744 10317.000147842 10317.000147881 10317.000147891 10317.000147894 
    ##               2               2               2               2               2 
    ## 10317.000147897 10317.000147906 10317.000147920 10317.000147945 10317.000147956 
    ##               2               2               2               2               2 
    ## 10317.000147960 10317.000147963 10317.000147992 10317.000148025 10317.000148107 
    ##               2               2               2               2               2 
    ## 10317.000148136 10317.000148278 10317.000148280 10317.000148316 10317.000148407 
    ##               2               2               2               2               2 
    ## 10317.000148416 10317.000148418 10317.000148425 10317.000148725 10317.000148793 
    ##               2               2               2               2               2 
    ## 10317.000148800 10317.000148808 10317.000148824 10317.000148938 10317.000148974 
    ##               2               2               2               2               2 
    ## 10317.000148980 10317.000149018 10317.000149036 10317.000149039 10317.000149216 
    ##               2               2               2               2               2 
    ## 10317.000149220 10317.000149226 10317.000149315 10317.000149343 10317.000149348 
    ##               2               2               2               2               2 
    ## 10317.000149417 10317.000149424 10317.000149658 10317.000149660 10317.000149664 
    ##               2               2               2               2               2 
    ## 10317.000149882 10317.000149883 10317.000149975 10317.000149984 10317.000150020 
    ##               2               2               2               2               2 
    ## 10317.000150066 10317.000150090 10317.000150201 10317.000150222 10317.000150249 
    ##               2               2               2               2               2 
    ## 10317.000150300 10317.000150492 10317.000150501 10317.000150510 10317.000150516 
    ##               2               2               2               2               2 
    ## 10317.000150518 10317.000150521 10317.000150531 10317.000150626 10317.000150644 
    ##               2               2               2               2               2 
    ## 10317.000150704 10317.000150990 10317.000151037 10317.000151055 10317.000151118 
    ##               2               2               2               2               2 
    ## 10317.000151124 10317.000151145 10317.000151218 10317.000151220 10317.000151239 
    ##               2               2               2               2               2 
    ## 10317.000151256 10317.000151259 10317.000151275 10317.000151317 10317.000151319 
    ##               2               2               2               2               2 
    ## 10317.000151412 10317.000151445 10317.000151812 10317.000151827 10317.000151851 
    ##               2               2               2               2               3 
    ## 10317.000151865 10317.000151905 10317.000152253 10317.000152364 10317.000152373 
    ##               2               2               2               2               2 
    ## 10317.000152418 10317.000152426 10317.000152430 10317.000152474 10317.000152478 
    ##               2               2               2               2               2 
    ## 10317.000152765 10317.000152951 10317.000153062 10317.000153093 10317.000153168 
    ##               2               2               2               2               2 
    ## 10317.000154308 10317.000154338 10317.000154383 10317.000154386 10317.000154397 
    ##               2               2               2               2               2 
    ## 10317.000154415 10317.000154427 10317.000154443 10317.000154457 10317.000154460 
    ##               2               2               2               2               2 
    ## 10317.000154482 10317.000154509 10317.000154512 10317.000154554 10317.000154559 
    ##               2               2               2               2               2 
    ## 10317.000154571 10317.000154575 10317.000154577 10317.000154596 10317.000154605 
    ##               2               2               2               2               2 
    ## 10317.000154614 10317.000154631 10317.000154656 10317.000154662 10317.000154679 
    ##               2               2               2               2               2 
    ## 10317.000154724 10317.000154742 10317.000154785 10317.000154791 10317.000154893 
    ##               2               2               2               2               2 
    ## 10317.000154896 10317.000154923 10317.000155016 10317.000155039 10317.000155043 
    ##               2               2               2               2               2 
    ## 10317.000155045 10317.000155058 10317.000155060 10317.000155064 10317.000155102 
    ##               2               2               2               3               2 
    ## 10317.000155109 10317.000155121 10317.000155132 10317.000155141 10317.000155150 
    ##               2               2               2               2               2 
    ## 10317.000155184 10317.000155196 10317.000155205 10317.000155208 10317.000155235 
    ##               2               2               2               2               2 
    ## 10317.000155238 10317.000155252 10317.000155268 10317.000155295 10317.000155298 
    ##               2               2               2               2               2 
    ## 10317.000155318 10317.000155372 10317.000155388 10317.000155396 10317.000155403 
    ##               2               2               2               2               2 
    ## 10317.000155408 10317.000155436 10317.000155439 10317.000155595 10317.000155600 
    ##               2               2               2               2               2 
    ## 10317.000155633 10317.000155660 10317.000155852 10317.000155859 10317.000155862 
    ##               2               2               2               2               2 
    ## 10317.000155864 10317.000155871 10317.000155907 10317.000155915 10317.000155922 
    ##               2               2               2               2               2 
    ## 10317.000155930 10317.000155963 10317.000155991 10317.000155993 10317.000156009 
    ##               2               2               2               2               2 
    ## 10317.000156011 10317.000156035 10317.000156080 10317.000156083 10317.000156095 
    ##               2               2               2               2               2 
    ## 10317.000156107 10317.000156123 10317.000156125 10317.000157472 10317.000157527 
    ##               2               2               2               2               2 
    ## 10317.000157551 10317.000157566 10317.000157593 10317.000157617 10317.000157680 
    ##               2               2               2               2               2 
    ## 10317.000157778 10317.000157812 10317.000157824 10317.000157835 10317.000157844 
    ##               2               2               2               2               2 
    ## 10317.000157929 10317.000157935 10317.000157986 10317.000157994 10317.000158018 
    ##               2               2               2               2               2 
    ## 10317.000158031 10317.000158036 10317.000158042 10317.000158064 10317.000158073 
    ##               2               2               2               2               2 
    ## 10317.000158090 10317.X00160990 10317.X00161003 10317.X00161177 10317.X00161471 
    ##               2               2               2               2               2 
    ## 10317.X00161645 10317.X00161728 10317.X00162145 10317.X00162295 10317.X00162298 
    ##               2               2               2               2               2 
    ## 10317.X00162439 10317.X00162795 10317.X00162916 10317.X00163024 10317.X00163410 
    ##               2               2               2               2               2 
    ## 10317.X00163859 10317.X00164100 10317.X00164607 10317.X00164793 10317.X00164869 
    ##               2               2               2               2               2 
    ## 10317.X00164945 10317.X00165204 10317.X00166293 10317.X00166568 10317.X00166575 
    ##               2               2               3               2               2 
    ## 10317.X00166596 10317.X00166778 10317.X00167426 10317.X00167703 10317.X00167867 
    ##               2               3               2               2               2 
    ## 10317.X00167935 10317.X00168378 10317.X00168400 10317.X00168784 10317.X00168816 
    ##               2               2               2               2               2 
    ## 10317.X00168905 10317.X00168933 10317.X00169174 10317.X00169446 10317.X00169586 
    ##               2               2               2               2               2 
    ## 10317.X00170092 10317.X00170432 10317.X00170491 10317.X00170646 10317.X00170832 
    ##               2               2               2               2               2 
    ## 10317.X00170918 10317.X00171056 10317.X00171376 10317.X00171452 10317.X00172216 
    ##               2               2               2               2               2 
    ## 10317.X00172577 10317.X00173492 10317.X00173511 10317.X00173550 10317.X00173879 
    ##               2               2               2               2               3 
    ## 10317.X00174046 10317.X00174145 10317.X00174573 10317.X00175241 10317.X00175353 
    ##               2               2               2               2               2

``` r
table(metadata_clean$Site)
```

    ## 
    ## human oral metagenome 
    ##                   905

``` r
table(metadata_clean$matched_samples)
```

    ## 
    ## YES 
    ## 905

``` r
# Use cleaned dataset for subsequent analysis
metadata <- metadata_clean
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 905"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 415"

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

    ##     BioProject        Run     sample_name Age     BMI         Sex Region
    ## 46  PRJEB11419 ERR1077608 10317.000014147  NA 18.5-30           F     V4
    ## 52  PRJEB11419 ERR1077622 10317.000015056  NA 18.5-30 unspecified     V4
    ## 58  PRJEB11419 ERR1077632 10317.000015189  NA 18.5-30           M     V4
    ## 79  PRJEB11419 ERR1077679 10317.000017141  NA 18.5-30           F     V4
    ## 134 PRJEB11419 ERR1073499 10317.000002688  NA 18.5-30           M     V4
    ## 147 PRJEB11419 ERR1073522 10317.000003224  NA 18.5-30           F     V4
    ##           Platform       Country                  Site          public_time
    ## 46  Illumina MiSeq  Not provided human oral metagenome 2015-10-28T00:00:00Z
    ## 52  Illumina MiSeq United States human oral metagenome 2015-10-28T00:00:00Z
    ## 58  Illumina MiSeq United States human oral metagenome 2015-10-28T00:00:00Z
    ## 79  Illumina MiSeq United States human oral metagenome 2015-10-28T00:00:00Z
    ## 134 Illumina MiSeq United States human oral metagenome 2015-10-28T00:00:00Z
    ## 147 Illumina MiSeq United States human oral metagenome 2015-10-28T00:00:00Z
    ##     PCR_primers
    ## 46    515F-806R
    ## 52    515F-806R
    ## 58    515F-806R
    ## 79    515F-806R
    ## 134   515F-806R
    ## 147   515F-806R
    ##                                                                      title
    ## 46  American Gut: an Open Platform for Citizen Science Microbiome Research
    ## 52  American Gut: an Open Platform for Citizen Science Microbiome Research
    ## 58  American Gut: an Open Platform for Citizen Science Microbiome Research
    ## 79  American Gut: an Open Platform for Citizen Science Microbiome Research
    ## 134 American Gut: an Open Platform for Citizen Science Microbiome Research
    ## 147 American Gut: an Open Platform for Citizen Science Microbiome Research
    ##     Disease_study Status matched_samples
    ## 46             NO     HC             YES
    ## 52             NO     HC             YES
    ## 58             NO     HC             YES
    ## 79             NO     HC             YES
    ## 134            NO     HC             YES
    ## 147            NO     HC             YES

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
