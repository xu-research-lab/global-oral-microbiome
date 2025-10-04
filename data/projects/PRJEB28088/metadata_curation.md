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
setwd("D:/oral_project/PRJEB28088")

# Read the initial metadata files
initial_meta<-read.csv("SraRunTable.csv",header = T)

# Combine the two metadata files
#meta <- rbind(initial_meta, initial_meta1)

#Remove entries where Organism is "gut metagenome" (non-oral samples)
#initial_meta <- subset(initial_meta,isolation_source!= "feces")

# Check the number of samples containing "HT" in Sample.Name (e.g., healthy controls)
#nrow(subset(meta, grepl("HT", Sample.Name)))

# Display the first few rows of the metadata
head(initial_meta)
```

    ##      abi_left1    abi_left2    abi_left3   abi_right1   abi_right2   abi_right3
    ## 1 not provided not provided not provided not provided not provided not provided
    ## 2 not provided not provided not provided not provided not provided not provided
    ## 3 not provided not provided not provided not provided not provided not provided
    ## 4 not provided not provided not provided not provided not provided not provided
    ## 5 not provided not provided not provided not provided not provided not provided
    ## 6 not provided not provided not provided not provided not provided not provided
    ##          Run AGE ancylostomide_eggs anonymized_name arm_circumference
    ## 1 ERR2728158  55       not provided            15-1      not provided
    ## 2 ERR2728159  34       not provided         15-1053      not provided
    ## 3 ERR2728160  34       not provided         15-1055      not provided
    ## 4 ERR2728161  44       not provided         15-1062      not provided
    ## 5 ERR2728162   3       not provided         15-1081      not provided
    ## 6 ERR2728163  34       not provided         15-1083      not provided
    ##     arm_length ascaris_lumbricoides_eggs Assay.Type AvgSpotLen barcode..exp.
    ## 1 not provided              not provided      OTHER        150  GTCATGCTCCAG
    ## 2 not provided              not provided      OTHER        150  CAATCGGCTTGC
    ## 3 not provided              not provided      OTHER        150  TACAGTTACGCG
    ## 4 not provided              not provided      OTHER        150  CGGCCTAAGTTC
    ## 5 not provided              not provided      OTHER        150  TAGCGCGAACTT
    ## 6 not provided              not provided      OTHER        150  TCGAGCCGATCT
    ##     Bases BioProject    BioSample   Birth_date blastocystis_spp          BMI
    ## 1 5415214 PRJEB28088 SAMEA4820302 not provided     not provided not provided
    ## 2 4130375 PRJEB28088 SAMEA4820303 not provided     not provided not provided
    ## 3 4814351 PRJEB28088 SAMEA4820304 not provided     not provided not provided
    ## 4 2887503 PRJEB28088 SAMEA4820305 not provided     not provided not provided
    ## 5 3352410 PRJEB28088 SAMEA4820306 not provided     not provided not provided
    ## 6 5911123 PRJEB28088 SAMEA4820307 not provided     not provided not provided
    ##             bp   Bytes          c_i   cavi_left1   cavi_left2   cavi_left3
    ## 1 not provided 3283154 not provided not provided not provided not provided
    ## 2 not provided 2504788 not provided not provided not provided not provided
    ## 3 not provided 2667913 not provided not provided not provided not provided
    ## 4 not provided 1773606 not provided not provided not provided not provided
    ## 5 not provided 2028762 not provided not provided not provided not provided
    ## 6 not provided 3638991 not provided not provided not provided not provided
    ##    cavi_right1  cavi_right2  cavi_right3 center_name..exp.
    ## 1 not provided not provided not provided            NYULMC
    ## 2 not provided not provided not provided            NYULMC
    ## 3 not provided not provided not provided            NYULMC
    ## 4 not provided not provided not provided            NYULMC
    ## 5 not provided not provided not provided            NYULMC
    ## 6 not provided not provided not provided            NYULMC
    ##                                          Center.Name center_project_name..exp.
    ## 1 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT              Bolivar-Diet
    ## 2 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT              Bolivar-Diet
    ## 3 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT              Bolivar-Diet
    ## 4 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT              Bolivar-Diet
    ## 5 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT              Bolivar-Diet
    ## 6 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT              Bolivar-Diet
    ##   chilomastix_mesnili_cysts  cholesterol collection_timestamp Consent
    ## 1              not provided not provided           2015-09-30  public
    ## 2              not provided not provided           2015-10-12  public
    ## 3              not provided not provided           2015-10-12  public
    ## 4              not provided not provided           2015-10-13  public
    ## 5              not provided not provided           2015-10-15  public
    ## 6              not provided not provided           2015-10-15  public
    ##             COUNTRY    cuff_size DATASTORE.filetype DATASTORE.provider
    ## 1 Restricted access not provided   fastq,run.zq,sra     ena,gs,ncbi,s3
    ## 2 Restricted access not provided   fastq,run.zq,sra     ena,gs,ncbi,s3
    ## 3 Restricted access not provided   fastq,run.zq,sra     ena,gs,ncbi,s3
    ## 4 Restricted access not provided   fastq,run.zq,sra     ena,gs,ncbi,s3
    ## 5 Restricted access not provided   fastq,run.zq,sra     ena,gs,ncbi,s3
    ## 6 Restricted access not provided   fastq,run.zq,sra     ena,gs,ncbi,s3
    ##                           DATASTORE.region days_of_study
    ## 1 ena,gs.us-east1,ncbi.public,s3.us-east-1             1
    ## 2 ena,gs.us-east1,ncbi.public,s3.us-east-1            13
    ## 3 ena,gs.us-east1,ncbi.public,s3.us-east-1            13
    ## 4 ena,gs.us-east1,ncbi.public,s3.us-east-1            14
    ## 5 ena,gs.us-east1,ncbi.public,s3.us-east-1            16
    ## 6 ena,gs.us-east1,ncbi.public,s3.us-east-1            16
    ##   days_relative_to_expedition  diastolic_1  diastolic_2   diatolic_3
    ## 1                  Kanarakuni not provided not provided not provided
    ## 2                  Kanarakuni not provided not provided not provided
    ## 3                  Kanarakuni not provided not provided not provided
    ## 4               +1 Kanarakuni not provided not provided not provided
    ## 5               +4 Kanarakuni not provided not provided not provided
    ## 6               +3 Kanarakuni not provided not provided not provided
    ##   dientamoeba_fragilis_cyst     diet_24h
    ## 1              not provided not provided
    ## 2              not provided not provided
    ## 3              not provided not provided
    ## 4              not provided not provided
    ## 5              not provided not provided
    ## 6              not provided not provided
    ##                                                                                                                                                                                                                                                              Diet
    ## 1 Breakfast- yukuta/ aji/ casabe/ boiled fish -coporo-/ fried fish/ peccary stew; Morning snack- guava/ sugar cane; Lunch- yukuta/ aji/ casabe/ boiled fish. -coporo-/ roasted peccary; Afternoon snack- banana; Dinner- yukuta/ aji/ casabe/ boiled fish -cop...
    ## 2                                                                     Breakfast- arepa/ aji/ cereal/ tuna/ milk/ coffee/ pumpkin soup/ fried and boiled lowland paca; Lunch- aji/ casabe/ pastas/ tuna/ nestea; Dinner- rice/ aji/ casabe/ tuna/ sardine/ nestea.
    ## 3                                                                     Breakfast- arepa/ aji/ cereal/ tuna/ milk/ coffee/ pumpkin soup/ fried and boiled lowland paca; Lunch- aji/ casabe/ pastas/ tuna/ nestea; Dinner- rice/ aji/ casabe/ tuna/ sardine/ nestea.
    ## 4                                                                                                                                                                                                                                                    not provided
    ## 5                                                                                                                                                                                                                                                    not provided
    ## 6                                                                                                                                                                                                                                                    not provided
    ##   dna_extracted elevation ENA.FIRST.PUBLIC..run. ENA.LAST.UPDATE..run.
    ## 1          TRUE        10             2018-08-11            2018-11-16
    ## 2          TRUE        10             2018-08-11            2018-11-16
    ## 3          TRUE        10             2018-08-11            2018-11-16
    ## 4          TRUE        10             2018-08-11            2018-11-16
    ## 5          TRUE        10             2018-08-11            2018-11-16
    ## 6          TRUE        10             2018-08-11            2018-11-16
    ##   endolimax_nana_cysts entamoeba_coli_cysts entamoeba_coli_trofozoites
    ## 1         not provided         not provided               not provided
    ## 2         not provided         not provided               not provided
    ## 3         not provided         not provided               not provided
    ## 4         not provided         not provided               not provided
    ## 5         not provided         not provided               not provided
    ## 6         not provided         not provided               not provided
    ##   entamoeba_histolytica_dispar entamoeba_histolytica_trofozoites
    ## 1                 not provided                      not provided
    ## 2                 not provided                      not provided
    ## 3                 not provided                      not provided
    ## 4                 not provided                      not provided
    ## 5                 not provided                      not provided
    ## 6                 not provided                      not provided
    ##   enterobius_vermicularis_eggs              env_biome              env_feature
    ## 1                 not provided dense settlement biome human-associated habitat
    ## 2                 not provided dense settlement biome human-associated habitat
    ## 3                 not provided dense settlement biome human-associated habitat
    ## 4                 not provided dense settlement biome human-associated habitat
    ## 5                 not provided dense settlement biome human-associated habitat
    ## 6                 not provided dense settlement biome human-associated habitat
    ##   env_material env_package    ethnicity Experiment  external_id    Family_Id
    ## 1        feces   human-gut not provided ERX2741403 SAMEA4820302 not provided
    ## 2        feces   human-gut not provided ERX2741404 SAMEA4820303 not provided
    ## 3        feces   human-gut not provided ERX2741405 SAMEA4820304 not provided
    ## 4        feces   human-gut not provided ERX2741406 SAMEA4820305 not provided
    ## 5        feces   human-gut not provided ERX2741407 SAMEA4820306 not provided
    ## 6        feces   human-gut not provided ERX2741408 SAMEA4820307 not provided
    ##   family_relations fitbit_bit_no fitbit_highly_active_minutes
    ## 1     not provided  not provided                 not provided
    ## 2     not provided  not provided                 not provided
    ## 3     not provided  not provided                 not provided
    ## 4     not provided  not provided                 not provided
    ## 5     not provided  not provided                 not provided
    ## 6     not provided  not provided                 not provided
    ##   fitbit_miles_walked fitbit_received fitbit_returned fitbit_steps gender
    ## 1        not provided    not provided    not provided not provided female
    ## 2        not provided    not provided    not provided not provided female
    ## 3        not provided    not provided    not provided not provided   male
    ## 4        not provided    not provided    not provided not provided female
    ## 5        not provided    not provided    not provided not provided   male
    ## 6        not provided    not provided    not provided not provided female
    ##   geo_loc_name_country geo_loc_name_country_continent   geo_loc_name
    ## 1         uncalculated                   uncalculated Sourth America
    ## 2         uncalculated                   uncalculated Sourth America
    ## 3         uncalculated                   uncalculated Sourth America
    ## 4         uncalculated                   uncalculated Sourth America
    ## 5         uncalculated                   uncalculated Sourth America
    ## 6         uncalculated                   uncalculated Sourth America
    ##   giardia_lamblia_cysts giardia_lamblia_trofozoites      Glucose   hemoglobin
    ## 1          not provided                not provided not provided not provided
    ## 2          not provided                not provided not provided not provided
    ## 3          not provided                not provided not provided not provided
    ## 4          not provided                not provided not provided not provided
    ## 5          not provided                not provided not provided not provided
    ## 6          not provided                not provided not provided not provided
    ##   Host_Age host_age_units host_body_mass_index host_common_name  host_height
    ## 1       55          years         not provided            human not provided
    ## 2       34          years         not provided            human not provided
    ## 3       34          years         not provided            human not provided
    ## 4       44          years         not provided            human not provided
    ## 5        3          years         not provided            human not provided
    ## 6       34          years         not provided            human not provided
    ##   host_height_units host_scientific_name host_subject_id host_taxid
    ## 1                cm         Homo sapiens              A7       9606
    ## 2                cm         Homo sapiens              A5       9606
    ## 3                cm         Homo sapiens              A6       9606
    ## 4                cm         Homo sapiens              A2       9606
    ## 5                cm         Homo sapiens              A3       9606
    ## 6                cm         Homo sapiens              A5       9606
    ##    host_weight host_weight_units     house_id hymenolepis_nana_eggs
    ## 1 not provided                kg not provided          not provided
    ## 2 not provided                kg not provided          not provided
    ## 3 not provided                kg not provided          not provided
    ## 4 not provided                kg not provided          not provided
    ## 5 not provided                kg not provided          not provided
    ## 6 not provided                kg not provided          not provided
    ##   INSDC_center_name   INSDC_first_public    INSDC_last_update INSDC_status
    ## 1            UCSDMI 2018-08-11T17:03:01Z 2018-08-06T22:41:45Z       public
    ## 2            UCSDMI 2018-08-11T17:03:01Z 2018-08-06T22:41:45Z       public
    ## 3            UCSDMI 2018-08-11T17:03:01Z 2018-08-06T22:41:45Z       public
    ## 4            UCSDMI 2018-08-11T17:03:01Z 2018-08-06T22:41:45Z       public
    ## 5            UCSDMI 2018-08-11T17:03:01Z 2018-08-06T22:41:45Z       public
    ## 6            UCSDMI 2018-08-11T17:03:01Z 2018-08-06T22:41:45Z       public
    ##       Instrument iodamoeba_butschlii_cysts latitude   Library.Name
    ## 1 Illumina MiSeq              not provided   40.737    11874.X15.1
    ## 2 Illumina MiSeq              not provided   40.737 11874.X15.1053
    ## 3 Illumina MiSeq              not provided   40.737 11874.X15.1055
    ## 4 Illumina MiSeq              not provided   40.737 11874.X15.1062
    ## 5 Illumina MiSeq              not provided   40.737 11874.X15.1081
    ## 6 Illumina MiSeq              not provided   40.737 11874.X15.1083
    ##   LibraryLayout LibrarySelection LibrarySource linker..exp. longitude
    ## 1        SINGLE              PCR   METAGENOMIC           CC  -73.9767
    ## 2        SINGLE              PCR   METAGENOMIC           CC  -73.9767
    ## 3        SINGLE              PCR   METAGENOMIC           CC  -73.9767
    ## 4        SINGLE              PCR   METAGENOMIC           CC  -73.9767
    ## 5        SINGLE              PCR   METAGENOMIC           CC  -73.9767
    ## 6        SINGLE              PCR   METAGENOMIC           CC  -73.9767
    ##                name observations   onchocerca             Organism
    ## 1 Restricted access not provided not provided human gut metagenome
    ## 2 Restricted access not provided not provided human gut metagenome
    ## 3 Restricted access not provided not provided human gut metagenome
    ## 4 Restricted access not provided not provided human gut metagenome
    ## 5 Restricted access not provided not provided human gut metagenome
    ## 6 Restricted access not provided not provided human gut metagenome
    ##                                   pcr_primers..exp. peak_inflation_level
    ## 1 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT         not provided
    ## 2 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT         not provided
    ## 3 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT         not provided
    ## 4 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT         not provided
    ## 5 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT         not provided
    ## 6 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT         not provided
    ##   physical_specimen_location physical_specimen_remaining Platform
    ## 1            Rutgers_MGDBLab                       FALSE ILLUMINA
    ## 2            Rutgers_MGDBLab                       FALSE ILLUMINA
    ## 3            Rutgers_MGDBLab                       FALSE ILLUMINA
    ## 4            Rutgers_MGDBLab                       FALSE ILLUMINA
    ## 5            Rutgers_MGDBLab                       FALSE ILLUMINA
    ## 6            Rutgers_MGDBLab                       FALSE ILLUMINA
    ##           primer..exp. pulse_obliteration_pressure        pulse
    ## 1 GGACTACHVGGGTWTCTAAT                not provided not provided
    ## 2 GGACTACHVGGGTWTCTAAT                not provided not provided
    ## 3 GGACTACHVGGGTWTCTAAT                not provided not provided
    ## 4 GGACTACHVGGGTWTCTAAT                not provided not provided
    ## 5 GGACTACHVGGGTWTCTAAT                not provided not provided
    ## 6 GGACTACHVGGGTWTCTAAT                not provided not provided
    ##            ReleaseDate run_center..exp. run_date..exp.          create_date
    ## 1 2018-08-12T00:00:00Z           NYULMC        2016-09 2018-08-12T17:27:00Z
    ## 2 2018-08-12T00:00:00Z           NYULMC        2016-09 2018-08-12T17:27:00Z
    ## 3 2018-08-12T00:00:00Z           NYULMC        2016-09 2018-08-12T17:27:00Z
    ## 4 2018-08-12T00:00:00Z           NYULMC        2016-09 2018-08-12T17:28:00Z
    ## 5 2018-08-12T00:00:00Z           NYULMC        2016-09 2018-08-12T17:28:00Z
    ## 6 2018-08-12T00:00:00Z           NYULMC        2016-09 2018-08-12T17:28:00Z
    ##   version              run_prefix..exp.  Sample.Name
    ## 1       1 lane1_NoIndex_L001_R2_001-13A SAMEA4820302
    ## 2       1 lane1_NoIndex_L001_R2_001-13A SAMEA4820303
    ## 3       1 lane1_NoIndex_L001_R2_001-13A SAMEA4820304
    ## 4       1 lane1_NoIndex_L001_R2_001-13A SAMEA4820305
    ## 5       1 lane1_NoIndex_L001_R2_001-13A SAMEA4820306
    ## 6       1 lane1_NoIndex_L001_R2_001-13A SAMEA4820307
    ##                      Sample_name sample_type samplegroup  seqs
    ## 1    qiita_sid_11874:11874.X15.1       feces    Visitors 31972
    ## 2 qiita_sid_11874:11874.X15.1053       feces    Visitors 24311
    ## 3 qiita_sid_11874:11874.X15.1055       feces    Visitors 29369
    ## 4 qiita_sid_11874:11874.X15.1062       feces    Visitors 16980
    ## 5 qiita_sid_11874:11874.X15.1081       feces    Visitors 20288
    ## 6 qiita_sid_11874:11874.X15.1083       feces    Visitors 35556
    ##     sequencing_meth..exp. sequencing_run    sex SRA.Study
    ## 1 sequencing by synthesis        RUN_13A female ERP110254
    ## 2 sequencing by synthesis        RUN_13A female ERP110254
    ## 3 sequencing by synthesis        RUN_13A   male ERP110254
    ## 4 sequencing by synthesis        RUN_13A female ERP110254
    ## 5 sequencing by synthesis        RUN_13A   male ERP110254
    ## 6 sequencing by synthesis        RUN_13A female ERP110254
    ##   strongyloides_stercolaris_larvas subject_id                   Submitter_Id
    ## 1                     not provided         A7    qiita_sid_11874:11874.X15.1
    ## 2                     not provided         A5 qiita_sid_11874:11874.X15.1053
    ## 3                     not provided         A6 qiita_sid_11874:11874.X15.1055
    ## 4                     not provided         A2 qiita_sid_11874:11874.X15.1062
    ## 5                     not provided         A3 qiita_sid_11874:11874.X15.1081
    ## 6                     not provided         A5 qiita_sid_11874:11874.X15.1083
    ##    supine_dbp1  supine_dbp2  supine_dbp3  supine_sbp1  supine_sbp2  supine_sbp3
    ## 1 not provided not provided not provided not provided not provided not provided
    ## 2 not provided not provided not provided not provided not provided not provided
    ## 3 not provided not provided not provided not provided not provided not provided
    ## 4 not provided not provided not provided not provided not provided not provided
    ## 5 not provided not provided not provided not provided not provided not provided
    ## 6 not provided not provided not provided not provided not provided not provided
    ##     systolic_1
    ## 1 not provided
    ## 2 not provided
    ## 3 not provided
    ## 4 not provided
    ## 5 not provided
    ## 6 not provided

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
  sample_name = initial_meta$host_subject_id,
  Age = initial_meta$Host_Age,
  BMI = initial_meta$host_body_mass_index,
  Sex=initial_meta$gender,
  Region = "V4",
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  Site = initial_meta$env_material,
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
metadata$title <- "Changes in the Gut Microbiota of Urban Subjects during an Immersion in the Traditional Diet and Lifestyle of a Rainforest Village"
#metadata<-subset(uniform_meta,Site!="feces")
metadata$Disease_study <- "NO"
# Assign disease status based on sample_name patterns
metadata$Status <- NA
#extract subject id by remove prefix of sample_name
#metadata$sample_name <- gsub("[FSP].*", "", metadata$sample_name)
#get subjects id with Saliva available
keep_subjects<-subset(metadata,Site %in% c("saliva"))
#remove subjects without Saliva samples
metadata<-subset(metadata,sample_name %in% keep_subjects$sample_name)
#get subjects sampling counts
sample_counts <- table(metadata$sample_name)
# Create matched_samples column, any subjects have >=2 samples will be label "YES", otherwise "NO"
metadata$matched_samples <- ifelse(metadata$sample_name %in% names(sample_counts[sample_counts >= 2]), "YES", "NO")
metadata$BMI<-gsub("not provided",NA,metadata$BMI)
metadata$BMI<-as.numeric(metadata$BMI)
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 572"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 86"

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

    ##   BioProject        Run sample_name Age BMI    Sex Region       Platform
    ## 1 PRJEB28088 ERR2728158          A7  55  NA female     V4 Illumina MiSeq
    ## 2 PRJEB28088 ERR2728159          A5  34  NA female     V4 Illumina MiSeq
    ## 3 PRJEB28088 ERR2728160          A6  34  NA   male     V4 Illumina MiSeq
    ## 4 PRJEB28088 ERR2728161          A2  44  NA female     V4 Illumina MiSeq
    ## 5 PRJEB28088 ERR2728162          A3   3  NA   male     V4 Illumina MiSeq
    ## 6 PRJEB28088 ERR2728163          A5  34  NA female     V4 Illumina MiSeq
    ##        Country  Site          public_time PCR_primers
    ## 1 uncalculated feces 2018-08-12T00:00:00Z          NA
    ## 2 uncalculated feces 2018-08-12T00:00:00Z          NA
    ## 3 uncalculated feces 2018-08-12T00:00:00Z          NA
    ## 4 uncalculated feces 2018-08-12T00:00:00Z          NA
    ## 5 uncalculated feces 2018-08-12T00:00:00Z          NA
    ## 6 uncalculated feces 2018-08-12T00:00:00Z          NA
    ##                                                                                                                               title
    ## 1 Changes in the Gut Microbiota of Urban Subjects during an Immersion in the Traditional Diet and Lifestyle of a Rainforest Village
    ## 2 Changes in the Gut Microbiota of Urban Subjects during an Immersion in the Traditional Diet and Lifestyle of a Rainforest Village
    ## 3 Changes in the Gut Microbiota of Urban Subjects during an Immersion in the Traditional Diet and Lifestyle of a Rainforest Village
    ## 4 Changes in the Gut Microbiota of Urban Subjects during an Immersion in the Traditional Diet and Lifestyle of a Rainforest Village
    ## 5 Changes in the Gut Microbiota of Urban Subjects during an Immersion in the Traditional Diet and Lifestyle of a Rainforest Village
    ## 6 Changes in the Gut Microbiota of Urban Subjects during an Immersion in the Traditional Diet and Lifestyle of a Rainforest Village
    ##   Disease_study Status matched_samples
    ## 1            NO     NA             YES
    ## 2            NO     NA             YES
    ## 3            NO     NA             YES
    ## 4            NO     NA             YES
    ## 5            NO     NA             YES
    ## 6            NO     NA             YES

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
