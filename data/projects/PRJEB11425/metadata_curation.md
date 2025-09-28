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
setwd("D:/oral_project/PRJEB11425")

# Read the initial metadata files
initial_meta<-read.csv("SraRunTable.csv",header = T)

# Combine the two metadata files
#meta <- rbind(initial_meta, initial_meta1)

#Remove entries where Organism is "gut metagenome" (non-oral samples)
#initial_meta <- subset(initial_meta,env_material== "saliva")

# Check the number of samples containing "HT" in Sample.Name (e.g., healthy controls)
#nrow(subset(meta, grepl("HT", Sample.Name)))

# Display the first few rows of the metadata
head(initial_meta)
```

    ##          Run age_unit age_units anonymized_name Assay.Type AvgSpotLen
    ## 1 ERR1076083    years     years            6143      OTHER        151
    ## 2 ERR1076084    years     years            6144      OTHER        151
    ## 3 ERR1076085    years     years            6145      OTHER        151
    ## 4 ERR1076086    years     years            6146      OTHER        151
    ## 5 ERR1076087    years     years            6147      OTHER        151
    ## 6 ERR1076088    years     years            6149      OTHER        151
    ##   barcode..exp.   Bases BioProject    BioSample   Bytes center_name..exp.
    ## 1  TACCGCTTCTTC 1065909 PRJEB11425 SAMEA3611562  838881              CCME
    ## 2  ACGAGACTGATT 5142456 PRJEB11425 SAMEA3611563 3238501              CCME
    ## 3  GTCGAATTTGCG 3913316 PRJEB11425 SAMEA3611564 2439267              CCME
    ## 4  GTGCAACCAATC 2813432 PRJEB11425 SAMEA3611565 1967312              CCME
    ## 5  ATGGGTTCCGTC 5658272 PRJEB11425 SAMEA3611566 3375152              CCME
    ## 6  GATCACGAGAGG 4977564 PRJEB11425 SAMEA3611567 3195861              CCME
    ##                                          Center.Name center_project_name..exp.
    ## 1 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT            ICU microbiome
    ## 2 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT            ICU microbiome
    ## 3 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT            ICU microbiome
    ## 4 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT            ICU microbiome
    ## 5 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT            ICU microbiome
    ## 6 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT            ICU microbiome
    ##   collection_timestamp           common_name Consent COUNTRY DATASTORE.filetype
    ## 1               1/1/59 human oral metagenome  public     USA   fastq,run.zq,sra
    ## 2               1/1/59 human skin metagenome  public     USA   fastq,run.zq,sra
    ## 3               1/1/59  human gut metagenome  public     USA   fastq,run.zq,sra
    ## 4               1/1/59 human oral metagenome  public     USA   fastq,run.zq,sra
    ## 5               1/1/59 human skin metagenome  public     USA   fastq,run.zq,sra
    ## 6               1/1/59 human oral metagenome  public     USA   fastq,run.zq,sra
    ##   DATASTORE.provider                         DATASTORE.region Depth
    ## 1     ena,gs,ncbi,s3 ena,gs.us-east1,ncbi.public,s3.us-east-1     0
    ## 2     ena,gs,ncbi,s3 ena,gs.us-east1,ncbi.public,s3.us-east-1     0
    ## 3     ena,gs,ncbi,s3 ena,gs.us-east1,ncbi.public,s3.us-east-1     0
    ## 4     ena,gs,ncbi,s3 ena,gs.us-east1,ncbi.public,s3.us-east-1     0
    ## 5     ena,gs,ncbi,s3 ena,gs.us-east1,ncbi.public,s3.us-east-1     0
    ## 6     ena,gs,ncbi,s3 ena,gs.us-east1,ncbi.public,s3.us-east-1     0
    ##   dna_extracted elevation emp_status..exp. ENA.FIRST.PUBLIC..run.
    ## 1          True       193              EMP             2015-10-15
    ## 2          True       193              EMP             2015-10-15
    ## 3          True       193              EMP             2015-10-15
    ## 4          True       193              EMP             2015-10-15
    ## 5          True       193              EMP             2015-10-15
    ## 6          True       193              EMP             2015-10-15
    ##   ENA_first_public ENA.LAST.UPDATE..run. ENA_last_update   env_biome
    ## 1       2015-10-15            2018-11-16      2016-05-11 urban biome
    ## 2       2015-10-15            2018-11-16      2016-05-11 urban biome
    ## 3       2015-10-15            2018-11-16      2016-05-11 urban biome
    ## 4       2015-10-15            2018-11-16      2016-05-11 urban biome
    ## 5       2015-10-15            2018-11-16      2016-05-11 urban biome
    ## 6       2015-10-15            2018-11-16      2016-05-11 urban biome
    ##                env_feature env_material  env_matter Env_package
    ## 1 human-associated habitat       saliva ENVO:saliva  human-oral
    ## 2 human-associated habitat        sebum  ENVO:sebum  human-skin
    ## 3 human-associated habitat        feces  ENVO:feces   human-gut
    ## 4 human-associated habitat       saliva ENVO:saliva  human-oral
    ## 5 human-associated habitat        sebum  ENVO:sebum  human-skin
    ## 6 human-associated habitat       saliva ENVO:saliva  human-oral
    ##   experiment_center..exp. Experiment experiment_title..exp.  external_id
    ## 1               UC Denver ERX1156148     AGP ICU microbiome SAMEA3611562
    ## 2               UC Denver ERX1156149     AGP ICU microbiome SAMEA3611563
    ## 3               UC Denver ERX1156150     AGP ICU microbiome SAMEA3611564
    ## 4               UC Denver ERX1156151     AGP ICU microbiome SAMEA3611565
    ## 5               UC Denver ERX1156152     AGP ICU microbiome SAMEA3611566
    ## 6               UC Denver ERX1156153     AGP ICU microbiome SAMEA3611567
    ##   extraction_robot..exp. extractionkit_lot..exp. has_physical_specimen
    ## 1                   HOWE                 PS13I16                  True
    ## 2                   HOWE                 PS13I16                  True
    ## 3                   NEWE                 PS13G10                  True
    ## 4                   HOWE                 PS13I16                  True
    ## 5                   HOWE                 PS13I16                  True
    ## 6                   HOWE                 PS13I16                  True
    ##   host_common_name host_scientific_name host_subject_id host_taxid
    ## 1            human         Homo sapiens           980-1       9606
    ## 2            human         Homo sapiens           980-1       9606
    ## 3            human         Homo sapiens           980-1       9606
    ## 4            human         Homo sapiens           980-1       9606
    ## 5            human         Homo sapiens           980-1       9606
    ## 6            human         Homo sapiens           980-2       9606
    ##   INSDC_center_alias                                        INSDC_center_name
    ## 1             UCSDMI University of California San Diego Microbiome Initiative
    ## 2             UCSDMI University of California San Diego Microbiome Initiative
    ## 3             UCSDMI University of California San Diego Microbiome Initiative
    ## 4             UCSDMI University of California San Diego Microbiome Initiative
    ## 5             UCSDMI University of California San Diego Microbiome Initiative
    ## 6             UCSDMI University of California San Diego Microbiome Initiative
    ##     INSDC_first_public    INSDC_last_update INSDC_status     Instrument
    ## 1 2015-10-15T17:01:22Z 2016-05-11T15:47:54Z       public Illumina MiSeq
    ## 2 2015-10-15T17:01:24Z 2016-05-11T15:47:54Z       public Illumina MiSeq
    ## 3 2015-10-15T17:01:25Z 2016-05-11T15:47:54Z       public Illumina MiSeq
    ## 4 2015-10-15T17:01:25Z 2016-05-11T15:47:54Z       public Illumina MiSeq
    ## 5 2015-10-15T17:01:19Z 2016-05-11T15:47:54Z       public Illumina MiSeq
    ## 6 2015-10-15T17:01:25Z 2016-05-11T15:47:54Z       public Illumina MiSeq
    ##   latitude   Library.Name LibraryLayout LibrarySelection LibrarySource
    ## 1     32.8 2136.000006143        SINGLE              PCR   METAGENOMIC
    ## 2     32.8 2136.000006144        SINGLE              PCR   METAGENOMIC
    ## 3     32.8 2136.000006145        SINGLE              PCR   METAGENOMIC
    ## 4     32.8 2136.000006146        SINGLE              PCR   METAGENOMIC
    ## 5     32.8 2136.000006147        SINGLE              PCR   METAGENOMIC
    ## 6     32.8 2136.000006149        SINGLE              PCR   METAGENOMIC
    ##   linker..exp. longitude mastermix_lot..exp.              Organism
    ## 1           GT    -117.3            14503598 human oral metagenome
    ## 2           GT    -117.3            14503598 human skin metagenome
    ## 3           GT    -117.3            14503598  human gut metagenome
    ## 4           GT    -117.3            14503598 human oral metagenome
    ## 5           GT    -117.3            14503598 human skin metagenome
    ## 6           GT    -117.3            14503598 human oral metagenome
    ##   orig_name..exp.                                 pcr_primers..exp.
    ## 1            6143 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 2            6144 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 3            6145 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 4            6146 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 5            6147 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 6            6149 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ##   physical_specimen_location physical_specimen_remaining Platform plating..exp.
    ## 1                     UCSDMI                        True ILLUMINA            GG
    ## 2                     UCSDMI                        True ILLUMINA            JG
    ## 3                     UCSDMI                        True ILLUMINA            GG
    ## 4                     UCSDMI                        True ILLUMINA            GG
    ## 5                     UCSDMI                        True ILLUMINA            JG
    ## 6                     UCSDMI                        True ILLUMINA            GG
    ##   primer_date..exp.        primer..exp. primer_plate..exp.
    ## 1            101013 GTGCCAGCMGCCGCGGTAA                  2
    ## 2            101013 GTGCCAGCMGCCGCGGTAA                  1
    ## 3            101013 GTGCCAGCMGCCGCGGTAA                  3
    ## 4            101013 GTGCCAGCMGCCGCGGTAA                  2
    ## 5            101013 GTGCCAGCMGCCGCGGTAA                  1
    ## 6            101013 GTGCCAGCMGCCGCGGTAA                  2
    ##   processing_robot..exp. project_name..exp. public          ReleaseDate
    ## 1                   JERE        AmericanGut   True 2015-10-28T00:00:00Z
    ## 2                   JERE        AmericanGut   True 2015-10-28T00:00:00Z
    ## 3                   ROBE        AmericanGut   True 2015-10-28T00:00:00Z
    ## 4                   JERE        AmericanGut   True 2015-10-28T00:00:00Z
    ## 5                   JERE        AmericanGut   True 2015-10-28T00:00:00Z
    ## 6                   JERE        AmericanGut   True 2015-10-28T00:00:00Z
    ##   required_sample_info_status run_center..exp. run_date..exp.
    ## 1                   completed             CCME       12/18/13
    ## 2                   completed             CCME       12/18/13
    ## 3                   completed             CCME       12/18/13
    ## 4                   completed             CCME       12/18/13
    ## 5                   completed             CCME       12/18/13
    ## 6                   completed             CCME       12/18/13
    ##          run_prefix..exp. samp_size..exp.  Sample.Name
    ## 1 AGP6_121813_16S_NoIndex           0.1,g SAMEA3611562
    ## 2 AGP6_121813_16S_NoIndex           0.1,g SAMEA3611563
    ## 3 AGP6_121813_16S_NoIndex           0.1,g SAMEA3611564
    ## 4 AGP6_121813_16S_NoIndex           0.1,g SAMEA3611565
    ## 5 AGP6_121813_16S_NoIndex           0.1,g SAMEA3611566
    ## 6 AGP6_121813_16S_NoIndex           0.1,g SAMEA3611567
    ##                     Sample_name sample_plate..exp. sample_type
    ## 1 qiita_sid_2136:2136.000006143           ICUOral2      saliva
    ## 2 qiita_sid_2136:2136.000006144           ICUSkin1        skin
    ## 3 qiita_sid_2136:2136.000006145          ICUFecal3       stool
    ## 4 qiita_sid_2136:2136.000006146           ICUOral2      saliva
    ## 5 qiita_sid_2136:2136.000006147           ICUSkin1        skin
    ## 6 qiita_sid_2136:2136.000006149           ICUOral2      saliva
    ##     sequencing_meth..exp. SRA.Study                  Submitter_Id
    ## 1 Sequencing by synthesis ERP012810 qiita_sid_2136:2136.000006143
    ## 2 Sequencing by synthesis ERP012810 qiita_sid_2136:2136.000006144
    ## 3 Sequencing by synthesis ERP012810 qiita_sid_2136:2136.000006145
    ## 4 Sequencing by synthesis ERP012810 qiita_sid_2136:2136.000006146
    ## 5 Sequencing by synthesis ERP012810 qiita_sid_2136:2136.000006147
    ## 6 Sequencing by synthesis ERP012810 qiita_sid_2136:2136.000006149
    ##   target_gene..exp. target_subfragment..exp. tm1000_8_tool..exp.
    ## 1          16S rRNA                       V4             311169B
    ## 2          16S rRNA                       V4             311169B
    ## 3          16S rRNA                       V4             311172B
    ## 4          16S rRNA                       V4             311169B
    ## 5          16S rRNA                       V4             311169B
    ## 6          16S rRNA                       V4             311169B
    ##   tm300_8_tool..exp. tm50_8_tool..exp. water_lot..exp. well_description..exp.
    ## 1            111569B           311426B        RNBC2791       ICUOral2_6143_A1
    ## 2            111569B           311426B        RNBC2791       ICUSkin1_6144_A2
    ## 3            311313B           111075B        RNBC2791      ICUFecal3_6145_A1
    ## 4            111569B           311426B        RNBC2791       ICUOral2_6146_H1
    ## 5            111569B           311426B        RNBC2791       ICUSkin1_6147_G2
    ## 6            111569B           311426B        RNBC2791       ICUOral2_6149_F3
    ##   well_id..exp.          create_date version adequacycalsen_class
    ## 1            A1 2015-12-08T16:15:00Z       1                  mid
    ## 2            A2 2015-12-08T16:15:00Z       1                  mid
    ## 3            A1 2015-12-08T16:15:00Z       1                  mid
    ## 4            H1 2015-12-08T16:15:00Z       1                  mid
    ## 5            G2 2015-12-08T16:15:00Z       1                  mid
    ## 6            F3 2015-12-08T16:15:00Z       1                  low
    ##   adequacycalsen_median adequacyproten_class adequacyproten_median
    ## 1              40.77586                  low              27.48175
    ## 2              40.77586                  low              27.48175
    ## 3              40.77586                  low              27.48175
    ## 4              40.77586                  low              27.48175
    ## 5              40.77586                  low              27.48175
    ## 6               0.00000                  low               0.00000
    ##      admission_type age_category AGE apache_score ards_yn        bmi_category
    ## 1           Medical          60s  67           31      No Very severely obese
    ## 2           Medical          60s  67           31      No Very severely obese
    ## 3           Medical          60s  67           31      No Very severely obese
    ## 4           Medical          60s  67           31      No Very severely obese
    ## 5           Medical          60s  67           31      No Very severely obese
    ## 6 Surgical Elective          20s  29           21      No Very severely obese
    ##    BMI       body_habitat  body_product     body_site daysfromhosp
    ## 1 42.1 UBERON:oral cavity UBERON:saliva UBERON:tongue   0.09444445
    ## 2 42.1        UBERON:skin  UBERON:sebum   UBERON:skin   0.09444445
    ## 3 42.1       UBERON:feces  UBERON:feces  UBERON:feces   0.09444445
    ## 4 42.1 UBERON:oral cavity UBERON:saliva UBERON:tongue   0.09444445
    ## 5 42.1        UBERON:skin  UBERON:sebum   UBERON:skin   0.09444445
    ## 6 40.8 UBERON:oral cavity UBERON:saliva UBERON:tongue   0.29513890
    ##                                       first_pn highest_glucose_mmoll
    ## 1 PN not initiated during first 12 days in ICU             10.555555
    ## 2 PN not initiated during first 12 days in ICU             10.555555
    ## 3 PN not initiated during first 12 days in ICU             10.555555
    ## 4 PN not initiated during first 12 days in ICU             10.555555
    ## 5 PN not initiated during first 12 days in ICU             10.555555
    ## 6 PN not initiated during first 12 days in ICU              8.666667
    ##   highest_glucose hmp_site hospdays_lt_22  hospdays hospevent
    ## 1             190     ORAL              f 25.490278         1
    ## 2             190     SKIN              f 25.490278         1
    ## 3             190    FECAL              f 25.490278         1
    ## 4             190     ORAL              f 25.490278         1
    ## 5             190     SKIN              f 25.490278         1
    ## 6             156     ORAL              t  6.104166         1
    ##                                           icu icudays_lt_15   icudays icuevent
    ## 1 University of Colorado - Surgical - medical             f 25.490278        1
    ## 2 University of Colorado - Surgical - medical             f 25.490278        1
    ## 3 University of Colorado - Surgical - medical             f 25.490278        1
    ## 4 University of Colorado - Surgical - medical             f 25.490278        1
    ## 5 University of Colorado - Surgical - medical             f 25.490278        1
    ## 6 University of Colorado - Surgical - medical             t  4.177084        1
    ##   is_pgp lowest_glucose_mmoll lowest_glucose mv_start_yn nutric_age
    ## 1  FALSE             7.555555            136         Yes          1
    ## 2  FALSE             7.555555            136         Yes          1
    ## 3  FALSE             7.555555            136         Yes          1
    ## 4  FALSE             7.555555            136         Yes          1
    ## 5  FALSE             7.555555            136         Yes          1
    ## 6  FALSE             6.722222            121         Yes          0
    ##   nutric_apache nutric_comob nutric_score_lt_5 nutric_score nutric_sofa
    ## 1             3            1                 f            6           1
    ## 2             3            1                 f            6           1
    ## 3             3            1                 f            6           1
    ## 4             3            1                 f            6           1
    ## 5             3            1                 f            6           1
    ## 6             2            0                 t            3           1
    ##   patientdied        pn_gt_5 probiotics_yn pt_die_icu_yn    sex
    ## 1           0 NO_PN_RECEIVED            No            No female
    ## 2           0 NO_PN_RECEIVED            No            No female
    ## 3           0 NO_PN_RECEIVED            No            No female
    ## 4           0 NO_PN_RECEIVED            No            No female
    ## 5           0 NO_PN_RECEIVED            No            No female
    ## 6           0 NO_PN_RECEIVED            No            No female
    ##   simple_body_site sofa_score study_id  timepoint title_acronym title_body_site
    ## 1             ORAL          9    980-1  ADMISSION           ICU        ICU-ORAL
    ## 2             SKIN          9    980-1  ADMISSION           ICU        ICU-SKIN
    ## 3            FECAL          9    980-1 D10 or D/C           ICU       ICU-FECAL
    ## 4             ORAL          9    980-1 D10 or D/C           ICU        ICU-ORAL
    ## 5             SKIN          9    980-1 D10 or D/C           ICU        ICU-SKIN
    ## 6             ORAL          7    980-2  ADMISSION           ICU        ICU-ORAL
    ##   ventdays_lt_7  ventdays daily_kcal_presc_entered icu_day_60_yn pt_die_hosp_yn
    ## 1             f 25.490278                     2320            No             No
    ## 2             f 25.490278                     2320            No             No
    ## 3             f 25.490278                     2320            No             No
    ## 4             f 25.490278                     2320            No             No
    ## 5             f 25.490278                     2320            No             No
    ## 6             t  3.919445                     1899            No             No

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
  Age = initial_meta$AGE,
  BMI = initial_meta$BMI,
  Sex=initial_meta$sex,
  Region = initial_meta$target_subfragment..exp.,
  Platform = initial_meta$Instrument,
  Country = initial_meta$COUNTRY,
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
metadata$PCR_primers <- NA_character_
metadata$title <- "Extreme Dysbiosis of the Microbiome in Critical Illness"
metadata$Disease_study <- "YES"

# Assign disease status based on sample_name patterns
metadata$Status <- "ICU patients"
keep_subjects<-subset(metadata,Site %in% c("saliva"))
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

    ## [1] "Number of unique runs: 435"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 92"

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

    ##   BioProject        Run sample_name Age  BMI    Sex Region       Platform
    ## 1 PRJEB11425 ERR1076083       980-1  67 42.1 female     V4 Illumina MiSeq
    ## 2 PRJEB11425 ERR1076084       980-1  67 42.1 female     V4 Illumina MiSeq
    ## 3 PRJEB11425 ERR1076085       980-1  67 42.1 female     V4 Illumina MiSeq
    ## 4 PRJEB11425 ERR1076086       980-1  67 42.1 female     V4 Illumina MiSeq
    ## 5 PRJEB11425 ERR1076087       980-1  67 42.1 female     V4 Illumina MiSeq
    ## 6 PRJEB11425 ERR1076088       980-2  29 40.8 female     V4 Illumina MiSeq
    ##   Country   Site          public_time PCR_primers
    ## 1     USA saliva 2015-10-28T00:00:00Z        <NA>
    ## 2     USA  sebum 2015-10-28T00:00:00Z        <NA>
    ## 3     USA  feces 2015-10-28T00:00:00Z        <NA>
    ## 4     USA saliva 2015-10-28T00:00:00Z        <NA>
    ## 5     USA  sebum 2015-10-28T00:00:00Z        <NA>
    ## 6     USA saliva 2015-10-28T00:00:00Z        <NA>
    ##                                                     title Disease_study
    ## 1 Extreme Dysbiosis of the Microbiome in Critical Illness           YES
    ## 2 Extreme Dysbiosis of the Microbiome in Critical Illness           YES
    ## 3 Extreme Dysbiosis of the Microbiome in Critical Illness           YES
    ## 4 Extreme Dysbiosis of the Microbiome in Critical Illness           YES
    ## 5 Extreme Dysbiosis of the Microbiome in Critical Illness           YES
    ## 6 Extreme Dysbiosis of the Microbiome in Critical Illness           YES
    ##         Status matched_samples
    ## 1 ICU patients             YES
    ## 2 ICU patients             YES
    ## 3 ICU patients             YES
    ## 4 ICU patients             YES
    ## 5 ICU patients             YES
    ## 6 ICU patients             YES

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
