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
setwd("D:/oral_project/PRJEB11425")

# Read the initial metadata files
initial_meta<-read.csv("SraRunTable.csv",header = T)

# Combine the two metadata files
#meta <- rbind(initial_meta, initial_meta1)

#Remove entries where Organism is "gut metagenome" (non-oral samples)
initial_meta <- subset(initial_meta,env_material== "saliva")

# Check the number of samples containing "HT" in Sample.Name (e.g., healthy controls)
#nrow(subset(meta, grepl("HT", Sample.Name)))

# Display the first few rows of the metadata
head(initial_meta)
```

    ##           Run age_unit age_units anonymized_name Assay.Type AvgSpotLen
    ## 1  ERR1076083    years     years            6143      OTHER        151
    ## 4  ERR1076086    years     years            6146      OTHER        151
    ## 6  ERR1076088    years     years            6149      OTHER        151
    ## 8  ERR1076090    years     years            6152      OTHER        151
    ## 10 ERR1076092    years     years            6155      OTHER        151
    ## 12 ERR1076094    years     years            6158      OTHER        151
    ##    barcode..exp.   Bases BioProject    BioSample   Bytes center_name..exp.
    ## 1   TACCGCTTCTTC 1065909 PRJEB11425 SAMEA3611562  838881              CCME
    ## 4   GTGCAACCAATC 2813432 PRJEB11425 SAMEA3611565 1967312              CCME
    ## 6   GATCACGAGAGG 4977564 PRJEB11425 SAMEA3611567 3195861              CCME
    ## 8   TGCAGCAAGATT 5256763 PRJEB11425 SAMEA3611569 2839246              CCME
    ## 10  TTAGAGCCATGC 4285229 PRJEB11425 SAMEA3611571 2487414              CCME
    ## 12  ACAACACTCCGA 5036454 PRJEB11425 SAMEA3611573 3047079              CCME
    ##                                           Center.Name center_project_name..exp.
    ## 1  UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT            ICU microbiome
    ## 4  UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT            ICU microbiome
    ## 6  UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT            ICU microbiome
    ## 8  UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT            ICU microbiome
    ## 10 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT            ICU microbiome
    ## 12 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT            ICU microbiome
    ##    collection_timestamp           common_name Consent COUNTRY
    ## 1                1/1/59 human oral metagenome  public     USA
    ## 4                1/1/59 human oral metagenome  public     USA
    ## 6                1/1/59 human oral metagenome  public     USA
    ## 8                1/1/59 human oral metagenome  public     USA
    ## 10               1/1/59 human oral metagenome  public     USA
    ## 12               1/1/59 human oral metagenome  public     USA
    ##    DATASTORE.filetype DATASTORE.provider
    ## 1    fastq,run.zq,sra     ena,gs,ncbi,s3
    ## 4    fastq,run.zq,sra     ena,gs,ncbi,s3
    ## 6    fastq,run.zq,sra     ena,gs,ncbi,s3
    ## 8    fastq,run.zq,sra     ena,gs,ncbi,s3
    ## 10   fastq,run.zq,sra     ena,gs,ncbi,s3
    ## 12   fastq,run.zq,sra     ena,gs,ncbi,s3
    ##                            DATASTORE.region Depth dna_extracted elevation
    ## 1  ena,gs.us-east1,ncbi.public,s3.us-east-1     0          True       193
    ## 4  ena,gs.us-east1,ncbi.public,s3.us-east-1     0          True       193
    ## 6  ena,gs.us-east1,ncbi.public,s3.us-east-1     0          True       193
    ## 8  ena,gs.us-east1,ncbi.public,s3.us-east-1     0          True       193
    ## 10 ena,gs.us-east1,ncbi.public,s3.us-east-1     0          True       193
    ## 12 ena,gs.us-east1,ncbi.public,s3.us-east-1     0          True       193
    ##    emp_status..exp. ENA.FIRST.PUBLIC..run. ENA_first_public
    ## 1               EMP             2015-10-15       2015-10-15
    ## 4               EMP             2015-10-15       2015-10-15
    ## 6               EMP             2015-10-15       2015-10-15
    ## 8               EMP             2015-10-15       2015-10-15
    ## 10              EMP             2015-10-15       2015-10-15
    ## 12              EMP             2015-10-15       2015-10-15
    ##    ENA.LAST.UPDATE..run. ENA_last_update   env_biome              env_feature
    ## 1             2018-11-16      2016-05-11 urban biome human-associated habitat
    ## 4             2018-11-16      2016-05-11 urban biome human-associated habitat
    ## 6             2018-11-16      2016-05-11 urban biome human-associated habitat
    ## 8             2018-11-16      2016-05-11 urban biome human-associated habitat
    ## 10            2018-11-16      2016-05-11 urban biome human-associated habitat
    ## 12            2018-11-16      2016-05-11 urban biome human-associated habitat
    ##    env_material  env_matter Env_package experiment_center..exp. Experiment
    ## 1        saliva ENVO:saliva  human-oral               UC Denver ERX1156148
    ## 4        saliva ENVO:saliva  human-oral               UC Denver ERX1156151
    ## 6        saliva ENVO:saliva  human-oral               UC Denver ERX1156153
    ## 8        saliva ENVO:saliva  human-oral               UC Denver ERX1156155
    ## 10       saliva ENVO:saliva  human-oral               UC Denver ERX1156157
    ## 12       saliva ENVO:saliva  human-oral               UC Denver ERX1156159
    ##    experiment_title..exp.  external_id extraction_robot..exp.
    ## 1      AGP ICU microbiome SAMEA3611562                   HOWE
    ## 4      AGP ICU microbiome SAMEA3611565                   HOWE
    ## 6      AGP ICU microbiome SAMEA3611567                   HOWE
    ## 8      AGP ICU microbiome SAMEA3611569                   HOWE
    ## 10     AGP ICU microbiome SAMEA3611571                   HOWE
    ## 12     AGP ICU microbiome SAMEA3611573                   HOWE
    ##    extractionkit_lot..exp. has_physical_specimen host_common_name
    ## 1                  PS13I16                  True            human
    ## 4                  PS13I16                  True            human
    ## 6                  PS13I16                  True            human
    ## 8                  PS13I16                  True            human
    ## 10                 PS13I16                  True            human
    ## 12                 PS13I16                  True            human
    ##    host_scientific_name host_subject_id host_taxid INSDC_center_alias
    ## 1          Homo sapiens           980-1       9606             UCSDMI
    ## 4          Homo sapiens           980-1       9606             UCSDMI
    ## 6          Homo sapiens           980-2       9606             UCSDMI
    ## 8          Homo sapiens           980-2       9606             UCSDMI
    ## 10         Homo sapiens           980-3       9606             UCSDMI
    ## 12         Homo sapiens           980-3       9606             UCSDMI
    ##                                           INSDC_center_name
    ## 1  University of California San Diego Microbiome Initiative
    ## 4  University of California San Diego Microbiome Initiative
    ## 6  University of California San Diego Microbiome Initiative
    ## 8  University of California San Diego Microbiome Initiative
    ## 10 University of California San Diego Microbiome Initiative
    ## 12 University of California San Diego Microbiome Initiative
    ##      INSDC_first_public    INSDC_last_update INSDC_status     Instrument
    ## 1  2015-10-15T17:01:22Z 2016-05-11T15:47:54Z       public Illumina MiSeq
    ## 4  2015-10-15T17:01:25Z 2016-05-11T15:47:54Z       public Illumina MiSeq
    ## 6  2015-10-15T17:01:25Z 2016-05-11T15:47:54Z       public Illumina MiSeq
    ## 8  2015-10-15T17:01:22Z 2016-05-11T15:47:54Z       public Illumina MiSeq
    ## 10 2015-10-15T17:01:22Z 2016-05-11T15:47:54Z       public Illumina MiSeq
    ## 12 2015-10-15T17:01:19Z 2016-05-11T15:47:54Z       public Illumina MiSeq
    ##    latitude   Library.Name LibraryLayout LibrarySelection LibrarySource
    ## 1      32.8 2136.000006143        SINGLE              PCR   METAGENOMIC
    ## 4      32.8 2136.000006146        SINGLE              PCR   METAGENOMIC
    ## 6      32.8 2136.000006149        SINGLE              PCR   METAGENOMIC
    ## 8      32.8 2136.000006152        SINGLE              PCR   METAGENOMIC
    ## 10     32.8 2136.000006155        SINGLE              PCR   METAGENOMIC
    ## 12     32.8 2136.000006158        SINGLE              PCR   METAGENOMIC
    ##    linker..exp. longitude mastermix_lot..exp.              Organism
    ## 1            GT    -117.3            14503598 human oral metagenome
    ## 4            GT    -117.3            14503598 human oral metagenome
    ## 6            GT    -117.3            14503598 human oral metagenome
    ## 8            GT    -117.3            14503598 human oral metagenome
    ## 10           GT    -117.3            14503598 human oral metagenome
    ## 12           GT    -117.3            14503598 human oral metagenome
    ##    orig_name..exp.                                 pcr_primers..exp.
    ## 1             6143 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 4             6146 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 6             6149 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 8             6152 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 10            6155 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ## 12            6158 FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT
    ##    physical_specimen_location physical_specimen_remaining Platform
    ## 1                      UCSDMI                        True ILLUMINA
    ## 4                      UCSDMI                        True ILLUMINA
    ## 6                      UCSDMI                        True ILLUMINA
    ## 8                      UCSDMI                        True ILLUMINA
    ## 10                     UCSDMI                        True ILLUMINA
    ## 12                     UCSDMI                        True ILLUMINA
    ##    plating..exp. primer_date..exp.        primer..exp. primer_plate..exp.
    ## 1             GG            101013 GTGCCAGCMGCCGCGGTAA                  2
    ## 4             GG            101013 GTGCCAGCMGCCGCGGTAA                  2
    ## 6             GG            101013 GTGCCAGCMGCCGCGGTAA                  2
    ## 8             GG            101013 GTGCCAGCMGCCGCGGTAA                  2
    ## 10            GG            101013 GTGCCAGCMGCCGCGGTAA                  2
    ## 12            GG            101013 GTGCCAGCMGCCGCGGTAA                  2
    ##    processing_robot..exp. project_name..exp. public          ReleaseDate
    ## 1                    JERE        AmericanGut   True 2015-10-28T00:00:00Z
    ## 4                    JERE        AmericanGut   True 2015-10-28T00:00:00Z
    ## 6                    JERE        AmericanGut   True 2015-10-28T00:00:00Z
    ## 8                    JERE        AmericanGut   True 2015-10-28T00:00:00Z
    ## 10                   JERE        AmericanGut   True 2015-10-28T00:00:00Z
    ## 12                   JERE        AmericanGut   True 2015-10-28T00:00:00Z
    ##    required_sample_info_status run_center..exp. run_date..exp.
    ## 1                    completed             CCME       12/18/13
    ## 4                    completed             CCME       12/18/13
    ## 6                    completed             CCME       12/18/13
    ## 8                    completed             CCME       12/18/13
    ## 10                   completed             CCME       12/18/13
    ## 12                   completed             CCME       12/18/13
    ##           run_prefix..exp. samp_size..exp.  Sample.Name
    ## 1  AGP6_121813_16S_NoIndex           0.1,g SAMEA3611562
    ## 4  AGP6_121813_16S_NoIndex           0.1,g SAMEA3611565
    ## 6  AGP6_121813_16S_NoIndex           0.1,g SAMEA3611567
    ## 8  AGP6_121813_16S_NoIndex           0.1,g SAMEA3611569
    ## 10 AGP6_121813_16S_NoIndex           0.1,g SAMEA3611571
    ## 12 AGP6_121813_16S_NoIndex           0.1,g SAMEA3611573
    ##                      Sample_name sample_plate..exp. sample_type
    ## 1  qiita_sid_2136:2136.000006143           ICUOral2      saliva
    ## 4  qiita_sid_2136:2136.000006146           ICUOral2      saliva
    ## 6  qiita_sid_2136:2136.000006149           ICUOral2      saliva
    ## 8  qiita_sid_2136:2136.000006152           ICUOral2      saliva
    ## 10 qiita_sid_2136:2136.000006155           ICUOral2      saliva
    ## 12 qiita_sid_2136:2136.000006158           ICUOral2      saliva
    ##      sequencing_meth..exp. SRA.Study                  Submitter_Id
    ## 1  Sequencing by synthesis ERP012810 qiita_sid_2136:2136.000006143
    ## 4  Sequencing by synthesis ERP012810 qiita_sid_2136:2136.000006146
    ## 6  Sequencing by synthesis ERP012810 qiita_sid_2136:2136.000006149
    ## 8  Sequencing by synthesis ERP012810 qiita_sid_2136:2136.000006152
    ## 10 Sequencing by synthesis ERP012810 qiita_sid_2136:2136.000006155
    ## 12 Sequencing by synthesis ERP012810 qiita_sid_2136:2136.000006158
    ##    target_gene..exp. target_subfragment..exp. tm1000_8_tool..exp.
    ## 1           16S rRNA                       V4             311169B
    ## 4           16S rRNA                       V4             311169B
    ## 6           16S rRNA                       V4             311169B
    ## 8           16S rRNA                       V4             311169B
    ## 10          16S rRNA                       V4             311169B
    ## 12          16S rRNA                       V4             311169B
    ##    tm300_8_tool..exp. tm50_8_tool..exp. water_lot..exp. well_description..exp.
    ## 1             111569B           311426B        RNBC2791       ICUOral2_6143_A1
    ## 4             111569B           311426B        RNBC2791       ICUOral2_6146_H1
    ## 6             111569B           311426B        RNBC2791       ICUOral2_6149_F3
    ## 8             111569B           311426B        RNBC2791       ICUOral2_6152_D3
    ## 10            111569B           311426B        RNBC2791       ICUOral2_6155_E1
    ## 12            111569B           311426B        RNBC2791       ICUOral2_6158_E4
    ##    well_id..exp.          create_date version adequacycalsen_class
    ## 1             A1 2015-12-08T16:15:00Z       1                  mid
    ## 4             H1 2015-12-08T16:15:00Z       1                  mid
    ## 6             F3 2015-12-08T16:15:00Z       1                  low
    ## 8             D3 2015-12-08T16:15:00Z       1                  low
    ## 10            E1 2015-12-08T16:15:00Z       1               normal
    ## 12            E4 2015-12-08T16:15:00Z       1               normal
    ##    adequacycalsen_median adequacyproten_class adequacyproten_median
    ## 1               40.77586                  low              27.48175
    ## 4               40.77586                  low              27.48175
    ## 6                0.00000                  low               0.00000
    ## 8                0.00000                  low               0.00000
    ## 10              82.35904                 high             110.99290
    ## 12              82.35904                 high             110.99290
    ##       admission_type age_category AGE apache_score ards_yn        bmi_category
    ## 1            Medical          60s  67           31      No Very severely obese
    ## 4            Medical          60s  67           31      No Very severely obese
    ## 6  Surgical Elective          20s  29           21      No Very severely obese
    ## 8  Surgical Elective          20s  29           21      No Very severely obese
    ## 10           Medical          40s  47           22      No          Overweight
    ## 12           Medical          40s  47           22      No          Overweight
    ##     BMI       body_habitat  body_product     body_site daysfromhosp
    ## 1  42.1 UBERON:oral cavity UBERON:saliva UBERON:tongue   0.09444445
    ## 4  42.1 UBERON:oral cavity UBERON:saliva UBERON:tongue   0.09444445
    ## 6  40.8 UBERON:oral cavity UBERON:saliva UBERON:tongue   0.29513890
    ## 8  40.8 UBERON:oral cavity UBERON:saliva UBERON:tongue   0.29513890
    ## 10 28.4 UBERON:oral cavity UBERON:saliva UBERON:tongue   0.00000000
    ## 12 28.4 UBERON:oral cavity UBERON:saliva UBERON:tongue   0.00000000
    ##                                        first_pn highest_glucose_mmoll
    ## 1  PN not initiated during first 12 days in ICU             10.555555
    ## 4  PN not initiated during first 12 days in ICU             10.555555
    ## 6  PN not initiated during first 12 days in ICU              8.666667
    ## 8  PN not initiated during first 12 days in ICU              8.666667
    ## 10 PN not initiated during first 12 days in ICU             10.222222
    ## 12 PN not initiated during first 12 days in ICU             10.222222
    ##    highest_glucose hmp_site hospdays_lt_22  hospdays hospevent
    ## 1              190     ORAL              f 25.490278         1
    ## 4              190     ORAL              f 25.490278         1
    ## 6              156     ORAL              t  6.104166         1
    ## 8              156     ORAL              t  6.104166         1
    ## 10             184     ORAL              t 13.945833         1
    ## 12             184     ORAL              t 13.945833         1
    ##                                            icu icudays_lt_15   icudays icuevent
    ## 1  University of Colorado - Surgical - medical             f 25.490278        1
    ## 4  University of Colorado - Surgical - medical             f 25.490278        1
    ## 6  University of Colorado - Surgical - medical             t  4.177084        1
    ## 8  University of Colorado - Surgical - medical             t  4.177084        1
    ## 10 University of Colorado - Surgical - medical             t 13.945833        1
    ## 12 University of Colorado - Surgical - medical             t 13.945833        1
    ##    is_pgp lowest_glucose_mmoll lowest_glucose mv_start_yn nutric_age
    ## 1   FALSE             7.555555            136         Yes          1
    ## 4   FALSE             7.555555            136         Yes          1
    ## 6   FALSE             6.722222            121         Yes          0
    ## 8   FALSE             6.722222            121         Yes          0
    ## 10  FALSE             4.944445             89         Yes          0
    ## 12  FALSE             4.944445             89         Yes          0
    ##    nutric_apache nutric_comob nutric_score_lt_5 nutric_score nutric_sofa
    ## 1              3            1                 f            6           1
    ## 4              3            1                 f            6           1
    ## 6              2            0                 t            3           1
    ## 8              2            0                 t            3           1
    ## 10             2            0                 t            3           1
    ## 12             2            0                 t            3           1
    ##    patientdied        pn_gt_5 probiotics_yn pt_die_icu_yn    sex
    ## 1            0 NO_PN_RECEIVED            No            No female
    ## 4            0 NO_PN_RECEIVED            No            No female
    ## 6            0 NO_PN_RECEIVED            No            No female
    ## 8            0 NO_PN_RECEIVED            No            No female
    ## 10           0 NO_PN_RECEIVED            No            No   male
    ## 12           0 NO_PN_RECEIVED            No            No   male
    ##    simple_body_site sofa_score study_id  timepoint title_acronym
    ## 1              ORAL          9    980-1  ADMISSION           ICU
    ## 4              ORAL          9    980-1 D10 or D/C           ICU
    ## 6              ORAL          7    980-2  ADMISSION           ICU
    ## 8              ORAL          7    980-2 D10 or D/C           ICU
    ## 10             ORAL          9    980-3  ADMISSION           ICU
    ## 12             ORAL          9    980-3 D10 or D/C           ICU
    ##    title_body_site ventdays_lt_7   ventdays daily_kcal_presc_entered
    ## 1         ICU-ORAL             f 25.4902780                     2320
    ## 4         ICU-ORAL             f 25.4902780                     2320
    ## 6         ICU-ORAL             t  3.9194446                     1899
    ## 8         ICU-ORAL             t  3.9194446                     1899
    ## 10        ICU-ORAL             t  0.8895833                     2412
    ## 12        ICU-ORAL             t  0.8895833                     2412
    ##    icu_day_60_yn pt_die_hosp_yn
    ## 1             No             No
    ## 4             No             No
    ## 6             No             No
    ## 8             No             No
    ## 10            No             No
    ## 12            No             No

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
  sample_name = initial_meta$Sample_name,
  Age = initial_meta$AGE,
  BMI = initial_meta$BMI,
  Sex=initial_meta$sex,
  Region = initial_meta$target_subfragment..exp.,
  Platform = initial_meta$Instrument,
  Country = initial_meta$COUNTRY,
  Site = initial_meta$env_material,
  public_time = initial_meta$ENA.FIRST.PUBLIC..run.,
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
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 171"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 171"

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
