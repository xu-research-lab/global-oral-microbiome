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
setwd("D:/oral_project/PRJEB1799")

# Read the initial metadata files
initial_meta<-read.csv("SraRunTable.csv",header = T)

# Combine the two metadata files
#meta <- rbind(initial_meta, initial_meta1)

#Remove entries non-oral samples
#initial_meta <- subset(initial_meta,isolate %in% c("saliva","tongue coat"))
initial_meta <- initial_meta[initial_meta$Organism != "gut metagenome", ]
initial_meta <- initial_meta[initial_meta$Organism != "human gut metagenome", ]
initial_meta <- initial_meta[initial_meta$Organism != "organismal metagenomes", ]
initial_meta <- initial_meta[initial_meta$Organism != "Parnassius szechenyii incognitus", ]
# Check the number of samples containing "HT" in Sample.Name (e.g., healthy controls)
#nrow(subset(meta, grepl("HT", Sample.Name)))

# Display the first few rows of the metadata
head(initial_meta)
```

    ##          Run age_unit altitude anonymized_name Assay.Type assigned_from_geo
    ## 3  ERR261871       78        0          FS_174   AMPLICON                 n
    ## 7  ERR261875       78        0          FS_153   AMPLICON                 n
    ## 12 ERR261880       78        0           FS_10   AMPLICON                 n
    ## 15 ERR261883       78        0           FS_28   AMPLICON                 n
    ## 19 ERR261887       78        0          FS_141   AMPLICON                 n
    ## 25 ERR261893       78        0          FS_105   AMPLICON                 n
    ##    AvgSpotLen   Bases BioProject    BioSample       body_habitat  body_product
    ## 3         125  788791  PRJEB1799 SAMEA1881262 UBERON:oral cavity UBERON:saliva
    ## 7         126 2142990  PRJEB1799 SAMEA1878909 UBERON:oral cavity UBERON:saliva
    ## 12        126 6794075  PRJEB1799 SAMEA1879886 UBERON:oral cavity UBERON:saliva
    ## 15        125 4467465  PRJEB1799 SAMEA1880711 UBERON:oral cavity UBERON:saliva
    ## 19        128 4632400  PRJEB1799 SAMEA1880524 UBERON:oral cavity UBERON:saliva
    ## 25        125 3258282  PRJEB1799 SAMEA1878915 UBERON:oral cavity UBERON:saliva
    ##        body_site   Bytes   Center.Name           common_name Consent
    ## 3  UBERON:tongue  642184 CCME-COLORADO human oral metagenome  public
    ## 7  UBERON:tongue 1651774 CCME-COLORADO human oral metagenome  public
    ## 12 UBERON:tongue 5217211 CCME-COLORADO human oral metagenome  public
    ## 15 UBERON:tongue 3497919 CCME-COLORADO human oral metagenome  public
    ## 19 UBERON:tongue 3628760 CCME-COLORADO human oral metagenome  public
    ## 25 UBERON:tongue 2542006 CCME-COLORADO human oral metagenome  public
    ##    DATASTORE.filetype DATASTORE.provider
    ## 3    fastq,run.zq,sra     ena,gs,ncbi,s3
    ## 7    fastq,run.zq,sra     ena,gs,ncbi,s3
    ## 12   fastq,run.zq,sra     ena,gs,ncbi,s3
    ## 15   fastq,run.zq,sra     ena,gs,ncbi,s3
    ## 19   fastq,run.zq,sra     ena,gs,ncbi,s3
    ## 25   fastq,run.zq,sra     ena,gs,ncbi,s3
    ##                            DATASTORE.region Depth elevation
    ## 3  ena,gs.us-east1,ncbi.public,s3.us-east-1     0    1624.1
    ## 7  ena,gs.us-east1,ncbi.public,s3.us-east-1     0    1624.1
    ## 12 ena,gs.us-east1,ncbi.public,s3.us-east-1     0    1624.1
    ## 15 ena,gs.us-east1,ncbi.public,s3.us-east-1     0    1624.1
    ## 19 ena,gs.us-east1,ncbi.public,s3.us-east-1     0    1624.1
    ## 25 ena,gs.us-east1,ncbi.public,s3.us-east-1     0    1624.1
    ##    ENA.FIRST.PUBLIC..run. ena_first_public ENA.LAST.UPDATE..run.
    ## 3              2013-04-17       2013-04-17            2018-11-16
    ## 7              2013-04-17       2013-04-17            2018-11-16
    ## 12             2013-04-17       2013-04-17            2018-11-16
    ## 15             2013-04-17       2013-04-17            2018-11-16
    ## 19             2013-04-17       2013-04-17            2018-11-16
    ## 25             2013-04-17       2013-04-17            2018-11-16
    ##    ena_last_update                     env_biome                   env_feature
    ## 3       2013-04-16 ENVO:human-associated habitat ENVO:human-associated habitat
    ## 7       2013-04-16 ENVO:human-associated habitat ENVO:human-associated habitat
    ## 12      2013-04-16 ENVO:human-associated habitat ENVO:human-associated habitat
    ## 15      2013-04-16 ENVO:human-associated habitat ENVO:human-associated habitat
    ## 19      2013-04-16 ENVO:human-associated habitat ENVO:human-associated habitat
    ## 25      2013-04-16 ENVO:human-associated habitat ENVO:human-associated habitat
    ##     env_matter experiment_center..exp.
    ## 3  ENVO:saliva                    CCME
    ## 7  ENVO:saliva                    CCME
    ## 12 ENVO:saliva                    CCME
    ## 15 ENVO:saliva                    CCME
    ## 19 ENVO:saliva                    CCME
    ## 25 ENVO:saliva                    CCME
    ##                                               experiment_design_description..exp.
    ## 3  Explaining variability among and between the microbial communities of families
    ## 7  Explaining variability among and between the microbial communities of families
    ## 12 Explaining variability among and between the microbial communities of families
    ## 15 Explaining variability among and between the microbial communities of families
    ## 19 Explaining variability among and between the microbial communities of families
    ## 25 Explaining variability among and between the microbial communities of families
    ##    Experiment
    ## 3   ERX236389
    ## 7   ERX236393
    ## 12  ERX236398
    ## 15  ERX236401
    ## 19  ERX236405
    ## 25  ERX236411
    ##                                                            experiment_title..exp.
    ## 3  Explaining variability among and between the microbial communities of families
    ## 7  Explaining variability among and between the microbial communities of families
    ## 12 Explaining variability among and between the microbial communities of families
    ## 15 Explaining variability among and between the microbial communities of families
    ## 19 Explaining variability among and between the microbial communities of families
    ## 25 Explaining variability among and between the microbial communities of families
    ##     external_id family geo_loc_name_country geo_loc_name_country_continent
    ## 3  SAMEA1881262   1075                  USA                  North America
    ## 7  SAMEA1878909   1053                  USA                  North America
    ## 12 SAMEA1879886      5                  USA                  North America
    ## 15 SAMEA1880711     18                  USA                  North America
    ## 19 SAMEA1880524     68                  USA                  North America
    ## 25 SAMEA1878915     53                  USA                  North America
    ##    geographic_location_.country_and.or_sea. host_common_name host_subject_id
    ## 3                                  USA: GAZ            human      797:FS_174
    ## 7                                  USA: GAZ            human      797:FS_153
    ## 12                                 USA: GAZ            human       797:FS_10
    ## 15                                 USA: GAZ            human       797:FS_28
    ## 19                                 USA: GAZ            human      797:FS_141
    ## 25                                 USA: GAZ            human      797:FS_105
    ##    host_taxid
    ## 3        9606
    ## 7        9606
    ## 12       9606
    ## 15       9606
    ## 19       9606
    ## 25       9606
    ##                                                             INSDC_center_name
    ## 3  Colorado Center for Microbial Ecology\\, University of Colorado at Boulder
    ## 7  Colorado Center for Microbial Ecology\\, University of Colorado at Boulder
    ## 12 Colorado Center for Microbial Ecology\\, University of Colorado at Boulder
    ## 15 Colorado Center for Microbial Ecology\\, University of Colorado at Boulder
    ## 19 Colorado Center for Microbial Ecology\\, University of Colorado at Boulder
    ## 25 Colorado Center for Microbial Ecology\\, University of Colorado at Boulder
    ##      INSDC_first_public    INSDC_last_update INSDC_status  Instrument latitude
    ## 3  2013-04-17T17:01:49Z 2013-04-16T19:29:10Z       public unspecified    40.01
    ## 7  2013-04-17T17:01:50Z 2013-04-16T19:29:10Z       public unspecified    40.01
    ## 12 2013-04-17T17:01:50Z 2013-04-16T19:29:10Z       public unspecified    40.01
    ## 15 2013-04-17T17:01:50Z 2013-04-16T19:29:10Z       public unspecified    40.01
    ## 19 2013-04-17T17:01:50Z 2013-04-16T19:29:11Z       public unspecified    40.01
    ## 25 2013-04-17T17:01:50Z 2013-04-16T19:29:11Z       public unspecified    40.01
    ##    library_construction_protocol..exp.      Library.Name LibraryLayout
    ## 3                                 none   1075.Tongue:130        SINGLE
    ## 7                                 none   1053.Tongue:125        SINGLE
    ## 12                                none   05.HB.Tongue:56        SINGLE
    ## 15                                none  18.HB.Tongue:289        SINGLE
    ## 19                                none 68.HB.Tongue:1023        SINGLE
    ## 25                                none  53.HB.Tongue:796        SINGLE
    ##    LibrarySelection LibrarySource longitude num_sequences..exp.
    ## 3               PCR   METAGENOMIC   -105.27                6273
    ## 7               PCR   METAGENOMIC   -105.27               16875
    ## 12              PCR   METAGENOMIC   -105.27               53655
    ## 15              PCR   METAGENOMIC   -105.27               35654
    ## 19              PCR   METAGENOMIC   -105.27               36186
    ## 25              PCR   METAGENOMIC   -105.27               25889
    ##                 Organism
    ## 3  human oral metagenome
    ## 7  human oral metagenome
    ## 12 human oral metagenome
    ## 15 human oral metagenome
    ## 19 human oral metagenome
    ## 25 human oral metagenome
    ##                                                        pcr_primers..exp.
    ## 3  FWD:GCCTTGCCAGCCCGCTCAGTCAGAGTTTGATCCTGGCTCAG;REV:TGCTGCCTCCCGTAGGAGT
    ## 7  FWD:GCCTTGCCAGCCCGCTCAGTCAGAGTTTGATCCTGGCTCAG;REV:TGCTGCCTCCCGTAGGAGT
    ## 12 FWD:GCCTTGCCAGCCCGCTCAGTCAGAGTTTGATCCTGGCTCAG;REV:TGCTGCCTCCCGTAGGAGT
    ## 15 FWD:GCCTTGCCAGCCCGCTCAGTCAGAGTTTGATCCTGGCTCAG;REV:TGCTGCCTCCCGTAGGAGT
    ## 19 FWD:GCCTTGCCAGCCCGCTCAGTCAGAGTTTGATCCTGGCTCAG;REV:TGCTGCCTCCCGTAGGAGT
    ## 25 FWD:GCCTTGCCAGCCCGCTCAGTCAGAGTTTGATCCTGGCTCAG;REV:TGCTGCCTCCCGTAGGAGT
    ##    platform..exp. Platform          ReleaseDate    role row_number..exp.
    ## 3        Illumina ILLUMINA 2013-04-26T00:00:00Z  Single              130
    ## 7        Illumina ILLUMINA 2013-04-26T00:00:00Z  Single              125
    ## 12       Illumina ILLUMINA 2013-04-26T00:00:00Z Partner               56
    ## 15       Illumina ILLUMINA 2013-04-26T00:00:00Z Partner              289
    ## 19       Illumina ILLUMINA 2013-04-26T00:00:00Z Partner             1023
    ## 25       Illumina ILLUMINA 2013-04-26T00:00:00Z Partner              796
    ##    run_center..exp. run_date..exp. samp_size sample_center..exp. Sample_ID
    ## 3              CCME      4/17/2011    .1\\,g                CCME    321853
    ## 7              CCME      4/17/2011    .1\\,g                CCME    321857
    ## 12             CCME      4/17/2011    .1\\,g                CCME    321862
    ## 15             CCME      4/17/2011    .1\\,g                CCME    321865
    ## 19             CCME      4/17/2011    .1\\,g                CCME    321869
    ## 25             CCME      4/17/2011    .1\\,g                CCME    321875
    ##     Sample.Name                               Sample_name sequencing_meth..exp.
    ## 3  SAMEA1881262   1075.Tongue,qiime_study_797:1075.Tongue        pyrosequencing
    ## 7  SAMEA1878909   1053.Tongue,qiime_study_797:1053.Tongue        pyrosequencing
    ## 12 SAMEA1879886 05.HB.Tongue,qiime_study_797:05.HB.Tongue        pyrosequencing
    ## 15 SAMEA1880711 18.HB.Tongue,qiime_study_797:18.HB.Tongue        pyrosequencing
    ## 19 SAMEA1880524 68.HB.Tongue,qiime_study_797:68.HB.Tongue        pyrosequencing
    ## 25 SAMEA1878915 53.HB.Tongue,qiime_study_797:53.HB.Tongue        pyrosequencing
    ##    SRA.Study study_center..exp.                 Submitter_Id target_gene..exp.
    ## 3  ERP002478               CCME  qiime_study_797:1075.Tongue          16S rRNA
    ## 7  ERP002478               CCME  qiime_study_797:1053.Tongue          16S rRNA
    ## 12 ERP002478               CCME qiime_study_797:05.HB.Tongue          16S rRNA
    ## 15 ERP002478               CCME qiime_study_797:18.HB.Tongue          16S rRNA
    ## 19 ERP002478               CCME qiime_study_797:68.HB.Tongue          16S rRNA
    ## 25 ERP002478               CCME qiime_study_797:53.HB.Tongue          16S rRNA
    ##    target_subfragment..exp. taxon_id                 family_relationship
    ## 3                        V2   447426                                none
    ## 7                        V2   447426                                none
    ## 12                       V2   447426   father of FS_09, partner of FS_11
    ## 15                       V2   447426                    partner of FS_27
    ## 19                       V2   447426 father of FS_142, partner of FS_140
    ## 25                       V2   447426                   partner of FS_104
    ##    agegroup pet_farm_animal          create_date version PUBLIC
    ## 3    Senior            none 2013-04-17T21:57:00Z       1      n
    ## 7    Senior            none 2013-04-17T21:58:00Z       1      n
    ## 12    Adult            none 2013-04-17T21:58:00Z       1      n
    ## 15    Adult             dog 2013-04-17T21:58:00Z       1      n
    ## 19    Adult            none 2013-04-17T22:16:00Z       1      n
    ## 25    Adult             dog 2013-04-17T22:16:00Z       1      n
    ##    INSDC_center_alias weight_loss_3_month sex
    ## 3       CCME-COLORADO                   n  48
    ## 7       CCME-COLORADO                   y  48
    ## 12      CCME-COLORADO                   n  47
    ## 15      CCME-COLORADO                   n  47
    ## 19      CCME-COLORADO                   n  47
    ## 25      CCME-COLORADO                   n  47
    ##                                      drug_usage         Diet AGE
    ## 3                                          none regular.diet  69
    ## 7  Tums (occasionally), glucosamine/chondroitin         diet  78
    ## 12                                         none         diet  35
    ## 15                                            0 regular.diet  46
    ## 19                                         none regular.diet  29
    ## 25                                    Ibuprofen    no.gluten  40
    ##    Collection_Date tot_mass sub_species public scientific_name
    ## 3                        NA                                   
    ## 7                    60.702                                   
    ## 12         1/30/10   58.890                                   
    ## 15         2/21/10   95.130                                   
    ## 19         6/11/10   86.070                                   
    ## 25         8/12/10   67.950

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
  sample_name = initial_meta$anonymized_name,
  Age = initial_meta$AGE,
  BMI = NA,
  Sex = initial_meta$sex,
  Region = initial_meta$target_subfragment..exp.,
  Platform = initial_meta$platform..exp.,
  Country = initial_meta$geo_loc_name_country,
  Site = "tongue",
  public_time = "2013-04-26",
  PCR_primers = initial_meta$pcr_primers..exp.,
  stringsAsFactors = FALSE
)
```

``` r
metadata <- metadata %>%
  mutate(Sex = recode(Sex, "47" = "M", "48" = "F"))
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$title <- "Prevalence and Source of Fecal and Oral Bacteria on Infant, Child, and Adult Hands"
metadata$Disease_study <- "NO"
# Assign disease status based on sample_name patterns
metadata$Status <- "HC"

#extract subject id by remove prefix of sample_name


# Check initial data distribution
cat("Initial sample distribution:\n")
```

    ## Initial sample distribution:

``` r
table(metadata$sample_name)
```

    ## 
    ##  FS_01  FS_03  FS_04  FS_05  FS_06  FS_07  FS_08  FS_09  FS_10 FS_100 FS_101 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## FS_102 FS_103 FS_104 FS_105 FS_106 FS_107 FS_108 FS_109  FS_11 FS_110 FS_111 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## FS_112 FS_113 FS_114 FS_115 FS_116 FS_117 FS_118 FS_119  FS_12 FS_120 FS_121 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## FS_122 FS_123 FS_124 FS_125 FS_126 FS_127 FS_128 FS_129  FS_13 FS_130 FS_131 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## FS_132 FS_133 FS_134 FS_135 FS_136 FS_137 FS_138 FS_139  FS_14 FS_140 FS_141 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## FS_142 FS_143 FS_144 FS_145 FS_146 FS_147  FS_15 FS_151 FS_152 FS_153 FS_154 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## FS_155 FS_156 FS_157 FS_158 FS_159  FS_16 FS_160 FS_161 FS_162 FS_163 FS_164 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## FS_165 FS_166 FS_167 FS_168 FS_169  FS_17 FS_171 FS_172 FS_173 FS_174 FS_175 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  FS_18  FS_19  FS_20  FS_21  FS_22  FS_23  FS_24  FS_25  FS_26  FS_27  FS_28 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  FS_29  FS_30  FS_31  FS_32  FS_33  FS_34  FS_35  FS_36  FS_37  FS_38  FS_39 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  FS_40  FS_41  FS_42  FS_43  FS_44  FS_45  FS_46  FS_47  FS_48  FS_49  FS_50 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  FS_51  FS_52  FS_53  FS_54  FS_55  FS_56  FS_57  FS_58  FS_59  FS_60  FS_61 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  FS_62  FS_63  FS_64  FS_65  FS_66  FS_67  FS_68  FS_69  FS_70  FS_71  FS_72 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  FS_73  FS_74  FS_75  FS_76  FS_77  FS_78  FS_79  FS_80  FS_81  FS_82  FS_83 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  FS_84  FS_85  FS_86  FS_87  FS_88  FS_89  FS_90  FS_91  FS_92  FS_93  FS_94 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  FS_95  FS_96  FS_97  FS_98  FS_99 
    ##      1      1      1      1      1

``` r
table(metadata$Site)
```

    ## 
    ## tongue 
    ##    170

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
    ##  FS_01  FS_03  FS_04  FS_05  FS_06  FS_07  FS_08  FS_09  FS_10 FS_100 FS_101 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## FS_102 FS_103 FS_104 FS_105 FS_106 FS_107 FS_108 FS_109  FS_11 FS_110 FS_111 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## FS_112 FS_113 FS_114 FS_115 FS_116 FS_117 FS_118 FS_119  FS_12 FS_120 FS_121 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## FS_122 FS_123 FS_124 FS_125 FS_126 FS_127 FS_128 FS_129  FS_13 FS_130 FS_131 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## FS_132 FS_133 FS_134 FS_135 FS_136 FS_137 FS_138 FS_139  FS_14 FS_140 FS_141 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## FS_142 FS_143 FS_144 FS_145 FS_146 FS_147  FS_15 FS_151 FS_152 FS_153 FS_154 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## FS_155 FS_156 FS_157 FS_158 FS_159  FS_16 FS_160 FS_161 FS_162 FS_163 FS_164 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ## FS_165 FS_166 FS_167 FS_168 FS_169  FS_17 FS_171 FS_172 FS_173 FS_174 FS_175 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  FS_18  FS_19  FS_20  FS_21  FS_22  FS_23  FS_24  FS_25  FS_26  FS_27  FS_28 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  FS_29  FS_30  FS_31  FS_32  FS_33  FS_34  FS_35  FS_36  FS_37  FS_38  FS_39 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  FS_40  FS_41  FS_42  FS_43  FS_44  FS_45  FS_46  FS_47  FS_48  FS_49  FS_50 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  FS_51  FS_52  FS_53  FS_54  FS_55  FS_56  FS_57  FS_58  FS_59  FS_60  FS_61 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  FS_62  FS_63  FS_64  FS_65  FS_66  FS_67  FS_68  FS_69  FS_70  FS_71  FS_72 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  FS_73  FS_74  FS_75  FS_76  FS_77  FS_78  FS_79  FS_80  FS_81  FS_82  FS_83 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  FS_84  FS_85  FS_86  FS_87  FS_88  FS_89  FS_90  FS_91  FS_92  FS_93  FS_94 
    ##      1      1      1      1      1      1      1      1      1      1      1 
    ##  FS_95  FS_96  FS_97  FS_98  FS_99 
    ##      1      1      1      1      1

``` r
# Create matched_samples column
metadata$matched_samples <- ifelse(metadata$sample_name %in% names(sample_counts[sample_counts >= 2]), "YES", "NO")
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 170"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 170"

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
write.csv(metadata, file = "metadata.csv", quote = TRUE, row.names = FALSE)
head(metadata)
```

    ##   BioProject       Run sample_name Age BMI Sex Region Platform Country   Site
    ## 1  PRJEB1799 ERR261871      FS_174  69  NA   F     V2 Illumina     USA tongue
    ## 2  PRJEB1799 ERR261875      FS_153  78  NA   F     V2 Illumina     USA tongue
    ## 3  PRJEB1799 ERR261880       FS_10  35  NA   M     V2 Illumina     USA tongue
    ## 4  PRJEB1799 ERR261883       FS_28  46  NA   M     V2 Illumina     USA tongue
    ## 5  PRJEB1799 ERR261887      FS_141  29  NA   M     V2 Illumina     USA tongue
    ## 6  PRJEB1799 ERR261893      FS_105  40  NA   M     V2 Illumina     USA tongue
    ##   public_time
    ## 1  2013-04-26
    ## 2  2013-04-26
    ## 3  2013-04-26
    ## 4  2013-04-26
    ## 5  2013-04-26
    ## 6  2013-04-26
    ##                                                             PCR_primers
    ## 1 FWD:GCCTTGCCAGCCCGCTCAGTCAGAGTTTGATCCTGGCTCAG;REV:TGCTGCCTCCCGTAGGAGT
    ## 2 FWD:GCCTTGCCAGCCCGCTCAGTCAGAGTTTGATCCTGGCTCAG;REV:TGCTGCCTCCCGTAGGAGT
    ## 3 FWD:GCCTTGCCAGCCCGCTCAGTCAGAGTTTGATCCTGGCTCAG;REV:TGCTGCCTCCCGTAGGAGT
    ## 4 FWD:GCCTTGCCAGCCCGCTCAGTCAGAGTTTGATCCTGGCTCAG;REV:TGCTGCCTCCCGTAGGAGT
    ## 5 FWD:GCCTTGCCAGCCCGCTCAGTCAGAGTTTGATCCTGGCTCAG;REV:TGCTGCCTCCCGTAGGAGT
    ## 6 FWD:GCCTTGCCAGCCCGCTCAGTCAGAGTTTGATCCTGGCTCAG;REV:TGCTGCCTCCCGTAGGAGT
    ##                                                                                title
    ## 1 Prevalence and Source of Fecal and Oral Bacteria on Infant, Child, and Adult Hands
    ## 2 Prevalence and Source of Fecal and Oral Bacteria on Infant, Child, and Adult Hands
    ## 3 Prevalence and Source of Fecal and Oral Bacteria on Infant, Child, and Adult Hands
    ## 4 Prevalence and Source of Fecal and Oral Bacteria on Infant, Child, and Adult Hands
    ## 5 Prevalence and Source of Fecal and Oral Bacteria on Infant, Child, and Adult Hands
    ## 6 Prevalence and Source of Fecal and Oral Bacteria on Infant, Child, and Adult Hands
    ##   Disease_study Status matched_samples
    ## 1            NO     HC              NO
    ## 2            NO     HC              NO
    ## 3            NO     HC              NO
    ## 4            NO     HC              NO
    ## 5            NO     HC              NO
    ## 6            NO     HC              NO

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
