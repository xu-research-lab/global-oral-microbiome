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
library(stringr)

# Set working directory (update this path to your local directory)
setwd("/home/hby/oral_project/PRJNA416501")

# Read the initial metadata files
initial_meta <- read.csv("SraRunTable.csv", header = T)

# Combine the two metadata files
# meta <- rbind(initial_meta, initial_meta1)

# Remove entries non-oral samples
# initial_meta <- subset(initial_meta,isolate %in% c("saliva","tongue coat"))

# Check the number of samples containing "HT" in Sample.Name (e.g., healthy controls)
# nrow(subset(meta, grepl("HT", Sample.Name)))

# Display the first few rows of the metadata
head(initial_meta)
```

    ##          Run Assay.Type                                                                    assembly AvgSpotLen    Bases  BioProject    BioSample                       BioSampleModel
    ## 1 SRR6239859   AMPLICON Flash software;error rate of 1 in 100 bp;overlapping and merging read pairs        447 10859358 PRJNA416501 SAMN07958299 MIMS.me,MIGS/MIMS/MIMARKS.human-oral
    ## 2 SRR6239860   AMPLICON Flash software;error rate of 1 in 100 bp;overlapping and merging read pairs        446  9364855 PRJNA416501 SAMN07958300 MIMS.me,MIGS/MIMS/MIMARKS.human-oral
    ## 3 SRR6239861   AMPLICON Flash software;error rate of 1 in 100 bp;overlapping and merging read pairs        451 11344944 PRJNA416501 SAMN07958295 MIMS.me,MIGS/MIMS/MIMARKS.human-oral
    ## 4 SRR6239862   AMPLICON Flash software;error rate of 1 in 100 bp;overlapping and merging read pairs        448 19693298 PRJNA416501 SAMN07958304 MIMS.me,MIGS/MIMS/MIMARKS.human-oral
    ## 5 SRR6239863   AMPLICON Flash software;error rate of 1 in 100 bp;overlapping and merging read pairs        446 18969416 PRJNA416501 SAMN07958306 MIMS.me,MIGS/MIMS/MIMARKS.human-oral
    ## 6 SRR6239864   AMPLICON Flash software;error rate of 1 in 100 bp;overlapping and merging read pairs        452 22684384 PRJNA416501 SAMN07958293 MIMS.me,MIGS/MIMS/MIMARKS.human-oral
    ##     Bytes            Center.Name Collection_Date Consent DATASTORE.filetype DATASTORE.provider                     DATASTORE.region    env_biome env_feature     env_material env_package
    ## 1 2752848 UNIVERSITE DE RENNES 1       2011/5/19  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 Homo sapiens     biofilm biofilm material  human-oral
    ## 2 2394750 UNIVERSITE DE RENNES 1       2011/5/19  public       fastq,run.zq         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 Homo sapiens     biofilm biofilm material  human-oral
    ## 3 2899243 UNIVERSITE DE RENNES 1       2011/3/24  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 Homo sapiens     biofilm biofilm material  human-oral
    ## 4 5065298 UNIVERSITE DE RENNES 1       2011/7/28  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 Homo sapiens     biofilm biofilm material  human-oral
    ## 5 4787809 UNIVERSITE DE RENNES 1        2011/9/9  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 Homo sapiens     biofilm biofilm material  human-oral
    ## 6 5809392 UNIVERSITE DE RENNES 1       2011/3/22  public   fastq,run.zq,sra         gs,ncbi,s3 gs.us-east1,ncbi.public,s3.us-east-1 Homo sapiens     biofilm biofilm material  human-oral
    ##   ethnicity Experiment geo_loc_name_country geo_loc_name_country_continent             geo_loc_name Host_Age  host_body_site      host_disease_stat         HOST host_sex host_subject_id
    ## 1 Caucasian SRX3347299               France                         Europe France:Brittany\\,Rennes       54 gingival groove hemochromatosis type 1 Homo sapiens     male              76
    ## 2 Caucasian SRX3347298               France                         Europe France:Brittany\\,Rennes       43 gingival groove hemochromatosis type 1 Homo sapiens     male              77
    ## 3 Caucasian SRX3347297               France                         Europe France:Brittany\\,Rennes       59 gingival groove hemochromatosis type 1 Homo sapiens     male              69
    ## 4 Caucasian SRX3347296               France                         Europe France:Brittany\\,Rennes       43 gingival groove hemochromatosis type 1 Homo sapiens     male              82
    ## 5 Caucasian SRX3347295               France                         Europe France:Brittany\\,Rennes       63 gingival groove hemochromatosis type 1 Homo sapiens   female              86
    ## 6 Caucasian SRX3347294               France                         Europe France:Brittany\\,Rennes       41 gingival groove hemochromatosis type 1 Homo sapiens     male              67
    ##       Instrument investigation_type                lat_lon lib_const_meth lib_reads_seqd Library.Name LibraryLayout LibrarySelection LibrarySource mean_bleeding_on_probing mean_cal
    ## 1 Illumina MiSeq         metagenome 48.119547 N 1.694538 W     paired-end          24265    HP_16S_76        PAIRED              PCR   METAGENOMIC                     0.88     3.60
    ## 2 Illumina MiSeq         metagenome 48.119547 N 1.694538 W     paired-end          20984    HP_16S_77        PAIRED              PCR   METAGENOMIC                     0.30     2.15
    ## 3 Illumina MiSeq         metagenome 48.119547 N 1.694538 W     paired-end          25152    HP_16S_69        PAIRED              PCR   METAGENOMIC                     0.44     2.86
    ## 4 Illumina MiSeq         metagenome 48.119547 N 1.694538 W     paired-end          43950    HP_16S_82        PAIRED              PCR   METAGENOMIC                     0.38     2.90
    ## 5 Illumina MiSeq         metagenome 48.119547 N 1.694538 W     paired-end          42486    HP_16S_86        PAIRED              PCR   METAGENOMIC                     0.22     2.52
    ## 6 Illumina MiSeq         metagenome 48.119547 N 1.694538 W     paired-end          50181    HP_16S_67        PAIRED              PCR   METAGENOMIC                     0.52     2.88
    ##   mean_gingival_index mean_plaque_index mean_ppd number_of_teeth              Organism                                                                                       pcr_cond
    ## 1                0.57              0.38     2.99              24 human oral metagenome initial denaturation:98_0.5; annealing:45_0.5; elongation: 72_0.75; final elongation:72_10; 25
    ## 2                0.33              0.22     1.91              27 human oral metagenome initial denaturation:98_0.5; annealing:45_0.5; elongation: 72_0.75; final elongation:72_10; 25
    ## 3                0.43              0.22     2.45              27 human oral metagenome initial denaturation:98_0.5; annealing:45_0.5; elongation: 72_0.75; final elongation:72_10; 25
    ## 4                0.26              0.38     2.06              26 human oral metagenome initial denaturation:98_0.5; annealing:45_0.5; elongation: 72_0.75; final elongation:72_10; 25
    ## 5                0.19              0.00     2.02              27 human oral metagenome initial denaturation:98_0.5; annealing:45_0.5; elongation: 72_0.75; final elongation:72_10; 25
    ## 6                0.15              0.15     2.18              27 human oral metagenome initial denaturation:98_0.5; annealing:45_0.5; elongation: 72_0.75; final elongation:72_10; 25
    ##                                       PCR_primers periodontitis_classification_eke_2012 Platform project_name          ReleaseDate samp_collect_device samp_store_temp Sample_Name
    ## 1 FWD:ACTCCTACGGGAGGCAGCAG REV:TACNVGGGTATCTAATCC                                severe ILLUMINA     Hemoparo 2017-10-31T00:00:00Z sterile paper point             -80  Hemoparo76
    ## 2 FWD:ACTCCTACGGGAGGCAGCAG REV:TACNVGGGTATCTAATCC                              moderate ILLUMINA     Hemoparo 2017-10-31T00:00:00Z sterile paper point             -80  Hemoparo77
    ## 3 FWD:ACTCCTACGGGAGGCAGCAG REV:TACNVGGGTATCTAATCC                              moderate ILLUMINA     Hemoparo 2017-10-31T00:00:00Z sterile paper point             -80  Hemoparo69
    ## 4 FWD:ACTCCTACGGGAGGCAGCAG REV:TACNVGGGTATCTAATCC                              moderate ILLUMINA     Hemoparo 2017-10-31T00:00:00Z sterile paper point             -80  Hemoparo82
    ## 5 FWD:ACTCCTACGGGAGGCAGCAG REV:TACNVGGGTATCTAATCC                              moderate ILLUMINA     Hemoparo 2017-10-31T00:00:00Z sterile paper point             -80  Hemoparo86
    ## 6 FWD:ACTCCTACGGGAGGCAGCAG REV:TACNVGGGTATCTAATCC                                severe ILLUMINA     Hemoparo 2017-10-31T00:00:00Z sterile paper point             -80  Hemoparo67
    ##         seq_meth          seq_quality_check serum_ferritin_micrograms_per_liter sites_with_cal_greater_or_equal_3mm_percentage sites_with_cal_greater_or_equal_5mm_percentage
    ## 1 Illumina MiSeq manually edited using PhiX                                  37                                          75.69                                          20.14
    ## 2 Illumina MiSeq manually edited using PhiX                                  36                                          30.25                                           2.47
    ## 3 Illumina MiSeq manually edited using PhiX                                 240                                          61.11                                           4.94
    ## 4 Illumina MiSeq manually edited using PhiX                                1022                                          57.69                                           8.97
    ## 5 Illumina MiSeq manually edited using PhiX                                  43                                          44.44                                           8.02
    ## 6 Illumina MiSeq manually edited using PhiX                                 873                                          51.85                                          11.73
    ##   sites_with_ppd_greater_or_equal_4mm_percentage sites_with_ppd_greater_or_equal_6mm_percentage SRA.Study target_gene target_subfragment transferrin_saturation_percentage
    ## 1                                          34.72                                           2.08 SRP123221    16S rRNA               V3V4                              56.0
    ## 2                                           7.41                                           0.00 SRP123221    16S rRNA               V3V4                              30.9
    ## 3                                           8.64                                           0.00 SRP123221    16S rRNA               V3V4                              71.5
    ## 4                                           3.21                                           0.00 SRP123221    16S rRNA               V3V4                              66.4
    ## 5                                           8.64                                           0.00 SRP123221    16S rRNA               V3V4                              71.0
    ## 6                                          11.11                                           0.00 SRP123221    16S rRNA               V3V4                              61.5
    ##            create_date version host_body_mass_index host_height host_weight
    ## 1 2017-10-31T19:05:00Z       1                 20.0         176          62
    ## 2                           NA                 28.1         178          89
    ## 3 2017-10-31T19:06:00Z       1                   NA          NA          NA
    ## 4 2017-10-31T19:06:00Z       1                 29.5         185         101
    ## 5 2017-10-31T19:06:00Z       1                   NA          NA          NA
    ## 6 2017-10-31T19:06:00Z       1                   NA          NA          NA

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
  sample_name = initial_meta$Sample_Name,
  Platform = initial_meta$Instrument,
  Country = initial_meta$geo_loc_name_country,
  public_time = initial_meta$ReleaseDate,
  Site = initial_meta$host_body_site,
  Sex = initial_meta$host_sex,
  Age = initial_meta$Host_Age,
  BMI = initial_meta$host_body_mass_index
) %>%
  mutate(
    Sex = case_when(
      Sex == "female" ~ "F",
      Sex == "male" ~ "M",
    )
  )
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$Status <- "hemochromatosis、periodontitis"
metadata$Region <- "V3-V4"
metadata$PCR_primers <- "338F-802R"
metadata$title <- "Periodontal pathogens and clinical parameters in chronic periodontitis"
metadata$Disease_study <- "YES"
metadata$matched_samples <- "NO"
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 68"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 68"

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
  "PCR_primers", "title", "Disease_study", "Status", "matched_samples"
)]

# Export the uniform metadata to a CSV file
write.csv(metadata, file = "metadata.csv", quote = FALSE)
head(metadata)
```

    ##    BioProject        Run sample_name Age  BMI Sex Region       Platform Country            Site          public_time PCR_primers
    ## 1 PRJNA416501 SRR6239859  Hemoparo76  54 20.0   M  V3-V4 Illumina MiSeq  France gingival groove 2017-10-31T00:00:00Z   338F-802R
    ## 2 PRJNA416501 SRR6239860  Hemoparo77  43 28.1   M  V3-V4 Illumina MiSeq  France gingival groove 2017-10-31T00:00:00Z   338F-802R
    ## 3 PRJNA416501 SRR6239861  Hemoparo69  59   NA   M  V3-V4 Illumina MiSeq  France gingival groove 2017-10-31T00:00:00Z   338F-802R
    ## 4 PRJNA416501 SRR6239862  Hemoparo82  43 29.5   M  V3-V4 Illumina MiSeq  France gingival groove 2017-10-31T00:00:00Z   338F-802R
    ## 5 PRJNA416501 SRR6239863  Hemoparo86  63   NA   F  V3-V4 Illumina MiSeq  France gingival groove 2017-10-31T00:00:00Z   338F-802R
    ## 6 PRJNA416501 SRR6239864  Hemoparo67  41   NA   M  V3-V4 Illumina MiSeq  France gingival groove 2017-10-31T00:00:00Z   338F-802R
    ##                                                                    title Disease_study                         Status matched_samples
    ## 1 Periodontal pathogens and clinical parameters in chronic periodontitis           YES hemochromatosis、periodontitis              NO
    ## 2 Periodontal pathogens and clinical parameters in chronic periodontitis           YES hemochromatosis、periodontitis              NO
    ## 3 Periodontal pathogens and clinical parameters in chronic periodontitis           YES hemochromatosis、periodontitis              NO
    ## 4 Periodontal pathogens and clinical parameters in chronic periodontitis           YES hemochromatosis、periodontitis              NO
    ## 5 Periodontal pathogens and clinical parameters in chronic periodontitis           YES hemochromatosis、periodontitis              NO
    ## 6 Periodontal pathogens and clinical parameters in chronic periodontitis           YES hemochromatosis、periodontitis              NO

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
