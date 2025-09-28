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
setwd("E:/bioproject/PRJEB86033")

# Read the initial metadata files
initial_meta <- read.csv("SraRunTable.csv", header = T) 

# Display the first few rows of the metadata
head(initial_meta)
```

    ##           Run anonymized_name Assay.Type AvgSpotLen    Bases BioProject
    ## 1 ERR14648080        SAMPLE10   AMPLICON        272 14071648 PRJEB86033
    ## 2 ERR14648081       SAMPLE100   AMPLICON        272 12905040 PRJEB86033
    ## 3 ERR14648082       SAMPLE101   AMPLICON        272 14947488 PRJEB86033
    ## 4 ERR14648083       SAMPLE102   AMPLICON        272 13282576 PRJEB86033
    ## 5 ERR14648085       SAMPLE104   AMPLICON        272 12982288 PRJEB86033
    ## 6 ERR14648086       SAMPLE105   AMPLICON        272 12212528 PRJEB86033
    ##        BioSample    Bytes center_name..exp.
    ## 1 SAMEA117741011  9622457 HKU + QIBEBT, CAS
    ## 2 SAMEA117741012  8343713 HKU + QIBEBT, CAS
    ## 3 SAMEA117741013 11073876 HKU + QIBEBT, CAS
    ## 4 SAMEA117741014  8584287 HKU + QIBEBT, CAS
    ## 5 SAMEA117741016  9618330 HKU + QIBEBT, CAS
    ## 6 SAMEA117741017  9027042 HKU + QIBEBT, CAS
    ##                                          Center.Name Collection_Date
    ## 1 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT            2015
    ## 2 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT            2015
    ## 3 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT            2015
    ## 4 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT            2015
    ## 5 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT            2015
    ## 6 UNIVERSITY OF CALIFORNIA SAN DIEGO MICROBIOME INIT            2015
    ##   collection_timestamp Consent   COUNTRY DATASTORE.filetype DATASTORE.provider
    ## 1            2015-2016  public Hong Kong   fastq,run.zq,sra         gs,ncbi,s3
    ## 2            2015-2016  public Hong Kong   fastq,run.zq,sra         gs,ncbi,s3
    ## 3            2015-2016  public Hong Kong   fastq,run.zq,sra         gs,ncbi,s3
    ## 4            2015-2016  public Hong Kong   fastq,run.zq,sra         gs,ncbi,s3
    ## 5            2015-2016  public Hong Kong   fastq,run.zq,sra         gs,ncbi,s3
    ## 6            2015-2016  public Hong Kong   fastq,run.zq,sra         gs,ncbi,s3
    ##                       DATASTORE.region dmfs_host dmfs_tooth dna_extracted
    ## 1 gs.us-east1,ncbi.public,s3.us-east-1         0          0          TRUE
    ## 2 gs.us-east1,ncbi.public,s3.us-east-1         1          0          TRUE
    ## 3 gs.us-east1,ncbi.public,s3.us-east-1         4          0          TRUE
    ## 4 gs.us-east1,ncbi.public,s3.us-east-1         1          0          TRUE
    ## 5 gs.us-east1,ncbi.public,s3.us-east-1         1          0          TRUE
    ## 6 gs.us-east1,ncbi.public,s3.us-east-1         4          0          TRUE
    ##   elevation          empo_1                       empo_2              empo_3
    ## 1         0 Host-associated Host-associated (non-saline) Animal (non-saline)
    ## 2         0 Host-associated Host-associated (non-saline) Animal (non-saline)
    ## 3         0 Host-associated Host-associated (non-saline) Animal (non-saline)
    ## 4         0 Host-associated Host-associated (non-saline) Animal (non-saline)
    ## 5         0 Host-associated Host-associated (non-saline) Animal (non-saline)
    ## 6         0 Host-associated Host-associated (non-saline) Animal (non-saline)
    ##                          empo_4 ENA.FIRST.PUBLIC..run. ENA_first_public
    ## 1 Animal secretion (non-saline)              2025/2/21        2025/2/21
    ## 2 Animal secretion (non-saline)              2025/2/21        2025/2/21
    ## 3 Animal secretion (non-saline)              2025/2/21        2025/2/21
    ## 4 Animal secretion (non-saline)              2025/2/21        2025/2/21
    ## 5 Animal secretion (non-saline)              2025/2/21        2025/2/21
    ## 6 Animal secretion (non-saline)              2025/2/21        2025/2/21
    ##   ENA.LAST.UPDATE..run.   env_biome              env_feature     env_material
    ## 1             2025/2/21 urban biome human-associated habitat organic material
    ## 2             2025/2/21 urban biome human-associated habitat organic material
    ## 3             2025/2/21 urban biome human-associated habitat organic material
    ## 4             2025/2/21 urban biome human-associated habitat organic material
    ## 5             2025/2/21 urban biome human-associated habitat organic material
    ## 6             2025/2/21 urban biome human-associated habitat organic material
    ##   env_package  Experiment future_status_tooth gender geo_loc_name_country
    ## 1  human-oral ERX14049208          ConfidentH      B            Hong Kong
    ## 2  human-oral ERX14049209                 C_H      B            Hong Kong
    ## 3  human-oral ERX14049210                 C_H      B            Hong Kong
    ## 4  human-oral ERX14049211                 C_H      B            Hong Kong
    ## 5  human-oral ERX14049213                 C_H      B            Hong Kong
    ## 6  human-oral ERX14049214                 C_H      B            Hong Kong
    ##   geo_loc_name_country_continent geo_loc_name
    ## 1                           Asia    Hong Kong
    ## 2                           Asia    Hong Kong
    ## 3                           Asia    Hong Kong
    ## 4                           Asia    Hong Kong
    ## 5                           Asia    Hong Kong
    ## 6                           Asia    Hong Kong
    ##   geographic_location_.country_and.or_sea.     Host_age host_age_units
    ## 1                                Hong Kong not provided          years
    ## 2                                Hong Kong not provided          years
    ## 3                                Hong Kong not provided          years
    ## 4                                Hong Kong not provided          years
    ## 5                                Hong Kong not provided          years
    ## 6                                Hong Kong not provided          years
    ##    host_body_habitat host_body_mass_index    host_body_product host_body_site
    ## 1 UBERON:oral cavity       not applicable UBERON:dental plaque   UBERON:tooth
    ## 2 UBERON:oral cavity       not applicable UBERON:dental plaque   UBERON:tooth
    ## 3 UBERON:oral cavity       not applicable UBERON:dental plaque   UBERON:tooth
    ## 4 UBERON:oral cavity       not applicable UBERON:dental plaque   UBERON:tooth
    ## 5 UBERON:oral cavity       not applicable UBERON:dental plaque   UBERON:tooth
    ## 6 UBERON:oral cavity       not applicable UBERON:dental plaque   UBERON:tooth
    ##   host_common_name  host_height host_height_units host_life_stage
    ## 1            human not provided                cm        juvenile
    ## 2            human not provided                cm        juvenile
    ## 3            human not provided                cm        juvenile
    ## 4            human not provided                cm        juvenile
    ## 5            human not provided                cm        juvenile
    ## 6            human not provided                cm        juvenile
    ##   host_scientific_name host_sex host_subject_id host_taxid  host_weight
    ## 1         Homo sapiens     male            A100       9606 not provided
    ## 2         Homo sapiens     male            A154       9606 not provided
    ## 3         Homo sapiens     male            A154       9606 not provided
    ## 4         Homo sapiens     male            A154       9606 not provided
    ## 5         Homo sapiens     male            A154       9606 not provided
    ## 6         Homo sapiens     male            A154       9606 not provided
    ##   host_weight_units hostgroup HostID INSDC_center_name INSDC_status
    ## 1                kg       H2H   A100            UCSDMI       public
    ## 2                kg       C2C   A154            UCSDMI       public
    ## 3                kg       C2C   A154            UCSDMI       public
    ## 4                kg       C2C   A154            UCSDMI       public
    ## 5                kg       C2C   A154            UCSDMI       public
    ## 6                kg       C2C   A154            UCSDMI       public
    ##       Instrument LATITUDE    Library.Name LibraryLayout LibrarySelection
    ## 1 Illumina MiSeq     22.3  14341.SAMPLE10        SINGLE              PCR
    ## 2 Illumina MiSeq     22.3 14341.SAMPLE100        SINGLE              PCR
    ## 3 Illumina MiSeq     22.3 14341.SAMPLE101        SINGLE              PCR
    ## 4 Illumina MiSeq     22.3 14341.SAMPLE102        SINGLE              PCR
    ## 5 Illumina MiSeq     22.3 14341.SAMPLE104        SINGLE              PCR
    ## 6 Illumina MiSeq     22.3 14341.SAMPLE105        SINGLE              PCR
    ##   LibrarySource longitude Month   nicheap nicheicm nichelr nicheul
    ## 1   METAGENOMIC     114.2    67 Posterior    Molar    Left   Upper
    ## 2   METAGENOMIC     114.2    49 Posterior    Molar    Left   Upper
    ## 3   METAGENOMIC     114.2    56 Posterior    Molar    Left   Upper
    ## 4   METAGENOMIC     114.2    49 Posterior    Molar    Left   Upper
    ## 5   METAGENOMIC     114.2    49 Posterior    Molar    Left   Lower
    ## 6   METAGENOMIC     114.2    56 Posterior    Molar    Left   Lower
    ##                Organism physical_specimen_location physical_specimen_remaining
    ## 1 human oral metagenome                        HKU                       FALSE
    ## 2 human oral metagenome                        HKU                       FALSE
    ## 3 human oral metagenome                        HKU                       FALSE
    ## 4 human oral metagenome                        HKU                       FALSE
    ## 5 human oral metagenome                        HKU                       FALSE
    ## 6 human oral metagenome                        HKU                       FALSE
    ##   Platform position position_tooth qiita_sample_type          ReleaseDate
    ## 1 ILLUMINA      T64             T4     dental plaque 2025-02-22T00:00:00Z
    ## 2 ILLUMINA      T64             T4     dental plaque 2025-02-22T00:00:00Z
    ## 3 ILLUMINA      T64             T4     dental plaque 2025-02-22T00:00:00Z
    ## 4 ILLUMINA      T65             T5     dental plaque 2025-02-22T00:00:00Z
    ## 5 ILLUMINA      T74             T4     dental plaque 2025-02-22T00:00:00Z
    ## 6 ILLUMINA      T74             T4     dental plaque 2025-02-22T00:00:00Z
    ##   run_center..exp.          create_date version    run_prefix..exp.
    ## 1        OE_biotec 2025-03-07T00:12:00Z       1  AWGB16SAA05317-5-7
    ## 2        OE_biotec 2025-03-03T18:38:00Z       1 AWGB16SAA05443-7-11
    ## 3        OE_biotec 2025-03-06T10:19:00Z       1  AWGB16SAA05452-8-7
    ## 4        OE_biotec 2025-03-06T10:18:00Z       1 AWGB16SAA05444-7-12
    ## 5        OE_biotec 2025-03-06T10:21:00Z       1  AWGB16SAA05445-8-1
    ## 6        OE_biotec 2025-03-06T10:22:00Z       1  AWGB16SAA05454-8-9
    ##      Sample.Name                     Sample_name   sample_type
    ## 1 SAMEA117741011  qiita_sid_14341:14341.SAMPLE10 dental plaque
    ## 2 SAMEA117741012 qiita_sid_14341:14341.SAMPLE100 dental plaque
    ## 3 SAMEA117741013 qiita_sid_14341:14341.SAMPLE101 dental plaque
    ## 4 SAMEA117741014 qiita_sid_14341:14341.SAMPLE102 dental plaque
    ## 5 SAMEA117741016 qiita_sid_14341:14341.SAMPLE104 dental plaque
    ## 6 SAMEA117741017 qiita_sid_14341:14341.SAMPLE105 dental plaque
    ##         Scientific_Name   sequencing_meth..exp.    sex SRA.Study status_host
    ## 1 human oral metagenome Sequencing by synthesis female ERP169403           H
    ## 2 human oral metagenome Sequencing by synthesis female ERP169403           C
    ## 3 human oral metagenome Sequencing by synthesis female ERP169403           C
    ## 4 human oral metagenome Sequencing by synthesis female ERP169403           C
    ## 5 human oral metagenome Sequencing by synthesis female ERP169403           C
    ## 6 human oral metagenome Sequencing by synthesis female ERP169403           C
    ##   status_tooth                    Submitter_Id timepoint         x         y
    ## 1            H  qiita_sid_14341:14341.SAMPLE10        T4 0.6200585 0.5288878
    ## 2            H qiita_sid_14341:14341.SAMPLE100        T3 0.6200585 0.5288878
    ## 3            H qiita_sid_14341:14341.SAMPLE101        T5 0.6200585 0.5288878
    ## 4            H qiita_sid_14341:14341.SAMPLE102        T3 0.6950439 0.2516267
    ## 5            H qiita_sid_14341:14341.SAMPLE104        T3 0.6250731 0.5235239
    ## 6            H qiita_sid_14341:14341.SAMPLE105        T5 0.6250731 0.5235239
    ##       z barcode..exp.
    ## 1  0.12         AAAAA
    ## 2  0.12         AAAAA
    ## 3  0.12         AAAAA
    ## 4  0.10         AAAAA
    ## 5 -0.12         AAAAA
    ## 6 -0.12         AAAAA
    ##                                                                        pcr_primers..exp.
    ## 1 Forward primer: 341F (CCTACGGGNGGCWGCAG), Reverse primer: 805R (GACTACHVGGGTATCTAATCC)
    ## 2 Forward primer: 341F (CCTACGGGNGGCWGCAG), Reverse primer: 805R (GACTACHVGGGTATCTAATCC)
    ## 3 Forward primer: 341F (CCTACGGGNGGCWGCAG), Reverse primer: 805R (GACTACHVGGGTATCTAATCC)
    ## 4 Forward primer: 341F (CCTACGGGNGGCWGCAG), Reverse primer: 805R (GACTACHVGGGTATCTAATCC)
    ## 5 Forward primer: 341F (CCTACGGGNGGCWGCAG), Reverse primer: 805R (GACTACHVGGGTATCTAATCC)
    ## 6 Forward primer: 341F (CCTACGGGNGGCWGCAG), Reverse primer: 805R (GACTACHVGGGTATCTAATCC)
    ##   target_gene..exp. target_subfragment..exp. center_project_name..exp.
    ## 1          16S rRNA             V3-V4 region                          
    ## 2          16S rRNA             V3-V4 region                          
    ## 3          16S rRNA             V3-V4 region                          
    ## 4          16S rRNA             V3-V4 region                          
    ## 5          16S rRNA             V3-V4 region                          
    ## 6          16S rRNA             V3-V4 region

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
  BMI = NA,
  Sex = initial_meta$sex,
  Region = "V3-V4",
  Platform = initial_meta$Instrument,
  Site = initial_meta$qiita_sample_type,
  public_time = initial_meta$ReleaseDate,
  Status = initial_meta$status_host,
  matched_samples = NA,
  stringsAsFactors = FALSE
)

# Remove all occurrences of "region" at the end of the Region column
# Use abbreviations to indicate gender
metadata <- metadata %>%
  mutate(Sex = recode(Sex, "female" = "F", "male" = "M"))
```

# Step 3: Manually Add Information from the Publication

Some information (e.g. PCR primers, study title, disease status) must be
manually extracted from the publication and added to the metadata.

``` r
# Add information obtained from the publication
metadata$Country <- "China"
metadata$PCR_primers <- "341F-805R"
metadata$title <- "Single-tooth resolved, whole-mouth prediction of early childhood caries
via spatiotemporal variations of plaque microbiota"
metadata$Disease_study <- "YES"

# Check initial data distribution
cat("Initial sample distribution:\n")
```

    ## Initial sample distribution:

``` r
table(metadata$sample_name)
```

    ## 
    ##      14341.A0001      14341.A0002      14341.A0003      14341.A0004 
    ##                1                1                1                1 
    ##      14341.A0005      14341.A0006      14341.A0007      14341.A0008 
    ##                1                1                1                1 
    ##      14341.A0009      14341.A0010      14341.A0011      14341.A0012 
    ##                1                1                1                1 
    ##      14341.A0013      14341.A0014      14341.A0015      14341.A0016 
    ##                1                1                1                1 
    ##      14341.A0017      14341.A0018      14341.A0019      14341.A0020 
    ##                1                1                1                1 
    ##      14341.A0021      14341.A0022      14341.A0023      14341.A0024 
    ##                1                1                1                1 
    ##      14341.A0025      14341.A0026      14341.A0027      14341.A0028 
    ##                1                1                1                1 
    ##      14341.A0029      14341.A0030      14341.A0031      14341.A0032 
    ##                1                1                1                1 
    ##      14341.A0033      14341.A0034      14341.A0035      14341.A0036 
    ##                1                1                1                1 
    ##      14341.A0037      14341.A0038      14341.A0039      14341.A0040 
    ##                1                1                1                1 
    ##      14341.A0041      14341.A0042      14341.A0043      14341.A0044 
    ##                1                1                1                1 
    ##      14341.A0045      14341.A0046      14341.A0047      14341.A0048 
    ##                1                1                1                1 
    ##      14341.A0049      14341.A0050      14341.A0051      14341.A0052 
    ##                1                1                1                1 
    ##      14341.A0053      14341.A0054      14341.A0055      14341.A0056 
    ##                1                1                1                1 
    ##      14341.A0057      14341.A0058      14341.A0059      14341.A0060 
    ##                1                1                1                1 
    ##      14341.A0061      14341.A0062      14341.A0063      14341.A0064 
    ##                1                1                1                1 
    ##      14341.A0065      14341.A0066      14341.A0067      14341.A0068 
    ##                1                1                1                1 
    ##      14341.A0069      14341.A0070      14341.A0071      14341.A0072 
    ##                1                1                1                1 
    ##      14341.A0073      14341.A0074      14341.A0075      14341.A0076 
    ##                1                1                1                1 
    ##      14341.A0077      14341.A0078      14341.A0079      14341.A0080 
    ##                1                1                1                1 
    ##      14341.A0081      14341.A0082      14341.A0083      14341.A0084 
    ##                1                1                1                1 
    ##      14341.A0085      14341.A0086      14341.A0087      14341.A0088 
    ##                1                1                1                1 
    ##      14341.A0089      14341.A0090      14341.A0091      14341.A0092 
    ##                1                1                1                1 
    ##      14341.A0093      14341.A0094      14341.A0095      14341.A0096 
    ##                1                1                1                1 
    ##      14341.A0097      14341.A0098      14341.A0099      14341.A0100 
    ##                1                1                1                1 
    ##      14341.A0101      14341.A0102      14341.A0103      14341.A0104 
    ##                1                1                1                1 
    ##      14341.A0105      14341.A0106      14341.A0107      14341.A0108 
    ##                1                1                1                1 
    ##      14341.A0109      14341.A0110      14341.A0111      14341.A0112 
    ##                1                1                1                1 
    ##      14341.A0113      14341.A0114      14341.A0115      14341.A0116 
    ##                1                1                1                1 
    ##      14341.A0117      14341.A0118      14341.A0119      14341.A0120 
    ##                1                1                1                1 
    ##      14341.A0121      14341.A0122      14341.A0123      14341.A0124 
    ##                1                1                1                1 
    ##      14341.A0125      14341.A0126      14341.A0127      14341.A0128 
    ##                1                1                1                1 
    ##      14341.A0129      14341.A0130      14341.A0131      14341.A0132 
    ##                1                1                1                1 
    ##      14341.A0133      14341.A0134      14341.A0135      14341.A0136 
    ##                1                1                1                1 
    ##      14341.A0137      14341.A0138      14341.A0139      14341.A0140 
    ##                1                1                1                1 
    ##      14341.A0141      14341.A0142      14341.A0143      14341.A0144 
    ##                1                1                1                1 
    ##      14341.A0145      14341.A0146      14341.A0147      14341.A0148 
    ##                1                1                1                1 
    ##      14341.A0149      14341.A0150      14341.A0151      14341.A0152 
    ##                1                1                1                1 
    ##      14341.A0153      14341.A0154      14341.A0155      14341.A0156 
    ##                1                1                1                1 
    ##      14341.A0157      14341.A0158      14341.A0159      14341.A0160 
    ##                1                1                1                1 
    ##      14341.A0161      14341.A0162      14341.A0163      14341.A0164 
    ##                1                1                1                1 
    ##      14341.A0165      14341.A0166      14341.A0167      14341.A0168 
    ##                1                1                1                1 
    ##      14341.A0169      14341.A0170      14341.A0171      14341.A0172 
    ##                1                1                1                1 
    ##      14341.A0173      14341.A0174      14341.A0175      14341.A0176 
    ##                1                1                1                1 
    ##      14341.A0177      14341.A0178      14341.A0179      14341.A0180 
    ##                1                1                1                1 
    ##      14341.A0181      14341.A0182      14341.A0183      14341.A0184 
    ##                1                1                1                1 
    ##      14341.A0185      14341.A0186      14341.A0187      14341.A0188 
    ##                1                1                1                1 
    ##      14341.A0189      14341.A0190      14341.A0191      14341.A0192 
    ##                1                1                1                1 
    ##      14341.A0193      14341.A0194      14341.A0195      14341.A0196 
    ##                1                1                1                1 
    ##      14341.A0197      14341.A0198      14341.A0199      14341.A0200 
    ##                1                1                1                1 
    ##      14341.A0201      14341.A0202      14341.A0203      14341.A0204 
    ##                1                1                1                1 
    ##      14341.A0205      14341.A0206      14341.A0207      14341.A0208 
    ##                1                1                1                1 
    ##      14341.A0209      14341.A0210      14341.A0211      14341.A0212 
    ##                1                1                1                1 
    ##      14341.A0213      14341.A0214      14341.A0215      14341.A0216 
    ##                1                1                1                1 
    ##      14341.A0217      14341.A0218      14341.A0219      14341.A0220 
    ##                1                1                1                1 
    ##      14341.A0221      14341.A0222      14341.A0223      14341.A0224 
    ##                1                1                1                1 
    ##      14341.A0225      14341.A0226      14341.A0227      14341.A0228 
    ##                1                1                1                1 
    ##      14341.A0229      14341.A0230      14341.A0231      14341.A0232 
    ##                1                1                1                1 
    ##      14341.A0233      14341.A0234      14341.A0235      14341.A0236 
    ##                1                1                1                1 
    ##      14341.A0237      14341.A0238      14341.A0239      14341.A0240 
    ##                1                1                1                1 
    ##      14341.A0241      14341.A0242      14341.A0243      14341.A0244 
    ##                1                1                1                1 
    ##      14341.A0245      14341.A0246      14341.A0247      14341.A0248 
    ##                1                1                1                1 
    ##      14341.A0249      14341.A0250      14341.A0251      14341.A0252 
    ##                1                1                1                1 
    ##      14341.A0253      14341.A0254      14341.A0255      14341.A0256 
    ##                1                1                1                1 
    ##      14341.A0257      14341.A0258      14341.A0259      14341.A0260 
    ##                1                1                1                1 
    ##      14341.A0261      14341.A0262      14341.A0263      14341.A0264 
    ##                1                1                1                1 
    ##      14341.A0265      14341.A0266      14341.A0267      14341.A0268 
    ##                1                1                1                1 
    ##      14341.A0269      14341.A0270      14341.A0271      14341.A0272 
    ##                1                1                1                1 
    ##      14341.A0273      14341.A0274      14341.A0275      14341.A0276 
    ##                1                1                1                1 
    ##      14341.A0277      14341.A0278      14341.A0279      14341.A0280 
    ##                1                1                1                1 
    ##      14341.A0281      14341.A0282      14341.A0283      14341.A0284 
    ##                1                1                1                1 
    ##      14341.A0285      14341.A0286      14341.A0287      14341.A0288 
    ##                1                1                1                1 
    ##      14341.A0289      14341.A0290      14341.A0291      14341.A0292 
    ##                1                1                1                1 
    ##      14341.A0293      14341.A0294      14341.A0295      14341.A0296 
    ##                1                1                1                1 
    ##      14341.A0297      14341.A0298      14341.A0299      14341.A0300 
    ##                1                1                1                1 
    ##      14341.A0301      14341.A0302      14341.A0303      14341.A0304 
    ##                1                1                1                1 
    ##      14341.A0305      14341.A0306      14341.A0307      14341.A0308 
    ##                1                1                1                1 
    ##      14341.A0309      14341.A0310      14341.A0311      14341.A0312 
    ##                1                1                1                1 
    ##      14341.A0313      14341.A0314      14341.A0315      14341.A0316 
    ##                1                1                1                1 
    ##      14341.A0317      14341.A0318      14341.A0319      14341.A0320 
    ##                1                1                1                1 
    ##      14341.A0321      14341.A0322      14341.A0323      14341.A0324 
    ##                1                1                1                1 
    ##      14341.A0325      14341.A0326      14341.A0327      14341.A0328 
    ##                1                1                1                1 
    ##      14341.A0329      14341.A0330      14341.A0331      14341.A0332 
    ##                1                1                1                1 
    ##      14341.A0333      14341.A0334      14341.A0335      14341.A0336 
    ##                1                1                1                1 
    ##      14341.A0337      14341.A0338      14341.A0339      14341.A0340 
    ##                1                1                1                1 
    ##      14341.A0341      14341.A0342      14341.A0343      14341.A0344 
    ##                1                1                1                1 
    ##      14341.A0345      14341.A0346      14341.A0347      14341.A0348 
    ##                1                1                1                1 
    ##      14341.A0349      14341.A0350      14341.A0351      14341.A0352 
    ##                1                1                1                1 
    ##      14341.A0353      14341.A0354      14341.A0355      14341.A0356 
    ##                1                1                1                1 
    ##      14341.A0357      14341.A0358      14341.A0359      14341.A0360 
    ##                1                1                1                1 
    ##      14341.A0361      14341.A0362      14341.A0363      14341.A0364 
    ##                1                1                1                1 
    ##      14341.A0365      14341.A0366      14341.A0367      14341.A0368 
    ##                1                1                1                1 
    ##      14341.A0369      14341.A0370      14341.A0371      14341.A0372 
    ##                1                1                1                1 
    ##      14341.A0373      14341.A0374      14341.A0375      14341.A0376 
    ##                1                1                1                1 
    ##      14341.A0377      14341.A0378      14341.A0379      14341.A0380 
    ##                1                1                1                1 
    ##      14341.A0381      14341.A0382      14341.A0383      14341.A0384 
    ##                1                1                1                1 
    ##      14341.A0385      14341.A0386      14341.A0387      14341.A0388 
    ##                1                1                1                1 
    ##      14341.A0389      14341.A0390      14341.A0391      14341.A0392 
    ##                1                1                1                1 
    ##      14341.A0393      14341.A0394      14341.A0395      14341.A0396 
    ##                1                1                1                1 
    ##      14341.A0397      14341.A0398      14341.A0399      14341.A0400 
    ##                1                1                1                1 
    ##      14341.A0401      14341.A0402      14341.A0403      14341.A0404 
    ##                1                1                1                1 
    ##      14341.A0405      14341.A0406      14341.A0407      14341.A0408 
    ##                1                1                1                1 
    ##      14341.A0409      14341.A0410      14341.A0411      14341.A0412 
    ##                1                1                1                1 
    ##      14341.A0413      14341.A0414      14341.A0415      14341.A0416 
    ##                1                1                1                1 
    ##      14341.A0417      14341.A0418      14341.A0419      14341.A0420 
    ##                1                1                1                1 
    ##      14341.A0421      14341.A0422      14341.A0423      14341.A0424 
    ##                1                1                1                1 
    ##      14341.A0425      14341.A0426      14341.A0427      14341.A0428 
    ##                1                1                1                1 
    ##      14341.A0429      14341.A0430      14341.A0431      14341.A0432 
    ##                1                1                1                1 
    ##      14341.A0433      14341.A0434      14341.A0435      14341.A0436 
    ##                1                1                1                1 
    ##      14341.A0437      14341.A0438      14341.A0439      14341.A0440 
    ##                1                1                1                1 
    ##      14341.A0441      14341.A0442      14341.A0443      14341.A0444 
    ##                1                1                1                1 
    ##      14341.A0445      14341.A0446      14341.A0447      14341.A0448 
    ##                1                1                1                1 
    ##      14341.A0449      14341.A0450      14341.A0451      14341.A0452 
    ##                1                1                1                1 
    ##      14341.A0453      14341.A0454      14341.A0455      14341.A0456 
    ##                1                1                1                1 
    ##      14341.A0457      14341.A0458      14341.A0459      14341.A0460 
    ##                1                1                1                1 
    ##      14341.A0461      14341.A0462      14341.A0463      14341.A0464 
    ##                1                1                1                1 
    ##      14341.A0465      14341.A0466      14341.A0467      14341.A0468 
    ##                1                1                1                1 
    ##      14341.A0469      14341.A0470      14341.A0471      14341.A0472 
    ##                1                1                1                1 
    ##      14341.A0473      14341.A0474      14341.A0475      14341.A0476 
    ##                1                1                1                1 
    ##      14341.A0477      14341.A0478      14341.A0479      14341.A0480 
    ##                1                1                1                1 
    ##      14341.A0481      14341.A0482      14341.A0483      14341.A0484 
    ##                1                1                1                1 
    ##      14341.A0485      14341.A0486      14341.A0487      14341.A0488 
    ##                1                1                1                1 
    ##      14341.A0489      14341.A0490      14341.A0491      14341.A0492 
    ##                1                1                1                1 
    ##      14341.A0493      14341.A0494      14341.A0495      14341.A0496 
    ##                1                1                1                1 
    ##      14341.A0497      14341.A0498      14341.A0499      14341.A0500 
    ##                1                1                1                1 
    ##      14341.A0501      14341.A0502      14341.A0503      14341.A0504 
    ##                1                1                1                1 
    ##      14341.A0505      14341.A0506      14341.A0507      14341.A0508 
    ##                1                1                1                1 
    ##      14341.A0509      14341.A0510      14341.A0511      14341.A0512 
    ##                1                1                1                1 
    ##      14341.A0513      14341.A0514      14341.A0515      14341.A0516 
    ##                1                1                1                1 
    ##      14341.A0517      14341.A0518      14341.A0519      14341.A0520 
    ##                1                1                1                1 
    ##      14341.A0521      14341.A0522      14341.A0523      14341.A0524 
    ##                1                1                1                1 
    ##      14341.A0525      14341.A0526      14341.A0527      14341.A0528 
    ##                1                1                1                1 
    ##      14341.A0529      14341.A0530      14341.A0531      14341.A0532 
    ##                1                1                1                1 
    ##      14341.A0533      14341.A0534      14341.A0535      14341.A0536 
    ##                1                1                1                1 
    ##      14341.A0537      14341.A0538      14341.A0539      14341.A0540 
    ##                1                1                1                1 
    ##      14341.A0541      14341.A0542      14341.A0543      14341.A0544 
    ##                1                1                1                1 
    ##      14341.A0545      14341.A0546      14341.A0547      14341.A0548 
    ##                1                1                1                1 
    ##      14341.A0549      14341.A0550      14341.A0551      14341.A0552 
    ##                1                1                1                1 
    ##      14341.A0553      14341.A0554      14341.A0555      14341.A0556 
    ##                1                1                1                1 
    ##      14341.A0557      14341.A0558      14341.A0559      14341.A0560 
    ##                1                1                1                1 
    ##      14341.A0561      14341.A0562      14341.A0563      14341.A0564 
    ##                1                1                1                1 
    ##      14341.A0565      14341.A0566      14341.A0567      14341.A0568 
    ##                1                1                1                1 
    ##      14341.A0569      14341.A0570      14341.A0571      14341.A0572 
    ##                1                1                1                1 
    ##      14341.A0573      14341.A0574      14341.A0575      14341.A0576 
    ##                1                1                1                1 
    ##      14341.A0577      14341.A0578      14341.A0579      14341.A0580 
    ##                1                1                1                1 
    ##      14341.A0581      14341.A0582      14341.A0583      14341.A0584 
    ##                1                1                1                1 
    ##      14341.A0585      14341.A0586      14341.A0587      14341.A0588 
    ##                1                1                1                1 
    ##      14341.A0589      14341.A0590      14341.A0591      14341.A0592 
    ##                1                1                1                1 
    ##      14341.A0593      14341.A0594      14341.A0595      14341.A0596 
    ##                1                1                1                1 
    ##      14341.A0597      14341.A0598      14341.A0599      14341.A0600 
    ##                1                1                1                1 
    ##      14341.A0601      14341.A0602      14341.A0603      14341.A0604 
    ##                1                1                1                1 
    ##      14341.A0605      14341.A0606      14341.A0607      14341.A0608 
    ##                1                1                1                1 
    ##      14341.A0609      14341.A0610      14341.A0611      14341.A0612 
    ##                1                1                1                1 
    ##      14341.A0613      14341.A0614      14341.A0615      14341.A0616 
    ##                1                1                1                1 
    ##      14341.A0617      14341.A0618      14341.A0619      14341.A0620 
    ##                1                1                1                1 
    ##      14341.A0621      14341.A0622      14341.A0623      14341.A0624 
    ##                1                1                1                1 
    ##      14341.A0625      14341.A0626      14341.A0627      14341.A0628 
    ##                1                1                1                1 
    ##      14341.A0629      14341.A0630      14341.A0631      14341.A0632 
    ##                1                1                1                1 
    ##      14341.A0633      14341.A0634      14341.A0635      14341.A0636 
    ##                1                1                1                1 
    ##      14341.A0637         14341.S1        14341.S10        14341.S11 
    ##                1                1                1                1 
    ##        14341.S12        14341.S13        14341.S14        14341.S15 
    ##                1                1                1                1 
    ##        14341.S16        14341.S17        14341.S18        14341.S19 
    ##                1                1                1                1 
    ##         14341.S2        14341.S20        14341.S21        14341.S22 
    ##                1                1                1                1 
    ##        14341.S23        14341.S24         14341.S3         14341.S4 
    ##                1                1                1                1 
    ##         14341.S5         14341.S6         14341.S7         14341.S8 
    ##                1                1                1                1 
    ##         14341.S9    14341.SAMPLE1   14341.SAMPLE10  14341.SAMPLE100 
    ##                1                1                1                1 
    ## 14341.SAMPLE1000 14341.SAMPLE1001 14341.SAMPLE1002 14341.SAMPLE1003 
    ##                1                1                1                1 
    ## 14341.SAMPLE1004 14341.SAMPLE1005 14341.SAMPLE1006 14341.SAMPLE1007 
    ##                1                1                1                1 
    ## 14341.SAMPLE1008 14341.SAMPLE1009  14341.SAMPLE101 14341.SAMPLE1010 
    ##                1                1                1                1 
    ## 14341.SAMPLE1011 14341.SAMPLE1012 14341.SAMPLE1013 14341.SAMPLE1014 
    ##                1                1                1                1 
    ## 14341.SAMPLE1015 14341.SAMPLE1016 14341.SAMPLE1017 14341.SAMPLE1018 
    ##                1                1                1                1 
    ## 14341.SAMPLE1019  14341.SAMPLE102 14341.SAMPLE1020 14341.SAMPLE1021 
    ##                1                1                1                1 
    ## 14341.SAMPLE1022 14341.SAMPLE1023 14341.SAMPLE1024 14341.SAMPLE1025 
    ##                1                1                1                1 
    ## 14341.SAMPLE1026 14341.SAMPLE1027 14341.SAMPLE1028 14341.SAMPLE1029 
    ##                1                1                1                1 
    ##  14341.SAMPLE103 14341.SAMPLE1030 14341.SAMPLE1031 14341.SAMPLE1032 
    ##                1                1                1                1 
    ## 14341.SAMPLE1033 14341.SAMPLE1034 14341.SAMPLE1035 14341.SAMPLE1036 
    ##                1                1                1                1 
    ## 14341.SAMPLE1037 14341.SAMPLE1038 14341.SAMPLE1039  14341.SAMPLE104 
    ##                1                1                1                1 
    ## 14341.SAMPLE1040 14341.SAMPLE1041 14341.SAMPLE1042 14341.SAMPLE1043 
    ##                1                1                1                1 
    ## 14341.SAMPLE1044 14341.SAMPLE1045 14341.SAMPLE1046 14341.SAMPLE1047 
    ##                1                1                1                1 
    ## 14341.SAMPLE1048 14341.SAMPLE1049  14341.SAMPLE105 14341.SAMPLE1050 
    ##                1                1                1                1 
    ## 14341.SAMPLE1051 14341.SAMPLE1052 14341.SAMPLE1053 14341.SAMPLE1054 
    ##                1                1                1                1 
    ## 14341.SAMPLE1055 14341.SAMPLE1056 14341.SAMPLE1057 14341.SAMPLE1058 
    ##                1                1                1                1 
    ## 14341.SAMPLE1059  14341.SAMPLE106 14341.SAMPLE1060 14341.SAMPLE1061 
    ##                1                1                1                1 
    ## 14341.SAMPLE1062 14341.SAMPLE1063 14341.SAMPLE1064 14341.SAMPLE1065 
    ##                1                1                1                1 
    ## 14341.SAMPLE1066 14341.SAMPLE1067 14341.SAMPLE1068 14341.SAMPLE1069 
    ##                1                1                1                1 
    ##  14341.SAMPLE107 14341.SAMPLE1070 14341.SAMPLE1071 14341.SAMPLE1072 
    ##                1                1                1                1 
    ## 14341.SAMPLE1073 14341.SAMPLE1074 14341.SAMPLE1075 14341.SAMPLE1076 
    ##                1                1                1                1 
    ## 14341.SAMPLE1077 14341.SAMPLE1078 14341.SAMPLE1079  14341.SAMPLE108 
    ##                1                1                1                1 
    ## 14341.SAMPLE1080 14341.SAMPLE1081 14341.SAMPLE1082 14341.SAMPLE1083 
    ##                1                1                1                1 
    ## 14341.SAMPLE1084 14341.SAMPLE1085 14341.SAMPLE1086 14341.SAMPLE1087 
    ##                1                1                1                1 
    ## 14341.SAMPLE1088 14341.SAMPLE1089  14341.SAMPLE109 14341.SAMPLE1090 
    ##                1                1                1                1 
    ## 14341.SAMPLE1091 14341.SAMPLE1092 14341.SAMPLE1093 14341.SAMPLE1094 
    ##                1                1                1                1 
    ## 14341.SAMPLE1095 14341.SAMPLE1096 14341.SAMPLE1097 14341.SAMPLE1098 
    ##                1                1                1                1 
    ## 14341.SAMPLE1099   14341.SAMPLE11  14341.SAMPLE110 14341.SAMPLE1100 
    ##                1                1                1                1 
    ## 14341.SAMPLE1101 14341.SAMPLE1102 14341.SAMPLE1103 14341.SAMPLE1104 
    ##                1                1                1                1 
    ## 14341.SAMPLE1105 14341.SAMPLE1106 14341.SAMPLE1107 14341.SAMPLE1108 
    ##                1                1                1                1 
    ## 14341.SAMPLE1109  14341.SAMPLE111 14341.SAMPLE1110 14341.SAMPLE1111 
    ##                1                1                1                1 
    ## 14341.SAMPLE1112 14341.SAMPLE1113 14341.SAMPLE1114 14341.SAMPLE1115 
    ##                1                1                1                1 
    ## 14341.SAMPLE1116 14341.SAMPLE1117 14341.SAMPLE1118 14341.SAMPLE1119 
    ##                1                1                1                1 
    ##  14341.SAMPLE112 14341.SAMPLE1120 14341.SAMPLE1121 14341.SAMPLE1122 
    ##                1                1                1                1 
    ## 14341.SAMPLE1123 14341.SAMPLE1124 14341.SAMPLE1125 14341.SAMPLE1126 
    ##                1                1                1                1 
    ## 14341.SAMPLE1127 14341.SAMPLE1128 14341.SAMPLE1129  14341.SAMPLE113 
    ##                1                1                1                1 
    ## 14341.SAMPLE1130 14341.SAMPLE1131 14341.SAMPLE1132 14341.SAMPLE1133 
    ##                1                1                1                1 
    ## 14341.SAMPLE1134 14341.SAMPLE1135 14341.SAMPLE1136 14341.SAMPLE1137 
    ##                1                1                1                1 
    ## 14341.SAMPLE1138 14341.SAMPLE1139  14341.SAMPLE114 14341.SAMPLE1140 
    ##                1                1                1                1 
    ## 14341.SAMPLE1141 14341.SAMPLE1142 14341.SAMPLE1143 14341.SAMPLE1144 
    ##                1                1                1                1 
    ## 14341.SAMPLE1145 14341.SAMPLE1146 14341.SAMPLE1147 14341.SAMPLE1148 
    ##                1                1                1                1 
    ## 14341.SAMPLE1149  14341.SAMPLE115 14341.SAMPLE1150 14341.SAMPLE1151 
    ##                1                1                1                1 
    ## 14341.SAMPLE1152 14341.SAMPLE1153 14341.SAMPLE1154 14341.SAMPLE1155 
    ##                1                1                1                1 
    ## 14341.SAMPLE1156 14341.SAMPLE1157 14341.SAMPLE1158 14341.SAMPLE1159 
    ##                1                1                1                1 
    ##  14341.SAMPLE116 14341.SAMPLE1160 14341.SAMPLE1161 14341.SAMPLE1162 
    ##                1                1                1                1 
    ## 14341.SAMPLE1163 14341.SAMPLE1164 14341.SAMPLE1165 14341.SAMPLE1166 
    ##                1                1                1                1 
    ## 14341.SAMPLE1167 14341.SAMPLE1168 14341.SAMPLE1169  14341.SAMPLE117 
    ##                1                1                1                1 
    ## 14341.SAMPLE1170 14341.SAMPLE1171 14341.SAMPLE1172 14341.SAMPLE1173 
    ##                1                1                1                1 
    ## 14341.SAMPLE1174 14341.SAMPLE1175 14341.SAMPLE1176 14341.SAMPLE1177 
    ##                1                1                1                1 
    ## 14341.SAMPLE1178 14341.SAMPLE1179  14341.SAMPLE118 14341.SAMPLE1180 
    ##                1                1                1                1 
    ## 14341.SAMPLE1181 14341.SAMPLE1182 14341.SAMPLE1183 14341.SAMPLE1184 
    ##                1                1                1                1 
    ## 14341.SAMPLE1185 14341.SAMPLE1186 14341.SAMPLE1187 14341.SAMPLE1188 
    ##                1                1                1                1 
    ## 14341.SAMPLE1189  14341.SAMPLE119 14341.SAMPLE1190 14341.SAMPLE1191 
    ##                1                1                1                1 
    ## 14341.SAMPLE1192 14341.SAMPLE1193 14341.SAMPLE1194 14341.SAMPLE1195 
    ##                1                1                1                1 
    ## 14341.SAMPLE1196 14341.SAMPLE1197 14341.SAMPLE1198 14341.SAMPLE1199 
    ##                1                1                1                1 
    ##   14341.SAMPLE12  14341.SAMPLE120 14341.SAMPLE1200 14341.SAMPLE1201 
    ##                1                1                1                1 
    ## 14341.SAMPLE1202 14341.SAMPLE1203 14341.SAMPLE1204 14341.SAMPLE1205 
    ##                1                1                1                1 
    ## 14341.SAMPLE1206 14341.SAMPLE1207 14341.SAMPLE1208 14341.SAMPLE1209 
    ##                1                1                1                1 
    ##  14341.SAMPLE121 14341.SAMPLE1210 14341.SAMPLE1211 14341.SAMPLE1212 
    ##                1                1                1                1 
    ## 14341.SAMPLE1213 14341.SAMPLE1214 14341.SAMPLE1215 14341.SAMPLE1216 
    ##                1                1                1                1 
    ## 14341.SAMPLE1217 14341.SAMPLE1218 14341.SAMPLE1219  14341.SAMPLE122 
    ##                1                1                1                1 
    ## 14341.SAMPLE1220 14341.SAMPLE1221 14341.SAMPLE1222 14341.SAMPLE1223 
    ##                1                1                1                1 
    ## 14341.SAMPLE1224 14341.SAMPLE1225 14341.SAMPLE1226 14341.SAMPLE1227 
    ##                1                1                1                1 
    ## 14341.SAMPLE1228 14341.SAMPLE1229  14341.SAMPLE123 14341.SAMPLE1230 
    ##                1                1                1                1 
    ## 14341.SAMPLE1231 14341.SAMPLE1232 14341.SAMPLE1233 14341.SAMPLE1234 
    ##                1                1                1                1 
    ## 14341.SAMPLE1235 14341.SAMPLE1236 14341.SAMPLE1237 14341.SAMPLE1238 
    ##                1                1                1                1 
    ## 14341.SAMPLE1239  14341.SAMPLE124 14341.SAMPLE1240 14341.SAMPLE1241 
    ##                1                1                1                1 
    ## 14341.SAMPLE1242 14341.SAMPLE1243 14341.SAMPLE1244 14341.SAMPLE1245 
    ##                1                1                1                1 
    ## 14341.SAMPLE1246 14341.SAMPLE1247 14341.SAMPLE1248 14341.SAMPLE1249 
    ##                1                1                1                1 
    ##  14341.SAMPLE125 14341.SAMPLE1250 14341.SAMPLE1251 14341.SAMPLE1252 
    ##                1                1                1                1 
    ## 14341.SAMPLE1253 14341.SAMPLE1254 14341.SAMPLE1255 14341.SAMPLE1256 
    ##                1                1                1                1 
    ## 14341.SAMPLE1257 14341.SAMPLE1258 14341.SAMPLE1259  14341.SAMPLE126 
    ##                1                1                1                1 
    ## 14341.SAMPLE1260 14341.SAMPLE1261 14341.SAMPLE1262 14341.SAMPLE1263 
    ##                1                1                1                1 
    ## 14341.SAMPLE1264 14341.SAMPLE1265 14341.SAMPLE1266 14341.SAMPLE1267 
    ##                1                1                1                1 
    ## 14341.SAMPLE1268 14341.SAMPLE1269  14341.SAMPLE127 14341.SAMPLE1270 
    ##                1                1                1                1 
    ## 14341.SAMPLE1271 14341.SAMPLE1272 14341.SAMPLE1273 14341.SAMPLE1274 
    ##                1                1                1                1 
    ## 14341.SAMPLE1275 14341.SAMPLE1276 14341.SAMPLE1277 14341.SAMPLE1278 
    ##                1                1                1                1 
    ## 14341.SAMPLE1279  14341.SAMPLE128 14341.SAMPLE1280 14341.SAMPLE1281 
    ##                1                1                1                1 
    ## 14341.SAMPLE1282 14341.SAMPLE1283 14341.SAMPLE1284 14341.SAMPLE1285 
    ##                1                1                1                1 
    ## 14341.SAMPLE1286 14341.SAMPLE1287 14341.SAMPLE1288 14341.SAMPLE1289 
    ##                1                1                1                1 
    ##  14341.SAMPLE129 14341.SAMPLE1290 14341.SAMPLE1291 14341.SAMPLE1292 
    ##                1                1                1                1 
    ## 14341.SAMPLE1293 14341.SAMPLE1294 14341.SAMPLE1295 14341.SAMPLE1296 
    ##                1                1                1                1 
    ## 14341.SAMPLE1297 14341.SAMPLE1298 14341.SAMPLE1299   14341.SAMPLE13 
    ##                1                1                1                1 
    ##  14341.SAMPLE130 14341.SAMPLE1300 14341.SAMPLE1301 14341.SAMPLE1302 
    ##                1                1                1                1 
    ## 14341.SAMPLE1303 14341.SAMPLE1304 14341.SAMPLE1305 14341.SAMPLE1306 
    ##                1                1                1                1 
    ## 14341.SAMPLE1307 14341.SAMPLE1308 14341.SAMPLE1309  14341.SAMPLE131 
    ##                1                1                1                1 
    ## 14341.SAMPLE1310 14341.SAMPLE1311 14341.SAMPLE1312 14341.SAMPLE1313 
    ##                1                1                1                1 
    ## 14341.SAMPLE1314 14341.SAMPLE1315 14341.SAMPLE1316 14341.SAMPLE1317 
    ##                1                1                1                1 
    ## 14341.SAMPLE1318 14341.SAMPLE1319  14341.SAMPLE132 14341.SAMPLE1320 
    ##                1                1                1                1 
    ## 14341.SAMPLE1321 14341.SAMPLE1322 14341.SAMPLE1323 14341.SAMPLE1324 
    ##                1                1                1                1 
    ## 14341.SAMPLE1325 14341.SAMPLE1326 14341.SAMPLE1327 14341.SAMPLE1328 
    ##                1                1                1                1 
    ## 14341.SAMPLE1329  14341.SAMPLE133 14341.SAMPLE1330 14341.SAMPLE1331 
    ##                1                1                1                1 
    ## 14341.SAMPLE1332 14341.SAMPLE1333 14341.SAMPLE1334 14341.SAMPLE1335 
    ##                1                1                1                1 
    ## 14341.SAMPLE1336 14341.SAMPLE1337 14341.SAMPLE1338 14341.SAMPLE1339 
    ##                1                1                1                1 
    ##  14341.SAMPLE134 14341.SAMPLE1340 14341.SAMPLE1341 14341.SAMPLE1342 
    ##                1                1                1                1 
    ## 14341.SAMPLE1343 14341.SAMPLE1344 14341.SAMPLE1345 14341.SAMPLE1346 
    ##                1                1                1                1 
    ## 14341.SAMPLE1347 14341.SAMPLE1348 14341.SAMPLE1349  14341.SAMPLE135 
    ##                1                1                1                1 
    ## 14341.SAMPLE1350 14341.SAMPLE1351 14341.SAMPLE1352 14341.SAMPLE1353 
    ##                1                1                1                1 
    ## 14341.SAMPLE1354 14341.SAMPLE1355 14341.SAMPLE1356 14341.SAMPLE1357 
    ##                1                1                1                1 
    ## 14341.SAMPLE1358 14341.SAMPLE1359  14341.SAMPLE136 14341.SAMPLE1360 
    ##                1                1                1                1 
    ## 14341.SAMPLE1361 14341.SAMPLE1362 14341.SAMPLE1363 14341.SAMPLE1364 
    ##                1                1                1                1 
    ## 14341.SAMPLE1365 14341.SAMPLE1366 14341.SAMPLE1367 14341.SAMPLE1368 
    ##                1                1                1                1 
    ## 14341.SAMPLE1369  14341.SAMPLE137 14341.SAMPLE1370 14341.SAMPLE1371 
    ##                1                1                1                1 
    ## 14341.SAMPLE1372 14341.SAMPLE1373 14341.SAMPLE1374 14341.SAMPLE1375 
    ##                1                1                1                1 
    ## 14341.SAMPLE1376 14341.SAMPLE1377 14341.SAMPLE1378 14341.SAMPLE1379 
    ##                1                1                1                1 
    ##  14341.SAMPLE138 14341.SAMPLE1380 14341.SAMPLE1381 14341.SAMPLE1382 
    ##                1                1                1                1 
    ## 14341.SAMPLE1383 14341.SAMPLE1384 14341.SAMPLE1385 14341.SAMPLE1386 
    ##                1                1                1                1 
    ## 14341.SAMPLE1387 14341.SAMPLE1388 14341.SAMPLE1389  14341.SAMPLE139 
    ##                1                1                1                1 
    ## 14341.SAMPLE1390 14341.SAMPLE1391 14341.SAMPLE1392 14341.SAMPLE1393 
    ##                1                1                1                1 
    ## 14341.SAMPLE1394 14341.SAMPLE1395 14341.SAMPLE1396 14341.SAMPLE1397 
    ##                1                1                1                1 
    ## 14341.SAMPLE1398 14341.SAMPLE1399   14341.SAMPLE14  14341.SAMPLE140 
    ##                1                1                1                1 
    ## 14341.SAMPLE1400 14341.SAMPLE1401 14341.SAMPLE1402 14341.SAMPLE1403 
    ##                1                1                1                1 
    ## 14341.SAMPLE1404 14341.SAMPLE1405 14341.SAMPLE1406 14341.SAMPLE1407 
    ##                1                1                1                1 
    ## 14341.SAMPLE1408 14341.SAMPLE1409  14341.SAMPLE141 14341.SAMPLE1410 
    ##                1                1                1                1 
    ## 14341.SAMPLE1411 14341.SAMPLE1412 14341.SAMPLE1413 14341.SAMPLE1414 
    ##                1                1                1                1 
    ## 14341.SAMPLE1415 14341.SAMPLE1416 14341.SAMPLE1417 14341.SAMPLE1418 
    ##                1                1                1                1 
    ## 14341.SAMPLE1419  14341.SAMPLE142 14341.SAMPLE1420 14341.SAMPLE1421 
    ##                1                1                1                1 
    ## 14341.SAMPLE1422 14341.SAMPLE1423 14341.SAMPLE1424 14341.SAMPLE1425 
    ##                1                1                1                1 
    ## 14341.SAMPLE1426 14341.SAMPLE1427 14341.SAMPLE1428 14341.SAMPLE1429 
    ##                1                1                1                1 
    ##  14341.SAMPLE143 14341.SAMPLE1430 14341.SAMPLE1431 14341.SAMPLE1432 
    ##                1                1                1                1 
    ## 14341.SAMPLE1433 14341.SAMPLE1434 14341.SAMPLE1435 14341.SAMPLE1436 
    ##                1                1                1                1 
    ## 14341.SAMPLE1437 14341.SAMPLE1438 14341.SAMPLE1439  14341.SAMPLE144 
    ##                1                1                1                1 
    ## 14341.SAMPLE1440 14341.SAMPLE1441 14341.SAMPLE1442 14341.SAMPLE1443 
    ##                1                1                1                1 
    ## 14341.SAMPLE1444 14341.SAMPLE1445 14341.SAMPLE1446 14341.SAMPLE1447 
    ##                1                1                1                1 
    ## 14341.SAMPLE1448 14341.SAMPLE1449  14341.SAMPLE145 14341.SAMPLE1450 
    ##                1                1                1                1 
    ## 14341.SAMPLE1451 14341.SAMPLE1452 14341.SAMPLE1453 14341.SAMPLE1454 
    ##                1                1                1                1 
    ## 14341.SAMPLE1455 14341.SAMPLE1456 14341.SAMPLE1457 14341.SAMPLE1458 
    ##                1                1                1                1 
    ## 14341.SAMPLE1459  14341.SAMPLE146 14341.SAMPLE1460 14341.SAMPLE1461 
    ##                1                1                1                1 
    ## 14341.SAMPLE1462 14341.SAMPLE1463 14341.SAMPLE1464 14341.SAMPLE1465 
    ##                1                1                1                1 
    ## 14341.SAMPLE1466 14341.SAMPLE1467 14341.SAMPLE1468 14341.SAMPLE1469 
    ##                1                1                1                1 
    ##  14341.SAMPLE147 14341.SAMPLE1470 14341.SAMPLE1471 14341.SAMPLE1472 
    ##                1                1                1                1 
    ## 14341.SAMPLE1473 14341.SAMPLE1474 14341.SAMPLE1475 14341.SAMPLE1476 
    ##                1                1                1                1 
    ## 14341.SAMPLE1477 14341.SAMPLE1478 14341.SAMPLE1479  14341.SAMPLE148 
    ##                1                1                1                1 
    ## 14341.SAMPLE1480 14341.SAMPLE1481 14341.SAMPLE1482 14341.SAMPLE1483 
    ##                1                1                1                1 
    ## 14341.SAMPLE1484 14341.SAMPLE1485 14341.SAMPLE1486 14341.SAMPLE1487 
    ##                1                1                1                1 
    ## 14341.SAMPLE1488 14341.SAMPLE1489  14341.SAMPLE149 14341.SAMPLE1490 
    ##                1                1                1                1 
    ## 14341.SAMPLE1491 14341.SAMPLE1492 14341.SAMPLE1493 14341.SAMPLE1494 
    ##                1                1                1                1 
    ## 14341.SAMPLE1495 14341.SAMPLE1496 14341.SAMPLE1497 14341.SAMPLE1498 
    ##                1                1                1                1 
    ## 14341.SAMPLE1499   14341.SAMPLE15  14341.SAMPLE150 14341.SAMPLE1500 
    ##                1                1                1                1 
    ## 14341.SAMPLE1501 14341.SAMPLE1502 14341.SAMPLE1503 14341.SAMPLE1504 
    ##                1                1                1                1 
    ## 14341.SAMPLE1505 14341.SAMPLE1506 14341.SAMPLE1507 14341.SAMPLE1508 
    ##                1                1                1                1 
    ## 14341.SAMPLE1509  14341.SAMPLE151 14341.SAMPLE1510 14341.SAMPLE1511 
    ##                1                1                1                1 
    ## 14341.SAMPLE1512 14341.SAMPLE1513 14341.SAMPLE1514 14341.SAMPLE1515 
    ##                1                1                1                1 
    ## 14341.SAMPLE1516 14341.SAMPLE1517 14341.SAMPLE1518 14341.SAMPLE1519 
    ##                1                1                1                1 
    ##  14341.SAMPLE152 14341.SAMPLE1520 14341.SAMPLE1521 14341.SAMPLE1522 
    ##                1                1                1                1 
    ## 14341.SAMPLE1523 14341.SAMPLE1524 14341.SAMPLE1525 14341.SAMPLE1526 
    ##                1                1                1                1 
    ## 14341.SAMPLE1527 14341.SAMPLE1528 14341.SAMPLE1529  14341.SAMPLE153 
    ##                1                1                1                1 
    ## 14341.SAMPLE1530 14341.SAMPLE1531 14341.SAMPLE1532 14341.SAMPLE1533 
    ##                1                1                1                1 
    ## 14341.SAMPLE1534 14341.SAMPLE1535 14341.SAMPLE1536 14341.SAMPLE1537 
    ##                1                1                1                1 
    ## 14341.SAMPLE1538 14341.SAMPLE1539  14341.SAMPLE154 14341.SAMPLE1540 
    ##                1                1                1                1 
    ## 14341.SAMPLE1541 14341.SAMPLE1542 14341.SAMPLE1543 14341.SAMPLE1544 
    ##                1                1                1                1 
    ## 14341.SAMPLE1545 14341.SAMPLE1546 14341.SAMPLE1547 14341.SAMPLE1548 
    ##                1                1                1                1 
    ## 14341.SAMPLE1549  14341.SAMPLE155 14341.SAMPLE1550 14341.SAMPLE1551 
    ##                1                1                1                1 
    ## 14341.SAMPLE1552 14341.SAMPLE1553 14341.SAMPLE1554 14341.SAMPLE1555 
    ##                1                1                1                1 
    ## 14341.SAMPLE1556 14341.SAMPLE1557 14341.SAMPLE1558 14341.SAMPLE1559 
    ##                1                1                1                1 
    ##  14341.SAMPLE156 14341.SAMPLE1560 14341.SAMPLE1561 14341.SAMPLE1562 
    ##                1                1                1                1 
    ## 14341.SAMPLE1563 14341.SAMPLE1564 14341.SAMPLE1565 14341.SAMPLE1566 
    ##                1                1                1                1 
    ## 14341.SAMPLE1567 14341.SAMPLE1568 14341.SAMPLE1569  14341.SAMPLE157 
    ##                1                1                1                1 
    ## 14341.SAMPLE1570 14341.SAMPLE1571 14341.SAMPLE1572 14341.SAMPLE1573 
    ##                1                1                1                1 
    ## 14341.SAMPLE1574 14341.SAMPLE1575 14341.SAMPLE1576 14341.SAMPLE1577 
    ##                1                1                1                1 
    ## 14341.SAMPLE1578 14341.SAMPLE1579  14341.SAMPLE158 14341.SAMPLE1580 
    ##                1                1                1                1 
    ## 14341.SAMPLE1581 14341.SAMPLE1582 14341.SAMPLE1583 14341.SAMPLE1584 
    ##                1                1                1                1 
    ## 14341.SAMPLE1585 14341.SAMPLE1586 14341.SAMPLE1587 14341.SAMPLE1588 
    ##                1                1                1                1 
    ## 14341.SAMPLE1589  14341.SAMPLE159 14341.SAMPLE1590 14341.SAMPLE1591 
    ##                1                1                1                1 
    ## 14341.SAMPLE1592 14341.SAMPLE1593 14341.SAMPLE1594 14341.SAMPLE1595 
    ##                1                1                1                1 
    ## 14341.SAMPLE1596 14341.SAMPLE1597 14341.SAMPLE1598 14341.SAMPLE1599 
    ##                1                1                1                1 
    ##   14341.SAMPLE16  14341.SAMPLE160 14341.SAMPLE1600 14341.SAMPLE1601 
    ##                1                1                1                1 
    ## 14341.SAMPLE1602 14341.SAMPLE1603 14341.SAMPLE1604 14341.SAMPLE1605 
    ##                1                1                1                1 
    ## 14341.SAMPLE1606 14341.SAMPLE1607 14341.SAMPLE1608 14341.SAMPLE1609 
    ##                1                1                1                1 
    ##  14341.SAMPLE161 14341.SAMPLE1610 14341.SAMPLE1611 14341.SAMPLE1612 
    ##                1                1                1                1 
    ## 14341.SAMPLE1613 14341.SAMPLE1614 14341.SAMPLE1615 14341.SAMPLE1616 
    ##                1                1                1                1 
    ## 14341.SAMPLE1617 14341.SAMPLE1618 14341.SAMPLE1619  14341.SAMPLE162 
    ##                1                1                1                1 
    ## 14341.SAMPLE1620 14341.SAMPLE1621 14341.SAMPLE1622 14341.SAMPLE1623 
    ##                1                1                1                1 
    ## 14341.SAMPLE1624 14341.SAMPLE1625 14341.SAMPLE1626 14341.SAMPLE1627 
    ##                1                1                1                1 
    ## 14341.SAMPLE1628 14341.SAMPLE1629  14341.SAMPLE163 14341.SAMPLE1630 
    ##                1                1                1                1 
    ## 14341.SAMPLE1631 14341.SAMPLE1632 14341.SAMPLE1633 14341.SAMPLE1634 
    ##                1                1                1                1 
    ## 14341.SAMPLE1635 14341.SAMPLE1636 14341.SAMPLE1637 14341.SAMPLE1638 
    ##                1                1                1                1 
    ## 14341.SAMPLE1639  14341.SAMPLE164 14341.SAMPLE1640 14341.SAMPLE1641 
    ##                1                1                1                1 
    ## 14341.SAMPLE1642 14341.SAMPLE1643 14341.SAMPLE1644 14341.SAMPLE1645 
    ##                1                1                1                1 
    ## 14341.SAMPLE1646 14341.SAMPLE1647 14341.SAMPLE1648 14341.SAMPLE1649 
    ##                1                1                1                1 
    ##  14341.SAMPLE165 14341.SAMPLE1650 14341.SAMPLE1651 14341.SAMPLE1652 
    ##                1                1                1                1 
    ## 14341.SAMPLE1653 14341.SAMPLE1654 14341.SAMPLE1655 14341.SAMPLE1656 
    ##                1                1                1                1 
    ## 14341.SAMPLE1657 14341.SAMPLE1658 14341.SAMPLE1659  14341.SAMPLE166 
    ##                1                1                1                1 
    ## 14341.SAMPLE1660 14341.SAMPLE1661 14341.SAMPLE1662 14341.SAMPLE1663 
    ##                1                1                1                1 
    ## 14341.SAMPLE1664 14341.SAMPLE1665 14341.SAMPLE1666 14341.SAMPLE1667 
    ##                1                1                1                1 
    ## 14341.SAMPLE1668 14341.SAMPLE1669  14341.SAMPLE167 14341.SAMPLE1670 
    ##                1                1                1                1 
    ## 14341.SAMPLE1671 14341.SAMPLE1672 14341.SAMPLE1673 14341.SAMPLE1674 
    ##                1                1                1                1 
    ## 14341.SAMPLE1675 14341.SAMPLE1676 14341.SAMPLE1677 14341.SAMPLE1678 
    ##                1                1                1                1 
    ## 14341.SAMPLE1679  14341.SAMPLE168 14341.SAMPLE1680 14341.SAMPLE1681 
    ##                1                1                1                1 
    ## 14341.SAMPLE1682 14341.SAMPLE1683 14341.SAMPLE1684 14341.SAMPLE1685 
    ##                1                1                1                1 
    ## 14341.SAMPLE1686 14341.SAMPLE1687 14341.SAMPLE1688 14341.SAMPLE1689 
    ##                1                1                1                1 
    ##  14341.SAMPLE169 14341.SAMPLE1690 14341.SAMPLE1691 14341.SAMPLE1692 
    ##                1                1                1                1 
    ## 14341.SAMPLE1693 14341.SAMPLE1694 14341.SAMPLE1695 14341.SAMPLE1696 
    ##                1                1                1                1 
    ## 14341.SAMPLE1697 14341.SAMPLE1698 14341.SAMPLE1699   14341.SAMPLE17 
    ##                1                1                1                1 
    ##  14341.SAMPLE170 14341.SAMPLE1700 14341.SAMPLE1701 14341.SAMPLE1702 
    ##                1                1                1                1 
    ## 14341.SAMPLE1703 14341.SAMPLE1704 14341.SAMPLE1705 14341.SAMPLE1706 
    ##                1                1                1                1 
    ## 14341.SAMPLE1707 14341.SAMPLE1708 14341.SAMPLE1709  14341.SAMPLE171 
    ##                1                1                1                1 
    ## 14341.SAMPLE1710 14341.SAMPLE1711 14341.SAMPLE1712 14341.SAMPLE1713 
    ##                1                1                1                1 
    ## 14341.SAMPLE1714 14341.SAMPLE1715 14341.SAMPLE1716 14341.SAMPLE1717 
    ##                1                1                1                1 
    ## 14341.SAMPLE1718 14341.SAMPLE1719  14341.SAMPLE172 14341.SAMPLE1720 
    ##                1                1                1                1 
    ## 14341.SAMPLE1721 14341.SAMPLE1722 14341.SAMPLE1723 14341.SAMPLE1724 
    ##                1                1                1                1 
    ## 14341.SAMPLE1725 14341.SAMPLE1726 14341.SAMPLE1727 14341.SAMPLE1728 
    ##                1                1                1                1 
    ## 14341.SAMPLE1729  14341.SAMPLE173 14341.SAMPLE1730 14341.SAMPLE1731 
    ##                1                1                1                1 
    ## 14341.SAMPLE1732 14341.SAMPLE1733 14341.SAMPLE1734 14341.SAMPLE1735 
    ##                1                1                1                1 
    ## 14341.SAMPLE1736 14341.SAMPLE1737 14341.SAMPLE1738 14341.SAMPLE1739 
    ##                1                1                1                1 
    ##  14341.SAMPLE174 14341.SAMPLE1740 14341.SAMPLE1741 14341.SAMPLE1742 
    ##                1                1                1                1 
    ## 14341.SAMPLE1743 14341.SAMPLE1744 14341.SAMPLE1745 14341.SAMPLE1746 
    ##                1                1                1                1 
    ## 14341.SAMPLE1747 14341.SAMPLE1748 14341.SAMPLE1749  14341.SAMPLE175 
    ##                1                1                1                1 
    ## 14341.SAMPLE1750 14341.SAMPLE1751 14341.SAMPLE1752 14341.SAMPLE1753 
    ##                1                1                1                1 
    ## 14341.SAMPLE1754 14341.SAMPLE1755 14341.SAMPLE1756 14341.SAMPLE1757 
    ##                1                1                1                1 
    ## 14341.SAMPLE1758 14341.SAMPLE1759  14341.SAMPLE176 14341.SAMPLE1760 
    ##                1                1                1                1 
    ## 14341.SAMPLE1761 14341.SAMPLE1762 14341.SAMPLE1763 14341.SAMPLE1764 
    ##                1                1                1                1 
    ## 14341.SAMPLE1765 14341.SAMPLE1766 14341.SAMPLE1767 14341.SAMPLE1768 
    ##                1                1                1                1 
    ## 14341.SAMPLE1769  14341.SAMPLE177 14341.SAMPLE1770 14341.SAMPLE1771 
    ##                1                1                1                1 
    ## 14341.SAMPLE1772 14341.SAMPLE1773 14341.SAMPLE1774 14341.SAMPLE1775 
    ##                1                1                1                1 
    ## 14341.SAMPLE1776 14341.SAMPLE1777 14341.SAMPLE1778 14341.SAMPLE1779 
    ##                1                1                1                1 
    ##  14341.SAMPLE178 14341.SAMPLE1780 14341.SAMPLE1781 14341.SAMPLE1782 
    ##                1                1                1                1 
    ## 14341.SAMPLE1783 14341.SAMPLE1784 14341.SAMPLE1785 14341.SAMPLE1786 
    ##                1                1                1                1 
    ## 14341.SAMPLE1787 14341.SAMPLE1788 14341.SAMPLE1789  14341.SAMPLE179 
    ##                1                1                1                1 
    ## 14341.SAMPLE1790 14341.SAMPLE1791 14341.SAMPLE1792 14341.SAMPLE1793 
    ##                1                1                1                1 
    ## 14341.SAMPLE1794 14341.SAMPLE1795 14341.SAMPLE1796 14341.SAMPLE1797 
    ##                1                1                1                1 
    ## 14341.SAMPLE1798 14341.SAMPLE1799   14341.SAMPLE18  14341.SAMPLE180 
    ##                1                1                1                1 
    ## 14341.SAMPLE1800 14341.SAMPLE1801 14341.SAMPLE1802 14341.SAMPLE1803 
    ##                1                1                1                1 
    ## 14341.SAMPLE1804 14341.SAMPLE1805 14341.SAMPLE1806 14341.SAMPLE1807 
    ##                1                1                1                1 
    ## 14341.SAMPLE1808 14341.SAMPLE1809  14341.SAMPLE181 14341.SAMPLE1810 
    ##                1                1                1                1 
    ## 14341.SAMPLE1811 14341.SAMPLE1812 14341.SAMPLE1813 14341.SAMPLE1814 
    ##                1                1                1                1 
    ## 14341.SAMPLE1815 14341.SAMPLE1816 14341.SAMPLE1817 14341.SAMPLE1818 
    ##                1                1                1                1 
    ## 14341.SAMPLE1819  14341.SAMPLE182 14341.SAMPLE1820 14341.SAMPLE1821 
    ##                1                1                1                1 
    ## 14341.SAMPLE1822 14341.SAMPLE1823 14341.SAMPLE1824 14341.SAMPLE1825 
    ##                1                1                1                1 
    ## 14341.SAMPLE1826 14341.SAMPLE1827 14341.SAMPLE1828 14341.SAMPLE1829 
    ##                1                1                1                1 
    ##  14341.SAMPLE183 14341.SAMPLE1830 14341.SAMPLE1831 14341.SAMPLE1832 
    ##                1                1                1                1 
    ## 14341.SAMPLE1833 14341.SAMPLE1834 14341.SAMPLE1835 14341.SAMPLE1836 
    ##                1                1                1                1 
    ## 14341.SAMPLE1837 14341.SAMPLE1838 14341.SAMPLE1839  14341.SAMPLE184 
    ##                1                1                1                1 
    ## 14341.SAMPLE1840 14341.SAMPLE1841 14341.SAMPLE1842 14341.SAMPLE1843 
    ##                1                1                1                1 
    ## 14341.SAMPLE1844 14341.SAMPLE1845 14341.SAMPLE1846 14341.SAMPLE1847 
    ##                1                1                1                1 
    ## 14341.SAMPLE1848 14341.SAMPLE1849  14341.SAMPLE185 14341.SAMPLE1850 
    ##                1                1                1                1 
    ## 14341.SAMPLE1851 14341.SAMPLE1852 14341.SAMPLE1853 14341.SAMPLE1854 
    ##                1                1                1                1 
    ## 14341.SAMPLE1855 14341.SAMPLE1856 14341.SAMPLE1857 14341.SAMPLE1858 
    ##                1                1                1                1 
    ## 14341.SAMPLE1859  14341.SAMPLE186 14341.SAMPLE1860 14341.SAMPLE1861 
    ##                1                1                1                1 
    ## 14341.SAMPLE1862 14341.SAMPLE1863 14341.SAMPLE1864 14341.SAMPLE1865 
    ##                1                1                1                1 
    ## 14341.SAMPLE1866 14341.SAMPLE1867 14341.SAMPLE1868 14341.SAMPLE1869 
    ##                1                1                1                1 
    ##  14341.SAMPLE187 14341.SAMPLE1870 14341.SAMPLE1871 14341.SAMPLE1872 
    ##                1                1                1                1 
    ##  14341.SAMPLE188  14341.SAMPLE189   14341.SAMPLE19  14341.SAMPLE190 
    ##                1                1                1                1 
    ##  14341.SAMPLE191  14341.SAMPLE192  14341.SAMPLE193  14341.SAMPLE194 
    ##                1                1                1                1 
    ##  14341.SAMPLE195  14341.SAMPLE196  14341.SAMPLE197  14341.SAMPLE198 
    ##                1                1                1                1 
    ##  14341.SAMPLE199    14341.SAMPLE2   14341.SAMPLE20  14341.SAMPLE200 
    ##                1                1                1                1 
    ##  14341.SAMPLE201  14341.SAMPLE202  14341.SAMPLE203  14341.SAMPLE204 
    ##                1                1                1                1 
    ##  14341.SAMPLE205  14341.SAMPLE206  14341.SAMPLE207  14341.SAMPLE208 
    ##                1                1                1                1 
    ##  14341.SAMPLE209   14341.SAMPLE21  14341.SAMPLE210  14341.SAMPLE211 
    ##                1                1                1                1 
    ##  14341.SAMPLE212  14341.SAMPLE213  14341.SAMPLE214  14341.SAMPLE215 
    ##                1                1                1                1 
    ##  14341.SAMPLE216  14341.SAMPLE217  14341.SAMPLE218  14341.SAMPLE219 
    ##                1                1                1                1 
    ##   14341.SAMPLE22  14341.SAMPLE220  14341.SAMPLE221  14341.SAMPLE222 
    ##                1                1                1                1 
    ##  14341.SAMPLE223  14341.SAMPLE224  14341.SAMPLE225  14341.SAMPLE226 
    ##                1                1                1                1 
    ##  14341.SAMPLE227  14341.SAMPLE228  14341.SAMPLE229   14341.SAMPLE23 
    ##                1                1                1                1 
    ##  14341.SAMPLE230  14341.SAMPLE231  14341.SAMPLE232  14341.SAMPLE233 
    ##                1                1                1                1 
    ##  14341.SAMPLE234  14341.SAMPLE235  14341.SAMPLE236  14341.SAMPLE237 
    ##                1                1                1                1 
    ##  14341.SAMPLE238  14341.SAMPLE239   14341.SAMPLE24  14341.SAMPLE240 
    ##                1                1                1                1 
    ##  14341.SAMPLE241  14341.SAMPLE242  14341.SAMPLE243  14341.SAMPLE244 
    ##                1                1                1                1 
    ##  14341.SAMPLE245  14341.SAMPLE246  14341.SAMPLE247  14341.SAMPLE248 
    ##                1                1                1                1 
    ##  14341.SAMPLE249   14341.SAMPLE25  14341.SAMPLE250  14341.SAMPLE251 
    ##                1                1                1                1 
    ##  14341.SAMPLE252  14341.SAMPLE253  14341.SAMPLE254  14341.SAMPLE255 
    ##                1                1                1                1 
    ##  14341.SAMPLE256  14341.SAMPLE257  14341.SAMPLE258  14341.SAMPLE259 
    ##                1                1                1                1 
    ##   14341.SAMPLE26  14341.SAMPLE260  14341.SAMPLE261  14341.SAMPLE262 
    ##                1                1                1                1 
    ##  14341.SAMPLE263  14341.SAMPLE264  14341.SAMPLE265  14341.SAMPLE266 
    ##                1                1                1                1 
    ##  14341.SAMPLE267  14341.SAMPLE268  14341.SAMPLE269   14341.SAMPLE27 
    ##                1                1                1                1 
    ##  14341.SAMPLE270  14341.SAMPLE271  14341.SAMPLE272  14341.SAMPLE273 
    ##                1                1                1                1 
    ##  14341.SAMPLE274  14341.SAMPLE275  14341.SAMPLE276  14341.SAMPLE277 
    ##                1                1                1                1 
    ##  14341.SAMPLE278  14341.SAMPLE279   14341.SAMPLE28  14341.SAMPLE280 
    ##                1                1                1                1 
    ##  14341.SAMPLE281  14341.SAMPLE282  14341.SAMPLE283  14341.SAMPLE284 
    ##                1                1                1                1 
    ##  14341.SAMPLE285  14341.SAMPLE286  14341.SAMPLE287  14341.SAMPLE288 
    ##                1                1                1                1 
    ##  14341.SAMPLE289   14341.SAMPLE29  14341.SAMPLE290  14341.SAMPLE291 
    ##                1                1                1                1 
    ##  14341.SAMPLE292  14341.SAMPLE293  14341.SAMPLE294  14341.SAMPLE295 
    ##                1                1                1                1 
    ##  14341.SAMPLE296  14341.SAMPLE297  14341.SAMPLE298  14341.SAMPLE299 
    ##                1                1                1                1 
    ##    14341.SAMPLE3   14341.SAMPLE30  14341.SAMPLE300  14341.SAMPLE301 
    ##                1                1                1                1 
    ##  14341.SAMPLE302  14341.SAMPLE303  14341.SAMPLE304  14341.SAMPLE305 
    ##                1                1                1                1 
    ##  14341.SAMPLE306  14341.SAMPLE307  14341.SAMPLE308  14341.SAMPLE309 
    ##                1                1                1                1 
    ##   14341.SAMPLE31  14341.SAMPLE310  14341.SAMPLE311  14341.SAMPLE312 
    ##                1                1                1                1 
    ##  14341.SAMPLE313  14341.SAMPLE314  14341.SAMPLE315  14341.SAMPLE316 
    ##                1                1                1                1 
    ##  14341.SAMPLE318  14341.SAMPLE319   14341.SAMPLE32  14341.SAMPLE321 
    ##                1                1                1                1 
    ##  14341.SAMPLE322  14341.SAMPLE324  14341.SAMPLE325  14341.SAMPLE327 
    ##                1                1                1                1 
    ##  14341.SAMPLE329   14341.SAMPLE33  14341.SAMPLE330  14341.SAMPLE331 
    ##                1                1                1                1 
    ##  14341.SAMPLE332  14341.SAMPLE333  14341.SAMPLE334  14341.SAMPLE335 
    ##                1                1                1                1 
    ##  14341.SAMPLE336  14341.SAMPLE337  14341.SAMPLE338  14341.SAMPLE339 
    ##                1                1                1                1 
    ##   14341.SAMPLE34  14341.SAMPLE340  14341.SAMPLE341  14341.SAMPLE342 
    ##                1                1                1                1 
    ##  14341.SAMPLE343  14341.SAMPLE344  14341.SAMPLE345  14341.SAMPLE346 
    ##                1                1                1                1 
    ##  14341.SAMPLE347  14341.SAMPLE348  14341.SAMPLE349   14341.SAMPLE35 
    ##                1                1                1                1 
    ##  14341.SAMPLE350  14341.SAMPLE351  14341.SAMPLE352  14341.SAMPLE353 
    ##                1                1                1                1 
    ##  14341.SAMPLE354  14341.SAMPLE355  14341.SAMPLE356  14341.SAMPLE357 
    ##                1                1                1                1 
    ##  14341.SAMPLE358  14341.SAMPLE359   14341.SAMPLE36  14341.SAMPLE360 
    ##                1                1                1                1 
    ##  14341.SAMPLE361  14341.SAMPLE362  14341.SAMPLE363  14341.SAMPLE364 
    ##                1                1                1                1 
    ##  14341.SAMPLE365  14341.SAMPLE366  14341.SAMPLE367  14341.SAMPLE368 
    ##                1                1                1                1 
    ##  14341.SAMPLE369   14341.SAMPLE37  14341.SAMPLE370  14341.SAMPLE371 
    ##                1                1                1                1 
    ##  14341.SAMPLE372  14341.SAMPLE373  14341.SAMPLE374  14341.SAMPLE375 
    ##                1                1                1                1 
    ##  14341.SAMPLE376  14341.SAMPLE377  14341.SAMPLE378  14341.SAMPLE379 
    ##                1                1                1                1 
    ##   14341.SAMPLE38  14341.SAMPLE380  14341.SAMPLE381  14341.SAMPLE382 
    ##                1                1                1                1 
    ##  14341.SAMPLE383  14341.SAMPLE384  14341.SAMPLE385  14341.SAMPLE386 
    ##                1                1                1                1 
    ##  14341.SAMPLE387  14341.SAMPLE388  14341.SAMPLE389   14341.SAMPLE39 
    ##                1                1                1                1 
    ##  14341.SAMPLE390  14341.SAMPLE391  14341.SAMPLE392  14341.SAMPLE393 
    ##                1                1                1                1 
    ##  14341.SAMPLE394  14341.SAMPLE395  14341.SAMPLE396  14341.SAMPLE397 
    ##                1                1                1                1 
    ##  14341.SAMPLE398  14341.SAMPLE399    14341.SAMPLE4   14341.SAMPLE40 
    ##                1                1                1                1 
    ##  14341.SAMPLE400  14341.SAMPLE401  14341.SAMPLE402  14341.SAMPLE403 
    ##                1                1                1                1 
    ##  14341.SAMPLE404  14341.SAMPLE405  14341.SAMPLE406  14341.SAMPLE407 
    ##                1                1                1                1 
    ##  14341.SAMPLE408  14341.SAMPLE409   14341.SAMPLE41  14341.SAMPLE410 
    ##                1                1                1                1 
    ##  14341.SAMPLE411  14341.SAMPLE412  14341.SAMPLE413  14341.SAMPLE414 
    ##                1                1                1                1 
    ##  14341.SAMPLE415  14341.SAMPLE416  14341.SAMPLE417  14341.SAMPLE418 
    ##                1                1                1                1 
    ##  14341.SAMPLE419   14341.SAMPLE42  14341.SAMPLE420  14341.SAMPLE421 
    ##                1                1                1                1 
    ##  14341.SAMPLE422  14341.SAMPLE423  14341.SAMPLE424  14341.SAMPLE425 
    ##                1                1                1                1 
    ##  14341.SAMPLE426  14341.SAMPLE427  14341.SAMPLE428  14341.SAMPLE429 
    ##                1                1                1                1 
    ##   14341.SAMPLE43  14341.SAMPLE430  14341.SAMPLE431  14341.SAMPLE432 
    ##                1                1                1                1 
    ##  14341.SAMPLE433  14341.SAMPLE434  14341.SAMPLE435  14341.SAMPLE436 
    ##                1                1                1                1 
    ##  14341.SAMPLE437  14341.SAMPLE438  14341.SAMPLE439   14341.SAMPLE44 
    ##                1                1                1                1 
    ##  14341.SAMPLE440  14341.SAMPLE441  14341.SAMPLE442  14341.SAMPLE443 
    ##                1                1                1                1 
    ##  14341.SAMPLE444  14341.SAMPLE445  14341.SAMPLE446  14341.SAMPLE447 
    ##                1                1                1                1 
    ##  14341.SAMPLE448  14341.SAMPLE449   14341.SAMPLE45  14341.SAMPLE450 
    ##                1                1                1                1 
    ##  14341.SAMPLE451  14341.SAMPLE452  14341.SAMPLE453  14341.SAMPLE454 
    ##                1                1                1                1 
    ##  14341.SAMPLE455  14341.SAMPLE456  14341.SAMPLE457  14341.SAMPLE458 
    ##                1                1                1                1 
    ##  14341.SAMPLE459   14341.SAMPLE46  14341.SAMPLE460  14341.SAMPLE461 
    ##                1                1                1                1 
    ##  14341.SAMPLE462  14341.SAMPLE463  14341.SAMPLE464  14341.SAMPLE465 
    ##                1                1                1                1 
    ##  14341.SAMPLE466  14341.SAMPLE467  14341.SAMPLE468  14341.SAMPLE469 
    ##                1                1                1                1 
    ##   14341.SAMPLE47  14341.SAMPLE470  14341.SAMPLE471  14341.SAMPLE472 
    ##                1                1                1                1 
    ##  14341.SAMPLE473  14341.SAMPLE474  14341.SAMPLE475  14341.SAMPLE476 
    ##                1                1                1                1 
    ##  14341.SAMPLE477  14341.SAMPLE478  14341.SAMPLE479   14341.SAMPLE48 
    ##                1                1                1                1 
    ##  14341.SAMPLE480  14341.SAMPLE481  14341.SAMPLE482  14341.SAMPLE483 
    ##                1                1                1                1 
    ##  14341.SAMPLE484  14341.SAMPLE485  14341.SAMPLE486  14341.SAMPLE487 
    ##                1                1                1                1 
    ##  14341.SAMPLE488  14341.SAMPLE489   14341.SAMPLE49  14341.SAMPLE490 
    ##                1                1                1                1 
    ##  14341.SAMPLE491  14341.SAMPLE492  14341.SAMPLE493  14341.SAMPLE494 
    ##                1                1                1                1 
    ##  14341.SAMPLE495  14341.SAMPLE496  14341.SAMPLE497  14341.SAMPLE498 
    ##                1                1                1                1 
    ##  14341.SAMPLE499    14341.SAMPLE5   14341.SAMPLE50  14341.SAMPLE500 
    ##                1                1                1                1 
    ##  14341.SAMPLE501  14341.SAMPLE502  14341.SAMPLE503  14341.SAMPLE504 
    ##                1                1                1                1 
    ##  14341.SAMPLE505  14341.SAMPLE506  14341.SAMPLE507  14341.SAMPLE508 
    ##                1                1                1                1 
    ##  14341.SAMPLE509   14341.SAMPLE51  14341.SAMPLE510  14341.SAMPLE511 
    ##                1                1                1                1 
    ##  14341.SAMPLE512  14341.SAMPLE513  14341.SAMPLE514  14341.SAMPLE515 
    ##                1                1                1                1 
    ##  14341.SAMPLE516  14341.SAMPLE517  14341.SAMPLE518  14341.SAMPLE519 
    ##                1                1                1                1 
    ##   14341.SAMPLE52  14341.SAMPLE520  14341.SAMPLE521  14341.SAMPLE522 
    ##                1                1                1                1 
    ##  14341.SAMPLE523  14341.SAMPLE524  14341.SAMPLE525  14341.SAMPLE526 
    ##                1                1                1                1 
    ##  14341.SAMPLE527  14341.SAMPLE528  14341.SAMPLE529   14341.SAMPLE53 
    ##                1                1                1                1 
    ##  14341.SAMPLE530  14341.SAMPLE531  14341.SAMPLE532  14341.SAMPLE533 
    ##                1                1                1                1 
    ##  14341.SAMPLE534  14341.SAMPLE535  14341.SAMPLE536  14341.SAMPLE537 
    ##                1                1                1                1 
    ##  14341.SAMPLE538  14341.SAMPLE539   14341.SAMPLE54  14341.SAMPLE540 
    ##                1                1                1                1 
    ##  14341.SAMPLE541  14341.SAMPLE542  14341.SAMPLE543  14341.SAMPLE544 
    ##                1                1                1                1 
    ##  14341.SAMPLE545  14341.SAMPLE546  14341.SAMPLE547  14341.SAMPLE548 
    ##                1                1                1                1 
    ##  14341.SAMPLE549   14341.SAMPLE55  14341.SAMPLE550  14341.SAMPLE551 
    ##                1                1                1                1 
    ##  14341.SAMPLE552  14341.SAMPLE553  14341.SAMPLE554  14341.SAMPLE555 
    ##                1                1                1                1 
    ##  14341.SAMPLE556  14341.SAMPLE557  14341.SAMPLE558  14341.SAMPLE559 
    ##                1                1                1                1 
    ##   14341.SAMPLE56  14341.SAMPLE560  14341.SAMPLE561  14341.SAMPLE562 
    ##                1                1                1                1 
    ##  14341.SAMPLE563  14341.SAMPLE564  14341.SAMPLE565  14341.SAMPLE566 
    ##                1                1                1                1 
    ##  14341.SAMPLE567  14341.SAMPLE568  14341.SAMPLE569   14341.SAMPLE57 
    ##                1                1                1                1 
    ##  14341.SAMPLE570  14341.SAMPLE571  14341.SAMPLE572  14341.SAMPLE573 
    ##                1                1                1                1 
    ##  14341.SAMPLE574  14341.SAMPLE575  14341.SAMPLE576  14341.SAMPLE577 
    ##                1                1                1                1 
    ##  14341.SAMPLE578  14341.SAMPLE579   14341.SAMPLE58  14341.SAMPLE580 
    ##                1                1                1                1 
    ##  14341.SAMPLE581  14341.SAMPLE582  14341.SAMPLE583  14341.SAMPLE584 
    ##                1                1                1                1 
    ##  14341.SAMPLE585  14341.SAMPLE586  14341.SAMPLE587  14341.SAMPLE588 
    ##                1                1                1                1 
    ##  14341.SAMPLE589   14341.SAMPLE59  14341.SAMPLE590  14341.SAMPLE591 
    ##                1                1                1                1 
    ##  14341.SAMPLE592  14341.SAMPLE593  14341.SAMPLE594  14341.SAMPLE595 
    ##                1                1                1                1 
    ##  14341.SAMPLE596  14341.SAMPLE597  14341.SAMPLE598  14341.SAMPLE599 
    ##                1                1                1                1 
    ##    14341.SAMPLE6   14341.SAMPLE60  14341.SAMPLE600  14341.SAMPLE601 
    ##                1                1                1                1 
    ##  14341.SAMPLE602  14341.SAMPLE603  14341.SAMPLE604  14341.SAMPLE605 
    ##                1                1                1                1 
    ##  14341.SAMPLE606  14341.SAMPLE607  14341.SAMPLE608  14341.SAMPLE609 
    ##                1                1                1                1 
    ##   14341.SAMPLE61  14341.SAMPLE610  14341.SAMPLE611  14341.SAMPLE612 
    ##                1                1                1                1 
    ##  14341.SAMPLE613  14341.SAMPLE614  14341.SAMPLE615  14341.SAMPLE616 
    ##                1                1                1                1 
    ##  14341.SAMPLE617  14341.SAMPLE618  14341.SAMPLE619   14341.SAMPLE62 
    ##                1                1                1                1 
    ##  14341.SAMPLE620  14341.SAMPLE621  14341.SAMPLE622  14341.SAMPLE623 
    ##                1                1                1                1 
    ##  14341.SAMPLE624  14341.SAMPLE625  14341.SAMPLE626  14341.SAMPLE627 
    ##                1                1                1                1 
    ##  14341.SAMPLE628  14341.SAMPLE629   14341.SAMPLE63  14341.SAMPLE630 
    ##                1                1                1                1 
    ##  14341.SAMPLE631  14341.SAMPLE632  14341.SAMPLE633  14341.SAMPLE634 
    ##                1                1                1                1 
    ##  14341.SAMPLE635  14341.SAMPLE636  14341.SAMPLE637  14341.SAMPLE638 
    ##                1                1                1                1 
    ##  14341.SAMPLE639   14341.SAMPLE64  14341.SAMPLE640  14341.SAMPLE641 
    ##                1                1                1                1 
    ##  14341.SAMPLE642  14341.SAMPLE643  14341.SAMPLE644  14341.SAMPLE645 
    ##                1                1                1                1 
    ##  14341.SAMPLE646  14341.SAMPLE647  14341.SAMPLE648  14341.SAMPLE649 
    ##                1                1                1                1 
    ##   14341.SAMPLE65  14341.SAMPLE650  14341.SAMPLE651  14341.SAMPLE652 
    ##                1                1                1                1 
    ##  14341.SAMPLE653  14341.SAMPLE654  14341.SAMPLE655  14341.SAMPLE656 
    ##                1                1                1                1 
    ##  14341.SAMPLE657  14341.SAMPLE658  14341.SAMPLE659   14341.SAMPLE66 
    ##                1                1                1                1 
    ##  14341.SAMPLE660  14341.SAMPLE661  14341.SAMPLE662  14341.SAMPLE663 
    ##                1                1                1                1 
    ##  14341.SAMPLE664  14341.SAMPLE665  14341.SAMPLE666  14341.SAMPLE667 
    ##                1                1                1                1 
    ##  14341.SAMPLE668  14341.SAMPLE669   14341.SAMPLE67  14341.SAMPLE670 
    ##                1                1                1                1 
    ##  14341.SAMPLE671  14341.SAMPLE672  14341.SAMPLE673  14341.SAMPLE674 
    ##                1                1                1                1 
    ##  14341.SAMPLE675  14341.SAMPLE676  14341.SAMPLE677  14341.SAMPLE678 
    ##                1                1                1                1 
    ##  14341.SAMPLE679   14341.SAMPLE68  14341.SAMPLE680  14341.SAMPLE681 
    ##                1                1                1                1 
    ##  14341.SAMPLE682  14341.SAMPLE683  14341.SAMPLE684  14341.SAMPLE685 
    ##                1                1                1                1 
    ##  14341.SAMPLE686  14341.SAMPLE687  14341.SAMPLE688  14341.SAMPLE689 
    ##                1                1                1                1 
    ##   14341.SAMPLE69  14341.SAMPLE690  14341.SAMPLE691  14341.SAMPLE692 
    ##                1                1                1                1 
    ##  14341.SAMPLE693  14341.SAMPLE694  14341.SAMPLE695  14341.SAMPLE696 
    ##                1                1                1                1 
    ##  14341.SAMPLE697  14341.SAMPLE698  14341.SAMPLE699    14341.SAMPLE7 
    ##                1                1                1                1 
    ##   14341.SAMPLE70  14341.SAMPLE700  14341.SAMPLE701  14341.SAMPLE702 
    ##                1                1                1                1 
    ##  14341.SAMPLE703  14341.SAMPLE704  14341.SAMPLE705  14341.SAMPLE706 
    ##                1                1                1                1 
    ##  14341.SAMPLE707  14341.SAMPLE708  14341.SAMPLE709   14341.SAMPLE71 
    ##                1                1                1                1 
    ##  14341.SAMPLE710  14341.SAMPLE711  14341.SAMPLE712  14341.SAMPLE713 
    ##                1                1                1                1 
    ##  14341.SAMPLE714  14341.SAMPLE715  14341.SAMPLE716  14341.SAMPLE717 
    ##                1                1                1                1 
    ##  14341.SAMPLE718  14341.SAMPLE719   14341.SAMPLE72  14341.SAMPLE720 
    ##                1                1                1                1 
    ##  14341.SAMPLE721  14341.SAMPLE722  14341.SAMPLE723  14341.SAMPLE724 
    ##                1                1                1                1 
    ##  14341.SAMPLE725  14341.SAMPLE726  14341.SAMPLE727  14341.SAMPLE728 
    ##                1                1                1                1 
    ##  14341.SAMPLE729   14341.SAMPLE73  14341.SAMPLE730  14341.SAMPLE731 
    ##                1                1                1                1 
    ##  14341.SAMPLE732  14341.SAMPLE733  14341.SAMPLE734  14341.SAMPLE735 
    ##                1                1                1                1 
    ##  14341.SAMPLE736  14341.SAMPLE737  14341.SAMPLE738  14341.SAMPLE739 
    ##                1                1                1                1 
    ##   14341.SAMPLE74  14341.SAMPLE740  14341.SAMPLE741  14341.SAMPLE742 
    ##                1                1                1                1 
    ##  14341.SAMPLE743  14341.SAMPLE744  14341.SAMPLE745  14341.SAMPLE746 
    ##                1                1                1                1 
    ##  14341.SAMPLE747  14341.SAMPLE748  14341.SAMPLE749   14341.SAMPLE75 
    ##                1                1                1                1 
    ##  14341.SAMPLE750  14341.SAMPLE751  14341.SAMPLE752  14341.SAMPLE753 
    ##                1                1                1                1 
    ##  14341.SAMPLE754  14341.SAMPLE755  14341.SAMPLE756  14341.SAMPLE757 
    ##                1                1                1                1 
    ##  14341.SAMPLE758  14341.SAMPLE759   14341.SAMPLE76  14341.SAMPLE760 
    ##                1                1                1                1 
    ##  14341.SAMPLE761  14341.SAMPLE762  14341.SAMPLE763  14341.SAMPLE764 
    ##                1                1                1                1 
    ##  14341.SAMPLE765  14341.SAMPLE766  14341.SAMPLE767  14341.SAMPLE768 
    ##                1                1                1                1 
    ##  14341.SAMPLE769   14341.SAMPLE77  14341.SAMPLE770  14341.SAMPLE771 
    ##                1                1                1                1 
    ##  14341.SAMPLE772  14341.SAMPLE773  14341.SAMPLE774  14341.SAMPLE775 
    ##                1                1                1                1 
    ##  14341.SAMPLE776  14341.SAMPLE777  14341.SAMPLE778  14341.SAMPLE779 
    ##                1                1                1                1 
    ##   14341.SAMPLE78  14341.SAMPLE780  14341.SAMPLE781  14341.SAMPLE782 
    ##                1                1                1                1 
    ##  14341.SAMPLE783  14341.SAMPLE784  14341.SAMPLE785  14341.SAMPLE786 
    ##                1                1                1                1 
    ##  14341.SAMPLE787  14341.SAMPLE788  14341.SAMPLE789   14341.SAMPLE79 
    ##                1                1                1                1 
    ##  14341.SAMPLE790  14341.SAMPLE791  14341.SAMPLE792  14341.SAMPLE793 
    ##                1                1                1                1 
    ##  14341.SAMPLE794  14341.SAMPLE795  14341.SAMPLE796  14341.SAMPLE797 
    ##                1                1                1                1 
    ##  14341.SAMPLE798  14341.SAMPLE799    14341.SAMPLE8   14341.SAMPLE80 
    ##                1                1                1                1 
    ##  14341.SAMPLE800  14341.SAMPLE801  14341.SAMPLE802  14341.SAMPLE803 
    ##                1                1                1                1 
    ##  14341.SAMPLE804  14341.SAMPLE805  14341.SAMPLE806  14341.SAMPLE807 
    ##                1                1                1                1 
    ##  14341.SAMPLE808  14341.SAMPLE809   14341.SAMPLE81  14341.SAMPLE810 
    ##                1                1                1                1 
    ##  14341.SAMPLE811  14341.SAMPLE812  14341.SAMPLE813  14341.SAMPLE814 
    ##                1                1                1                1 
    ##  14341.SAMPLE815  14341.SAMPLE816  14341.SAMPLE817  14341.SAMPLE818 
    ##                1                1                1                1 
    ##  14341.SAMPLE819   14341.SAMPLE82  14341.SAMPLE820  14341.SAMPLE821 
    ##                1                1                1                1 
    ##  14341.SAMPLE822  14341.SAMPLE823  14341.SAMPLE824  14341.SAMPLE825 
    ##                1                1                1                1 
    ##  14341.SAMPLE826  14341.SAMPLE827  14341.SAMPLE828  14341.SAMPLE829 
    ##                1                1                1                1 
    ##   14341.SAMPLE83  14341.SAMPLE830  14341.SAMPLE831  14341.SAMPLE832 
    ##                1                1                1                1 
    ##  14341.SAMPLE833  14341.SAMPLE834  14341.SAMPLE835  14341.SAMPLE836 
    ##                1                1                1                1 
    ##  14341.SAMPLE837  14341.SAMPLE838  14341.SAMPLE839   14341.SAMPLE84 
    ##                1                1                1                1 
    ##  14341.SAMPLE840  14341.SAMPLE841  14341.SAMPLE842  14341.SAMPLE843 
    ##                1                1                1                1 
    ##  14341.SAMPLE844  14341.SAMPLE845  14341.SAMPLE846  14341.SAMPLE847 
    ##                1                1                1                1 
    ##  14341.SAMPLE848  14341.SAMPLE849   14341.SAMPLE85  14341.SAMPLE850 
    ##                1                1                1                1 
    ##  14341.SAMPLE851  14341.SAMPLE852  14341.SAMPLE853  14341.SAMPLE854 
    ##                1                1                1                1 
    ##  14341.SAMPLE855  14341.SAMPLE856  14341.SAMPLE857  14341.SAMPLE858 
    ##                1                1                1                1 
    ##  14341.SAMPLE859   14341.SAMPLE86  14341.SAMPLE860  14341.SAMPLE861 
    ##                1                1                1                1 
    ##  14341.SAMPLE862  14341.SAMPLE863  14341.SAMPLE864  14341.SAMPLE865 
    ##                1                1                1                1 
    ##  14341.SAMPLE866  14341.SAMPLE867  14341.SAMPLE868  14341.SAMPLE869 
    ##                1                1                1                1 
    ##   14341.SAMPLE87  14341.SAMPLE870  14341.SAMPLE871  14341.SAMPLE872 
    ##                1                1                1                1 
    ##  14341.SAMPLE873  14341.SAMPLE874  14341.SAMPLE875  14341.SAMPLE876 
    ##                1                1                1                1 
    ##  14341.SAMPLE877  14341.SAMPLE878  14341.SAMPLE879   14341.SAMPLE88 
    ##                1                1                1                1 
    ##  14341.SAMPLE880  14341.SAMPLE881  14341.SAMPLE882  14341.SAMPLE883 
    ##                1                1                1                1 
    ##  14341.SAMPLE884  14341.SAMPLE885  14341.SAMPLE886  14341.SAMPLE887 
    ##                1                1                1                1 
    ##  14341.SAMPLE888  14341.SAMPLE889   14341.SAMPLE89  14341.SAMPLE890 
    ##                1                1                1                1 
    ##  14341.SAMPLE891  14341.SAMPLE892  14341.SAMPLE893  14341.SAMPLE894 
    ##                1                1                1                1 
    ##  14341.SAMPLE895  14341.SAMPLE896  14341.SAMPLE897  14341.SAMPLE898 
    ##                1                1                1                1 
    ##  14341.SAMPLE899    14341.SAMPLE9   14341.SAMPLE90  14341.SAMPLE900 
    ##                1                1                1                1 
    ##  14341.SAMPLE901  14341.SAMPLE902  14341.SAMPLE903  14341.SAMPLE904 
    ##                1                1                1                1 
    ##  14341.SAMPLE905  14341.SAMPLE906  14341.SAMPLE907  14341.SAMPLE908 
    ##                1                1                1                1 
    ##  14341.SAMPLE909   14341.SAMPLE91  14341.SAMPLE910  14341.SAMPLE911 
    ##                1                1                1                1 
    ##  14341.SAMPLE912  14341.SAMPLE913  14341.SAMPLE914  14341.SAMPLE915 
    ##                1                1                1                1 
    ##  14341.SAMPLE916  14341.SAMPLE917  14341.SAMPLE918  14341.SAMPLE919 
    ##                1                1                1                1 
    ##   14341.SAMPLE92  14341.SAMPLE920  14341.SAMPLE921  14341.SAMPLE922 
    ##                1                1                1                1 
    ##  14341.SAMPLE923  14341.SAMPLE924  14341.SAMPLE925  14341.SAMPLE926 
    ##                1                1                1                1 
    ##  14341.SAMPLE927  14341.SAMPLE928  14341.SAMPLE929   14341.SAMPLE93 
    ##                1                1                1                1 
    ##  14341.SAMPLE930  14341.SAMPLE931  14341.SAMPLE932  14341.SAMPLE933 
    ##                1                1                1                1 
    ##  14341.SAMPLE934  14341.SAMPLE935  14341.SAMPLE936  14341.SAMPLE937 
    ##                1                1                1                1 
    ##  14341.SAMPLE938  14341.SAMPLE939   14341.SAMPLE94  14341.SAMPLE940 
    ##                1                1                1                1 
    ##  14341.SAMPLE941  14341.SAMPLE942  14341.SAMPLE943  14341.SAMPLE944 
    ##                1                1                1                1 
    ##  14341.SAMPLE945  14341.SAMPLE946  14341.SAMPLE947  14341.SAMPLE948 
    ##                1                1                1                1 
    ##  14341.SAMPLE949   14341.SAMPLE95  14341.SAMPLE950  14341.SAMPLE951 
    ##                1                1                1                1 
    ##  14341.SAMPLE952  14341.SAMPLE953  14341.SAMPLE954  14341.SAMPLE955 
    ##                1                1                1                1 
    ##  14341.SAMPLE956  14341.SAMPLE957  14341.SAMPLE958  14341.SAMPLE959 
    ##                1                1                1                1 
    ##   14341.SAMPLE96  14341.SAMPLE960  14341.SAMPLE961  14341.SAMPLE962 
    ##                1                1                1                1 
    ##  14341.SAMPLE963  14341.SAMPLE964  14341.SAMPLE965  14341.SAMPLE966 
    ##                1                1                1                1 
    ##  14341.SAMPLE967  14341.SAMPLE968  14341.SAMPLE969   14341.SAMPLE97 
    ##                1                1                1                1 
    ##  14341.SAMPLE970  14341.SAMPLE971  14341.SAMPLE972  14341.SAMPLE973 
    ##                1                1                1                1 
    ##  14341.SAMPLE974  14341.SAMPLE975  14341.SAMPLE976  14341.SAMPLE977 
    ##                1                1                1                1 
    ##  14341.SAMPLE978  14341.SAMPLE979   14341.SAMPLE98  14341.SAMPLE980 
    ##                1                1                1                1 
    ##  14341.SAMPLE981  14341.SAMPLE982  14341.SAMPLE983  14341.SAMPLE984 
    ##                1                1                1                1 
    ##  14341.SAMPLE985  14341.SAMPLE986  14341.SAMPLE987  14341.SAMPLE988 
    ##                1                1                1                1 
    ##  14341.SAMPLE989   14341.SAMPLE99  14341.SAMPLE990  14341.SAMPLE991 
    ##                1                1                1                1 
    ##  14341.SAMPLE992  14341.SAMPLE993  14341.SAMPLE994  14341.SAMPLE995 
    ##                1                1                1                1 
    ##  14341.SAMPLE996  14341.SAMPLE997  14341.SAMPLE998  14341.SAMPLE999 
    ##                1                1                1                1

``` r
table(metadata$Site)
```

    ## 
    ## dental plaque 
    ##          2528

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
    ##      14341.A0001      14341.A0002      14341.A0003      14341.A0004 
    ##                1                1                1                1 
    ##      14341.A0005      14341.A0006      14341.A0007      14341.A0008 
    ##                1                1                1                1 
    ##      14341.A0009      14341.A0010      14341.A0011      14341.A0012 
    ##                1                1                1                1 
    ##      14341.A0013      14341.A0014      14341.A0015      14341.A0016 
    ##                1                1                1                1 
    ##      14341.A0017      14341.A0018      14341.A0019      14341.A0020 
    ##                1                1                1                1 
    ##      14341.A0021      14341.A0022      14341.A0023      14341.A0024 
    ##                1                1                1                1 
    ##      14341.A0025      14341.A0026      14341.A0027      14341.A0028 
    ##                1                1                1                1 
    ##      14341.A0029      14341.A0030      14341.A0031      14341.A0032 
    ##                1                1                1                1 
    ##      14341.A0033      14341.A0034      14341.A0035      14341.A0036 
    ##                1                1                1                1 
    ##      14341.A0037      14341.A0038      14341.A0039      14341.A0040 
    ##                1                1                1                1 
    ##      14341.A0041      14341.A0042      14341.A0043      14341.A0044 
    ##                1                1                1                1 
    ##      14341.A0045      14341.A0046      14341.A0047      14341.A0048 
    ##                1                1                1                1 
    ##      14341.A0049      14341.A0050      14341.A0051      14341.A0052 
    ##                1                1                1                1 
    ##      14341.A0053      14341.A0054      14341.A0055      14341.A0056 
    ##                1                1                1                1 
    ##      14341.A0057      14341.A0058      14341.A0059      14341.A0060 
    ##                1                1                1                1 
    ##      14341.A0061      14341.A0062      14341.A0063      14341.A0064 
    ##                1                1                1                1 
    ##      14341.A0065      14341.A0066      14341.A0067      14341.A0068 
    ##                1                1                1                1 
    ##      14341.A0069      14341.A0070      14341.A0071      14341.A0072 
    ##                1                1                1                1 
    ##      14341.A0073      14341.A0074      14341.A0075      14341.A0076 
    ##                1                1                1                1 
    ##      14341.A0077      14341.A0078      14341.A0079      14341.A0080 
    ##                1                1                1                1 
    ##      14341.A0081      14341.A0082      14341.A0083      14341.A0084 
    ##                1                1                1                1 
    ##      14341.A0085      14341.A0086      14341.A0087      14341.A0088 
    ##                1                1                1                1 
    ##      14341.A0089      14341.A0090      14341.A0091      14341.A0092 
    ##                1                1                1                1 
    ##      14341.A0093      14341.A0094      14341.A0095      14341.A0096 
    ##                1                1                1                1 
    ##      14341.A0097      14341.A0098      14341.A0099      14341.A0100 
    ##                1                1                1                1 
    ##      14341.A0101      14341.A0102      14341.A0103      14341.A0104 
    ##                1                1                1                1 
    ##      14341.A0105      14341.A0106      14341.A0107      14341.A0108 
    ##                1                1                1                1 
    ##      14341.A0109      14341.A0110      14341.A0111      14341.A0112 
    ##                1                1                1                1 
    ##      14341.A0113      14341.A0114      14341.A0115      14341.A0116 
    ##                1                1                1                1 
    ##      14341.A0117      14341.A0118      14341.A0119      14341.A0120 
    ##                1                1                1                1 
    ##      14341.A0121      14341.A0122      14341.A0123      14341.A0124 
    ##                1                1                1                1 
    ##      14341.A0125      14341.A0126      14341.A0127      14341.A0128 
    ##                1                1                1                1 
    ##      14341.A0129      14341.A0130      14341.A0131      14341.A0132 
    ##                1                1                1                1 
    ##      14341.A0133      14341.A0134      14341.A0135      14341.A0136 
    ##                1                1                1                1 
    ##      14341.A0137      14341.A0138      14341.A0139      14341.A0140 
    ##                1                1                1                1 
    ##      14341.A0141      14341.A0142      14341.A0143      14341.A0144 
    ##                1                1                1                1 
    ##      14341.A0145      14341.A0146      14341.A0147      14341.A0148 
    ##                1                1                1                1 
    ##      14341.A0149      14341.A0150      14341.A0151      14341.A0152 
    ##                1                1                1                1 
    ##      14341.A0153      14341.A0154      14341.A0155      14341.A0156 
    ##                1                1                1                1 
    ##      14341.A0157      14341.A0158      14341.A0159      14341.A0160 
    ##                1                1                1                1 
    ##      14341.A0161      14341.A0162      14341.A0163      14341.A0164 
    ##                1                1                1                1 
    ##      14341.A0165      14341.A0166      14341.A0167      14341.A0168 
    ##                1                1                1                1 
    ##      14341.A0169      14341.A0170      14341.A0171      14341.A0172 
    ##                1                1                1                1 
    ##      14341.A0173      14341.A0174      14341.A0175      14341.A0176 
    ##                1                1                1                1 
    ##      14341.A0177      14341.A0178      14341.A0179      14341.A0180 
    ##                1                1                1                1 
    ##      14341.A0181      14341.A0182      14341.A0183      14341.A0184 
    ##                1                1                1                1 
    ##      14341.A0185      14341.A0186      14341.A0187      14341.A0188 
    ##                1                1                1                1 
    ##      14341.A0189      14341.A0190      14341.A0191      14341.A0192 
    ##                1                1                1                1 
    ##      14341.A0193      14341.A0194      14341.A0195      14341.A0196 
    ##                1                1                1                1 
    ##      14341.A0197      14341.A0198      14341.A0199      14341.A0200 
    ##                1                1                1                1 
    ##      14341.A0201      14341.A0202      14341.A0203      14341.A0204 
    ##                1                1                1                1 
    ##      14341.A0205      14341.A0206      14341.A0207      14341.A0208 
    ##                1                1                1                1 
    ##      14341.A0209      14341.A0210      14341.A0211      14341.A0212 
    ##                1                1                1                1 
    ##      14341.A0213      14341.A0214      14341.A0215      14341.A0216 
    ##                1                1                1                1 
    ##      14341.A0217      14341.A0218      14341.A0219      14341.A0220 
    ##                1                1                1                1 
    ##      14341.A0221      14341.A0222      14341.A0223      14341.A0224 
    ##                1                1                1                1 
    ##      14341.A0225      14341.A0226      14341.A0227      14341.A0228 
    ##                1                1                1                1 
    ##      14341.A0229      14341.A0230      14341.A0231      14341.A0232 
    ##                1                1                1                1 
    ##      14341.A0233      14341.A0234      14341.A0235      14341.A0236 
    ##                1                1                1                1 
    ##      14341.A0237      14341.A0238      14341.A0239      14341.A0240 
    ##                1                1                1                1 
    ##      14341.A0241      14341.A0242      14341.A0243      14341.A0244 
    ##                1                1                1                1 
    ##      14341.A0245      14341.A0246      14341.A0247      14341.A0248 
    ##                1                1                1                1 
    ##      14341.A0249      14341.A0250      14341.A0251      14341.A0252 
    ##                1                1                1                1 
    ##      14341.A0253      14341.A0254      14341.A0255      14341.A0256 
    ##                1                1                1                1 
    ##      14341.A0257      14341.A0258      14341.A0259      14341.A0260 
    ##                1                1                1                1 
    ##      14341.A0261      14341.A0262      14341.A0263      14341.A0264 
    ##                1                1                1                1 
    ##      14341.A0265      14341.A0266      14341.A0267      14341.A0268 
    ##                1                1                1                1 
    ##      14341.A0269      14341.A0270      14341.A0271      14341.A0272 
    ##                1                1                1                1 
    ##      14341.A0273      14341.A0274      14341.A0275      14341.A0276 
    ##                1                1                1                1 
    ##      14341.A0277      14341.A0278      14341.A0279      14341.A0280 
    ##                1                1                1                1 
    ##      14341.A0281      14341.A0282      14341.A0283      14341.A0284 
    ##                1                1                1                1 
    ##      14341.A0285      14341.A0286      14341.A0287      14341.A0288 
    ##                1                1                1                1 
    ##      14341.A0289      14341.A0290      14341.A0291      14341.A0292 
    ##                1                1                1                1 
    ##      14341.A0293      14341.A0294      14341.A0295      14341.A0296 
    ##                1                1                1                1 
    ##      14341.A0297      14341.A0298      14341.A0299      14341.A0300 
    ##                1                1                1                1 
    ##      14341.A0301      14341.A0302      14341.A0303      14341.A0304 
    ##                1                1                1                1 
    ##      14341.A0305      14341.A0306      14341.A0307      14341.A0308 
    ##                1                1                1                1 
    ##      14341.A0309      14341.A0310      14341.A0311      14341.A0312 
    ##                1                1                1                1 
    ##      14341.A0313      14341.A0314      14341.A0315      14341.A0316 
    ##                1                1                1                1 
    ##      14341.A0317      14341.A0318      14341.A0319      14341.A0320 
    ##                1                1                1                1 
    ##      14341.A0321      14341.A0322      14341.A0323      14341.A0324 
    ##                1                1                1                1 
    ##      14341.A0325      14341.A0326      14341.A0327      14341.A0328 
    ##                1                1                1                1 
    ##      14341.A0329      14341.A0330      14341.A0331      14341.A0332 
    ##                1                1                1                1 
    ##      14341.A0333      14341.A0334      14341.A0335      14341.A0336 
    ##                1                1                1                1 
    ##      14341.A0337      14341.A0338      14341.A0339      14341.A0340 
    ##                1                1                1                1 
    ##      14341.A0341      14341.A0342      14341.A0343      14341.A0344 
    ##                1                1                1                1 
    ##      14341.A0345      14341.A0346      14341.A0347      14341.A0348 
    ##                1                1                1                1 
    ##      14341.A0349      14341.A0350      14341.A0351      14341.A0352 
    ##                1                1                1                1 
    ##      14341.A0353      14341.A0354      14341.A0355      14341.A0356 
    ##                1                1                1                1 
    ##      14341.A0357      14341.A0358      14341.A0359      14341.A0360 
    ##                1                1                1                1 
    ##      14341.A0361      14341.A0362      14341.A0363      14341.A0364 
    ##                1                1                1                1 
    ##      14341.A0365      14341.A0366      14341.A0367      14341.A0368 
    ##                1                1                1                1 
    ##      14341.A0369      14341.A0370      14341.A0371      14341.A0372 
    ##                1                1                1                1 
    ##      14341.A0373      14341.A0374      14341.A0375      14341.A0376 
    ##                1                1                1                1 
    ##      14341.A0377      14341.A0378      14341.A0379      14341.A0380 
    ##                1                1                1                1 
    ##      14341.A0381      14341.A0382      14341.A0383      14341.A0384 
    ##                1                1                1                1 
    ##      14341.A0385      14341.A0386      14341.A0387      14341.A0388 
    ##                1                1                1                1 
    ##      14341.A0389      14341.A0390      14341.A0391      14341.A0392 
    ##                1                1                1                1 
    ##      14341.A0393      14341.A0394      14341.A0395      14341.A0396 
    ##                1                1                1                1 
    ##      14341.A0397      14341.A0398      14341.A0399      14341.A0400 
    ##                1                1                1                1 
    ##      14341.A0401      14341.A0402      14341.A0403      14341.A0404 
    ##                1                1                1                1 
    ##      14341.A0405      14341.A0406      14341.A0407      14341.A0408 
    ##                1                1                1                1 
    ##      14341.A0409      14341.A0410      14341.A0411      14341.A0412 
    ##                1                1                1                1 
    ##      14341.A0413      14341.A0414      14341.A0415      14341.A0416 
    ##                1                1                1                1 
    ##      14341.A0417      14341.A0418      14341.A0419      14341.A0420 
    ##                1                1                1                1 
    ##      14341.A0421      14341.A0422      14341.A0423      14341.A0424 
    ##                1                1                1                1 
    ##      14341.A0425      14341.A0426      14341.A0427      14341.A0428 
    ##                1                1                1                1 
    ##      14341.A0429      14341.A0430      14341.A0431      14341.A0432 
    ##                1                1                1                1 
    ##      14341.A0433      14341.A0434      14341.A0435      14341.A0436 
    ##                1                1                1                1 
    ##      14341.A0437      14341.A0438      14341.A0439      14341.A0440 
    ##                1                1                1                1 
    ##      14341.A0441      14341.A0442      14341.A0443      14341.A0444 
    ##                1                1                1                1 
    ##      14341.A0445      14341.A0446      14341.A0447      14341.A0448 
    ##                1                1                1                1 
    ##      14341.A0449      14341.A0450      14341.A0451      14341.A0452 
    ##                1                1                1                1 
    ##      14341.A0453      14341.A0454      14341.A0455      14341.A0456 
    ##                1                1                1                1 
    ##      14341.A0457      14341.A0458      14341.A0459      14341.A0460 
    ##                1                1                1                1 
    ##      14341.A0461      14341.A0462      14341.A0463      14341.A0464 
    ##                1                1                1                1 
    ##      14341.A0465      14341.A0466      14341.A0467      14341.A0468 
    ##                1                1                1                1 
    ##      14341.A0469      14341.A0470      14341.A0471      14341.A0472 
    ##                1                1                1                1 
    ##      14341.A0473      14341.A0474      14341.A0475      14341.A0476 
    ##                1                1                1                1 
    ##      14341.A0477      14341.A0478      14341.A0479      14341.A0480 
    ##                1                1                1                1 
    ##      14341.A0481      14341.A0482      14341.A0483      14341.A0484 
    ##                1                1                1                1 
    ##      14341.A0485      14341.A0486      14341.A0487      14341.A0488 
    ##                1                1                1                1 
    ##      14341.A0489      14341.A0490      14341.A0491      14341.A0492 
    ##                1                1                1                1 
    ##      14341.A0493      14341.A0494      14341.A0495      14341.A0496 
    ##                1                1                1                1 
    ##      14341.A0497      14341.A0498      14341.A0499      14341.A0500 
    ##                1                1                1                1 
    ##      14341.A0501      14341.A0502      14341.A0503      14341.A0504 
    ##                1                1                1                1 
    ##      14341.A0505      14341.A0506      14341.A0507      14341.A0508 
    ##                1                1                1                1 
    ##      14341.A0509      14341.A0510      14341.A0511      14341.A0512 
    ##                1                1                1                1 
    ##      14341.A0513      14341.A0514      14341.A0515      14341.A0516 
    ##                1                1                1                1 
    ##      14341.A0517      14341.A0518      14341.A0519      14341.A0520 
    ##                1                1                1                1 
    ##      14341.A0521      14341.A0522      14341.A0523      14341.A0524 
    ##                1                1                1                1 
    ##      14341.A0525      14341.A0526      14341.A0527      14341.A0528 
    ##                1                1                1                1 
    ##      14341.A0529      14341.A0530      14341.A0531      14341.A0532 
    ##                1                1                1                1 
    ##      14341.A0533      14341.A0534      14341.A0535      14341.A0536 
    ##                1                1                1                1 
    ##      14341.A0537      14341.A0538      14341.A0539      14341.A0540 
    ##                1                1                1                1 
    ##      14341.A0541      14341.A0542      14341.A0543      14341.A0544 
    ##                1                1                1                1 
    ##      14341.A0545      14341.A0546      14341.A0547      14341.A0548 
    ##                1                1                1                1 
    ##      14341.A0549      14341.A0550      14341.A0551      14341.A0552 
    ##                1                1                1                1 
    ##      14341.A0553      14341.A0554      14341.A0555      14341.A0556 
    ##                1                1                1                1 
    ##      14341.A0557      14341.A0558      14341.A0559      14341.A0560 
    ##                1                1                1                1 
    ##      14341.A0561      14341.A0562      14341.A0563      14341.A0564 
    ##                1                1                1                1 
    ##      14341.A0565      14341.A0566      14341.A0567      14341.A0568 
    ##                1                1                1                1 
    ##      14341.A0569      14341.A0570      14341.A0571      14341.A0572 
    ##                1                1                1                1 
    ##      14341.A0573      14341.A0574      14341.A0575      14341.A0576 
    ##                1                1                1                1 
    ##      14341.A0577      14341.A0578      14341.A0579      14341.A0580 
    ##                1                1                1                1 
    ##      14341.A0581      14341.A0582      14341.A0583      14341.A0584 
    ##                1                1                1                1 
    ##      14341.A0585      14341.A0586      14341.A0587      14341.A0588 
    ##                1                1                1                1 
    ##      14341.A0589      14341.A0590      14341.A0591      14341.A0592 
    ##                1                1                1                1 
    ##      14341.A0593      14341.A0594      14341.A0595      14341.A0596 
    ##                1                1                1                1 
    ##      14341.A0597      14341.A0598      14341.A0599      14341.A0600 
    ##                1                1                1                1 
    ##      14341.A0601      14341.A0602      14341.A0603      14341.A0604 
    ##                1                1                1                1 
    ##      14341.A0605      14341.A0606      14341.A0607      14341.A0608 
    ##                1                1                1                1 
    ##      14341.A0609      14341.A0610      14341.A0611      14341.A0612 
    ##                1                1                1                1 
    ##      14341.A0613      14341.A0614      14341.A0615      14341.A0616 
    ##                1                1                1                1 
    ##      14341.A0617      14341.A0618      14341.A0619      14341.A0620 
    ##                1                1                1                1 
    ##      14341.A0621      14341.A0622      14341.A0623      14341.A0624 
    ##                1                1                1                1 
    ##      14341.A0625      14341.A0626      14341.A0627      14341.A0628 
    ##                1                1                1                1 
    ##      14341.A0629      14341.A0630      14341.A0631      14341.A0632 
    ##                1                1                1                1 
    ##      14341.A0633      14341.A0634      14341.A0635      14341.A0636 
    ##                1                1                1                1 
    ##      14341.A0637         14341.S1        14341.S10        14341.S11 
    ##                1                1                1                1 
    ##        14341.S12        14341.S13        14341.S14        14341.S15 
    ##                1                1                1                1 
    ##        14341.S16        14341.S17        14341.S18        14341.S19 
    ##                1                1                1                1 
    ##         14341.S2        14341.S20        14341.S21        14341.S22 
    ##                1                1                1                1 
    ##        14341.S23        14341.S24         14341.S3         14341.S4 
    ##                1                1                1                1 
    ##         14341.S5         14341.S6         14341.S7         14341.S8 
    ##                1                1                1                1 
    ##         14341.S9    14341.SAMPLE1   14341.SAMPLE10  14341.SAMPLE100 
    ##                1                1                1                1 
    ## 14341.SAMPLE1000 14341.SAMPLE1001 14341.SAMPLE1002 14341.SAMPLE1003 
    ##                1                1                1                1 
    ## 14341.SAMPLE1004 14341.SAMPLE1005 14341.SAMPLE1006 14341.SAMPLE1007 
    ##                1                1                1                1 
    ## 14341.SAMPLE1008 14341.SAMPLE1009  14341.SAMPLE101 14341.SAMPLE1010 
    ##                1                1                1                1 
    ## 14341.SAMPLE1011 14341.SAMPLE1012 14341.SAMPLE1013 14341.SAMPLE1014 
    ##                1                1                1                1 
    ## 14341.SAMPLE1015 14341.SAMPLE1016 14341.SAMPLE1017 14341.SAMPLE1018 
    ##                1                1                1                1 
    ## 14341.SAMPLE1019  14341.SAMPLE102 14341.SAMPLE1020 14341.SAMPLE1021 
    ##                1                1                1                1 
    ## 14341.SAMPLE1022 14341.SAMPLE1023 14341.SAMPLE1024 14341.SAMPLE1025 
    ##                1                1                1                1 
    ## 14341.SAMPLE1026 14341.SAMPLE1027 14341.SAMPLE1028 14341.SAMPLE1029 
    ##                1                1                1                1 
    ##  14341.SAMPLE103 14341.SAMPLE1030 14341.SAMPLE1031 14341.SAMPLE1032 
    ##                1                1                1                1 
    ## 14341.SAMPLE1033 14341.SAMPLE1034 14341.SAMPLE1035 14341.SAMPLE1036 
    ##                1                1                1                1 
    ## 14341.SAMPLE1037 14341.SAMPLE1038 14341.SAMPLE1039  14341.SAMPLE104 
    ##                1                1                1                1 
    ## 14341.SAMPLE1040 14341.SAMPLE1041 14341.SAMPLE1042 14341.SAMPLE1043 
    ##                1                1                1                1 
    ## 14341.SAMPLE1044 14341.SAMPLE1045 14341.SAMPLE1046 14341.SAMPLE1047 
    ##                1                1                1                1 
    ## 14341.SAMPLE1048 14341.SAMPLE1049  14341.SAMPLE105 14341.SAMPLE1050 
    ##                1                1                1                1 
    ## 14341.SAMPLE1051 14341.SAMPLE1052 14341.SAMPLE1053 14341.SAMPLE1054 
    ##                1                1                1                1 
    ## 14341.SAMPLE1055 14341.SAMPLE1056 14341.SAMPLE1057 14341.SAMPLE1058 
    ##                1                1                1                1 
    ## 14341.SAMPLE1059  14341.SAMPLE106 14341.SAMPLE1060 14341.SAMPLE1061 
    ##                1                1                1                1 
    ## 14341.SAMPLE1062 14341.SAMPLE1063 14341.SAMPLE1064 14341.SAMPLE1065 
    ##                1                1                1                1 
    ## 14341.SAMPLE1066 14341.SAMPLE1067 14341.SAMPLE1068 14341.SAMPLE1069 
    ##                1                1                1                1 
    ##  14341.SAMPLE107 14341.SAMPLE1070 14341.SAMPLE1071 14341.SAMPLE1072 
    ##                1                1                1                1 
    ## 14341.SAMPLE1073 14341.SAMPLE1074 14341.SAMPLE1075 14341.SAMPLE1076 
    ##                1                1                1                1 
    ## 14341.SAMPLE1077 14341.SAMPLE1078 14341.SAMPLE1079  14341.SAMPLE108 
    ##                1                1                1                1 
    ## 14341.SAMPLE1080 14341.SAMPLE1081 14341.SAMPLE1082 14341.SAMPLE1083 
    ##                1                1                1                1 
    ## 14341.SAMPLE1084 14341.SAMPLE1085 14341.SAMPLE1086 14341.SAMPLE1087 
    ##                1                1                1                1 
    ## 14341.SAMPLE1088 14341.SAMPLE1089  14341.SAMPLE109 14341.SAMPLE1090 
    ##                1                1                1                1 
    ## 14341.SAMPLE1091 14341.SAMPLE1092 14341.SAMPLE1093 14341.SAMPLE1094 
    ##                1                1                1                1 
    ## 14341.SAMPLE1095 14341.SAMPLE1096 14341.SAMPLE1097 14341.SAMPLE1098 
    ##                1                1                1                1 
    ## 14341.SAMPLE1099   14341.SAMPLE11  14341.SAMPLE110 14341.SAMPLE1100 
    ##                1                1                1                1 
    ## 14341.SAMPLE1101 14341.SAMPLE1102 14341.SAMPLE1103 14341.SAMPLE1104 
    ##                1                1                1                1 
    ## 14341.SAMPLE1105 14341.SAMPLE1106 14341.SAMPLE1107 14341.SAMPLE1108 
    ##                1                1                1                1 
    ## 14341.SAMPLE1109  14341.SAMPLE111 14341.SAMPLE1110 14341.SAMPLE1111 
    ##                1                1                1                1 
    ## 14341.SAMPLE1112 14341.SAMPLE1113 14341.SAMPLE1114 14341.SAMPLE1115 
    ##                1                1                1                1 
    ## 14341.SAMPLE1116 14341.SAMPLE1117 14341.SAMPLE1118 14341.SAMPLE1119 
    ##                1                1                1                1 
    ##  14341.SAMPLE112 14341.SAMPLE1120 14341.SAMPLE1121 14341.SAMPLE1122 
    ##                1                1                1                1 
    ## 14341.SAMPLE1123 14341.SAMPLE1124 14341.SAMPLE1125 14341.SAMPLE1126 
    ##                1                1                1                1 
    ## 14341.SAMPLE1127 14341.SAMPLE1128 14341.SAMPLE1129  14341.SAMPLE113 
    ##                1                1                1                1 
    ## 14341.SAMPLE1130 14341.SAMPLE1131 14341.SAMPLE1132 14341.SAMPLE1133 
    ##                1                1                1                1 
    ## 14341.SAMPLE1134 14341.SAMPLE1135 14341.SAMPLE1136 14341.SAMPLE1137 
    ##                1                1                1                1 
    ## 14341.SAMPLE1138 14341.SAMPLE1139  14341.SAMPLE114 14341.SAMPLE1140 
    ##                1                1                1                1 
    ## 14341.SAMPLE1141 14341.SAMPLE1142 14341.SAMPLE1143 14341.SAMPLE1144 
    ##                1                1                1                1 
    ## 14341.SAMPLE1145 14341.SAMPLE1146 14341.SAMPLE1147 14341.SAMPLE1148 
    ##                1                1                1                1 
    ## 14341.SAMPLE1149  14341.SAMPLE115 14341.SAMPLE1150 14341.SAMPLE1151 
    ##                1                1                1                1 
    ## 14341.SAMPLE1152 14341.SAMPLE1153 14341.SAMPLE1154 14341.SAMPLE1155 
    ##                1                1                1                1 
    ## 14341.SAMPLE1156 14341.SAMPLE1157 14341.SAMPLE1158 14341.SAMPLE1159 
    ##                1                1                1                1 
    ##  14341.SAMPLE116 14341.SAMPLE1160 14341.SAMPLE1161 14341.SAMPLE1162 
    ##                1                1                1                1 
    ## 14341.SAMPLE1163 14341.SAMPLE1164 14341.SAMPLE1165 14341.SAMPLE1166 
    ##                1                1                1                1 
    ## 14341.SAMPLE1167 14341.SAMPLE1168 14341.SAMPLE1169  14341.SAMPLE117 
    ##                1                1                1                1 
    ## 14341.SAMPLE1170 14341.SAMPLE1171 14341.SAMPLE1172 14341.SAMPLE1173 
    ##                1                1                1                1 
    ## 14341.SAMPLE1174 14341.SAMPLE1175 14341.SAMPLE1176 14341.SAMPLE1177 
    ##                1                1                1                1 
    ## 14341.SAMPLE1178 14341.SAMPLE1179  14341.SAMPLE118 14341.SAMPLE1180 
    ##                1                1                1                1 
    ## 14341.SAMPLE1181 14341.SAMPLE1182 14341.SAMPLE1183 14341.SAMPLE1184 
    ##                1                1                1                1 
    ## 14341.SAMPLE1185 14341.SAMPLE1186 14341.SAMPLE1187 14341.SAMPLE1188 
    ##                1                1                1                1 
    ## 14341.SAMPLE1189  14341.SAMPLE119 14341.SAMPLE1190 14341.SAMPLE1191 
    ##                1                1                1                1 
    ## 14341.SAMPLE1192 14341.SAMPLE1193 14341.SAMPLE1194 14341.SAMPLE1195 
    ##                1                1                1                1 
    ## 14341.SAMPLE1196 14341.SAMPLE1197 14341.SAMPLE1198 14341.SAMPLE1199 
    ##                1                1                1                1 
    ##   14341.SAMPLE12  14341.SAMPLE120 14341.SAMPLE1200 14341.SAMPLE1201 
    ##                1                1                1                1 
    ## 14341.SAMPLE1202 14341.SAMPLE1203 14341.SAMPLE1204 14341.SAMPLE1205 
    ##                1                1                1                1 
    ## 14341.SAMPLE1206 14341.SAMPLE1207 14341.SAMPLE1208 14341.SAMPLE1209 
    ##                1                1                1                1 
    ##  14341.SAMPLE121 14341.SAMPLE1210 14341.SAMPLE1211 14341.SAMPLE1212 
    ##                1                1                1                1 
    ## 14341.SAMPLE1213 14341.SAMPLE1214 14341.SAMPLE1215 14341.SAMPLE1216 
    ##                1                1                1                1 
    ## 14341.SAMPLE1217 14341.SAMPLE1218 14341.SAMPLE1219  14341.SAMPLE122 
    ##                1                1                1                1 
    ## 14341.SAMPLE1220 14341.SAMPLE1221 14341.SAMPLE1222 14341.SAMPLE1223 
    ##                1                1                1                1 
    ## 14341.SAMPLE1224 14341.SAMPLE1225 14341.SAMPLE1226 14341.SAMPLE1227 
    ##                1                1                1                1 
    ## 14341.SAMPLE1228 14341.SAMPLE1229  14341.SAMPLE123 14341.SAMPLE1230 
    ##                1                1                1                1 
    ## 14341.SAMPLE1231 14341.SAMPLE1232 14341.SAMPLE1233 14341.SAMPLE1234 
    ##                1                1                1                1 
    ## 14341.SAMPLE1235 14341.SAMPLE1236 14341.SAMPLE1237 14341.SAMPLE1238 
    ##                1                1                1                1 
    ## 14341.SAMPLE1239  14341.SAMPLE124 14341.SAMPLE1240 14341.SAMPLE1241 
    ##                1                1                1                1 
    ## 14341.SAMPLE1242 14341.SAMPLE1243 14341.SAMPLE1244 14341.SAMPLE1245 
    ##                1                1                1                1 
    ## 14341.SAMPLE1246 14341.SAMPLE1247 14341.SAMPLE1248 14341.SAMPLE1249 
    ##                1                1                1                1 
    ##  14341.SAMPLE125 14341.SAMPLE1250 14341.SAMPLE1251 14341.SAMPLE1252 
    ##                1                1                1                1 
    ## 14341.SAMPLE1253 14341.SAMPLE1254 14341.SAMPLE1255 14341.SAMPLE1256 
    ##                1                1                1                1 
    ## 14341.SAMPLE1257 14341.SAMPLE1258 14341.SAMPLE1259  14341.SAMPLE126 
    ##                1                1                1                1 
    ## 14341.SAMPLE1260 14341.SAMPLE1261 14341.SAMPLE1262 14341.SAMPLE1263 
    ##                1                1                1                1 
    ## 14341.SAMPLE1264 14341.SAMPLE1265 14341.SAMPLE1266 14341.SAMPLE1267 
    ##                1                1                1                1 
    ## 14341.SAMPLE1268 14341.SAMPLE1269  14341.SAMPLE127 14341.SAMPLE1270 
    ##                1                1                1                1 
    ## 14341.SAMPLE1271 14341.SAMPLE1272 14341.SAMPLE1273 14341.SAMPLE1274 
    ##                1                1                1                1 
    ## 14341.SAMPLE1275 14341.SAMPLE1276 14341.SAMPLE1277 14341.SAMPLE1278 
    ##                1                1                1                1 
    ## 14341.SAMPLE1279  14341.SAMPLE128 14341.SAMPLE1280 14341.SAMPLE1281 
    ##                1                1                1                1 
    ## 14341.SAMPLE1282 14341.SAMPLE1283 14341.SAMPLE1284 14341.SAMPLE1285 
    ##                1                1                1                1 
    ## 14341.SAMPLE1286 14341.SAMPLE1287 14341.SAMPLE1288 14341.SAMPLE1289 
    ##                1                1                1                1 
    ##  14341.SAMPLE129 14341.SAMPLE1290 14341.SAMPLE1291 14341.SAMPLE1292 
    ##                1                1                1                1 
    ## 14341.SAMPLE1293 14341.SAMPLE1294 14341.SAMPLE1295 14341.SAMPLE1296 
    ##                1                1                1                1 
    ## 14341.SAMPLE1297 14341.SAMPLE1298 14341.SAMPLE1299   14341.SAMPLE13 
    ##                1                1                1                1 
    ##  14341.SAMPLE130 14341.SAMPLE1300 14341.SAMPLE1301 14341.SAMPLE1302 
    ##                1                1                1                1 
    ## 14341.SAMPLE1303 14341.SAMPLE1304 14341.SAMPLE1305 14341.SAMPLE1306 
    ##                1                1                1                1 
    ## 14341.SAMPLE1307 14341.SAMPLE1308 14341.SAMPLE1309  14341.SAMPLE131 
    ##                1                1                1                1 
    ## 14341.SAMPLE1310 14341.SAMPLE1311 14341.SAMPLE1312 14341.SAMPLE1313 
    ##                1                1                1                1 
    ## 14341.SAMPLE1314 14341.SAMPLE1315 14341.SAMPLE1316 14341.SAMPLE1317 
    ##                1                1                1                1 
    ## 14341.SAMPLE1318 14341.SAMPLE1319  14341.SAMPLE132 14341.SAMPLE1320 
    ##                1                1                1                1 
    ## 14341.SAMPLE1321 14341.SAMPLE1322 14341.SAMPLE1323 14341.SAMPLE1324 
    ##                1                1                1                1 
    ## 14341.SAMPLE1325 14341.SAMPLE1326 14341.SAMPLE1327 14341.SAMPLE1328 
    ##                1                1                1                1 
    ## 14341.SAMPLE1329  14341.SAMPLE133 14341.SAMPLE1330 14341.SAMPLE1331 
    ##                1                1                1                1 
    ## 14341.SAMPLE1332 14341.SAMPLE1333 14341.SAMPLE1334 14341.SAMPLE1335 
    ##                1                1                1                1 
    ## 14341.SAMPLE1336 14341.SAMPLE1337 14341.SAMPLE1338 14341.SAMPLE1339 
    ##                1                1                1                1 
    ##  14341.SAMPLE134 14341.SAMPLE1340 14341.SAMPLE1341 14341.SAMPLE1342 
    ##                1                1                1                1 
    ## 14341.SAMPLE1343 14341.SAMPLE1344 14341.SAMPLE1345 14341.SAMPLE1346 
    ##                1                1                1                1 
    ## 14341.SAMPLE1347 14341.SAMPLE1348 14341.SAMPLE1349  14341.SAMPLE135 
    ##                1                1                1                1 
    ## 14341.SAMPLE1350 14341.SAMPLE1351 14341.SAMPLE1352 14341.SAMPLE1353 
    ##                1                1                1                1 
    ## 14341.SAMPLE1354 14341.SAMPLE1355 14341.SAMPLE1356 14341.SAMPLE1357 
    ##                1                1                1                1 
    ## 14341.SAMPLE1358 14341.SAMPLE1359  14341.SAMPLE136 14341.SAMPLE1360 
    ##                1                1                1                1 
    ## 14341.SAMPLE1361 14341.SAMPLE1362 14341.SAMPLE1363 14341.SAMPLE1364 
    ##                1                1                1                1 
    ## 14341.SAMPLE1365 14341.SAMPLE1366 14341.SAMPLE1367 14341.SAMPLE1368 
    ##                1                1                1                1 
    ## 14341.SAMPLE1369  14341.SAMPLE137 14341.SAMPLE1370 14341.SAMPLE1371 
    ##                1                1                1                1 
    ## 14341.SAMPLE1372 14341.SAMPLE1373 14341.SAMPLE1374 14341.SAMPLE1375 
    ##                1                1                1                1 
    ## 14341.SAMPLE1376 14341.SAMPLE1377 14341.SAMPLE1378 14341.SAMPLE1379 
    ##                1                1                1                1 
    ##  14341.SAMPLE138 14341.SAMPLE1380 14341.SAMPLE1381 14341.SAMPLE1382 
    ##                1                1                1                1 
    ## 14341.SAMPLE1383 14341.SAMPLE1384 14341.SAMPLE1385 14341.SAMPLE1386 
    ##                1                1                1                1 
    ## 14341.SAMPLE1387 14341.SAMPLE1388 14341.SAMPLE1389  14341.SAMPLE139 
    ##                1                1                1                1 
    ## 14341.SAMPLE1390 14341.SAMPLE1391 14341.SAMPLE1392 14341.SAMPLE1393 
    ##                1                1                1                1 
    ## 14341.SAMPLE1394 14341.SAMPLE1395 14341.SAMPLE1396 14341.SAMPLE1397 
    ##                1                1                1                1 
    ## 14341.SAMPLE1398 14341.SAMPLE1399   14341.SAMPLE14  14341.SAMPLE140 
    ##                1                1                1                1 
    ## 14341.SAMPLE1400 14341.SAMPLE1401 14341.SAMPLE1402 14341.SAMPLE1403 
    ##                1                1                1                1 
    ## 14341.SAMPLE1404 14341.SAMPLE1405 14341.SAMPLE1406 14341.SAMPLE1407 
    ##                1                1                1                1 
    ## 14341.SAMPLE1408 14341.SAMPLE1409  14341.SAMPLE141 14341.SAMPLE1410 
    ##                1                1                1                1 
    ## 14341.SAMPLE1411 14341.SAMPLE1412 14341.SAMPLE1413 14341.SAMPLE1414 
    ##                1                1                1                1 
    ## 14341.SAMPLE1415 14341.SAMPLE1416 14341.SAMPLE1417 14341.SAMPLE1418 
    ##                1                1                1                1 
    ## 14341.SAMPLE1419  14341.SAMPLE142 14341.SAMPLE1420 14341.SAMPLE1421 
    ##                1                1                1                1 
    ## 14341.SAMPLE1422 14341.SAMPLE1423 14341.SAMPLE1424 14341.SAMPLE1425 
    ##                1                1                1                1 
    ## 14341.SAMPLE1426 14341.SAMPLE1427 14341.SAMPLE1428 14341.SAMPLE1429 
    ##                1                1                1                1 
    ##  14341.SAMPLE143 14341.SAMPLE1430 14341.SAMPLE1431 14341.SAMPLE1432 
    ##                1                1                1                1 
    ## 14341.SAMPLE1433 14341.SAMPLE1434 14341.SAMPLE1435 14341.SAMPLE1436 
    ##                1                1                1                1 
    ## 14341.SAMPLE1437 14341.SAMPLE1438 14341.SAMPLE1439  14341.SAMPLE144 
    ##                1                1                1                1 
    ## 14341.SAMPLE1440 14341.SAMPLE1441 14341.SAMPLE1442 14341.SAMPLE1443 
    ##                1                1                1                1 
    ## 14341.SAMPLE1444 14341.SAMPLE1445 14341.SAMPLE1446 14341.SAMPLE1447 
    ##                1                1                1                1 
    ## 14341.SAMPLE1448 14341.SAMPLE1449  14341.SAMPLE145 14341.SAMPLE1450 
    ##                1                1                1                1 
    ## 14341.SAMPLE1451 14341.SAMPLE1452 14341.SAMPLE1453 14341.SAMPLE1454 
    ##                1                1                1                1 
    ## 14341.SAMPLE1455 14341.SAMPLE1456 14341.SAMPLE1457 14341.SAMPLE1458 
    ##                1                1                1                1 
    ## 14341.SAMPLE1459  14341.SAMPLE146 14341.SAMPLE1460 14341.SAMPLE1461 
    ##                1                1                1                1 
    ## 14341.SAMPLE1462 14341.SAMPLE1463 14341.SAMPLE1464 14341.SAMPLE1465 
    ##                1                1                1                1 
    ## 14341.SAMPLE1466 14341.SAMPLE1467 14341.SAMPLE1468 14341.SAMPLE1469 
    ##                1                1                1                1 
    ##  14341.SAMPLE147 14341.SAMPLE1470 14341.SAMPLE1471 14341.SAMPLE1472 
    ##                1                1                1                1 
    ## 14341.SAMPLE1473 14341.SAMPLE1474 14341.SAMPLE1475 14341.SAMPLE1476 
    ##                1                1                1                1 
    ## 14341.SAMPLE1477 14341.SAMPLE1478 14341.SAMPLE1479  14341.SAMPLE148 
    ##                1                1                1                1 
    ## 14341.SAMPLE1480 14341.SAMPLE1481 14341.SAMPLE1482 14341.SAMPLE1483 
    ##                1                1                1                1 
    ## 14341.SAMPLE1484 14341.SAMPLE1485 14341.SAMPLE1486 14341.SAMPLE1487 
    ##                1                1                1                1 
    ## 14341.SAMPLE1488 14341.SAMPLE1489  14341.SAMPLE149 14341.SAMPLE1490 
    ##                1                1                1                1 
    ## 14341.SAMPLE1491 14341.SAMPLE1492 14341.SAMPLE1493 14341.SAMPLE1494 
    ##                1                1                1                1 
    ## 14341.SAMPLE1495 14341.SAMPLE1496 14341.SAMPLE1497 14341.SAMPLE1498 
    ##                1                1                1                1 
    ## 14341.SAMPLE1499   14341.SAMPLE15  14341.SAMPLE150 14341.SAMPLE1500 
    ##                1                1                1                1 
    ## 14341.SAMPLE1501 14341.SAMPLE1502 14341.SAMPLE1503 14341.SAMPLE1504 
    ##                1                1                1                1 
    ## 14341.SAMPLE1505 14341.SAMPLE1506 14341.SAMPLE1507 14341.SAMPLE1508 
    ##                1                1                1                1 
    ## 14341.SAMPLE1509  14341.SAMPLE151 14341.SAMPLE1510 14341.SAMPLE1511 
    ##                1                1                1                1 
    ## 14341.SAMPLE1512 14341.SAMPLE1513 14341.SAMPLE1514 14341.SAMPLE1515 
    ##                1                1                1                1 
    ## 14341.SAMPLE1516 14341.SAMPLE1517 14341.SAMPLE1518 14341.SAMPLE1519 
    ##                1                1                1                1 
    ##  14341.SAMPLE152 14341.SAMPLE1520 14341.SAMPLE1521 14341.SAMPLE1522 
    ##                1                1                1                1 
    ## 14341.SAMPLE1523 14341.SAMPLE1524 14341.SAMPLE1525 14341.SAMPLE1526 
    ##                1                1                1                1 
    ## 14341.SAMPLE1527 14341.SAMPLE1528 14341.SAMPLE1529  14341.SAMPLE153 
    ##                1                1                1                1 
    ## 14341.SAMPLE1530 14341.SAMPLE1531 14341.SAMPLE1532 14341.SAMPLE1533 
    ##                1                1                1                1 
    ## 14341.SAMPLE1534 14341.SAMPLE1535 14341.SAMPLE1536 14341.SAMPLE1537 
    ##                1                1                1                1 
    ## 14341.SAMPLE1538 14341.SAMPLE1539  14341.SAMPLE154 14341.SAMPLE1540 
    ##                1                1                1                1 
    ## 14341.SAMPLE1541 14341.SAMPLE1542 14341.SAMPLE1543 14341.SAMPLE1544 
    ##                1                1                1                1 
    ## 14341.SAMPLE1545 14341.SAMPLE1546 14341.SAMPLE1547 14341.SAMPLE1548 
    ##                1                1                1                1 
    ## 14341.SAMPLE1549  14341.SAMPLE155 14341.SAMPLE1550 14341.SAMPLE1551 
    ##                1                1                1                1 
    ## 14341.SAMPLE1552 14341.SAMPLE1553 14341.SAMPLE1554 14341.SAMPLE1555 
    ##                1                1                1                1 
    ## 14341.SAMPLE1556 14341.SAMPLE1557 14341.SAMPLE1558 14341.SAMPLE1559 
    ##                1                1                1                1 
    ##  14341.SAMPLE156 14341.SAMPLE1560 14341.SAMPLE1561 14341.SAMPLE1562 
    ##                1                1                1                1 
    ## 14341.SAMPLE1563 14341.SAMPLE1564 14341.SAMPLE1565 14341.SAMPLE1566 
    ##                1                1                1                1 
    ## 14341.SAMPLE1567 14341.SAMPLE1568 14341.SAMPLE1569  14341.SAMPLE157 
    ##                1                1                1                1 
    ## 14341.SAMPLE1570 14341.SAMPLE1571 14341.SAMPLE1572 14341.SAMPLE1573 
    ##                1                1                1                1 
    ## 14341.SAMPLE1574 14341.SAMPLE1575 14341.SAMPLE1576 14341.SAMPLE1577 
    ##                1                1                1                1 
    ## 14341.SAMPLE1578 14341.SAMPLE1579  14341.SAMPLE158 14341.SAMPLE1580 
    ##                1                1                1                1 
    ## 14341.SAMPLE1581 14341.SAMPLE1582 14341.SAMPLE1583 14341.SAMPLE1584 
    ##                1                1                1                1 
    ## 14341.SAMPLE1585 14341.SAMPLE1586 14341.SAMPLE1587 14341.SAMPLE1588 
    ##                1                1                1                1 
    ## 14341.SAMPLE1589  14341.SAMPLE159 14341.SAMPLE1590 14341.SAMPLE1591 
    ##                1                1                1                1 
    ## 14341.SAMPLE1592 14341.SAMPLE1593 14341.SAMPLE1594 14341.SAMPLE1595 
    ##                1                1                1                1 
    ## 14341.SAMPLE1596 14341.SAMPLE1597 14341.SAMPLE1598 14341.SAMPLE1599 
    ##                1                1                1                1 
    ##   14341.SAMPLE16  14341.SAMPLE160 14341.SAMPLE1600 14341.SAMPLE1601 
    ##                1                1                1                1 
    ## 14341.SAMPLE1602 14341.SAMPLE1603 14341.SAMPLE1604 14341.SAMPLE1605 
    ##                1                1                1                1 
    ## 14341.SAMPLE1606 14341.SAMPLE1607 14341.SAMPLE1608 14341.SAMPLE1609 
    ##                1                1                1                1 
    ##  14341.SAMPLE161 14341.SAMPLE1610 14341.SAMPLE1611 14341.SAMPLE1612 
    ##                1                1                1                1 
    ## 14341.SAMPLE1613 14341.SAMPLE1614 14341.SAMPLE1615 14341.SAMPLE1616 
    ##                1                1                1                1 
    ## 14341.SAMPLE1617 14341.SAMPLE1618 14341.SAMPLE1619  14341.SAMPLE162 
    ##                1                1                1                1 
    ## 14341.SAMPLE1620 14341.SAMPLE1621 14341.SAMPLE1622 14341.SAMPLE1623 
    ##                1                1                1                1 
    ## 14341.SAMPLE1624 14341.SAMPLE1625 14341.SAMPLE1626 14341.SAMPLE1627 
    ##                1                1                1                1 
    ## 14341.SAMPLE1628 14341.SAMPLE1629  14341.SAMPLE163 14341.SAMPLE1630 
    ##                1                1                1                1 
    ## 14341.SAMPLE1631 14341.SAMPLE1632 14341.SAMPLE1633 14341.SAMPLE1634 
    ##                1                1                1                1 
    ## 14341.SAMPLE1635 14341.SAMPLE1636 14341.SAMPLE1637 14341.SAMPLE1638 
    ##                1                1                1                1 
    ## 14341.SAMPLE1639  14341.SAMPLE164 14341.SAMPLE1640 14341.SAMPLE1641 
    ##                1                1                1                1 
    ## 14341.SAMPLE1642 14341.SAMPLE1643 14341.SAMPLE1644 14341.SAMPLE1645 
    ##                1                1                1                1 
    ## 14341.SAMPLE1646 14341.SAMPLE1647 14341.SAMPLE1648 14341.SAMPLE1649 
    ##                1                1                1                1 
    ##  14341.SAMPLE165 14341.SAMPLE1650 14341.SAMPLE1651 14341.SAMPLE1652 
    ##                1                1                1                1 
    ## 14341.SAMPLE1653 14341.SAMPLE1654 14341.SAMPLE1655 14341.SAMPLE1656 
    ##                1                1                1                1 
    ## 14341.SAMPLE1657 14341.SAMPLE1658 14341.SAMPLE1659  14341.SAMPLE166 
    ##                1                1                1                1 
    ## 14341.SAMPLE1660 14341.SAMPLE1661 14341.SAMPLE1662 14341.SAMPLE1663 
    ##                1                1                1                1 
    ## 14341.SAMPLE1664 14341.SAMPLE1665 14341.SAMPLE1666 14341.SAMPLE1667 
    ##                1                1                1                1 
    ## 14341.SAMPLE1668 14341.SAMPLE1669  14341.SAMPLE167 14341.SAMPLE1670 
    ##                1                1                1                1 
    ## 14341.SAMPLE1671 14341.SAMPLE1672 14341.SAMPLE1673 14341.SAMPLE1674 
    ##                1                1                1                1 
    ## 14341.SAMPLE1675 14341.SAMPLE1676 14341.SAMPLE1677 14341.SAMPLE1678 
    ##                1                1                1                1 
    ## 14341.SAMPLE1679  14341.SAMPLE168 14341.SAMPLE1680 14341.SAMPLE1681 
    ##                1                1                1                1 
    ## 14341.SAMPLE1682 14341.SAMPLE1683 14341.SAMPLE1684 14341.SAMPLE1685 
    ##                1                1                1                1 
    ## 14341.SAMPLE1686 14341.SAMPLE1687 14341.SAMPLE1688 14341.SAMPLE1689 
    ##                1                1                1                1 
    ##  14341.SAMPLE169 14341.SAMPLE1690 14341.SAMPLE1691 14341.SAMPLE1692 
    ##                1                1                1                1 
    ## 14341.SAMPLE1693 14341.SAMPLE1694 14341.SAMPLE1695 14341.SAMPLE1696 
    ##                1                1                1                1 
    ## 14341.SAMPLE1697 14341.SAMPLE1698 14341.SAMPLE1699   14341.SAMPLE17 
    ##                1                1                1                1 
    ##  14341.SAMPLE170 14341.SAMPLE1700 14341.SAMPLE1701 14341.SAMPLE1702 
    ##                1                1                1                1 
    ## 14341.SAMPLE1703 14341.SAMPLE1704 14341.SAMPLE1705 14341.SAMPLE1706 
    ##                1                1                1                1 
    ## 14341.SAMPLE1707 14341.SAMPLE1708 14341.SAMPLE1709  14341.SAMPLE171 
    ##                1                1                1                1 
    ## 14341.SAMPLE1710 14341.SAMPLE1711 14341.SAMPLE1712 14341.SAMPLE1713 
    ##                1                1                1                1 
    ## 14341.SAMPLE1714 14341.SAMPLE1715 14341.SAMPLE1716 14341.SAMPLE1717 
    ##                1                1                1                1 
    ## 14341.SAMPLE1718 14341.SAMPLE1719  14341.SAMPLE172 14341.SAMPLE1720 
    ##                1                1                1                1 
    ## 14341.SAMPLE1721 14341.SAMPLE1722 14341.SAMPLE1723 14341.SAMPLE1724 
    ##                1                1                1                1 
    ## 14341.SAMPLE1725 14341.SAMPLE1726 14341.SAMPLE1727 14341.SAMPLE1728 
    ##                1                1                1                1 
    ## 14341.SAMPLE1729  14341.SAMPLE173 14341.SAMPLE1730 14341.SAMPLE1731 
    ##                1                1                1                1 
    ## 14341.SAMPLE1732 14341.SAMPLE1733 14341.SAMPLE1734 14341.SAMPLE1735 
    ##                1                1                1                1 
    ## 14341.SAMPLE1736 14341.SAMPLE1737 14341.SAMPLE1738 14341.SAMPLE1739 
    ##                1                1                1                1 
    ##  14341.SAMPLE174 14341.SAMPLE1740 14341.SAMPLE1741 14341.SAMPLE1742 
    ##                1                1                1                1 
    ## 14341.SAMPLE1743 14341.SAMPLE1744 14341.SAMPLE1745 14341.SAMPLE1746 
    ##                1                1                1                1 
    ## 14341.SAMPLE1747 14341.SAMPLE1748 14341.SAMPLE1749  14341.SAMPLE175 
    ##                1                1                1                1 
    ## 14341.SAMPLE1750 14341.SAMPLE1751 14341.SAMPLE1752 14341.SAMPLE1753 
    ##                1                1                1                1 
    ## 14341.SAMPLE1754 14341.SAMPLE1755 14341.SAMPLE1756 14341.SAMPLE1757 
    ##                1                1                1                1 
    ## 14341.SAMPLE1758 14341.SAMPLE1759  14341.SAMPLE176 14341.SAMPLE1760 
    ##                1                1                1                1 
    ## 14341.SAMPLE1761 14341.SAMPLE1762 14341.SAMPLE1763 14341.SAMPLE1764 
    ##                1                1                1                1 
    ## 14341.SAMPLE1765 14341.SAMPLE1766 14341.SAMPLE1767 14341.SAMPLE1768 
    ##                1                1                1                1 
    ## 14341.SAMPLE1769  14341.SAMPLE177 14341.SAMPLE1770 14341.SAMPLE1771 
    ##                1                1                1                1 
    ## 14341.SAMPLE1772 14341.SAMPLE1773 14341.SAMPLE1774 14341.SAMPLE1775 
    ##                1                1                1                1 
    ## 14341.SAMPLE1776 14341.SAMPLE1777 14341.SAMPLE1778 14341.SAMPLE1779 
    ##                1                1                1                1 
    ##  14341.SAMPLE178 14341.SAMPLE1780 14341.SAMPLE1781 14341.SAMPLE1782 
    ##                1                1                1                1 
    ## 14341.SAMPLE1783 14341.SAMPLE1784 14341.SAMPLE1785 14341.SAMPLE1786 
    ##                1                1                1                1 
    ## 14341.SAMPLE1787 14341.SAMPLE1788 14341.SAMPLE1789  14341.SAMPLE179 
    ##                1                1                1                1 
    ## 14341.SAMPLE1790 14341.SAMPLE1791 14341.SAMPLE1792 14341.SAMPLE1793 
    ##                1                1                1                1 
    ## 14341.SAMPLE1794 14341.SAMPLE1795 14341.SAMPLE1796 14341.SAMPLE1797 
    ##                1                1                1                1 
    ## 14341.SAMPLE1798 14341.SAMPLE1799   14341.SAMPLE18  14341.SAMPLE180 
    ##                1                1                1                1 
    ## 14341.SAMPLE1800 14341.SAMPLE1801 14341.SAMPLE1802 14341.SAMPLE1803 
    ##                1                1                1                1 
    ## 14341.SAMPLE1804 14341.SAMPLE1805 14341.SAMPLE1806 14341.SAMPLE1807 
    ##                1                1                1                1 
    ## 14341.SAMPLE1808 14341.SAMPLE1809  14341.SAMPLE181 14341.SAMPLE1810 
    ##                1                1                1                1 
    ## 14341.SAMPLE1811 14341.SAMPLE1812 14341.SAMPLE1813 14341.SAMPLE1814 
    ##                1                1                1                1 
    ## 14341.SAMPLE1815 14341.SAMPLE1816 14341.SAMPLE1817 14341.SAMPLE1818 
    ##                1                1                1                1 
    ## 14341.SAMPLE1819  14341.SAMPLE182 14341.SAMPLE1820 14341.SAMPLE1821 
    ##                1                1                1                1 
    ## 14341.SAMPLE1822 14341.SAMPLE1823 14341.SAMPLE1824 14341.SAMPLE1825 
    ##                1                1                1                1 
    ## 14341.SAMPLE1826 14341.SAMPLE1827 14341.SAMPLE1828 14341.SAMPLE1829 
    ##                1                1                1                1 
    ##  14341.SAMPLE183 14341.SAMPLE1830 14341.SAMPLE1831 14341.SAMPLE1832 
    ##                1                1                1                1 
    ## 14341.SAMPLE1833 14341.SAMPLE1834 14341.SAMPLE1835 14341.SAMPLE1836 
    ##                1                1                1                1 
    ## 14341.SAMPLE1837 14341.SAMPLE1838 14341.SAMPLE1839  14341.SAMPLE184 
    ##                1                1                1                1 
    ## 14341.SAMPLE1840 14341.SAMPLE1841 14341.SAMPLE1842 14341.SAMPLE1843 
    ##                1                1                1                1 
    ## 14341.SAMPLE1844 14341.SAMPLE1845 14341.SAMPLE1846 14341.SAMPLE1847 
    ##                1                1                1                1 
    ## 14341.SAMPLE1848 14341.SAMPLE1849  14341.SAMPLE185 14341.SAMPLE1850 
    ##                1                1                1                1 
    ## 14341.SAMPLE1851 14341.SAMPLE1852 14341.SAMPLE1853 14341.SAMPLE1854 
    ##                1                1                1                1 
    ## 14341.SAMPLE1855 14341.SAMPLE1856 14341.SAMPLE1857 14341.SAMPLE1858 
    ##                1                1                1                1 
    ## 14341.SAMPLE1859  14341.SAMPLE186 14341.SAMPLE1860 14341.SAMPLE1861 
    ##                1                1                1                1 
    ## 14341.SAMPLE1862 14341.SAMPLE1863 14341.SAMPLE1864 14341.SAMPLE1865 
    ##                1                1                1                1 
    ## 14341.SAMPLE1866 14341.SAMPLE1867 14341.SAMPLE1868 14341.SAMPLE1869 
    ##                1                1                1                1 
    ##  14341.SAMPLE187 14341.SAMPLE1870 14341.SAMPLE1871 14341.SAMPLE1872 
    ##                1                1                1                1 
    ##  14341.SAMPLE188  14341.SAMPLE189   14341.SAMPLE19  14341.SAMPLE190 
    ##                1                1                1                1 
    ##  14341.SAMPLE191  14341.SAMPLE192  14341.SAMPLE193  14341.SAMPLE194 
    ##                1                1                1                1 
    ##  14341.SAMPLE195  14341.SAMPLE196  14341.SAMPLE197  14341.SAMPLE198 
    ##                1                1                1                1 
    ##  14341.SAMPLE199    14341.SAMPLE2   14341.SAMPLE20  14341.SAMPLE200 
    ##                1                1                1                1 
    ##  14341.SAMPLE201  14341.SAMPLE202  14341.SAMPLE203  14341.SAMPLE204 
    ##                1                1                1                1 
    ##  14341.SAMPLE205  14341.SAMPLE206  14341.SAMPLE207  14341.SAMPLE208 
    ##                1                1                1                1 
    ##  14341.SAMPLE209   14341.SAMPLE21  14341.SAMPLE210  14341.SAMPLE211 
    ##                1                1                1                1 
    ##  14341.SAMPLE212  14341.SAMPLE213  14341.SAMPLE214  14341.SAMPLE215 
    ##                1                1                1                1 
    ##  14341.SAMPLE216  14341.SAMPLE217  14341.SAMPLE218  14341.SAMPLE219 
    ##                1                1                1                1 
    ##   14341.SAMPLE22  14341.SAMPLE220  14341.SAMPLE221  14341.SAMPLE222 
    ##                1                1                1                1 
    ##  14341.SAMPLE223  14341.SAMPLE224  14341.SAMPLE225  14341.SAMPLE226 
    ##                1                1                1                1 
    ##  14341.SAMPLE227  14341.SAMPLE228  14341.SAMPLE229   14341.SAMPLE23 
    ##                1                1                1                1 
    ##  14341.SAMPLE230  14341.SAMPLE231  14341.SAMPLE232  14341.SAMPLE233 
    ##                1                1                1                1 
    ##  14341.SAMPLE234  14341.SAMPLE235  14341.SAMPLE236  14341.SAMPLE237 
    ##                1                1                1                1 
    ##  14341.SAMPLE238  14341.SAMPLE239   14341.SAMPLE24  14341.SAMPLE240 
    ##                1                1                1                1 
    ##  14341.SAMPLE241  14341.SAMPLE242  14341.SAMPLE243  14341.SAMPLE244 
    ##                1                1                1                1 
    ##  14341.SAMPLE245  14341.SAMPLE246  14341.SAMPLE247  14341.SAMPLE248 
    ##                1                1                1                1 
    ##  14341.SAMPLE249   14341.SAMPLE25  14341.SAMPLE250  14341.SAMPLE251 
    ##                1                1                1                1 
    ##  14341.SAMPLE252  14341.SAMPLE253  14341.SAMPLE254  14341.SAMPLE255 
    ##                1                1                1                1 
    ##  14341.SAMPLE256  14341.SAMPLE257  14341.SAMPLE258  14341.SAMPLE259 
    ##                1                1                1                1 
    ##   14341.SAMPLE26  14341.SAMPLE260  14341.SAMPLE261  14341.SAMPLE262 
    ##                1                1                1                1 
    ##  14341.SAMPLE263  14341.SAMPLE264  14341.SAMPLE265  14341.SAMPLE266 
    ##                1                1                1                1 
    ##  14341.SAMPLE267  14341.SAMPLE268  14341.SAMPLE269   14341.SAMPLE27 
    ##                1                1                1                1 
    ##  14341.SAMPLE270  14341.SAMPLE271  14341.SAMPLE272  14341.SAMPLE273 
    ##                1                1                1                1 
    ##  14341.SAMPLE274  14341.SAMPLE275  14341.SAMPLE276  14341.SAMPLE277 
    ##                1                1                1                1 
    ##  14341.SAMPLE278  14341.SAMPLE279   14341.SAMPLE28  14341.SAMPLE280 
    ##                1                1                1                1 
    ##  14341.SAMPLE281  14341.SAMPLE282  14341.SAMPLE283  14341.SAMPLE284 
    ##                1                1                1                1 
    ##  14341.SAMPLE285  14341.SAMPLE286  14341.SAMPLE287  14341.SAMPLE288 
    ##                1                1                1                1 
    ##  14341.SAMPLE289   14341.SAMPLE29  14341.SAMPLE290  14341.SAMPLE291 
    ##                1                1                1                1 
    ##  14341.SAMPLE292  14341.SAMPLE293  14341.SAMPLE294  14341.SAMPLE295 
    ##                1                1                1                1 
    ##  14341.SAMPLE296  14341.SAMPLE297  14341.SAMPLE298  14341.SAMPLE299 
    ##                1                1                1                1 
    ##    14341.SAMPLE3   14341.SAMPLE30  14341.SAMPLE300  14341.SAMPLE301 
    ##                1                1                1                1 
    ##  14341.SAMPLE302  14341.SAMPLE303  14341.SAMPLE304  14341.SAMPLE305 
    ##                1                1                1                1 
    ##  14341.SAMPLE306  14341.SAMPLE307  14341.SAMPLE308  14341.SAMPLE309 
    ##                1                1                1                1 
    ##   14341.SAMPLE31  14341.SAMPLE310  14341.SAMPLE311  14341.SAMPLE312 
    ##                1                1                1                1 
    ##  14341.SAMPLE313  14341.SAMPLE314  14341.SAMPLE315  14341.SAMPLE316 
    ##                1                1                1                1 
    ##  14341.SAMPLE318  14341.SAMPLE319   14341.SAMPLE32  14341.SAMPLE321 
    ##                1                1                1                1 
    ##  14341.SAMPLE322  14341.SAMPLE324  14341.SAMPLE325  14341.SAMPLE327 
    ##                1                1                1                1 
    ##  14341.SAMPLE329   14341.SAMPLE33  14341.SAMPLE330  14341.SAMPLE331 
    ##                1                1                1                1 
    ##  14341.SAMPLE332  14341.SAMPLE333  14341.SAMPLE334  14341.SAMPLE335 
    ##                1                1                1                1 
    ##  14341.SAMPLE336  14341.SAMPLE337  14341.SAMPLE338  14341.SAMPLE339 
    ##                1                1                1                1 
    ##   14341.SAMPLE34  14341.SAMPLE340  14341.SAMPLE341  14341.SAMPLE342 
    ##                1                1                1                1 
    ##  14341.SAMPLE343  14341.SAMPLE344  14341.SAMPLE345  14341.SAMPLE346 
    ##                1                1                1                1 
    ##  14341.SAMPLE347  14341.SAMPLE348  14341.SAMPLE349   14341.SAMPLE35 
    ##                1                1                1                1 
    ##  14341.SAMPLE350  14341.SAMPLE351  14341.SAMPLE352  14341.SAMPLE353 
    ##                1                1                1                1 
    ##  14341.SAMPLE354  14341.SAMPLE355  14341.SAMPLE356  14341.SAMPLE357 
    ##                1                1                1                1 
    ##  14341.SAMPLE358  14341.SAMPLE359   14341.SAMPLE36  14341.SAMPLE360 
    ##                1                1                1                1 
    ##  14341.SAMPLE361  14341.SAMPLE362  14341.SAMPLE363  14341.SAMPLE364 
    ##                1                1                1                1 
    ##  14341.SAMPLE365  14341.SAMPLE366  14341.SAMPLE367  14341.SAMPLE368 
    ##                1                1                1                1 
    ##  14341.SAMPLE369   14341.SAMPLE37  14341.SAMPLE370  14341.SAMPLE371 
    ##                1                1                1                1 
    ##  14341.SAMPLE372  14341.SAMPLE373  14341.SAMPLE374  14341.SAMPLE375 
    ##                1                1                1                1 
    ##  14341.SAMPLE376  14341.SAMPLE377  14341.SAMPLE378  14341.SAMPLE379 
    ##                1                1                1                1 
    ##   14341.SAMPLE38  14341.SAMPLE380  14341.SAMPLE381  14341.SAMPLE382 
    ##                1                1                1                1 
    ##  14341.SAMPLE383  14341.SAMPLE384  14341.SAMPLE385  14341.SAMPLE386 
    ##                1                1                1                1 
    ##  14341.SAMPLE387  14341.SAMPLE388  14341.SAMPLE389   14341.SAMPLE39 
    ##                1                1                1                1 
    ##  14341.SAMPLE390  14341.SAMPLE391  14341.SAMPLE392  14341.SAMPLE393 
    ##                1                1                1                1 
    ##  14341.SAMPLE394  14341.SAMPLE395  14341.SAMPLE396  14341.SAMPLE397 
    ##                1                1                1                1 
    ##  14341.SAMPLE398  14341.SAMPLE399    14341.SAMPLE4   14341.SAMPLE40 
    ##                1                1                1                1 
    ##  14341.SAMPLE400  14341.SAMPLE401  14341.SAMPLE402  14341.SAMPLE403 
    ##                1                1                1                1 
    ##  14341.SAMPLE404  14341.SAMPLE405  14341.SAMPLE406  14341.SAMPLE407 
    ##                1                1                1                1 
    ##  14341.SAMPLE408  14341.SAMPLE409   14341.SAMPLE41  14341.SAMPLE410 
    ##                1                1                1                1 
    ##  14341.SAMPLE411  14341.SAMPLE412  14341.SAMPLE413  14341.SAMPLE414 
    ##                1                1                1                1 
    ##  14341.SAMPLE415  14341.SAMPLE416  14341.SAMPLE417  14341.SAMPLE418 
    ##                1                1                1                1 
    ##  14341.SAMPLE419   14341.SAMPLE42  14341.SAMPLE420  14341.SAMPLE421 
    ##                1                1                1                1 
    ##  14341.SAMPLE422  14341.SAMPLE423  14341.SAMPLE424  14341.SAMPLE425 
    ##                1                1                1                1 
    ##  14341.SAMPLE426  14341.SAMPLE427  14341.SAMPLE428  14341.SAMPLE429 
    ##                1                1                1                1 
    ##   14341.SAMPLE43  14341.SAMPLE430  14341.SAMPLE431  14341.SAMPLE432 
    ##                1                1                1                1 
    ##  14341.SAMPLE433  14341.SAMPLE434  14341.SAMPLE435  14341.SAMPLE436 
    ##                1                1                1                1 
    ##  14341.SAMPLE437  14341.SAMPLE438  14341.SAMPLE439   14341.SAMPLE44 
    ##                1                1                1                1 
    ##  14341.SAMPLE440  14341.SAMPLE441  14341.SAMPLE442  14341.SAMPLE443 
    ##                1                1                1                1 
    ##  14341.SAMPLE444  14341.SAMPLE445  14341.SAMPLE446  14341.SAMPLE447 
    ##                1                1                1                1 
    ##  14341.SAMPLE448  14341.SAMPLE449   14341.SAMPLE45  14341.SAMPLE450 
    ##                1                1                1                1 
    ##  14341.SAMPLE451  14341.SAMPLE452  14341.SAMPLE453  14341.SAMPLE454 
    ##                1                1                1                1 
    ##  14341.SAMPLE455  14341.SAMPLE456  14341.SAMPLE457  14341.SAMPLE458 
    ##                1                1                1                1 
    ##  14341.SAMPLE459   14341.SAMPLE46  14341.SAMPLE460  14341.SAMPLE461 
    ##                1                1                1                1 
    ##  14341.SAMPLE462  14341.SAMPLE463  14341.SAMPLE464  14341.SAMPLE465 
    ##                1                1                1                1 
    ##  14341.SAMPLE466  14341.SAMPLE467  14341.SAMPLE468  14341.SAMPLE469 
    ##                1                1                1                1 
    ##   14341.SAMPLE47  14341.SAMPLE470  14341.SAMPLE471  14341.SAMPLE472 
    ##                1                1                1                1 
    ##  14341.SAMPLE473  14341.SAMPLE474  14341.SAMPLE475  14341.SAMPLE476 
    ##                1                1                1                1 
    ##  14341.SAMPLE477  14341.SAMPLE478  14341.SAMPLE479   14341.SAMPLE48 
    ##                1                1                1                1 
    ##  14341.SAMPLE480  14341.SAMPLE481  14341.SAMPLE482  14341.SAMPLE483 
    ##                1                1                1                1 
    ##  14341.SAMPLE484  14341.SAMPLE485  14341.SAMPLE486  14341.SAMPLE487 
    ##                1                1                1                1 
    ##  14341.SAMPLE488  14341.SAMPLE489   14341.SAMPLE49  14341.SAMPLE490 
    ##                1                1                1                1 
    ##  14341.SAMPLE491  14341.SAMPLE492  14341.SAMPLE493  14341.SAMPLE494 
    ##                1                1                1                1 
    ##  14341.SAMPLE495  14341.SAMPLE496  14341.SAMPLE497  14341.SAMPLE498 
    ##                1                1                1                1 
    ##  14341.SAMPLE499    14341.SAMPLE5   14341.SAMPLE50  14341.SAMPLE500 
    ##                1                1                1                1 
    ##  14341.SAMPLE501  14341.SAMPLE502  14341.SAMPLE503  14341.SAMPLE504 
    ##                1                1                1                1 
    ##  14341.SAMPLE505  14341.SAMPLE506  14341.SAMPLE507  14341.SAMPLE508 
    ##                1                1                1                1 
    ##  14341.SAMPLE509   14341.SAMPLE51  14341.SAMPLE510  14341.SAMPLE511 
    ##                1                1                1                1 
    ##  14341.SAMPLE512  14341.SAMPLE513  14341.SAMPLE514  14341.SAMPLE515 
    ##                1                1                1                1 
    ##  14341.SAMPLE516  14341.SAMPLE517  14341.SAMPLE518  14341.SAMPLE519 
    ##                1                1                1                1 
    ##   14341.SAMPLE52  14341.SAMPLE520  14341.SAMPLE521  14341.SAMPLE522 
    ##                1                1                1                1 
    ##  14341.SAMPLE523  14341.SAMPLE524  14341.SAMPLE525  14341.SAMPLE526 
    ##                1                1                1                1 
    ##  14341.SAMPLE527  14341.SAMPLE528  14341.SAMPLE529   14341.SAMPLE53 
    ##                1                1                1                1 
    ##  14341.SAMPLE530  14341.SAMPLE531  14341.SAMPLE532  14341.SAMPLE533 
    ##                1                1                1                1 
    ##  14341.SAMPLE534  14341.SAMPLE535  14341.SAMPLE536  14341.SAMPLE537 
    ##                1                1                1                1 
    ##  14341.SAMPLE538  14341.SAMPLE539   14341.SAMPLE54  14341.SAMPLE540 
    ##                1                1                1                1 
    ##  14341.SAMPLE541  14341.SAMPLE542  14341.SAMPLE543  14341.SAMPLE544 
    ##                1                1                1                1 
    ##  14341.SAMPLE545  14341.SAMPLE546  14341.SAMPLE547  14341.SAMPLE548 
    ##                1                1                1                1 
    ##  14341.SAMPLE549   14341.SAMPLE55  14341.SAMPLE550  14341.SAMPLE551 
    ##                1                1                1                1 
    ##  14341.SAMPLE552  14341.SAMPLE553  14341.SAMPLE554  14341.SAMPLE555 
    ##                1                1                1                1 
    ##  14341.SAMPLE556  14341.SAMPLE557  14341.SAMPLE558  14341.SAMPLE559 
    ##                1                1                1                1 
    ##   14341.SAMPLE56  14341.SAMPLE560  14341.SAMPLE561  14341.SAMPLE562 
    ##                1                1                1                1 
    ##  14341.SAMPLE563  14341.SAMPLE564  14341.SAMPLE565  14341.SAMPLE566 
    ##                1                1                1                1 
    ##  14341.SAMPLE567  14341.SAMPLE568  14341.SAMPLE569   14341.SAMPLE57 
    ##                1                1                1                1 
    ##  14341.SAMPLE570  14341.SAMPLE571  14341.SAMPLE572  14341.SAMPLE573 
    ##                1                1                1                1 
    ##  14341.SAMPLE574  14341.SAMPLE575  14341.SAMPLE576  14341.SAMPLE577 
    ##                1                1                1                1 
    ##  14341.SAMPLE578  14341.SAMPLE579   14341.SAMPLE58  14341.SAMPLE580 
    ##                1                1                1                1 
    ##  14341.SAMPLE581  14341.SAMPLE582  14341.SAMPLE583  14341.SAMPLE584 
    ##                1                1                1                1 
    ##  14341.SAMPLE585  14341.SAMPLE586  14341.SAMPLE587  14341.SAMPLE588 
    ##                1                1                1                1 
    ##  14341.SAMPLE589   14341.SAMPLE59  14341.SAMPLE590  14341.SAMPLE591 
    ##                1                1                1                1 
    ##  14341.SAMPLE592  14341.SAMPLE593  14341.SAMPLE594  14341.SAMPLE595 
    ##                1                1                1                1 
    ##  14341.SAMPLE596  14341.SAMPLE597  14341.SAMPLE598  14341.SAMPLE599 
    ##                1                1                1                1 
    ##    14341.SAMPLE6   14341.SAMPLE60  14341.SAMPLE600  14341.SAMPLE601 
    ##                1                1                1                1 
    ##  14341.SAMPLE602  14341.SAMPLE603  14341.SAMPLE604  14341.SAMPLE605 
    ##                1                1                1                1 
    ##  14341.SAMPLE606  14341.SAMPLE607  14341.SAMPLE608  14341.SAMPLE609 
    ##                1                1                1                1 
    ##   14341.SAMPLE61  14341.SAMPLE610  14341.SAMPLE611  14341.SAMPLE612 
    ##                1                1                1                1 
    ##  14341.SAMPLE613  14341.SAMPLE614  14341.SAMPLE615  14341.SAMPLE616 
    ##                1                1                1                1 
    ##  14341.SAMPLE617  14341.SAMPLE618  14341.SAMPLE619   14341.SAMPLE62 
    ##                1                1                1                1 
    ##  14341.SAMPLE620  14341.SAMPLE621  14341.SAMPLE622  14341.SAMPLE623 
    ##                1                1                1                1 
    ##  14341.SAMPLE624  14341.SAMPLE625  14341.SAMPLE626  14341.SAMPLE627 
    ##                1                1                1                1 
    ##  14341.SAMPLE628  14341.SAMPLE629   14341.SAMPLE63  14341.SAMPLE630 
    ##                1                1                1                1 
    ##  14341.SAMPLE631  14341.SAMPLE632  14341.SAMPLE633  14341.SAMPLE634 
    ##                1                1                1                1 
    ##  14341.SAMPLE635  14341.SAMPLE636  14341.SAMPLE637  14341.SAMPLE638 
    ##                1                1                1                1 
    ##  14341.SAMPLE639   14341.SAMPLE64  14341.SAMPLE640  14341.SAMPLE641 
    ##                1                1                1                1 
    ##  14341.SAMPLE642  14341.SAMPLE643  14341.SAMPLE644  14341.SAMPLE645 
    ##                1                1                1                1 
    ##  14341.SAMPLE646  14341.SAMPLE647  14341.SAMPLE648  14341.SAMPLE649 
    ##                1                1                1                1 
    ##   14341.SAMPLE65  14341.SAMPLE650  14341.SAMPLE651  14341.SAMPLE652 
    ##                1                1                1                1 
    ##  14341.SAMPLE653  14341.SAMPLE654  14341.SAMPLE655  14341.SAMPLE656 
    ##                1                1                1                1 
    ##  14341.SAMPLE657  14341.SAMPLE658  14341.SAMPLE659   14341.SAMPLE66 
    ##                1                1                1                1 
    ##  14341.SAMPLE660  14341.SAMPLE661  14341.SAMPLE662  14341.SAMPLE663 
    ##                1                1                1                1 
    ##  14341.SAMPLE664  14341.SAMPLE665  14341.SAMPLE666  14341.SAMPLE667 
    ##                1                1                1                1 
    ##  14341.SAMPLE668  14341.SAMPLE669   14341.SAMPLE67  14341.SAMPLE670 
    ##                1                1                1                1 
    ##  14341.SAMPLE671  14341.SAMPLE672  14341.SAMPLE673  14341.SAMPLE674 
    ##                1                1                1                1 
    ##  14341.SAMPLE675  14341.SAMPLE676  14341.SAMPLE677  14341.SAMPLE678 
    ##                1                1                1                1 
    ##  14341.SAMPLE679   14341.SAMPLE68  14341.SAMPLE680  14341.SAMPLE681 
    ##                1                1                1                1 
    ##  14341.SAMPLE682  14341.SAMPLE683  14341.SAMPLE684  14341.SAMPLE685 
    ##                1                1                1                1 
    ##  14341.SAMPLE686  14341.SAMPLE687  14341.SAMPLE688  14341.SAMPLE689 
    ##                1                1                1                1 
    ##   14341.SAMPLE69  14341.SAMPLE690  14341.SAMPLE691  14341.SAMPLE692 
    ##                1                1                1                1 
    ##  14341.SAMPLE693  14341.SAMPLE694  14341.SAMPLE695  14341.SAMPLE696 
    ##                1                1                1                1 
    ##  14341.SAMPLE697  14341.SAMPLE698  14341.SAMPLE699    14341.SAMPLE7 
    ##                1                1                1                1 
    ##   14341.SAMPLE70  14341.SAMPLE700  14341.SAMPLE701  14341.SAMPLE702 
    ##                1                1                1                1 
    ##  14341.SAMPLE703  14341.SAMPLE704  14341.SAMPLE705  14341.SAMPLE706 
    ##                1                1                1                1 
    ##  14341.SAMPLE707  14341.SAMPLE708  14341.SAMPLE709   14341.SAMPLE71 
    ##                1                1                1                1 
    ##  14341.SAMPLE710  14341.SAMPLE711  14341.SAMPLE712  14341.SAMPLE713 
    ##                1                1                1                1 
    ##  14341.SAMPLE714  14341.SAMPLE715  14341.SAMPLE716  14341.SAMPLE717 
    ##                1                1                1                1 
    ##  14341.SAMPLE718  14341.SAMPLE719   14341.SAMPLE72  14341.SAMPLE720 
    ##                1                1                1                1 
    ##  14341.SAMPLE721  14341.SAMPLE722  14341.SAMPLE723  14341.SAMPLE724 
    ##                1                1                1                1 
    ##  14341.SAMPLE725  14341.SAMPLE726  14341.SAMPLE727  14341.SAMPLE728 
    ##                1                1                1                1 
    ##  14341.SAMPLE729   14341.SAMPLE73  14341.SAMPLE730  14341.SAMPLE731 
    ##                1                1                1                1 
    ##  14341.SAMPLE732  14341.SAMPLE733  14341.SAMPLE734  14341.SAMPLE735 
    ##                1                1                1                1 
    ##  14341.SAMPLE736  14341.SAMPLE737  14341.SAMPLE738  14341.SAMPLE739 
    ##                1                1                1                1 
    ##   14341.SAMPLE74  14341.SAMPLE740  14341.SAMPLE741  14341.SAMPLE742 
    ##                1                1                1                1 
    ##  14341.SAMPLE743  14341.SAMPLE744  14341.SAMPLE745  14341.SAMPLE746 
    ##                1                1                1                1 
    ##  14341.SAMPLE747  14341.SAMPLE748  14341.SAMPLE749   14341.SAMPLE75 
    ##                1                1                1                1 
    ##  14341.SAMPLE750  14341.SAMPLE751  14341.SAMPLE752  14341.SAMPLE753 
    ##                1                1                1                1 
    ##  14341.SAMPLE754  14341.SAMPLE755  14341.SAMPLE756  14341.SAMPLE757 
    ##                1                1                1                1 
    ##  14341.SAMPLE758  14341.SAMPLE759   14341.SAMPLE76  14341.SAMPLE760 
    ##                1                1                1                1 
    ##  14341.SAMPLE761  14341.SAMPLE762  14341.SAMPLE763  14341.SAMPLE764 
    ##                1                1                1                1 
    ##  14341.SAMPLE765  14341.SAMPLE766  14341.SAMPLE767  14341.SAMPLE768 
    ##                1                1                1                1 
    ##  14341.SAMPLE769   14341.SAMPLE77  14341.SAMPLE770  14341.SAMPLE771 
    ##                1                1                1                1 
    ##  14341.SAMPLE772  14341.SAMPLE773  14341.SAMPLE774  14341.SAMPLE775 
    ##                1                1                1                1 
    ##  14341.SAMPLE776  14341.SAMPLE777  14341.SAMPLE778  14341.SAMPLE779 
    ##                1                1                1                1 
    ##   14341.SAMPLE78  14341.SAMPLE780  14341.SAMPLE781  14341.SAMPLE782 
    ##                1                1                1                1 
    ##  14341.SAMPLE783  14341.SAMPLE784  14341.SAMPLE785  14341.SAMPLE786 
    ##                1                1                1                1 
    ##  14341.SAMPLE787  14341.SAMPLE788  14341.SAMPLE789   14341.SAMPLE79 
    ##                1                1                1                1 
    ##  14341.SAMPLE790  14341.SAMPLE791  14341.SAMPLE792  14341.SAMPLE793 
    ##                1                1                1                1 
    ##  14341.SAMPLE794  14341.SAMPLE795  14341.SAMPLE796  14341.SAMPLE797 
    ##                1                1                1                1 
    ##  14341.SAMPLE798  14341.SAMPLE799    14341.SAMPLE8   14341.SAMPLE80 
    ##                1                1                1                1 
    ##  14341.SAMPLE800  14341.SAMPLE801  14341.SAMPLE802  14341.SAMPLE803 
    ##                1                1                1                1 
    ##  14341.SAMPLE804  14341.SAMPLE805  14341.SAMPLE806  14341.SAMPLE807 
    ##                1                1                1                1 
    ##  14341.SAMPLE808  14341.SAMPLE809   14341.SAMPLE81  14341.SAMPLE810 
    ##                1                1                1                1 
    ##  14341.SAMPLE811  14341.SAMPLE812  14341.SAMPLE813  14341.SAMPLE814 
    ##                1                1                1                1 
    ##  14341.SAMPLE815  14341.SAMPLE816  14341.SAMPLE817  14341.SAMPLE818 
    ##                1                1                1                1 
    ##  14341.SAMPLE819   14341.SAMPLE82  14341.SAMPLE820  14341.SAMPLE821 
    ##                1                1                1                1 
    ##  14341.SAMPLE822  14341.SAMPLE823  14341.SAMPLE824  14341.SAMPLE825 
    ##                1                1                1                1 
    ##  14341.SAMPLE826  14341.SAMPLE827  14341.SAMPLE828  14341.SAMPLE829 
    ##                1                1                1                1 
    ##   14341.SAMPLE83  14341.SAMPLE830  14341.SAMPLE831  14341.SAMPLE832 
    ##                1                1                1                1 
    ##  14341.SAMPLE833  14341.SAMPLE834  14341.SAMPLE835  14341.SAMPLE836 
    ##                1                1                1                1 
    ##  14341.SAMPLE837  14341.SAMPLE838  14341.SAMPLE839   14341.SAMPLE84 
    ##                1                1                1                1 
    ##  14341.SAMPLE840  14341.SAMPLE841  14341.SAMPLE842  14341.SAMPLE843 
    ##                1                1                1                1 
    ##  14341.SAMPLE844  14341.SAMPLE845  14341.SAMPLE846  14341.SAMPLE847 
    ##                1                1                1                1 
    ##  14341.SAMPLE848  14341.SAMPLE849   14341.SAMPLE85  14341.SAMPLE850 
    ##                1                1                1                1 
    ##  14341.SAMPLE851  14341.SAMPLE852  14341.SAMPLE853  14341.SAMPLE854 
    ##                1                1                1                1 
    ##  14341.SAMPLE855  14341.SAMPLE856  14341.SAMPLE857  14341.SAMPLE858 
    ##                1                1                1                1 
    ##  14341.SAMPLE859   14341.SAMPLE86  14341.SAMPLE860  14341.SAMPLE861 
    ##                1                1                1                1 
    ##  14341.SAMPLE862  14341.SAMPLE863  14341.SAMPLE864  14341.SAMPLE865 
    ##                1                1                1                1 
    ##  14341.SAMPLE866  14341.SAMPLE867  14341.SAMPLE868  14341.SAMPLE869 
    ##                1                1                1                1 
    ##   14341.SAMPLE87  14341.SAMPLE870  14341.SAMPLE871  14341.SAMPLE872 
    ##                1                1                1                1 
    ##  14341.SAMPLE873  14341.SAMPLE874  14341.SAMPLE875  14341.SAMPLE876 
    ##                1                1                1                1 
    ##  14341.SAMPLE877  14341.SAMPLE878  14341.SAMPLE879   14341.SAMPLE88 
    ##                1                1                1                1 
    ##  14341.SAMPLE880  14341.SAMPLE881  14341.SAMPLE882  14341.SAMPLE883 
    ##                1                1                1                1 
    ##  14341.SAMPLE884  14341.SAMPLE885  14341.SAMPLE886  14341.SAMPLE887 
    ##                1                1                1                1 
    ##  14341.SAMPLE888  14341.SAMPLE889   14341.SAMPLE89  14341.SAMPLE890 
    ##                1                1                1                1 
    ##  14341.SAMPLE891  14341.SAMPLE892  14341.SAMPLE893  14341.SAMPLE894 
    ##                1                1                1                1 
    ##  14341.SAMPLE895  14341.SAMPLE896  14341.SAMPLE897  14341.SAMPLE898 
    ##                1                1                1                1 
    ##  14341.SAMPLE899    14341.SAMPLE9   14341.SAMPLE90  14341.SAMPLE900 
    ##                1                1                1                1 
    ##  14341.SAMPLE901  14341.SAMPLE902  14341.SAMPLE903  14341.SAMPLE904 
    ##                1                1                1                1 
    ##  14341.SAMPLE905  14341.SAMPLE906  14341.SAMPLE907  14341.SAMPLE908 
    ##                1                1                1                1 
    ##  14341.SAMPLE909   14341.SAMPLE91  14341.SAMPLE910  14341.SAMPLE911 
    ##                1                1                1                1 
    ##  14341.SAMPLE912  14341.SAMPLE913  14341.SAMPLE914  14341.SAMPLE915 
    ##                1                1                1                1 
    ##  14341.SAMPLE916  14341.SAMPLE917  14341.SAMPLE918  14341.SAMPLE919 
    ##                1                1                1                1 
    ##   14341.SAMPLE92  14341.SAMPLE920  14341.SAMPLE921  14341.SAMPLE922 
    ##                1                1                1                1 
    ##  14341.SAMPLE923  14341.SAMPLE924  14341.SAMPLE925  14341.SAMPLE926 
    ##                1                1                1                1 
    ##  14341.SAMPLE927  14341.SAMPLE928  14341.SAMPLE929   14341.SAMPLE93 
    ##                1                1                1                1 
    ##  14341.SAMPLE930  14341.SAMPLE931  14341.SAMPLE932  14341.SAMPLE933 
    ##                1                1                1                1 
    ##  14341.SAMPLE934  14341.SAMPLE935  14341.SAMPLE936  14341.SAMPLE937 
    ##                1                1                1                1 
    ##  14341.SAMPLE938  14341.SAMPLE939   14341.SAMPLE94  14341.SAMPLE940 
    ##                1                1                1                1 
    ##  14341.SAMPLE941  14341.SAMPLE942  14341.SAMPLE943  14341.SAMPLE944 
    ##                1                1                1                1 
    ##  14341.SAMPLE945  14341.SAMPLE946  14341.SAMPLE947  14341.SAMPLE948 
    ##                1                1                1                1 
    ##  14341.SAMPLE949   14341.SAMPLE95  14341.SAMPLE950  14341.SAMPLE951 
    ##                1                1                1                1 
    ##  14341.SAMPLE952  14341.SAMPLE953  14341.SAMPLE954  14341.SAMPLE955 
    ##                1                1                1                1 
    ##  14341.SAMPLE956  14341.SAMPLE957  14341.SAMPLE958  14341.SAMPLE959 
    ##                1                1                1                1 
    ##   14341.SAMPLE96  14341.SAMPLE960  14341.SAMPLE961  14341.SAMPLE962 
    ##                1                1                1                1 
    ##  14341.SAMPLE963  14341.SAMPLE964  14341.SAMPLE965  14341.SAMPLE966 
    ##                1                1                1                1 
    ##  14341.SAMPLE967  14341.SAMPLE968  14341.SAMPLE969   14341.SAMPLE97 
    ##                1                1                1                1 
    ##  14341.SAMPLE970  14341.SAMPLE971  14341.SAMPLE972  14341.SAMPLE973 
    ##                1                1                1                1 
    ##  14341.SAMPLE974  14341.SAMPLE975  14341.SAMPLE976  14341.SAMPLE977 
    ##                1                1                1                1 
    ##  14341.SAMPLE978  14341.SAMPLE979   14341.SAMPLE98  14341.SAMPLE980 
    ##                1                1                1                1 
    ##  14341.SAMPLE981  14341.SAMPLE982  14341.SAMPLE983  14341.SAMPLE984 
    ##                1                1                1                1 
    ##  14341.SAMPLE985  14341.SAMPLE986  14341.SAMPLE987  14341.SAMPLE988 
    ##                1                1                1                1 
    ##  14341.SAMPLE989   14341.SAMPLE99  14341.SAMPLE990  14341.SAMPLE991 
    ##                1                1                1                1 
    ##  14341.SAMPLE992  14341.SAMPLE993  14341.SAMPLE994  14341.SAMPLE995 
    ##                1                1                1                1 
    ##  14341.SAMPLE996  14341.SAMPLE997  14341.SAMPLE998  14341.SAMPLE999 
    ##                1                1                1                1

``` r
# Create matched_samples column
metadata$matched_samples <- ifelse(metadata$sample_name %in% names(sample_counts[sample_counts >= 2]), "YES", "NO")
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(metadata$Run))))
```

    ## [1] "Number of unique runs: 2528"

``` r
print(paste("Number of unique sample names:", length(unique(metadata$sample_name))))
```

    ## [1] "Number of unique sample names: 2528"

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
  "PCR_primers", "title", "Disease_study", "Status", "matched_samples"
)]

# Export the uniform metadata to a CSV file
write.csv(metadata, file = "metadata.csv", quote = FALSE)
head(metadata)
```

    ##   BioProject         Run     sample_name Age BMI Sex Region       Platform
    ## 1 PRJEB86033 ERR14648080  14341.SAMPLE10  NA  NA   F  V3-V4 Illumina MiSeq
    ## 2 PRJEB86033 ERR14648081 14341.SAMPLE100  NA  NA   F  V3-V4 Illumina MiSeq
    ## 3 PRJEB86033 ERR14648082 14341.SAMPLE101  NA  NA   F  V3-V4 Illumina MiSeq
    ## 4 PRJEB86033 ERR14648083 14341.SAMPLE102  NA  NA   F  V3-V4 Illumina MiSeq
    ## 5 PRJEB86033 ERR14648085 14341.SAMPLE104  NA  NA   F  V3-V4 Illumina MiSeq
    ## 6 PRJEB86033 ERR14648086 14341.SAMPLE105  NA  NA   F  V3-V4 Illumina MiSeq
    ##   Country          Site          public_time PCR_primers
    ## 1   China dental plaque 2025-02-22T00:00:00Z   341F-805R
    ## 2   China dental plaque 2025-02-22T00:00:00Z   341F-805R
    ## 3   China dental plaque 2025-02-22T00:00:00Z   341F-805R
    ## 4   China dental plaque 2025-02-22T00:00:00Z   341F-805R
    ## 5   China dental plaque 2025-02-22T00:00:00Z   341F-805R
    ## 6   China dental plaque 2025-02-22T00:00:00Z   341F-805R
    ##                                                                                                                         title
    ## 1 Single-tooth resolved, whole-mouth prediction of early childhood caries\nvia spatiotemporal variations of plaque microbiota
    ## 2 Single-tooth resolved, whole-mouth prediction of early childhood caries\nvia spatiotemporal variations of plaque microbiota
    ## 3 Single-tooth resolved, whole-mouth prediction of early childhood caries\nvia spatiotemporal variations of plaque microbiota
    ## 4 Single-tooth resolved, whole-mouth prediction of early childhood caries\nvia spatiotemporal variations of plaque microbiota
    ## 5 Single-tooth resolved, whole-mouth prediction of early childhood caries\nvia spatiotemporal variations of plaque microbiota
    ## 6 Single-tooth resolved, whole-mouth prediction of early childhood caries\nvia spatiotemporal variations of plaque microbiota
    ##   Disease_study Status matched_samples
    ## 1           YES      H              NO
    ## 2           YES      C              NO
    ## 3           YES      C              NO
    ## 4           YES      C              NO
    ## 5           YES      C              NO
    ## 6           YES      C              NO

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
