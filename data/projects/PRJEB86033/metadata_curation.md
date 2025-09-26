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
  Region = initial_meta$target_subfragment..exp.,
  Platform = initial_meta$Instrument,
  Site = initial_meta$env_material,
  public_time = initial_meta$ReleaseDate,
  Status = initial_meta$status_host,
  matched_samples = NA,
  stringsAsFactors = FALSE
)

# Remove all occurrences of "region" at the end of the Region column
# Use abbreviations to indicate gender
metadata <- metadata %>%
  mutate(Region = str_replace(Region, "region$", ""),
         Sex = recode(Sex, "female" = "F", "male" = "M"))
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

#extract subject id by remove prefix of sample_name
metadata$sample_name <- gsub("14341", "", metadata$sample_name)

# Check initial data distribution
cat("Initial sample distribution:\n")
```

    ## Initial sample distribution:

``` r
table(metadata$sample_name)
```

    ## 
    ##      .A0001      .A0002      .A0003      .A0004      .A0005      .A0006 
    ##           1           1           1           1           1           1 
    ##      .A0007      .A0008      .A0009      .A0010      .A0011      .A0012 
    ##           1           1           1           1           1           1 
    ##      .A0013      .A0014      .A0015      .A0016      .A0017      .A0018 
    ##           1           1           1           1           1           1 
    ##      .A0019      .A0020      .A0021      .A0022      .A0023      .A0024 
    ##           1           1           1           1           1           1 
    ##      .A0025      .A0026      .A0027      .A0028      .A0029      .A0030 
    ##           1           1           1           1           1           1 
    ##      .A0031      .A0032      .A0033      .A0034      .A0035      .A0036 
    ##           1           1           1           1           1           1 
    ##      .A0037      .A0038      .A0039      .A0040      .A0041      .A0042 
    ##           1           1           1           1           1           1 
    ##      .A0043      .A0044      .A0045      .A0046      .A0047      .A0048 
    ##           1           1           1           1           1           1 
    ##      .A0049      .A0050      .A0051      .A0052      .A0053      .A0054 
    ##           1           1           1           1           1           1 
    ##      .A0055      .A0056      .A0057      .A0058      .A0059      .A0060 
    ##           1           1           1           1           1           1 
    ##      .A0061      .A0062      .A0063      .A0064      .A0065      .A0066 
    ##           1           1           1           1           1           1 
    ##      .A0067      .A0068      .A0069      .A0070      .A0071      .A0072 
    ##           1           1           1           1           1           1 
    ##      .A0073      .A0074      .A0075      .A0076      .A0077      .A0078 
    ##           1           1           1           1           1           1 
    ##      .A0079      .A0080      .A0081      .A0082      .A0083      .A0084 
    ##           1           1           1           1           1           1 
    ##      .A0085      .A0086      .A0087      .A0088      .A0089      .A0090 
    ##           1           1           1           1           1           1 
    ##      .A0091      .A0092      .A0093      .A0094      .A0095      .A0096 
    ##           1           1           1           1           1           1 
    ##      .A0097      .A0098      .A0099      .A0100      .A0101      .A0102 
    ##           1           1           1           1           1           1 
    ##      .A0103      .A0104      .A0105      .A0106      .A0107      .A0108 
    ##           1           1           1           1           1           1 
    ##      .A0109      .A0110      .A0111      .A0112      .A0113      .A0114 
    ##           1           1           1           1           1           1 
    ##      .A0115      .A0116      .A0117      .A0118      .A0119      .A0120 
    ##           1           1           1           1           1           1 
    ##      .A0121      .A0122      .A0123      .A0124      .A0125      .A0126 
    ##           1           1           1           1           1           1 
    ##      .A0127      .A0128      .A0129      .A0130      .A0131      .A0132 
    ##           1           1           1           1           1           1 
    ##      .A0133      .A0134      .A0135      .A0136      .A0137      .A0138 
    ##           1           1           1           1           1           1 
    ##      .A0139      .A0140      .A0141      .A0142      .A0143      .A0144 
    ##           1           1           1           1           1           1 
    ##      .A0145      .A0146      .A0147      .A0148      .A0149      .A0150 
    ##           1           1           1           1           1           1 
    ##      .A0151      .A0152      .A0153      .A0154      .A0155      .A0156 
    ##           1           1           1           1           1           1 
    ##      .A0157      .A0158      .A0159      .A0160      .A0161      .A0162 
    ##           1           1           1           1           1           1 
    ##      .A0163      .A0164      .A0165      .A0166      .A0167      .A0168 
    ##           1           1           1           1           1           1 
    ##      .A0169      .A0170      .A0171      .A0172      .A0173      .A0174 
    ##           1           1           1           1           1           1 
    ##      .A0175      .A0176      .A0177      .A0178      .A0179      .A0180 
    ##           1           1           1           1           1           1 
    ##      .A0181      .A0182      .A0183      .A0184      .A0185      .A0186 
    ##           1           1           1           1           1           1 
    ##      .A0187      .A0188      .A0189      .A0190      .A0191      .A0192 
    ##           1           1           1           1           1           1 
    ##      .A0193      .A0194      .A0195      .A0196      .A0197      .A0198 
    ##           1           1           1           1           1           1 
    ##      .A0199      .A0200      .A0201      .A0202      .A0203      .A0204 
    ##           1           1           1           1           1           1 
    ##      .A0205      .A0206      .A0207      .A0208      .A0209      .A0210 
    ##           1           1           1           1           1           1 
    ##      .A0211      .A0212      .A0213      .A0214      .A0215      .A0216 
    ##           1           1           1           1           1           1 
    ##      .A0217      .A0218      .A0219      .A0220      .A0221      .A0222 
    ##           1           1           1           1           1           1 
    ##      .A0223      .A0224      .A0225      .A0226      .A0227      .A0228 
    ##           1           1           1           1           1           1 
    ##      .A0229      .A0230      .A0231      .A0232      .A0233      .A0234 
    ##           1           1           1           1           1           1 
    ##      .A0235      .A0236      .A0237      .A0238      .A0239      .A0240 
    ##           1           1           1           1           1           1 
    ##      .A0241      .A0242      .A0243      .A0244      .A0245      .A0246 
    ##           1           1           1           1           1           1 
    ##      .A0247      .A0248      .A0249      .A0250      .A0251      .A0252 
    ##           1           1           1           1           1           1 
    ##      .A0253      .A0254      .A0255      .A0256      .A0257      .A0258 
    ##           1           1           1           1           1           1 
    ##      .A0259      .A0260      .A0261      .A0262      .A0263      .A0264 
    ##           1           1           1           1           1           1 
    ##      .A0265      .A0266      .A0267      .A0268      .A0269      .A0270 
    ##           1           1           1           1           1           1 
    ##      .A0271      .A0272      .A0273      .A0274      .A0275      .A0276 
    ##           1           1           1           1           1           1 
    ##      .A0277      .A0278      .A0279      .A0280      .A0281      .A0282 
    ##           1           1           1           1           1           1 
    ##      .A0283      .A0284      .A0285      .A0286      .A0287      .A0288 
    ##           1           1           1           1           1           1 
    ##      .A0289      .A0290      .A0291      .A0292      .A0293      .A0294 
    ##           1           1           1           1           1           1 
    ##      .A0295      .A0296      .A0297      .A0298      .A0299      .A0300 
    ##           1           1           1           1           1           1 
    ##      .A0301      .A0302      .A0303      .A0304      .A0305      .A0306 
    ##           1           1           1           1           1           1 
    ##      .A0307      .A0308      .A0309      .A0310      .A0311      .A0312 
    ##           1           1           1           1           1           1 
    ##      .A0313      .A0314      .A0315      .A0316      .A0317      .A0318 
    ##           1           1           1           1           1           1 
    ##      .A0319      .A0320      .A0321      .A0322      .A0323      .A0324 
    ##           1           1           1           1           1           1 
    ##      .A0325      .A0326      .A0327      .A0328      .A0329      .A0330 
    ##           1           1           1           1           1           1 
    ##      .A0331      .A0332      .A0333      .A0334      .A0335      .A0336 
    ##           1           1           1           1           1           1 
    ##      .A0337      .A0338      .A0339      .A0340      .A0341      .A0342 
    ##           1           1           1           1           1           1 
    ##      .A0343      .A0344      .A0345      .A0346      .A0347      .A0348 
    ##           1           1           1           1           1           1 
    ##      .A0349      .A0350      .A0351      .A0352      .A0353      .A0354 
    ##           1           1           1           1           1           1 
    ##      .A0355      .A0356      .A0357      .A0358      .A0359      .A0360 
    ##           1           1           1           1           1           1 
    ##      .A0361      .A0362      .A0363      .A0364      .A0365      .A0366 
    ##           1           1           1           1           1           1 
    ##      .A0367      .A0368      .A0369      .A0370      .A0371      .A0372 
    ##           1           1           1           1           1           1 
    ##      .A0373      .A0374      .A0375      .A0376      .A0377      .A0378 
    ##           1           1           1           1           1           1 
    ##      .A0379      .A0380      .A0381      .A0382      .A0383      .A0384 
    ##           1           1           1           1           1           1 
    ##      .A0385      .A0386      .A0387      .A0388      .A0389      .A0390 
    ##           1           1           1           1           1           1 
    ##      .A0391      .A0392      .A0393      .A0394      .A0395      .A0396 
    ##           1           1           1           1           1           1 
    ##      .A0397      .A0398      .A0399      .A0400      .A0401      .A0402 
    ##           1           1           1           1           1           1 
    ##      .A0403      .A0404      .A0405      .A0406      .A0407      .A0408 
    ##           1           1           1           1           1           1 
    ##      .A0409      .A0410      .A0411      .A0412      .A0413      .A0414 
    ##           1           1           1           1           1           1 
    ##      .A0415      .A0416      .A0417      .A0418      .A0419      .A0420 
    ##           1           1           1           1           1           1 
    ##      .A0421      .A0422      .A0423      .A0424      .A0425      .A0426 
    ##           1           1           1           1           1           1 
    ##      .A0427      .A0428      .A0429      .A0430      .A0431      .A0432 
    ##           1           1           1           1           1           1 
    ##      .A0433      .A0434      .A0435      .A0436      .A0437      .A0438 
    ##           1           1           1           1           1           1 
    ##      .A0439      .A0440      .A0441      .A0442      .A0443      .A0444 
    ##           1           1           1           1           1           1 
    ##      .A0445      .A0446      .A0447      .A0448      .A0449      .A0450 
    ##           1           1           1           1           1           1 
    ##      .A0451      .A0452      .A0453      .A0454      .A0455      .A0456 
    ##           1           1           1           1           1           1 
    ##      .A0457      .A0458      .A0459      .A0460      .A0461      .A0462 
    ##           1           1           1           1           1           1 
    ##      .A0463      .A0464      .A0465      .A0466      .A0467      .A0468 
    ##           1           1           1           1           1           1 
    ##      .A0469      .A0470      .A0471      .A0472      .A0473      .A0474 
    ##           1           1           1           1           1           1 
    ##      .A0475      .A0476      .A0477      .A0478      .A0479      .A0480 
    ##           1           1           1           1           1           1 
    ##      .A0481      .A0482      .A0483      .A0484      .A0485      .A0486 
    ##           1           1           1           1           1           1 
    ##      .A0487      .A0488      .A0489      .A0490      .A0491      .A0492 
    ##           1           1           1           1           1           1 
    ##      .A0493      .A0494      .A0495      .A0496      .A0497      .A0498 
    ##           1           1           1           1           1           1 
    ##      .A0499      .A0500      .A0501      .A0502      .A0503      .A0504 
    ##           1           1           1           1           1           1 
    ##      .A0505      .A0506      .A0507      .A0508      .A0509      .A0510 
    ##           1           1           1           1           1           1 
    ##      .A0511      .A0512      .A0513      .A0514      .A0515      .A0516 
    ##           1           1           1           1           1           1 
    ##      .A0517      .A0518      .A0519      .A0520      .A0521      .A0522 
    ##           1           1           1           1           1           1 
    ##      .A0523      .A0524      .A0525      .A0526      .A0527      .A0528 
    ##           1           1           1           1           1           1 
    ##      .A0529      .A0530      .A0531      .A0532      .A0533      .A0534 
    ##           1           1           1           1           1           1 
    ##      .A0535      .A0536      .A0537      .A0538      .A0539      .A0540 
    ##           1           1           1           1           1           1 
    ##      .A0541      .A0542      .A0543      .A0544      .A0545      .A0546 
    ##           1           1           1           1           1           1 
    ##      .A0547      .A0548      .A0549      .A0550      .A0551      .A0552 
    ##           1           1           1           1           1           1 
    ##      .A0553      .A0554      .A0555      .A0556      .A0557      .A0558 
    ##           1           1           1           1           1           1 
    ##      .A0559      .A0560      .A0561      .A0562      .A0563      .A0564 
    ##           1           1           1           1           1           1 
    ##      .A0565      .A0566      .A0567      .A0568      .A0569      .A0570 
    ##           1           1           1           1           1           1 
    ##      .A0571      .A0572      .A0573      .A0574      .A0575      .A0576 
    ##           1           1           1           1           1           1 
    ##      .A0577      .A0578      .A0579      .A0580      .A0581      .A0582 
    ##           1           1           1           1           1           1 
    ##      .A0583      .A0584      .A0585      .A0586      .A0587      .A0588 
    ##           1           1           1           1           1           1 
    ##      .A0589      .A0590      .A0591      .A0592      .A0593      .A0594 
    ##           1           1           1           1           1           1 
    ##      .A0595      .A0596      .A0597      .A0598      .A0599      .A0600 
    ##           1           1           1           1           1           1 
    ##      .A0601      .A0602      .A0603      .A0604      .A0605      .A0606 
    ##           1           1           1           1           1           1 
    ##      .A0607      .A0608      .A0609      .A0610      .A0611      .A0612 
    ##           1           1           1           1           1           1 
    ##      .A0613      .A0614      .A0615      .A0616      .A0617      .A0618 
    ##           1           1           1           1           1           1 
    ##      .A0619      .A0620      .A0621      .A0622      .A0623      .A0624 
    ##           1           1           1           1           1           1 
    ##      .A0625      .A0626      .A0627      .A0628      .A0629      .A0630 
    ##           1           1           1           1           1           1 
    ##      .A0631      .A0632      .A0633      .A0634      .A0635      .A0636 
    ##           1           1           1           1           1           1 
    ##      .A0637         .S1        .S10        .S11        .S12        .S13 
    ##           1           1           1           1           1           1 
    ##        .S14        .S15        .S16        .S17        .S18        .S19 
    ##           1           1           1           1           1           1 
    ##         .S2        .S20        .S21        .S22        .S23        .S24 
    ##           1           1           1           1           1           1 
    ##         .S3         .S4         .S5         .S6         .S7         .S8 
    ##           1           1           1           1           1           1 
    ##         .S9    .SAMPLE1   .SAMPLE10  .SAMPLE100 .SAMPLE1000 .SAMPLE1001 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1002 .SAMPLE1003 .SAMPLE1004 .SAMPLE1005 .SAMPLE1006 .SAMPLE1007 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1008 .SAMPLE1009  .SAMPLE101 .SAMPLE1010 .SAMPLE1011 .SAMPLE1012 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1013 .SAMPLE1014 .SAMPLE1015 .SAMPLE1016 .SAMPLE1017 .SAMPLE1018 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1019  .SAMPLE102 .SAMPLE1020 .SAMPLE1021 .SAMPLE1022 .SAMPLE1023 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1024 .SAMPLE1025 .SAMPLE1026 .SAMPLE1027 .SAMPLE1028 .SAMPLE1029 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE103 .SAMPLE1030 .SAMPLE1031 .SAMPLE1032 .SAMPLE1033 .SAMPLE1034 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1035 .SAMPLE1036 .SAMPLE1037 .SAMPLE1038 .SAMPLE1039  .SAMPLE104 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1040 .SAMPLE1041 .SAMPLE1042 .SAMPLE1043 .SAMPLE1044 .SAMPLE1045 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1046 .SAMPLE1047 .SAMPLE1048 .SAMPLE1049  .SAMPLE105 .SAMPLE1050 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1051 .SAMPLE1052 .SAMPLE1053 .SAMPLE1054 .SAMPLE1055 .SAMPLE1056 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1057 .SAMPLE1058 .SAMPLE1059  .SAMPLE106 .SAMPLE1060 .SAMPLE1061 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1062 .SAMPLE1063 .SAMPLE1064 .SAMPLE1065 .SAMPLE1066 .SAMPLE1067 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1068 .SAMPLE1069  .SAMPLE107 .SAMPLE1070 .SAMPLE1071 .SAMPLE1072 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1073 .SAMPLE1074 .SAMPLE1075 .SAMPLE1076 .SAMPLE1077 .SAMPLE1078 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1079  .SAMPLE108 .SAMPLE1080 .SAMPLE1081 .SAMPLE1082 .SAMPLE1083 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1084 .SAMPLE1085 .SAMPLE1086 .SAMPLE1087 .SAMPLE1088 .SAMPLE1089 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE109 .SAMPLE1090 .SAMPLE1091 .SAMPLE1092 .SAMPLE1093 .SAMPLE1094 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1095 .SAMPLE1096 .SAMPLE1097 .SAMPLE1098 .SAMPLE1099   .SAMPLE11 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE110 .SAMPLE1100 .SAMPLE1101 .SAMPLE1102 .SAMPLE1103 .SAMPLE1104 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1105 .SAMPLE1106 .SAMPLE1107 .SAMPLE1108 .SAMPLE1109  .SAMPLE111 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1110 .SAMPLE1111 .SAMPLE1112 .SAMPLE1113 .SAMPLE1114 .SAMPLE1115 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1116 .SAMPLE1117 .SAMPLE1118 .SAMPLE1119  .SAMPLE112 .SAMPLE1120 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1121 .SAMPLE1122 .SAMPLE1123 .SAMPLE1124 .SAMPLE1125 .SAMPLE1126 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1127 .SAMPLE1128 .SAMPLE1129  .SAMPLE113 .SAMPLE1130 .SAMPLE1131 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1132 .SAMPLE1133 .SAMPLE1134 .SAMPLE1135 .SAMPLE1136 .SAMPLE1137 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1138 .SAMPLE1139  .SAMPLE114 .SAMPLE1140 .SAMPLE1141 .SAMPLE1142 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1143 .SAMPLE1144 .SAMPLE1145 .SAMPLE1146 .SAMPLE1147 .SAMPLE1148 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1149  .SAMPLE115 .SAMPLE1150 .SAMPLE1151 .SAMPLE1152 .SAMPLE1153 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1154 .SAMPLE1155 .SAMPLE1156 .SAMPLE1157 .SAMPLE1158 .SAMPLE1159 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE116 .SAMPLE1160 .SAMPLE1161 .SAMPLE1162 .SAMPLE1163 .SAMPLE1164 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1165 .SAMPLE1166 .SAMPLE1167 .SAMPLE1168 .SAMPLE1169  .SAMPLE117 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1170 .SAMPLE1171 .SAMPLE1172 .SAMPLE1173 .SAMPLE1174 .SAMPLE1175 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1176 .SAMPLE1177 .SAMPLE1178 .SAMPLE1179  .SAMPLE118 .SAMPLE1180 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1181 .SAMPLE1182 .SAMPLE1183 .SAMPLE1184 .SAMPLE1185 .SAMPLE1186 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1187 .SAMPLE1188 .SAMPLE1189  .SAMPLE119 .SAMPLE1190 .SAMPLE1191 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1192 .SAMPLE1193 .SAMPLE1194 .SAMPLE1195 .SAMPLE1196 .SAMPLE1197 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1198 .SAMPLE1199   .SAMPLE12  .SAMPLE120 .SAMPLE1200 .SAMPLE1201 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1202 .SAMPLE1203 .SAMPLE1204 .SAMPLE1205 .SAMPLE1206 .SAMPLE1207 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1208 .SAMPLE1209  .SAMPLE121 .SAMPLE1210 .SAMPLE1211 .SAMPLE1212 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1213 .SAMPLE1214 .SAMPLE1215 .SAMPLE1216 .SAMPLE1217 .SAMPLE1218 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1219  .SAMPLE122 .SAMPLE1220 .SAMPLE1221 .SAMPLE1222 .SAMPLE1223 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1224 .SAMPLE1225 .SAMPLE1226 .SAMPLE1227 .SAMPLE1228 .SAMPLE1229 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE123 .SAMPLE1230 .SAMPLE1231 .SAMPLE1232 .SAMPLE1233 .SAMPLE1234 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1235 .SAMPLE1236 .SAMPLE1237 .SAMPLE1238 .SAMPLE1239  .SAMPLE124 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1240 .SAMPLE1241 .SAMPLE1242 .SAMPLE1243 .SAMPLE1244 .SAMPLE1245 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1246 .SAMPLE1247 .SAMPLE1248 .SAMPLE1249  .SAMPLE125 .SAMPLE1250 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1251 .SAMPLE1252 .SAMPLE1253 .SAMPLE1254 .SAMPLE1255 .SAMPLE1256 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1257 .SAMPLE1258 .SAMPLE1259  .SAMPLE126 .SAMPLE1260 .SAMPLE1261 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1262 .SAMPLE1263 .SAMPLE1264 .SAMPLE1265 .SAMPLE1266 .SAMPLE1267 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1268 .SAMPLE1269  .SAMPLE127 .SAMPLE1270 .SAMPLE1271 .SAMPLE1272 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1273 .SAMPLE1274 .SAMPLE1275 .SAMPLE1276 .SAMPLE1277 .SAMPLE1278 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1279  .SAMPLE128 .SAMPLE1280 .SAMPLE1281 .SAMPLE1282 .SAMPLE1283 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1284 .SAMPLE1285 .SAMPLE1286 .SAMPLE1287 .SAMPLE1288 .SAMPLE1289 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE129 .SAMPLE1290 .SAMPLE1291 .SAMPLE1292 .SAMPLE1293 .SAMPLE1294 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1295 .SAMPLE1296 .SAMPLE1297 .SAMPLE1298 .SAMPLE1299   .SAMPLE13 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE130 .SAMPLE1300 .SAMPLE1301 .SAMPLE1302 .SAMPLE1303 .SAMPLE1304 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1305 .SAMPLE1306 .SAMPLE1307 .SAMPLE1308 .SAMPLE1309  .SAMPLE131 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1310 .SAMPLE1311 .SAMPLE1312 .SAMPLE1313 .SAMPLE1314 .SAMPLE1315 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1316 .SAMPLE1317 .SAMPLE1318 .SAMPLE1319  .SAMPLE132 .SAMPLE1320 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1321 .SAMPLE1322 .SAMPLE1323 .SAMPLE1324 .SAMPLE1325 .SAMPLE1326 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1327 .SAMPLE1328 .SAMPLE1329  .SAMPLE133 .SAMPLE1330 .SAMPLE1331 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1332 .SAMPLE1333 .SAMPLE1334 .SAMPLE1335 .SAMPLE1336 .SAMPLE1337 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1338 .SAMPLE1339  .SAMPLE134 .SAMPLE1340 .SAMPLE1341 .SAMPLE1342 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1343 .SAMPLE1344 .SAMPLE1345 .SAMPLE1346 .SAMPLE1347 .SAMPLE1348 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1349  .SAMPLE135 .SAMPLE1350 .SAMPLE1351 .SAMPLE1352 .SAMPLE1353 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1354 .SAMPLE1355 .SAMPLE1356 .SAMPLE1357 .SAMPLE1358 .SAMPLE1359 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE136 .SAMPLE1360 .SAMPLE1361 .SAMPLE1362 .SAMPLE1363 .SAMPLE1364 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1365 .SAMPLE1366 .SAMPLE1367 .SAMPLE1368 .SAMPLE1369  .SAMPLE137 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1370 .SAMPLE1371 .SAMPLE1372 .SAMPLE1373 .SAMPLE1374 .SAMPLE1375 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1376 .SAMPLE1377 .SAMPLE1378 .SAMPLE1379  .SAMPLE138 .SAMPLE1380 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1381 .SAMPLE1382 .SAMPLE1383 .SAMPLE1384 .SAMPLE1385 .SAMPLE1386 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1387 .SAMPLE1388 .SAMPLE1389  .SAMPLE139 .SAMPLE1390 .SAMPLE1391 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1392 .SAMPLE1393 .SAMPLE1394 .SAMPLE1395 .SAMPLE1396 .SAMPLE1397 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1398 .SAMPLE1399   .SAMPLE14  .SAMPLE140 .SAMPLE1400 .SAMPLE1401 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1402 .SAMPLE1403 .SAMPLE1404 .SAMPLE1405 .SAMPLE1406 .SAMPLE1407 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1408 .SAMPLE1409  .SAMPLE141 .SAMPLE1410 .SAMPLE1411 .SAMPLE1412 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1413 .SAMPLE1414 .SAMPLE1415 .SAMPLE1416 .SAMPLE1417 .SAMPLE1418 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1419  .SAMPLE142 .SAMPLE1420 .SAMPLE1421 .SAMPLE1422 .SAMPLE1423 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1424 .SAMPLE1425 .SAMPLE1426 .SAMPLE1427 .SAMPLE1428 .SAMPLE1429 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE143 .SAMPLE1430 .SAMPLE1431 .SAMPLE1432 .SAMPLE1433 .SAMPLE1434 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1435 .SAMPLE1436 .SAMPLE1437 .SAMPLE1438 .SAMPLE1439  .SAMPLE144 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1440 .SAMPLE1441 .SAMPLE1442 .SAMPLE1443 .SAMPLE1444 .SAMPLE1445 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1446 .SAMPLE1447 .SAMPLE1448 .SAMPLE1449  .SAMPLE145 .SAMPLE1450 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1451 .SAMPLE1452 .SAMPLE1453 .SAMPLE1454 .SAMPLE1455 .SAMPLE1456 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1457 .SAMPLE1458 .SAMPLE1459  .SAMPLE146 .SAMPLE1460 .SAMPLE1461 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1462 .SAMPLE1463 .SAMPLE1464 .SAMPLE1465 .SAMPLE1466 .SAMPLE1467 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1468 .SAMPLE1469  .SAMPLE147 .SAMPLE1470 .SAMPLE1471 .SAMPLE1472 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1473 .SAMPLE1474 .SAMPLE1475 .SAMPLE1476 .SAMPLE1477 .SAMPLE1478 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1479  .SAMPLE148 .SAMPLE1480 .SAMPLE1481 .SAMPLE1482 .SAMPLE1483 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1484 .SAMPLE1485 .SAMPLE1486 .SAMPLE1487 .SAMPLE1488 .SAMPLE1489 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE149 .SAMPLE1490 .SAMPLE1491 .SAMPLE1492 .SAMPLE1493 .SAMPLE1494 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1495 .SAMPLE1496 .SAMPLE1497 .SAMPLE1498 .SAMPLE1499   .SAMPLE15 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE150 .SAMPLE1500 .SAMPLE1501 .SAMPLE1502 .SAMPLE1503 .SAMPLE1504 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1505 .SAMPLE1506 .SAMPLE1507 .SAMPLE1508 .SAMPLE1509  .SAMPLE151 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1510 .SAMPLE1511 .SAMPLE1512 .SAMPLE1513 .SAMPLE1514 .SAMPLE1515 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1516 .SAMPLE1517 .SAMPLE1518 .SAMPLE1519  .SAMPLE152 .SAMPLE1520 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1521 .SAMPLE1522 .SAMPLE1523 .SAMPLE1524 .SAMPLE1525 .SAMPLE1526 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1527 .SAMPLE1528 .SAMPLE1529  .SAMPLE153 .SAMPLE1530 .SAMPLE1531 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1532 .SAMPLE1533 .SAMPLE1534 .SAMPLE1535 .SAMPLE1536 .SAMPLE1537 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1538 .SAMPLE1539  .SAMPLE154 .SAMPLE1540 .SAMPLE1541 .SAMPLE1542 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1543 .SAMPLE1544 .SAMPLE1545 .SAMPLE1546 .SAMPLE1547 .SAMPLE1548 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1549  .SAMPLE155 .SAMPLE1550 .SAMPLE1551 .SAMPLE1552 .SAMPLE1553 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1554 .SAMPLE1555 .SAMPLE1556 .SAMPLE1557 .SAMPLE1558 .SAMPLE1559 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE156 .SAMPLE1560 .SAMPLE1561 .SAMPLE1562 .SAMPLE1563 .SAMPLE1564 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1565 .SAMPLE1566 .SAMPLE1567 .SAMPLE1568 .SAMPLE1569  .SAMPLE157 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1570 .SAMPLE1571 .SAMPLE1572 .SAMPLE1573 .SAMPLE1574 .SAMPLE1575 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1576 .SAMPLE1577 .SAMPLE1578 .SAMPLE1579  .SAMPLE158 .SAMPLE1580 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1581 .SAMPLE1582 .SAMPLE1583 .SAMPLE1584 .SAMPLE1585 .SAMPLE1586 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1587 .SAMPLE1588 .SAMPLE1589  .SAMPLE159 .SAMPLE1590 .SAMPLE1591 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1592 .SAMPLE1593 .SAMPLE1594 .SAMPLE1595 .SAMPLE1596 .SAMPLE1597 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1598 .SAMPLE1599   .SAMPLE16  .SAMPLE160 .SAMPLE1600 .SAMPLE1601 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1602 .SAMPLE1603 .SAMPLE1604 .SAMPLE1605 .SAMPLE1606 .SAMPLE1607 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1608 .SAMPLE1609  .SAMPLE161 .SAMPLE1610 .SAMPLE1611 .SAMPLE1612 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1613 .SAMPLE1614 .SAMPLE1615 .SAMPLE1616 .SAMPLE1617 .SAMPLE1618 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1619  .SAMPLE162 .SAMPLE1620 .SAMPLE1621 .SAMPLE1622 .SAMPLE1623 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1624 .SAMPLE1625 .SAMPLE1626 .SAMPLE1627 .SAMPLE1628 .SAMPLE1629 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE163 .SAMPLE1630 .SAMPLE1631 .SAMPLE1632 .SAMPLE1633 .SAMPLE1634 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1635 .SAMPLE1636 .SAMPLE1637 .SAMPLE1638 .SAMPLE1639  .SAMPLE164 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1640 .SAMPLE1641 .SAMPLE1642 .SAMPLE1643 .SAMPLE1644 .SAMPLE1645 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1646 .SAMPLE1647 .SAMPLE1648 .SAMPLE1649  .SAMPLE165 .SAMPLE1650 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1651 .SAMPLE1652 .SAMPLE1653 .SAMPLE1654 .SAMPLE1655 .SAMPLE1656 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1657 .SAMPLE1658 .SAMPLE1659  .SAMPLE166 .SAMPLE1660 .SAMPLE1661 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1662 .SAMPLE1663 .SAMPLE1664 .SAMPLE1665 .SAMPLE1666 .SAMPLE1667 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1668 .SAMPLE1669  .SAMPLE167 .SAMPLE1670 .SAMPLE1671 .SAMPLE1672 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1673 .SAMPLE1674 .SAMPLE1675 .SAMPLE1676 .SAMPLE1677 .SAMPLE1678 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1679  .SAMPLE168 .SAMPLE1680 .SAMPLE1681 .SAMPLE1682 .SAMPLE1683 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1684 .SAMPLE1685 .SAMPLE1686 .SAMPLE1687 .SAMPLE1688 .SAMPLE1689 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE169 .SAMPLE1690 .SAMPLE1691 .SAMPLE1692 .SAMPLE1693 .SAMPLE1694 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1695 .SAMPLE1696 .SAMPLE1697 .SAMPLE1698 .SAMPLE1699   .SAMPLE17 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE170 .SAMPLE1700 .SAMPLE1701 .SAMPLE1702 .SAMPLE1703 .SAMPLE1704 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1705 .SAMPLE1706 .SAMPLE1707 .SAMPLE1708 .SAMPLE1709  .SAMPLE171 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1710 .SAMPLE1711 .SAMPLE1712 .SAMPLE1713 .SAMPLE1714 .SAMPLE1715 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1716 .SAMPLE1717 .SAMPLE1718 .SAMPLE1719  .SAMPLE172 .SAMPLE1720 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1721 .SAMPLE1722 .SAMPLE1723 .SAMPLE1724 .SAMPLE1725 .SAMPLE1726 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1727 .SAMPLE1728 .SAMPLE1729  .SAMPLE173 .SAMPLE1730 .SAMPLE1731 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1732 .SAMPLE1733 .SAMPLE1734 .SAMPLE1735 .SAMPLE1736 .SAMPLE1737 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1738 .SAMPLE1739  .SAMPLE174 .SAMPLE1740 .SAMPLE1741 .SAMPLE1742 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1743 .SAMPLE1744 .SAMPLE1745 .SAMPLE1746 .SAMPLE1747 .SAMPLE1748 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1749  .SAMPLE175 .SAMPLE1750 .SAMPLE1751 .SAMPLE1752 .SAMPLE1753 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1754 .SAMPLE1755 .SAMPLE1756 .SAMPLE1757 .SAMPLE1758 .SAMPLE1759 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE176 .SAMPLE1760 .SAMPLE1761 .SAMPLE1762 .SAMPLE1763 .SAMPLE1764 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1765 .SAMPLE1766 .SAMPLE1767 .SAMPLE1768 .SAMPLE1769  .SAMPLE177 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1770 .SAMPLE1771 .SAMPLE1772 .SAMPLE1773 .SAMPLE1774 .SAMPLE1775 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1776 .SAMPLE1777 .SAMPLE1778 .SAMPLE1779  .SAMPLE178 .SAMPLE1780 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1781 .SAMPLE1782 .SAMPLE1783 .SAMPLE1784 .SAMPLE1785 .SAMPLE1786 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1787 .SAMPLE1788 .SAMPLE1789  .SAMPLE179 .SAMPLE1790 .SAMPLE1791 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1792 .SAMPLE1793 .SAMPLE1794 .SAMPLE1795 .SAMPLE1796 .SAMPLE1797 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1798 .SAMPLE1799   .SAMPLE18  .SAMPLE180 .SAMPLE1800 .SAMPLE1801 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1802 .SAMPLE1803 .SAMPLE1804 .SAMPLE1805 .SAMPLE1806 .SAMPLE1807 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1808 .SAMPLE1809  .SAMPLE181 .SAMPLE1810 .SAMPLE1811 .SAMPLE1812 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1813 .SAMPLE1814 .SAMPLE1815 .SAMPLE1816 .SAMPLE1817 .SAMPLE1818 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1819  .SAMPLE182 .SAMPLE1820 .SAMPLE1821 .SAMPLE1822 .SAMPLE1823 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1824 .SAMPLE1825 .SAMPLE1826 .SAMPLE1827 .SAMPLE1828 .SAMPLE1829 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE183 .SAMPLE1830 .SAMPLE1831 .SAMPLE1832 .SAMPLE1833 .SAMPLE1834 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1835 .SAMPLE1836 .SAMPLE1837 .SAMPLE1838 .SAMPLE1839  .SAMPLE184 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1840 .SAMPLE1841 .SAMPLE1842 .SAMPLE1843 .SAMPLE1844 .SAMPLE1845 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1846 .SAMPLE1847 .SAMPLE1848 .SAMPLE1849  .SAMPLE185 .SAMPLE1850 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1851 .SAMPLE1852 .SAMPLE1853 .SAMPLE1854 .SAMPLE1855 .SAMPLE1856 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1857 .SAMPLE1858 .SAMPLE1859  .SAMPLE186 .SAMPLE1860 .SAMPLE1861 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1862 .SAMPLE1863 .SAMPLE1864 .SAMPLE1865 .SAMPLE1866 .SAMPLE1867 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1868 .SAMPLE1869  .SAMPLE187 .SAMPLE1870 .SAMPLE1871 .SAMPLE1872 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE188  .SAMPLE189   .SAMPLE19  .SAMPLE190  .SAMPLE191  .SAMPLE192 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE193  .SAMPLE194  .SAMPLE195  .SAMPLE196  .SAMPLE197  .SAMPLE198 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE199    .SAMPLE2   .SAMPLE20  .SAMPLE200  .SAMPLE201  .SAMPLE202 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE203  .SAMPLE204  .SAMPLE205  .SAMPLE206  .SAMPLE207  .SAMPLE208 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE209   .SAMPLE21  .SAMPLE210  .SAMPLE211  .SAMPLE212  .SAMPLE213 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE214  .SAMPLE215  .SAMPLE216  .SAMPLE217  .SAMPLE218  .SAMPLE219 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE22  .SAMPLE220  .SAMPLE221  .SAMPLE222  .SAMPLE223  .SAMPLE224 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE225  .SAMPLE226  .SAMPLE227  .SAMPLE228  .SAMPLE229   .SAMPLE23 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE230  .SAMPLE231  .SAMPLE232  .SAMPLE233  .SAMPLE234  .SAMPLE235 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE236  .SAMPLE237  .SAMPLE238  .SAMPLE239   .SAMPLE24  .SAMPLE240 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE241  .SAMPLE242  .SAMPLE243  .SAMPLE244  .SAMPLE245  .SAMPLE246 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE247  .SAMPLE248  .SAMPLE249   .SAMPLE25  .SAMPLE250  .SAMPLE251 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE252  .SAMPLE253  .SAMPLE254  .SAMPLE255  .SAMPLE256  .SAMPLE257 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE258  .SAMPLE259   .SAMPLE26  .SAMPLE260  .SAMPLE261  .SAMPLE262 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE263  .SAMPLE264  .SAMPLE265  .SAMPLE266  .SAMPLE267  .SAMPLE268 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE269   .SAMPLE27  .SAMPLE270  .SAMPLE271  .SAMPLE272  .SAMPLE273 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE274  .SAMPLE275  .SAMPLE276  .SAMPLE277  .SAMPLE278  .SAMPLE279 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE28  .SAMPLE280  .SAMPLE281  .SAMPLE282  .SAMPLE283  .SAMPLE284 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE285  .SAMPLE286  .SAMPLE287  .SAMPLE288  .SAMPLE289   .SAMPLE29 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE290  .SAMPLE291  .SAMPLE292  .SAMPLE293  .SAMPLE294  .SAMPLE295 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE296  .SAMPLE297  .SAMPLE298  .SAMPLE299    .SAMPLE3   .SAMPLE30 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE300  .SAMPLE301  .SAMPLE302  .SAMPLE303  .SAMPLE304  .SAMPLE305 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE306  .SAMPLE307  .SAMPLE308  .SAMPLE309   .SAMPLE31  .SAMPLE310 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE311  .SAMPLE312  .SAMPLE313  .SAMPLE314  .SAMPLE315  .SAMPLE316 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE318  .SAMPLE319   .SAMPLE32  .SAMPLE321  .SAMPLE322  .SAMPLE324 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE325  .SAMPLE327  .SAMPLE329   .SAMPLE33  .SAMPLE330  .SAMPLE331 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE332  .SAMPLE333  .SAMPLE334  .SAMPLE335  .SAMPLE336  .SAMPLE337 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE338  .SAMPLE339   .SAMPLE34  .SAMPLE340  .SAMPLE341  .SAMPLE342 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE343  .SAMPLE344  .SAMPLE345  .SAMPLE346  .SAMPLE347  .SAMPLE348 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE349   .SAMPLE35  .SAMPLE350  .SAMPLE351  .SAMPLE352  .SAMPLE353 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE354  .SAMPLE355  .SAMPLE356  .SAMPLE357  .SAMPLE358  .SAMPLE359 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE36  .SAMPLE360  .SAMPLE361  .SAMPLE362  .SAMPLE363  .SAMPLE364 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE365  .SAMPLE366  .SAMPLE367  .SAMPLE368  .SAMPLE369   .SAMPLE37 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE370  .SAMPLE371  .SAMPLE372  .SAMPLE373  .SAMPLE374  .SAMPLE375 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE376  .SAMPLE377  .SAMPLE378  .SAMPLE379   .SAMPLE38  .SAMPLE380 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE381  .SAMPLE382  .SAMPLE383  .SAMPLE384  .SAMPLE385  .SAMPLE386 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE387  .SAMPLE388  .SAMPLE389   .SAMPLE39  .SAMPLE390  .SAMPLE391 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE392  .SAMPLE393  .SAMPLE394  .SAMPLE395  .SAMPLE396  .SAMPLE397 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE398  .SAMPLE399    .SAMPLE4   .SAMPLE40  .SAMPLE400  .SAMPLE401 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE402  .SAMPLE403  .SAMPLE404  .SAMPLE405  .SAMPLE406  .SAMPLE407 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE408  .SAMPLE409   .SAMPLE41  .SAMPLE410  .SAMPLE411  .SAMPLE412 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE413  .SAMPLE414  .SAMPLE415  .SAMPLE416  .SAMPLE417  .SAMPLE418 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE419   .SAMPLE42  .SAMPLE420  .SAMPLE421  .SAMPLE422  .SAMPLE423 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE424  .SAMPLE425  .SAMPLE426  .SAMPLE427  .SAMPLE428  .SAMPLE429 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE43  .SAMPLE430  .SAMPLE431  .SAMPLE432  .SAMPLE433  .SAMPLE434 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE435  .SAMPLE436  .SAMPLE437  .SAMPLE438  .SAMPLE439   .SAMPLE44 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE440  .SAMPLE441  .SAMPLE442  .SAMPLE443  .SAMPLE444  .SAMPLE445 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE446  .SAMPLE447  .SAMPLE448  .SAMPLE449   .SAMPLE45  .SAMPLE450 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE451  .SAMPLE452  .SAMPLE453  .SAMPLE454  .SAMPLE455  .SAMPLE456 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE457  .SAMPLE458  .SAMPLE459   .SAMPLE46  .SAMPLE460  .SAMPLE461 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE462  .SAMPLE463  .SAMPLE464  .SAMPLE465  .SAMPLE466  .SAMPLE467 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE468  .SAMPLE469   .SAMPLE47  .SAMPLE470  .SAMPLE471  .SAMPLE472 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE473  .SAMPLE474  .SAMPLE475  .SAMPLE476  .SAMPLE477  .SAMPLE478 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE479   .SAMPLE48  .SAMPLE480  .SAMPLE481  .SAMPLE482  .SAMPLE483 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE484  .SAMPLE485  .SAMPLE486  .SAMPLE487  .SAMPLE488  .SAMPLE489 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE49  .SAMPLE490  .SAMPLE491  .SAMPLE492  .SAMPLE493  .SAMPLE494 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE495  .SAMPLE496  .SAMPLE497  .SAMPLE498  .SAMPLE499    .SAMPLE5 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE50  .SAMPLE500  .SAMPLE501  .SAMPLE502  .SAMPLE503  .SAMPLE504 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE505  .SAMPLE506  .SAMPLE507  .SAMPLE508  .SAMPLE509   .SAMPLE51 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE510  .SAMPLE511  .SAMPLE512  .SAMPLE513  .SAMPLE514  .SAMPLE515 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE516  .SAMPLE517  .SAMPLE518  .SAMPLE519   .SAMPLE52  .SAMPLE520 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE521  .SAMPLE522  .SAMPLE523  .SAMPLE524  .SAMPLE525  .SAMPLE526 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE527  .SAMPLE528  .SAMPLE529   .SAMPLE53  .SAMPLE530  .SAMPLE531 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE532  .SAMPLE533  .SAMPLE534  .SAMPLE535  .SAMPLE536  .SAMPLE537 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE538  .SAMPLE539   .SAMPLE54  .SAMPLE540  .SAMPLE541  .SAMPLE542 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE543  .SAMPLE544  .SAMPLE545  .SAMPLE546  .SAMPLE547  .SAMPLE548 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE549   .SAMPLE55  .SAMPLE550  .SAMPLE551  .SAMPLE552  .SAMPLE553 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE554  .SAMPLE555  .SAMPLE556  .SAMPLE557  .SAMPLE558  .SAMPLE559 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE56  .SAMPLE560  .SAMPLE561  .SAMPLE562  .SAMPLE563  .SAMPLE564 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE565  .SAMPLE566  .SAMPLE567  .SAMPLE568  .SAMPLE569   .SAMPLE57 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE570  .SAMPLE571  .SAMPLE572  .SAMPLE573  .SAMPLE574  .SAMPLE575 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE576  .SAMPLE577  .SAMPLE578  .SAMPLE579   .SAMPLE58  .SAMPLE580 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE581  .SAMPLE582  .SAMPLE583  .SAMPLE584  .SAMPLE585  .SAMPLE586 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE587  .SAMPLE588  .SAMPLE589   .SAMPLE59  .SAMPLE590  .SAMPLE591 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE592  .SAMPLE593  .SAMPLE594  .SAMPLE595  .SAMPLE596  .SAMPLE597 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE598  .SAMPLE599    .SAMPLE6   .SAMPLE60  .SAMPLE600  .SAMPLE601 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE602  .SAMPLE603  .SAMPLE604  .SAMPLE605  .SAMPLE606  .SAMPLE607 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE608  .SAMPLE609   .SAMPLE61  .SAMPLE610  .SAMPLE611  .SAMPLE612 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE613  .SAMPLE614  .SAMPLE615  .SAMPLE616  .SAMPLE617  .SAMPLE618 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE619   .SAMPLE62  .SAMPLE620  .SAMPLE621  .SAMPLE622  .SAMPLE623 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE624  .SAMPLE625  .SAMPLE626  .SAMPLE627  .SAMPLE628  .SAMPLE629 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE63  .SAMPLE630  .SAMPLE631  .SAMPLE632  .SAMPLE633  .SAMPLE634 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE635  .SAMPLE636  .SAMPLE637  .SAMPLE638  .SAMPLE639   .SAMPLE64 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE640  .SAMPLE641  .SAMPLE642  .SAMPLE643  .SAMPLE644  .SAMPLE645 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE646  .SAMPLE647  .SAMPLE648  .SAMPLE649   .SAMPLE65  .SAMPLE650 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE651  .SAMPLE652  .SAMPLE653  .SAMPLE654  .SAMPLE655  .SAMPLE656 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE657  .SAMPLE658  .SAMPLE659   .SAMPLE66  .SAMPLE660  .SAMPLE661 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE662  .SAMPLE663  .SAMPLE664  .SAMPLE665  .SAMPLE666  .SAMPLE667 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE668  .SAMPLE669   .SAMPLE67  .SAMPLE670  .SAMPLE671  .SAMPLE672 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE673  .SAMPLE674  .SAMPLE675  .SAMPLE676  .SAMPLE677  .SAMPLE678 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE679   .SAMPLE68  .SAMPLE680  .SAMPLE681  .SAMPLE682  .SAMPLE683 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE684  .SAMPLE685  .SAMPLE686  .SAMPLE687  .SAMPLE688  .SAMPLE689 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE69  .SAMPLE690  .SAMPLE691  .SAMPLE692  .SAMPLE693  .SAMPLE694 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE695  .SAMPLE696  .SAMPLE697  .SAMPLE698  .SAMPLE699    .SAMPLE7 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE70  .SAMPLE700  .SAMPLE701  .SAMPLE702  .SAMPLE703  .SAMPLE704 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE705  .SAMPLE706  .SAMPLE707  .SAMPLE708  .SAMPLE709   .SAMPLE71 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE710  .SAMPLE711  .SAMPLE712  .SAMPLE713  .SAMPLE714  .SAMPLE715 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE716  .SAMPLE717  .SAMPLE718  .SAMPLE719   .SAMPLE72  .SAMPLE720 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE721  .SAMPLE722  .SAMPLE723  .SAMPLE724  .SAMPLE725  .SAMPLE726 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE727  .SAMPLE728  .SAMPLE729   .SAMPLE73  .SAMPLE730  .SAMPLE731 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE732  .SAMPLE733  .SAMPLE734  .SAMPLE735  .SAMPLE736  .SAMPLE737 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE738  .SAMPLE739   .SAMPLE74  .SAMPLE740  .SAMPLE741  .SAMPLE742 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE743  .SAMPLE744  .SAMPLE745  .SAMPLE746  .SAMPLE747  .SAMPLE748 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE749   .SAMPLE75  .SAMPLE750  .SAMPLE751  .SAMPLE752  .SAMPLE753 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE754  .SAMPLE755  .SAMPLE756  .SAMPLE757  .SAMPLE758  .SAMPLE759 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE76  .SAMPLE760  .SAMPLE761  .SAMPLE762  .SAMPLE763  .SAMPLE764 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE765  .SAMPLE766  .SAMPLE767  .SAMPLE768  .SAMPLE769   .SAMPLE77 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE770  .SAMPLE771  .SAMPLE772  .SAMPLE773  .SAMPLE774  .SAMPLE775 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE776  .SAMPLE777  .SAMPLE778  .SAMPLE779   .SAMPLE78  .SAMPLE780 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE781  .SAMPLE782  .SAMPLE783  .SAMPLE784  .SAMPLE785  .SAMPLE786 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE787  .SAMPLE788  .SAMPLE789   .SAMPLE79  .SAMPLE790  .SAMPLE791 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE792  .SAMPLE793  .SAMPLE794  .SAMPLE795  .SAMPLE796  .SAMPLE797 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE798  .SAMPLE799    .SAMPLE8   .SAMPLE80  .SAMPLE800  .SAMPLE801 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE802  .SAMPLE803  .SAMPLE804  .SAMPLE805  .SAMPLE806  .SAMPLE807 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE808  .SAMPLE809   .SAMPLE81  .SAMPLE810  .SAMPLE811  .SAMPLE812 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE813  .SAMPLE814  .SAMPLE815  .SAMPLE816  .SAMPLE817  .SAMPLE818 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE819   .SAMPLE82  .SAMPLE820  .SAMPLE821  .SAMPLE822  .SAMPLE823 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE824  .SAMPLE825  .SAMPLE826  .SAMPLE827  .SAMPLE828  .SAMPLE829 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE83  .SAMPLE830  .SAMPLE831  .SAMPLE832  .SAMPLE833  .SAMPLE834 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE835  .SAMPLE836  .SAMPLE837  .SAMPLE838  .SAMPLE839   .SAMPLE84 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE840  .SAMPLE841  .SAMPLE842  .SAMPLE843  .SAMPLE844  .SAMPLE845 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE846  .SAMPLE847  .SAMPLE848  .SAMPLE849   .SAMPLE85  .SAMPLE850 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE851  .SAMPLE852  .SAMPLE853  .SAMPLE854  .SAMPLE855  .SAMPLE856 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE857  .SAMPLE858  .SAMPLE859   .SAMPLE86  .SAMPLE860  .SAMPLE861 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE862  .SAMPLE863  .SAMPLE864  .SAMPLE865  .SAMPLE866  .SAMPLE867 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE868  .SAMPLE869   .SAMPLE87  .SAMPLE870  .SAMPLE871  .SAMPLE872 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE873  .SAMPLE874  .SAMPLE875  .SAMPLE876  .SAMPLE877  .SAMPLE878 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE879   .SAMPLE88  .SAMPLE880  .SAMPLE881  .SAMPLE882  .SAMPLE883 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE884  .SAMPLE885  .SAMPLE886  .SAMPLE887  .SAMPLE888  .SAMPLE889 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE89  .SAMPLE890  .SAMPLE891  .SAMPLE892  .SAMPLE893  .SAMPLE894 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE895  .SAMPLE896  .SAMPLE897  .SAMPLE898  .SAMPLE899    .SAMPLE9 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE90  .SAMPLE900  .SAMPLE901  .SAMPLE902  .SAMPLE903  .SAMPLE904 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE905  .SAMPLE906  .SAMPLE907  .SAMPLE908  .SAMPLE909   .SAMPLE91 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE910  .SAMPLE911  .SAMPLE912  .SAMPLE913  .SAMPLE914  .SAMPLE915 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE916  .SAMPLE917  .SAMPLE918  .SAMPLE919   .SAMPLE92  .SAMPLE920 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE921  .SAMPLE922  .SAMPLE923  .SAMPLE924  .SAMPLE925  .SAMPLE926 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE927  .SAMPLE928  .SAMPLE929   .SAMPLE93  .SAMPLE930  .SAMPLE931 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE932  .SAMPLE933  .SAMPLE934  .SAMPLE935  .SAMPLE936  .SAMPLE937 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE938  .SAMPLE939   .SAMPLE94  .SAMPLE940  .SAMPLE941  .SAMPLE942 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE943  .SAMPLE944  .SAMPLE945  .SAMPLE946  .SAMPLE947  .SAMPLE948 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE949   .SAMPLE95  .SAMPLE950  .SAMPLE951  .SAMPLE952  .SAMPLE953 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE954  .SAMPLE955  .SAMPLE956  .SAMPLE957  .SAMPLE958  .SAMPLE959 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE96  .SAMPLE960  .SAMPLE961  .SAMPLE962  .SAMPLE963  .SAMPLE964 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE965  .SAMPLE966  .SAMPLE967  .SAMPLE968  .SAMPLE969   .SAMPLE97 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE970  .SAMPLE971  .SAMPLE972  .SAMPLE973  .SAMPLE974  .SAMPLE975 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE976  .SAMPLE977  .SAMPLE978  .SAMPLE979   .SAMPLE98  .SAMPLE980 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE981  .SAMPLE982  .SAMPLE983  .SAMPLE984  .SAMPLE985  .SAMPLE986 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE987  .SAMPLE988  .SAMPLE989   .SAMPLE99  .SAMPLE990  .SAMPLE991 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE992  .SAMPLE993  .SAMPLE994  .SAMPLE995  .SAMPLE996  .SAMPLE997 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE998  .SAMPLE999 
    ##           1           1

``` r
table(metadata$Site)
```

    ## 
    ## organic material 
    ##             2528

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
    ##      .A0001      .A0002      .A0003      .A0004      .A0005      .A0006 
    ##           1           1           1           1           1           1 
    ##      .A0007      .A0008      .A0009      .A0010      .A0011      .A0012 
    ##           1           1           1           1           1           1 
    ##      .A0013      .A0014      .A0015      .A0016      .A0017      .A0018 
    ##           1           1           1           1           1           1 
    ##      .A0019      .A0020      .A0021      .A0022      .A0023      .A0024 
    ##           1           1           1           1           1           1 
    ##      .A0025      .A0026      .A0027      .A0028      .A0029      .A0030 
    ##           1           1           1           1           1           1 
    ##      .A0031      .A0032      .A0033      .A0034      .A0035      .A0036 
    ##           1           1           1           1           1           1 
    ##      .A0037      .A0038      .A0039      .A0040      .A0041      .A0042 
    ##           1           1           1           1           1           1 
    ##      .A0043      .A0044      .A0045      .A0046      .A0047      .A0048 
    ##           1           1           1           1           1           1 
    ##      .A0049      .A0050      .A0051      .A0052      .A0053      .A0054 
    ##           1           1           1           1           1           1 
    ##      .A0055      .A0056      .A0057      .A0058      .A0059      .A0060 
    ##           1           1           1           1           1           1 
    ##      .A0061      .A0062      .A0063      .A0064      .A0065      .A0066 
    ##           1           1           1           1           1           1 
    ##      .A0067      .A0068      .A0069      .A0070      .A0071      .A0072 
    ##           1           1           1           1           1           1 
    ##      .A0073      .A0074      .A0075      .A0076      .A0077      .A0078 
    ##           1           1           1           1           1           1 
    ##      .A0079      .A0080      .A0081      .A0082      .A0083      .A0084 
    ##           1           1           1           1           1           1 
    ##      .A0085      .A0086      .A0087      .A0088      .A0089      .A0090 
    ##           1           1           1           1           1           1 
    ##      .A0091      .A0092      .A0093      .A0094      .A0095      .A0096 
    ##           1           1           1           1           1           1 
    ##      .A0097      .A0098      .A0099      .A0100      .A0101      .A0102 
    ##           1           1           1           1           1           1 
    ##      .A0103      .A0104      .A0105      .A0106      .A0107      .A0108 
    ##           1           1           1           1           1           1 
    ##      .A0109      .A0110      .A0111      .A0112      .A0113      .A0114 
    ##           1           1           1           1           1           1 
    ##      .A0115      .A0116      .A0117      .A0118      .A0119      .A0120 
    ##           1           1           1           1           1           1 
    ##      .A0121      .A0122      .A0123      .A0124      .A0125      .A0126 
    ##           1           1           1           1           1           1 
    ##      .A0127      .A0128      .A0129      .A0130      .A0131      .A0132 
    ##           1           1           1           1           1           1 
    ##      .A0133      .A0134      .A0135      .A0136      .A0137      .A0138 
    ##           1           1           1           1           1           1 
    ##      .A0139      .A0140      .A0141      .A0142      .A0143      .A0144 
    ##           1           1           1           1           1           1 
    ##      .A0145      .A0146      .A0147      .A0148      .A0149      .A0150 
    ##           1           1           1           1           1           1 
    ##      .A0151      .A0152      .A0153      .A0154      .A0155      .A0156 
    ##           1           1           1           1           1           1 
    ##      .A0157      .A0158      .A0159      .A0160      .A0161      .A0162 
    ##           1           1           1           1           1           1 
    ##      .A0163      .A0164      .A0165      .A0166      .A0167      .A0168 
    ##           1           1           1           1           1           1 
    ##      .A0169      .A0170      .A0171      .A0172      .A0173      .A0174 
    ##           1           1           1           1           1           1 
    ##      .A0175      .A0176      .A0177      .A0178      .A0179      .A0180 
    ##           1           1           1           1           1           1 
    ##      .A0181      .A0182      .A0183      .A0184      .A0185      .A0186 
    ##           1           1           1           1           1           1 
    ##      .A0187      .A0188      .A0189      .A0190      .A0191      .A0192 
    ##           1           1           1           1           1           1 
    ##      .A0193      .A0194      .A0195      .A0196      .A0197      .A0198 
    ##           1           1           1           1           1           1 
    ##      .A0199      .A0200      .A0201      .A0202      .A0203      .A0204 
    ##           1           1           1           1           1           1 
    ##      .A0205      .A0206      .A0207      .A0208      .A0209      .A0210 
    ##           1           1           1           1           1           1 
    ##      .A0211      .A0212      .A0213      .A0214      .A0215      .A0216 
    ##           1           1           1           1           1           1 
    ##      .A0217      .A0218      .A0219      .A0220      .A0221      .A0222 
    ##           1           1           1           1           1           1 
    ##      .A0223      .A0224      .A0225      .A0226      .A0227      .A0228 
    ##           1           1           1           1           1           1 
    ##      .A0229      .A0230      .A0231      .A0232      .A0233      .A0234 
    ##           1           1           1           1           1           1 
    ##      .A0235      .A0236      .A0237      .A0238      .A0239      .A0240 
    ##           1           1           1           1           1           1 
    ##      .A0241      .A0242      .A0243      .A0244      .A0245      .A0246 
    ##           1           1           1           1           1           1 
    ##      .A0247      .A0248      .A0249      .A0250      .A0251      .A0252 
    ##           1           1           1           1           1           1 
    ##      .A0253      .A0254      .A0255      .A0256      .A0257      .A0258 
    ##           1           1           1           1           1           1 
    ##      .A0259      .A0260      .A0261      .A0262      .A0263      .A0264 
    ##           1           1           1           1           1           1 
    ##      .A0265      .A0266      .A0267      .A0268      .A0269      .A0270 
    ##           1           1           1           1           1           1 
    ##      .A0271      .A0272      .A0273      .A0274      .A0275      .A0276 
    ##           1           1           1           1           1           1 
    ##      .A0277      .A0278      .A0279      .A0280      .A0281      .A0282 
    ##           1           1           1           1           1           1 
    ##      .A0283      .A0284      .A0285      .A0286      .A0287      .A0288 
    ##           1           1           1           1           1           1 
    ##      .A0289      .A0290      .A0291      .A0292      .A0293      .A0294 
    ##           1           1           1           1           1           1 
    ##      .A0295      .A0296      .A0297      .A0298      .A0299      .A0300 
    ##           1           1           1           1           1           1 
    ##      .A0301      .A0302      .A0303      .A0304      .A0305      .A0306 
    ##           1           1           1           1           1           1 
    ##      .A0307      .A0308      .A0309      .A0310      .A0311      .A0312 
    ##           1           1           1           1           1           1 
    ##      .A0313      .A0314      .A0315      .A0316      .A0317      .A0318 
    ##           1           1           1           1           1           1 
    ##      .A0319      .A0320      .A0321      .A0322      .A0323      .A0324 
    ##           1           1           1           1           1           1 
    ##      .A0325      .A0326      .A0327      .A0328      .A0329      .A0330 
    ##           1           1           1           1           1           1 
    ##      .A0331      .A0332      .A0333      .A0334      .A0335      .A0336 
    ##           1           1           1           1           1           1 
    ##      .A0337      .A0338      .A0339      .A0340      .A0341      .A0342 
    ##           1           1           1           1           1           1 
    ##      .A0343      .A0344      .A0345      .A0346      .A0347      .A0348 
    ##           1           1           1           1           1           1 
    ##      .A0349      .A0350      .A0351      .A0352      .A0353      .A0354 
    ##           1           1           1           1           1           1 
    ##      .A0355      .A0356      .A0357      .A0358      .A0359      .A0360 
    ##           1           1           1           1           1           1 
    ##      .A0361      .A0362      .A0363      .A0364      .A0365      .A0366 
    ##           1           1           1           1           1           1 
    ##      .A0367      .A0368      .A0369      .A0370      .A0371      .A0372 
    ##           1           1           1           1           1           1 
    ##      .A0373      .A0374      .A0375      .A0376      .A0377      .A0378 
    ##           1           1           1           1           1           1 
    ##      .A0379      .A0380      .A0381      .A0382      .A0383      .A0384 
    ##           1           1           1           1           1           1 
    ##      .A0385      .A0386      .A0387      .A0388      .A0389      .A0390 
    ##           1           1           1           1           1           1 
    ##      .A0391      .A0392      .A0393      .A0394      .A0395      .A0396 
    ##           1           1           1           1           1           1 
    ##      .A0397      .A0398      .A0399      .A0400      .A0401      .A0402 
    ##           1           1           1           1           1           1 
    ##      .A0403      .A0404      .A0405      .A0406      .A0407      .A0408 
    ##           1           1           1           1           1           1 
    ##      .A0409      .A0410      .A0411      .A0412      .A0413      .A0414 
    ##           1           1           1           1           1           1 
    ##      .A0415      .A0416      .A0417      .A0418      .A0419      .A0420 
    ##           1           1           1           1           1           1 
    ##      .A0421      .A0422      .A0423      .A0424      .A0425      .A0426 
    ##           1           1           1           1           1           1 
    ##      .A0427      .A0428      .A0429      .A0430      .A0431      .A0432 
    ##           1           1           1           1           1           1 
    ##      .A0433      .A0434      .A0435      .A0436      .A0437      .A0438 
    ##           1           1           1           1           1           1 
    ##      .A0439      .A0440      .A0441      .A0442      .A0443      .A0444 
    ##           1           1           1           1           1           1 
    ##      .A0445      .A0446      .A0447      .A0448      .A0449      .A0450 
    ##           1           1           1           1           1           1 
    ##      .A0451      .A0452      .A0453      .A0454      .A0455      .A0456 
    ##           1           1           1           1           1           1 
    ##      .A0457      .A0458      .A0459      .A0460      .A0461      .A0462 
    ##           1           1           1           1           1           1 
    ##      .A0463      .A0464      .A0465      .A0466      .A0467      .A0468 
    ##           1           1           1           1           1           1 
    ##      .A0469      .A0470      .A0471      .A0472      .A0473      .A0474 
    ##           1           1           1           1           1           1 
    ##      .A0475      .A0476      .A0477      .A0478      .A0479      .A0480 
    ##           1           1           1           1           1           1 
    ##      .A0481      .A0482      .A0483      .A0484      .A0485      .A0486 
    ##           1           1           1           1           1           1 
    ##      .A0487      .A0488      .A0489      .A0490      .A0491      .A0492 
    ##           1           1           1           1           1           1 
    ##      .A0493      .A0494      .A0495      .A0496      .A0497      .A0498 
    ##           1           1           1           1           1           1 
    ##      .A0499      .A0500      .A0501      .A0502      .A0503      .A0504 
    ##           1           1           1           1           1           1 
    ##      .A0505      .A0506      .A0507      .A0508      .A0509      .A0510 
    ##           1           1           1           1           1           1 
    ##      .A0511      .A0512      .A0513      .A0514      .A0515      .A0516 
    ##           1           1           1           1           1           1 
    ##      .A0517      .A0518      .A0519      .A0520      .A0521      .A0522 
    ##           1           1           1           1           1           1 
    ##      .A0523      .A0524      .A0525      .A0526      .A0527      .A0528 
    ##           1           1           1           1           1           1 
    ##      .A0529      .A0530      .A0531      .A0532      .A0533      .A0534 
    ##           1           1           1           1           1           1 
    ##      .A0535      .A0536      .A0537      .A0538      .A0539      .A0540 
    ##           1           1           1           1           1           1 
    ##      .A0541      .A0542      .A0543      .A0544      .A0545      .A0546 
    ##           1           1           1           1           1           1 
    ##      .A0547      .A0548      .A0549      .A0550      .A0551      .A0552 
    ##           1           1           1           1           1           1 
    ##      .A0553      .A0554      .A0555      .A0556      .A0557      .A0558 
    ##           1           1           1           1           1           1 
    ##      .A0559      .A0560      .A0561      .A0562      .A0563      .A0564 
    ##           1           1           1           1           1           1 
    ##      .A0565      .A0566      .A0567      .A0568      .A0569      .A0570 
    ##           1           1           1           1           1           1 
    ##      .A0571      .A0572      .A0573      .A0574      .A0575      .A0576 
    ##           1           1           1           1           1           1 
    ##      .A0577      .A0578      .A0579      .A0580      .A0581      .A0582 
    ##           1           1           1           1           1           1 
    ##      .A0583      .A0584      .A0585      .A0586      .A0587      .A0588 
    ##           1           1           1           1           1           1 
    ##      .A0589      .A0590      .A0591      .A0592      .A0593      .A0594 
    ##           1           1           1           1           1           1 
    ##      .A0595      .A0596      .A0597      .A0598      .A0599      .A0600 
    ##           1           1           1           1           1           1 
    ##      .A0601      .A0602      .A0603      .A0604      .A0605      .A0606 
    ##           1           1           1           1           1           1 
    ##      .A0607      .A0608      .A0609      .A0610      .A0611      .A0612 
    ##           1           1           1           1           1           1 
    ##      .A0613      .A0614      .A0615      .A0616      .A0617      .A0618 
    ##           1           1           1           1           1           1 
    ##      .A0619      .A0620      .A0621      .A0622      .A0623      .A0624 
    ##           1           1           1           1           1           1 
    ##      .A0625      .A0626      .A0627      .A0628      .A0629      .A0630 
    ##           1           1           1           1           1           1 
    ##      .A0631      .A0632      .A0633      .A0634      .A0635      .A0636 
    ##           1           1           1           1           1           1 
    ##      .A0637         .S1        .S10        .S11        .S12        .S13 
    ##           1           1           1           1           1           1 
    ##        .S14        .S15        .S16        .S17        .S18        .S19 
    ##           1           1           1           1           1           1 
    ##         .S2        .S20        .S21        .S22        .S23        .S24 
    ##           1           1           1           1           1           1 
    ##         .S3         .S4         .S5         .S6         .S7         .S8 
    ##           1           1           1           1           1           1 
    ##         .S9    .SAMPLE1   .SAMPLE10  .SAMPLE100 .SAMPLE1000 .SAMPLE1001 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1002 .SAMPLE1003 .SAMPLE1004 .SAMPLE1005 .SAMPLE1006 .SAMPLE1007 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1008 .SAMPLE1009  .SAMPLE101 .SAMPLE1010 .SAMPLE1011 .SAMPLE1012 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1013 .SAMPLE1014 .SAMPLE1015 .SAMPLE1016 .SAMPLE1017 .SAMPLE1018 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1019  .SAMPLE102 .SAMPLE1020 .SAMPLE1021 .SAMPLE1022 .SAMPLE1023 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1024 .SAMPLE1025 .SAMPLE1026 .SAMPLE1027 .SAMPLE1028 .SAMPLE1029 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE103 .SAMPLE1030 .SAMPLE1031 .SAMPLE1032 .SAMPLE1033 .SAMPLE1034 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1035 .SAMPLE1036 .SAMPLE1037 .SAMPLE1038 .SAMPLE1039  .SAMPLE104 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1040 .SAMPLE1041 .SAMPLE1042 .SAMPLE1043 .SAMPLE1044 .SAMPLE1045 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1046 .SAMPLE1047 .SAMPLE1048 .SAMPLE1049  .SAMPLE105 .SAMPLE1050 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1051 .SAMPLE1052 .SAMPLE1053 .SAMPLE1054 .SAMPLE1055 .SAMPLE1056 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1057 .SAMPLE1058 .SAMPLE1059  .SAMPLE106 .SAMPLE1060 .SAMPLE1061 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1062 .SAMPLE1063 .SAMPLE1064 .SAMPLE1065 .SAMPLE1066 .SAMPLE1067 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1068 .SAMPLE1069  .SAMPLE107 .SAMPLE1070 .SAMPLE1071 .SAMPLE1072 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1073 .SAMPLE1074 .SAMPLE1075 .SAMPLE1076 .SAMPLE1077 .SAMPLE1078 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1079  .SAMPLE108 .SAMPLE1080 .SAMPLE1081 .SAMPLE1082 .SAMPLE1083 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1084 .SAMPLE1085 .SAMPLE1086 .SAMPLE1087 .SAMPLE1088 .SAMPLE1089 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE109 .SAMPLE1090 .SAMPLE1091 .SAMPLE1092 .SAMPLE1093 .SAMPLE1094 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1095 .SAMPLE1096 .SAMPLE1097 .SAMPLE1098 .SAMPLE1099   .SAMPLE11 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE110 .SAMPLE1100 .SAMPLE1101 .SAMPLE1102 .SAMPLE1103 .SAMPLE1104 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1105 .SAMPLE1106 .SAMPLE1107 .SAMPLE1108 .SAMPLE1109  .SAMPLE111 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1110 .SAMPLE1111 .SAMPLE1112 .SAMPLE1113 .SAMPLE1114 .SAMPLE1115 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1116 .SAMPLE1117 .SAMPLE1118 .SAMPLE1119  .SAMPLE112 .SAMPLE1120 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1121 .SAMPLE1122 .SAMPLE1123 .SAMPLE1124 .SAMPLE1125 .SAMPLE1126 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1127 .SAMPLE1128 .SAMPLE1129  .SAMPLE113 .SAMPLE1130 .SAMPLE1131 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1132 .SAMPLE1133 .SAMPLE1134 .SAMPLE1135 .SAMPLE1136 .SAMPLE1137 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1138 .SAMPLE1139  .SAMPLE114 .SAMPLE1140 .SAMPLE1141 .SAMPLE1142 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1143 .SAMPLE1144 .SAMPLE1145 .SAMPLE1146 .SAMPLE1147 .SAMPLE1148 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1149  .SAMPLE115 .SAMPLE1150 .SAMPLE1151 .SAMPLE1152 .SAMPLE1153 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1154 .SAMPLE1155 .SAMPLE1156 .SAMPLE1157 .SAMPLE1158 .SAMPLE1159 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE116 .SAMPLE1160 .SAMPLE1161 .SAMPLE1162 .SAMPLE1163 .SAMPLE1164 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1165 .SAMPLE1166 .SAMPLE1167 .SAMPLE1168 .SAMPLE1169  .SAMPLE117 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1170 .SAMPLE1171 .SAMPLE1172 .SAMPLE1173 .SAMPLE1174 .SAMPLE1175 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1176 .SAMPLE1177 .SAMPLE1178 .SAMPLE1179  .SAMPLE118 .SAMPLE1180 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1181 .SAMPLE1182 .SAMPLE1183 .SAMPLE1184 .SAMPLE1185 .SAMPLE1186 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1187 .SAMPLE1188 .SAMPLE1189  .SAMPLE119 .SAMPLE1190 .SAMPLE1191 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1192 .SAMPLE1193 .SAMPLE1194 .SAMPLE1195 .SAMPLE1196 .SAMPLE1197 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1198 .SAMPLE1199   .SAMPLE12  .SAMPLE120 .SAMPLE1200 .SAMPLE1201 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1202 .SAMPLE1203 .SAMPLE1204 .SAMPLE1205 .SAMPLE1206 .SAMPLE1207 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1208 .SAMPLE1209  .SAMPLE121 .SAMPLE1210 .SAMPLE1211 .SAMPLE1212 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1213 .SAMPLE1214 .SAMPLE1215 .SAMPLE1216 .SAMPLE1217 .SAMPLE1218 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1219  .SAMPLE122 .SAMPLE1220 .SAMPLE1221 .SAMPLE1222 .SAMPLE1223 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1224 .SAMPLE1225 .SAMPLE1226 .SAMPLE1227 .SAMPLE1228 .SAMPLE1229 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE123 .SAMPLE1230 .SAMPLE1231 .SAMPLE1232 .SAMPLE1233 .SAMPLE1234 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1235 .SAMPLE1236 .SAMPLE1237 .SAMPLE1238 .SAMPLE1239  .SAMPLE124 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1240 .SAMPLE1241 .SAMPLE1242 .SAMPLE1243 .SAMPLE1244 .SAMPLE1245 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1246 .SAMPLE1247 .SAMPLE1248 .SAMPLE1249  .SAMPLE125 .SAMPLE1250 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1251 .SAMPLE1252 .SAMPLE1253 .SAMPLE1254 .SAMPLE1255 .SAMPLE1256 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1257 .SAMPLE1258 .SAMPLE1259  .SAMPLE126 .SAMPLE1260 .SAMPLE1261 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1262 .SAMPLE1263 .SAMPLE1264 .SAMPLE1265 .SAMPLE1266 .SAMPLE1267 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1268 .SAMPLE1269  .SAMPLE127 .SAMPLE1270 .SAMPLE1271 .SAMPLE1272 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1273 .SAMPLE1274 .SAMPLE1275 .SAMPLE1276 .SAMPLE1277 .SAMPLE1278 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1279  .SAMPLE128 .SAMPLE1280 .SAMPLE1281 .SAMPLE1282 .SAMPLE1283 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1284 .SAMPLE1285 .SAMPLE1286 .SAMPLE1287 .SAMPLE1288 .SAMPLE1289 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE129 .SAMPLE1290 .SAMPLE1291 .SAMPLE1292 .SAMPLE1293 .SAMPLE1294 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1295 .SAMPLE1296 .SAMPLE1297 .SAMPLE1298 .SAMPLE1299   .SAMPLE13 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE130 .SAMPLE1300 .SAMPLE1301 .SAMPLE1302 .SAMPLE1303 .SAMPLE1304 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1305 .SAMPLE1306 .SAMPLE1307 .SAMPLE1308 .SAMPLE1309  .SAMPLE131 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1310 .SAMPLE1311 .SAMPLE1312 .SAMPLE1313 .SAMPLE1314 .SAMPLE1315 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1316 .SAMPLE1317 .SAMPLE1318 .SAMPLE1319  .SAMPLE132 .SAMPLE1320 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1321 .SAMPLE1322 .SAMPLE1323 .SAMPLE1324 .SAMPLE1325 .SAMPLE1326 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1327 .SAMPLE1328 .SAMPLE1329  .SAMPLE133 .SAMPLE1330 .SAMPLE1331 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1332 .SAMPLE1333 .SAMPLE1334 .SAMPLE1335 .SAMPLE1336 .SAMPLE1337 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1338 .SAMPLE1339  .SAMPLE134 .SAMPLE1340 .SAMPLE1341 .SAMPLE1342 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1343 .SAMPLE1344 .SAMPLE1345 .SAMPLE1346 .SAMPLE1347 .SAMPLE1348 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1349  .SAMPLE135 .SAMPLE1350 .SAMPLE1351 .SAMPLE1352 .SAMPLE1353 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1354 .SAMPLE1355 .SAMPLE1356 .SAMPLE1357 .SAMPLE1358 .SAMPLE1359 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE136 .SAMPLE1360 .SAMPLE1361 .SAMPLE1362 .SAMPLE1363 .SAMPLE1364 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1365 .SAMPLE1366 .SAMPLE1367 .SAMPLE1368 .SAMPLE1369  .SAMPLE137 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1370 .SAMPLE1371 .SAMPLE1372 .SAMPLE1373 .SAMPLE1374 .SAMPLE1375 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1376 .SAMPLE1377 .SAMPLE1378 .SAMPLE1379  .SAMPLE138 .SAMPLE1380 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1381 .SAMPLE1382 .SAMPLE1383 .SAMPLE1384 .SAMPLE1385 .SAMPLE1386 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1387 .SAMPLE1388 .SAMPLE1389  .SAMPLE139 .SAMPLE1390 .SAMPLE1391 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1392 .SAMPLE1393 .SAMPLE1394 .SAMPLE1395 .SAMPLE1396 .SAMPLE1397 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1398 .SAMPLE1399   .SAMPLE14  .SAMPLE140 .SAMPLE1400 .SAMPLE1401 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1402 .SAMPLE1403 .SAMPLE1404 .SAMPLE1405 .SAMPLE1406 .SAMPLE1407 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1408 .SAMPLE1409  .SAMPLE141 .SAMPLE1410 .SAMPLE1411 .SAMPLE1412 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1413 .SAMPLE1414 .SAMPLE1415 .SAMPLE1416 .SAMPLE1417 .SAMPLE1418 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1419  .SAMPLE142 .SAMPLE1420 .SAMPLE1421 .SAMPLE1422 .SAMPLE1423 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1424 .SAMPLE1425 .SAMPLE1426 .SAMPLE1427 .SAMPLE1428 .SAMPLE1429 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE143 .SAMPLE1430 .SAMPLE1431 .SAMPLE1432 .SAMPLE1433 .SAMPLE1434 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1435 .SAMPLE1436 .SAMPLE1437 .SAMPLE1438 .SAMPLE1439  .SAMPLE144 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1440 .SAMPLE1441 .SAMPLE1442 .SAMPLE1443 .SAMPLE1444 .SAMPLE1445 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1446 .SAMPLE1447 .SAMPLE1448 .SAMPLE1449  .SAMPLE145 .SAMPLE1450 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1451 .SAMPLE1452 .SAMPLE1453 .SAMPLE1454 .SAMPLE1455 .SAMPLE1456 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1457 .SAMPLE1458 .SAMPLE1459  .SAMPLE146 .SAMPLE1460 .SAMPLE1461 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1462 .SAMPLE1463 .SAMPLE1464 .SAMPLE1465 .SAMPLE1466 .SAMPLE1467 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1468 .SAMPLE1469  .SAMPLE147 .SAMPLE1470 .SAMPLE1471 .SAMPLE1472 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1473 .SAMPLE1474 .SAMPLE1475 .SAMPLE1476 .SAMPLE1477 .SAMPLE1478 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1479  .SAMPLE148 .SAMPLE1480 .SAMPLE1481 .SAMPLE1482 .SAMPLE1483 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1484 .SAMPLE1485 .SAMPLE1486 .SAMPLE1487 .SAMPLE1488 .SAMPLE1489 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE149 .SAMPLE1490 .SAMPLE1491 .SAMPLE1492 .SAMPLE1493 .SAMPLE1494 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1495 .SAMPLE1496 .SAMPLE1497 .SAMPLE1498 .SAMPLE1499   .SAMPLE15 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE150 .SAMPLE1500 .SAMPLE1501 .SAMPLE1502 .SAMPLE1503 .SAMPLE1504 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1505 .SAMPLE1506 .SAMPLE1507 .SAMPLE1508 .SAMPLE1509  .SAMPLE151 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1510 .SAMPLE1511 .SAMPLE1512 .SAMPLE1513 .SAMPLE1514 .SAMPLE1515 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1516 .SAMPLE1517 .SAMPLE1518 .SAMPLE1519  .SAMPLE152 .SAMPLE1520 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1521 .SAMPLE1522 .SAMPLE1523 .SAMPLE1524 .SAMPLE1525 .SAMPLE1526 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1527 .SAMPLE1528 .SAMPLE1529  .SAMPLE153 .SAMPLE1530 .SAMPLE1531 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1532 .SAMPLE1533 .SAMPLE1534 .SAMPLE1535 .SAMPLE1536 .SAMPLE1537 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1538 .SAMPLE1539  .SAMPLE154 .SAMPLE1540 .SAMPLE1541 .SAMPLE1542 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1543 .SAMPLE1544 .SAMPLE1545 .SAMPLE1546 .SAMPLE1547 .SAMPLE1548 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1549  .SAMPLE155 .SAMPLE1550 .SAMPLE1551 .SAMPLE1552 .SAMPLE1553 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1554 .SAMPLE1555 .SAMPLE1556 .SAMPLE1557 .SAMPLE1558 .SAMPLE1559 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE156 .SAMPLE1560 .SAMPLE1561 .SAMPLE1562 .SAMPLE1563 .SAMPLE1564 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1565 .SAMPLE1566 .SAMPLE1567 .SAMPLE1568 .SAMPLE1569  .SAMPLE157 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1570 .SAMPLE1571 .SAMPLE1572 .SAMPLE1573 .SAMPLE1574 .SAMPLE1575 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1576 .SAMPLE1577 .SAMPLE1578 .SAMPLE1579  .SAMPLE158 .SAMPLE1580 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1581 .SAMPLE1582 .SAMPLE1583 .SAMPLE1584 .SAMPLE1585 .SAMPLE1586 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1587 .SAMPLE1588 .SAMPLE1589  .SAMPLE159 .SAMPLE1590 .SAMPLE1591 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1592 .SAMPLE1593 .SAMPLE1594 .SAMPLE1595 .SAMPLE1596 .SAMPLE1597 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1598 .SAMPLE1599   .SAMPLE16  .SAMPLE160 .SAMPLE1600 .SAMPLE1601 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1602 .SAMPLE1603 .SAMPLE1604 .SAMPLE1605 .SAMPLE1606 .SAMPLE1607 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1608 .SAMPLE1609  .SAMPLE161 .SAMPLE1610 .SAMPLE1611 .SAMPLE1612 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1613 .SAMPLE1614 .SAMPLE1615 .SAMPLE1616 .SAMPLE1617 .SAMPLE1618 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1619  .SAMPLE162 .SAMPLE1620 .SAMPLE1621 .SAMPLE1622 .SAMPLE1623 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1624 .SAMPLE1625 .SAMPLE1626 .SAMPLE1627 .SAMPLE1628 .SAMPLE1629 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE163 .SAMPLE1630 .SAMPLE1631 .SAMPLE1632 .SAMPLE1633 .SAMPLE1634 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1635 .SAMPLE1636 .SAMPLE1637 .SAMPLE1638 .SAMPLE1639  .SAMPLE164 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1640 .SAMPLE1641 .SAMPLE1642 .SAMPLE1643 .SAMPLE1644 .SAMPLE1645 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1646 .SAMPLE1647 .SAMPLE1648 .SAMPLE1649  .SAMPLE165 .SAMPLE1650 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1651 .SAMPLE1652 .SAMPLE1653 .SAMPLE1654 .SAMPLE1655 .SAMPLE1656 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1657 .SAMPLE1658 .SAMPLE1659  .SAMPLE166 .SAMPLE1660 .SAMPLE1661 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1662 .SAMPLE1663 .SAMPLE1664 .SAMPLE1665 .SAMPLE1666 .SAMPLE1667 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1668 .SAMPLE1669  .SAMPLE167 .SAMPLE1670 .SAMPLE1671 .SAMPLE1672 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1673 .SAMPLE1674 .SAMPLE1675 .SAMPLE1676 .SAMPLE1677 .SAMPLE1678 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1679  .SAMPLE168 .SAMPLE1680 .SAMPLE1681 .SAMPLE1682 .SAMPLE1683 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1684 .SAMPLE1685 .SAMPLE1686 .SAMPLE1687 .SAMPLE1688 .SAMPLE1689 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE169 .SAMPLE1690 .SAMPLE1691 .SAMPLE1692 .SAMPLE1693 .SAMPLE1694 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1695 .SAMPLE1696 .SAMPLE1697 .SAMPLE1698 .SAMPLE1699   .SAMPLE17 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE170 .SAMPLE1700 .SAMPLE1701 .SAMPLE1702 .SAMPLE1703 .SAMPLE1704 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1705 .SAMPLE1706 .SAMPLE1707 .SAMPLE1708 .SAMPLE1709  .SAMPLE171 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1710 .SAMPLE1711 .SAMPLE1712 .SAMPLE1713 .SAMPLE1714 .SAMPLE1715 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1716 .SAMPLE1717 .SAMPLE1718 .SAMPLE1719  .SAMPLE172 .SAMPLE1720 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1721 .SAMPLE1722 .SAMPLE1723 .SAMPLE1724 .SAMPLE1725 .SAMPLE1726 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1727 .SAMPLE1728 .SAMPLE1729  .SAMPLE173 .SAMPLE1730 .SAMPLE1731 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1732 .SAMPLE1733 .SAMPLE1734 .SAMPLE1735 .SAMPLE1736 .SAMPLE1737 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1738 .SAMPLE1739  .SAMPLE174 .SAMPLE1740 .SAMPLE1741 .SAMPLE1742 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1743 .SAMPLE1744 .SAMPLE1745 .SAMPLE1746 .SAMPLE1747 .SAMPLE1748 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1749  .SAMPLE175 .SAMPLE1750 .SAMPLE1751 .SAMPLE1752 .SAMPLE1753 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1754 .SAMPLE1755 .SAMPLE1756 .SAMPLE1757 .SAMPLE1758 .SAMPLE1759 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE176 .SAMPLE1760 .SAMPLE1761 .SAMPLE1762 .SAMPLE1763 .SAMPLE1764 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1765 .SAMPLE1766 .SAMPLE1767 .SAMPLE1768 .SAMPLE1769  .SAMPLE177 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1770 .SAMPLE1771 .SAMPLE1772 .SAMPLE1773 .SAMPLE1774 .SAMPLE1775 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1776 .SAMPLE1777 .SAMPLE1778 .SAMPLE1779  .SAMPLE178 .SAMPLE1780 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1781 .SAMPLE1782 .SAMPLE1783 .SAMPLE1784 .SAMPLE1785 .SAMPLE1786 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1787 .SAMPLE1788 .SAMPLE1789  .SAMPLE179 .SAMPLE1790 .SAMPLE1791 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1792 .SAMPLE1793 .SAMPLE1794 .SAMPLE1795 .SAMPLE1796 .SAMPLE1797 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1798 .SAMPLE1799   .SAMPLE18  .SAMPLE180 .SAMPLE1800 .SAMPLE1801 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1802 .SAMPLE1803 .SAMPLE1804 .SAMPLE1805 .SAMPLE1806 .SAMPLE1807 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1808 .SAMPLE1809  .SAMPLE181 .SAMPLE1810 .SAMPLE1811 .SAMPLE1812 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1813 .SAMPLE1814 .SAMPLE1815 .SAMPLE1816 .SAMPLE1817 .SAMPLE1818 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1819  .SAMPLE182 .SAMPLE1820 .SAMPLE1821 .SAMPLE1822 .SAMPLE1823 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1824 .SAMPLE1825 .SAMPLE1826 .SAMPLE1827 .SAMPLE1828 .SAMPLE1829 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE183 .SAMPLE1830 .SAMPLE1831 .SAMPLE1832 .SAMPLE1833 .SAMPLE1834 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1835 .SAMPLE1836 .SAMPLE1837 .SAMPLE1838 .SAMPLE1839  .SAMPLE184 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1840 .SAMPLE1841 .SAMPLE1842 .SAMPLE1843 .SAMPLE1844 .SAMPLE1845 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1846 .SAMPLE1847 .SAMPLE1848 .SAMPLE1849  .SAMPLE185 .SAMPLE1850 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1851 .SAMPLE1852 .SAMPLE1853 .SAMPLE1854 .SAMPLE1855 .SAMPLE1856 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1857 .SAMPLE1858 .SAMPLE1859  .SAMPLE186 .SAMPLE1860 .SAMPLE1861 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1862 .SAMPLE1863 .SAMPLE1864 .SAMPLE1865 .SAMPLE1866 .SAMPLE1867 
    ##           1           1           1           1           1           1 
    ## .SAMPLE1868 .SAMPLE1869  .SAMPLE187 .SAMPLE1870 .SAMPLE1871 .SAMPLE1872 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE188  .SAMPLE189   .SAMPLE19  .SAMPLE190  .SAMPLE191  .SAMPLE192 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE193  .SAMPLE194  .SAMPLE195  .SAMPLE196  .SAMPLE197  .SAMPLE198 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE199    .SAMPLE2   .SAMPLE20  .SAMPLE200  .SAMPLE201  .SAMPLE202 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE203  .SAMPLE204  .SAMPLE205  .SAMPLE206  .SAMPLE207  .SAMPLE208 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE209   .SAMPLE21  .SAMPLE210  .SAMPLE211  .SAMPLE212  .SAMPLE213 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE214  .SAMPLE215  .SAMPLE216  .SAMPLE217  .SAMPLE218  .SAMPLE219 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE22  .SAMPLE220  .SAMPLE221  .SAMPLE222  .SAMPLE223  .SAMPLE224 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE225  .SAMPLE226  .SAMPLE227  .SAMPLE228  .SAMPLE229   .SAMPLE23 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE230  .SAMPLE231  .SAMPLE232  .SAMPLE233  .SAMPLE234  .SAMPLE235 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE236  .SAMPLE237  .SAMPLE238  .SAMPLE239   .SAMPLE24  .SAMPLE240 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE241  .SAMPLE242  .SAMPLE243  .SAMPLE244  .SAMPLE245  .SAMPLE246 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE247  .SAMPLE248  .SAMPLE249   .SAMPLE25  .SAMPLE250  .SAMPLE251 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE252  .SAMPLE253  .SAMPLE254  .SAMPLE255  .SAMPLE256  .SAMPLE257 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE258  .SAMPLE259   .SAMPLE26  .SAMPLE260  .SAMPLE261  .SAMPLE262 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE263  .SAMPLE264  .SAMPLE265  .SAMPLE266  .SAMPLE267  .SAMPLE268 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE269   .SAMPLE27  .SAMPLE270  .SAMPLE271  .SAMPLE272  .SAMPLE273 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE274  .SAMPLE275  .SAMPLE276  .SAMPLE277  .SAMPLE278  .SAMPLE279 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE28  .SAMPLE280  .SAMPLE281  .SAMPLE282  .SAMPLE283  .SAMPLE284 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE285  .SAMPLE286  .SAMPLE287  .SAMPLE288  .SAMPLE289   .SAMPLE29 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE290  .SAMPLE291  .SAMPLE292  .SAMPLE293  .SAMPLE294  .SAMPLE295 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE296  .SAMPLE297  .SAMPLE298  .SAMPLE299    .SAMPLE3   .SAMPLE30 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE300  .SAMPLE301  .SAMPLE302  .SAMPLE303  .SAMPLE304  .SAMPLE305 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE306  .SAMPLE307  .SAMPLE308  .SAMPLE309   .SAMPLE31  .SAMPLE310 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE311  .SAMPLE312  .SAMPLE313  .SAMPLE314  .SAMPLE315  .SAMPLE316 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE318  .SAMPLE319   .SAMPLE32  .SAMPLE321  .SAMPLE322  .SAMPLE324 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE325  .SAMPLE327  .SAMPLE329   .SAMPLE33  .SAMPLE330  .SAMPLE331 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE332  .SAMPLE333  .SAMPLE334  .SAMPLE335  .SAMPLE336  .SAMPLE337 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE338  .SAMPLE339   .SAMPLE34  .SAMPLE340  .SAMPLE341  .SAMPLE342 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE343  .SAMPLE344  .SAMPLE345  .SAMPLE346  .SAMPLE347  .SAMPLE348 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE349   .SAMPLE35  .SAMPLE350  .SAMPLE351  .SAMPLE352  .SAMPLE353 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE354  .SAMPLE355  .SAMPLE356  .SAMPLE357  .SAMPLE358  .SAMPLE359 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE36  .SAMPLE360  .SAMPLE361  .SAMPLE362  .SAMPLE363  .SAMPLE364 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE365  .SAMPLE366  .SAMPLE367  .SAMPLE368  .SAMPLE369   .SAMPLE37 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE370  .SAMPLE371  .SAMPLE372  .SAMPLE373  .SAMPLE374  .SAMPLE375 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE376  .SAMPLE377  .SAMPLE378  .SAMPLE379   .SAMPLE38  .SAMPLE380 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE381  .SAMPLE382  .SAMPLE383  .SAMPLE384  .SAMPLE385  .SAMPLE386 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE387  .SAMPLE388  .SAMPLE389   .SAMPLE39  .SAMPLE390  .SAMPLE391 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE392  .SAMPLE393  .SAMPLE394  .SAMPLE395  .SAMPLE396  .SAMPLE397 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE398  .SAMPLE399    .SAMPLE4   .SAMPLE40  .SAMPLE400  .SAMPLE401 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE402  .SAMPLE403  .SAMPLE404  .SAMPLE405  .SAMPLE406  .SAMPLE407 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE408  .SAMPLE409   .SAMPLE41  .SAMPLE410  .SAMPLE411  .SAMPLE412 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE413  .SAMPLE414  .SAMPLE415  .SAMPLE416  .SAMPLE417  .SAMPLE418 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE419   .SAMPLE42  .SAMPLE420  .SAMPLE421  .SAMPLE422  .SAMPLE423 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE424  .SAMPLE425  .SAMPLE426  .SAMPLE427  .SAMPLE428  .SAMPLE429 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE43  .SAMPLE430  .SAMPLE431  .SAMPLE432  .SAMPLE433  .SAMPLE434 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE435  .SAMPLE436  .SAMPLE437  .SAMPLE438  .SAMPLE439   .SAMPLE44 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE440  .SAMPLE441  .SAMPLE442  .SAMPLE443  .SAMPLE444  .SAMPLE445 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE446  .SAMPLE447  .SAMPLE448  .SAMPLE449   .SAMPLE45  .SAMPLE450 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE451  .SAMPLE452  .SAMPLE453  .SAMPLE454  .SAMPLE455  .SAMPLE456 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE457  .SAMPLE458  .SAMPLE459   .SAMPLE46  .SAMPLE460  .SAMPLE461 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE462  .SAMPLE463  .SAMPLE464  .SAMPLE465  .SAMPLE466  .SAMPLE467 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE468  .SAMPLE469   .SAMPLE47  .SAMPLE470  .SAMPLE471  .SAMPLE472 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE473  .SAMPLE474  .SAMPLE475  .SAMPLE476  .SAMPLE477  .SAMPLE478 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE479   .SAMPLE48  .SAMPLE480  .SAMPLE481  .SAMPLE482  .SAMPLE483 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE484  .SAMPLE485  .SAMPLE486  .SAMPLE487  .SAMPLE488  .SAMPLE489 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE49  .SAMPLE490  .SAMPLE491  .SAMPLE492  .SAMPLE493  .SAMPLE494 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE495  .SAMPLE496  .SAMPLE497  .SAMPLE498  .SAMPLE499    .SAMPLE5 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE50  .SAMPLE500  .SAMPLE501  .SAMPLE502  .SAMPLE503  .SAMPLE504 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE505  .SAMPLE506  .SAMPLE507  .SAMPLE508  .SAMPLE509   .SAMPLE51 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE510  .SAMPLE511  .SAMPLE512  .SAMPLE513  .SAMPLE514  .SAMPLE515 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE516  .SAMPLE517  .SAMPLE518  .SAMPLE519   .SAMPLE52  .SAMPLE520 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE521  .SAMPLE522  .SAMPLE523  .SAMPLE524  .SAMPLE525  .SAMPLE526 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE527  .SAMPLE528  .SAMPLE529   .SAMPLE53  .SAMPLE530  .SAMPLE531 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE532  .SAMPLE533  .SAMPLE534  .SAMPLE535  .SAMPLE536  .SAMPLE537 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE538  .SAMPLE539   .SAMPLE54  .SAMPLE540  .SAMPLE541  .SAMPLE542 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE543  .SAMPLE544  .SAMPLE545  .SAMPLE546  .SAMPLE547  .SAMPLE548 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE549   .SAMPLE55  .SAMPLE550  .SAMPLE551  .SAMPLE552  .SAMPLE553 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE554  .SAMPLE555  .SAMPLE556  .SAMPLE557  .SAMPLE558  .SAMPLE559 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE56  .SAMPLE560  .SAMPLE561  .SAMPLE562  .SAMPLE563  .SAMPLE564 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE565  .SAMPLE566  .SAMPLE567  .SAMPLE568  .SAMPLE569   .SAMPLE57 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE570  .SAMPLE571  .SAMPLE572  .SAMPLE573  .SAMPLE574  .SAMPLE575 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE576  .SAMPLE577  .SAMPLE578  .SAMPLE579   .SAMPLE58  .SAMPLE580 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE581  .SAMPLE582  .SAMPLE583  .SAMPLE584  .SAMPLE585  .SAMPLE586 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE587  .SAMPLE588  .SAMPLE589   .SAMPLE59  .SAMPLE590  .SAMPLE591 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE592  .SAMPLE593  .SAMPLE594  .SAMPLE595  .SAMPLE596  .SAMPLE597 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE598  .SAMPLE599    .SAMPLE6   .SAMPLE60  .SAMPLE600  .SAMPLE601 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE602  .SAMPLE603  .SAMPLE604  .SAMPLE605  .SAMPLE606  .SAMPLE607 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE608  .SAMPLE609   .SAMPLE61  .SAMPLE610  .SAMPLE611  .SAMPLE612 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE613  .SAMPLE614  .SAMPLE615  .SAMPLE616  .SAMPLE617  .SAMPLE618 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE619   .SAMPLE62  .SAMPLE620  .SAMPLE621  .SAMPLE622  .SAMPLE623 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE624  .SAMPLE625  .SAMPLE626  .SAMPLE627  .SAMPLE628  .SAMPLE629 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE63  .SAMPLE630  .SAMPLE631  .SAMPLE632  .SAMPLE633  .SAMPLE634 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE635  .SAMPLE636  .SAMPLE637  .SAMPLE638  .SAMPLE639   .SAMPLE64 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE640  .SAMPLE641  .SAMPLE642  .SAMPLE643  .SAMPLE644  .SAMPLE645 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE646  .SAMPLE647  .SAMPLE648  .SAMPLE649   .SAMPLE65  .SAMPLE650 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE651  .SAMPLE652  .SAMPLE653  .SAMPLE654  .SAMPLE655  .SAMPLE656 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE657  .SAMPLE658  .SAMPLE659   .SAMPLE66  .SAMPLE660  .SAMPLE661 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE662  .SAMPLE663  .SAMPLE664  .SAMPLE665  .SAMPLE666  .SAMPLE667 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE668  .SAMPLE669   .SAMPLE67  .SAMPLE670  .SAMPLE671  .SAMPLE672 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE673  .SAMPLE674  .SAMPLE675  .SAMPLE676  .SAMPLE677  .SAMPLE678 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE679   .SAMPLE68  .SAMPLE680  .SAMPLE681  .SAMPLE682  .SAMPLE683 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE684  .SAMPLE685  .SAMPLE686  .SAMPLE687  .SAMPLE688  .SAMPLE689 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE69  .SAMPLE690  .SAMPLE691  .SAMPLE692  .SAMPLE693  .SAMPLE694 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE695  .SAMPLE696  .SAMPLE697  .SAMPLE698  .SAMPLE699    .SAMPLE7 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE70  .SAMPLE700  .SAMPLE701  .SAMPLE702  .SAMPLE703  .SAMPLE704 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE705  .SAMPLE706  .SAMPLE707  .SAMPLE708  .SAMPLE709   .SAMPLE71 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE710  .SAMPLE711  .SAMPLE712  .SAMPLE713  .SAMPLE714  .SAMPLE715 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE716  .SAMPLE717  .SAMPLE718  .SAMPLE719   .SAMPLE72  .SAMPLE720 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE721  .SAMPLE722  .SAMPLE723  .SAMPLE724  .SAMPLE725  .SAMPLE726 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE727  .SAMPLE728  .SAMPLE729   .SAMPLE73  .SAMPLE730  .SAMPLE731 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE732  .SAMPLE733  .SAMPLE734  .SAMPLE735  .SAMPLE736  .SAMPLE737 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE738  .SAMPLE739   .SAMPLE74  .SAMPLE740  .SAMPLE741  .SAMPLE742 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE743  .SAMPLE744  .SAMPLE745  .SAMPLE746  .SAMPLE747  .SAMPLE748 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE749   .SAMPLE75  .SAMPLE750  .SAMPLE751  .SAMPLE752  .SAMPLE753 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE754  .SAMPLE755  .SAMPLE756  .SAMPLE757  .SAMPLE758  .SAMPLE759 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE76  .SAMPLE760  .SAMPLE761  .SAMPLE762  .SAMPLE763  .SAMPLE764 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE765  .SAMPLE766  .SAMPLE767  .SAMPLE768  .SAMPLE769   .SAMPLE77 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE770  .SAMPLE771  .SAMPLE772  .SAMPLE773  .SAMPLE774  .SAMPLE775 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE776  .SAMPLE777  .SAMPLE778  .SAMPLE779   .SAMPLE78  .SAMPLE780 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE781  .SAMPLE782  .SAMPLE783  .SAMPLE784  .SAMPLE785  .SAMPLE786 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE787  .SAMPLE788  .SAMPLE789   .SAMPLE79  .SAMPLE790  .SAMPLE791 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE792  .SAMPLE793  .SAMPLE794  .SAMPLE795  .SAMPLE796  .SAMPLE797 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE798  .SAMPLE799    .SAMPLE8   .SAMPLE80  .SAMPLE800  .SAMPLE801 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE802  .SAMPLE803  .SAMPLE804  .SAMPLE805  .SAMPLE806  .SAMPLE807 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE808  .SAMPLE809   .SAMPLE81  .SAMPLE810  .SAMPLE811  .SAMPLE812 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE813  .SAMPLE814  .SAMPLE815  .SAMPLE816  .SAMPLE817  .SAMPLE818 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE819   .SAMPLE82  .SAMPLE820  .SAMPLE821  .SAMPLE822  .SAMPLE823 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE824  .SAMPLE825  .SAMPLE826  .SAMPLE827  .SAMPLE828  .SAMPLE829 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE83  .SAMPLE830  .SAMPLE831  .SAMPLE832  .SAMPLE833  .SAMPLE834 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE835  .SAMPLE836  .SAMPLE837  .SAMPLE838  .SAMPLE839   .SAMPLE84 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE840  .SAMPLE841  .SAMPLE842  .SAMPLE843  .SAMPLE844  .SAMPLE845 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE846  .SAMPLE847  .SAMPLE848  .SAMPLE849   .SAMPLE85  .SAMPLE850 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE851  .SAMPLE852  .SAMPLE853  .SAMPLE854  .SAMPLE855  .SAMPLE856 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE857  .SAMPLE858  .SAMPLE859   .SAMPLE86  .SAMPLE860  .SAMPLE861 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE862  .SAMPLE863  .SAMPLE864  .SAMPLE865  .SAMPLE866  .SAMPLE867 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE868  .SAMPLE869   .SAMPLE87  .SAMPLE870  .SAMPLE871  .SAMPLE872 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE873  .SAMPLE874  .SAMPLE875  .SAMPLE876  .SAMPLE877  .SAMPLE878 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE879   .SAMPLE88  .SAMPLE880  .SAMPLE881  .SAMPLE882  .SAMPLE883 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE884  .SAMPLE885  .SAMPLE886  .SAMPLE887  .SAMPLE888  .SAMPLE889 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE89  .SAMPLE890  .SAMPLE891  .SAMPLE892  .SAMPLE893  .SAMPLE894 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE895  .SAMPLE896  .SAMPLE897  .SAMPLE898  .SAMPLE899    .SAMPLE9 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE90  .SAMPLE900  .SAMPLE901  .SAMPLE902  .SAMPLE903  .SAMPLE904 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE905  .SAMPLE906  .SAMPLE907  .SAMPLE908  .SAMPLE909   .SAMPLE91 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE910  .SAMPLE911  .SAMPLE912  .SAMPLE913  .SAMPLE914  .SAMPLE915 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE916  .SAMPLE917  .SAMPLE918  .SAMPLE919   .SAMPLE92  .SAMPLE920 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE921  .SAMPLE922  .SAMPLE923  .SAMPLE924  .SAMPLE925  .SAMPLE926 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE927  .SAMPLE928  .SAMPLE929   .SAMPLE93  .SAMPLE930  .SAMPLE931 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE932  .SAMPLE933  .SAMPLE934  .SAMPLE935  .SAMPLE936  .SAMPLE937 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE938  .SAMPLE939   .SAMPLE94  .SAMPLE940  .SAMPLE941  .SAMPLE942 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE943  .SAMPLE944  .SAMPLE945  .SAMPLE946  .SAMPLE947  .SAMPLE948 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE949   .SAMPLE95  .SAMPLE950  .SAMPLE951  .SAMPLE952  .SAMPLE953 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE954  .SAMPLE955  .SAMPLE956  .SAMPLE957  .SAMPLE958  .SAMPLE959 
    ##           1           1           1           1           1           1 
    ##   .SAMPLE96  .SAMPLE960  .SAMPLE961  .SAMPLE962  .SAMPLE963  .SAMPLE964 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE965  .SAMPLE966  .SAMPLE967  .SAMPLE968  .SAMPLE969   .SAMPLE97 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE970  .SAMPLE971  .SAMPLE972  .SAMPLE973  .SAMPLE974  .SAMPLE975 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE976  .SAMPLE977  .SAMPLE978  .SAMPLE979   .SAMPLE98  .SAMPLE980 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE981  .SAMPLE982  .SAMPLE983  .SAMPLE984  .SAMPLE985  .SAMPLE986 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE987  .SAMPLE988  .SAMPLE989   .SAMPLE99  .SAMPLE990  .SAMPLE991 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE992  .SAMPLE993  .SAMPLE994  .SAMPLE995  .SAMPLE996  .SAMPLE997 
    ##           1           1           1           1           1           1 
    ##  .SAMPLE998  .SAMPLE999 
    ##           1           1

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

    ##   BioProject         Run sample_name Age BMI Sex Region       Platform Country
    ## 1 PRJEB86033 ERR14648080   .SAMPLE10  NA  NA   F V3-V4  Illumina MiSeq   China
    ## 2 PRJEB86033 ERR14648081  .SAMPLE100  NA  NA   F V3-V4  Illumina MiSeq   China
    ## 3 PRJEB86033 ERR14648082  .SAMPLE101  NA  NA   F V3-V4  Illumina MiSeq   China
    ## 4 PRJEB86033 ERR14648083  .SAMPLE102  NA  NA   F V3-V4  Illumina MiSeq   China
    ## 5 PRJEB86033 ERR14648085  .SAMPLE104  NA  NA   F V3-V4  Illumina MiSeq   China
    ## 6 PRJEB86033 ERR14648086  .SAMPLE105  NA  NA   F V3-V4  Illumina MiSeq   China
    ##               Site          public_time PCR_primers
    ## 1 organic material 2025-02-22T00:00:00Z   341F-805R
    ## 2 organic material 2025-02-22T00:00:00Z   341F-805R
    ## 3 organic material 2025-02-22T00:00:00Z   341F-805R
    ## 4 organic material 2025-02-22T00:00:00Z   341F-805R
    ## 5 organic material 2025-02-22T00:00:00Z   341F-805R
    ## 6 organic material 2025-02-22T00:00:00Z   341F-805R
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
