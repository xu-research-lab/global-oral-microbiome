R Notebook: Metadata Processing Workflow
================

# Step 1: Obtain Initial Metadata from Each Project

This step sets the project directory, loads the raw SRA metadata, and
optionally merges supplement tables.

``` r
# Load libraries
library(dplyr)
library(ggplot2)
library(cowplot)

# (1) Setting working directory
# EDIT HERE: change project_id for each study
project_id <- "PRJNA744870" 

base_dir <- "E:/oral/26/metadata-2604/Oral-Metadata-YL-03"
setwd(file.path(base_dir, project_id))

# (2) Loading raw SRA metadata file
sra_file <- "SraRunTable.csv"
initial_meta <- read.csv(sra_file, header = TRUE)

# (3) OPTIONAL: If extra metadata (e.g., Table.x or supplemental materials) are available from the article, you may need to merge the file with initial_meta manually. Below is example code for merging. Depending on the specific case of each study, you can customize the code as needed.
# supp_table <- read.csv("supplemental_table.csv", header = TRUE)
# initial_meta <- left_join(initial_meta, supp_table, by = merge_col_name)
```

# Step 2: Data Overview

First, you need to get an overview of the data, which includes two
steps: reading the article and reading the initial metadata. This step
is very IMPORTANT!

``` r
# (1) Read the article in the following order: abstract, the beginning of the results section, and the methods. Understand the sample size, experimental design, and disease group overview. Clarify the sampling sites; for oral and gut samples, specify the exact sites like supragingival plaque, tongue, etc.

# (2) IMPORTANT: MUST go through each column of initial_meta
head(initial_meta, 10)
```

    ##            Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample
    ## 1  SRR15719926   AMPLICON        602 23042754 PRJNA744870 SAMN21215907
    ## 2  SRR15719937   AMPLICON        602 23526160 PRJNA744870 SAMN21215897
    ## 3  SRR15719951   AMPLICON        602 28818944 PRJNA744870 SAMN21215884
    ## 4  SRR15719961   AMPLICON        602 22373932 PRJNA744870 SAMN21215816
    ## 5  SRR15719966   AMPLICON        602  4212796 PRJNA744870 SAMN21215815
    ## 6  SRR15719973   AMPLICON        602 36488424 PRJNA744870 SAMN21215951
    ## 7  SRR15719989   AMPLICON        602  5851440 PRJNA744870 SAMN21215830
    ## 8  SRR15720004   AMPLICON        602  4814194 PRJNA744870 SAMN21215826
    ## 9  SRR15720010   AMPLICON        602  7121058 PRJNA744870 SAMN21215825
    ## 10 SRR15720015   AMPLICON        602 21083244 PRJNA744870 SAMN21215837
    ##                 BioSampleModel    Bytes                         Center.Name
    ## 1  Metagenome or environmental 15416009 THE CHINESE UNIVERSITY OF HONG KONG
    ## 2  Metagenome or environmental 15389328 THE CHINESE UNIVERSITY OF HONG KONG
    ## 3  Metagenome or environmental 17930945 THE CHINESE UNIVERSITY OF HONG KONG
    ## 4  Metagenome or environmental 13334903 THE CHINESE UNIVERSITY OF HONG KONG
    ## 5  Metagenome or environmental  2130157 THE CHINESE UNIVERSITY OF HONG KONG
    ## 6  Metagenome or environmental 22725477 THE CHINESE UNIVERSITY OF HONG KONG
    ## 7  Metagenome or environmental  2882275 THE CHINESE UNIVERSITY OF HONG KONG
    ## 8  Metagenome or environmental  2357889 THE CHINESE UNIVERSITY OF HONG KONG
    ## 9  Metagenome or environmental  3421802 THE CHINESE UNIVERSITY OF HONG KONG
    ## 10 Metagenome or environmental 12511086 THE CHINESE UNIVERSITY OF HONG KONG
    ##    Collection_Date Consent DATASTORE.filetype DATASTORE.provider
    ## 1             2016  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 2             2016  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 3             2016  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 4             2016  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5             2016  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 6             2016  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 7             2016  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 8             2016  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 9             2016  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 10            2016  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                        DATASTORE.region  Experiment geo_loc_name_country
    ## 1  gs.us-east1,ncbi.public,s3.us-east-1 SRX12015501            Hong Kong
    ## 2  gs.us-east1,ncbi.public,s3.us-east-1 SRX12015490            Hong Kong
    ## 3  gs.us-east1,ncbi.public,s3.us-east-1 SRX12015476            Hong Kong
    ## 4  gs.us-east1,ncbi.public,s3.us-east-1 SRX12015466            Hong Kong
    ## 5  gs.us-east1,ncbi.public,s3.us-east-1 SRX12015461            Hong Kong
    ## 6  gs.us-east1,ncbi.public,s3.us-east-1 SRX12015834            Hong Kong
    ## 7  gs.us-east1,ncbi.public,s3.us-east-1 SRX12015818            Hong Kong
    ## 8  gs.us-east1,ncbi.public,s3.us-east-1 SRX12015803            Hong Kong
    ## 9  gs.us-east1,ncbi.public,s3.us-east-1 SRX12015797            Hong Kong
    ## 10 gs.us-east1,ncbi.public,s3.us-east-1 SRX12015792            Hong Kong
    ##    geo_loc_name_country_continent geo_loc_name         HOST     Instrument
    ## 1                            Asia    Hong Kong Homo sapiens Illumina MiSeq
    ## 2                            Asia    Hong Kong Homo sapiens Illumina MiSeq
    ## 3                            Asia    Hong Kong Homo sapiens Illumina MiSeq
    ## 4                            Asia    Hong Kong Homo sapiens Illumina MiSeq
    ## 5                            Asia    Hong Kong Homo sapiens Illumina MiSeq
    ## 6                            Asia    Hong Kong Homo sapiens Illumina MiSeq
    ## 7                            Asia    Hong Kong Homo sapiens Illumina MiSeq
    ## 8                            Asia    Hong Kong Homo sapiens Illumina MiSeq
    ## 9                            Asia    Hong Kong Homo sapiens Illumina MiSeq
    ## 10                           Asia    Hong Kong Homo sapiens Illumina MiSeq
    ##    lat_lon             Library.Name LibraryLayout LibrarySelection
    ## 1  missing    Health26-Oral-16SV3V4        PAIRED              PCR
    ## 2  missing    Health16-Oral-16SV3V4        PAIRED              PCR
    ## 3  missing    Health03-Oral-16SV3V4        PAIRED              PCR
    ## 4  missing Patient11-PostM3-16SV3V4        PAIRED              PCR
    ## 5  missing Patient10-PostM3-16SV3V4        PAIRED              PCR
    ## 6  missing    Health70-Oral-16SV3V4        PAIRED              PCR
    ## 7  missing Patient25-PostM1-16SV3V4        PAIRED              PCR
    ## 8  missing Patient21-PostM6-16SV3V4        PAIRED              PCR
    ## 9  missing Patient20-PostM3-16SV3V4        PAIRED              PCR
    ## 10 missing Patient32-PostM3-16SV3V4        PAIRED              PCR
    ##    LibrarySource        Organism Platform          ReleaseDate
    ## 1        GENOMIC oral metagenome ILLUMINA 2021-09-28T00:00:00Z
    ## 2        GENOMIC oral metagenome ILLUMINA 2021-09-28T00:00:00Z
    ## 3        GENOMIC oral metagenome ILLUMINA 2021-09-28T00:00:00Z
    ## 4        GENOMIC oral metagenome ILLUMINA 2021-09-28T00:00:00Z
    ## 5        GENOMIC oral metagenome ILLUMINA 2021-09-28T00:00:00Z
    ## 6        GENOMIC oral metagenome ILLUMINA 2021-09-28T00:00:00Z
    ## 7        GENOMIC oral metagenome ILLUMINA 2021-09-28T00:00:00Z
    ## 8        GENOMIC oral metagenome ILLUMINA 2021-09-28T00:00:00Z
    ## 9        GENOMIC oral metagenome ILLUMINA 2021-09-28T00:00:00Z
    ## 10       GENOMIC oral metagenome ILLUMINA 2021-09-28T00:00:00Z
    ##             create_date version Sample.Name SRA.Study   isolate sample_type
    ## 1  2021-09-04T12:36:00Z       1    Health26 SRP328065  Health26  oral rinse
    ## 2  2021-09-04T12:36:00Z       1    Health16 SRP328065  Health16  oral rinse
    ## 3  2021-09-04T12:37:00Z       1    Health03 SRP328065  Health03  oral rinse
    ## 4  2021-09-04T12:36:00Z       1   Patient11 SRP328065 Patient11  oral rinse
    ## 5  2021-09-04T12:36:00Z       1   Patient10 SRP328065 Patient10  oral rinse
    ## 6  2021-09-04T12:36:00Z       1    Health70 SRP328065  Health70  oral rinse
    ## 7  2021-09-04T12:36:00Z       1   Patient25 SRP328065 Patient25  oral rinse
    ## 8  2021-09-04T12:36:00Z       1   Patient21 SRP328065 Patient21  oral rinse
    ## 9  2021-09-04T12:36:00Z       1   Patient20 SRP328065 Patient20  oral rinse
    ## 10 2021-09-04T12:37:00Z       1   Patient32 SRP328065 Patient32  oral rinse
    ##      strain isolation_source
    ## 1  Bacteria               NA
    ## 2  Bacteria               NA
    ## 3  Bacteria               NA
    ## 4  Bacteria               NA
    ## 5  Bacteria               NA
    ## 6  Bacteria               NA
    ## 7  Bacteria               NA
    ## 8  Bacteria               NA
    ## 9  Bacteria               NA
    ## 10 Bacteria               NA

``` r
# IMPORTANT: Check the unique values in the column of interest or query column.
lapply(initial_meta[c('sample_type')], table) # site
```

    ## $sample_type
    ## 
    ## oral rinse 
    ##        384

``` r
lapply(initial_meta[c('Organism')], table) # organism
```

    ## $Organism
    ## 
    ## human oral metagenome       oral metagenome 
    ##                     4                   380

``` r
table(initial_meta$Sample.Name) # subject id
```

    ## 
    ##   Health01   Health02   Health03   Health04   Health05   Health06   Health07 
    ##          1          1          1          1          1          1          1 
    ##   Health08   Health09   Health10   Health11   Health12   Health13   Health14 
    ##          1          1          1          1          1          1          1 
    ##   Health15   Health16   Health17   Health18   Health19   Health20   Health21 
    ##          1          1          1          1          1          1          1 
    ##   Health22   Health23   Health24   Health25   Health26   Health27   Health28 
    ##          1          1          1          1          1          1          1 
    ##   Health29   Health30   Health31   Health32   Health33   Health34   Health35 
    ##          1          1          1          1          1          1          1 
    ##   Health36   Health37   Health38   Health39   Health40   Health41   Health42 
    ##          1          1          1          1          1          1          1 
    ##   Health43   Health44   Health45   Health46   Health47   Health48   Health49 
    ##          1          1          1          1          1          1          1 
    ##   Health50   Health51   Health52   Health53   Health54   Health55   Health56 
    ##          1          1          1          1          1          1          1 
    ##   Health57   Health58   Health59   Health60   Health61   Health62   Health63 
    ##          1          1          1          1          1          1          1 
    ##   Health64   Health65   Health66   Health67   Health68   Health69   Health70 
    ##          1          1          1          1          1          1          1 
    ##   Health71   Health72   Health73   Health74   Health75   Health76  Patient01 
    ##          1          1          1          1          1          1          4 
    ## Patient01R  Patient02  Patient03  Patient04  Patient05  Patient06  Patient07 
    ##          4          4          4          4          4          4          4 
    ##  Patient08  Patient09  Patient10  Patient11  Patient12  Patient13  Patient14 
    ##          4          4          4          4          4          4          4 
    ##  Patient15  Patient16  Patient17  Patient18  Patient19  Patient20  Patient21 
    ##          4          4          4          4          4          4          4 
    ##  Patient22  Patient23  Patient24  Patient25  Patient26  Patient27  Patient28 
    ##          4          4          4          4          4          4          4 
    ##  Patient29  Patient30  Patient31  Patient32  Patient33  Patient34  Patient35 
    ##          4          4          4          4          4          4          4 
    ##  Patient36  Patient37  Patient38  Patient39  Patient40  Patient41  Patient42 
    ##          4          4          4          4          4          4          4 
    ##  Patient43  Patient44  Patient45  Patient46  Patient47  Patient48  Patient49 
    ##          4          4          4          4          4          4          4 
    ##  Patient50  Patient51  Patient52  Patient53  Patient54  Patient55  Patient56 
    ##          4          4          4          4          4          4          4 
    ##  Patient57  Patient58  Patient59  Patient60  Patient61  Patient62  Patient63 
    ##          4          4          4          4          4          4          4 
    ##  Patient64  Patient65  Patient66  Patient67  Patient68  Patient69  Patient70 
    ##          4          4          4          4          4          4          4 
    ##  Patient71  Patient72  Patient73  Patient74  Patient75  Patient76 
    ##          4          4          4          4          4          4

``` r
#table(initial_meta$sequencing_meth..exp.) # output is Sequencing by synthesis
table(initial_meta$Platform) # output is Illumina. according to the methods,  Illumina is correct.
```

    ## 
    ## ILLUMINA 
    ##      384

``` r
#table(initial_meta$target_subfragment) # output is V4
```

# Step 3: Removing Control and Non-human Samples

First, we need to remove control samples(including both negative
controls and positive controls) and non-human samples.

``` r
# (1) EDIT HERE: Parameter setup
control_column <- "" # column name in initial_meta
# unique(initial_meta[[control_column]])
control_keywords <- c() # Enter exact keywords

species_column <- "" # column name in initial_meta
unique(initial_meta[[species_column]])
```

    ## NULL

``` r
non_human_keywords <- c() # Enter exact keywords

# (2) Remove controls and non-human samples. DO NOT need to modify code in this part, just run it.
n_samples <- nrow(initial_meta)
# Control evaluation
if (length(control_keywords) > 0) {
  if (is.null(control_column) || control_column == "") {
    stop("Error: control_keywords provided but control_column is missing.")
  }
  
  if (!(control_column %in% colnames(initial_meta))) {
    stop(sprintf("Error: Column '%s' not found in initial_meta.", control_column))
  }
  
  control_pattern <- paste(control_keywords, collapse = "|")
  control_text <- initial_meta[[control_column]]
  is_control <- grepl(control_pattern, control_text)
  
} else {
  is_control <- rep(FALSE, n_samples)
}

# Non-human evaluation
if (length(non_human_keywords) > 0) {
  if (is.null(species_column) || species_column == "") {
    stop("Error: non_human_keywords provided but species_column is missing.")
  }
  
  if (!(species_column %in% colnames(initial_meta))) {
    stop(sprintf("Error: Column '%s' not found in initial_meta.", species_column))
  }
  
  non_human_pattern <- paste(non_human_keywords, collapse = "|")
  species_text <- initial_meta[[species_column]]
  is_non_human <- grepl(non_human_pattern, species_text)
} else {
  is_non_human <- rep(FALSE, n_samples)
}

# Execute removal
samples_to_remove <- is_control | is_non_human
initial_meta <- initial_meta[!samples_to_remove, ]

cat("Initial samples counts:",n_samples,"\n")
```

    ## Initial samples counts: 384

``` r
cat("Removed control samples:", sum(is_control), "\n")
```

    ## Removed control samples: 0

``` r
cat("Removed non-human samples:", sum(is_non_human), "\n")
```

    ## Removed non-human samples: 0

``` r
cat("Total samples removed ", sum(samples_to_remove), "\n")
```

    ## Total samples removed  0

``` r
cat("Remaining samples after Step 2:", nrow(initial_meta), "\n")
```

    ## Remaining samples after Step 2: 384

# Step 4: Creating a Uniform Metadata Template

This step converts the raw SRA table into the standard metadata format
used in our meta-analysis. The uniform metadata should include the
following columns: 1. **bioproject**: Project ID from SRA or other
database e.g. `PRJNA1234`, `CNP1234`, ‘OEP1234’, ‘qiita_1234’ 2.
**sample_id**: SRA sample/run ID. Unique id for each sample!
e.g. `DRR1234567`, `ERR1234567` and `SRR1234567`. If no such column
exists, can use column like sample_name. 3. **subject_id**: subject id
for each person if available. e.g. `patient_01`, `Host_03_before`,
`subject12` 4. **age**: Age of the subject if available. e.g. `35`, `42`
\| not `adult`, `chilren` 5. **age_info**: If the age is not explicitly
stated but approximate age range information is available, mark it in
this column. e.g., `18-55`, `children`, `infant`, `elder` 6. **bmi**:
Body Mass Index if available. e.g. `22.5`, `18.9` \| not `normal`,
`overweight` 7. **sex**: Sex of the subject. Mandatory to write as M/F
\| not `female`, `woman` 8. **sequencing_region**: Targeted 16S rRNA
region. e.g. `V3-V4`, `V4` 9. **platform**: Sequencing platform.
e.g. `Illumina MiSeq`, `454 GS FLX Titanium`, `Ion Torrent` 10.
**country**: Country of sample origin. e.g. `China`, `Japan` \| not city
name like `Beijing` 11. **site**: Sampling site. For oral and gut sites,
specify the exact sites as detailed as possible. For example, if it is
known to be supragingival plaque, do not record it as ‘dental plaque’,
unless the original text provides only a vague description.
e.g. `saliva`, `tongue coating` \| not broad term like `oral` 12.
**site_group**: only accepts one of following three values: ‘oral’,
‘gut’, or ‘other’. 13. **public_time**: Data release date.
e.g. `2021-06-22` \| not `June 2021` 14. **pcr_primers**: PCR primers
used. e.g. `515F-806R`, `27F-338R` \| not primer sequences 15.
**title**: Title of the associated publication. 16. **disease_study**:
Whether the study involves disease samples. If it’s not clear, can
record as NA. Only three values accepted: ‘yes’, ‘no’ or NA 17.
**status**: Disease status for each sample if available, else record as
NA. Only accept following abbreviations, Write full names for all other
diseases. HC: healthy control CD: Crohn’s disease CRC: colorectal cancer
COPD: chronic obstructive pulmonary disease HIV: human immunodeficiency
virus ICU: intensive care unit OSCC: oral squamous cell carcinoma RA:
rheumatoid arthritis UC: ulcerative colitis 18. **matched_samples**:
Whether the subject has matched samples. Only accept ‘yes’ or ‘no’. 19.
**is_time_series**: Whether the study is a time series study. Only
accept ‘yes’ or ‘no’.

``` r
# Extract columns from initial_meta as follows: (1) No modification needed → extract directly. (2) Modification needed → create a new column with the prefix raw_ and mark the original column as “pending”. (3) Need to find or determine based on the article.
metadata <- data.frame(
  bioproject = initial_meta$BioProject,
  sample_id = initial_meta$Run,
  subject_id = initial_meta$Sample.Name,
  age = NA, 
  age_info = "above 18 years old", # extract from methods section 
  bmi = NA,
  sex = NA,
  sequencing_region = "V3-V4", # extract from methods section
  platform = "Illumina MiSeq", # extract from methods section, NOT initial_meta
  country = "China",
  site = initial_meta$sample_type, 
  site_group = "pending", # modify later
  public_time = initial_meta$ReleaseDate,
  pcr_primers = "341F-806R", # Format: xxF-xxR
  title = "Restoration of the oral microbiota after surgery for head and neck squamous cell carcinoma is associated with patient outcomes", # extract from article
  disease_study = "yes", # determine according to the article
  raw_status = initial_meta$Sample.Name,
  status = "pending", # modify later
  matched_samples = "pending", # modify later
  is_time_series = "no", 
  stringsAsFactors = FALSE
)
metadata$public_time <- as.Date(metadata$public_time)
```

# Step 5: Manually solving pending issues

This step involves two parts: (1) handle the columns with the `raw_`
prefix, and next (2) handle the `site_group` and `matched_samples`
columns.

``` r
# (1) Handle the columns with the `raw_` prefix.
# In this part, you need to customize the code based on the situation of the raw_ columns to obtain the content of the target columns. IMPORTANT: Add comments if the modification is not straightforward!!
metadata <- metadata %>%
  mutate(status = case_when(
    grepl("Patient", raw_status) ~ "head and neck squamous cell carcinoma",
    grepl("Health", raw_status) ~ "HC",
    TRUE ~ NA_character_
  )) # according to the 'Sample.Name' column of initial_meta
```

``` r
# (2) Filling the 'site_group' and 'matched_samples' columns
# editing site_group
unique(metadata$site)
```

    ## [1] "oral rinse"

``` r
# EDIT HERE: add oral sites
Oral_site = unique(metadata$site)
# EDIT HERE: add gut sites
Gut_site = c()
# just run
metadata <- metadata %>%
  mutate(site_group = case_when(
    site %in% Oral_site ~ "oral",
    site %in% Gut_site ~ "gut",
    TRUE ~ "other" 
  ))

# editing matched_samples, just run.
metadata$matched_samples <- ifelse(
  ave(metadata$subject_id, metadata$subject_id, FUN = length) > 1, 
  "yes", 
  "no"
)
```

``` r
# (3) check if there is pending value in metadata. If there is pending value, go on modify the metadata.
# Function to check if a column contains "pending"
has_pending <- function(col) {
  if (is.character(col) || is.factor(col)) {
    any(col == "pending", na.rm = TRUE)
  } else {
    FALSE
  }
}

# Check for any pending value
any_pending_value <- any(sapply(metadata, has_pending))

if (!any_pending_value) {
  message("Good! No cell contains the value 'pending'.")
} else {
  pending_cols <- names(metadata)[sapply(metadata, has_pending)]
  message("Columns containing 'pending': ", paste(pending_cols, collapse = ", "))
  message("Please continue to resolve these pending columns!")
}
```

    ## Good! No cell contains the value 'pending'.

# Step 6: Visualize the Metadata

This step makes a quick set of plots for checking the distribution of
key metadata fields.

``` r
# Select columns for visualization, just run
plot_data <- metadata %>%
  select(-bioproject, -sample_id, -subject_id, -title)

create_simple_plots <- function(data) {
  plots <- list()
  continuous_vars <- c("age", "bmi")

  for (col in colnames(data)) {
    if (all(is.na(data[[col]]))) {
      p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "this column is NA", size = 6) +
        labs(title = col) +
        theme_void()
    } else if (col %in% continuous_vars && is.numeric(data[[col]])) {
      p <- ggplot(data, aes(x = .data[[col]])) +
        geom_histogram(fill = "blue", alpha = 0.6, bins = 20, na.rm = TRUE) +
        labs(title = col, x = col, y = "Count") +
        theme_minimal()
    } else {
      count_data <- as.data.frame(table(data[[col]], useNA = "ifany"))
      p <- ggplot(count_data, aes(x = reorder(Var1, -Freq), y = Freq)) +
        geom_col(fill = "red", alpha = 0.6) +
        labs(title = col, x = NULL, y = "Count") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 15, hjust = 1))
    }
    plots[[col]] <- p
  }

  return(plots)
}

# Generate and arrange plots
plots <- create_simple_plots(plot_data)
p <- plot_grid(plots$age,
               plots$age_info,
               plots$bmi,
               plots$sex,
               plots$site,
               plots$site_group,
               plots$sequencing_region,
               plots$platform,
               plots$country,
               plots$pcr_primers,
               plots$status,
               plots$matched_samples,
               ncol = 3, align = "hv", axis = "lb")
ggsave("Metadata_Summary_v2.png", plot = p, width = 8, height = 6, dpi = 150)
```

``` r
# Insert Metadata_Summary.png into the current Rmd output
knitr::include_graphics("Metadata_Summary_v2.png")
```

<img src="Metadata_Summary_v2.png" alt="" width="1200" />

# Step 7: Sort and export the metadata

This step saves the cleaned metadata table as the final `metadata.csv`
file.

``` r
# Reorder the data frame
fixed_cols <- c("bioproject", "sample_id", "subject_id", "age", "age_info",
                "bmi", "sex", "platform", "sequencing_region", "pcr_primers", 
                "country", "site", "site_group", "disease_study", "status", 
                "title", "public_time", "matched_samples", "is_time_series")  
metadata <- metadata %>%
  select(all_of(fixed_cols), everything())

# save metadata 
# head(metadata)
write.csv(metadata, file = "metadata_v2.csv", quote = TRUE)
```
