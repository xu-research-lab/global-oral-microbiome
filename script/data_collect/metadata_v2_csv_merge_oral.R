library(tidyverse)
library(fs)
library(countrycode)
library(tools) 

# ------------------------------------------------------------------------------
# 1. Path Setup & Initialization
# ------------------------------------------------------------------------------
base_dir <- "D:/YL2025/2026/global-oral-github-YL-v2/global-oral-microbiome/data/projects"
output_dir <- "E:/oral/26/metadata-2604/Oral-Metadata-merge-YL/results-260515"
output_raw <- file.path(output_dir, "metadata-v2-merge-oral.csv")

# Create output directory if it doesn't exist
dir_create(output_dir)

required_cols <- c(
  "bioproject", "sample_id", "subject_id", "age", "age_info", 
  "bmi", "sex", "platform", "sequencing_region", "pcr_primers", 
  "country", "site", "site_group", "disease_study", "status", 
  "title", "public_time", "matched_samples", "is_time_series"
)

# ------------------------------------------------------------------------------
# 2. Scan, Read, and Merge Data
# ------------------------------------------------------------------------------
message("Scanning and merging metadata files...")
subdirs <- dir_ls(base_dir, type = "directory")
expected_files <- file.path(subdirs, "metadata_v2.csv")
valid_files <- expected_files[file.exists(expected_files)]

df_list <- list()
for (file_path in valid_files) {
  folder_name <- basename(dirname(file_path))
  temp_df <- read_csv(file_path, col_types = cols(.default = "c"), show_col_types = FALSE)
  df_list[[folder_name]] <- temp_df
}

# Diagnostic: Identify folders containing the unnamed "...1" column
culprit_folders <- names(df_list) %>% keep(~ "...1" %in% colnames(df_list[[.x]]))
cat("\n[Diagnostic] Folders containing '...1' column:\n")
if (length(culprit_folders) > 0) print(culprit_folders) else cat("None found.\n")

# Merge dataframes and automatically drop unnamed columns
merged_df <- bind_rows(df_list) %>% select(-starts_with("..."))

# Ensure all 19 required columns exist (fill missing with NA)
for (col in required_cols) {
  if (!col %in% colnames(merged_df)) {
    merged_df[[col]] <- NA_character_
  }
}

# Standardize base layout: Required columns first, then others
other_cols <- setdiff(colnames(merged_df), required_cols)
df <- merged_df %>% select(all_of(required_cols), all_of(other_cols))

# ------------------------------------------------------------------------------
# 3. Data Filtering
# ------------------------------------------------------------------------------
message("Filtering valid site groups...")
df_filtered <- df %>%
  filter((site_group == "oral" | bioproject == "PRJNA659632") & 
           !(bioproject == "PRJNA471556" & site == "stool" & site_group == "oral"))

# ------------------------------------------------------------------------------
# 4. Core Cleaning & Generating New Columns
# ------------------------------------------------------------------------------
message("Cleaning data and mapping values...")

df_new <- df_filtered %>%
  mutate(
    # Basic age cleaning
    age_cleaned = case_when(
      str_detect(age, "(?i)week") ~ as.character(as.numeric(str_replace_all(age, "[^0-9.]", "")) / 52),
      TRUE ~ age
    ),
    age_numeric = as.numeric(age_cleaned),
    
    # Base logic to evaluate preliminary age_cleaned_group
    age_cleaned_group_base = case_when(
      !is.na(age_numeric) & (age_numeric < 0 | age_numeric > 150) ~ "Anomaly",
      !is.na(age_numeric) & age_numeric <= 0.5 ~ "0-6 months",
      !is.na(age_numeric) & age_numeric > 0.5 & age_numeric <= 3 ~ "6-36 months",
      !is.na(age_numeric) & age_numeric > 3 & age_numeric <= 12 ~ "3-12 years",
      !is.na(age_numeric) & age_numeric > 12 & age_numeric <= 18 ~ "12-18 years",
      !is.na(age_numeric) & age_numeric > 18 & age_numeric <= 60 ~ "18-60 years",
      !is.na(age_numeric) & age_numeric > 60 ~ "＞60 years",
      TRUE ~ NA_character_
    ),

    # Temporary lowercase standardization for age info mapping
    age_info_lower = tolower(trimws(age_info)),

    # Evaluate age_info_cleaned mapping rules
    age_info_cleaned = case_when(
      is.na(age_info_lower) ~ NA_character_,

      # Explicit NA matches
      age_info_lower %in% c(">=21", "≥21", "1-3 years at enrollment; longitudinal follow-up saliva samples in sra span approximately 1.34-5.57 years", "50+", "8 months-25 years", "age_group_1", "baby", "children", "child", "chilren", "elderly adults", "elder", "full cohort median age 41.5 years; healthy 38 (20-66), gingivitis 43 (18-71), chronic periodontitis 43 (21-60)", "infant", "not provided", "over 18 years", "toddlers 1-3 years; children 4-12 years; teenager-adults 13-59 years", "unspecified") ~ NA_character_,

      # Exact String Matches -> 18-60 years
      age_info_lower %in% c(
        "adult male betel-quid users", 
        "adult male non-smokers; cohort mean age 44 years (range 31-60)", 
        "adult participants; eligibility age >=18 years; total cohort median age 40.0 years (iqr 31.0-50.0)", 
        "adults; group means reported separately by disease severity", 
        "control 30.50±5.74; osa 40.11±11.25; htn 36.22±10.20; osa+htn 42.17±9.52; male adults >18 years", 
        "control non-cp: 34.7 ± 15.9 years; control cp: 46.3 ± 13.2 years; sle non-cp: 31.2 ± 5.2 years; sle cp: 40.5 ± 10.1 years", 
        "healthy dental students; mean age 26.36 +/- 1.58 years"
      ) ~ "18-60 years",

      # Range Mapping: (A, B)
      str_detect(age_info_lower, "^\\(\\s*\\d+\\.?\\d*\\s*,\\s*\\d+\\.?\\d*\\s*\\)$") ~ case_when(
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 0.5 ~ "0-6 months",
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 3 ~ "6-36 months",
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 12 ~ "3-12 years",
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 18 ~ "12-18 years",
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 60 ~ "18-60 years",
        TRUE ~ "＞60 years"
      ),

      # Month-based mapping
      str_detect(age_info_lower, "month") & str_detect(age_info_lower, "\\d+") ~ case_when(
        as.numeric(str_extract(age_info_lower, "\\d+")) <= 6 ~ "0-6 months",
        as.numeric(str_extract(age_info_lower, "\\d+")) <= 36 ~ "6-36 months",
        TRUE ~ NA_character_
      ),

      # Range Mapping: A-B (Calculates Mean)
      str_detect(age_info_lower, "^\\s*\\d+\\.?\\d*\\s*-\\s*\\d+\\.?\\d*\\s*$") ~ case_when(
        (as.numeric(str_extract(age_info_lower, "^\\s*\\d+\\.?\\d*")) + as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*\\s*$"))) / 2 <= 0.5 ~ "0-6 months",
        (as.numeric(str_extract(age_info_lower, "^\\s*\\d+\\.?\\d*")) + as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*\\s*$"))) / 2 <= 3 ~ "6-36 months",
        (as.numeric(str_extract(age_info_lower, "^\\s*\\d+\\.?\\d*")) + as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*\\s*$"))) / 2 <= 12 ~ "3-12 years",
        (as.numeric(str_extract(age_info_lower, "^\\s*\\d+\\.?\\d*")) + as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*\\s*$"))) / 2 <= 18 ~ "12-18 years",
        (as.numeric(str_extract(age_info_lower, "^\\s*\\d+\\.?\\d*")) + as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*\\s*$"))) / 2 <= 60 ~ "18-60 years",
        TRUE ~ "＞60 years"
      ),

      # Average/Mean mapping
      str_detect(age_info_lower, "average|mean") & !str_detect(age_info_lower, "±|\\+-|;|migraine") & str_detect(age_info_lower, "\\d+") ~ case_when(
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 0.5 ~ "0-6 months",
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 3 ~ "6-36 months",
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 12 ~ "3-12 years",
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 18 ~ "12-18 years",
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 60 ~ "18-60 years",
        TRUE ~ "＞60 years"
      ),

      # Decade strings (20s, 30s) mapping
      str_detect(age_info_lower, "\\b[1-9]0s\\b") ~ case_when(
        as.numeric(str_extract(age_info_lower, "\\d{2}")) < 12 ~ "3-12 years",
        as.numeric(str_extract(age_info_lower, "\\d{2}")) < 18 ~ "12-18 years",
        as.numeric(str_extract(age_info_lower, "\\d{2}")) < 60 ~ "18-60 years",
        TRUE ~ "＞60 years"
      ),

      # Variance: A ± B mapping
      str_detect(age_info_lower, "±") & str_count(age_info_lower, "±") == 1 & str_detect(age_info_lower, "\\d+") ~ case_when(
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 0.5 ~ "0-6 months",
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 3 ~ "6-36 months",
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 12 ~ "3-12 years",
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 18 ~ "12-18 years",
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 60 ~ "18-60 years",
        TRUE ~ "＞60 years"
      ),

      # Complex ';' statements
      str_detect(age_info_lower, ";") & str_detect(age_info_lower, "25-58 years; migraine mean 39\\.5") ~ "18-60 years",

      # 60+ mappings
      str_detect(age_info_lower, "60\\+|above 60") ~ "＞60 years",

      # Text-based adult/senior matches
      age_info_lower %in% c("adult", "adults") ~ "18-60 years",
      str_detect(age_info_lower, "centenarians") ~ "＞60 years",
      str_count(age_info_lower, "±") >= 2 ~ "18-60 years",

      # Child matching
      str_detect(age_info_lower, "child \\(\\d+-year-olds\\)") ~ case_when(
        as.numeric(str_extract(age_info_lower, "\\d+")) <= 0.5 ~ "0-6 months",
        as.numeric(str_extract(age_info_lower, "\\d+")) <= 3 ~ "6-36 months",
        as.numeric(str_extract(age_info_lower, "\\d+")) <= 12 ~ "3-12 years",
        as.numeric(str_extract(age_info_lower, "\\d+")) <= 18 ~ "12-18 years",
        TRUE ~ NA_character_
      ),

      # Pre-adult mapping
      str_detect(age_info_lower, "newborn|preterm infant") ~ "0-6 months",
      str_detect(age_info_lower, "older adults.*±") ~ "＞60 years",
      age_info_lower %in% c("teen", "teenager") ~ "12-18 years",
      age_info_lower == "young adults" ~ "18-60 years",

      # Variance: mean A +- B
      str_detect(age_info_lower, "mean.*\\+-") & str_detect(age_info_lower, "\\d+") ~ case_when(
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 18 ~ "12-18 years",
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 60 ~ "18-60 years",
        TRUE ~ "＞60 years"
      ),

      TRUE ~ NA_character_
    ),

    # Final combined logic for age_cleaned_group
    age_cleaned_group = case_when(
      !is.na(age_cleaned_group_base) & age_cleaned_group_base != "Anomaly" ~ age_cleaned_group_base,
      !is.na(age_info_cleaned) ~ age_info_cleaned,
      (!is.na(raw_host) & tolower(trimws(raw_host)) == "pregnant women") ~ "18-60 years",
      !is.na(age_cleaned_group_base) & age_cleaned_group_base == "Anomaly" ~ "Anomaly",
      TRUE ~ NA_character_
    ),

    # BMI Standardization
    bmi_calc = case_when(
      bioproject == "PRJEB34272" & bmi == "not collected" ~ NA_character_,
      TRUE ~ bmi
    ),
    bmi_raw_numeric = as.numeric(bmi_calc),
    bmi_cleaned = case_when(
      bmi_raw_numeric > 50 ~ NA_character_,
      TRUE ~ bmi_calc
    ),
    bmi_numeric = as.numeric(bmi_cleaned),
    bmi_cleaned_group = case_when(
      is.na(bmi_numeric) & bmi_raw_numeric > 50 ~ "Anomaly_High",
      is.na(bmi_numeric) ~ NA_character_,
      bmi_numeric < 0 ~ "Anomaly_High",
      bmi_numeric < 18.5 ~ "＜18.5",
      bmi_numeric >= 18.5 & bmi_numeric < 24.0 ~ "18.5-24.0",
      bmi_numeric >= 24.0 & bmi_numeric < 28.0 ~ "24.0-28.0",
      bmi_numeric >= 28.0 ~ "≥28.0",
      TRUE ~ NA_character_
    ),
    
    # Metadata string standardization
    sex_cleaned = case_when(
      bioproject %in% c("PRJDB12092", "PRJNA720269", "PRJDB10426") & sex == "female" ~ "F",
      bioproject %in% c("PRJDB12092", "PRJNA720269", "PRJDB10426") & sex == "male" ~ "M",
      TRUE ~ sex
    ),
    platform_cleaned = case_when(
      is.na(platform) ~ NA_character_,
      bioproject == "PRJNA659632" & platform == "ION TORRENT" ~ "Ion Torrent",
      bioproject == "PRJNA745286" & platform == "Illumina iSeq 100" ~ "Illumina MiSeq",
      bioproject == "PRJNA362794" & platform == "Illumina PE250" ~ "Illumina MiSeq",
      platform == "Illumina Miseq" ~ "Illumina MiSeq",
      str_detect(platform, "(?i)NextSeq") ~ "Illumina NextSeq",
      str_detect(platform, "(?i)454") ~ "454 GS FLX Titanium",
      str_detect(platform, "(?i)HiSeq") ~ "Illumina HiSeq",
      str_detect(platform, "(?i)NovaSeq") ~ "Illumina NovaSeq",
      str_detect(platform, "(?i)Ion Torrent") ~ "Ion Torrent",
      str_detect(platform, "(?i)GAIIx|Genome Analyzer IIx") ~ "Illumina GAIIx",
      TRUE ~ platform
    ),
    sequencing_region_cleaned = case_when(
      is.na(sequencing_region) ~ NA_character_,
      bioproject == "PRJNA659632" & sequencing_region == "V2, V3, V4, V6-7, V8, V9" ~ "V2,V3,V4,V6-V7,V8,V9",
      bioproject == "PRJNA391982" & sequencing_region == "V1-V3/V7-V9" ~ "V1-V3,V7-V9",
      str_detect(sequencing_region, "[–—_]") ~ str_replace_all(sequencing_region, "[–—_]+", "-"),
      TRUE ~ sequencing_region
    ),
    country_cleaned = case_when(
      bioproject == "PRJEB12953" & country == "Korean" ~ "Republic of Korea",
      bioproject == "PRJNA489752" & country == "Korea" ~ "Republic of Korea",
      bioproject == "PRJEB11419" & country == "Korea, Republic of" ~ "Republic of Korea",
      bioproject == "PRJEB11419" & country %in% c("Not provided", "Unspecified") ~ NA_character_,
      country %in% c("South Korea", "South Korean") ~ "Republic of Korea",
      country %in% c("United States of America", "USA", "United States Minor Outlying Islands") ~ "United States",
      TRUE ~ country
    ),
    continent_temp = countrycode(country_cleaned, origin = "country.name", destination = "continent", warn = FALSE),
    continent_group = case_when(
      bioproject %in% c("PRJEB12953", "PRJNA489752") ~ "Asia",
      is.na(country_cleaned) ~ NA_character_,
      continent_temp == "Africa" ~ "Africa",
      continent_temp == "Asia" ~ "Asia",
      continent_temp == "Europe" ~ "Europe",
      continent_temp == "Oceania" ~ "Oceania",
      continent_temp == "Americas" & countrycode(country_cleaned, origin="country.name", destination="un.regionsub.name", warn=FALSE) %in% c("Northern America") ~ "North America",
      continent_temp == "Americas" ~ "South America",
      str_detect(tolower(country_cleaned), "antarctica") ~ "Antarctica",
      TRUE ~ NA_character_ 
    ),
    site_lower = tolower(site),
    site_cleaned = case_when(
      is.na(site) ~ NA_character_,
      bioproject == "PRJEB4953" & site %in% c("Right retromolar trigone", "Right post maxilla ging", "Right floor of mouth", "Right alveolar ridge", "Midline hard palate", "Left retromolar trigone", "Left floor of mouth", "Left alveolar ridge", "floor of Mouth", "floor of mouth") ~ "oral mucosa",
      bioproject == "PRJNA542898" & site == "buccal surfaces" ~ "oral mucosa",
      bioproject == "PRJNA412445" & site == "mucous" ~ "oral mucosa",
      bioproject == "PRJNA820027" & site == "oral_cavity" ~ "other", 
      bioproject == "PRJNA362794" & site == "paracancerous tissues" ~ "other",
      bioproject == "PRJNA595346" & site == "posterior oropharyngeal swab" ~ "oral swab",
      bioproject == "PRJNA48333" & site == "throat" ~ "oral swab",
      bioproject == "PRJNA841771" & site == "oral scraping" ~ "other",
      site == "Left post mandible ging" ~ "gingiva",
      str_detect(site_lower, "hard palate") ~ "oral mucosa", 
      str_detect(site_lower, "oropharyngeal swab|throat") ~ "oral swab",
      site_lower %in% tolower(c("paracancerous tissues", "paracancerous tissue", "tumor tissue", "tumor", "oral cavity")) ~ "other",
      site %in% c("Right post maxilla ging", "implant sulcus") ~ "gingiva",
      site == "faucial pillars" ~ "oral mucosa",
      site %in% c("oral rinse", "oral washing", "oral lavage") ~ "saliva",
      site %in% c("cancer surface", "cancer tissue", "cancerous lesions", "oral", "normal tissue", "oral tissue (location unspecified)", "swab_nasopharynx_or_oropharynx", "cleft") ~ "other",
      str_detect(site_lower, "plaque|supragingival|subgingival|dental|carious dentin tissue|periodontal|carious|dentition") ~ "tooth",
      str_detect(site_lower, "oral wash|mouthwash|mouth wash") ~ "saliva",
      str_detect(site_lower, "mucosa|mucosal|buccal swab|buccal mucosa swabs") ~ "oral mucosa",
      str_detect(site_lower, "tongue") ~ "tongue",
      str_detect(site_lower, "gingiva|gingival") ~ "gingiva",
      str_detect(site_lower, "tonsil") ~ "Other",
      str_detect(site_lower, "oral swab") ~ "oral swab",
      str_detect(site_lower, "saliva") ~ "saliva",
      TRUE ~ site
    ),
    site_group_cleaned = case_when(
      bioproject == "PRJNA659632" & site_group == "other" ~ "oral",
      TRUE ~ site_group
    ),
    disease_study_cleaned = case_when(
      bioproject %in% c("PRJEB38610", "PRJEB39420", "PRJEB39436", "PRJNA417321", "PRJNA451534", "PRJNA876796") & tolower(disease_study) == "no" ~ NA_character_,
      bioproject %in% c("PRJNA380250", "PRJNA685840", "PRJNA847970", "PRJNA762545") ~ "yes",
      bioproject %in% c("PRJNA518756", "PRJNA613363", "PRJNA782768") ~ "no",
      TRUE ~ disease_study
    ),
    status_lower = tolower(status),
    status_temp = case_when(
      is.na(status) ~ NA_character_,
      status_lower %in% tolower(c("mucositis healed", "pre-conditioning", "atopic non-asthmatic controls", "pre-cancer", "full-term", "allogeneic hematopoietic stem cell transplantation recipient", "immunological ART responder", "dentally compromised", "preterm newborns", "preterm", "pulmonary nodule", "thyroid nodules", "black dental stain", "thyroid nodules patient", "epithelial precursor lesion")) ~ NA_character_,
      str_detect(status_lower, "breastfed|full-term|orthodontic|third molar extraction|non-responder|high-risk for ra") ~ NA_character_,
      bioproject == "PRJNA417321" & status == "dental implant" ~ NA_character_,
      status_lower %in% tolower(c("halitosis", "dysphagia", "watery diarrhea NOS", "truncal ataxia NOS", "oral benign tumors", "cancer patient", "Carcinoma in situ", "Chronic mucosal inflammation", "acute pancreatitis", "aortic thrombosis", "tibial fracture", "hip fracture", "hypotension", "STEMI", "NSTEMI", "temporal lobe epilepsy", "Fanconi anemia with leukoplakia", "Fanconi anemia with OSCC", "systemic lupus erythematosus with chronic periodontitis", "systemic lupus erythematosus", "acute lymphoblastic leukemia", "helminth infection", "pre-eclampsia", "stomatitis", "early childhood caries and congenital heart disease", "proliferative verrucous leukoplakia", "Refractory Celiac disease", "IBD+Diabetes")) ~ "Other",
      status_lower == tolower("Lung Cancer") ~ "Lung Cancer",
      status_lower == tolower("Breast Cancer") ~ "Breast Cancer",
      status_lower == tolower("Pancreatic Cancer") ~ "Pancreatic Cancer",
      status_lower == tolower("Oral Cancer") ~ "Oral Cancer",
      status_lower == tolower("Digestive Tract Cancer") ~ "Digestive Tract Cancer",
      status_lower == tolower("Prostate Cancer") ~ "Prostate Cancer",
      status_lower == tolower("Esophageal Cancer") ~ "Esophageal Cancer",
      status_lower == tolower("Colon Cancer") ~ "Colon Cancer",
      status_lower == tolower("viremic ART-untreated") ~ "HIV",
      status_lower %in% tolower(c("functional dyspepsia and gastroesophageal reflux disease", "GERD", "Reflux Esophagitis", "reflux esophagitis")) ~ "Gastroesophageal Reflux Disease",
      status_lower == tolower("RAU") ~ "Aphthous",
      status_lower == tolower("chronic periodontitis associated with type 2 diabetes mellitus") ~ "Diabetes and Periodontitis",
      bioproject == "PRJNA316469" & disease_study == "no" & status %in% c("never smoker", "current smoker") ~ "HC",
      bioproject == "PRJEB40804" & status == "certified saturation divers with valid health certificates" ~ "HC",
      bioproject == "PRJNA784772" & status == "DTC" ~ "differentiated thyroid carcinoma",
      bioproject == "PRJNA656387" & status == "first-degree relatives" ~ "HC",
      bioproject == "PRJNA853196" & status == "ESCC" ~ "esophageal squamous cell carcinoma",
      bioproject == "PRJDB8324" & status == "healthy" ~ "HC",
      bioproject == "PRJNA873789" & status == "axSpA" ~ "axial spondyloarthritis",
      bioproject == "PRJEB11419" & status == "C.diff" ~ "Clostridium difficile infection",
      bioproject == "PRJEB11419" & status == "IBD+C.diff" ~ "IBD and Clostridium difficile infection",
      bioproject == "PRJNA362643" & status == "Sjögren’s syndrome" ~ "Sjögren's syndrome",
      bioproject == "PRJNA664107" & status == "chronic periodontitis associated with T2DM" ~ "Diabetes and Periodontitis",
      bioproject == "PRJNA663926" & status %in% c("hand, foot and mouth disease symptomatic", "hand, foot and mouth disease asymptomatic") ~ "hand, foot and mouth disease",
      bioproject == "PRJNA544600" & status == "obstructive sleep apnea syndrome" ~ "obstructive sleep apnea",
      status == "caries" ~ "dental caries",
      status == "recurrent aphthous ulcer" ~ "aphthous",
      status_lower == tolower("CD and UC") ~ "CD and UC",
      status_lower == tolower("no-OM") ~ "HC",
      str_detect(status_lower, "periodontitis|gingivitis|periodontal|peri-implantitis") ~ "periodontal disease",
      str_detect(status_lower, "caries") ~ "dental caries",
      str_detect(status_lower, "coronavirus disease of 2019|covid") ~ "COVID-19",
      str_detect(status_lower, "chagas") ~ "chagas disease",
      str_detect(status_lower, "hiv") ~ "HIV",
      str_detect(status_lower, "obese|obesity") ~ "obesity",
      str_detect(status_lower, "aphthous") ~ "aphthous",
      str_detect(status_lower, "type 2 diabetes mellitus|t2dm|type-2 diabetes|type 1 diabetes|diabetes") ~ "diabetes",
      str_detect(status_lower, "atherosclerosis") ~ "cardiovascular disease",
      str_detect(status_lower, "oral mucositis|severe-om|grade 1 om") ~ "oral mucositis",
      str_detect(status_lower, "olp") ~ "oral lichen planus",
      str_detect(status_lower, "helicobacter pylori") ~ "helicobacter pylori infection",
      str_detect(status_lower, "thyroid cancer") ~ "thyroid cancer",
      str_detect(status_lower, "cancer") ~ "cancer",
      str_detect(status_lower, "reflux esophagitis") ~ "Gastroesophageal Reflux Disease",
      str_detect(status_lower, "axspa") ~ "axial spondyloarthritis",
      str_detect(status_lower, "nasopharyngeal carcinoma") ~ "nasopharyngeal carcinoma",
      str_detect(status_lower, "oscc|oral squamous cell carcinoma") ~ "OSCC",
      str_detect(status_lower, "head and neck squamous cell carcinoma|hnscc|head and neck cancer|hnc") ~ "head and neck cancer",
      str_detect(status_lower, "cleft lip and palate|clp") ~ "cleft lip and palate",
      str_detect(status_lower, "ad|alzheimer") ~ "Alzheimer's disease",
      str_detect(status_lower, "parkinson") ~ "Parkinson's disease",
      str_detect(status_lower, "post-traumatic stress disorder|ptsd") ~ "PTSD",
      TRUE ~ status
    ),
    status_cleaned = ifelse(!is.na(status_temp) & tolower(status_temp) != "na", tools::toTitleCase(status_temp), NA_character_)
  ) %>%
  select(-age_cleaned_group_base, -age_info_lower)

# Uniform Title Case processing 
df_new <- df_new %>%
  mutate(across(c(sex_cleaned, platform_cleaned, sequencing_region_cleaned, 
                  country_cleaned, continent_group, site_cleaned, 
                  site_group_cleaned, status_cleaned),
                ~ ifelse(!is.na(.), tools::toTitleCase(as.character(.)), NA_character_)))

# Final mapping adjustments and condition to convert NA to HC for specific bioprojects
df_new <- df_new %>%
    mutate(status_cleaned = case_when(
        is.na(status) & bioproject %in% c("PRJEB1799", "PRJEB20889", "PRJEB20890", "PRJEB20894", "PRJEB27304", "PRJEB28088", "PRJEB28769", "PRJEB33119", "PRJEB42056", "PRJEB6518", "PRJNA255332", "PRJNA395436", "PRJNA427101", "PRJNA634162", "PRJNA689848", "PRJNA721325", "PRJNA747639", "PRJNA778070", "PRJNA828039", "PRJNA841067") ~ "HC",
        status_cleaned == "Hiv" ~ "HIV",
        status_cleaned == "Oscc" ~ "OSCC",
        status_cleaned == "Hc" ~ "HC",
        status_cleaned == "Cd And Uc" ~ "CD and UC",
        status_cleaned == "Covid-19" ~ "COVID-19",
        status_cleaned == "Gastroesophageal Reflux Disease" ~ "Gastroesophageal Reflux Disease",
        TRUE ~ status_cleaned
    ))

# Standardize NA and empty strings and resolve specific anomalies to NAs
df_new <- df_new %>%
  mutate(
    across(c(age_cleaned, bmi_cleaned, age_cleaned_group, bmi_cleaned_group, sex_cleaned, platform_cleaned, 
             sequencing_region_cleaned, country_cleaned, continent_group, 
             site_cleaned, site_group_cleaned, disease_study_cleaned, status_cleaned),
           ~ ifelse(trimws(.) %in% c("", "NA"), NA_character_, .)),
    age_cleaned_group = ifelse(age_cleaned_group == "Anomaly", NA_character_, age_cleaned_group),
    bmi_cleaned_group = ifelse(bmi_cleaned_group == "Anomaly_High", NA_character_, bmi_cleaned_group)
  )

# ------------------------------------------------------------------------------
# 5. Reorder Columns and Save Output
# ------------------------------------------------------------------------------
message("Reordering columns and saving output dataset...")

idx_time_series <- which(names(df_new) == "is_time_series")

# Insert newly generated columns exactly after the core fields (around idx_time_series bounds)
new_cols_order <- c("age_cleaned", "bmi_cleaned", "age_cleaned_group", "age_info_cleaned", "bmi_cleaned_group", 
                    "sex_cleaned", "platform_cleaned", "sequencing_region_cleaned", 
                    "country_cleaned", "continent_group", "site_cleaned", 
                    "site_group_cleaned", "disease_study_cleaned", "status_cleaned")

# Fetch preceding and trailing columns natively based on dynamic indices
cols_before <- names(df_new)[1:idx_time_series]
cols_before <- setdiff(cols_before, new_cols_order) 

cols_after  <- names(df_new)[(idx_time_series + 1):length(names(df_new))]
cols_after  <- setdiff(cols_after, new_cols_order)

final_cols  <- c(cols_before, new_cols_order, cols_after)

df_final_raw <- df_new %>% select(all_of(final_cols))
write_csv(df_final_raw, output_raw)

message(sprintf("\nPipeline execution complete! Output saved to: %s\n", output_raw))