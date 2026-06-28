library(tidyverse)
library(fs)
library(countrycode)
library(tools)

# ------------------------------------------------------------------------------
# 1. Path Setup & Initialization
# ------------------------------------------------------------------------------
input_file <- "/home/yanglun/YSD/global_oral/metadata_v2_merge/oral-v2/results/metadata-v2-merge.csv"
output_file <- "/home/yanglun/YSD/global_oral/metadata_v2_merge/oral-v2/results/metadata-v2-merge-cleaned.csv"

required_cols <- c(
  "bioproject", "sample_id", "subject_id", "age", "age_info",
  "bmi", "sex", "platform", "sequencing_region", "pcr_primers",
  "country", "site", "site_group", "disease_study", "status",
  "title", "public_time", "matched_samples", "is_time_series"
)

message("Reading merged metadata...")
df_raw <- read_csv(input_file, col_types = cols(.default = "c"), show_col_types = FALSE)

# Keep only the required 19 columns plus raw_host (used in one cleaning rule)
# raw_* columns will be dropped from the final output
keep_cols <- c(required_cols, "raw_host")
keep_cols <- keep_cols[keep_cols %in% names(df_raw)]
df <- df_raw %>% select(all_of(keep_cols))

# Ensure all 19 required columns exist (fill missing with NA)
for (col in required_cols) {
  if (!col %in% colnames(df)) {
    df[[col]] <- NA_character_
  }
}

# Standardize base layout: Required columns first, then others
other_cols <- setdiff(colnames(df), required_cols)
df <- df %>% select(all_of(required_cols), all_of(other_cols))

# ------------------------------------------------------------------------------
# 2. Core Cleaning & Generating New Columns
# ------------------------------------------------------------------------------
message("Cleaning data and mapping values...")

df_new <- df %>%
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
      age_info_lower %in% c(
        "1-3 years at enrollment; longitudinal follow-up saliva samples in sra span approximately 1.34-5.57 years",
        "not provided", "unspecified"
      ) ~ NA_character_,

      # age ranges with "years" suffix
      str_detect(age_info_lower, "\\d+\\.?\\d*\\s*[-–—]\\s*\\d+\\.?\\d*\\s*years") ~ case_when(
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 0.5 ~ "0-6 months",
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 3 ~ "6-36 months",
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 12 ~ "3-12 years",
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 18 ~ "12-18 years",
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 60 ~ "18-60 years",
        TRUE ~ "＞60 years"
      ),

      # >= / ≥ thresholds
      str_detect(age_info_lower, "^\\s*[≥>]=?\\s*\\d+") ~ case_when(
        as.numeric(str_extract(age_info_lower, "\\d+")) >= 60 ~ "＞60 years",
        as.numeric(str_extract(age_info_lower, "\\d+")) >= 18 ~ "18-60 years",
        as.numeric(str_extract(age_info_lower, "\\d+")) > 12 ~ "12-18 years",
        as.numeric(str_extract(age_info_lower, "\\d+")) > 3 ~ "3-12 years",
        as.numeric(str_extract(age_info_lower, "\\d+")) > 0.5 ~ "6-36 months",
        TRUE ~ "0-6 months"
      ),

      # above / over X years
      str_detect(age_info_lower, "above\\s+\\d+\\s*years?\\s*old|over\\s+\\d+\\s*years?|above\\s+\\d+\\s*years?") ~ case_when(
        as.numeric(str_extract(age_info_lower, "\\d+")) >= 60 ~ "＞60 years",
        as.numeric(str_extract(age_info_lower, "\\d+")) >= 18 ~ "18-60 years",
        as.numeric(str_extract(age_info_lower, "\\d+")) > 12 ~ "12-18 years",
        as.numeric(str_extract(age_info_lower, "\\d+")) > 3 ~ "3-12 years",
        as.numeric(str_extract(age_info_lower, "\\d+")) > 0.5 ~ "6-36 months",
        TRUE ~ "0-6 months"
      ),

      # median age
      str_detect(age_info_lower, "median age.*\\d+") ~ case_when(
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 0.5 ~ "0-6 months",
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 3 ~ "6-36 months",
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 12 ~ "3-12 years",
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 18 ~ "12-18 years",
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 60 ~ "18-60 years",
        TRUE ~ "＞60 years"
      ),

      # average/mean ... years
      str_detect(age_info_lower, "(average|mean).*\\d+\\.?\\d*\\s*years") & !str_detect(age_info_lower, "±|\\+-") ~ case_when(
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 0.5 ~ "0-6 months",
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 3 ~ "6-36 months",
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 12 ~ "3-12 years",
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 18 ~ "12-18 years",
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 60 ~ "18-60 years",
        TRUE ~ "＞60 years"
      ),

      # simple numeric range without "years"
      str_detect(age_info_lower, "^\\s*\\d+\\.?\\d*\\s*[-–—]\\s*\\d+\\.?\\d*\\s*$") ~ case_when(
        (as.numeric(str_extract(age_info_lower, "^\\s*\\d+\\.?\\d*")) + as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*\\s*$"))) / 2 <= 0.5 ~ "0-6 months",
        (as.numeric(str_extract(age_info_lower, "^\\s*\\d+\\.?\\d*")) + as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*\\s*$"))) / 2 <= 3 ~ "6-36 months",
        (as.numeric(str_extract(age_info_lower, "^\\s*\\d+\\.?\\d*")) + as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*\\s*$"))) / 2 <= 12 ~ "3-12 years",
        (as.numeric(str_extract(age_info_lower, "^\\s*\\d+\\.?\\d*")) + as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*\\s*$"))) / 2 <= 18 ~ "12-18 years",
        (as.numeric(str_extract(age_info_lower, "^\\s*\\d+\\.?\\d*")) + as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*\\s*$"))) / 2 <= 60 ~ "18-60 years",
        TRUE ~ "＞60 years"
      ),

      # Exact String Matches -> 18-60 years
      age_info_lower %in% c(
        "adult male betel-quid users",
        "adult male non-smokers; cohort mean age 44 years (range 31-60)",
        "adult participants; eligibility age >=18 years; total cohort median age 40.0 years (iqr 31.0-50.0)",
        "adults; group means reported separately by disease severity",
        "control 30.50±5.74; osa 40.11±11.25; htn 36.22±10.20; osa+htn 42.17±9.52; male adults >18 years",
        "control non-cp: 34.7 ± 15.9 years; control cp: 46.3 ± 13.2 years; sle non-cp: 31.2 ± 5.2 years; sle cp: 40.5 ± 10.1 years",
        "healthy dental students; mean age 26.36 +/- 1.58 years",
        "adult participants; eligibility age >=18 years; total cohort median age 40.0 years (IQR 31.0-50.0)"
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

      # Average/Mean mapping (without "years" suffix, no ±/;)
      str_detect(age_info_lower, "average|mean") & !str_detect(age_info_lower, "±|\\+-|;") & str_detect(age_info_lower, "\\d+") ~ case_when(
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

      # 60+ / 50+ / general X+ mappings
      str_detect(age_info_lower, "\\d+\\+") ~ case_when(
        as.numeric(str_extract(age_info_lower, "\\d+")) >= 60 ~ "＞60 years",
        as.numeric(str_extract(age_info_lower, "\\d+")) >= 18 ~ "18-60 years",
        as.numeric(str_extract(age_info_lower, "\\d+")) > 12 ~ "12-18 years",
        as.numeric(str_extract(age_info_lower, "\\d+")) > 3 ~ "3-12 years",
        as.numeric(str_extract(age_info_lower, "\\d+")) > 0.5 ~ "6-36 months",
        TRUE ~ "0-6 months"
      ),
      str_detect(age_info_lower, "above 60") ~ "＞60 years",

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
      age_info_lower == "infant" ~ "0-6 months",
      age_info_lower %in% c("children", "child", "chilren") ~ "3-12 years",
      age_info_lower %in% c("elderly adults", "elderly people and adults", "elder") ~ "＞60 years",
      age_info_lower == "baby" ~ "0-6 months",
      age_info_lower %in% c("teen", "teenager") ~ "12-18 years",
      age_info_lower == "young adults" ~ "18-60 years",
      age_info_lower == "age_group_1" ~ NA_character_,
      age_info_lower == "full cohort median age 41.5 years; healthy 38 (20-66), gingivitis 43 (18-71), chronic periodontitis 43 (21-60)" ~ NA_character_,
      age_info_lower == "toddlers 1-3 years; children 4-12 years; teenager-adults 13-59 years" ~ NA_character_,

      # Variance: mean A +- B
      str_detect(age_info_lower, "mean.*\\+-") & str_detect(age_info_lower, "\\d+") ~ case_when(
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 18 ~ "12-18 years",
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 60 ~ "18-60 years",
        TRUE ~ "＞60 years"
      ),

      # between X and Y years
      str_detect(age_info_lower, "between\\s+\\d+\\.?\\d*\\s+and\\s+\\d+\\.?\\d*\\s*years") ~ case_when(
        (as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) + as.numeric(str_extract(age_info_lower, "(?<=and\\s)\\d+\\.?\\d*"))) / 2 <= 0.5 ~ "0-6 months",
        (as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) + as.numeric(str_extract(age_info_lower, "(?<=and\\s)\\d+\\.?\\d*"))) / 2 <= 3 ~ "6-36 months",
        (as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) + as.numeric(str_extract(age_info_lower, "(?<=and\\s)\\d+\\.?\\d*"))) / 2 <= 12 ~ "3-12 years",
        (as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) + as.numeric(str_extract(age_info_lower, "(?<=and\\s)\\d+\\.?\\d*"))) / 2 <= 18 ~ "12-18 years",
        (as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) + as.numeric(str_extract(age_info_lower, "(?<=and\\s)\\d+\\.?\\d*"))) / 2 <= 60 ~ "18-60 years",
        TRUE ~ "＞60 years"
      ),

      # catch-all for remaining text with numbers
      str_detect(age_info_lower, "\\d+\\.?\\d*\\s*[-–—]\\s*\\d+\\.?\\d*") ~ case_when(
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 0.5 ~ "0-6 months",
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 3 ~ "6-36 months",
        as.numeric(str_extract(age_info_lower, "\\d+\\.?\\d*")) <= 12 ~ "3-12 years",
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
      country == "Korean" ~ "Republic of Korea",
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
      site_lower == tolower("right maxillary dentition") ~ "tooth",
      site_lower %in% tolower(c("gingival surface", "Right post maxilla ging")) ~ "gingiva",
      site_lower %in% tolower(c("buccal surfaces", "mucosal surface of the soft palate", "oral mucosal surfaces")) ~ "oral mucosa",
      site_lower %in% tolower(c("the anterior dorsal tongue surface", "tongue tumor surface")) ~ "tongue",
      site_lower %in% tolower(c("cancer surface", "swab_nasopharynx_or_oropharynx")) ~ "other",
      site_lower == tolower("tonsil surface") ~ "Other",

      str_detect(site_lower, "feces|faeces|fecal|stool|rectal swab|rectum|^gut$") ~ "feces",
      str_detect(site_lower, "skin|face|forehead|palm|scalp|axilla|cheek-nose crease|umbilicus|antecubital fossa|retroauricular crease|sebum|periumbilical|forearm|areola") ~ "skin",
      str_detect(site_lower, "nasal cavity|anterior nares|nasal swab|nasal mucosa|nasal mucus|nostril mucus|nasopharyngeal|nasopharynx") ~ "nasal",
      str_detect(site_lower, "vagina|cervical mucus|posterior fornix|mid vagina") ~ "vaginal",
      str_detect(site_lower, "^blood$|bronchoalveolar lavage|sputum|endotracheal aspirate") ~ "respiratory/blood",
      str_detect(site_lower, "tears|milk|ear wax|gastric fluid") ~ "other",

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
      site_lower %in% tolower(c(
        "swab_nasopharynx_or_oropharynx", "buccal surfaces", "cancer surface",
        "gingival surface", "mucosal surface of the soft palate", "oral mucosal surfaces",
        "right maxillary dentition", "Right post maxilla ging",
        "the anterior dorsal tongue surface", "tongue tumor surface", "tonsil surface"
      )) ~ "oral",
      TRUE ~ site_group
    ),

    status_lower = tolower(status),

    disease_study_cleaned = case_when(
      bioproject %in% c("PRJEB38610", "PRJEB39420", "PRJEB39436", "PRJNA417321", "PRJNA451534", "PRJNA876796") & tolower(disease_study) == "no" ~ NA_character_,
      bioproject %in% c("PRJNA380250", "PRJNA685840", "PRJNA847970", "PRJNA762545") ~ "yes",
      bioproject %in% c("PRJNA518756", "PRJNA613363", "PRJNA782768") ~ "no",
      status_lower %in% tolower(c("preterm", "preterm newborns")) ~ "yes",
      TRUE ~ disease_study
    ),
    
    status_temp = case_when(
     is.na(status) & bioproject %in% c(
       "PRJEB1799", "PRJEB20889", "PRJEB20890", "PRJEB20894", "PRJEB27304",
       "PRJEB28088", "PRJEB28769", "PRJEB33119", "PRJEB42056", "PRJEB6518",
       "PRJNA255332", "PRJNA395436", "PRJNA427101", "PRJNA634162", "PRJNA689848",
       "PRJNA721325", "PRJNA747639", "PRJNA778070", "PRJNA828039", "PRJNA841067"
      ) ~ "HC",

      is.na(status) ~ NA_character_,

      # NA / Control categories
      status_lower %in% tolower(c(
        "mucositis healed", "pre-conditioning", "atopic non-asthmatic controls",
        "allogeneic hematopoietic stem cell transplantation recipient",
        "immunological ART responder", "dentally compromised",
        "black dental stain", "myasthenia gravis",
        "functional gastrointestinal complaints", "missing"
      )) ~ NA_character_,

      status_lower == "icu discharge or on icu day 10" ~ "ICU discharge or on ICU day 10",
      status_lower == "full-term" ~ "HC",
      str_detect(status_lower, "tonsillar hyperplasia|tonsillar hypertrophy") ~ NA_character_,
      str_detect(status_lower, "breastfed|never breastfed|full-term|orthodontic|non-responder|high-risk for ra|third molar extraction") ~ NA_character_,
      status_lower %in% tolower(c("never smoker", "current smoker")) ~ NA_character_,
      bioproject == "PRJNA417321" & status == "dental implant" ~ NA_character_,
      bioproject == "PRJNA841771" & status == "oral benign tumors" ~ NA_character_,

      # Other category (rare / non-target diseases)
      status_lower %in% tolower(c(
        "halitosis", "dysphagia", "watery diarrhea NOS", "truncal ataxia NOS",
        "oral benign tumors", "cancer patient", "Carcinoma in situ",
        "Chronic mucosal inflammation", "acute pancreatitis", "aortic thrombosis",
        "tibial fracture", "hip fracture", "hypotension", "STEMI", "NSTEMI",
        "temporal lobe epilepsy", "Fanconi anemia with leukoplakia",
        "Fanconi anemia with OSCC", "systemic lupus erythematosus with chronic periodontitis",
        "systemic lupus erythematosus", "acute lymphoblastic leukemia",
        "helminth infection", "stomatitis",
        "early childhood caries and congenital heart disease",
        "proliferative verrucous leukoplakia", "Refractory Celiac disease",
        "IBD+Diabetes", "acute diarrhea", "bradycardia NOS", "juvenile idiopathic arthritis",
        "chronic tonsillitis", "liver cirrhosis", "chronic kidney disease",
        "immunoglobulin A nephropathy", "membranous nephropathy",
        "liver carcinoma", "gliomas", "sub-acute stage of stroke with dysphagia"
      )) ~ "Other",

      # Kept as NA: conditions restored from Other back to NA
      status_lower %in% tolower(c(
        "pre-cancer", "epithelial precursor lesion", "pulmonary nodule",
        "thyroid nodules", "thyroid nodules patient",
        "primary biliary cholangitis", "autism spectrum disorder",
        "chronic fatigue syndrome", "chronic recurrent multifocal osteomyelitis",
        "severe congenital neutropenia", "methamphetamine addiction",
        "pre-eclampsia", "down syndrome", "migraine",
        "polycystic ovary syndrome"
      )) ~ NA_character_,

      # Preterm
      status_lower %in% tolower(c("preterm", "preterm infant", "preterm newborns")) ~ "Preterm",

      # ==============================================================================
      # CRITICAL FIX: Cancers (Correct Order & Capitalization)
      # ==============================================================================
      str_detect(status_lower, "head and neck squamous cell carcinoma|hnscc|head and neck cancer|hnc|nasopharyngeal carcinoma") ~ "Head and Neck Cancer",
      str_detect(status_lower, "oscc|oral squamous cell carcinoma|oral cancer") ~ "OSCC",
      str_detect(status_lower, "esophageal squamous cell carcinoma|escc|esophageal cancer") ~ "Esophageal Cancer",
      str_detect(status_lower, "thyroid cancer|differentiated thyroid carcinoma|dtc") ~ "Thyroid Cancer",
      str_detect(status_lower, "colon cancer|colorectal cancer|crc") ~ "Colorectal Cancer",
      str_detect(status_lower, "lung cancer") ~ "Lung Cancer",
      str_detect(status_lower, "breast cancer") ~ "Breast Cancer",
      str_detect(status_lower, "pancreatic cancer") ~ "Pancreatic Cancer",
      str_detect(status_lower, "digestive tract cancer") ~ "Digestive Tract Cancer",
      str_detect(status_lower, "prostate cancer") ~ "Prostate Cancer",
      str_detect(status_lower, "cancer") ~ "Other Cancer", 

      # ==============================================================================
      # CRITICAL FIX: OPMD & Local Diseases (No merging, exact capitalization)
      # ==============================================================================
      str_detect(status_lower, "oral verrucous hyperplasia") ~ "Oral Verrucous Hyperplasia",
      str_detect(status_lower, "oral submucous fibrosis") ~ "Oral Submucous Fibrosis",
      str_detect(status_lower, "olp|oral lichen planus") ~ "Oral Lichen Planus",
      status_lower == "chronic periodontitis associated with type 2 diabetes mellitus" ~ "Periodontal Disease",
      str_detect(status_lower, "periodontitis|gingivitis|periodontal|peri-implantitis") ~ "Periodontal Disease",
      str_detect(status_lower, "caries") ~ "Dental Caries",
      str_detect(status_lower, "aphthous|rau") ~ "Aphthous",
      str_detect(status_lower, "oral mucositis|severe-om|grade 1 om") ~ "Oral Mucositis",
      str_detect(status_lower, "cleft lip and palate|clp") ~ "Cleft Lip and Palate",
      
      # Cerebral Palsy Caries mapping
      status_lower == "cerebral palsy caries-free" ~ "Other",
      status_lower %in% tolower(c("cerebral palsy with mild-moderate caries", "cerebral palsy with severe caries")) ~ "Dental Caries",

      # ==============================================================================
      # CRITICAL FIX: Systemic / Autoimmune / GI (Preserving exact acronyms)
      # ==============================================================================
      str_detect(status_lower, "functional dyspepsia and gastroesophageal reflux disease|gerd|reflux esophagitis") ~ "Gastroesophageal Reflux Disease",
      str_detect(status_lower, "cd and uc|uc|ulcerative colitis|cd|crohn's disease|ibd|inflammatory bowel disease") ~ "IBD",
      status_lower == "ibd+c.diff" ~ "IBD", 
      status_lower == "c.diff" ~ "Clostridium Difficile Infection",
      str_detect(status_lower, "type 2 diabetes mellitus|t2dm|type-2 diabetes|type 1 diabetes|diabetes") ~ "Diabetes",
      str_detect(status_lower, "obese|obesity") ~ "Obesity",
      str_detect(status_lower, "atherosclerosis") ~ "Cardiovascular Disease",
      str_detect(status_lower, "coronavirus disease of 2019|covid") & !str_detect(status_lower, "discharge|convalescent") ~ "COVID-19",
      status_lower == "covid-19 discharge" ~ "COVID-19 discharge",
      status_lower == "1-year covid-19 convalescent" ~ "1-year COVID-19 convalescent",
      str_detect(status_lower, "hiv|viremic art-untreated") ~ "HIV",
      str_detect(status_lower, "axspa|axial spondyloarthritis") ~ "Axial Spondyloarthritis",
      str_detect(status_lower, "sjögren's syndrome|sjögren’s syndrome") ~ "Sjögren's Syndrome",
      str_detect(status_lower, "ad|alzheimer") ~ "Alzheimer's Disease",
      str_detect(status_lower, "parkinson") ~ "Parkinson's Disease",
      str_detect(status_lower, "copd") ~ "COPD",
      str_detect(status_lower, "obstructive sleep apnea") ~ "Obstructive Sleep Apnea",
      str_detect(status_lower, "post-traumatic stress disorder|ptsd") ~ "PTSD",
      str_detect(status_lower, "celiac disease") ~ "Celiac Disease",
      str_detect(status_lower, "alcohol use disorder") ~ "Alcohol Use Disorder",
      str_detect(status_lower, "atopic dermatitis") ~ "Atopic Dermatitis",
      str_detect(status_lower, "atopic asthma|asthma") ~ "Asthma",
      str_detect(status_lower, "chagas") ~ "Chagas Disease",
      str_detect(status_lower, "helicobacter pylori") ~ "Helicobacter pylori Infection",

      # Bioproject-specific overrides & direct HC matches
      status_lower %in% c("healthy", "no-om") ~ "HC",
      bioproject == "PRJNA316469" & disease_study == "no" & status %in% c("never smoker", "current smoker") ~ "HC",
      bioproject == "PRJEB40804" & status == "certified saturation divers with valid health certificates" ~ "HC",
      bioproject == "PRJNA656387" & status == "first-degree relatives" ~ "HC",
      bioproject == "PRJDB8324" & status == "healthy" ~ "HC",
      bioproject == "PRJNA663926" & status %in% c("hand, foot and mouth disease symptomatic", "hand, foot and mouth disease asymptomatic") ~ "Hand, Foot and Mouth Disease",

      TRUE ~ tools::toTitleCase(status) # Fallback for anything not caught above
    ),
    
    # Directly assign status_temp without doing another toTitleCase which destroys Acronyms
    status_cleaned = ifelse(!is.na(status_temp) & tolower(status_temp) != "na", status_temp, NA_character_)
  ) %>%
  select(-age_cleaned_group_base, -age_info_lower, -status_temp) # Clean up temp variables

# Apply Uniform Title Case processing, BUT exclude status_cleaned to protect acronyms
df_new <- df_new %>%
  mutate(across(c(sex_cleaned, platform_cleaned, sequencing_region_cleaned,
                  country_cleaned, continent_group, site_group_cleaned,
                  disease_study_cleaned),
                ~ ifelse(!is.na(.), tools::toTitleCase(as.character(.)), NA_character_)))

# Apply TitleCase to site_cleaned BUT keep standard categories consistent
df_new <- df_new %>%
  mutate(site_cleaned = ifelse(!is.na(site_cleaned), tools::toTitleCase(as.character(site_cleaned)), NA_character_))

# Fix specific categories that should not be TitleCase-d or need correction
fix_site_case <- function(x) {
  x <- case_when(
    tolower(x) == "feces" ~ "feces",
    tolower(x) == "skin" ~ "skin",
    tolower(x) == "nasal" ~ "nasal",
    tolower(x) == "vaginal" ~ "vaginal",
    tolower(x) == "respiratory/blood" ~ "respiratory/blood",
    tolower(x) == "other" ~ "other",
    tolower(x) == "oral mucosa" ~ "oral mucosa",
    tolower(x) == "oral swab" ~ "oral swab",
    tolower(x) == "saliva" ~ "saliva",
    tolower(x) == "tooth" ~ "tooth",
    tolower(x) == "tongue" ~ "tongue",
    tolower(x) == "gingiva" ~ "gingiva",
    tolower(x) == "other" ~ "Other",  # preserve original "Other" for tonsil etc.
    TRUE ~ x
  )
  x
}

df_new <- df_new %>%
  mutate(site_cleaned = fix_site_case(site_cleaned))

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

# --- NEW: fill subject_id NA with sample_id for specific bioprojects ---
subject_id_fill_bioprojects <- c(
  "PRJEB35824", "PRJNA779698", "PRJNA876796", "PRJEB36081",
  "PRJNA820027", "PRJNA545410", "PRJNA770171", "PRJNA663437"
)
df_new <- df_new %>%
  mutate(subject_id = ifelse(
    bioproject %in% subject_id_fill_bioprojects & is.na(subject_id),
    sample_id,
    subject_id
  ))

# ------------------------------------------------------------------------------
# 3. Reorder Columns and Save Output
# ------------------------------------------------------------------------------
message("Reordering columns and saving output dataset...")

# Define cleaned columns in a logical order
cleaned_cols <- c(
  "age_cleaned", "bmi_cleaned", "age_cleaned_group", "age_info_cleaned",
  "bmi_cleaned_group", "sex_cleaned", "platform_cleaned", "sequencing_region_cleaned",
  "country_cleaned", "continent_group", "site_cleaned", "site_group_cleaned",
  "disease_study_cleaned", "status_cleaned"
)

# Final output: required_cols + cleaned_cols ONLY (drop all raw_* columns)
final_cols <- c(required_cols, cleaned_cols)
final_cols <- final_cols[final_cols %in% names(df_new)]

df_final <- df_new %>% select(all_of(final_cols))
write_csv(df_final, output_file)

message(sprintf("\nPipeline execution complete! Output saved to: %s\n", output_file))
message(sprintf("Rows: %d | Columns: %d\n", nrow(df_final), ncol(df_final)))

# ------------------------------------------------------------------------------
# 4. Export Oral-Only Subset (matching original oral script filtering logic)
# ------------------------------------------------------------------------------
oral_output_file <- "/home/yanglun/YSD/global_oral/metadata_v2_merge/oral-v2/results/metadata-v2-merge-oral.csv"

df_oral <- df_final %>%
  filter(
    tolower(site_group_cleaned) == "oral" &
      !(bioproject == "PRJNA471556" & tolower(site) == "stool" & tolower(site_group_cleaned) == "oral")
  )

write_csv(df_oral, oral_output_file)
message(sprintf("Oral-only subset saved to: %s\n", oral_output_file))
message(sprintf("Oral subset rows: %d | Columns: %d\n", nrow(df_oral), ncol(df_oral)))