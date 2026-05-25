library(dplyr)
library(readr)
library(purrr)
library(stringr)

# Define paths
project_list_input <- "D:/YL2025/2026/global-oral-github-YL-v2/global-oral-microbiome/data/project_oral_gut_pair_list.txt"
base_project_dir <- "D:/YL2025/2026/global-oral-github-YL-v2/global-oral-microbiome/data/projects"

output_dir <- "E:/oral/26/metadata-2604/Oral-Metadata-merge-YL/pair/260518/pair-v2/results-0520/merge-pair-0525"
output_file <- "metadata_oral_gut_pair_v2_merge.csv"

output_additional_file <- "E:/oral/26/metadata-2604/Oral-Metadata-merge-YL/pair/260518/pair-v2/results-0520/merge-pair-0525/metadata_pair_v2_additional_oral_site_merge.csv"

# Create output directories if they don't exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(dirname(output_additional_file))) dir.create(dirname(output_additional_file), recursive = TRUE)

# Read project list and locate valid files
project_ids <- readLines(project_list_input, warn = FALSE)
project_ids <- project_ids[trimws(project_ids) != ""]

all_files <- file.path(base_project_dir, project_ids, "metadata_oral_gut_pair_v2.csv")
valid_files <- all_files[file.exists(all_files)]

# Define required columns in correct order
target_cols <- c(
  "pair_id", "pair_time_point", "pair_bioproject", "pair_type", "pair_sample_number", "pair_subject_number", 
  "bioproject", "sample_id", "subject_id", "age", "age_info", "bmi", "sex", "platform", 
  "sequencing_region", "pcr_primers", "country", "site", "site_group", "disease_study", "status", 
  "title", "public_time", "matched_samples", "is_time_series", "age_cleaned", "bmi_cleaned", 
  "age_cleaned_group", "age_info_cleaned", "bmi_cleaned_group", "sex_cleaned", "platform_cleaned", 
  "sequencing_region_cleaned", "country_cleaned", "continent_group", "site_cleaned", 
  "site_group_cleaned", "disease_study_cleaned", "status_cleaned"
)

merged_data_list <- list()
additional_data_list <- list()

# Process each valid file
for (f in valid_files) {
  df <- read_csv(f, col_types = cols(.default = "c"), show_col_types = FALSE)
  
  if(nrow(df) > 0) {
    proj_name <- unique(df$pair_bioproject)[1]
    if(is.na(proj_name)) proj_name <- unique(df$bioproject)[1]
    
    needs_recalc <- FALSE 
    
    # Filter specific projects
    if(!is.na(proj_name) && proj_name %in% c("PRJEB69530", "PRJEB11425")) {
      df <- df %>% filter(!is.na(status_cleaned))
      needs_recalc <- TRUE
    }
    
    # Isolate 4-part pair_ids into an additional dataframe
    oral_sites_pattern <- "(saliva|tooth|oral[ _]mucosa|oral[ _]swab|tongue|gingiva|other)$"
    is_4part <- grepl(oral_sites_pattern, trimws(df$pair_id), ignore.case = TRUE)
    
    if (sum(is_4part, na.rm = TRUE) > 0) {
       df_4part <- df[is_4part, ]
       df <- df[!is_4part, ]
       needs_recalc <- TRUE
       
       available_cols_4part <- target_cols[target_cols %in% names(df_4part)]
       df_4part_filtered <- df_4part %>% select(all_of(available_cols_4part))
       for(c in setdiff(target_cols, names(df_4part_filtered))) {
          df_4part_filtered[[c]] <- NA_character_
       }
       additional_data_list[[f]] <- df_4part_filtered %>% select(all_of(target_cols))
    }
    
    # Recalculate pairing metrics if filtered
    if(needs_recalc && nrow(df) > 0) {
      valid_pairs <- df %>% group_by(pair_id) %>% summarise(count = n(), .groups='drop') %>% filter(count == 2) %>% pull(pair_id)
      df <- df %>% filter(pair_id %in% valid_pairs)
      
      if(nrow(df) > 0) {
        new_pair_sample_number <- length(unique(df$pair_id))
        new_pair_subject_number <- length(unique(df$subject_id))
        
        subj_pairs_stats <- df %>% group_by(subject_id) %>% summarise(p_count = n_distinct(pair_id), .groups='drop')
        new_pair_type <- ifelse(max(subj_pairs_stats$p_count, na.rm = TRUE) > 1, "many pairs", "single pair")
        
        df <- df %>% mutate(
          pair_sample_number = as.character(new_pair_sample_number),
          pair_subject_number = as.character(new_pair_subject_number),
          pair_type = new_pair_type
        )
      }
    }
    
    # Align columns and store main data
    if(nrow(df) > 0) {
      available_cols <- target_cols[target_cols %in% names(df)]
      df_filtered <- df %>% select(all_of(available_cols))
      
      for(c in setdiff(target_cols, names(df_filtered))) {
         df_filtered[[c]] <- NA_character_
      }
      merged_data_list[[f]] <- df_filtered %>% select(all_of(target_cols))
    }
  } 
}

# Bind lists into final dataframes
merged_data <- bind_rows(merged_data_list)
additional_data <- bind_rows(additional_data_list)

# Format target columns and type conversion
final_data <- merged_data %>% 
  select(all_of(target_cols)) %>% 
  mutate(
    age                 = as.numeric(age),
    bmi                 = as.numeric(bmi),
    age_cleaned         = as.numeric(age_cleaned),
    bmi_cleaned         = as.numeric(bmi_cleaned),
    pair_sample_number  = as.numeric(pair_sample_number),
    pair_subject_number = as.numeric(pair_subject_number)
  )

# Process and save additional 4-part data
if(nrow(additional_data) > 0) {
  final_additional_data <- additional_data %>%
    mutate(
      age                 = as.numeric(age),
      bmi                 = as.numeric(bmi),
      age_cleaned         = as.numeric(age_cleaned),
      bmi_cleaned         = as.numeric(bmi_cleaned),
      pair_sample_number  = NA_real_,
      pair_subject_number = NA_real_
    )
  write_csv(final_additional_data, output_additional_file)
}

# Save main merged CSV
write_csv(final_data, file.path(output_dir, output_file))