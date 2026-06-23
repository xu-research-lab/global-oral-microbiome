library(tidyverse)
library(fs)

# Define paths and required columns
base_dir <- "/home/yanglun/YSD/global_oral/metadata_v2_merge/oral-v2/global-oral-microbiome/data/projects/"
output_dir <- "/home/yanglun/YSD/global_oral/metadata_v2_merge/oral-v2/results/"

dir_create(output_dir)

output_csv <- file.path(output_dir, "metadata-v2-merge.csv")
log_file <- file.path(output_dir, "merge_report.txt")

required_cols <- c(
  "bioproject", "sample_id", "subject_id", "age", "age_info", 
  "bmi", "sex", "platform", "sequencing_region", "pcr_primers", 
  "country", "site", "site_group", "disease_study", "status", 
  "title", "public_time", "matched_samples", "is_time_series"
)

# Scan directories for metadata files
subdirs <- dir_ls(base_dir, type = "directory")
expected_files <- file.path(subdirs, "metadata_v2.csv")

file_exists_log <- file.exists(expected_files)
valid_files <- expected_files[file_exists_log]
missing_files_dirs <- basename(subdirs[!file_exists_log])

# Read data and track missing columns
df_list <- list()
missing_cols_report <- list()

for (file_path in valid_files) {
  folder_name <- basename(dirname(file_path))
  
  # Read as character to prevent type mismatch during merge; explicitly disable progress bar
  temp_df <- read_csv(file_path, col_types = cols(.default = "c"), show_col_types = FALSE)
  
  missing_cols <- setdiff(required_cols, colnames(temp_df))
  
  if (length(missing_cols) > 0) {
    missing_cols_report[[folder_name]] <- missing_cols
  }
  
  df_list[[folder_name]] <- temp_df
}

# Merge dataframes and ensure all 19 required columns exist
merged_df <- bind_rows(df_list)

for (col in required_cols) {
  if (!col %in% colnames(merged_df)) {
    merged_df[[col]] <- NA
  }
}

# Reorder columns: required first, others after
other_cols <- setdiff(colnames(merged_df), required_cols)
final_df <- merged_df %>% select(all_of(required_cols), all_of(other_cols))

write_csv(final_df, output_csv)

# Generate log report
report_lines <- c(
  "================ Merge Report ================",
  sprintf("Execution Time: %s", Sys.time()),
  sprintf("Total subdirectories checked: %d", length(subdirs)),
  sprintf("Successfully merged metadata_v2.csv files: %d", length(valid_files)),
  "",
  "[1. Subdirectories missing metadata_v2.csv]",
  if (length(missing_files_dirs) > 0) paste(" -", missing_files_dirs) else "None. All subdirectories contain the file.",
  "",
  "[2. Files missing required columns]"
)

if (length(missing_cols_report) > 0) {
  for (folder in names(missing_cols_report)) {
    report_lines <- c(report_lines, sprintf(" - Folder '%s' is missing: %s", folder, paste(missing_cols_report[[folder]], collapse = ", ")))
  }
} else {
  report_lines <- c(report_lines, "None. All read files contain the 19 required columns.")
}

cat(paste(report_lines, collapse = "\n"), "\n")
writeLines(report_lines, log_file)

cat(sprintf("\nMerge complete! Data saved to: %s\nReport saved to: %s\n", output_csv, log_file))