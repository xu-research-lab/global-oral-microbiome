# ==============================================================================
# 01_analysis_oral_continent_maaslin3_array.R
# MaAsLin3 pairwise continent comparisons - Array Job Version
#
# This script is designed for Slurm array jobs. Each task processes ONE
# continent pair comparison.
#
# Usage: Rscript script.R <task_id>
#   task_id: 1-15 (one for each pairwise comparison)
# ==============================================================================

suppressPackageStartupMessages({
  library(qiime2R)
  library(maaslin3)
  library(dplyr)
  library(data.table)
  library(Matrix)
})

# ==========================================================================
# 0. Get task ID from command line
# ==========================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript script.R <task_id>")
}
task_id <- as.integer(args[1])

set.seed(20260615)

# ==========================================================================
# 1. Configuration
# ==========================================================================
input_dir    <- "input"
results_dir  <- "results/maaslin3_pairwise"

genus_qza    <- file.path(input_dir, "genus-table.qza")
metadata_file <- "/home/yanglun/YSD/global_oral/global_oral_v2/metadata-v2/metadata_oral_qc.tsv"

# --- ANALYSIS PARAMETERS ---
target_site <- "saliva"

# 修正：将 sample_id_col 设置为 metadata 中的 sample_id 列
sample_id_col <- "sample_id"
continent_col <- "continent_group"
site_col      <- "site_cleaned"
age_col       <- "age_cleaned"
sex_col       <- "sex_cleaned"
bmi_col       <- "bmi_cleaned"
project_col   <- "bioproject_id"
platform_col  <- "platform_cleaned"
seq_region_col <- "sequencing_region_cleaned"
pcr_col       <- "pcr_primers"

continent_order <- c("Africa", "Asia", "Europe",
                     "North America", "Oceania", "South America")

min_samples_per_continent <- 1

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
alpha <- 0.05

# ==========================================================================
# 2. Load and prepare data
# ==========================================================================
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

cat("==========================================\n")
cat("MaAsLin3 Array Task: ", task_id, "\n")
cat("Start: ", format(Sys.time()), "\n")
cat("==========================================\n\n")

cat("[1/4] Loading data...\n")

# Load genus abundance table
qza_obj <- read_qza(genus_qza)
genus_sparse <- qza_obj$data

# Load metadata
metadata_raw <- fread(metadata_file, sep = "\t", header = TRUE,
                      quote = "", fill = TRUE, check.names = FALSE,
                      data.table = FALSE)

# ==========================================================================
# 3. Filter and clean metadata
# ==========================================================================
cat("\n[2/4] Filtering metadata...\n")

# Only require structural columns; covariate columns are handled with "Unknown" fill
required_cols <- c(sample_id_col, continent_col, site_col, project_col)
missing_cols <- setdiff(required_cols, names(metadata_raw))
if (length(missing_cols) > 0) {
  stop("Metadata missing required columns: ", paste(missing_cols, collapse = ", "))
}

metadata_clean <- metadata_raw %>%
  filter(!is.na(.data[[sample_id_col]]) & .data[[sample_id_col]] != "") %>%
  distinct(.data[[sample_id_col]], .keep_all = TRUE) %>%
  mutate(
    # Extract continuous variables to temp columns to prevent data loss
    temp_age = suppressWarnings(as.numeric(as.character(.data[[age_col]]))),
    temp_bmi = suppressWarnings(as.numeric(as.character(.data[[bmi_col]]))),
    !!site_col := trimws(as.character(.data[[site_col]])),
    !!sex_col := trimws(as.character(.data[[sex_col]])),
    !!continent_col := trimws(as.character(.data[[continent_col]])),
    !!platform_col := trimws(as.character(.data[[platform_col]])),
    !!seq_region_col := trimws(as.character(.data[[seq_region_col]])),
    !!pcr_col := trimws(as.character(.data[[pcr_col]]))
  ) %>%
  # Categorize BMI (Chinese Standard) and handle NAs
  mutate(
    !!bmi_col := case_when(
      is.na(temp_bmi) ~ "Unknown",
      temp_bmi < 18.5 ~ "Underweight",
      temp_bmi >= 18.5 & temp_bmi < 24.0 ~ "Normal",
      temp_bmi >= 24.0 & temp_bmi < 28.0 ~ "Overweight",
      temp_bmi >= 28.0 ~ "Obese",
      TRUE ~ "Unknown"
    ),
    !!age_col := case_when(
      is.na(temp_age) ~ "Unknown",
      temp_age < 18 ~ "Minor",
      temp_age >= 18 & temp_age < 45 ~ "Young_Adult",
      temp_age >= 45 & temp_age < 65 ~ "Middle_Aged",
      temp_age >= 65 ~ "Elderly",
      TRUE ~ "Unknown"
    )
  ) %>%
  # Convert categorized variables to factors
  mutate(
    !!bmi_col := factor(.data[[bmi_col]], levels = c("Normal", "Underweight", "Overweight", "Obese", "Unknown")),
    !!age_col := factor(.data[[age_col]], levels = c("Young_Adult", "Minor", "Middle_Aged", "Elderly", "Unknown"))
  ) %>%
  # Replace missing categorical values with "Unknown"
  mutate(
    !!sex_col := ifelse(.data[[sex_col]] %in% c("", "NA"), "Unknown", .data[[sex_col]]),
    !!platform_col := ifelse(.data[[platform_col]] %in% c("", "NA"), "Unknown", .data[[platform_col]]),
    !!seq_region_col := ifelse(.data[[seq_region_col]] %in% c("", "NA"), "Unknown", .data[[seq_region_col]]),
    !!pcr_col := ifelse(.data[[pcr_col]] %in% c("", "NA"), "Unknown", .data[[pcr_col]])
  ) %>%
  filter(
    !is.na(.data[[continent_col]]) &
    .data[[continent_col]] != "" & .data[[continent_col]] != "NA" &
    .data[[site_col]] == target_site
  ) %>%
  select(-temp_age, -temp_bmi)

# Critical: dplyr pipelines return tibbles which DO NOT support rownames.
# MaAsLin3 requires a base data.frame with valid rownames.
metadata_clean <- as.data.frame(metadata_clean)
rownames(metadata_clean) <- metadata_clean[[sample_id_col]]

# Find common samples
common_samples <- intersect(rownames(metadata_clean), colnames(genus_sparse))
if (length(common_samples) == 0) {
  stop("No common samples between metadata and genus table!")
}

metadata_clean <- metadata_clean[common_samples, ]
genus_sparse <- genus_sparse[, common_samples]

# Identify valid continents
continent_counts <- table(metadata_clean[[continent_col]])
valid_continents <- names(continent_counts)[continent_counts >= min_samples_per_continent]
valid_continents <- intersect(continent_order, valid_continents)

if (length(valid_continents) < 2) {
  stop("Need at least 2 continents with >= ", min_samples_per_continent, " samples!")
}

cat("    Valid continents: ", paste(valid_continents, collapse = ", "), "\n", sep = "")

# ==========================================================================
# 4. Generate all pairs and select the one for this task
# ==========================================================================
cat("\n[3/4] Identifying comparison for task ", task_id, "...\n", sep = "")

continent_pairs <- combn(valid_continents, 2, simplify = FALSE)
n_pairs_total <- length(continent_pairs)

if (task_id > n_pairs_total) {
  cat("Task ID", task_id, "exceeds total pairs (", n_pairs_total, ") - nothing to do. Exiting.\n", sep = " ")
  quit(save = "no", status = 0)
}

# Select the pair for this task
pair <- continent_pairs[[task_id]]
continent1 <- pair[1]
continent2 <- pair[2]

comparison_name <- paste0(continent1, "_vs_", continent2)
output_subdir <- file.path(results_dir, comparison_name)

cat("    Task ", task_id, "/", n_pairs_total, ": ", comparison_name, "\n", sep = "")

# Check if already completed (for resumability)
done_file <- file.path(output_subdir, "analysis.done")
if (file.exists(done_file)) {
  cat("    Already completed. Skipping.\n")
  quit(save = "no", status = 0)
}

# ==========================================================================
# 5. Run MaAsLin3 for this pair
# ==========================================================================
cat("\n[4/4] Running MaAsLin3 analysis...\n")

# Subset samples for this pair
pair_samples <- metadata_clean[[continent_col]] %in% c(continent1, continent2)
pair_meta <- metadata_clean[pair_samples, ]
pair_genus <- genus_sparse[, pair_meta[[sample_id_col]]]

# Set continent as factor with reference level
pair_meta[[continent_col]] <- factor(pair_meta[[continent_col]],
                                     levels = c(continent1, continent2))

# Prepare input for MaAsLin3
input_data <- as.data.frame(as.matrix(t(pair_genus)))
input_meta <- as.data.frame(pair_meta)
rownames(input_meta) <- input_meta[[sample_id_col]]
# Ensure input_data and input_meta have matching sample order
input_meta <- input_meta[rownames(input_data), ]

# --- Dynamic Covariate Filtering to Prevent Singular Fit ---
all_covars <- c(age_col, sex_col, platform_col, seq_region_col, bmi_col)
valid_covars <- character(0)

for (cov in all_covars) {
  vals <- as.character(input_meta[[cov]])
  tab <- sort(table(vals), decreasing = TRUE)
  
  if (length(tab) < 2) {
    cat("    WARNING: ", cov, " has only 1 level, dropping from formula.\n", sep = "")
    next
  }
  
  cross_tab <- table(input_meta[[continent_col]], input_meta[[cov]])
  if (any(colSums(cross_tab > 0) < 2)) {
      cat("    WARNING: ", cov, " has non-overlapping levels between continents, dropping to prevent Singular Fit.\n", sep = "")
      next
  }
  
  ref_level <- names(tab)[1]
  if (cov == bmi_col && "Normal" %in% names(tab)) ref_level <- "Normal"
  if (cov == age_col && "Young_Adult" %in% names(tab)) ref_level <- "Young_Adult"
  
  other_levels <- setdiff(names(tab), ref_level)
  input_meta[[cov]] <- factor(vals, levels = c(ref_level, other_levels))
  valid_covars <- c(valid_covars, cov)
  cat("    [OK] Kept ", cov, " (Ref: ", ref_level, ")\n", sep = "")
}

if (length(valid_covars) > 0) {
  base_formula <- paste0("~ ", continent_col, " + ", paste(valid_covars, collapse = " + "))
} else {
  base_formula <- paste0("~ ", continent_col)
}

project_tab <- table(input_meta[[continent_col]], input_meta[[project_col]])
shared_projects <- sum(colSums(project_tab > 0) == 2)

if (shared_projects > 0) {
  formula_str <- paste0(base_formula, " + (1|", project_col, ")")
  cat("    Using Random Effect for BioProject (", shared_projects, " shared projects found).\n", sep="")
} else {
  formula_str <- base_formula
  cat("    WARNING: No shared BioProjects between continents. Falling back to fixed-effects only (LM).\n", sep="")
}

cat("    Formula: ", formula_str, "\n", sep = "")
cat("    Samples: ", continent1, " (", sum(pair_meta[[continent_col]] == continent1), "), ",
    continent2, " (", sum(pair_meta[[continent_col]] == continent2), ")\n", sep = "")

# Run MaAsLin3
if (!dir.exists(output_subdir)) dir.create(output_subdir, recursive = TRUE)

tryCatch({
  fit_out <- maaslin3(
    input_data     = input_data,
    input_metadata = input_meta,
    output         = output_subdir,
    formula        = formula_str,
    reference      = paste0(continent_col, ",", continent1),
    normalization  = "TSS",
    transform      = "LOG",
    standardize    = TRUE,
    cores          = n_cores,
    min_prevalence = 0.05,
    plot_associations = TRUE,
    plot_summary      = TRUE
  )

  # Read and annotate results
  res_file <- file.path(output_subdir, "all_results.tsv")
  if (file.exists(res_file)) {
    res_df <- fread(res_file, sep = "\t", header = TRUE, data.table = FALSE)

    # Add comparison metadata
    res_df$Continent1 <- continent1
    res_df$Continent2 <- continent2
    res_df$Comparison <- comparison_name
    res_df$Target_Site <- target_site

    # Filter significant results
    sig_df <- res_df[!is.na(res_df[["qval"]]) & res_df[["qval"]] < alpha, ]

    if (nrow(sig_df) > 0) {
      sig_file <- file.path(output_subdir, "significant_results.tsv")
      fwrite(sig_df, sig_file, sep = "\t", quote = FALSE)
    }

    n_sig <- nrow(sig_df)
    cat("    Significant taxa (qval < ", alpha, "): ", n_sig, "\n", sep = "")

    # Save summary
    summary_file <- file.path(output_subdir, "summary.tsv")
    summary_df <- data.frame(
      Continent1 = continent1,
      Continent2 = continent2,
      Comparison = comparison_name,
      N_Significant = n_sig,
      stringsAsFactors = FALSE
    )
    fwrite(summary_df, summary_file, sep = "\t", quote = FALSE)
  }

  # Mark as completed
  writeLines("done", done_file)
  cat("\n    SUCCESS!\n")

}, error = function(e) {
  cat("    ERROR: ", e$message, "\n", sep = "")
  quit(save = "no", status = 1)
})

cat("\n==========================================\n")
cat("Task ", task_id, " Complete\n", sep = "")
cat("End: ", format(Sys.time()), "\n")
cat("==========================================\n")