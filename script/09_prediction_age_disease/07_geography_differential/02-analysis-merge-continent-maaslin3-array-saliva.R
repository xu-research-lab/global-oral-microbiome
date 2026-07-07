# ==============================================================================
# 01_merge_maaslin3_results.R
# Merge results from all MaAsLin3 array tasks
#
# This script collects and combines results from all pairwise comparisons
# after the array job completes.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# ==========================================================================
# Configuration
# ==========================================================================
results_dir <- "results/maaslin3_pairwise"
alpha <- 0.05
continent_col <- "continent_group"

cat("==========================================\n")
cat("Merging MaAsLin3 Array Results\n")
cat("Start: ", format(Sys.time()), "\n")
cat("==========================================\n\n")

# ==========================================================================
# 1. Find all completed comparisons
# ==========================================================================
cat("[1/3] Scanning for completed comparisons...\n")

comparison_dirs <- list.dirs(results_dir, recursive = FALSE, full.names = TRUE)
comparison_names <- basename(comparison_dirs)

# Filter to directories with completed analysis
completed <- sapply(comparison_dirs, function(d) {
  file.exists(file.path(d, "analysis.done"))
})

comparison_dirs <- comparison_dirs[completed]
comparison_names <- comparison_names[completed]

n_completed <- length(comparison_dirs)
cat("    Found ", n_completed, " completed comparisons\n", sep = "")

if (n_completed == 0) {
  stop("No completed comparisons found!")
}

# ==========================================================================
# 2. Merge all results
# ==========================================================================
cat("\n[2/3] Merging results...\n")

all_results_list <- list()

for (i in seq_along(comparison_dirs)) {
  comp_dir <- comparison_dirs[i]
  comp_name <- comparison_names[i]

  cat("    [", i, "/", n_completed, "] ", comp_name, "\n", sep = "")

  # Read full results (MaAsLin3 writes with quote=TRUE; fread strips quotes
  # automatically so column names come back clean: qval_individual / qval_joint)
  res_file <- file.path(comp_dir, "all_results.tsv")
  if (file.exists(res_file)) {
    res_df <- fread(res_file, sep = "\t", header = TRUE, data.table = FALSE)

    # If the array R script forgot to attach Continent1/2/Comparison, attach now
    # by parsing the directory name "<C1>_vs_<C2>".
    if (!"Comparison" %in% names(res_df)) res_df$Comparison <- comp_name
    if (!"Continent1" %in% names(res_df) || !"Continent2" %in% names(res_df)) {
      parts <- strsplit(comp_name, "_vs_", fixed = TRUE)[[1]]
      if (length(parts) == 2) {
        res_df$Continent1 <- parts[1]
        res_df$Continent2 <- parts[2]
      }
    }
    all_results_list[[comp_name]] <- res_df
  }
}

# Combine all results
if (length(all_results_list) == 0) stop("No all_results.tsv files found!")

full_results <- bind_rows(all_results_list)
rownames(full_results) <- NULL

all_results_file <- file.path(results_dir, "all_pairwise_results.tsv")
fwrite(full_results, all_results_file, sep = "\t", quote = FALSE)
cat("\n    Saved: ", all_results_file, "\n", sep = "")

# ==========================================================================
# 2b. Recompute summary counts directly from full_results
# ----------------------------------------------------------------------------
# We DO NOT trust per-pair summary.tsv files because the original array R
# script filtered on a non-existent `qval` column (MaAsLin3 outputs
# qval_individual / qval_joint) and did not restrict to continent_group rows.
# Compute the corrected N_Significant here instead.
# ==========================================================================
qcol <- if ("qval_individual" %in% names(full_results)) "qval_individual" else "qval"
cat("\n    Using q-value column: ", qcol, "\n", sep = "")

is_continent_row <- if ("metadata" %in% names(full_results)) {
  full_results[["metadata"]] == continent_col
} else {
  rep(TRUE, nrow(full_results))
}

sig_all <- full_results[is_continent_row &
                          !is.na(full_results[[qcol]]) &
                          full_results[[qcol]] < alpha, ]

# Save the corrected significant-only table for downstream plots
sig_file <- file.path(results_dir, "significant_pairwise_results.tsv")
fwrite(sig_all, sig_file, sep = "\t", quote = FALSE)
cat("    Saved: ", sig_file, " (", nrow(sig_all), " significant rows)\n", sep = "")

summary_counts <- sig_all %>%
  group_by(Continent1, Continent2, Comparison) %>%
  summarise(N_Significant = n_distinct(feature), .groups = "drop") %>%
  as.data.frame()

# Make sure every comparison appears, even if it has 0 significant taxa
all_comparisons <- unique(data.frame(
  Comparison = sapply(all_results_list, function(d) d$Comparison[1]),
  Continent1 = sapply(all_results_list, function(d) d$Continent1[1]),
  Continent2 = sapply(all_results_list, function(d) d$Continent2[1]),
  stringsAsFactors = FALSE
))
summary_counts <- merge(all_comparisons, summary_counts,
                        by = c("Continent1", "Continent2", "Comparison"),
                        all.x = TRUE)
summary_counts$N_Significant[is.na(summary_counts$N_Significant)] <- 0

# Symmetric form (needed by Figure 5A bubble matrix)
summary_symmetric <- summary_counts
for (i in seq_len(nrow(summary_counts))) {
  summary_symmetric <- rbind(summary_symmetric, data.frame(
    Continent1 = summary_counts$Continent2[i],
    Continent2 = summary_counts$Continent1[i],
    Comparison = paste0(summary_counts$Continent2[i], "_vs_", summary_counts$Continent1[i]),
    N_Significant = summary_counts$N_Significant[i],
    stringsAsFactors = FALSE
  ))
}

summary_file_out <- file.path(results_dir, "summary_diff_counts.tsv")
fwrite(summary_symmetric, summary_file_out, sep = "\t", quote = FALSE)
cat("    Saved: ", summary_file_out, "\n", sep = "")

# ==========================================================================
# 3. Print summary
# ==========================================================================
cat("\n[3/3] Summary statistics...\n")

cat("\n==========================================\n")
cat("Pairwise Significant Taxa Counts (", qcol, " < ", alpha, ", metadata=", continent_col, "):\n", sep = "")
cat("==========================================\n")
print(summary_counts)

cat("\n    Total comparisons: ", nrow(summary_counts), "\n", sep = "")
cat("    Total associations tested: ", nrow(full_results), "\n", sep = "")
cat("    Mean diff taxa per comparison: ", round(mean(summary_counts$N_Significant), 1), "\n", sep = "")
cat("    Range: ", min(summary_counts$N_Significant), " - ", max(summary_counts$N_Significant), "\n", sep = "")

cat("\n==========================================\n")
cat("Merge Complete\n")
cat("End: ", format(Sys.time()), "\n")
cat("==========================================\n")
