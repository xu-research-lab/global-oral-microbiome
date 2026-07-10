#!/usr/bin/env Rscript
# ============================================================================
# Merge site pairwise MaAsLin3 results
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

setwd("/beegfs/huangxiaochang/global_oral/new_analysis/maaslin3_v2/site_pairwise")

results_dir <- "results"
alpha <- 0.05
site_col <- "site_cleaned"

cat("=== Merging site pairwise results ===\n")

comparison_dirs <- list.dirs(results_dir, recursive = FALSE, full.names = TRUE)
comparison_names <- basename(comparison_dirs)

completed <- sapply(comparison_dirs, function(d) file.exists(file.path(d, "analysis.done")))
comparison_dirs <- comparison_dirs[completed]
comparison_names <- comparison_names[completed]
cat(sprintf("Found %d completed comparisons\n", length(comparison_dirs)))

if (length(comparison_dirs) == 0) stop("No completed comparisons!")

# Merge all results
all_results_list <- list()
for (i in seq_along(comparison_dirs)) {
  comp_dir <- comparison_dirs[i]
  comp_name <- comparison_names[i]
  res_file <- file.path(comp_dir, "all_results.tsv")
  if (!file.exists(res_file)) next

  res_df <- read.delim(res_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"Comparison" %in% names(res_df)) res_df$Comparison <- comp_name
  if (!"Site1" %in% names(res_df)) {
    parts <- strsplit(comp_name, "_vs_", fixed = TRUE)[[1]]
    if (length(parts) == 2) { res_df$Site1 <- parts[1]; res_df$Site2 <- parts[2] }
  }
  all_results_list[[comp_name]] <- res_df
}

full_results <- bind_rows(all_results_list)

# Save full merged results
write.table(full_results, file.path(results_dir, "all_pairwise_results.tsv"),
            sep = "\t", quote = TRUE, row.names = FALSE)
cat(sprintf("Saved: all_pairwise_results.tsv (%d rows)\n", nrow(full_results)))

# Filter significant (only site_cleaned rows)
qcol <- "qval_individual"
sig_results <- full_results %>%
  filter(metadata == site_col, !is.na(.data[[qcol]]), .data[[qcol]] < alpha)

write.table(sig_results, file.path(results_dir, "significant_pairwise_results.tsv"),
            sep = "\t", quote = TRUE, row.names = FALSE)
cat(sprintf("Significant results (q<%.2f): %d rows\n", alpha, nrow(sig_results)))

# Summary counts per comparison
summary_counts <- sig_results %>%
  group_by(Site1, Site2, Comparison) %>%
  summarise(N_Significant = n_distinct(feature), .groups = "drop")

# Make symmetric
summary_symmetric <- summary_counts
for (i in seq_len(nrow(summary_counts))) {
  summary_symmetric <- bind_rows(summary_symmetric, data.frame(
    Site1 = summary_counts$Site2[i],
    Site2 = summary_counts$Site1[i],
    Comparison = paste0(summary_counts$Site2[i], "_vs_", summary_counts$Site1[i]),
    N_Significant = summary_counts$N_Significant[i],
    stringsAsFactors = FALSE
  ))
}

write.table(summary_symmetric, file.path(results_dir, "summary_diff_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n=== Significant taxa per comparison ===\n")
print(summary_counts %>% arrange(desc(N_Significant)))
cat("\n=== Merge done ===\n")
