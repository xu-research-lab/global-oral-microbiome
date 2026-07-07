#!/usr/bin/env Rscript
# ============================================================================
# MaAsLin3 Stratified Analysis (adapted for new data)
# Runs MaAsLin3 for each Site_Cleaned x Geography combination
# Fixed effects: Age_Group, Sex, Region_custom
# Random effect: BioProject
# ============================================================================

suppressPackageStartupMessages({
  library(maaslin3)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

setwd("/beegfs/huangxiaochang/global_oral/new_analysis/maaslin3")
source("input_helpers.R")

cat("=== Loading data ===\n")
genus_data <- load_genus_relabund()
metadata <- load_standard_metadata(dedup = TRUE)

sample_ids <- trimws(as.character(metadata$SampleID))
keep_ids <- !is.na(sample_ids) & sample_ids != "" & !duplicated(sample_ids)
metadata <- metadata[keep_ids, , drop = FALSE]
sample_ids <- sample_ids[keep_ids]

common_samples <- intersect(sample_ids, colnames(genus_data))
cat(sprintf("Common samples: %d\n", length(common_samples)))

abund <- genus_data[, common_samples, drop = FALSE]
meta <- metadata[match(common_samples, sample_ids), , drop = FALSE]
rownames(meta) <- common_samples

for (col in c("Age_Group", "Sex", "Region_custom", "BioProject", "Site_Cleaned", "Geography")) {
  if (col %in% colnames(meta)) meta[[col]] <- trimws(meta[[col]])
}

valid <- complete.cases(meta[, c("Age_Group", "Sex", "Region_custom", "BioProject", "Site_Cleaned", "Geography")]) &
  meta$Age_Group != "" & meta$Sex != "" & meta$Region_custom != "" & meta$BioProject != "" &
  meta$Age_Group != "NA" & meta$Sex != "NA"
meta <- meta[valid, , drop = FALSE]
abund <- abund[, rownames(meta), drop = FALSE]
cat(sprintf("Samples with complete metadata: %d\n", nrow(meta)))

meta$Age_Group <- factor(meta$Age_Group)
meta$Sex <- factor(meta$Sex)
meta$BioProject <- factor(meta$BioProject)

group_table <- meta %>%
  group_by(Site_Cleaned, Geography) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n))

cat("\n=== Group sizes ===\n")
print(group_table)

MIN_N <- 50
group_table <- group_table %>% filter(n >= MIN_N)
cat(sprintf("\nGroups with >= %d samples: %d\n", MIN_N, nrow(group_table)))
write.csv(group_table, "stratified_group_sizes.csv", row.names = FALSE)

results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

all_results <- list()
successful_groups <- character()

for (i in seq_len(nrow(group_table))) {
  site <- group_table$Site_Cleaned[i]
  geo <- group_table$Geography[i]
  n_samp <- group_table$n[i]

  group_key <- paste0(gsub(" ", "_", site), "_", gsub(" ", "_", geo))
  out_dir <- file.path(results_dir, group_key)

  cat(sprintf("\n[%d/%d] %s (n=%d)\n", i, nrow(group_table), group_key, n_samp))

  res_file_existing <- file.path(out_dir, "all_results.tsv")
  if (file.exists(res_file_existing)) {
    cat("  SKIP: existing completed result found\n")
    res <- read.delim(res_file_existing, stringsAsFactors = FALSE)
    res$Site_Cleaned <- site
    res$Geography <- geo
    all_results[[group_key]] <- res
    successful_groups <- c(successful_groups, group_key)
    next
  }

  grp_idx <- meta$Site_Cleaned == site & meta$Geography == geo
  sub_meta <- meta[grp_idx, , drop = FALSE]
  sub_abund <- abund[, rownames(sub_meta), drop = FALSE]

  sub_meta$Age_Group <- droplevels(sub_meta$Age_Group)
  sub_meta$Sex <- droplevels(sub_meta$Sex)
  sub_meta$Region_custom <- droplevels(factor(sub_meta$Region_custom))
  sub_meta$BioProject <- droplevels(sub_meta$BioProject)

  if ("18-60 years" %in% levels(sub_meta$Age_Group)) {
    sub_meta$Age_Group <- relevel(sub_meta$Age_Group, ref = "18-60 years")
  }
  if ("F" %in% levels(sub_meta$Sex)) sub_meta$Sex <- relevel(sub_meta$Sex, ref = "F")
  if ("V4" %in% levels(sub_meta$Region_custom)) sub_meta$Region_custom <- relevel(sub_meta$Region_custom, ref = "V4")

  n_age <- length(unique(sub_meta$Age_Group))
  n_sex <- length(unique(sub_meta$Sex))
  n_region <- length(unique(sub_meta$Region_custom))
  n_bp <- length(unique(sub_meta$BioProject))

  cat(sprintf("  Levels: Age=%d, Sex=%d, Region=%d, BioProject=%d\n", n_age, n_sex, n_region, n_bp))

  if (n_age < 2 || n_sex < 2 || n_region < 2 || n_bp < 2) {
    cat("  SKIP: insufficient factor levels\n")
    next
  }

  result <- tryCatch({
    maaslin3(
      input_data = t(sub_abund),
      input_metadata = sub_meta,
      output = out_dir,
      fixed_effects = c("Age_Group", "Sex", "Region_custom"),
      random_effects = c("BioProject"),
      min_abundance = 0,
      min_prevalence = 0.1,
      max_significance = 0.25,
      normalization = "TSS",
      transform = "LOG",
      correction = "BH",
      warn_prevalence = FALSE,
      cores = max(1, min(8, parallel::detectCores() - 1)),
      verbosity = "INFO"
    )

    res_file <- file.path(out_dir, "all_results.tsv")
    if (file.exists(res_file)) {
      res <- read.delim(res_file, stringsAsFactors = FALSE)
      res$Site_Cleaned <- site
      res$Geography <- geo
      n_sig <- sum(res$qval_individual < 0.1, na.rm = TRUE)
      cat(sprintf("  Done: %d associations, %d significant (q<0.1)\n", nrow(res), n_sig))
      res
    } else {
      cat("  WARN: no results file\n")
      NULL
    }
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
    NULL
  })

  if (!is.null(result)) {
    all_results[[group_key]] <- result
    successful_groups <- c(successful_groups, group_key)
  }
  gc(verbose = FALSE)
}

cat(sprintf("\n=== Completed %d/%d groups ===\n", length(all_results), nrow(group_table)))

saveRDS(all_results, file = "stratified_results.rds")
saveRDS(successful_groups, file = "successful_groups.rds")
if (length(all_results) > 0) {
  combined <- bind_rows(all_results, .id = "Group")
  write.table(combined, "stratified_results_all.tsv",
              sep = "\t", quote = TRUE, row.names = FALSE)
}

extract_genus <- function(x) {
  g <- regmatches(x, regexpr("g__[^;]*", x))
  if (length(g) == 0 || identical(g, "g__")) return(x)
  gsub("g__", "", g)
}

build_hm <- function(results_list, variable_name, sig_threshold = 0.1) {
  all_rows <- list()
  for (gk in names(results_list)) {
    res <- results_list[[gk]]
    sig <- res[res$metadata == variable_name & res$qval_individual < sig_threshold & !is.na(res$qval_individual), ]
    if (nrow(sig) == 0) next
    sig$Genus <- sapply(sig$feature, extract_genus)
    all_rows[[gk]] <- sig
  }
  if (length(all_rows) == 0) return(data.frame())
  bind_rows(all_rows)
}

age_hm <- build_hm(all_results, "Age_Group")
sex_hm <- build_hm(all_results, "Sex")
region_hm <- build_hm(all_results, "Region_custom")

cat(sprintf("Age significant: %d\n", nrow(age_hm)))
cat(sprintf("Sex significant: %d\n", nrow(sex_hm)))
cat(sprintf("Region significant: %d\n", nrow(region_hm)))

saveRDS(list(age = age_hm, sex = sex_hm, region = region_hm), file = "heatmap_data.rds")

plot_heatmap <- function(hm_data, x_var, x_label, title, nrow_facets = 1) {
  if (nrow(hm_data) == 0) {
    cat(sprintf("  No data for %s\n", title))
    return(invisible(NULL))
  }

  hm_data$sig_label <- ifelse(hm_data$coef > 0, "+", "-")
  hm_data$facet_label <- paste0(hm_data$Site_Cleaned, " | ", hm_data$Geography)

  genus_ord <- hm_data %>%
    group_by(Genus) %>%
    summarise(min_q = min(qval_individual, na.rm = TRUE),
              mean_coef = mean(coef, na.rm = TRUE), .groups = "drop") %>%
    arrange(min_q, -abs(mean_coef))

  hm_data$Genus <- factor(hm_data$Genus, levels = genus_ord$Genus)

  p <- ggplot(hm_data, aes_string(x = x_var, y = "Genus", fill = "sig_label")) +
    geom_tile(color = "grey80", linewidth = 0.2) +
    geom_text(aes(label = sig_label), size = 3, fontface = "bold") +
    scale_fill_manual(values = c("+" = "#d73027", "-" = "#4575b4"), name = "Direction") +
    facet_wrap(~ facet_label, scales = "free_x", nrow = nrow_facets) +
    labs(title = title, x = x_label, y = "Genus") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 7),
      strip.background = element_rect(fill = "grey90", colour = "grey70"),
      strip.text = element_text(size = 7, face = "bold"),
      legend.position = "bottom",
      panel.spacing = unit(0.3, "lines"),
      panel.grid = element_blank()
    )

  h <- max(6, length(genus_ord$Genus) * 0.25)
  n_f <- length(unique(hm_data$facet_label))
  w <- max(10, n_f * 2)
  fn <- gsub(" ", "_", tolower(title))
  ggsave(paste0(fn, ".pdf"), p, width = w, height = h, limitsize = FALSE)
  ggsave(paste0(fn, ".png"), p, width = w, height = h, dpi = 300, limitsize = FALSE)
  cat(sprintf("  Saved: %s.pdf (%d rows, %d facets)\n", fn, length(genus_ord$Genus), n_f))
}

cat("\n=== Generating heatmaps ===\n")
if (nrow(age_hm) > 0) {
  age_hm$value <- factor(age_hm$value)
  plot_heatmap(age_hm, "value", "Age Group", "Age-Associated Microbiota (q < 0.1)")
}
if (nrow(sex_hm) > 0) {
  sex_hm$value <- factor(sex_hm$value)
  plot_heatmap(sex_hm, "value", "Sex Contrast", "Sex-Associated Microbiota (q < 0.1)")
}
if (nrow(region_hm) > 0) {
  reg_order <- c("V1-V2", "V1-V3", "V2", "V3", "V3-V4", "V3-V5", "V4", "V4-V5", "V5-V7", "V6-V8")
  reg_order <- reg_order[reg_order %in% unique(region_hm$value)]
  region_hm$value <- factor(region_hm$value, levels = reg_order)
  plot_heatmap(region_hm, "value", "16S Region", "Primer-Region-Associated Microbiota (q < 0.1)")
}

cat("\n=== Done ===\n")

