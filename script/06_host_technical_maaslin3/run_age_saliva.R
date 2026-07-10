#!/usr/bin/env Rscript
# ============================================================================
# MaAsLin3 v2 — Analysis 1: Age-associated microbiota in Saliva
# Sample scope: Saliva only
# Model: age_cleaned_group + sex_cleaned + sequencing_region_cleaned + continent_group
# Random effect: bioproject_id
# ============================================================================

suppressPackageStartupMessages({
  library(maaslin3)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
})

setwd("/beegfs/huangxiaochang/global_oral/new_analysis/maaslin3_v2")

OUT_DIR <- "age_saliva_results"
MAASLIN3_CORES <- 1
Q_THRESHOLD <- 0.1
TAXONOMY <- "/beegfs/huangxiaochang/global_oral/taxonomy.tsv"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 0. Helper functions
# ============================================================================
clean_genus <- function(x) {
  x <- trimws(as.character(x))
  has_genus <- grepl("g__", x)
  extracted <- regmatches(x, regexpr("g__[^;]*", x))
  x[has_genus] <- sub("^g__", "", extracted)
  x <- sub("_.+$", "", x)
  ifelse(is.na(x) | x == "" | x == "g__", NA_character_, x)
}
extract_rank <- function(taxon, prefix) {
  parts <- strsplit(as.character(taxon), ";", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  hit <- parts[startsWith(parts, prefix)]
  if (length(hit) == 0 || hit[1] == prefix) return(NA_character_)
  hit[1]
}
build_family_map <- function() {
  tax <- read.delim(TAXONOMY, stringsAsFactors = FALSE, check.names = FALSE)
  tax_col <- if ("Taxon" %in% colnames(tax)) "Taxon" else colnames(tax)[2]
  fam <- vapply(tax[[tax_col]], extract_rank, character(1), prefix = "f__")
  gen <- vapply(tax[[tax_col]], extract_rank, character(1), prefix = "g__")
  gen_clean <- clean_genus(gen)
  df <- data.frame(Genus = gen_clean, Family = fam, stringsAsFactors = FALSE) %>%
    filter(!is.na(Genus), !is.na(Family), Family != "f__") %>%
    count(Genus, Family, name = "n") %>%
    group_by(Genus) %>%
    arrange(desc(n), Family, .by_group = TRUE) %>%
    slice(1) %>% ungroup()
  setNames(df$Family, df$Genus)
}
family_map <- build_family_map()

# ============================================================================
# 1. Load data
# ============================================================================
cat("=== Loading data ===\n")
genus_data <- read.table("genus_table_relabund.tsv", header = TRUE, sep = "\t",
                         row.names = 1, check.names = FALSE,
                         comment.char = "", stringsAsFactors = FALSE)
metadata <- read.table("../metadata_oral_qc_dedup.tsv", header = TRUE, sep = "\t",
                       fill = TRUE, comment.char = "", stringsAsFactors = FALSE, quote = "")

# Column mappings
site_map <- c("saliva" = "Saliva", "tongue" = "Tongue", "tooth" = "Tooth",
              "oral swab" = "Oral_swab", "oral mucosa" = "Oral_Mucosa",
              "gingiva" = "Gingiva", "other" = "Others")
metadata$site_cleaned <- ifelse(metadata$site_cleaned %in% names(site_map),
                                site_map[metadata$site_cleaned], metadata$site_cleaned)
metadata$status_cleaned[metadata$status_cleaned %in% c("HC", "Healthy")] <- "Healthy Control"
metadata$age_cleaned_group <- gsub("＞", ">", metadata$age_cleaned_group)

# ============================================================================
# 2. Match samples — Saliva only
# ============================================================================
sample_ids <- trimws(as.character(metadata$sample_id))
keep_ids <- !is.na(sample_ids) & sample_ids != "" & !duplicated(sample_ids)
metadata <- metadata[keep_ids, , drop = FALSE]
sample_ids <- sample_ids[keep_ids]

common_samples <- intersect(sample_ids, colnames(genus_data))
cat(sprintf("Common samples overall: %d\n", length(common_samples)))

abund <- genus_data[, common_samples, drop = FALSE]
meta <- metadata[match(common_samples, sample_ids), , drop = FALSE]
rownames(meta) <- common_samples

# Filter to Saliva
needed_cols <- c("age_cleaned_group", "sex_cleaned", "sequencing_region_cleaned",
                 "bioproject_id", "continent_group", "site_cleaned")
for (col in needed_cols) {
  meta[[col]] <- trimws(as.character(meta[[col]]))
  meta[[col]][is.na(meta[[col]]) | meta[[col]] == "" | meta[[col]] == "NA"] <- NA_character_
}

meta <- meta[meta$site_cleaned == "Saliva" & !is.na(meta$site_cleaned), , drop = FALSE]
abund <- abund[, rownames(meta), drop = FALSE]
cat(sprintf("Saliva samples before complete-case filter: %d\n", nrow(meta)))

valid <- complete.cases(meta[, needed_cols]) &
  meta$age_cleaned_group != "" & meta$sex_cleaned != "" &
  meta$sequencing_region_cleaned != "" & meta$bioproject_id != ""
meta <- meta[valid, , drop = FALSE]
abund <- abund[, rownames(meta), drop = FALSE]
cat(sprintf("Saliva samples with complete metadata: %d\n", nrow(meta)))

# Factorize
meta$age_cleaned_group <- factor(meta$age_cleaned_group)
meta$sex_cleaned <- factor(meta$sex_cleaned)
meta$sequencing_region_cleaned <- factor(meta$sequencing_region_cleaned)
meta$continent_group <- factor(meta$continent_group)
meta$bioproject_id <- factor(meta$bioproject_id)

# Reference levels
if ("18-60 years" %in% levels(meta$age_cleaned_group))
  meta$age_cleaned_group <- relevel(meta$age_cleaned_group, ref = "18-60 years")
if ("F" %in% levels(meta$sex_cleaned))
  meta$sex_cleaned <- relevel(meta$sex_cleaned, ref = "F")
if ("V4" %in% levels(meta$sequencing_region_cleaned))
  meta$sequencing_region_cleaned <- relevel(meta$sequencing_region_cleaned, ref = "V4")

cat(sprintf("Age levels (%d): %s\n", nlevels(meta$age_cleaned_group),
            paste(levels(meta$age_cleaned_group), collapse=", ")))
cat(sprintf("Sex levels (%d): %s\n", nlevels(meta$sex_cleaned),
            paste(levels(meta$sex_cleaned), collapse=", ")))
cat(sprintf("Region levels (%d): %s\n", nlevels(meta$sequencing_region_cleaned),
            paste(levels(meta$sequencing_region_cleaned), collapse=", ")))
cat(sprintf("BioProject levels: %d\n", nlevels(meta$bioproject_id)))

# ============================================================================
# 3. Run MaAsLin3
# ============================================================================
res_file <- file.path(OUT_DIR, "all_results.tsv")
if (!file.exists(res_file)) {
  cat("\n=== Running MaAsLin3 for Saliva age analysis ===\n")
  maaslin3(
    input_data = t(abund),
    input_metadata = meta,
    output = OUT_DIR,
    fixed_effects = c("age_cleaned_group", "sex_cleaned",
                      "sequencing_region_cleaned", "continent_group"),
    random_effects = c("bioproject_id"),
    evaluate_only = "abundance",
    min_abundance = 0, min_prevalence = 0.1, max_significance = 0.25,
    normalization = "TSS", transform = "LOG", correction = "BH",
    warn_prevalence = FALSE, cores = MAASLIN3_CORES, verbosity = "INFO"
  )
} else {
  cat("Existing results found; skipping fit.\n")
}

# ============================================================================
# 4. Extract age results
# ============================================================================
results <- read.delim(res_file, stringsAsFactors = FALSE)
age_results <- results %>%
  filter(metadata == "age_cleaned_group") %>%
  mutate(Genus = clean_genus(feature)) %>%
  filter(!is.na(Genus)) %>%
  arrange(qval_individual, desc(abs(coef)))

write.table(age_results, file.path(OUT_DIR, "age_results_all.tsv"),
            sep = "\t", quote = TRUE, row.names = FALSE)

age_sig <- age_results %>% filter(!is.na(qval_individual), qval_individual < Q_THRESHOLD)
write.table(age_sig, file.path(OUT_DIR, "age_results_significant_q0.1.tsv"),
            sep = "\t", quote = TRUE, row.names = FALSE)

cat(sprintf("\nTotal age associations: %d\n", nrow(age_results)))
cat(sprintf("Significant (q<%.2f): %d (%d unique genera)\n",
            Q_THRESHOLD, nrow(age_sig), n_distinct(age_sig$Genus)))

# ============================================================================
# 5. Age heatmap
# ============================================================================
if (nrow(age_sig) > 0) {
  hm_data <- age_sig %>%
    mutate(
      sig_label = ifelse(coef > 0, "+", "-"),
      q_plot = ifelse(qval_individual <= 0, 1e-300, qval_individual),
      signed_logq = pmin(-log10(q_plot), 10) * sign(coef)
    )

  genus_order <- hm_data %>%
    group_by(Genus) %>%
    summarise(min_q = min(qval_individual), .groups = "drop") %>%
    arrange(min_q) %>% pull(Genus)

  age_values <- unique(hm_data$value)
  hm_data <- hm_data %>%
    mutate(
      Genus = factor(Genus, levels = rev(genus_order)),
      value = factor(value, levels = age_values)
    )

  # Axis labels with family prefix
  axis_labels <- sapply(genus_order, function(g) {
    fam <- family_map[g]; if (is.na(fam)) fam <- "f__"
    paste0(str_pad(fam, width = max(nchar(family_map), na.rm=TRUE)+2, side="right"),
           "g__", g)
  })
  names(axis_labels) <- genus_order

  w <- max(8, length(age_values) * 1.2)
  h <- max(6, length(genus_order) * 0.25)

  p <- ggplot(hm_data, aes(x = value, y = Genus, fill = signed_logq)) +
    geom_tile(color = "grey85", linewidth = 0.15) +
    geom_text(aes(label = sig_label), size = 3, fontface = "bold") +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#D73027",
                         midpoint = 0, limits = c(-10, 10), name = "signed\n-log10(q)") +
    scale_y_discrete(labels = axis_labels) +
    labs(title = "Saliva Age-Associated Microbiota",
         subtitle = paste0("Ref: 18-60 years; q < ", Q_THRESHOLD,
                          "; +/- indicates direction vs reference"),
         x = "Age Group", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          axis.text.y = element_text(size = 7, family = "mono"),
          panel.grid = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(file.path(OUT_DIR, "age_associated_microbiota.pdf"), p,
         width = w, height = h, limitsize = FALSE)
  ggsave(file.path(OUT_DIR, "age_associated_microbiota.png"), p,
         width = w, height = h, limitsize = FALSE, dpi = 300)
  cat(sprintf("Saved age heatmap: %d genera x %d age levels (%.1f x %.1f in)\n",
              length(genus_order), length(age_values), w, h))
}

cat("=== Age saliva analysis done ===\n")
