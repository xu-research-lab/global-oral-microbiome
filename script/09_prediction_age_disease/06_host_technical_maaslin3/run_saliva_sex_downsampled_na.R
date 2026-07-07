#!/usr/bin/env Rscript
# ============================================================================
# Stratified downsampling of North America saliva samples, then MaAsLin3 for Sex
# within each Geography x Site_Cleaned stratum. Compares with other geographies.
#
# Design:
#   - Target North America sample size = 2500 (matches Europe ~2680, Asia ~2493)
#   - Stratify by BioProject x Sex x Age_Group x Region_custom
#   - Keep strata proportionate to original NA composition
#   - Run MaAsLin3: fixed = Age_Group + Sex + Region_custom, random = BioProject
#   - Plot heatmap comparing all geographies after downsampling
# ============================================================================

suppressPackageStartupMessages({
  library(maaslin3)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
})

set.seed(42)

BASE <- "/beegfs/huangxiaochang/global_oral/new_analysis/maaslin3"
setwd(BASE)
source("input_helpers.R")

OUT_ROOT <- file.path(BASE, "saliva_sex_downsampled_na")
RESULTS_DIR <- file.path(OUT_ROOT, "results")
PLOT_DIR <- file.path(OUT_ROOT, "plots")
TARGET_NA_N <- 2500
Q_THRESHOLD <- 0.1
CORES <- max(1, min(8, parallel::detectCores() - 1))

dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------------------
# Load data
# ----------------------------------------------------------------------------
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

# Clean columns
for (col in c("Age_Group", "Sex", "Region_custom", "BioProject", "Site_Cleaned", "Geography")) {
  if (col %in% colnames(meta)) {
    meta[[col]] <- trimws(as.character(meta[[col]]))
    meta[[col]][is.na(meta[[col]]) | meta[[col]] == "" | meta[[col]] == "NA"] <- NA_character_
  }
}

# Keep only saliva samples with complete metadata
valid <- meta$Site_Cleaned == "Saliva" & !is.na(meta$Site_Cleaned) &
  complete.cases(meta[, c("Age_Group", "Sex", "Region_custom", "BioProject", "Site_Cleaned", "Geography")])
meta <- meta[valid, , drop = FALSE]
abund <- abund[, rownames(meta), drop = FALSE]
cat(sprintf("Saliva samples with complete metadata: %d\n", nrow(meta)))

meta$Age_Group <- factor(meta$Age_Group)
meta$Sex <- factor(meta$Sex)
meta$Region_custom <- factor(meta$Region_custom)
meta$BioProject <- factor(meta$BioProject)
meta$Geography <- factor(meta$Geography)

# ----------------------------------------------------------------------------
# Stratified downsampling for North America
# ----------------------------------------------------------------------------
na_samples <- rownames(meta)[meta$Geography == "North America"]
na_meta <- meta[na_samples, , drop = FALSE]
na_abund <- abund[, na_samples, drop = FALSE]

cat(sprintf("\nOriginal North America saliva samples: %d\n", length(na_samples)))

# Create strata; collapse very small strata to avoid sampling failures
na_meta$stratum <- interaction(
  na_meta$BioProject,
  na_meta$Sex,
  na_meta$Age_Group,
  na_meta$Region_custom,
  drop = TRUE,
  sep = "__"
)

stratum_counts <- table(na_meta$stratum)
small_strata <- names(stratum_counts[stratum_counts < 3])
if (length(small_strata) > 0) {
  cat(sprintf("Collapsing %d very small strata (<3 samples)\n", length(small_strata)))
  for (s in small_strata) {
    parts <- strsplit(s, "__", fixed = TRUE)[[1]]
    na_meta$stratum[na_meta$stratum == s] <- paste(parts[1], parts[2], "__pooled_age__", parts[4], sep = "__")
  }
  na_meta$stratum <- factor(na_meta$stratum)
}

# Proportional allocation across strata, with minimum of 1 per non-empty stratum
stratum_n <- table(na_meta$stratum)
prop_n <- floor(TARGET_NA_N * stratum_n / sum(stratum_n))
prop_n[prop_n == 0] <- 1
# Adjust if sum exceeds target
if (sum(prop_n) > TARGET_NA_N) {
  extra <- sum(prop_n) - TARGET_NA_N
  sorted <- order(prop_n, decreasing = TRUE)
  for (i in seq_len(extra)) {
    k <- sorted[((i - 1) %% length(sorted)) + 1]
    if (prop_n[k] > 1) prop_n[k] <- prop_n[k] - 1
  }
}
cat(sprintf("Target strata sample allocation sum: %d\n", sum(prop_n)))

selected_in_na <- unlist(lapply(names(stratum_n), function(s) {
  rows <- which(na_meta$stratum == s)
  n_take <- min(prop_n[s], length(rows))
  sample(rows, n_take)
}))
na_selected_samples <- na_samples[selected_in_na]

cat(sprintf("Selected North America samples after stratified downsampling: %d\n", length(na_selected_samples)))

# Combine non-NA samples with downsampled NA samples
non_na_samples <- rownames(meta)[meta$Geography != "North America"]
keep_samples <- c(non_na_samples, na_selected_samples)
meta_ds <- meta[keep_samples, , drop = FALSE]
abund_ds <- abund[, keep_samples, drop = FALSE]

# Relevel factors
meta_ds$Age_Group <- droplevels(meta_ds$Age_Group)
meta_ds$Sex <- droplevels(meta_ds$Sex)
meta_ds$Region_custom <- droplevels(meta_ds$Region_custom)
meta_ds$BioProject <- droplevels(meta_ds$BioProject)
meta_ds$Geography <- droplevels(meta_ds$Geography)

if ("18-60 years" %in% levels(meta_ds$Age_Group)) {
  meta_ds$Age_Group <- relevel(meta_ds$Age_Group, ref = "18-60 years")
}
if ("F" %in% levels(meta_ds$Sex)) meta_ds$Sex <- relevel(meta_ds$Sex, ref = "F")
if ("V4" %in% levels(meta_ds$Region_custom)) meta_ds$Region_custom <- relevel(meta_ds$Region_custom, ref = "V4")

# Write downsampled metadata / abundance
write.table(meta_ds, file.path(OUT_ROOT, "downsampled_metadata.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(abund_ds, file.path(OUT_ROOT, "downsampled_genus_relabund.tsv"),
            sep = "\t", quote = FALSE, row.names = TRUE)

# Summary
summary_df <- meta_ds %>%
  group_by(Geography, BioProject, Sex, Age_Group, Region_custom) %>%
  summarise(n = n(), .groups = "drop")
write.csv(summary_df, file.path(OUT_ROOT, "downsampled_stratum_summary.csv"), row.names = FALSE)

geo_counts <- meta_ds %>% count(Geography)
write.csv(geo_counts, file.path(OUT_ROOT, "downsampled_geography_counts.csv"), row.names = FALSE)
cat("\nGeography counts after downsampling:\n")
print(geo_counts)

# ----------------------------------------------------------------------------
# Run MaAsLin3 per geography
# ----------------------------------------------------------------------------
geographies <- as.character(unique(meta_ds$Geography))
geographies <- geographies[geographies != "Africa"]  # exclude Africa as requested

all_results <- list()

for (geo in geographies) {
  cat(sprintf("\n[%s] Running MaAsLin3...\n", geo))
  geo_samples <- rownames(meta_ds)[meta_ds$Geography == geo]
  sub_meta <- meta_ds[geo_samples, , drop = FALSE]
  sub_abund <- abund_ds[, geo_samples, drop = FALSE]
  
  sub_meta$Age_Group <- droplevels(sub_meta$Age_Group)
  sub_meta$Sex <- droplevels(sub_meta$Sex)
  sub_meta$Region_custom <- droplevels(sub_meta$Region_custom)
  sub_meta$BioProject <- droplevels(sub_meta$BioProject)
  
  n_age <- nlevels(sub_meta$Age_Group)
  n_sex <- nlevels(sub_meta$Sex)
  n_region <- nlevels(sub_meta$Region_custom)
  n_bp <- nlevels(sub_meta$BioProject)
  cat(sprintf("  n=%d; Age=%d; Sex=%d; Region=%d; BioProject=%d\n",
              nrow(sub_meta), n_age, n_sex, n_region, n_bp))
  
  if (n_age < 2 || n_sex < 2 || n_region < 2 || n_bp < 2) {
    cat("  SKIP: insufficient factor levels\n")
    next
  }
  
  out_dir <- file.path(RESULTS_DIR, paste0("Saliva_", gsub(" ", "_", geo)))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  res_file <- file.path(out_dir, "all_results.tsv")
  if (!file.exists(res_file)) {
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
      cores = CORES,
      verbosity = "INFO"
    )
  }
  
  res <- read.delim(res_file, stringsAsFactors = FALSE)
  res$Geography <- geo
  res$Site <- "Saliva"
  all_results[[geo]] <- res
}

combined <- bind_rows(all_results, .id = "Geography")
write.table(combined, file.path(OUT_ROOT, "all_results_downsampled.tsv"),
            sep = "\t", quote = TRUE, row.names = FALSE)

# ----------------------------------------------------------------------------
# Summarize significant genera per geography
# ----------------------------------------------------------------------------
clean_genus <- function(x) {
  x <- trimws(as.character(x))
  has_genus <- grepl("g__", x)
  extracted <- regmatches(x, regexpr("g__[^;]*", x))
  x[has_genus] <- sub("^g__", "", extracted)
  x <- sub("_.+$", "", x)
  ifelse(is.na(x) | x == "" | x == "g__", NA_character_, x)
}

combined$Genus <- clean_genus(combined$feature)
sex_sig <- combined %>%
  filter(metadata == "Sex", !is.na(qval_individual), qval_individual < Q_THRESHOLD)

sig_counts <- sex_sig %>%
  group_by(Geography) %>%
  summarise(SignificantGenera = n_distinct(Genus), .groups = "drop")
write.csv(sig_counts, file.path(OUT_ROOT, "significant_genera_counts_downsampled.csv"), row.names = FALSE)
cat("\nSignificant genera per geography (downsampled):\n")
print(sig_counts)

# ----------------------------------------------------------------------------
# Heatmap
# ----------------------------------------------------------------------------
geography_order <- c("North America", "Asia", "Europe", "South America", "Oceania")
TAXONOMY <- "/beegfs/huangxiaochang/global_oral/taxonomy.tsv"

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
    slice(1) %>%
    ungroup()
  setNames(df$Family, df$Genus)
}

family_map <- build_family_map()

sex_data <- combined %>%
  filter(metadata == "Sex", !is.na(qval_individual)) %>%
  mutate(
    q_plot = ifelse(qval_individual <= 0, 1e-300, qval_individual),
    signed_logq = pmin(-log10(q_plot), 10) * sign(coef),
    sig_label = ifelse(qval_individual < Q_THRESHOLD, ifelse(coef > 0, "+", "-"), "")
  )

sig_rows <- sex_data %>% filter(qval_individual < Q_THRESHOLD)
kept_pairs <- sig_rows %>% distinct(value, Geography)
kept_genera <- sig_rows %>% distinct(value, Genus)
grid <- sex_data %>%
  inner_join(kept_pairs, by = c("value", "Geography")) %>%
  inner_join(kept_genera, by = c("value", "Genus"))

# Order genera
sig_first <- grid %>%
  filter(qval_individual < Q_THRESHOLD) %>%
  group_by(Genus) %>%
  summarise(
    min_q = min(qval_individual, na.rm = TRUE),
    first_geo = min(match(Geography, geography_order), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(first_geo, min_q)

all_genera <- unique(grid$Genus)
genus_order <- c(
  intersect(sig_first$Genus, all_genera),
  setdiff(all_genera, sig_first$Genus)
)

axis_labels <- sapply(genus_order, function(g) {
  fam <- family_map[g]
  if (is.na(fam) || fam == "f__") fam <- "f__"
  genus <- paste0("g__", g)
  fam_width <- max(nchar(family_map), na.rm = TRUE) + 2
  paste0(str_pad(fam, width = fam_width, side = "right"), genus)
})
names(axis_labels) <- genus_order

grid <- grid %>%
  mutate(
    sig_label = ifelse(is.na(sig_label), "", sig_label),
    Genus = factor(Genus, levels = rev(genus_order)),
    Geography = factor(Geography, levels = geography_order),
    value = factor(value)
  )

p <- ggplot(grid, aes(x = Geography, y = Genus, fill = signed_logq)) +
  geom_tile(color = "grey85", linewidth = 0.14) +
  geom_text(aes(label = sig_label), size = 2.2, fontface = "bold") +
  facet_wrap(~ value, nrow = 1, scales = "free_x") +
  scale_y_discrete(labels = axis_labels[genus_order]) +
  scale_x_discrete(drop = TRUE) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#D73027",
    midpoint = 0, limits = c(-10, 10), name = "signed\n-log10(q)"
  ) +
  labs(
    title = "Saliva Sex-Associated Microbiota (North America Downsampled)",
    subtitle = paste0("Stratified downsampling to n=", TARGET_NA_N,
                      "; q < ", Q_THRESHOLD, "; +/- indicates male coefficient direction"),
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7.5, face = "bold"),
    axis.text.y = element_text(size = 8.2, family = "mono", hjust = 0),
    strip.text = element_text(size = 8.2, face = "bold"),
    panel.grid = element_blank(),
    panel.spacing = unit(0.2, "lines"),
    legend.position = "bottom",
    legend.key.width = unit(0.75, "cm"),
    legend.key.height = unit(0.25, "cm"),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6.5),
    plot.title = element_text(hjust = 0.5, size = 11.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 8.0),
    plot.margin = margin(6, 6, 6, 6)
  )

width <- max(6.2, length(unique(grid$value)) * 2.5)
height <- max(5.2, length(genus_order) * 0.18)
ggsave(file.path(PLOT_DIR, "saliva_Sex_heatmap_downsampled_na.pdf"),
       p, width = width, height = height, limitsize = FALSE)
ggsave(file.path(PLOT_DIR, "saliva_Sex_heatmap_downsampled_na.png"),
       p, width = width, height = height, dpi = 300, limitsize = FALSE)

cat(sprintf("\nSaved heatmap to %s (%d genera, %.1f x %.1f in)\n",
            PLOT_DIR, length(genus_order), width, height))
cat("Done.\n")
