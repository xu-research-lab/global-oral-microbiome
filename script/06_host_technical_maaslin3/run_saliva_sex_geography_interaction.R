#!/usr/bin/env Rscript
# ============================================================================
# MaAsLin3 Saliva Sex x Geography Interaction Model
# Tests whether male-vs-female microbial associations differ by geography.
#
# Model: ~ Sex * Geography + Age_Group + Region_custom + (1 | BioProject)
# Reference levels: Sex = F, Geography = North America when available,
#                   Age_Group = 18-60 years when available,
#                   Region_custom = V4 when available.
# ============================================================================

suppressPackageStartupMessages({
  library(maaslin3)
  library(dplyr)
  library(tidyr)
})

setwd("/beegfs/huangxiaochang/global_oral/new_analysis/maaslin3")
source("input_helpers.R")

MAASLIN3_CORES <- 4
Q_THRESHOLDS <- c(0.05, 0.1)
MODEL_FORMULA <- "~ Sex * Geography + Age_Group + Region_custom + (1 | BioProject)"

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}

EVALUATE_ONLY <- get_arg("evaluate_only", "abundance")
if (!EVALUATE_ONLY %in% c("abundance", "prevalence", "all")) {
  stop("--evaluate_only must be one of: abundance, prevalence, all")
}
EVALUATE_ONLY_ARG <- if (EVALUATE_ONLY == "all") NULL else EVALUATE_ONLY

OUT_ROOT <- get_arg(
  "out",
  if (EVALUATE_ONLY == "all") {
    "saliva_sex_geography_interaction_model"
  } else {
    paste0("saliva_sex_geography_interaction_", EVALUATE_ONLY)
  }
)

dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)

clean_missing <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == "" | x %in% c("NA", "nan", "None", "Unknown")] <- NA_character_
  x
}

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

needed_cols <- c("Age_Group", "Sex", "Region_custom", "BioProject", "Site_Cleaned", "Geography")
missing_cols <- setdiff(needed_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing required metadata columns: %s", paste(missing_cols, collapse = ", ")))
}

for (col in needed_cols) {
  meta[[col]] <- clean_missing(meta[[col]])
}

valid <- meta$Site_Cleaned == "Saliva" & complete.cases(meta[, needed_cols])
meta <- meta[valid, , drop = FALSE]
abund <- abund[, rownames(meta), drop = FALSE]

cat(sprintf("Saliva samples with complete metadata: %d\n", nrow(meta)))
if (nrow(meta) < 50) stop("Too few complete Saliva samples for interaction model.")

meta$Age_Group <- droplevels(factor(meta$Age_Group))
meta$Sex <- droplevels(factor(meta$Sex))
meta$Region_custom <- droplevels(factor(meta$Region_custom))
meta$Geography <- droplevels(factor(meta$Geography))
meta$BioProject <- droplevels(factor(meta$BioProject))

if ("18-60 years" %in% levels(meta$Age_Group)) {
  meta$Age_Group <- relevel(meta$Age_Group, ref = "18-60 years")
}
if ("F" %in% levels(meta$Sex)) {
  meta$Sex <- relevel(meta$Sex, ref = "F")
}
if ("V4" %in% levels(meta$Region_custom)) {
  meta$Region_custom <- relevel(meta$Region_custom, ref = "V4")
}
if ("North America" %in% levels(meta$Geography)) {
  meta$Geography <- relevel(meta$Geography, ref = "North America")
}

factor_counts <- bind_rows(lapply(needed_cols, function(col) {
  tab <- table(meta[[col]])
  data.frame(
    Variable = col,
    Level = names(tab),
    N = as.integer(tab),
    stringsAsFactors = FALSE
  )
}))
write.csv(factor_counts, file.path(OUT_ROOT, "saliva_interaction_factor_counts.csv"), row.names = FALSE)

sex_by_geography <- as.data.frame.matrix(table(meta$Geography, meta$Sex))
sex_by_geography <- tibble::rownames_to_column(sex_by_geography, "Geography")
write.csv(sex_by_geography, file.path(OUT_ROOT, "saliva_sex_by_geography_counts.csv"), row.names = FALSE)

cat("\n=== Model factors ===\n")
cat(sprintf("Sex levels: %s\n", paste(levels(meta$Sex), collapse = ", ")))
cat(sprintf("Geography levels: %s\n", paste(levels(meta$Geography), collapse = ", ")))
cat(sprintf("Age_Group levels: %s\n", paste(levels(meta$Age_Group), collapse = ", ")))
cat(sprintf("Region_custom levels: %s\n", paste(levels(meta$Region_custom), collapse = ", ")))
cat(sprintf("BioProject levels: %d\n", nlevels(meta$BioProject)))
cat(sprintf("Formula: %s\n", MODEL_FORMULA))

if (nlevels(meta$Sex) < 2) stop("Sex has fewer than 2 levels after filtering.")
if (nlevels(meta$Geography) < 2) stop("Geography has fewer than 2 levels after filtering.")
if (nlevels(meta$Age_Group) < 2) stop("Age_Group has fewer than 2 levels after filtering.")
if (nlevels(meta$Region_custom) < 2) stop("Region_custom has fewer than 2 levels after filtering.")
if (nlevels(meta$BioProject) < 2) stop("BioProject has fewer than 2 levels after filtering.")

write.csv(data.frame(
  Samples = nrow(meta),
  Features = nrow(abund),
  BioProjects = nlevels(meta$BioProject),
  Formula = MODEL_FORMULA,
  Reference_Sex = levels(meta$Sex)[1],
  Reference_Geography = levels(meta$Geography)[1],
  Reference_Age_Group = levels(meta$Age_Group)[1],
  Reference_Region_custom = levels(meta$Region_custom)[1],
  Evaluate_Only = EVALUATE_ONLY,
  Cores = MAASLIN3_CORES,
  stringsAsFactors = FALSE
), file.path(OUT_ROOT, "saliva_interaction_run_summary.csv"), row.names = FALSE)

res_file <- file.path(OUT_ROOT, "all_results.tsv")
if (!file.exists(res_file)) {
  cat("\n=== Running MaAsLin3 ===\n")
  maaslin3(
    input_data = t(abund),
    input_metadata = meta,
    output = OUT_ROOT,
    formula = MODEL_FORMULA,
    min_abundance = 0,
    min_prevalence = 0.1,
    max_significance = 0.25,
    normalization = "TSS",
    transform = "LOG",
    correction = "BH",
    warn_prevalence = FALSE,
    evaluate_only = EVALUATE_ONLY_ARG,
    cores = MAASLIN3_CORES,
    verbosity = "INFO"
  )
} else {
  cat("\nExisting all_results.tsv found; skipping MaAsLin3 fit.\n")
}

cat("\n=== Summarizing results ===\n")
results <- read.delim(res_file, stringsAsFactors = FALSE, check.names = FALSE)

is_sex_main <- results$metadata == "Sex" | grepl("^Sex", results$name)
is_interaction <- grepl(":", results$name) | grepl(":", results$metadata)
is_sex_geo_interaction <- is_interaction &
  grepl("Sex", results$name) & grepl("Geography", results$name)

sex_main <- results %>%
  filter(is_sex_main, !is_interaction) %>%
  arrange(qval_individual, desc(abs(coef)))

sex_geo_interactions <- results %>%
  filter(is_sex_geo_interaction) %>%
  arrange(qval_individual, desc(abs(coef)))

write.table(sex_main, file.path(OUT_ROOT, "sex_main_effect_results_all.tsv"),
            sep = "\t", quote = TRUE, row.names = FALSE)
write.table(sex_geo_interactions, file.path(OUT_ROOT, "sex_geography_interaction_results_all.tsv"),
            sep = "\t", quote = TRUE, row.names = FALSE)

for (q in Q_THRESHOLDS) {
  q_label <- gsub("\\.", "p", as.character(q))
  write.table(sex_main %>% filter(!is.na(qval_individual), qval_individual <= q),
              file.path(OUT_ROOT, paste0("sex_main_effect_significant_q", q_label, ".tsv")),
              sep = "\t", quote = TRUE, row.names = FALSE)
  write.table(sex_geo_interactions %>% filter(!is.na(qval_individual), qval_individual <= q),
              file.path(OUT_ROOT, paste0("sex_geography_interaction_significant_q", q_label, ".tsv")),
              sep = "\t", quote = TRUE, row.names = FALSE)
}

extract_interaction_geography <- function(term_name, geography_levels) {
  non_ref <- geography_levels[geography_levels != levels(meta$Geography)[1]]
  hits <- non_ref[vapply(non_ref, function(g) grepl(g, term_name, fixed = TRUE), logical(1))]
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

sex_main_simple <- sex_main %>%
  filter(value == "M" | name == "SexM") %>%
  select(feature, model, sex_coef_ref = coef, sex_qval_ref = qval_individual)

interaction_simple <- sex_geo_interactions %>%
  mutate(Geography = vapply(name, extract_interaction_geography, character(1),
                            geography_levels = levels(meta$Geography))) %>%
  filter(!is.na(Geography)) %>%
  select(feature, model, Geography, interaction_coef = coef,
         interaction_qval = qval_individual, interaction_name = name)

ref_geo <- levels(meta$Geography)[1]
ref_effects <- sex_main_simple %>%
  mutate(Geography = ref_geo,
         estimated_sex_coef_in_geography = sex_coef_ref,
         interaction_coef = 0,
         interaction_qval = NA_real_,
         interaction_name = NA_character_)

nonref_effects <- interaction_simple %>%
  left_join(sex_main_simple, by = c("feature", "model")) %>%
  mutate(estimated_sex_coef_in_geography = sex_coef_ref + interaction_coef)

sex_effect_by_geography <- bind_rows(ref_effects, nonref_effects) %>%
  arrange(feature, model, Geography)

write.table(sex_effect_by_geography,
            file.path(OUT_ROOT, "estimated_sex_effect_by_geography.tsv"),
            sep = "\t", quote = TRUE, row.names = FALSE)

sig_summary <- bind_rows(lapply(Q_THRESHOLDS, function(q) {
  data.frame(
    q_threshold = q,
    significant_sex_main_associations = sum(!is.na(sex_main$qval_individual) & sex_main$qval_individual <= q),
    significant_sex_main_features = length(unique(sex_main$feature[!is.na(sex_main$qval_individual) & sex_main$qval_individual <= q])),
    significant_interaction_associations = sum(!is.na(sex_geo_interactions$qval_individual) & sex_geo_interactions$qval_individual <= q),
    significant_interaction_features = length(unique(sex_geo_interactions$feature[!is.na(sex_geo_interactions$qval_individual) & sex_geo_interactions$qval_individual <= q])),
    stringsAsFactors = FALSE
  )
}))
write.csv(sig_summary, file.path(OUT_ROOT, "saliva_interaction_significant_summary.csv"), row.names = FALSE)

cat("\n=== Significant summary ===\n")
print(sig_summary)
cat("\nOutputs written to: ", normalizePath(OUT_ROOT), "\n", sep = "")
cat("\n=== Done ===\n")
