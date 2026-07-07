#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(maaslin3)
  library(dplyr)
  library(tidyr)
})

setwd("/beegfs/huangxiaochang/global_oral/new_analysis/maaslin3")
source("input_helpers.R")

out_dir <- "disease_results/Dental_Caries_PRJNA383868"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Loading data ===\n")
genus_data <- load_genus_relabund()
metadata <- load_standard_metadata(dedup = TRUE)

sample_ids <- trimws(as.character(metadata$SampleID))
keep_ids <- !is.na(sample_ids) & sample_ids != "" & !duplicated(sample_ids)
metadata <- metadata[keep_ids, , drop = FALSE]
sample_ids <- sample_ids[keep_ids]

common_samples <- intersect(sample_ids, colnames(genus_data))
meta <- metadata[match(common_samples, sample_ids), , drop = FALSE]
rownames(meta) <- common_samples
abund <- genus_data[, common_samples, drop = FALSE]

for (col in c("Age_Group", "Sex", "Region_custom", "BioProject", "Site_Cleaned", "Geography", "Status_Cleaned")) {
  meta[[col]] <- trimws(as.character(meta[[col]]))
}

keep <- meta$BioProject == "PRJNA383868" & meta$Status_Cleaned %in% c("Healthy Control", "Dental Caries")
sub_meta <- meta[keep, , drop = FALSE]
sub_abund <- abund[, rownames(sub_meta), drop = FALSE]

sub_meta$Disease_Status <- factor(
  ifelse(sub_meta$Status_Cleaned == "Dental Caries", "Dental Caries", "Healthy Control"),
  levels = c("Healthy Control", "Dental Caries")
)

optional_effects <- c("Age_Group", "Sex", "Site_Cleaned", "Region_custom", "Geography")
for (eff in optional_effects) {
  sub_meta[[eff]][is.na(sub_meta[[eff]]) | sub_meta[[eff]] == ""] <- "Unknown"
  sub_meta[[eff]] <- droplevels(factor(sub_meta[[eff]]))
}

fixed_effects <- c("Disease_Status")
for (eff in optional_effects) {
  if (nlevels(sub_meta[[eff]]) >= 2) fixed_effects <- c(fixed_effects, eff)
}

run_info <- data.frame(
  BioProject = "PRJNA383868",
  Disease = "Dental Caries",
  Samples = nrow(sub_meta),
  Cases = sum(sub_meta$Disease_Status == "Dental Caries"),
  Controls = sum(sub_meta$Disease_Status == "Healthy Control"),
  Site = paste(names(table(sub_meta$Site_Cleaned)), table(sub_meta$Site_Cleaned), sep = "=", collapse = "; "),
  Age_Group = paste(names(table(sub_meta$Age_Group)), table(sub_meta$Age_Group), sep = "=", collapse = "; "),
  Sex = paste(names(table(sub_meta$Sex)), table(sub_meta$Sex), sep = "=", collapse = "; "),
  Region_custom = paste(names(table(sub_meta$Region_custom)), table(sub_meta$Region_custom), sep = "=", collapse = "; "),
  Geography = paste(names(table(sub_meta$Geography)), table(sub_meta$Geography), sep = "=", collapse = "; "),
  Fixed_Effects = paste(fixed_effects, collapse = ",")
)
write.csv(run_info, file.path(out_dir, "run_info.csv"), row.names = FALSE)
print(run_info)

if (sum(sub_meta$Disease_Status == "Dental Caries") < 10 || sum(sub_meta$Disease_Status == "Healthy Control") < 10) {
  stop("Insufficient case/control samples")
}

cat("=== Running MaAsLin3 ===\n")
cat("Fixed effects:", paste(fixed_effects, collapse = ", "), "\n")

maaslin3(
  input_data = t(sub_abund),
  input_metadata = sub_meta,
  output = out_dir,
  fixed_effects = fixed_effects,
  random_effects = NULL,
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

res <- read.delim(file.path(out_dir, "all_results.tsv"), stringsAsFactors = FALSE)
disease_res <- res %>%
  filter(metadata == "Disease_Status") %>%
  mutate(Direction = ifelse(coef > 0, "Enriched in Dental Caries", "Depleted in Dental Caries")) %>%
  arrange(qval_individual, pval_individual)
write.table(disease_res, file.path(out_dir, "disease_status_results.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)

sig <- disease_res %>% filter(!is.na(qval_individual), qval_individual < 0.1)
write.table(sig, file.path(out_dir, "significant_disease_status_q0.1.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)

summary <- data.frame(
  BioProject = "PRJNA383868",
  Disease = "Dental Caries",
  Cases = sum(sub_meta$Disease_Status == "Dental Caries"),
  Controls = sum(sub_meta$Disease_Status == "Healthy Control"),
  Fixed_Effects = paste(fixed_effects, collapse = ","),
  Tested_Features = n_distinct(disease_res$feature),
  Significant_q0.1 = nrow(sig),
  Enriched_q0.1 = sum(sig$coef > 0),
  Depleted_q0.1 = sum(sig$coef < 0)
)
write.csv(summary, file.path(out_dir, "summary.csv"), row.names = FALSE)
print(summary)

if (nrow(sig) > 0) {
  print(sig %>% select(feature, coef, pval_individual, qval_individual, Direction) %>% head(50))
}
cat("=== Done ===\n")
