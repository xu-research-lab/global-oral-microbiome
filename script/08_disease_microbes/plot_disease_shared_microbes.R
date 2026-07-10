#!/usr/bin/env Rscript
# ============================================================================
# MaAsLin3 v2 — Disease Shared Microbes Plots
# Heatmap + UpSet of genera shared across diseases (all sites).
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
})

setwd("/beegfs/huangxiaochang/global_oral/new_analysis/maaslin3_v2")

RESULTS_DIR <- "disease_results"
OUT_DIR <- "disease_shared_plots"
Q_THRESHOLD <- 0.05
TAXONOMY <- "/beegfs/huangxiaochang/global_oral/taxonomy.tsv"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Helpers
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

make_axis_labels <- function(genera, family_map) {
  family <- family_map[genera]
  family[is.na(family)] <- "f__"
  genus <- paste0("g__", genera)
  family_width <- max(nchar(family), na.rm = TRUE) + 2
  labels <- paste0(str_pad(family, width = family_width, side = "right"), genus)
  names(labels) <- genera
  labels
}

disease_labels <- c(
  "Dental_Caries" = "Dental Caries", "Head_and_Neck_Cancer" = "H&N Cancer",
  "HIV" = "HIV", "OSCC" = "OSCC", "Periodontal_Disease" = "Periodontal Dz",
  "RA" = "RA", "Breast_Cancer" = "Breast Cancer", "CRC" = "CRC",
  "IBD" = "IBD", "Lung_Cancer" = "Lung Cancer"
)

# ============================================================================
# Count BioProjects per disease (matched: BP has both Case + HC)
# ============================================================================
cat("Counting matched BioProjects per disease...\n")
metadata <- read.table("../metadata_oral_qc_dedup.tsv", header = TRUE, sep = "\t",
                       fill = TRUE, comment.char = "", stringsAsFactors = FALSE, quote = "")
site_map <- c("saliva" = "Saliva", "tongue" = "Tongue", "tooth" = "Tooth",
              "oral swab" = "Oral_swab", "oral mucosa" = "Oral_Mucosa",
              "gingiva" = "Gingiva", "other" = "Others")
metadata$site_cleaned <- ifelse(metadata$site_cleaned %in% names(site_map),
                                site_map[metadata$site_cleaned], metadata$site_cleaned)
metadata$status_cleaned[metadata$status_cleaned %in% c("HC", "Healthy")] <- "Healthy Control"

bp_counts <- list()
for (dis_dir in names(disease_labels)) {
  dis_original <- gsub("_", " ", dis_dir)
  if (dis_dir == "Head_and_Neck_Cancer") dis_original <- "Head and Neck Cancer"
  if (dis_dir == "Periodontal_Disease") dis_original <- "Periodontal Disease"
  if (dis_dir == "Breast_Cancer") dis_original <- "Breast Cancer"
  if (dis_dir == "Lung_Cancer") dis_original <- "Lung Cancer"
  if (dis_dir == "Dental_Caries") dis_original <- "Dental Caries"

  disease_bp <- unique(metadata$bioproject_id[metadata$status_cleaned == dis_original])
  bp_tbl <- table(
    metadata$bioproject_id[metadata$bioproject_id %in% disease_bp &
                           metadata$status_cleaned %in% c("Healthy Control", dis_original)],
    metadata$status_cleaned[metadata$bioproject_id %in% disease_bp &
                            metadata$status_cleaned %in% c("Healthy Control", dis_original)]
  )
  if ("Healthy Control" %in% colnames(bp_tbl) && dis_original %in% colnames(bp_tbl)) {
    n_bp <- sum(bp_tbl[, "Healthy Control"] > 0 & bp_tbl[, dis_original] > 0)
  } else {
    n_bp <- 0
  }
  bp_counts[[disease_labels[dis_dir]]] <- n_bp
  cat(sprintf("  %s: %d matched BPs\n", disease_labels[dis_dir], n_bp))
}

# Build disease label with BP count suffix
make_disease_label <- function(disease_name) {
  n_bp <- bp_counts[[disease_name]]
  if (is.null(n_bp) || is.na(n_bp)) return(disease_name)
  sprintf("%s [%dBP]", disease_name, n_bp)
}

disease_labels_with_bp <- sapply(names(disease_labels), function(d) {
  setNames(make_disease_label(disease_labels[d]), disease_labels[d])
})
# Create a lookup: Disease -> Disease_with_BP
disease_to_label_bp <- setNames(
  sapply(disease_labels, make_disease_label),
  disease_labels
)

# ============================================================================
# Collect disease results
# ============================================================================
result_files <- list.files(RESULTS_DIR, pattern = "^all_results.tsv$",
                           recursive = TRUE, full.names = TRUE)
if (length(result_files) == 0) stop("No all_results.tsv files under ", RESULTS_DIR)

cat(sprintf("Found %d disease result files\n", length(result_files)))

sig_data <- list()
summary_rows <- list()

for (f in sort(result_files)) {
  disease_dir <- basename(dirname(f))
  disease <- ifelse(disease_dir %in% names(disease_labels),
                    disease_labels[[disease_dir]], gsub("_", " ", disease_dir))
  res <- read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)

  sig <- res %>%
    filter(metadata == "Disease_Status", !is.na(qval_individual),
           qval_individual <= Q_THRESHOLD) %>%
    mutate(DiseaseDir = disease_dir, Disease = disease,
           Genus = clean_genus(feature)) %>%
    filter(!is.na(Genus))

  if (nrow(sig) == 0) {
    summary_rows[[disease_dir]] <- data.frame(
      Disease = disease, SignificantGenera = 0, stringsAsFactors = FALSE)
    next
  }

  sig_best <- sig %>%
    group_by(DiseaseDir, Disease, Genus) %>%
    arrange(qval_individual, desc(abs(coef)), .by_group = TRUE) %>%
    slice(1) %>% ungroup() %>%
    mutate(Direction = ifelse(coef > 0, "+", "-"))

  sig_data[[disease_dir]] <- sig_best
  summary_rows[[disease_dir]] <- data.frame(
    Disease = disease, SignificantGenera = n_distinct(sig_best$Genus),
    stringsAsFactors = FALSE)
}

sig_all <- bind_rows(sig_data)
if (nrow(sig_all) == 0) {
  cat("No significant disease-associated genera at q <= ", Q_THRESHOLD, "\n")
  quit(save = "no", status = 0)
}

summary_df <- bind_rows(summary_rows) %>% arrange(desc(SignificantGenera), Disease)
write.csv(summary_df, file.path(OUT_DIR, "disease_sig_genera_counts.csv"), row.names = FALSE)
write.table(sig_all, file.path(OUT_DIR, "disease_significant_genera_q0.05.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n=== Significant genera per disease ===\n")
print(summary_df)

# ============================================================================
# Shared genera (>=2 diseases)
# ============================================================================
disease_sets <- lapply(split(sig_all$Genus, sig_all$Disease), unique)
disease_sets <- disease_sets[sapply(disease_sets, length) > 0]

genus_disease_counts <- sig_all %>%
  distinct(Genus, Disease) %>% count(Genus, name = "DiseaseCount")
shared_genera <- genus_disease_counts %>% filter(DiseaseCount >= 2) %>% pull(Genus)
cat(sprintf("\nGenera shared by >= 2 diseases: %d\n", length(shared_genera)))

if (length(shared_genera) == 0) {
  cat("No shared genera found.\n"); quit(save = "no", status = 0)
}

# ============================================================================
# Heatmap
# ============================================================================
cat("Building shared genera heatmap...\n")

shared_df <- sig_all %>%
  filter(Genus %in% shared_genera) %>%
  mutate(
    q_plot = ifelse(qval_individual <= 0, 1e-300, qval_individual),
    signed_logq = pmin(-log10(q_plot), 10) * ifelse(coef > 0, 1, -1))

disease_order_orig <- shared_df %>%
  distinct(Disease, Genus) %>% count(Disease, name = "n") %>%
  arrange(desc(n), Disease) %>% pull(Disease)

# Map to BP-annotated labels
disease_order <- disease_to_label_bp[disease_order_orig]
names(disease_order) <- disease_order_orig

genus_order_df <- shared_df %>%
  distinct(Genus, Disease) %>%
  mutate(dr = match(Disease, disease_order_orig)) %>%
  group_by(Genus) %>%
  summarise(n_diseases = n(), first_disease = min(dr),
            pattern = paste(sort(dr), collapse = "."), .groups = "drop") %>%
  arrange(desc(n_diseases), first_disease, pattern, Genus)
genus_order <- genus_order_df$Genus

plot_df <- expand.grid(Genus = genus_order, Disease = disease_order_orig,
                       stringsAsFactors = FALSE) %>%
  left_join(
    shared_df %>% group_by(Genus, Disease) %>%
      arrange(qval_individual, desc(abs(coef)), .by_group = TRUE) %>%
      slice(1) %>% ungroup() %>%
      select(Genus, Disease, Direction, signed_logq, qval_individual, coef),
    by = c("Genus", "Disease")) %>%
  mutate(
    Direction = ifelse(is.na(Direction), "", Direction),
    signed_logq = ifelse(is.na(signed_logq), 0, signed_logq),
    Genus = factor(Genus, levels = rev(genus_order)),
    Disease = factor(disease_to_label_bp[Disease], levels = disease_order))

write.table(plot_df, file.path(OUT_DIR, "shared_genera_heatmap_data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Build family-prefixed axis labels (matching reference style)
axis_labels <- make_axis_labels(genus_order, family_map)

w <- max(6.8, length(disease_order) * 0.78)
h <- max(4.8, length(genus_order) * 0.145)
p <- ggplot(plot_df, aes(x = Disease, y = Genus, fill = signed_logq)) +
  geom_tile(color = "grey82", linewidth = 0.18) +
  geom_text(aes(label = Direction), size = 3.2, fontface = "bold") +
  scale_y_discrete(labels = axis_labels) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, limits = c(-10, 10), name = "signed -log10(q)") +
  labs(title = "Shared Disease-associated Genera",
       subtitle = paste0("Genera significant in at least two diseases (q <= ", Q_THRESHOLD, ")"),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9.5, face = "bold"),
        axis.text.y = element_text(size = 6.3, family = "mono", hjust = 0),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "grey55", fill = NA, linewidth = 0.4),
        legend.position = "bottom",
        legend.key.width = unit(1.6, "cm"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        plot.subtitle = element_text(hjust = 0.5, size = 9.5),
        plot.margin = margin(6, 6, 6, 6))

ggsave(file.path(OUT_DIR, "shared_genera_heatmap.pdf"), p, width = w, height = h, limitsize = FALSE)
ggsave(file.path(OUT_DIR, "shared_genera_heatmap.png"), p, width = w, height = h, limitsize = FALSE, dpi = 300)

# ============================================================================
# UpSet plot
# ============================================================================
cat("Building UpSet plot...\n")

library(UpSetR)
upset_order_orig <- names(sort(sapply(disease_sets, length), decreasing = TRUE))

# Rename disease_sets with BP counts for display, then reorder
disease_sets_bp <- disease_sets
names(disease_sets_bp) <- disease_to_label_bp[names(disease_sets_bp)]
# Reorder by original order (now using BP-annotated names)
upset_order_bp <- disease_to_label_bp[upset_order_orig]
disease_sets_bp <- disease_sets_bp[upset_order_bp]

write.csv(data.frame(Disease = upset_order_orig,
                     DiseaseLabel = disease_to_label_bp[upset_order_orig],
                     SignificantGenera = sapply(disease_sets, length)),
          file.path(OUT_DIR, "disease_set_sizes.csv"), row.names = FALSE)

upset_input <- fromList(disease_sets_bp)

pdf(file.path(OUT_DIR, "disease_genera_upset.pdf"), width = 13, height = 7.2)
upset(upset_input, nsets = length(upset_order_bp), nintersects = 80,
      order.by = "freq", mb.ratio = c(0.55, 0.45),
      matrix.color = "#2C7FB8", main.bar.color = "#2C7FB8", sets.bar.color = "#41B6C4",
      mainbar.y.label = "Shared genera count", sets.x.label = "Significant genera per disease",
      point.size = 3.6, line.size = 1.1,
      text.scale = c(1.45, 1.25, 1.1, 1.1, 0.95, 0.95))
dev.off()

png(file.path(OUT_DIR, "disease_genera_upset.png"), width = 13, height = 7.2,
    units = "in", res = 300)
upset(upset_input, nsets = length(upset_order_bp), nintersects = 80,
      order.by = "freq", mb.ratio = c(0.55, 0.45),
      matrix.color = "#2C7FB8", main.bar.color = "#2C7FB8", sets.bar.color = "#41B6C4",
      mainbar.y.label = "Shared genera count", sets.x.label = "Significant genera per disease",
      point.size = 3.6, line.size = 1.1,
      text.scale = c(1.45, 1.25, 1.1, 1.1, 0.95, 0.95))
dev.off()

cat("\nSaved outputs in ", OUT_DIR, "\n", sep = "")
cat("Shared genera in >=2 diseases: ", length(genus_order), "\n", sep = "")
cat("=== Disease shared plots done ===\n")
