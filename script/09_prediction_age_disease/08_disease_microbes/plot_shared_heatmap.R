#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

BASE <- "/beegfs/huangxiaochang/global_oral/new_analysis/maaslin3"
RESULTS_DIR <- file.path(BASE, "disease_results")
OUT_DIR <- file.path(BASE, "disease_shared_plots")
TAXONOMY <- "/beegfs/huangxiaochang/global_oral/taxonomy.tsv"
Q_THRESHOLD <- 0.05

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

disease_labels <- c(
  "Dental_Caries" = "Dental Caries",
  "Head_and_Neck_Cancer" = "H&N Cancer",
  "HIV" = "HIV",
  "OSCC" = "OSCC",
  "Periodontal_Disease" = "Periodontal Dz",
  "RA" = "RA",
  "Breast_Cancer" = "Breast Cancer",
  "CRC" = "CRC",
  "IBD" = "IBD",
  "Lung_Cancer" = "Lung Cancer"
)

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
    slice(1) %>%
    ungroup()
  setNames(df$Family, df$Genus)
}

collect_significant <- function() {
  files <- list.files(RESULTS_DIR, pattern = "^all_results.tsv$",
                      recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) stop("No all_results.tsv files found: ", RESULTS_DIR)

  rows <- list()
  for (f in sort(files)) {
    disease_dir <- basename(dirname(f))
    disease <- ifelse(disease_dir %in% names(disease_labels),
                      disease_labels[[disease_dir]], gsub("_", " ", disease_dir))
    res <- read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
    sig <- res %>%
      filter(metadata == "Disease_Status",
             !is.na(qval_individual),
             qval_individual <= Q_THRESHOLD) %>%
      mutate(
        DiseaseDir = disease_dir,
        Disease = disease,
        Genus = clean_genus(feature),
        Direction = ifelse(coef > 0, "+", "-")
      ) %>%
      filter(!is.na(Genus)) %>%
      group_by(DiseaseDir, Disease, Genus) %>%
      arrange(qval_individual, desc(abs(coef)), .by_group = TRUE) %>%
      slice(1) %>%
      ungroup()
    if (nrow(sig) > 0) rows[[disease_dir]] <- sig
  }
  bind_rows(rows)
}

sig_all <- collect_significant()
if (nrow(sig_all) == 0) stop("No significant disease-associated genera found.")

shared_genera <- sig_all %>%
  distinct(Genus, Disease) %>%
  count(Genus, name = "DiseaseCount") %>%
  filter(DiseaseCount >= 2) %>%
  pull(Genus)
if (length(shared_genera) == 0) stop("No genera shared by >=2 diseases.")

shared_df <- sig_all %>%
  filter(Genus %in% shared_genera) %>%
  mutate(
    q_plot = ifelse(qval_individual <= 0, 1e-300, qval_individual),
    signed_logq = pmin(-log10(q_plot), 10) * ifelse(coef > 0, 1, -1)
  )

disease_order <- shared_df %>%
  distinct(Disease, Genus) %>%
  count(Disease, name = "SharedGenera") %>%
  arrange(desc(SharedGenera), Disease) %>%
  pull(Disease)

genus_order <- shared_df %>%
  distinct(Genus, Disease) %>%
  mutate(DiseaseRank = match(Disease, disease_order)) %>%
  group_by(Genus) %>%
  summarize(
    DiseaseCount = n(),
    FirstDisease = min(DiseaseRank),
    Pattern = paste(sort(DiseaseRank), collapse = "."),
    .groups = "drop"
  ) %>%
  arrange(desc(DiseaseCount), FirstDisease, Pattern, Genus) %>%
  pull(Genus)

family_map <- build_family_map()
family_labels <- family_map[genus_order]
family_labels[is.na(family_labels)] <- "f__"
genus_labels <- paste0("g__", genus_order)
family_width <- max(nchar(family_labels), na.rm = TRUE) + 2
axis_labels <- paste0(str_pad(family_labels, width = family_width, side = "right"),
                      genus_labels)
names(axis_labels) <- genus_order

plot_df <- expand.grid(
  Genus = genus_order,
  Disease = disease_order,
  stringsAsFactors = FALSE
) %>%
  left_join(
    shared_df %>%
      group_by(Genus, Disease) %>%
      arrange(qval_individual, desc(abs(coef)), .by_group = TRUE) %>%
      slice(1) %>%
      ungroup() %>%
      select(Genus, Disease, Direction, signed_logq, qval_individual, coef),
    by = c("Genus", "Disease")
  ) %>%
  mutate(
    Direction = ifelse(is.na(Direction), "", Direction),
    signed_logq = ifelse(is.na(signed_logq), 0, signed_logq),
    Genus = factor(Genus, levels = rev(genus_order)),
    Disease = factor(Disease, levels = disease_order)
  )

write.table(plot_df, file.path(OUT_DIR, "shared_genera_heatmap_data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(
  data.frame(Genus = genus_order,
             Family = family_labels,
             AxisLabel = axis_labels,
             stringsAsFactors = FALSE),
  file.path(OUT_DIR, "shared_genera_family_labels.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Compressed dimensions relative to the initial plot.
heat_w <- max(6.8, length(disease_order) * 0.78)
heat_h <- max(4.8, length(genus_order) * 0.145)

p <- ggplot(plot_df, aes(x = Disease, y = Genus, fill = signed_logq)) +
  geom_tile(color = "grey82", linewidth = 0.18) +
  geom_text(aes(label = Direction), size = 3.2, fontface = "bold") +
  scale_y_discrete(labels = axis_labels) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0, limits = c(-10, 10),
    name = "signed -log10(q)"
  ) +
  labs(
    title = "Shared Disease-associated Genera",
    subtitle = paste0("Genera significant in at least two diseases (q <= ", Q_THRESHOLD, ")"),
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                               size = 9.5, face = "bold"),
    axis.text.y = element_text(size = 6.3, family = "mono", hjust = 0),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey55", fill = NA, linewidth = 0.4),
    legend.position = "bottom",
    legend.key.width = unit(1.6, "cm"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 9.5),
    plot.margin = margin(6, 6, 6, 6)
  )

ggsave(file.path(OUT_DIR, "shared_genera_heatmap.pdf"), p,
       width = heat_w, height = heat_h, limitsize = FALSE)
ggsave(file.path(OUT_DIR, "shared_genera_heatmap.png"), p,
       width = heat_w, height = heat_h, limitsize = FALSE, dpi = 300)

cat("Saved heatmap to ", OUT_DIR, "\n", sep = "")
cat("Shared genera in >=2 diseases: ", length(genus_order), "\n", sep = "")
cat(sprintf("Dimensions: %.2f x %.2f inches\n", heat_w, heat_h))
