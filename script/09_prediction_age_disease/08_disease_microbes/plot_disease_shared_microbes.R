#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(UpSetR)
})

setwd("/beegfs/huangxiaochang/global_oral/new_analysis/maaslin3")

RESULTS_DIR <- "disease_results"
OUT_DIR <- "disease_shared_plots"
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

extract_genus <- function(x) {
  x <- as.character(x)
  if (grepl("g__", x)) {
    g <- regmatches(x, regexpr("g__[^;]*", x))
    x <- sub("^g__", "", g)
  }
  x <- trimws(x)
  if (is.na(x) || x == "" || x == "g__") return(NA_character_)
  # Match the previous analysis: collapse GTDB-style suffixed names to genus.
  x <- sub("_.+$", "", x)
  if (x == "") return(NA_character_)
  x
}

result_files <- list.files(RESULTS_DIR, pattern = "^all_results.tsv$",
                           recursive = TRUE, full.names = TRUE)
if (length(result_files) == 0) {
  stop("No all_results.tsv files found under ", RESULTS_DIR)
}

sig_data <- list()
summary_rows <- list()

for (f in sort(result_files)) {
  disease_dir <- basename(dirname(f))
  disease <- ifelse(disease_dir %in% names(disease_labels),
                    disease_labels[[disease_dir]],
                    gsub("_", " ", disease_dir))
  res <- read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
  sig <- res %>%
    filter(
      metadata == "Disease_Status",
      !is.na(qval_individual),
      qval_individual <= Q_THRESHOLD
    ) %>%
    mutate(
      DiseaseDir = disease_dir,
      Disease = disease,
      Genus = vapply(feature, extract_genus, character(1))
    ) %>%
    filter(!is.na(Genus))

  if (nrow(sig) == 0) {
    summary_rows[[disease_dir]] <- data.frame(
      Disease = disease, SignificantGenera = 0, stringsAsFactors = FALSE
    )
    next
  }

  # If multiple MaAsLin3 rows collapse to the same genus, keep the strongest.
  sig_best <- sig %>%
    group_by(DiseaseDir, Disease, Genus) %>%
    arrange(qval_individual, desc(abs(coef)), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(Direction = ifelse(coef > 0, "+", "-"))

  sig_data[[disease_dir]] <- sig_best
  summary_rows[[disease_dir]] <- data.frame(
    Disease = disease,
    SignificantGenera = n_distinct(sig_best$Genus),
    stringsAsFactors = FALSE
  )
}

sig_all <- bind_rows(sig_data)
if (nrow(sig_all) == 0) {
  stop("No significant disease-associated genera found at q <= ", Q_THRESHOLD)
}

summary_df <- bind_rows(summary_rows) %>%
  arrange(desc(SignificantGenera), Disease)
write.csv(summary_df, file.path(OUT_DIR, "disease_significant_genera_counts.csv"),
          row.names = FALSE)
write.table(sig_all, file.path(OUT_DIR, "disease_significant_genera_q0.05.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

disease_sets <- split(sig_all$Genus, sig_all$Disease)
disease_sets <- lapply(disease_sets, unique)
disease_sets <- disease_sets[sapply(disease_sets, length) > 0]

genus_disease_counts <- sig_all %>%
  distinct(Genus, Disease) %>%
  count(Genus, name = "DiseaseCount")
shared_genera <- genus_disease_counts %>%
  filter(DiseaseCount >= 2) %>%
  pull(Genus)

if (length(shared_genera) == 0) {
  stop("No genera shared by two or more diseases at q <= ", Q_THRESHOLD)
}

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

genus_order_df <- shared_df %>%
  distinct(Genus, Disease) %>%
  mutate(DiseaseRank = match(Disease, disease_order)) %>%
  group_by(Genus) %>%
  summarize(
    DiseaseCount = n(),
    FirstDisease = min(DiseaseRank),
    Pattern = paste(sort(DiseaseRank), collapse = "."),
    .groups = "drop"
  ) %>%
  arrange(desc(DiseaseCount), FirstDisease, Pattern, Genus)
genus_order <- genus_order_df$Genus

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

heat_w <- max(8, length(disease_order) * 1.05)
heat_h <- max(6, length(genus_order) * 0.24)

p <- ggplot(plot_df, aes(x = Disease, y = Genus, fill = signed_logq)) +
  geom_tile(color = "grey82", linewidth = 0.25) +
  geom_text(aes(label = Direction), size = 4.2, fontface = "bold") +
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
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                               size = 11, face = "bold"),
    axis.text.y = element_text(size = 10, face = "italic"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.5),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )

ggsave(file.path(OUT_DIR, "shared_genera_heatmap.pdf"), p,
       width = heat_w, height = heat_h, limitsize = FALSE)
ggsave(file.path(OUT_DIR, "shared_genera_heatmap.png"), p,
       width = heat_w, height = heat_h, limitsize = FALSE, dpi = 300)

upset_order <- names(sort(sapply(disease_sets, length), decreasing = TRUE))
upset_input <- fromList(disease_sets[upset_order])
write.csv(
  data.frame(Disease = upset_order,
             SignificantGenera = sapply(disease_sets[upset_order], length),
             row.names = NULL),
  file.path(OUT_DIR, "upset_set_sizes.csv"),
  row.names = FALSE
)

pdf(file.path(OUT_DIR, "disease_genera_upset.pdf"), width = 14, height = 8)
upset(
  upset_input,
  nsets = length(upset_order),
  nintersects = 80,
  order.by = "freq",
  mb.ratio = c(0.55, 0.45),
  matrix.color = "#2C7FB8",
  main.bar.color = "#2C7FB8",
  sets.bar.color = "#41B6C4",
  mainbar.y.label = "Shared genera count",
  sets.x.label = "Significant genera per disease",
  point.size = 3.8,
  line.size = 1.2,
  text.scale = c(1.6, 1.4, 1.2, 1.2, 1.0, 1.0)
)
dev.off()

png(file.path(OUT_DIR, "disease_genera_upset.png"),
    width = 14, height = 8, units = "in", res = 300)
upset(
  upset_input,
  nsets = length(upset_order),
  nintersects = 80,
  order.by = "freq",
  mb.ratio = c(0.55, 0.45),
  matrix.color = "#2C7FB8",
  main.bar.color = "#2C7FB8",
  sets.bar.color = "#41B6C4",
  mainbar.y.label = "Shared genera count",
  sets.x.label = "Significant genera per disease",
  point.size = 3.8,
  line.size = 1.2,
  text.scale = c(1.6, 1.4, 1.2, 1.2, 1.0, 1.0)
)
dev.off()

cat("Saved outputs in ", OUT_DIR, "\n", sep = "")
cat("Shared genera in >=2 diseases: ", length(genus_order), "\n", sep = "")
print(summary_df)
