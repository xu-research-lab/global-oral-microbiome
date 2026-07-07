#!/usr/bin/env Rscript
# ============================================================================
# Combined heatmap comparing original vs North-America-downsampled results
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
})

BASE <- "/beegfs/huangxiaochang/global_oral/new_analysis/maaslin3"
ORIG_DIR <- file.path(BASE, "saliva_sex_filtered_no_geo_interaction")
DS_DIR <- file.path(BASE, "saliva_sex_downsampled_na")
OUT_DIR <- file.path(BASE, "saliva_sex_downsampled_na", "plots")
TAXONOMY <- "/beegfs/huangxiaochang/global_oral/taxonomy.tsv"
Q_THRESHOLD <- 0.1

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------------
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

family_map <- build_family_map()

# ----------------------------------------------------------------------------
# Load original heatmap data
# ----------------------------------------------------------------------------
orig <- read.delim(file.path(ORIG_DIR, "saliva_Sex_heatmap_data.tsv"),
                   stringsAsFactors = FALSE, check.names = FALSE) %>%
  filter(!is.na(Genus)) %>%
  mutate(
    Method = "Original (n=7596 NA)",
    sig = qval_individual < Q_THRESHOLD & !is.na(qval_individual)
  ) %>%
  select(Genus, Geography, value, coef, qval_individual, signed_logq, sig, Method)

# ----------------------------------------------------------------------------
# Load downsampled results
# ----------------------------------------------------------------------------
ds <- read.delim(file.path(DS_DIR, "all_results_downsampled.tsv"),
                 stringsAsFactors = FALSE, check.names = FALSE) %>%
  filter(metadata == "Sex", !is.na(qval_individual)) %>%
  mutate(
    Genus = clean_genus(feature),
    Method = "Downsampled NA (n=2470)",
    q_plot = ifelse(qval_individual <= 0, 1e-300, qval_individual),
    signed_logq = pmin(-log10(q_plot), 10) * sign(coef),
    sig = qval_individual < Q_THRESHOLD
  ) %>%
  filter(!is.na(Genus)) %>%
  select(Genus, Geography, value, coef, qval_individual, signed_logq, sig, Method)

# Combine
combined <- bind_rows(orig, ds)

# Keep only genera significant in at least one method x geography
sig_genera <- combined %>%
  filter(sig) %>%
  distinct(Genus) %>%
  pull(Genus)

plot_data <- combined %>% filter(Genus %in% sig_genera)

# Order genera by total number of significant geographies across both methods
# and by mean effect size
 genus_order_df <- plot_data %>%
  filter(sig) %>%
  group_by(Genus) %>%
  summarise(
    n_sig_regions = n_distinct(Geography),
    n_sig_total = sum(sig),
    min_q = min(qval_individual, na.rm = TRUE),
    mean_abs_coef = mean(abs(coef), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sig_regions), desc(n_sig_total), desc(mean_abs_coef), min_q)

genus_order <- genus_order_df$Genus

# Also require the genus to have data in each method x geography pair for plotting
# (some combinations may be NA if not tested)
plot_data <- plot_data %>%
  mutate(
    Genus = factor(Genus, levels = rev(genus_order)),
    Geography = factor(Geography, levels = c("North America", "Asia", "Europe", "South America", "Oceania")),
    Method = factor(Method, levels = c("Original (n=7596 NA)", "Downsampled NA (n=2470)")),
    sig_label = ifelse(sig, ifelse(coef > 0, "+", "-"), ""),
    signed_logq = ifelse(is.na(signed_logq), 0, signed_logq)
  )

# Axis labels with family prefix
axis_labels <- sapply(genus_order, function(g) {
  fam <- family_map[g]
  if (is.na(fam) || fam == "f__") fam <- "f__"
  genus <- paste0("g__", g)
  fam_width <- max(nchar(family_map), na.rm = TRUE) + 2
  paste0(str_pad(fam, width = fam_width, side = "right"), genus)
})
names(axis_labels) <- genus_order
plot_data$AxisLabel <- axis_labels[as.character(plot_data$Genus)]

# ----------------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------------
p <- ggplot(plot_data, aes(x = Geography, y = Genus, fill = signed_logq)) +
  geom_tile(color = "grey85", linewidth = 0.14) +
  geom_text(aes(label = sig_label), size = 2.2, fontface = "bold") +
  facet_wrap(~ Method, nrow = 2, scales = "free_x") +
  scale_y_discrete(labels = axis_labels[genus_order]) +
  scale_x_discrete(drop = TRUE) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#D73027",
    midpoint = 0, limits = c(-10, 10), name = "signed\n-log10(q)"
  ) +
  labs(
    title = "Saliva Sex-Associated Microbiota: Original vs North America Downsampled",
    subtitle = paste0("Stratified downsampling reduces NA from 7596 to 2470; q < ", Q_THRESHOLD,
                      "; +/- indicates male coefficient direction"),
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7.5, face = "bold"),
    axis.text.y = element_text(size = 8.2, family = "mono", hjust = 0),
    strip.text = element_text(size = 9, face = "bold"),
    panel.grid = element_blank(),
    panel.spacing = unit(0.3, "lines"),
    legend.position = "bottom",
    legend.key.width = unit(0.75, "cm"),
    legend.key.height = unit(0.25, "cm"),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6.5),
    plot.title = element_text(hjust = 0.5, size = 11.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 8.0),
    plot.margin = margin(6, 6, 6, 6)
  )

n_genera <- length(genus_order)
width <- 10
height <- max(6.5, n_genera * 0.22)

ggsave(file.path(OUT_DIR, "saliva_Sex_heatmap_original_vs_downsampled.pdf"),
       p, width = width, height = height, limitsize = FALSE)
ggsave(file.path(OUT_DIR, "saliva_Sex_heatmap_original_vs_downsampled.png"),
       p, width = width, height = height, dpi = 300, limitsize = FALSE)

# Write summary of comparison
comparison_summary <- plot_data %>%
  filter(sig) %>%
  group_by(Genus, Method) %>%
  summarise(n_sig_regions = n_distinct(Geography), .groups = "drop") %>%
  pivot_wider(names_from = Method, values_from = n_sig_regions, values_fill = 0)

write.csv(comparison_summary,
          file.path(OUT_DIR, "original_vs_downsampled_summary.csv"),
          row.names = FALSE)

cat(sprintf("Saved combined heatmap to %s (%d genera, %.1f x %.1f in)\n", OUT_DIR, n_genera, width, height))
cat("Done.\n")
