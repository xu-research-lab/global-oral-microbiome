#!/usr/bin/env Rscript
# ============================================================================
# Site pairwise plots: bubble matrix, heatmap, abundance heatmap, ridgeline
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggridges)
})

setwd("/beegfs/huangxiaochang/global_oral/new_analysis/maaslin3_v2/site_pairwise")

results_dir <- "results"
figures_dir <- "figures"
alpha <- 0.05
qcol <- "qval_individual"

dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

save_both <- function(name, plot, width, height) {
  ggsave(file.path(figures_dir, paste0(name, ".pdf")), plot, width = width, height = height)
  ggsave(file.path(figures_dir, paste0(name, ".png")), plot, width = width, height = height, dpi = 300)
}

# ============================================================================
# Load data
# ============================================================================
cat("Loading data...\n")
summary_counts <- read.delim(file.path(results_dir, "summary_diff_counts.tsv"), stringsAsFactors = FALSE)
all_results <- read.delim(file.path(results_dir, "all_pairwise_results.tsv"), stringsAsFactors = FALSE)

meta <- read.delim("../../metadata_oral_qc_dedup.tsv", header = TRUE, sep = "\t",
                    fill = TRUE, comment.char = "", stringsAsFactors = FALSE, quote = "")
genus_abund <- read.table("../genus_table_relabund.tsv", header = TRUE, sep = "\t",
                          row.names = 1, check.names = FALSE,
                          comment.char = "", stringsAsFactors = FALSE)

# Apply site map
site_map <- c("saliva" = "Saliva", "tongue" = "Tongue", "tooth" = "Tooth",
              "oral swab" = "Oral_swab", "oral mucosa" = "Oral_Mucosa",
              "gingiva" = "Gingiva", "other" = "Others")
meta$site_cleaned <- ifelse(meta$site_cleaned %in% names(site_map),
                            site_map[meta$site_cleaned], meta$site_cleaned)
meta$age_cleaned_group <- gsub("＞", ">", meta$age_cleaned_group)

TARGET_SITES <- c("Saliva", "Tooth", "Tongue", "Oral_swab", "Oral_Mucosa")
site_colors <- c("Saliva" = "#4C78A8", "Tooth" = "#F58518", "Tongue" = "#E45756",
                 "Oral_swab" = "#72B7B2", "Oral_Mucosa" = "#54A24B")

# Extract genus name
extract_genus <- function(x) {
  x <- as.character(x)
  if (grepl("g__", x)) {
    g <- regmatches(x, regexpr("g__[^;]*", x))
    x <- sub("^g__", "", g)
  }
  x <- sub("_.+$", "", x)
  ifelse(is.na(x) | x == "" | x == "g__", NA_character_, x)
}

# ============================================================================
# 1. Bubble matrix
# ============================================================================
cat("\n=== Figure 1: Bubble matrix ===\n")

bubble_data <- summary_counts %>%
  filter(Site1 != Site2) %>%
  select(Site1, Site2, N_Significant)

# Sort by total significant taxa
site_order <- TARGET_SITES
site_weights <- setNames(
  sapply(TARGET_SITES, function(s) sum(bubble_data$N_Significant[bubble_data$Site1 == s | bubble_data$Site2 == s])),
  TARGET_SITES
)
site_order <- names(sort(site_weights, decreasing = TRUE))

bubble_data <- bubble_data %>%
  mutate(Site1 = factor(Site1, levels = site_order),
         Site2 = factor(Site2, levels = rev(site_order)))

max_n <- max(bubble_data$N_Significant)

diag_labels <- data.frame(
  Site1 = factor(site_order, levels = site_order),
  Site2 = factor(rev(site_order), levels = rev(site_order)),
  label = gsub("_", "\n", site_order)
)

fig1 <- ggplot() +
  geom_point(data = bubble_data,
             aes(x = Site1, y = Site2, size = N_Significant),
             shape = 21, fill = "grey70", stroke = 0.5) +
  geom_text(data = bubble_data,
            aes(x = Site1, y = Site2, label = N_Significant),
            size = 5, vjust = 3.5, fontface = "bold") +
  geom_label(data = diag_labels,
             aes(x = Site1, y = Site2, label = label, fill = Site1),
             size = 5, fontface = "bold", label.size = 0.5) +
  scale_radius(range = c(4, 20), limits = c(0, max(max_n, 1) * 1.1)) +
  scale_fill_manual(values = site_colors) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.2))) +
  scale_y_discrete(expand = expansion(mult = c(0.2, 0.1))) +
  guides(size = guide_legend(title = "N Significant\nGenera")) +
  labs(x = NULL, y = NULL,
       title = "Pairwise Site-differential Genera",
       subtitle = paste0("MaAsLin3 abundance model, q < ", alpha)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        axis.text = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

save_both("bubble_matrix", fig1, width = 8, height = 7)

# ============================================================================
# 2. Heatmap: top 50 taxa q-values across all pairs
# ============================================================================
cat("\n=== Figure 2: Heatmap ===\n")

sig_results <- all_results %>%
  filter(metadata == "site_cleaned", !is.na(.data[[qcol]]), .data[[qcol]] < alpha)

top_taxa <- sig_results %>%
  count(feature, name = "n_sig") %>%
  filter(!is.na(feature)) %>%
  arrange(desc(n_sig)) %>%
  head(50)

if (nrow(top_taxa) > 0) {
  taxa_order <- top_taxa$feature
  taxa_labels <- sapply(taxa_order, extract_genus)
  names(taxa_labels) <- taxa_order

  # q-value matrix: taxon x comparison
  all_pairs <- unique(sig_results$Comparison)
  qval_mat <- sig_results %>%
    filter(feature %in% taxa_order) %>%
    group_by(feature, Comparison) %>%
    summarise(qval = min(.data[[qcol]], na.rm = TRUE), .groups = "drop")

  qval_wide <- expand.grid(
    feature = taxa_order,
    Comparison = all_pairs,
    stringsAsFactors = FALSE
  ) %>%
    left_join(qval_mat, by = c("feature", "Comparison")) %>%
    mutate(
      log_q = -log10(qval),
      log_q = ifelse(is.infinite(log_q) | log_q > 15, 15, log_q),
      log_q = ifelse(is.na(log_q), 0, log_q),
      feature = factor(feature, levels = rev(taxa_order)),
      Comparison = factor(Comparison, levels = all_pairs)
    )

  fig2 <- ggplot(qval_wide, aes(x = Comparison, y = feature, fill = log_q)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "#FDDBC7", high = "#B2182B", na.value = "grey90",
                        name = "-log10(q)") +
    scale_y_discrete(labels = taxa_labels) +
    labs(title = "Site-Differential Genera (top 50)",
         x = NULL, y = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
          axis.text.y = element_text(size = 9, face = "italic"),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

  h <- max(10, nrow(top_taxa) * 0.22)
  save_both("top50_heatmap", fig2, width = 12, height = h)
}

# ============================================================================
# 3. Mean abundance heatmap per site
# ============================================================================
cat("\n=== Figure 3: Mean abundance by site ===\n")

sample_ids <- trimws(as.character(meta$sample_id))
keep_ids <- !is.na(sample_ids) & sample_ids != "" & !duplicated(sample_ids)
meta <- meta[keep_ids, , drop = FALSE]
sample_ids <- sample_ids[keep_ids]

common <- intersect(sample_ids, colnames(genus_abund))
meta <- meta[match(common, sample_ids), , drop = FALSE]
rownames(meta) <- common
meta <- meta[meta$site_cleaned %in% TARGET_SITES, , drop = FALSE]

genus_rel <- sweep(genus_abund[, rownames(meta)], 2, colSums(genus_abund[, rownames(meta)]), "/")

# Mean abundance per site for top taxa
if (nrow(top_taxa) > 0) {
  taxa_for_abund <- intersect(taxa_order, rownames(genus_rel))
  mean_by_site <- sapply(TARGET_SITES, function(s) {
    s_samples <- rownames(meta)[meta$site_cleaned == s]
    rowMeans(genus_rel[taxa_for_abund, s_samples, drop = FALSE])
  })
  mean_df <- as.data.frame(mean_by_site) %>%
    mutate(feature = taxa_for_abund) %>%
    pivot_longer(-feature, names_to = "site", values_to = "mean_abund")

  mean_df <- mean_df %>%
    mutate(
      feature = factor(feature, levels = rev(taxa_order)),
      site = factor(site, levels = TARGET_SITES),
      log10_abund = log10(ifelse(mean_abund <= 0, 1e-10, mean_abund))
    )

  fig3 <- ggplot(mean_df, aes(x = site, y = feature, fill = log10_abund)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(option = "mako", direction = -1, name = "log10(mean abund)") +
    scale_y_discrete(labels = taxa_labels) +
    labs(title = "Mean Relative Abundance by Site", x = NULL, y = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
          axis.text.y = element_text(size = 9, face = "italic"),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

  save_both("mean_abundance_heatmap", fig3, width = 8, height = h)
}

# ============================================================================
# 4. Ridgeline: top 5 genera by abundance across sites
# ============================================================================
cat("\n=== Figure 4: Ridgeline ===\n")

top5_features <- taxa_order[1:min(5, length(taxa_order))]
top5_labels <- sapply(top5_features, extract_genus)
names(top5_labels) <- top5_features

ridge_data <- as.data.frame(t(genus_rel[top5_features, , drop = FALSE]))
ridge_data$site <- meta[rownames(ridge_data), "site_cleaned"]

ridge_long <- ridge_data %>%
  pivot_longer(-site, names_to = "feature", values_to = "abundance") %>%
  mutate(
    feature = factor(feature, levels = top5_features),
    site = factor(site, levels = TARGET_SITES),
    log10_abund = log10(ifelse(abundance <= 0, 1e-10, abundance))
  )

fig4 <- ggplot(ridge_long,
               aes(x = log10_abund, y = factor(site, levels = rev(TARGET_SITES)),
                   fill = site)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2, alpha = 0.7) +
  scale_fill_manual(values = site_colors) +
  facet_grid(~ feature, labeller = as_labeller(top5_labels)) +
  labs(x = "log10(Relative Abundance)", y = NULL,
       title = "Top 5 Genera Abundance Distribution by Site") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 14, face = "bold.italic"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

save_both("ridgeline_top5", fig4, width = 14, height = 6)

cat("\n=== Plotting done ===\n")
