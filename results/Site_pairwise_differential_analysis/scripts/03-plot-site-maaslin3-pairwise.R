suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(data.table)
  library(viridis); library(cowplot); library(patchwork); library(ggridges)
  library(qiime2R)
})

base_dir      <- "/home/yanglun/YSD/global_oral/global_oral_v2/oral-v2-02/03-site-differential-analysis"
input_dir     <- file.path(base_dir, "input")
output_dir    <- file.path(base_dir, "output")
results_dir   <- file.path(output_dir, "maaslin3_site_pairwise")
figures_dir   <- file.path(output_dir, "figures")

summary_file <- file.path(results_dir, "summary_diff_counts.tsv")
all_results_file <- file.path(results_dir, "all_pairwise_results.tsv")
genus_qza <- file.path(input_dir, "genus_table_all_qc.qza")
metadata_file <- file.path(input_dir, "metadata_oral_qc_update.tsv")

if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

# Helper: save both PDF and PNG for every figure
save_both <- function(path_pdf, plot, width, height, dpi = 300) {
  ggsave(path_pdf, plot, width = width, height = height)
  png_path <- sub("\\.pdf$", ".png", path_pdf)
  ggsave(png_path, plot, width = width, height = height, dpi = dpi, bg = "white")
}

# All oral sites (for analysis, must match 01 script)
site_order <- c("saliva", "tooth", "oral mucosa", "oral swab",
                "tongue", "gingiva")

# Only these 5 sites are shown in plots
plot_sites <- c("saliva", "oral swab", "oral mucosa", "tongue", "tooth")

site_colors <- c(
  "saliva"       = "#E41A1C",
  "oral swab"    = "#FF7F00",
  "oral mucosa"  = "#377EB8",
  "tongue"       = "#4DAF4A",
  "tooth"        = "#984EA3"
)

# Display names for plots
site_display <- c(
  "saliva"       = "Saliva",
  "oral swab"    = "Oral Swab",
  "oral mucosa"  = "Oral Mucosa",
  "tongue"       = "Tongue",
  "tooth"        = "Tooth",
  "gingiva"      = "Gingiva"
)

alpha <- 0.05
group_col <- "site_cleaned"
site_group_col <- "site_group_cleaned"

cat("==========================================\n")
cat("Site Differential Analysis Visualization\n")
cat("Start:", format(Sys.time()), "\n")
cat("==========================================\n\n")

# 1. Load data
cat("[1/7] Loading data...\n")
summary_counts <- fread(summary_file, data.table = FALSE)
all_results <- fread(all_results_file, data.table = FALSE, quote = "\"", fill = TRUE)

qza_obj <- read_qza(genus_qza)
genus_counts <- as.data.frame(as.matrix(t(qza_obj$data)))

metadata_raw <- fread(metadata_file, quote = "", fill = TRUE, check.names = FALSE, data.table = FALSE)
metadata_clean <- metadata_raw %>%
  filter(!is.na(sample_id), sample_id != "",
         site_group_cleaned == "Oral",
         site_cleaned != "", !is.na(site_cleaned),
         site_cleaned %in% site_order) %>%
  distinct(sample_id, .keep_all = TRUE)

common_samples <- intersect(metadata_clean$sample_id, rownames(genus_counts))
genus_counts <- genus_counts[common_samples, ]
metadata_clean <- metadata_clean[metadata_clean$sample_id %in% common_samples, ]
rownames(metadata_clean) <- metadata_clean$sample_id

# 2. Figure - Bubble matrix (only 5 plot sites)
cat("\n[2/7] Figure - Bubble matrix...\n")
bubble_data <- summary_counts %>%
  select(Site1, Site2, N_Significant) %>%
  rename(region1 = Site1, region2 = Site2, freq = N_Significant) %>%
  filter(region1 %in% plot_sites, region2 %in% plot_sites)

sites_in_data <- intersect(plot_sites, unique(c(bubble_data$region1, bubble_data$region2)))
bubble_data <- filter(bubble_data, region1 %in% sites_in_data, region2 %in% sites_in_data)

site_idx <- setNames(seq_along(sites_in_data), sites_in_data)
bubble_data <- bubble_data %>% filter(site_idx[region1] <= site_idx[region2])

diagonal_labels <- data.frame(
  region1 = sites_in_data,
  region2 = sites_in_data,
  label = site_display[sites_in_data]
)

max_freq <- max(bubble_data$freq, na.rm = TRUE)
size_limits <- c(0, ifelse(max_freq > 0, max_freq * 1.1, 1))

fig_bubble <- ggplot() +
  geom_point(data = bubble_data,
             aes(x = factor(region1, levels = sites_in_data),
                 y = factor(region2, levels = rev(sites_in_data)), size = freq),
             shape = 21, fill = "grey70", stroke = 0.5) +
  geom_text(data = bubble_data,
            aes(x = factor(region1, levels = sites_in_data),
                y = factor(region2, levels = rev(sites_in_data)), label = freq),
            size = 6, vjust = 3.5, fontface = "bold") +
  geom_label(data = diagonal_labels,
             aes(x = factor(region1, levels = sites_in_data),
                 y = factor(region2, levels = rev(sites_in_data)),
                 label = label, fill = region1),
             size = 6, fontface = "bold", label.size = 0.5) +
  scale_radius(range = c(4, 20), limits = size_limits) +
  scale_fill_manual(values = site_colors) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.2))) +
  scale_y_discrete(expand = expansion(mult = c(0.2, 0.1))) +
  guides(size = guide_legend(title.position = "top", title = "Number of Differently\nAbundant Genera",
                             label.position = "bottom", nrow = 1,
                             override.aes = list(fill = "grey70", stroke = 0.5)),
         fill = "none") +
  labs(x = NULL, y = NULL) +
  coord_fixed(ratio = 1, clip = "off") +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.75, 0.85),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.background = element_rect(fill = "white", colour = "grey60", linewidth = 0.3),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(20, 20, 20, 40)
  )

# External label strip appended to the bottom via patchwork
label_strip <- ggplot(data.frame(region = sites_in_data)) +
  geom_tile(aes(x = factor(region, levels = sites_in_data), y = 1, fill = region)) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.2))) +
  scale_fill_manual(values = site_colors) +
  theme_void() +
  theme(legend.position = "none")

fig_bubble_combined <- fig_bubble / label_strip + plot_layout(heights = c(20, 1))

save_both(file.path(figures_dir, "figureA_bubble_matrix.pdf"), fig_bubble_combined, width = 10, height = 10)

# 3. Prepare Data
cat("\n[3/7] Preparing data...\n")
qcol <- if ("qval_individual" %in% names(all_results)) "qval_individual" else "qval"

if ("metadata" %in% names(all_results)) {
  sig_results <- all_results[all_results$metadata == group_col & !is.na(all_results[[qcol]]) & all_results[[qcol]] < alpha, ]
} else {
  sig_results <- all_results[!is.na(all_results[[qcol]]) & all_results[[qcol]] < alpha, ]
}
# Filter to only the 5 plot sites
sig_results <- sig_results %>%
  filter(Site1 %in% plot_sites & Site2 %in% plot_sites)

top_taxa <- sig_results %>% count(feature, name = "n_sig") %>% arrange(desc(n_sig)) %>% head(50)

# 4. Mean Abundances
cat("\n[4/7] Calculating mean abundances...\n")
row_sums <- rowSums(genus_counts)
row_sums[row_sums == 0] <- 1  # prevent division by zero
genus_rel <- sweep(genus_counts, 1, row_sums, "/")

genus_rel_df <- as.data.frame(genus_rel)
genus_rel_df$site_cleaned <- metadata_clean[rownames(genus_rel_df), "site_cleaned"]

# Use non-zero mean abundance (only detected samples) for ordering and display
taxa_means <- genus_rel_df %>%
  select(-site_cleaned) %>%
  pivot_longer(everything(), names_to = "taxon", values_to = "abundance") %>%
  filter(abundance > 0) %>%
  group_by(taxon) %>%
  summarise(mean_abundance = mean(abundance), .groups = "drop") %>%
  arrange(desc(mean_abundance))

mean_by_site <- genus_rel_df %>%
  pivot_longer(cols = -site_cleaned, names_to = "taxon", values_to = "abundance") %>%
  filter(abundance > 0) %>%
  group_by(site_cleaned, taxon) %>% summarise(mean_abundance = mean(abundance), .groups = "drop")

# 5. Figure - Heatmaps + Barplot
cat("\n[5/7] Figure - Heatmaps + Barplot...\n")
taxa_order <- taxa_means %>% filter(taxon %in% top_taxa$feature) %>% pull(taxon)

extract_strict_genus_only <- function(taxon_name) {
  if (grepl("g__", taxon_name)) {
    extracted <- sub(".*g__([A-Za-z0-9_.-]+).*", "\\1", taxon_name)
    extracted <- gsub("_+$", "", extracted)
    if (nzchar(extracted)) return(extracted)
  }
  tail(strsplit(taxon_name, "[_;]")[[1]], 1)
}

taxa_display_labels <- vapply(taxa_order, extract_strict_genus_only, character(1))
names(taxa_display_labels) <- taxa_order

pval_for_plot <- sig_results %>%
  filter(feature %in% taxa_order) %>%
  mutate(pair = paste(gsub(" ", "_", Site1), gsub(" ", "_", Site2), sep = "_vs_")) %>%
  group_by(feature, pair) %>%
  summarise(qval = min(.data[[qcol]], na.rm = TRUE), .groups = "drop")

pair_levels <- unique(pval_for_plot$pair)

qval_long <- pval_for_plot %>%
  select(feature, pair, qval) %>%
  right_join(expand_grid(feature = taxa_order, pair = pair_levels), by = c("feature", "pair"))

qval_long <- qval_long %>%
  mutate(log_q = -log10(qval),
         log_q = ifelse(is.infinite(log_q) | log_q > 15, 15, log_q))

adj_p_heatmap <- ggplot(qval_long) +
  geom_tile(aes(x = factor(pair, levels = pair_levels), y = factor(feature, levels = rev(taxa_order)), fill = log_q), colour = "white") +
  scale_fill_gradient(low = "#FDDBC7", high = "#B2182B", na.value = "grey90",
                      name = "-log10\n(Adj P Value)", limits = c(0, 15), breaks = seq(0, 15, by = 5)) +
  theme_minimal() + coord_fixed(ratio = 1) +
  theme(axis.text.y = element_text(size = 20, face = "italic", hjust = 1),
        axis.text.x = element_text(size = 22, angle = 45, hjust = 1, face = "bold"),
        axis.title = element_blank(), axis.ticks = element_blank(), legend.position = "top",
        legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.8)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0), labels = taxa_display_labels)

mean_heatmap_data <- filter(mean_by_site, taxon %in% taxa_order, site_cleaned %in% sites_in_data)

mean_heatmap <- ggplot(mean_heatmap_data) +
  geom_tile(aes(x = factor(site_cleaned, levels = sites_in_data), y = factor(taxon, levels = rev(taxa_order)), fill = log10(mean_abundance)), colour = "white") +
  scale_fill_viridis(option = "mako", direction = -1, name = "log10(Mean\nAbundance)") +
  theme_minimal() + coord_fixed(ratio = 1) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 22, angle = 45, hjust = 1, face = "bold"),
        axis.title = element_blank(), axis.ticks = element_blank(), legend.position = "top",
        legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.8)) +
  scale_x_discrete(expand = c(0, 0), labels = site_display) + scale_y_discrete(expand = c(0, 0))

bar_df <- filter(taxa_means, taxon %in% taxa_order)
bar_breaks <- pretty(c(0, max(bar_df$mean_abundance, na.rm = TRUE)), n = 4)

mean_barplot <- ggplot(bar_df) +
  geom_col(aes(x = mean_abundance, y = factor(taxon, levels = rev(taxa_order))), fill = "steelblue", width = 0.8) +
  scale_x_continuous(trans = "sqrt", breaks = bar_breaks, expand = expansion(mult = c(0, 0.05))) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Mean Relative\nAbundance (sqrt)") +
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 20),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.8),
        plot.margin = margin(t = 0, r = 20, b = 0, l = 0))

fig_heatmaps <- adj_p_heatmap | mean_heatmap | mean_barplot
fig_heatmaps <- fig_heatmaps + plot_layout(widths = c(15, 6, 7))

save_both(file.path(figures_dir, "figureB_heatmaps.pdf"), fig_heatmaps, width = 22, height = 30)

# 6. Figure - Density ridgeline (top 5 by abundance)
cat("\n[6/7] Figure - Density ridgeline...\n")
top5_by_abundance <- taxa_order[1:5]

top5_genus_labels <- sapply(top5_by_abundance, extract_strict_genus_only)
names(top5_genus_labels) <- top5_by_abundance

genus_long <- genus_rel_df %>%
  select(all_of(c(top5_by_abundance, "site_cleaned"))) %>%
  pivot_longer(cols = -site_cleaned, names_to = "taxon", values_to = "abundance") %>%
  filter(abundance > 0) %>%
  mutate(taxon = factor(taxon, levels = top5_by_abundance))

fig_ridgeline <- ggplot(genus_long) +
  geom_density_ridges(aes(x = log10(abundance), y = factor(site_cleaned, levels = rev(sites_in_data)), fill = site_cleaned), quantile_lines = TRUE, quantiles = 2, alpha = 0.7) +
  scale_fill_manual(values = site_colors) +
  facet_grid(~ taxon, labeller = as_labeller(top5_genus_labels)) +
  labs(x = "log10(Relative Abundance)", y = NULL) + theme_minimal() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 16, face = "bold.italic"),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
        plot.margin = margin(10, 10, 10, 10))

save_both(file.path(figures_dir, "figureC_ridgeline.pdf"), fig_ridgeline, width = 14, height = 6)

# 7. Combined Figure
cat("\n[7/7] Combining figures...\n")
fig_combined <- wrap_elements(full = fig_bubble_combined) /
                 wrap_elements(full = fig_heatmaps) /
                 wrap_elements(full = fig_ridgeline) +
  plot_layout(heights = c(1.0, 3.0, 0.7)) +
  plot_annotation(tag_levels = 'A')

save_both(file.path(figures_dir, "figure_combined.pdf"), fig_combined, width = 26, height = 36)
cat("\n==========================================\nComplete!\n")