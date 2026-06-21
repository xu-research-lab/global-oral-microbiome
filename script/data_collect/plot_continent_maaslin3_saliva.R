suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(data.table)
  library(viridis); library(cowplot); library(patchwork); library(ggridges)
  library(qiime2R)
})

results_dir <- "results/maaslin3_pairwise"
figures_dir <- "results/figures"
input_dir <- "input"

summary_file <- file.path(results_dir, "summary_diff_counts.tsv")
all_results_file <- file.path(results_dir, "all_pairwise_results.tsv")
genus_qza <- file.path(input_dir, "genus-table.qza")
metadata_file <- "/home/yanglun/YSD/global_oral/global_oral_v2/metadata-v2/metadata_oral_qc.tsv"

if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

# Helper: save both PDF and PNG for every figure
save_both <- function(path_pdf, plot, width, height, dpi = 300) {
  ggsave(path_pdf, plot, width = width, height = height)
  png_path <- sub("\\.pdf$", ".png", path_pdf)
  ggsave(png_path, plot, width = width, height = height, dpi = dpi, bg = "white")
}

continent_order <- c("Africa", "Asia", "Europe", "North America", "Oceania", "South America")
continent_colors <- c("Africa" = "#E41A1C", "Asia" = "#FF7F00", "Europe" = "#377EB8",
                      "North America" = "#FFD700", "Oceania" = "#984EA3", "South America" = "#A65628")
alpha <- 0.05
continent_col <- "continent_group"

cat("==========================================\n")
cat("Figure 5 Visualization (Log Cap & Layout Fix)\n")
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
  filter(!is.na(sample_id), sample_id != "", site_cleaned == "saliva",
         !is.na(continent_group), continent_group != "", continent_group != "NA") %>%
  distinct(sample_id, .keep_all = TRUE)

common_samples <- intersect(metadata_clean$sample_id, rownames(genus_counts))
genus_counts <- genus_counts[common_samples, ]
metadata_clean <- metadata_clean[metadata_clean$sample_id %in% common_samples, ]
rownames(metadata_clean) <- metadata_clean$sample_id 

# 2. Figure 5A - Bubble matrix
cat("\n[2/7] Figure 5A - Bubble matrix...\n")
bubble_data <- summary_counts %>%
  select(Continent1, Continent2, N_Significant) %>%
  rename(region1 = Continent1, region2 = Continent2, freq = N_Significant)

continents_in_data <- intersect(continent_order, unique(c(bubble_data$region1, bubble_data$region2)))
bubble_data <- filter(bubble_data, region1 %in% continents_in_data, region2 %in% continents_in_data)

continent_idx <- setNames(seq_along(continents_in_data), continents_in_data)
bubble_data <- bubble_data %>% filter(continent_idx[region1] <= continent_idx[region2])

diagonal_labels <- data.frame(
  region1 = continents_in_data,
  region2 = continents_in_data,
  label = gsub(" ", "\n", continents_in_data)
)

max_freq <- max(bubble_data$freq, na.rm = TRUE)
size_limits <- c(0, ifelse(max_freq > 0, max_freq * 1.1, 1))

fig5A <- ggplot() +
  geom_point(data = bubble_data,
             aes(x = factor(region1, levels = continents_in_data),
                 y = factor(region2, levels = rev(continents_in_data)), size = freq),
             shape = 21, fill = "grey70", stroke = 0.5) +
  geom_text(data = bubble_data,
            aes(x = factor(region1, levels = continents_in_data),
                y = factor(region2, levels = rev(continents_in_data)), label = freq),
            size = 6, vjust = 3.5, fontface = "bold") +
  geom_label(data = diagonal_labels,
             aes(x = factor(region1, levels = continents_in_data),
                 y = factor(region2, levels = rev(continents_in_data)),
                 label = label, fill = region1),
             size = 6, fontface = "bold", label.size = 0.5) +
  scale_radius(range = c(4, 20), limits = size_limits) +
  scale_fill_manual(values = continent_colors) +
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
label_strip_A <- ggplot(data.frame(region = continents_in_data)) +
  geom_tile(aes(x = factor(region, levels = continents_in_data), y = 1, fill = region)) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.2))) +
  scale_fill_manual(values = continent_colors) +
  theme_void() +
  theme(legend.position = "none")

fig5A_combined <- fig5A / label_strip_A + plot_layout(heights = c(20, 1))

save_both(file.path(figures_dir, "figure5A_bubble_matrix.pdf"), fig5A_combined, width = 10, height = 10)

# 3. Prepare Data
cat("\n[3/7] Preparing data...\n")
qcol <- if ("qval_individual" %in% names(all_results)) "qval_individual" else "qval"

if ("metadata" %in% names(all_results)) {
  sig_results <- all_results[all_results$metadata == continent_col & !is.na(all_results[[qcol]]) & all_results[[qcol]] < alpha, ]
} else {
  sig_results <- all_results[!is.na(all_results[[qcol]]) & all_results[[qcol]] < alpha, ]
}
top_taxa <- sig_results %>% count(feature, name = "n_sig") %>% arrange(desc(n_sig)) %>% head(50)

# 4. Mean Abundances
cat("\n[4/7] Calculating mean abundances...\n")
genus_rel <- sweep(genus_counts, 1, rowSums(genus_counts), "/")
genus_rel[genus_rel == 0] <- 1e-10

genus_rel_df <- as.data.frame(genus_rel)
genus_rel_df$continent_group <- metadata_clean[rownames(genus_rel_df), "continent_group"]

taxa_means <- data.frame(
  taxon = colnames(genus_rel_df)[colnames(genus_rel_df) != "continent_group"],
  mean_abundance = colMeans(genus_rel_df[, colnames(genus_rel_df) != "continent_group"])
) %>% arrange(desc(mean_abundance))

mean_by_continent <- genus_rel_df %>%
  pivot_longer(cols = -continent_group, names_to = "taxon", values_to = "abundance") %>%
  group_by(continent_group, taxon) %>% summarise(mean_abundance = mean(abundance), .groups = "drop")

# 5. Figure 5B
cat("\n[5/7] Figure 5B - Heatmaps + Barplot...\n")
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
  mutate(pair = paste(Continent1, Continent2, sep = "_vs_")) %>%
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

mean_heatmap_data <- filter(mean_by_continent, taxon %in% taxa_order)

mean_heatmap <- ggplot(mean_heatmap_data) +
  geom_tile(aes(x = factor(continent_group, levels = continents_in_data), y = factor(taxon, levels = rev(taxa_order)), fill = log10(mean_abundance)), colour = "white") +
  scale_fill_viridis(option = "mako", direction = -1, name = "log10(Mean\nAbundance)") +
  theme_minimal() + coord_fixed(ratio = 1) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 22, angle = 45, hjust = 1, face = "bold"),
        axis.title = element_blank(), axis.ticks = element_blank(), legend.position = "top",
        legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.8)) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))

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

fig5B <- adj_p_heatmap | mean_heatmap | mean_barplot
fig5B <- fig5B + plot_layout(widths = c(15, 6, 7))

save_both(file.path(figures_dir, "figure5B_heatmaps.pdf"), fig5B, width = 22, height = 30)

# 6. Figure 5C - Density ridgeline
cat("\n[6/7] Figure 5C - Density ridgeline...\n")
top5_by_abundance <- taxa_order[1:5]

top5_genus_labels <- sapply(top5_by_abundance, extract_strict_genus_only)
names(top5_genus_labels) <- top5_by_abundance

genus_long <- genus_rel_df %>%
  select(all_of(c(top5_by_abundance, "continent_group"))) %>%
  pivot_longer(cols = -continent_group, names_to = "taxon", values_to = "abundance") %>%
  mutate(taxon = factor(taxon, levels = top5_by_abundance))

fig5C <- ggplot(genus_long) +
  geom_density_ridges(aes(x = log10(abundance), y = factor(continent_group, levels = rev(continents_in_data)), fill = continent_group), quantile_lines = TRUE, quantiles = 2, alpha = 0.7) +
  scale_fill_manual(values = continent_colors) +
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

save_both(file.path(figures_dir, "figure5C_ridgeline.pdf"), fig5C, width = 14, height = 6)

# 7. Combined Figure 5
cat("\n[7/7] Combining figures...\n")
fig5_combined <- wrap_elements(full = fig5A_combined) /
                 wrap_elements(full = fig5B) /
                 wrap_elements(full = fig5C) +
  plot_layout(heights = c(1.0, 3.0, 0.7)) +
  plot_annotation(tag_levels = 'A')

save_both(file.path(figures_dir, "figure5_combined.pdf"), fig5_combined, width = 26, height = 36)
cat("\n==========================================\nComplete!\n")