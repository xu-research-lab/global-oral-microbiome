#!/usr/bin/env Rscript
# pipeline_v3: ASV insertion hot spots and cold spots
# Figure 1: genus-level circular tree — tip colour = ASV insertion density,
#           tip size = number of ASVs, outer rings = HOMD reference counts
# Figure 2: distribution of ASV distance to the nearest HOMD reference
#           (against the background references)
# Figure 3: genus-level scatter — ASV density vs HOMD reference count
suppressMessages({
  library(ggtree)
  library(ape)
  library(tidyverse)
  library(ggnewscale)
  library(ggtreeExtra)
})

PIPE_V3 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
OUT <- file.path(PIPE_V3, "12h_genus_density")
FIGOUT <- file.path(PIPE_V3, "analysis_report", "figs")
dir.create(FIGOUT, showWarnings = FALSE)

# ============ Figure 1: genus-level circular tree ============
cat("== Figure 1: genus-level circular tree ==\n")
tree <- read.tree(file.path(OUT, "genus_tree_clean.nwk"))
cat(sprintf("  tree tips: %d\n", length(tree$tip.label)))

gs <- read.delim(file.path(OUT, "genus_stats.tsv"), sep = "\t", stringsAsFactors = FALSE)
gs$genus <- as.character(gs$genus)

tip_df <- data.frame(tip = tree$tip.label, stringsAsFactors = FALSE)
tip_df$cat <- ifelse(tip_df$tip %in% gs$genus, "genus",
              ifelse(grepl("^cold_", tip_df$tip), "cold", "unclassified_asv"))
tip_df <- merge(tip_df, gs, by.x = "tip", by.y = "genus", all.x = TRUE)
tip_df$n_asv[is.na(tip_df$n_asv)] <- 0
tip_df$n_homd_ref[is.na(tip_df$n_homd_ref)] <- 0
tip_df$log10_density[is.na(tip_df$log10_density)] <- 0
tip_df$median_asv_homd_dist[is.na(tip_df$median_asv_homd_dist)] <- NA
rownames(tip_df) <- tip_df$tip

p0 <- ggtree(tree, layout = "circular", size = 0.15) %<+% tip_df
dat <- p0$data[p0$data$isTip, ]

# colour: genus = density (red = hot spot), cold = grey, unclassified ASV = purple
p <- p0 +
  geom_tippoint(data = dat,
                aes(size = log2(n_asv + 1),
                    color = cat,
                    fill = ifelse(cat == "genus", log10_density, NA)),
                shape = 21, stroke = 0.15) +
  scale_color_manual(name = "Tip type",
                     values = c("genus" = "black", "cold" = "grey60",
                                "unclassified_asv" = "#7B3FA0"),
                     breaks = c("genus", "cold", "unclassified_asv"),
                     labels = c("genus clade", "cold clade (no ASV/HOMD)",
                                "unclassified ASV (no genus ref)")) +
  scale_fill_viridis_c(name = "ASV density\nlog10(ASV/ref+1)",
                       option = "A", direction = -1,
                       na.value = "grey92", limits = c(0, 1.2),
                       guide = guide_colorbar(order = 1)) +
  scale_size_continuous(name = "ASV count (log2)", range = c(0.2, 4),
                        guide = guide_legend(order = 2)) +
  theme(legend.position = "right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))

# outer ring 1: number of HOMD references
ring_homd <- data.frame(tip = tip_df$tip, value = tip_df$n_homd_ref)
p <- p + new_scale_fill()
p <- p + geom_fruit(data = ring_homd, geom = geom_tile,
                    mapping = aes(y = tip, fill = value), width = 0.4,
                    offset = 0.06, pwidth = 0.06,
                    axis.params = list(axis = "none")) +
  scale_fill_viridis_c(option = "C", name = "HOMD refs\nin genus", na.value = "grey92",
                       guide = guide_colorbar(order = 3))

# outer ring 2: median ASV-to-nearest-HOMD-reference distance within the genus
ring_med <- data.frame(tip = tip_df$tip, value = tip_df$median_asv_homd_dist)
p <- p + new_scale_fill()
p <- p + geom_fruit(data = ring_med, geom = geom_tile,
                    mapping = aes(y = tip, fill = value), width = 0.4,
                    offset = 0.06, pwidth = 0.06,
                    axis.params = list(axis = "none")) +
  scale_fill_viridis_c(option = "D", name = "median ASV-to-HOMD\ndistance", na.value = "grey92",
                       guide = guide_colorbar(order = 4))

p <- p + ggtitle("ASV insertion density across GG2 tree (genus level)") +
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))

ggsave(file.path(FIGOUT, "fig_v3_insertion_density_tree.png"), p,
       width = 15, height = 15, dpi = 300, limitsize = FALSE)
ggsave(file.path(FIGOUT, "fig_v3_insertion_density_tree.pdf"), p,
       width = 15, height = 15, limitsize = FALSE)
cat("  saved fig_v3_insertion_density_tree.png/pdf\n")

# ============ Figure 2: distance distribution ============
cat("== Figure 2: distance-to-HOMD distribution ==\n")
m <- read.delim(file.path(PIPE_V3, "07_master_asv_table.tsv"), sep = "\t",
                stringsAsFactors = FALSE)
asv_d <- m[, c("asv_id", "nearest_homd_ref_dist", "contamination_flag")]
asv_d$group <- ifelse(asv_d$contamination_flag, "ASV (organelle-flagged)", "ASV (clean)")
asv_d <- asv_d[, c("nearest_homd_ref_dist", "group")]

bg <- read.delim(file.path(OUT, "background_homd_dist.tsv"), sep = "\t",
                 stringsAsFactors = FALSE)
bg$group <- "background refs (non-HOMD)"
bg <- bg[, c("nearest_homd_ref_dist", "group")]

dist_all <- rbind(asv_d, bg)
dist_all$group <- factor(dist_all$group,
                         levels = c("ASV (clean)", "ASV (organelle-flagged)",
                                    "background refs (non-HOMD)"))

p2 <- ggplot(dist_all, aes(x = nearest_homd_ref_dist, fill = group, color = group)) +
  geom_density(alpha = 0.35, linewidth = 0.4) +
  scale_x_log10() +
  scale_fill_manual(values = c("#4E79A7", "#E15759", "grey50")) +
  scale_color_manual(values = c("#4E79A7", "#E15759", "grey50")) +
  labs(x = "distance to nearest HOMD reference (log10)",
       y = "density",
       title = "Distance from ASVs to nearest HOMD reference vs background") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", legend.title = element_blank())

ggsave(file.path(FIGOUT, "fig_v3_asv_homd_distance_distribution.png"), p2,
       width = 9, height = 6, dpi = 300, limitsize = FALSE)
ggsave(file.path(FIGOUT, "fig_v3_asv_homd_distance_distribution.pdf"), p2,
       width = 9, height = 6, limitsize = FALSE)
cat("  saved fig_v3_asv_homd_distance_distribution.png/pdf\n")

# ============ Figure 3: genus-level scatter ============
cat("== Figure 3: genus scatter (ASV density vs HOMD refs) ==\n")
gs3 <- gs[gs$n_asv > 0, ]
gs3$phylum[is.na(gs3$phylum) | gs3$phylum == ""] <- "unknown"
top_phy <- names(sort(table(gs3$phylum), decreasing = TRUE))[1:10]
gs3$phy_group <- ifelse(gs3$phylum %in% top_phy, gs3$phylum, "other")
gs3$phy_group <- factor(gs3$phy_group, levels = c(top_phy, "other"))

p3 <- ggplot(gs3, aes(x = n_homd_ref + 1, y = asv_density,
                      color = phy_group, size = n_asv)) +
  geom_point(alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c(hcl.colors(length(top_phy), "Spectral"), "grey70")) +
  scale_size_continuous(range = c(0.5, 6), name = "ASV count") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  labs(x = "HOMD refs in genus clade (+1, log10)",
       y = "ASV density (ASVs per ref, +1, log10)",
       title = "Per-genus ASV insertion density vs HOMD reference content",
       color = "Phylum") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")

ggsave(file.path(FIGOUT, "fig_v3_genus_density_vs_homd.png"), p3,
       width = 10, height = 7, dpi = 300, limitsize = FALSE)
ggsave(file.path(FIGOUT, "fig_v3_genus_density_vs_homd.pdf"), p3,
       width = 10, height = 7, limitsize = FALSE)
cat("  saved fig_v3_genus_density_vs_homd.png/pdf\n")

cat("Done!\n")
