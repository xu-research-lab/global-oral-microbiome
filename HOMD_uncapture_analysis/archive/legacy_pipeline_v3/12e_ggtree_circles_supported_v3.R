#!/usr/bin/env Rscript
# pipeline_v3: CIRCULAR cluster-circle tree for the SUPPORTED CATALOG (46).
# Each cluster = ONE circle: size = log member ASVs; fill = oral clade.
# Rings: Prevalence, BioProjects, Continents, Regions, Phylum.
suppressMessages({
  library(ggtree)
  library(ape)
  library(tidyverse)
  library(ggnewscale)
  library(ggtreeExtra)
})

PIPE_V3 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
REPORT  <- file.path(PIPE_V3, "analysis_report")
FIGOUT  <- file.path(REPORT, "figs")
dir.create(FIGOUT, showWarnings = FALSE)

INK <- "#2B2B2B"
COL_OUTSIDE <- "#E15759"
COL_IN_HOMD  <- "#4E79A7"

tree <- read.tree(file.path(REPORT, "tree_allcandidates_collapsed.nwk"))
cat(sprintf("Collapsed tree tips: %d\n", length(tree$tip.label)))

anchors <- read.delim(file.path(PIPE_V3, "12d_anchors_supported_v3.tsv"),
                      sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("  anchors: %d\n", nrow(anchors)))

keep <- as.character(anchors$anchor_tip)
missing <- setdiff(keep, tree$tip.label)
cat(sprintf("  anchors missing from tree: %d\n", length(missing)))
drop <- setdiff(tree$tip.label, keep)
tree <- ape::drop.tip(tree, drop)
cat(sprintf("  rep-only tree tips: %d\n", length(tree$tip.label)))

tip_annot <- data.frame(tip = tree$tip.label, stringsAsFactors = FALSE)
tip_annot <- merge(tip_annot, anchors, by.x = "tip", by.y = "anchor_tip", all.x = TRUE)
rownames(tip_annot) <- tip_annot$tip
for (col in c("max_subject_prevalence", "n_bioprojects_total", "n_continents",
              "n_regions", "n_members")) {
  tip_annot[[col]] <- as.numeric(tip_annot[[col]])
  tip_annot[[col]][is.na(tip_annot[[col]])] <- 0
}
tip_annot$oral_clade <- toupper(as.character(tip_annot$oral_clade))
tip_annot$phylum[is.na(tip_annot$phylum) | tip_annot$phylum == ""] <- "unknown"

p0 <- ggtree(tree, layout = "circular", size = 0.3) %<+% tip_annot
rep_data <- p0$data[p0$data$isTip, ]

p <- p0 +
  geom_tippoint(data = rep_data,
                aes(size = log2(n_members + 1), fill = oral_clade),
                shape = 21, color = "black", stroke = 0.3) +
  scale_fill_manual(name = "Oral clade",
                    values = c("TRUE" = COL_IN_HOMD, "FALSE" = COL_OUTSIDE),
                    breaks = c("TRUE", "FALSE"),
                    labels = c("within HOMD genus", "outside HOMD genera"),
                    na.value = "grey70") +
  scale_size_continuous(name = "member ASVs (log2)", range = c(0.3, 3)) +
  theme(legend.position = "right",
        legend.title = element_text(size = 9),
        legend.text  = element_text(size = 8))

ring_prev <- data.frame(tip = tip_annot$tip, value = tip_annot$max_subject_prevalence)
ring_bp   <- data.frame(tip = tip_annot$tip, value = tip_annot$n_bioprojects_total)
ring_cont <- data.frame(tip = tip_annot$tip, value = tip_annot$n_continents)
ring_reg  <- data.frame(tip = tip_annot$tip, value = tip_annot$n_regions)
p <- p + new_scale_fill()
p <- p + geom_fruit(data = ring_prev, geom = geom_tile,
                    mapping = aes(y = tip, fill = value), width = 0.35,
                    offset = 0.05, pwidth = 0.05,
                    axis.params = list(axis = "none")) +
  scale_fill_viridis_c(option = "C", name = "Prevalence", na.value = "grey92",
                       guide = guide_colorbar(order = 1))
p <- p + new_scale_fill()
p <- p + geom_fruit(data = ring_bp, geom = geom_tile,
                    mapping = aes(y = tip, fill = value), width = 0.35,
                    offset = 0.05, pwidth = 0.05,
                    axis.params = list(axis = "none")) +
  scale_fill_viridis_c(option = "D", name = "BioProjects", na.value = "grey92",
                       guide = guide_colorbar(order = 2))
p <- p + new_scale_fill()
p <- p + geom_fruit(data = ring_cont, geom = geom_tile,
                    mapping = aes(y = tip, fill = value), width = 0.35,
                    offset = 0.05, pwidth = 0.05,
                    axis.params = list(axis = "none")) +
  scale_fill_viridis_c(option = "B", name = "Continents", na.value = "grey92",
                       guide = guide_colorbar(order = 3))
p <- p + new_scale_fill()
p <- p + geom_fruit(data = ring_reg, geom = geom_tile,
                    mapping = aes(y = tip, fill = value), width = 0.35,
                    offset = 0.05, pwidth = 0.05,
                    axis.params = list(axis = "none")) +
  scale_fill_viridis_c(option = "A", name = "Regions", na.value = "grey92",
                       guide = guide_colorbar(order = 4))

ring_phy <- data.frame(tip = tip_annot$tip, phylum_group = tip_annot$phylum)
phy_levels <- names(sort(table(ring_phy$phylum_group), decreasing = TRUE))
ring_phy$phylum_group <- factor(ring_phy$phylum_group, levels = phy_levels)
phy_colors <- hcl.colors(max(3, length(phy_levels)), "Spectral")
names(phy_colors) <- phy_levels
phy_colors["unknown"] <- "white"
p <- p + new_scale_fill()
p <- p + geom_fruit(data = ring_phy, geom = geom_tile,
                    mapping = aes(y = tip, fill = phylum_group), width = 0.35,
                    offset = 0.05, pwidth = 0.05,
                    axis.params = list(axis = "none")) +
  scale_fill_manual(values = phy_colors,
                    breaks = setdiff(phy_levels, "unknown"),
                    name = "Phylum",
                    guide = guide_legend(ncol = 3, keyheight = 0.5, keywidth = 0.5,
                                         order = 5,
                                         title.theme = element_text(size = 8),
                                         label.theme = element_text(size = 6)))

p <- p + ggtitle(sprintf("Supported catalog (%d clean clusters) — circle size = member ASVs",
                         length(tree$tip.label))) +
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))

ggsave(file.path(FIGOUT, "fig_tree_circles_supported_v3.png"), p,
       width = 14, height = 14, dpi = 300, limitsize = FALSE)
ggsave(file.path(FIGOUT, "fig_tree_circles_supported_v3.pdf"), p,
       width = 14, height = 14, limitsize = FALSE)
cat("Saved:", file.path(FIGOUT, "fig_tree_circles_supported_v3.png / .pdf"), "\nDone!\n")
