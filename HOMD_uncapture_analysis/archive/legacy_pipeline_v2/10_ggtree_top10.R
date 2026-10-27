#!/usr/bin/env Rscript
# pipeline_v2: ggtree for the TOP-10 strongest strict candidates.
# - 10 phylotypes, each one a distinct colour; member ASVs share the colour
# - all 10 phylotypes labelled: phylotype_id (genus)
# - heatmap rings: subjects, BioProjects, regions, flagged
suppressMessages({
  library(ggtree)
  library(ape)
  library(tidyverse)
  library(ggnewscale)
})

PIPE_V2 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v2"
REPORT  <- file.path(PIPE_V2, "analysis_report")
FIGOUT  <- file.path(REPORT, "figs")
dir.create(FIGOUT, showWarnings = FALSE)

INK <- "#2B2B2B"
BG  <- "#E8E4E0"

# ---- Load tree ----
tree <- read.tree(file.path(REPORT, "phylotype_subtree_top10.nwk"))
cat(sprintf("Tree tips: %d\n", length(tree$tip.label)))
is_derep <- grepl("^DEREP_", tree$tip.label)

# ---- Load d2p ----
d2p <- read.delim(file.path(REPORT, "derep_to_phylotype_top10.tsv"),
                  sep = "\t", stringsAsFactors = FALSE)
d2p$tip <- as.character(d2p$tip)

# ---- Tip annotation ----
tip_annot <- data.frame(tip = tree$tip.label, stringsAsFactors = FALSE)
tip_annot$type <- ifelse(is_derep, "derep", "gg2")
tip_annot <- merge(tip_annot, d2p, by = "tip", all.x = TRUE)
rownames(tip_annot) <- tip_annot$tip

for (col in c("subjects", "bps", "continents", "regions", "n_flagged", "n_asvs")) {
  tip_annot[[col]] <- as.numeric(tip_annot[[col]])
  tip_annot[[col]][is.na(tip_annot[[col]])] <- 0
}

# ---- Per-phylotype colour (10 distinct hues) ----
strict_mask <- tip_annot$type == "derep" & !is.na(tip_annot$phylotype)
strict_idx <- which(strict_mask)
phy_ids <- unique(tip_annot$phylotype[strict_idx])
cat(sprintf("  phylotypes: %d\n", length(phy_ids)))
pal <- hcl(h = seq(0, 360, length.out = length(phy_ids) + 1)[1:length(phy_ids)],
           c = 80, l = 55)
names(pal) <- phy_ids
tip_annot$pt_color <- BG
tip_annot$pt_color[strict_idx] <- pal[tip_annot$phylotype[strict_idx]]
tip_annot$pt_color <- adjustcolor(tip_annot$pt_color, alpha.f = 0.85)

# ---- Labels: one representative tip per phylotype ----
# pick the first member tip per phylotype to carry the label
first_per_phy <- tapply(tip_annot$tip[strict_idx], tip_annot$phylotype[strict_idx], head, 1)
tip_annot$tip_label <- NA
for (phy_id in names(first_per_phy)) {
  t <- first_per_phy[[phy_id]]
  genus <- tip_annot$rep_gg2_genus[tip_annot$tip == t]
  label <- sprintf("%s (%s)", phy_id, genus)
  tip_annot$tip_label[tip_annot$tip == t] <- label
}

# ---- Base plot (rectangular) ----
p0 <- ggtree(tree, layout = "rectangular", size = 0.3) %<+% tip_annot
p <- p0 +
  geom_tippoint(aes(color = pt_color), size = 1.6) +
  scale_color_identity() +
  theme(legend.position = "none")

# ---- Heatmap rings ----
max_x <- max(node.depth.edgelength(tree))
ring_w <- 0.06 * max_x
off <- c(1.05, 1.14, 1.23, 1.32) * max_x

mk_ring <- function(vals, name) {
  d <- data.frame(vals, row.names = tip_annot$tip)
  colnames(d) <- name
  d
}
ring_subj <- mk_ring(log1p(tip_annot$subjects), "Subjects")
ring_bp   <- mk_ring(log1p(tip_annot$bps), "BioProjects")
ring_reg  <- mk_ring(tip_annot$regions, "Regions")
ring_flag <- mk_ring(tip_annot$n_flagged, "Flagged")

p <- p + new_scale_fill()
p <- gheatmap(p, ring_subj, offset = off[1], width = ring_w, colnames = TRUE,
              colnames_angle = 90, colnames_offset_y = 0.3, font.size = 2.6) +
  scale_fill_viridis_c(option = "C", name = "Subjects", na.value = "#F5F5F5")
p <- p + new_scale_fill()
p <- gheatmap(p, ring_bp, offset = off[2], width = ring_w, colnames = TRUE,
              colnames_angle = 90, colnames_offset_y = 0.3, font.size = 2.6) +
  scale_fill_viridis_c(option = "D", name = "BioProjects", na.value = "#F5F5F5")
p <- p + new_scale_fill()
p <- gheatmap(p, ring_reg, offset = off[3], width = ring_w, colnames = TRUE,
              colnames_angle = 90, colnames_offset_y = 0.3, font.size = 2.6) +
  scale_fill_viridis_c(option = "A", name = "Regions", na.value = "#F5F5F5")
p <- p + new_scale_fill()
p <- gheatmap(p, ring_flag, offset = off[4], width = ring_w, colnames = TRUE,
              colnames_angle = 90, colnames_offset_y = 0.3, font.size = 2.6) +
  scale_fill_viridis_c(option = "magma", name = "Flagged", na.value = "#F5F5F5")

# ---- Save ----
p_final <- p +
  geom_tiplab(aes(label = tip_label), size = 3.2, offset = 0.03 * max_x, color = INK) +
  xlim(0, max(off) + ring_w + 0.12 * max_x) +
  ggtitle("Top-10 strongest candidate novel oral phylotypes on the GG2 tree") +
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
        legend.position = "none")

ggsave(file.path(FIGOUT, "fig_tree_top10.png"), p_final,
       width = 12, height = 9, dpi = 300)
ggsave(file.path(FIGOUT, "fig_tree_top10.pdf"), p_final,
       width = 12, height = 9)
cat("Saved: fig_tree_top10.png / .pdf\n")
cat("Done!\n")
