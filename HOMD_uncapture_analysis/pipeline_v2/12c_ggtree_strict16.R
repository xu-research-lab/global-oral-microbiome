#!/usr/bin/env Rscript
# pipeline_v2: ggtree rectangular tree for the 16 STRICT candidates (new scheme)
# - each strict phylotype = one transparent hue (member ASVs share the colour)
# - background tips (GG2 refs, other DEREP) faint gray
# - heatmap rings (subjects, BioProjects, regions)
suppressMessages({
  library(ggtree)
  library(ape)
  library(tidyverse)
  library(ggnewscale)
})

PIPE_V2 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v2"
REPORT  <- file.path(PIPE_V2, "analysis_report")
FIGOUT  <- file.path(REPORT, "figs_v2")
dir.create(FIGOUT, showWarnings = FALSE)

INK <- "#2B2B2B"
BG  <- "#EDEAE6"
args <- commandArgs(trailingOnly = TRUE)
SUFFIX <- if (length(args) >= 1) args[1] else "headline53"

tree <- read.tree(file.path(REPORT, paste0("phylotype_subtree_", SUFFIX, ".nwk")))
cat(sprintf("Tree tips: %d\n", length(tree$tip.label)))
tip_labels <- tree$tip.label
is_derep <- grepl("^DEREP_", tip_labels)

d2p <- read.delim(file.path(REPORT, paste0("derep_to_phylotype_", SUFFIX, ".tsv")),
                  sep = "\t", stringsAsFactors = FALSE)
d2p$tip <- as.character(d2p$tip)
cat(sprintf("  d2p rows: %d\n", nrow(d2p)))

tip_annot <- data.frame(tip = tree$tip.label, stringsAsFactors = FALSE)
tip_annot$type <- ifelse(is_derep, "derep", "gg2")
tip_annot <- merge(tip_annot, d2p, by = "tip", all.x = TRUE)
rownames(tip_annot) <- tip_annot$tip

for (col in c("subjects", "bps", "regions", "n_asvs")) {
  tip_annot[[col]] <- as.numeric(tip_annot[[col]])
  tip_annot[[col]][is.na(tip_annot[[col]])] <- 0
}

strict_mask <- tip_annot$type == "derep" & !is.na(tip_annot$phylotype)
strict_idx <- which(strict_mask)
phy_ids <- unique(tip_annot$phylotype[strict_idx])
cat(sprintf("  strict phylotypes on tree: %d\n", length(phy_ids)))
pal <- hcl(h = seq(0, 360, length.out = length(phy_ids) + 1)[1:length(phy_ids)],
           c = 75, l = 55)
names(pal) <- phy_ids
tip_annot$pt_color <- BG
tip_annot$pt_color[strict_idx] <- pal[tip_annot$phylotype[strict_idx]]
tip_annot$pt_color <- adjustcolor(tip_annot$pt_color, alpha.f = 0.75)

tip_annot$pt_size <- ifelse(strict_mask, 0.9, 0.25)

# labels: species name when available, otherwise genus, otherwise phylotype id
lab <- ifelse(!is.na(tip_annot$rep_species) & tip_annot$rep_species != "",
              tip_annot$rep_species,
              ifelse(!is.na(tip_annot$rep_genus) & tip_annot$rep_genus != "",
                     tip_annot$rep_genus, tip_annot$phylotype))
tip_annot$tip_label <- ifelse(strict_mask, lab, NA)
tip_annot$tip_label[is.na(tip_annot$tip_label)] <- NA

p0 <- ggtree(tree, layout = "rectangular", size = 0.12) %<+% tip_annot
p <- p0 +
  geom_tippoint(aes(color = pt_color, size = pt_size)) +
  scale_color_identity() +
  scale_size_identity() +
  theme(legend.position = "none")

max_x <- max(node.depth.edgelength(tree))
ring_w <- 0.05 * max_x
off <- c(1.02, 1.09, 1.16) * max_x

mk_ring <- function(vals, name) {
  d <- data.frame(vals, row.names = tip_annot$tip)
  colnames(d) <- name
  d
}

p <- p + new_scale_fill()
p <- gheatmap(p, mk_ring(log1p(tip_annot$subjects), "Subjects"), offset = off[1],
              width = ring_w, colnames = TRUE, colnames_angle = 90,
              colnames_offset_y = 0.3, font.size = 2.5) +
  scale_fill_viridis_c(option = "C", name = "Subjects", na.value = "#F5F5F5")
p <- p + new_scale_fill()
p <- gheatmap(p, mk_ring(log1p(tip_annot$bps), "BioProjects"), offset = off[2],
              width = ring_w, colnames = TRUE, colnames_angle = 90,
              colnames_offset_y = 0.3, font.size = 2.5) +
  scale_fill_viridis_c(option = "D", name = "BioProjects", na.value = "#F5F5F5")
p <- p + new_scale_fill()
p <- gheatmap(p, mk_ring(tip_annot$regions, "Regions"), offset = off[3],
              width = ring_w, colnames = TRUE, colnames_angle = 90,
              colnames_offset_y = 0.3, font.size = 2.5) +
  scale_fill_viridis_c(option = "A", name = "Regions", na.value = "#F5F5F5")

p_final <- p +
  geom_tiplab(aes(label = tip_label), size = 2.4, offset = 0.02 * max_x, color = INK) +
  xlim(0, max(off) + ring_w + 0.06 * max_x) +
  ggtitle("16 strict candidates (GG2 2024.09 closed-reference) on the GG2 tree") +
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
        legend.position = "right")

ggsave(file.path(FIGOUT, paste0("fig_tree_", SUFFIX, ".png")), p_final,
       width = 10, height = 30, dpi = 300, limitsize = FALSE)
cat("Saved:", paste0("fig_tree_", SUFFIX, ".png"), "\nDone!\n")
