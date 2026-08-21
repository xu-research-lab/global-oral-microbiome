#!/usr/bin/env Rscript
# pipeline_v2: ggtree rectangular tree for the strict candidate phylotypes
# - each strict phylotype = one transparent hue (member ASVs share the colour)
# - background tips (GG2 refs, non-candidate DEREP) faint gray
# - heatmap rings (subjects, BioProjects, regions, flagged) with expanded xlim
# - top-N labels: phylotype_id (genus)
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

INK  <- "#2B2B2B"
BG   <- "#EDEAE6"   # background tip colour (GG2 refs / other DEREP)

# ---- Load tree ----
# suffix auto-derived from the strict-candidate count (09_build_strict_tree.py)
SUFFIX <- "strict71"
tree <- read.tree(file.path(REPORT, paste0("phylotype_subtree_", SUFFIX, ".nwk")))
cat(sprintf("Tree tips: %d\n", length(tree$tip.label)))
tip_labels <- tree$tip.label
is_derep <- grepl("^DEREP_", tip_labels)

# ---- Load DEREP → phylotype metadata (with genus) ----
d2p <- read.delim(file.path(REPORT, paste0("derep_to_phylotype_", SUFFIX, ".tsv")),
                  sep = "\t", stringsAsFactors = FALSE)
d2p$tip <- as.character(d2p$tip)
cat(sprintf("  d2p rows: %d\n", nrow(d2p)))

# ---- Tip annotation ----
tip_annot <- data.frame(tip = tree$tip.label, stringsAsFactors = FALSE)
tip_annot$type <- ifelse(is_derep, "derep", "gg2")
tip_annot <- merge(tip_annot, d2p, by = "tip", all.x = TRUE)
rownames(tip_annot) <- tip_annot$tip

# numeric columns
for (col in c("subjects", "bps", "continents", "regions", "n_flagged", "n_asvs")) {
  tip_annot[[col]] <- as.numeric(tip_annot[[col]])
  tip_annot[[col]][is.na(tip_annot[[col]])] <- 0
}

# ---- Per-phylotype transparent colour ----
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

# ---- Novelty (for reference, not colouring) ----
tip_annot$novelty[tip_annot$type == "gg2"] <- "GG2 ref"
tip_annot$novelty[is.na(tip_annot$novelty)] <- "other DEREP"

# point size: strict members bigger, background tiny
tip_annot$pt_size <- ifelse(strict_mask, 0.9, 0.25)

# ---- Top-N labels (phylotype_id + genus) — set before attaching data ----
phy_support <- d2p[!duplicated(d2p$phylotype), c("phylotype", "subjects", "rep_gg2_genus")]
phy_support <- phy_support[order(-phy_support$subjects), ]
top_ids <- head(phy_support$phylotype, 25)
tip_annot$tip_label <- ifelse(tip_annot$phylotype %in% top_ids,
                          sprintf("%s (%s)", tip_annot$phylotype, tip_annot$rep_gg2_genus),
                          NA)
tip_annot$tip_label[grepl("^\\S+ \\(NA\\)$", tip_annot$tip_label)] <- NA

# ---- Base plot (rectangular) ----
p0 <- ggtree(tree, layout = "rectangular", size = 0.12) %<+% tip_annot
p <- p0 +
  geom_tippoint(aes(color = pt_color, size = pt_size)) +
  scale_color_identity() +
  scale_size_identity() +
  theme(legend.position = "none")

# ---- Heatmap rings (rectangular) ----
max_x <- max(node.depth.edgelength(tree))
ring_w  <- 0.05 * max_x
off <- c(1.02, 1.09, 1.16, 1.23) * max_x

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
              colnames_angle = 90, colnames_offset_y = 0.3, font.size = 2.5) +
  scale_fill_viridis_c(option = "C", name = "Subjects", na.value = "#F5F5F5")

p <- p + new_scale_fill()
p <- gheatmap(p, ring_bp, offset = off[2], width = ring_w, colnames = TRUE,
              colnames_angle = 90, colnames_offset_y = 0.3, font.size = 2.5) +
  scale_fill_viridis_c(option = "D", name = "BioProjects", na.value = "#F5F5F5")

p <- p + new_scale_fill()
p <- gheatmap(p, ring_reg, offset = off[3], width = ring_w, colnames = TRUE,
              colnames_angle = 90, colnames_offset_y = 0.3, font.size = 2.5) +
  scale_fill_viridis_c(option = "A", name = "Regions", na.value = "#F5F5F5")

p <- p + new_scale_fill()
p <- gheatmap(p, ring_flag, offset = off[4], width = ring_w, colnames = TRUE,
              colnames_angle = 90, colnames_offset_y = 0.3, font.size = 2.5) +
  scale_fill_viridis_c(option = "magma", name = "Flagged", na.value = "#F5F5F5")

# ---- Save ----
p_final <- p +
  geom_tiplab(aes(label = tip_label), size = 2.4, offset = 0.02 * max_x, color = INK) +
  xlim(0, max(off) + ring_w + 0.06 * max_x) +
  ggtitle("71 strict-evidence candidate novel oral phylotypes on the GG2 tree") +
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
        legend.position = "right")

ggsave(file.path(FIGOUT, paste0("fig_tree_", SUFFIX, ".png")), p_final,
       width = 10, height = 30, dpi = 300, limitsize = FALSE)
ggsave(file.path(FIGOUT, paste0("fig_tree_", SUFFIX, ".pdf")), p_final,
       width = 10, height = 30, limitsize = FALSE)
cat("Saved:", paste0("fig_tree_", SUFFIX, ".png / .pdf"), "\n")
cat("Done!\n")
