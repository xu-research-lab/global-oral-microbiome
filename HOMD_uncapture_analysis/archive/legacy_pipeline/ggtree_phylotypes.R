#!/usr/bin/env Rscript
# ggtree circular tree with heatmap rings
suppressMessages({
  library(ggtree)
  library(ape)
  library(tidyverse)
  library(ggnewscale)
})

PIPE <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline"
ITOL <- file.path(PIPE, "analysis_report/itol")
OUT  <- file.path(PIPE, "analysis_report")

# ---- Load tree ----
tree <- read.tree(file.path(ITOL, "phylotype_subtree.nwk"))

# ---- Load phylotype metadata ----
phy <- read.delim(file.path(ITOL, "phylotype_info.tsv"), sep="\t", stringsAsFactors=FALSE)

# Get contamination info per phylotype
phy_sepp <- read.delim(file.path(PIPE, "06_phylotypes_sepp/phylotypes_candidate_novel_updated.tsv"), sep="\t", stringsAsFactors=FALSE)
flagged <- read.delim(file.path(PIPE, "05_integrated_candidates/candidates_flagged_updated.tsv"), sep="\t", stringsAsFactors=FALSE)

for (i in seq_len(nrow(phy))) {
  idx <- match(phy$id[i], phy_sepp$phylotype_id)
  if (!is.na(idx)) {
    members <- trimws(strsplit(as.character(phy_sepp$member_derep_ids[idx]), ";")[[1]])
    sub <- flagged[flagged$derep_id %in% members, ]
    n97 <- sum(sub$cluster97_contaminant == "true", na.rm=TRUE)
    ng  <- sum(sub$genus_contaminant == "true", na.rm=TRUE)
    phy$n97[i] <- n97
    phy$genusf[i] <- ng
    if (n97 > 0 && ng > 0) phy$concern[i] <- "BOTH"
    else if (n97 > 0) phy$concern[i] <- "cluster97"
    else if (ng > 0) phy$concern[i] <- "genus"
    else phy$concern[i] <- "clean"
  }
}
phy$n97[is.na(phy$n97)] <- 0
phy$genusf[is.na(phy$genusf)] <- 0
phy$concern[is.na(phy$concern)] <- "clean"

# ---- Prepare tip metadata: contamination + species annotation ----
tip_labels <- tree$tip.label
is_phylo <- grepl("^EDGE_", tip_labels)

# ---- Build heatmap data matrices (row = tip, col = value) ----
# Only for phylotype tips — use original tip labels (matching tree) as row names
phy_orig <- tip_labels[is_phylo]

# Map from original label back to phy$id
phy_map_idx <- match(phy_orig, phy$id)

# Matrix 1: Subjects (log-scale)
subj_mat <- matrix(NA, nrow=length(phy_orig), ncol=1, dimnames=list(phy_orig, "Subjects"))
subj_mat[,1] <- phy$subjects[phy_map_idx]

# Matrix 2: BioProjects
bps_mat <- matrix(NA, nrow=length(phy_orig), ncol=1, dimnames=list(phy_orig, "BPs"))
bps_mat[,1] <- phy$bps[phy_map_idx]

# Matrix 3: Continent presence
all_conts <- c("N.America","Europe","Asia","Oceania","S.America","Africa")
cont_mat <- matrix(0, nrow=length(phy_orig), ncol=length(all_conts),
                   dimnames=list(phy_orig, all_conts))
for (i in seq_along(phy_orig)) {
  idx <- phy_map_idx[i]
  if (!is.na(idx) && !is.na(phy$continents[idx]) && phy$continents[idx] != "") {
    # Continent field is Python-set literal: {'Asia', 'Europe', 'North America'}
    # Strip { } ' then split by comma
    conts_clean <- gsub("[{}']", "", as.character(phy$continents[idx]))
    conts <- trimws(strsplit(conts_clean, ",")[[1]])
    conts <- conts[conts != ""]
    for (ct in conts) {
      ct_short <- switch(ct,
                         "North America"="N.America", "South America"="S.America",
                         ct)
      if (ct_short %in% all_conts) cont_mat[phy_orig[i], ct_short] <- 1
    }
  }
}

# ---- Build tree ----
# Build metadata data frame with species annotation
d <- data.frame(label=tree$tip.label,
                is_phy=is_phylo,
                concern=NA_character_,
                species_label=NA_character_,
                stringsAsFactors=FALSE)
for (i in seq_len(nrow(d))) {
  if (is_phylo[i]) {
    idx2 <- match(tip_labels[i], phy$id)
    if (!is.na(idx2)) {
      d$concern[i] <- phy$concern[idx2]
      nm <- phy$gg2_name[idx2]
      if (is.na(nm) || nm == "") nm <- phy$gg2_ref[idx2]
      if (is.na(nm) || nm == "") nm <- "Unknown"
      d$species_label[i] <- paste0(tip_labels[i], "  ", nm, "  [", phy$asvs[idx2], "ASV]")
    }
  }
}

cat("=== Debug: first 30 phylotype labels ===\n")
phy_labels <- d$species_label[is_phylo]
for (j in seq_len(min(30, length(phy_labels)))) {
  cat(sprintf("  %s\n", phy_labels[j]))
}
cat(sprintf("Total phylotype tips with labels: %d\n", sum(!is.na(d$species_label))))

# Base tree
p <- ggtree(tree, layout="circular", size=0.2, ladderize=FALSE)

# Suppress default tip labels; use geom_tiplab2 with species annotation
p <- p %<+% d

# Color tip points by contamination status (instead of novelty ring)
p <- p + geom_tippoint(aes(color=concern, subset=is_phylo), size=1.3) +
  scale_color_manual(values=c(clean="#27ae60", cluster97="#f39c12", genus="#e67e22", BOTH="#e74c3c"),
                     name="Contamination",
                     labels=c(clean="Clean", cluster97="97% cluster", genus="Genus", BOTH="Both"),
                     na.translate=FALSE) +
  new_scale_color()

# ---- Add heatmap rings ----
# Ring 1: Subjects (log10 scale, red gradient)
p1 <- gheatmap(p, subj_mat, offset=0.06, width=0.06,
               colnames=TRUE, colnames_position="top",
               colnames_angle=90, colnames_offset_y=1.5,
               hjust=0, font.size=2.5, color=NA) +
  scale_fill_gradient(name="Subjects", na.value=NA,
                      low="#fee5d9", high="#cb181d",
                      trans="log10",
                      breaks=c(100, 500, 2000),
                      labels=c("100", "500", "2K"))

# Debug ring 1 data
cat(sprintf("\n--- Ring 1 (Subjects) data range: %.0f - %.0f, NA: %d ---\n",
            min(subj_mat, na.rm=TRUE), max(subj_mat, na.rm=TRUE), sum(is.na(subj_mat))))

# Ring 2: BioProjects (blue gradient)
p2 <- p1 + new_scale_fill()
p2 <- gheatmap(p2, bps_mat, offset=0.14, width=0.06,
               colnames=TRUE, colnames_position="top",
               colnames_angle=90, colnames_offset_y=1.5,
               hjust=0, font.size=2.5, color=NA) +
  scale_fill_gradient(name="BioProjects", na.value=NA,
                      low="#deebf7", high="#08519c",
                      breaks=c(5, 15, 30),
                      labels=c("5", "15", "30"))

cat(sprintf("--- Ring 2 (BPs) data range: %.0f - %.0f, NA: %d ---\n",
            min(bps_mat, na.rm=TRUE), max(bps_mat, na.rm=TRUE), sum(is.na(bps_mat))))

# Ring 3: Continents (binary presence, discrete colors)
# Convert 0/1 matrix to factor matrix for gheatmap
cont_mat_f <- matrix(as.character(cont_mat), nrow=nrow(cont_mat), ncol=ncol(cont_mat),
                     dimnames=dimnames(cont_mat))
cont_mat_f[cont_mat_f == "0"] <- "Absent"
cont_mat_f[cont_mat_f == "1"] <- "Present"

p3 <- p2 + new_scale_fill()
p3 <- gheatmap(p3, cont_mat_f, offset=0.22, width=0.06,
               colnames=TRUE, colnames_position="top",
               colnames_angle=90, colnames_offset_y=1.5,
               hjust=0, font.size=2, color=NA) +
  scale_fill_manual(name="Continents", na.value="grey90",
                    values=c(Absent="white", Present="#6a51a3"))

cat(sprintf("--- Ring 3 (Continents) ones: %d / %d cells ---\n",
            sum(cont_mat, na.rm=TRUE), length(cont_mat)))

# ---- Tip labels: tight to tips (offset=0), added LAST to render on top ----
# Use tree$tip.label directly (already has species annotation from d$species_label)
p4 <- p3 + geom_tiplab2(aes(label=species_label, subset=is_phylo),
                        size=1.6, offset=0, align=FALSE, linetype=NA, hjust=0)

# Title
p4 <- p4 + ggtitle("Candidate Novel Oral Phylotypes on Greengenes2 Backbone",
                   subtitle="Tip color: contamination flag | Ring 1: Subjects (red) | Ring 2: BioProjects (blue) | Ring 3: Continents (purple)") +
  theme(plot.title=element_text(hjust=0.5, size=13, face="bold"),
        plot.subtitle=element_text(hjust=0.5, size=7, color="grey50"),
        legend.position="right",
        legend.text=element_text(size=7),
        legend.title=element_text(size=8))

# Save
ggsave(file.path(OUT, "phylotypes_ggtree_circular.pdf"), p4,
       width=18, height=16, device=cairo_pdf, limitsize=FALSE)
ggsave(file.path(OUT, "phylotypes_ggtree_circular.png"), p4,
       width=18, height=16, dpi=250, bg="white", limitsize=FALSE)

cat("Saved:", file.path(OUT, "phylotypes_ggtree_circular.pdf/png"), "\n")
cat("Tips:", length(tree$tip.label), "| Phylotypes:", sum(is_phylo), "\n")
