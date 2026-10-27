#!/usr/bin/env Rscript
# ggtree rectangular tree with actual DEREP ASV tips + GG2 sister refs
suppressMessages({
  library(ggtree)
  library(ape)
  library(tidyverse)
  library(ggnewscale)
})

PIPE <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline"
ITOL <- file.path(PIPE, "analysis_report/itol")
OUT  <- file.path(PIPE, "analysis_report")

# ---- Load corrected tree (actual ASV tips + GG2 refs) ----
tree <- read.tree(file.path(ITOL, "phylotype_subtree_corrected.nwk"))

cat(sprintf("Tree tips: %d\n", length(tree$tip.label)))

# ---- Identify tip types ----
tip_labels <- tree$tip.label
is_derep <- grepl("^DEREP_", tip_labels)
is_gg2   <- !is_derep

cat(sprintf("  DEREP ASV tips: %d\n", sum(is_derep)))
cat(sprintf("  GG2 ref tips: %d\n", sum(is_gg2)))

# ---- Load DEREP → phylotype mapping ----
d2p <- read.delim(file.path(ITOL, "derep_to_phylotype.tsv"), sep="\t", stringsAsFactors=FALSE)

# ---- Load contamination data (per phylotype) ----
phy_sepp <- read.delim(file.path(PIPE, "06_phylotypes_sepp/phylotypes_candidate_novel_updated.tsv"), sep="\t", stringsAsFactors=FALSE)
flagged  <- read.delim(file.path(PIPE, "05_integrated_candidates/candidates_flagged_updated.tsv"), sep="\t", stringsAsFactors=FALSE)

# Build phylotype → concern mapping
compute_concern <- function(phylotype_id) {
  idx <- match(phylotype_id, phy_sepp$phylotype_id)
  if (is.na(idx)) return("clean")
  members <- trimws(strsplit(as.character(phy_sepp$member_derep_ids[idx]), ";")[[1]])
  sub <- flagged[flagged$derep_id %in% members, ]
  n97 <- sum(sub$cluster97_contaminant == "true", na.rm=TRUE)
  ng  <- sum(sub$genus_contaminant == "true", na.rm=TRUE)
  if (n97 > 0 && ng > 0) return("BOTH")
  if (n97 > 0) return("cluster97")
  if (ng > 0) return("genus")
  return("clean")
}

# Apply contamination per phylotype from d2p
d2p$concern <- sapply(d2p$phylotype, compute_concern, USE.NAMES=FALSE)

# Create a named vector for quick lookup: DEREP → concern
derep_concern <- setNames(d2p$concern, d2p$tip)
derep_subjects <- setNames(d2p$subjects, d2p$tip)
derep_bps <- setNames(d2p$bps, d2p$tip)
derep_novelty <- setNames(d2p$novelty, d2p$tip)
derep_gg2_name <- setNames(d2p$gg2_name, d2p$tip)
derep_phylo <- setNames(d2p$phylotype, d2p$tip)

# ---- Continent presence per DEREP ----
all_conts <- c("N.America","Europe","Asia","Oceania","S.America","Africa")
continent_set <- function(cont_str) {
  if (is.na(cont_str) || cont_str == "") return(character(0))
  clean <- gsub("[{}']", "", cont_str)
  parts <- trimws(strsplit(clean, ",")[[1]])
  parts <- parts[parts != ""]
  result <- c()
  for (ct in parts) {
    ct_short <- switch(ct,
      "North America"="N.America", "South America"="S.America", ct)
    if (ct_short %in% all_conts) result <- c(result, ct_short)
  }
  unique(result)
}

# ---- Build metadata data frame (row order = tree tip order) ----
d <- data.frame(
  label = tip_labels,
  is_der = is_derep,
  concern = NA_character_,
  subjects = NA_real_,
  bps = NA_real_,
  novelty = NA_character_,
  gg2_name = NA_character_,
  phylotype = NA_character_,
  stringsAsFactors = FALSE
)
rownames(d) <- tip_labels

for (i in seq_along(tip_labels)) {
  if (is_derep[i]) {
    nm <- tip_labels[i]
    d$concern[i] <- ifelse(nm %in% names(derep_concern), derep_concern[nm], "clean")
    d$subjects[i] <- ifelse(nm %in% names(derep_subjects), derep_subjects[nm], NA)
    d$bps[i] <- ifelse(nm %in% names(derep_bps), derep_bps[nm], NA)
    d$novelty[i] <- ifelse(nm %in% names(derep_novelty), derep_novelty[nm], "ref")
    d$gg2_name[i] <- ifelse(nm %in% names(derep_gg2_name), derep_gg2_name[nm], "")
    d$phylotype[i] <- ifelse(nm %in% names(derep_phylo), derep_phylo[nm], "")
  } else {
    # GG2 ref tips: short label from name (remove path-like prefix)
    nm <- tip_labels[i]
    if (nchar(nm) > 35) nm <- paste0(substr(nm, 1, 32), "...")
    d$gg2_name[i] <- nm
    d$novelty[i] <- "ref"
  }
}

# ---- Heatmap: Continents matrix (tips x continents) ----
cont_mat <- matrix(FALSE, nrow=length(tip_labels), ncol=length(all_conts),
                   dimnames=list(tip_labels, all_conts))
for (i in seq_along(tip_labels)) {
  if (is_derep[i]) {
    nm <- tip_labels[i]
    if (nm %in% d2p$tip) {
      cont_str <- d2p$continents[d2p$tip == nm]
      if (length(cont_str) > 0) {
        conts <- continent_set(cont_str[1])
        for (ct in conts) {
          if (ct %in% all_conts) cont_mat[nm, ct] <- TRUE
        }
      }
    }
  }
}
cont_mat_char <- matrix(ifelse(cont_mat, "Present", "Absent"),
                        nrow=nrow(cont_mat), ncol=ncol(cont_mat),
                        dimnames=dimnames(cont_mat))

# ---- Summary stats ----
cat(sprintf("\nMetadata:\n  DEREP with concern: %d\n", sum(!is.na(d$concern[d$is_der]))))
cat(sprintf("  DEREP with subjects: %d\n", sum(!is.na(d$subjects[d$is_der]))))
cat(sprintf("  Continents matrix: %d presents / %d cells\n",
            sum(cont_mat, na.rm=TRUE), prod(dim(cont_mat))))

# ============================================================
# Build rectangular tree
# ============================================================
p <- ggtree(tree, layout="rectangular", size=0.15, ladderize=FALSE)
p <- p %<+% d

# GG2 ref tips: show as small grey points with labels
p <- p + geom_tippoint(aes(subset=is_gg2), size=0.5, color="grey60")

# GG2 ref labels (important ones)
p <- p + geom_tiplab(aes(label=gg2_name, subset=is_gg2),
                     size=1.0, offset=0.001, color="grey40", linetype=NA)

# DEREP tips: colored by contamination, sized to show ASV clustering
p <- p + geom_tippoint(aes(color=concern, subset=is_der), size=1.2) +
  scale_color_manual(values=c(clean="#27ae60", cluster97="#f39c12", genus="#e67e22", BOTH="#e74c3c"),
                     name="Contamination",
                     labels=c(clean="Clean", cluster97="97% cluster", genus="Genus", BOTH="Both"),
                     na.translate=FALSE) +
  new_scale_color()

# ---- Heatmap columns ----
# Column 1: Subjects (per DEREP, log scale)
col1 <- d[, "subjects", drop=FALSE]
p1 <- gheatmap(p, col1, offset=0.005, width=0.03,
               colnames=TRUE, colnames_position="top",
               colnames_angle=0, colnames_offset_y=0.5,
               font.size=3, color=NA) +
  scale_fill_gradient(name="Subjects", na.value="#f0f0f0",
                      low="#fee5d9", high="#cb181d",
                      trans="log10",
                      breaks=c(100, 500, 2000),
                      labels=c("100", "500", "2K"))

# Column 2: BioProjects
col2 <- d[, "bps", drop=FALSE]
p2 <- p1 + new_scale_fill()
p2 <- gheatmap(p2, col2, offset=0.04, width=0.03,
               colnames=TRUE, colnames_position="top",
               colnames_angle=0, colnames_offset_y=0.5,
               font.size=3, color=NA) +
  scale_fill_gradient(name="BioProjects", na.value="#f0f0f0",
                      low="#deebf7", high="#08519c",
                      breaks=c(5, 15, 30),
                      labels=c("5", "15", "30"))

# Column 3-8: Continents (discrete)
p3 <- p2 + new_scale_fill()
p3 <- gheatmap(p3, cont_mat_char, offset=0.08, width=0.08,
               colnames=TRUE, colnames_position="top",
               colnames_angle=45, colnames_offset_y=1,
               font.size=2, color="grey90") +
  scale_fill_manual(name="Continents", na.value="grey90",
                    values=c(Absent="white", Present="#6a51a3"))

# Title
p3 <- p3 + ggtitle("Candidate Novel Oral Phylotypes on Greengenes2 Backbone",
                   subtitle="DEREP tips (points): contamination color | Heatmap: Subjects, BioProjects, Continents | Grey tips: GG2 references") +
  theme(plot.title=element_text(hjust=0.5, size=13, face="bold"),
        plot.subtitle=element_text(hjust=0.5, size=7, color="grey50"),
        legend.position="right",
        legend.text=element_text(size=7),
        legend.title=element_text(size=8))

# Save
ggsave(file.path(OUT, "phylotypes_ggtree_rect.pdf"), p3,
       width=24, height=30, device=cairo_pdf, limitsize=FALSE)
ggsave(file.path(OUT, "phylotypes_ggtree_rect.png"), p3,
       width=24, height=30, dpi=200, bg="white", limitsize=FALSE)

cat(sprintf("\nSaved: %s\n", file.path(OUT, "phylotypes_ggtree_rect.pdf/png")))
cat(sprintf("Tips: %d | DEREP (ASV): %d | GG2 ref: %d\n",
            length(tip_labels), sum(is_derep), sum(is_gg2)))
