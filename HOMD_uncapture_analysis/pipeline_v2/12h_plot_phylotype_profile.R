#!/usr/bin/env Rscript
# pipeline_v2: per-phylotype profile for ALL phylotypes (26,311), horizontal.
# y = phylotype ranked by max_subject_prevalence descending.
# Panels (left -> right): No. ASVs, GG2 id max, GG2 id min, HOMD id max,
# BioProjects, Continents, Subject prevalence (bars), Phylum,
# Species in HOMD, Contamination flag.
suppressPackageStartupMessages({
  library(ggplot2)
  library(aplot)
})

PIPE_V2 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v2"
FIGOUT  <- file.path(PIPE_V2, "analysis_report", "figs_v2")
dir.create(FIGOUT, showWarnings = FALSE, recursive = TRUE)

df <- read.delim(file.path(PIPE_V2, "12g_phylotype_profile.tsv"),
                 sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
message("[12h] loaded ", nrow(df), " phylotypes")

df$rank_f <- factor(df$rank, levels = rev(df$rank))   # highest prevalence on top
df$species_in_homd <- ifelse(df$species_in_homd, "in HOMD", "not in HOMD")
df$contam <- ifelse(df$clean_phylotype, "clean", "flagged")

theme_h <- theme_minimal(base_size = 7) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 5, color = "#222222"),
        axis.ticks.x = element_blank(), axis.title = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = margin(2, 4, 2, 2),
        legend.position = "bottom", legend.box = "horizontal",
        legend.title = element_text(size = 5.5),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.16, "cm"))

mk_tile <- function(d, fill_var, scale_fun, show_y = FALSE) {
  p <- ggplot(d, aes(y = rank_f, x = track, fill = !!sym(fill_var))) +
    geom_tile(width = 0.88, height = 0.78) +
    scale_fun + theme_h +
    theme(axis.text.y = if (show_y) element_text(size = 3.5) else element_blank())
  p
}

# ---- panel 1: No. ASVs ----
p_asv <- ggplot(df, aes(y = rank_f, x = "No. ASVs", fill = n_asvs)) +
  geom_tile(width = 0.88, height = 0.78) +
  scale_fill_viridis_c(option = "C", name = "No. ASVs") + theme_h

# ---- GG2 max / min ----
p_gg2max <- ggplot(df, aes(y = rank_f, x = "GG2 id max", fill = phylo_max_gg2_id)) +
  geom_tile(width = 0.88, height = 0.78) +
  scale_fill_viridis_c(option = "D", limits = c(0, 100), name = "GG2 id max") + theme_h
p_gg2min <- ggplot(df, aes(y = rank_f, x = "GG2 id min", fill = phylo_min_gg2_id)) +
  geom_tile(width = 0.88, height = 0.78) +
  scale_fill_viridis_c(option = "D", limits = c(0, 100), name = "GG2 id min") + theme_h

# ---- HOMD max ----
p_homd <- ggplot(df, aes(y = rank_f, x = "HOMD id max", fill = phylo_max_homd_id)) +
  geom_tile(width = 0.88, height = 0.78) +
  scale_fill_viridis_c(option = "B", limits = c(0, 100), name = "HOMD id max") + theme_h

# ---- BioProjects / Continents ----
p_bp <- ggplot(df, aes(y = rank_f, x = "BioProjects", fill = n_bioprojects_total)) +
  geom_tile(width = 0.88, height = 0.78) +
  scale_fill_viridis_c(option = "C", name = "BioProjects") + theme_h
p_cont <- ggplot(df, aes(y = rank_f, x = "Continents", fill = n_continents)) +
  geom_tile(width = 0.88, height = 0.78) +
  scale_fill_viridis_c(option = "B", name = "Continents") + theme_h

# ---- Subject prevalence (bars, far right) ----
p_prev <- ggplot(df, aes(y = rank_f, x = max_subject_prevalence)) +
  geom_col(fill = "#2F5D62", width = 0.78) +
  labs(x = "Subject prevalence") +
  theme_h +
  theme(axis.title.x = element_text(size = 6), legend.position = "none")

# ---- Phylum (categorical) ----
all_phy <- unique(df$phylum)
phy_colors <- hcl.colors(max(3, length(all_phy)), "Spectral")[seq_along(all_phy)]
names(phy_colors) <- all_phy
phy_colors[""] <- "#BBBBBB"
p_phylum <- ggplot(df, aes(y = rank_f, x = "Phylum", fill = phylum)) +
  geom_tile(width = 0.88, height = 0.78) +
  scale_fill_manual(values = phy_colors, name = "Phylum", drop = FALSE) + theme_h

# ---- Species in HOMD / contamination ----
p_homd_sp <- ggplot(df, aes(y = rank_f, x = "Species in HOMD", fill = species_in_homd)) +
  geom_tile(width = 0.88, height = 0.78) +
  scale_fill_manual(values = c("in HOMD" = "#4E79A7", "not in HOMD" = "#E15759"),
                    name = "Species in HOMD") + theme_h
p_contam <- ggplot(df, aes(y = rank_f, x = "Contamination", fill = contam)) +
  geom_tile(width = 0.88, height = 0.78) +
  scale_fill_manual(values = c("clean" = "#59A14F", "flagged" = "#8C8C8C"),
                    name = "Contamination") + theme_h

combined <- p_asv %>%
  insert_right(p_gg2max, width = 0.15) %>%
  insert_right(p_gg2min, width = 0.15) %>%
  insert_right(p_homd, width = 0.15) %>%
  insert_right(p_bp, width = 0.15) %>%
  insert_right(p_cont, width = 0.15) %>%
  insert_right(p_phylum, width = 0.40) %>%
  insert_right(p_homd_sp, width = 0.30) %>%
  insert_right(p_contam, width = 0.30) %>%
  insert_right(p_prev, width = 0.40)
combined$width <- c(0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.40, 0.30, 0.30, 0.40)

h_in <- max(8, nrow(df) / 350)   # ~75 in for 26k rows; vector PDF stays zoomable
out_base <- file.path(FIGOUT, "fig_v2_all_phylotype_profile")
ggsave(paste0(out_base, ".pdf"), combined, width = 9, height = h_in,
       device = cairo_pdf, limitsize = FALSE)
ggsave(paste0(out_base, ".png"), combined, width = 9, height = min(h_in, 60),
       dpi = 150, limitsize = FALSE)
cat("Saved:", paste0(out_base, ".pdf / .png"), "\nDone!\n")
