#!/usr/bin/env Rscript
# PREVIEW -- the identity-band composition panel on its own, at its own size.
#
# The panel itself is built by 13z2_panel_homd_band_composition.R, which panel B
# also sources; this script only renders the standalone copy.  It is the way to
# look at the panel without panel B around it.
#
# Run 13z2_homd_phylum_band_composition.py first.
source("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3/scripts/13z2_panel_homd_band_composition.R")

PIPE_V3 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIGOUT  <- file.path(PIPE_V3, "analysis_report", "figs")
dir.create(FIGOUT, showWarnings = FALSE, recursive = TRUE)

comp <- band_comp_read(file.path(PIPE_V3, "13z2_homd_phylum_band_composition.tsv"))
p <- band_comp_panel(comp)

FIG <- file.path(FIGOUT, "fig_v3_panelB_band_composition_preview")
# the legend sits under the bars on one line, so the canvas is a little wider
# than the facets need; the spare width is what the six keys sit in.
ggsave(paste0(FIG, ".png"), p, width = 7.4, height = 5.0, dpi = 300)
ggsave(paste0(FIG, ".pdf"), p, width = 7.4, height = 5.0, device = cairo_pdf)
message("[13z2] preview written: ", FIG, ".png")
