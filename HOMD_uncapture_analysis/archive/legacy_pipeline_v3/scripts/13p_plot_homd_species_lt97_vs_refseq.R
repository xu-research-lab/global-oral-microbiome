#!/usr/bin/env Rscript
# pipeline_v3: scatter of species-level coverage metric
#   x = Fraction of ASVs with best HOMD identity < 97%
#   y = No. HOMD reference sequences per species (log10, proxy for reference size)
# with linear trend line (lm on log10 y) and annotated Pearson R / Spearman rho + P.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

PIPE_V3 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIGOUT  <- file.path(PIPE_V3, "analysis_report", "figs")
dir.create(FIGOUT, showWarnings = FALSE, recursive = TRUE)

s <- read.delim(file.path(PIPE_V3, "13o2_homd_species_lt97.tsv"),
                sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
s$frac_lt97 <- pmin(pmax(s$frac_lt97, 0), 1)

# ---- statistics ----
pe <- cor.test(s$frac_lt97, log10(s$n_homd_refseq), method = "pearson")
sp <- cor.test(s$frac_lt97, s$n_homd_refseq, method = "spearman")
pe_r  <- unname(pe$estimate)
pe_p  <- pe$p.value
sp_r  <- unname(sp$estimate)
sp_p  <- sp$p.value
n_sp  <- nrow(s)

lab_txt <- sprintf("Pearson R = %.3f (P = %.2g)\nSpearman rho = %.3f (P = %.2g)\nn = %d species",
                   pe_r, pe_p, sp_r, sp_p, n_sp)
cat(lab_txt, "\n")

p <- ggplot(s, aes(x = frac_lt97, y = n_homd_refseq)) +
  geom_point(colour = "grey30", alpha = 0.55, size = 1.8) +
  geom_smooth(method = "lm", se = TRUE, colour = "#d73027", fill = "#d73027",
              alpha = 0.18, linewidth = 0.8) +
  annotate("text", x = 0.015, y = 250, label = lab_txt, hjust = 0, vjust = 1,
           size = 3.4, colour = "grey15", family = "mono", lineheight = 0.95) +
  scale_x_continuous(labels = percent_format(), limits = c(0, 1),
                     name = "Fraction of ASVs with identity < 97% to HOMD (per species)") +
  scale_y_log10(labels = label_number(), limits = c(0.5, 400),
                breaks = c(1, 2, 5, 10, 20, 50, 100, 200),
                name = "No. HOMD reference sequences per species (log10)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

fig <- file.path(FIGOUT, "fig_v3_homd_species_lt97_vs_refseq")
ggsave(paste0(fig, ".png"), p, width = 7.2, height = 5.6, dpi = 300)
ggsave(paste0(fig, ".pdf"), p, width = 7.2, height = 5.6)
message("[13p] scatter lt97 vs refseq written")
