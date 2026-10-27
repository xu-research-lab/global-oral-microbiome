#!/usr/bin/env Rscript
# pipeline_v3: species-level HOMD coverage scatter.
# x = fraction of ASVs <97% identity; y = No. derep ASVs (log10);
# point size = No. HOMD reference sequences per species;
# marginal density on top/right; label best-10 / worst-10 covered species.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(ggrepel)
  library(ggExtra)
})

PIPE_V3 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIGOUT  <- file.path(PIPE_V3, "analysis_report", "figs")
dir.create(FIGOUT, showWarnings = FALSE, recursive = TRUE)

s <- read.delim(file.path(PIPE_V3, "13j_homd_species_coverage_fullname.tsv"),
                sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
s$frac_ge97 <- pmin(pmax(s$frac_ge97, 0), 1)
s$frac_lt97 <- 1 - s$frac_ge97

# label best-10 / worst-10 by frac_ge97 among species with >=1000 ASVs
big <- s %>% filter(n_asvs_total >= 1000)
worst <- big %>% arrange(frac_ge97) %>% slice_head(n = 10)
best  <- big %>% arrange(desc(frac_ge97)) %>% slice_head(n = 10)
lab <- bind_rows(worst, best) %>% distinct(species, .keep_all = TRUE)

p <- ggplot(s, aes(x = frac_lt97, y = n_asvs_total)) +
  geom_point(aes(size = n_homd_refseq), colour = "grey35", alpha = 0.65) +
  geom_text_repel(data = lab, aes(label = species),
                  size = 2.3, max.overlaps = 30, segment.size = 0.2,
                  min.segment.length = 0, box.padding = 0.35,
                  colour = "black", fontface = "italic") +
  scale_x_continuous(labels = percent_format(),
                     name = "Fraction of ASVs <97% identity (HOMD coverage gap)") +
  scale_y_log10(labels = label_number(scale_cut = cut_short_scale()),
                name = "No. derep ASVs mapped (log10)") +
  scale_size_continuous(range = c(0.8, 7),
                        breaks = c(1, 5, 10, 25, 50, 100, 178),
                        name = "HOMD reference\nsequences per species") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right", panel.grid.minor = element_blank())

p <- ggMarginal(p, type = "density", margins = "both", size = 5,
                colour = "grey40", fill = "grey75", alpha = 0.5)

fig <- file.path(FIGOUT, "fig_v3_homd_species_coverage")
ggsave(paste0(fig, ".png"), p, width = 8.5, height = 6.5, dpi = 300)
ggsave(paste0(fig, ".pdf"), p, width = 8.5, height = 6.5)
message("[13k] species scatter written")
