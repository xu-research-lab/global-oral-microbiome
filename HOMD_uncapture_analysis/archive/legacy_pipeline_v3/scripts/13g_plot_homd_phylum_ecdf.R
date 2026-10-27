#!/usr/bin/env Rscript
# pipeline_v3: per-phylum cumulative identity curve.
# x = best HOMD identity; y = cumulative fraction of ASVs with identity >= x.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

PIPE_V3 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIGOUT  <- file.path(PIPE_V3, "analysis_report", "figs")
dir.create(FIGOUT, showWarnings = FALSE, recursive = TRUE)

d <- read.delim(file.path(PIPE_V3, "13f_homd_phylum_ecdf.tsv"),
                sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# total ASV count per phylum for legend labels
tot <- d %>% filter(identity == 0) %>% select(phylum, n_asvs = frac_ge)
tot$n_asvs <- tot$n_asvs * 0  # placeholder; real counts come from coverage table
cov <- read.delim(file.path(PIPE_V3, "13d_homd_phylum_coverage.tsv"),
                  sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
tot <- cov %>% select(phylum, n_asvs = n_asvs_total) %>%
  mutate(label = sprintf("%s (%.0fk)", phylum, n_asvs / 1000))
d <- d %>% left_join(tot %>% select(phylum, label), by = "phylum")
d$label <- factor(d$label, levels = tot$label)  # sorted by ASV count desc

# color palette (15 phyla, colorblind-friendly)
pal <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
         "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
         "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5")

p <- ggplot(d, aes(x = identity, y = frac_ge, colour = label)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = c(75, 80, 90, 95, 97), linetype = "dotted", colour = "grey75", linewidth = 0.35) +
  geom_vline(xintercept = 97, linetype = "dashed", colour = "grey30", linewidth = 0.55) +
  annotate("text", x = c(75.4, 80.4, 90.4, 95.4, 97.4), y = rep(1.04, 5),
           label = c("phylum\n~75%", "class\n~80%", "family\n~90%", "genus\n~95%", "97%"),
           size = 2.4, colour = "grey30", hjust = 0, vjust = 1, lineheight = 0.85) +
  scale_x_continuous(breaks = seq(70, 100, 2.5), limits = c(70, 100),
                     name = "Best identity to HOMD reference (%)") +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1.05),
                     name = "Cumulative fraction of ASVs (identity >= x)") +
  scale_colour_manual(values = pal, name = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right", legend.key.size = unit(0.45, "cm"),
        legend.text = element_text(size = 8), panel.grid.minor = element_blank())

fig <- file.path(FIGOUT, "fig_v3_homd_phylum_identity_curve")
ggsave(paste0(fig, ".png"), p, width = 9.5, height = 6, dpi = 300)
ggsave(paste0(fig, ".pdf"), p, width = 9.5, height = 6)
message("[13g] phylum identity curve written")
