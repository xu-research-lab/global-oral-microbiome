#!/usr/bin/env Rscript
# pipeline_v3: per-phylum identity-bin ASV counts.
# y = 1% identity bins (70-100); x = No. ASVs in bin (log10(n+1));
# one profile line per phylum. Peaks reveal where ASVs of each phylum sit
# relative to HOMD reference identity.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

PIPE_V3 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIGOUT  <- file.path(PIPE_V3, "analysis_report", "figs")
dir.create(FIGOUT, showWarnings = FALSE, recursive = TRUE)

d <- read.delim(file.path(PIPE_V3, "13n_homd_phylum_identity_bins.tsv"),
                sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
cov <- read.delim(file.path(PIPE_V3, "13d_homd_phylum_coverage.tsv"),
                  sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

meta <- cov %>% transmute(phylum = phylum, frac_ge97 = n_asvs_ge97 / n_asvs_total) %>%
  arrange(desc(frac_ge97))
meta$label <- sprintf("%s (%.1f%%)", meta$phylum, 100 * meta$frac_ge97)
meta$label <- factor(meta$label, levels = meta$label)
d <- d %>% left_join(meta, by = "phylum")
d$label <- factor(d$label, levels = levels(meta$label))
d <- d %>% mutate(x = log10(n_asvs + 1))

# axis break positions on log10(n+1) scale with true-count labels
cnt_breaks <- c(0, 1, 10, 100, 1000, 10000, 30000, 60000)
x_br <- log10(cnt_breaks + 1)
x_lb <- ifelse(cnt_breaks == 60000, "60k",
               ifelse(cnt_breaks == 30000, "30k",
                      ifelse(cnt_breaks == 10000, "10k",
                             ifelse(cnt_breaks == 1000, "1k", as.character(cnt_breaks)))))

p <- ggplot(d, aes(x = x, y = bin_mid, colour = label)) +
  geom_hline(yintercept = 97, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = log10(1 + 0), colour = "grey85", linewidth = 0.3) +  # zero reference
  geom_line(linewidth = 0.8) +
  geom_point(size = 0.9) +
  scale_x_continuous(breaks = x_br, labels = x_lb,
                     name = "No. ASVs in identity bin (log10)") +
  scale_y_continuous(breaks = seq(70, 100, 5), limits = c(70, 100),
                     name = "Best identity to HOMD reference (%)") +
  scale_colour_manual(values = hue_pal()(15), name = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right", legend.key.size = unit(0.45, "cm"),
        legend.text = element_text(size = 7.5), panel.grid.minor = element_blank())

fig <- file.path(FIGOUT, "fig_v3_homd_phylum_identity_bins")
ggsave(paste0(fig, ".png"), p, width = 8.5, height = 5.8, dpi = 300)
ggsave(paste0(fig, ".pdf"), p, width = 8.5, height = 5.8)
message("[13o] identity-bin profile written")
