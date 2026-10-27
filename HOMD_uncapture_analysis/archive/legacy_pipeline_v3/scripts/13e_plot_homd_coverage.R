#!/usr/bin/env Rscript
# pipeline_v3: HOMD coverage figures from full clean (non-chimeric, non-contaminant,
# non-organelle) derep ASVs mapped to HOMD.
#  - Fig A (phylum): stacked bar chart of ASV identity bins per phylum
#  - Fig B (species): scatter of ASV count vs fraction of ASVs <97% identity
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(ggrepel)
})

PIPE_V3 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIGOUT  <- file.path(PIPE_V3, "analysis_report", "figs")
dir.create(FIGOUT, showWarnings = FALSE, recursive = TRUE)

# ---------------- phylum ----------------
p <- read.delim(file.path(PIPE_V3, "13d_homd_phylum_coverage.tsv"),
                sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
p$phylum <- factor(p$phylum, levels = p$phylum)  # already sorted by n_asvs_total desc
p <- p %>%
  mutate(n_ge95_97 = n_asvs_ge95 - n_asvs_ge97,
         n_ge90_95 = n_asvs_ge90 - n_asvs_ge95,
         n_lt90    = n_asvs_lt90) %>%
  select(phylum, n_asvs_total, n_homd_refseq,
         n_ge97 = n_asvs_ge97, n_ge95_97, n_ge90_95, n_lt90)

p_long <- p %>%
  pivot_longer(cols = c(n_ge97, n_ge95_97, n_ge90_95, n_lt90),
               names_to = "bin", values_to = "n") %>%
  mutate(bin = factor(bin, levels = c("n_ge97", "n_ge95_97", "n_ge90_95", "n_lt90"),
                      labels = c(">=97%", "95-97%", "90-95%", "<90%")))

bin_cols <- c(">=97%" = "#1a9850", "95-97%" = "#a6d96a",
              "90-95%" = "#fee08b", "<90%" = "#d73027")

# Panel A: stacked counts (log10 y)
pa <- ggplot(p_long, aes(x = phylum, y = n, fill = bin)) +
  geom_col(width = 0.72, colour = "grey30", linewidth = 0.15) +
  scale_y_log10(labels = label_number(scale_cut = cut_short_scale()),
                breaks = c(1e2, 1e3, 1e4, 1e5, 1e6)) +
  scale_fill_manual(values = bin_cols, name = "HOMD identity") +
  labs(x = NULL, y = "No. derep ASVs (log10)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 9),
        legend.position = "top", legend.key.size = unit(0.35, "cm"),
        panel.grid.minor = element_blank())

# Panel B: 100% stacked proportions + HOMD refseq count annotation
pb <- ggplot(p_long, aes(x = phylum, y = n, fill = bin)) +
  geom_col(position = "fill", width = 0.72, colour = "grey30", linewidth = 0.15) +
  geom_text(data = p, aes(x = phylum, y = 1.02, label = n_homd_refseq),
            inherit.aes = FALSE, size = 2.6, angle = 40, hjust = 0) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.10))) +
  scale_fill_manual(values = bin_cols, name = "HOMD identity") +
  labs(x = NULL, y = "Fraction of ASVs") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 9),
        legend.position = "top", legend.key.size = unit(0.35, "cm"),
        panel.grid.minor = element_blank())

figA <- file.path(FIGOUT, "fig_v3_homd_phylum_coverage")
ggsave(paste0(figA, ".png"), pa, width = 9, height = 4.2, dpi = 300)
ggsave(paste0(figA, ".pdf"), pa, width = 9, height = 4.2)
ggsave(paste0(figA, "_frac.png"), pb, width = 9, height = 4.2, dpi = 300)
ggsave(paste0(figA, "_frac.pdf"), pb, width = 9, height = 4.2)
message("[13e] phylum figures written")

# ---------------- species ----------------
s <- read.delim(file.path(PIPE_V3, "13d_homd_species_coverage.tsv"),
                sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
s$frac_ge97 <- pmin(pmax(s$frac_ge97, 0), 1)
s$frac_lt97 <- 1 - s$frac_ge97

# label selection: top by ASV count; plus high-abundance blind spots (frac_lt97 high)
top_abund <- s %>% arrange(desc(n_asvs_total)) %>% slice_head(n = 12)
blind <- s %>% filter(n_asvs_total >= 3000, frac_lt97 >= 0.85) %>% arrange(desc(n_asvs_total))
well <- s %>% filter(n_asvs_total >= 3000, frac_lt97 <= 0.15)
lab <- bind_rows(top_abund, blind, well) %>% distinct(species, .keep_all = TRUE)

ps <- ggplot(s, aes(x = frac_lt97, y = n_asvs_total)) +
  geom_point(aes(colour = identity_median, size = n_asvs_total), alpha = 0.75) +
  geom_text_repel(data = lab, aes(label = species, colour = identity_median),
                  size = 2.4, max.overlaps = 30, segment.size = 0.2,
                  min.segment.length = 0, box.padding = 0.35) +
  scale_x_continuous(labels = percent_format(),
                     name = "Fraction of ASVs <97% identity (HOMD coverage gap)") +
  scale_y_log10(labels = label_number(scale_cut = cut_short_scale()),
                name = "No. derep ASVs mapped (log10)") +
  scale_colour_gradient2(low = "#d73027", mid = "#fee08b", high = "#1a9850",
                         midpoint = 90, name = "Median identity\n(%)") +
  scale_size_continuous(range = c(0.8, 6), guide = "none") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right", panel.grid.minor = element_blank())

figB <- file.path(FIGOUT, "fig_v3_homd_species_coverage")
ggsave(paste0(figB, ".png"), ps, width = 8.5, height = 6.5, dpi = 300)
ggsave(paste0(figB, ".pdf"), ps, width = 8.5, height = 6.5)
message("[13e] species figures written")
