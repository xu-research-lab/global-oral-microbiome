#!/usr/bin/env Rscript
options(ignore.negative.edge = TRUE)
suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
  library(ggtree)
  library(ggtreeExtra)
  library(ggnewscale)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  idx <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) args[[idx + 1]] else default
}

outdir <- get_arg('--outdir', '/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/region_aware_pipeline')
alignment_file <- get_arg('--alignment', '/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03_mafft_aligned.fasta')
tree_file <- get_arg('--tree', file.path(outdir, 'HOMD_16S_rRNA_RefSeq_V16.03_mafft_nj.nwk'))
annot_file <- get_arg('--annotation', file.path(outdir, 'homd_tree_tip_annotation.tsv'))
figdir <- file.path(outdir, 'figures')
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

message('[ggtree] loading annotation')
annot <- read_tsv(annot_file, show_col_types = FALSE)
annot <- annot %>%
  mutate(
    homd_target = as.character(homd_target),
    phylum_plot = ifelse(is.na(phylum_plot) | phylum_plot == '', 'NA', phylum_plot),
    log10_any_hit_derep = log10(any_hit_derep_sequences + 1),
    median_identity_for_plot = ifelse(any_hit_derep_sequences > 0, any_hit_median_identity, NA_real_),
    n_bioprojects_any = ifelse(is.na(n_bioprojects_any), 0, n_bioprojects_any),
    bioproject_bin = cut(n_bioprojects_any, breaks = c(-Inf, 0, 1, 3, 10, 30, Inf), labels = c('0', '1', '2-3', '4-10', '11-30', '>30')),
    burden_bin = cut(any_hit_derep_sequences, breaks = c(-Inf, 0, 10, 100, 1000, 5000, Inf), labels = c('0', '1-10', '11-100', '101-1k', '1k-5k', '>5k')),
    identity_bin = cut(median_identity_for_plot, breaks = c(-Inf, 80, 90, 95, 97, 99, 99.5, 100.0001), labels = c('<80/no hit', '80-90', '90-95', '95-97', '97-99', '99-99.5', '99.5-100'))
  )

if (!file.exists(tree_file)) {
  message('[ggtree] tree file not found; building NJ tree with ape::dist.dna')
  aln <- read.dna(alignment_file, format = 'fasta')
  d <- dist.dna(aln, model = 'JC69', pairwise.deletion = TRUE, as.matrix = FALSE)
  tr <- nj(d)
  tr <- ladderize(tr)
  write.tree(tr, tree_file)
  message('[ggtree] wrote tree: ', tree_file)
} else {
  message('[ggtree] using existing tree: ', tree_file)
  tr <- read.tree(tree_file)
  tr <- ladderize(tr)
}

# Match annotation to tree tips and keep plotting robust if a few tips are missing.
plot_annot <- annot %>% filter(homd_target %in% tr$tip.label) %>% mutate(ID = homd_target)
missing_annot <- setdiff(tr$tip.label, plot_annot$homd_target)
if (length(missing_annot) > 0) {
  plot_annot <- bind_rows(plot_annot, tibble(
    homd_target = missing_annot,
    ID = missing_annot,
    phylum_plot = 'NA',
    any_hit_derep_sequences = 0,
    log10_any_hit_derep = 0,
    median_identity_for_plot = NA_real_,
    n_bioprojects_any = 0,
    burden_bin = factor('0', levels = levels(annot$burden_bin)),
    identity_bin = factor('<80/no hit', levels = levels(annot$identity_bin)),
    bioproject_bin = factor('0', levels = levels(annot$bioproject_bin))
  ))
}

phyla <- sort(unique(plot_annot$phylum_plot))
base_cols <- c(
  'Bacillota' = '#4E79A7', 'Pseudomonadota' = '#F28E2B', 'Bacteroidota' = '#59A14F',
  'Actinomycetota' = '#E15759', 'Fusobacteriota' = '#B07AA1', 'Campylobacterota' = '#9C755F',
  'Spirochaetota' = '#76B7B2', 'Saccharibacteria' = '#EDC948', 'Synergistota' = '#FF9DA7',
  'Cyanobacteriota' = '#8CD17D', 'Methanobacteriota' = '#BAB0AC', 'Other' = '#999999', 'NA' = '#D0D0D0'
)
phyla_cols <- setNames(rep('#999999', length(phyla)), phyla)
phyla_cols[names(base_cols)[names(base_cols) %in% phyla]] <- base_cols[names(base_cols)[names(base_cols) %in% phyla]]


ring_burden <- plot_annot %>% select(ID, ring_value = log10_any_hit_derep)
ring_identity <- plot_annot %>% select(ID, ring_value = median_identity_for_plot)
ring_bioprojects <- plot_annot %>% select(ID, ring_value = n_bioprojects_any)
ring_fraction_id99 <- plot_annot %>% select(ID, ring_value = fraction_asv_id99_qcov99)
ring_gap_score <- plot_annot %>% select(ID, ring_value = reference_gap_score)

message('[ggtree] plotting circular tree')
p <- ggtree(tr, layout = 'circular', size = 0.08, open.angle = 8) %<+% plot_annot +
  geom_tippoint(aes(color = phylum_plot), size = 0.28, alpha = 0.9) +
  scale_color_manual(values = phyla_cols, name = 'Phylum', na.value = '#D0D0D0') +
  theme(legend.position = 'right')

p <- p + ggnewscale::new_scale_fill()
p <- p + geom_fruit(
  data = ring_burden,
  geom = geom_tile,
  mapping = aes(y = ID, x = 'ASV burden', fill = ring_value),
  offset = 0.020,
  width = 0.055,
  pwidth = 0.055,
  axis.params = list(axis = 'x', text.size = 2.0, text.angle = -45, hjust = 0),
  grid.params = list()
) + scale_fill_viridis_c(name = 'log10(unique ASV + 1)', option = 'magma', na.value = '#F0F0F0')

p <- p + ggnewscale::new_scale_fill()
p <- p + geom_fruit(
  data = ring_identity,
  geom = geom_tile,
  mapping = aes(y = ID, x = 'Median identity', fill = ring_value),
  offset = 0.080,
  width = 0.055,
  pwidth = 0.055,
  axis.params = list(axis = 'x', text.size = 2.0, text.angle = -45, hjust = 0),
  grid.params = list()
) + scale_fill_gradientn(
  colors = c('#762A83', '#9970AB', '#F7F7F7', '#5AAE61', '#1B7837'),
  values = rescale(c(80, 90, 95, 99, 100)),
  limits = c(80, 100),
  oob = squish,
  name = 'Median best identity (%)',
  na.value = '#E6E6E6'
)

p <- p + ggnewscale::new_scale_fill()
p <- p + geom_fruit(
  data = ring_bioprojects,
  geom = geom_tile,
  mapping = aes(y = ID, x = 'Bioprojects', fill = ring_value),
  offset = 0.140,
  width = 0.055,
  pwidth = 0.055,
  axis.params = list(axis = 'x', text.size = 2.0, text.angle = -45, hjust = 0),
  grid.params = list()
) + scale_fill_viridis_c(name = 'Bioprojects with best-hit ASVs', option = 'cividis', trans = 'sqrt', na.value = '#F0F0F0')

p <- p + ggnewscale::new_scale_fill()
p <- p + geom_fruit(
  data = ring_fraction_id99,
  geom = geom_tile,
  mapping = aes(y = ID, x = 'Fraction id99', fill = ring_value),
  offset = 0.200,
  width = 0.055,
  pwidth = 0.055,
  axis.params = list(axis = 'x', text.size = 2.0, text.angle = -45, hjust = 0),
  grid.params = list()
) + scale_fill_gradientn(
  colors = c('#8B0000', '#D95F02', '#FEE08B', '#91CF60', '#1A9850'),
  values = rescale(c(0, 0.25, 0.5, 0.75, 1)),
  limits = c(0, 1),
  oob = squish,
  name = 'Fraction ASVs id >=99%',
  na.value = '#E6E6E6'
)

p <- p + ggnewscale::new_scale_fill()
p <- p + geom_fruit(
  data = ring_gap_score,
  geom = geom_tile,
  mapping = aes(y = ID, x = 'Gap score', fill = ring_value),
  offset = 0.260,
  width = 0.055,
  pwidth = 0.055,
  axis.params = list(axis = 'x', text.size = 2.0, text.angle = -45, hjust = 0),
  grid.params = list()
) + scale_fill_viridis_c(name = 'Reference gap score', option = 'inferno', trans = 'sqrt', na.value = '#F0F0F0')

# Write top-priority references separately. Tip labels are omitted from the full tree
# to avoid overplotting and to keep the ring visualization robust.
top_labels <- plot_annot %>%
  filter(any_hit_derep_sequences >= 100, n_bioprojects_any >= 10, !is.na(median_identity_for_plot), median_identity_for_plot < 97) %>%
  arrange(desc(reference_gap_score)) %>%
  slice_head(n = 20) %>%
  mutate(label_text = paste0(homd_hmt, ' ', genus, ' ', species_epithet))

out_pdf <- file.path(figdir, 'fig6_homd_phylogeny_asv_mapping_5rings.pdf')
out_png <- file.path(figdir, 'fig6_homd_phylogeny_asv_mapping_5rings.png')
message('[ggtree] writing ', out_pdf)
ggsave(out_pdf, p, width = 13.5, height = 12.5, limitsize = FALSE, device = cairo_pdf)
ggsave(out_png, p, width = 13.5, height = 12.5, limitsize = FALSE, dpi = 350)

# Also write a compact top-priority table for interpretation.
top_out <- file.path(outdir, 'homd_tree_reference_gap_score_top20.tsv')
top_labels %>%
  select(homd_target, homd_hmt, homd_species, phylum, genus, species_epithet, any_hit_derep_sequences,
         any_hit_source_asvs, id99p0_qcov99_source_asvs, fraction_asv_id99_qcov99, any_hit_total_reads,
         any_hit_median_identity, n_bioprojects_any, n_countries_any, n_continents_any,
         clipped_identity_gap_to_99, reference_gap_score) %>%
  write_tsv(top_out)
message('[ggtree] wrote ', top_out)
