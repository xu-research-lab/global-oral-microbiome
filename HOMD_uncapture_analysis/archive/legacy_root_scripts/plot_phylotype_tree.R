library(ape)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)

BASE <- '/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis'

# ---- 1. Read data ----
cat('Reading data...\n')
phylo <- read_tsv(file.path(BASE, 'pipeline/07_novelty_sepp/phylotypes_novelty_classified.tsv'), show_col_types = FALSE)
tax <- read_tsv('/beegfs/dongbiao/greengene2/exported-taxonomy/taxonomy.tsv', show_col_types = FALSE)
names(tax) <- c('feature_id', 'taxonomy')
placements <- read_tsv(file.path(BASE, 'pipeline/05_sepp/sepp_placements_parsed.tsv'), show_col_types = FALSE)
annot <- read_tsv(file.path(BASE, 'homd_mapping/exact_derep_homd_annotated.tsv'), show_col_types = FALSE)
annot <- annot[, c('derep_id', 'n_bioprojects', 'total_reads', 'prevalence_samples', 'continents', 'subject_prevalence')]

# ---- 2. Build phylotype-level annotation ----
cat('Building phylotype annotations...\n')
phylo_long <- phylo %>%
  separate_rows(member_derep_ids, sep = ';') %>%
  left_join(annot, by = c('member_derep_ids' = 'derep_id'))

# Build aggregate annotations from ASV-level data
phylo_summary <- phylo_long %>%
  group_by(phylotype_id, novelty_class) %>%
  summarise(
    n_members = n(),
    n_bioprojects = max(n_bioprojects, na.rm = TRUE),
    prevalence = max(subject_prevalence, na.rm = TRUE),
    med_rel_abund = median(total_reads / prevalence_samples, na.rm = TRUE),
    continents = paste(unique(unlist(strsplit(paste(na.omit(continents), collapse = ';'), ';'))), collapse = ';'),
    .groups = 'drop'
  ) %>%
  mutate(
    n_bioprojects = ifelse(is.infinite(n_bioprojects) | is.nan(n_bioprojects), 0, n_bioprojects),
    prevalence = ifelse(is.infinite(prevalence) | is.nan(prevalence), 0, prevalence),
    med_rel_abund = ifelse(is.infinite(med_rel_abund) | is.nan(med_rel_abund), 0, med_rel_abund)
  )

# Join phylotype-level QC columns
phylo_qc <- phylo %>%
  select(phylotype_id, n_nonchimeric, gg2_best_id, silva_best_id, homd_best_id, candidate_novel)

phylo_annot <- phylo_summary %>%
  left_join(phylo_qc, by = 'phylotype_id')

# ---- 3. Assign phylum taxonomy ----
cat('Assigning phyla...\n')
phylo_phylum <- phylo %>%
  separate_rows(member_derep_ids, sep = ';') %>%
  left_join(placements, by = c('member_derep_ids' = 'query_id')) %>%
  filter(!is.na(nearest_ref), nearest_ref != '') %>%
  left_join(tax, by = c('nearest_ref' = 'feature_id')) %>%
  mutate(phylum = sub('.*p__([^;]+).*', '\\1', taxonomy)) %>%
  group_by(phylotype_id) %>%
  summarise(
    phylum = names(which.max(table(phylum))),
    full_tax = names(which.max(table(taxonomy))),
    .groups = 'drop'
  )

# Phylotypes without nearest_ref -> Unclassified
all_phylo_ids <- data.frame(phylotype_id = phylo[['phylotype_id']])
phylo_phylum <- all_phylo_ids %>%
  left_join(phylo_phylum, by = 'phylotype_id') %>%
  mutate(
    phylum = ifelse(is.na(phylum), 'Unclassified', phylum),
    full_tax = ifelse(is.na(full_tax), 'd__Bacteria; p__Unclassified', full_tax)
  )

phylo_annot <- phylo_annot %>%
  left_join(select(phylo_phylum, phylotype_id, phylum, full_tax), by = 'phylotype_id')

cat('Phyla:', sort(unique(phylo_annot[['phylum']])), '\n')

# ---- 4. Build taxonomy tree ----
cat('Building taxonomy tree...\n')
# For each phylotype, get the full taxonomy path
# Parse taxonomy into ranks
parse_tax_levels <- function(tax_str) {
  parts <- unlist(strsplit(tax_str, '; '))
  ranks <- c('d__' = 'domain', 'p__' = 'phylum', 'c__' = 'class', 
             'o__' = 'order', 'f__' = 'family', 'g__' = 'genus')
  result <- c()
  for (r in names(ranks)) {
    idx <- grep(r, parts)
    if (length(idx) > 0) {
      val <- sub(paste0('.*', r, '([^;]+).*'), '\\1', parts[idx[1]])
      result[ranks[r]] <- val
    }
  }
  result
}

# Build a taxonomic hierarchy tree
tax_hierarchy <- phylo_annot %>%
  rowwise() %>%
  mutate(
    tax_ranks = list(parse_tax_levels(full_tax))
  ) %>%
  ungroup()

# Create a tree from scratch using phylo object
# Structure: root -> domain -> phylum -> tip (phylotype)
# We'll build an edge matrix manually

# Get unique phyla (sorted)
phyla <- sort(unique(tax_hierarchy[['phylum']]))
n_phyla <- length(phyla)

# Sort phylotypes within each phylum by gg2_best_id (descending) then prevalence
tax_hierarchy <- tax_hierarchy %>%
  mutate(gg2_best_id = as.numeric(gg2_best_id)) %>%
  arrange(phylum, desc(gg2_best_id), desc(prevalence))

n_tips <- nrow(tax_hierarchy)

# Node layout:
# Root = n_tips + 1
# Domain node = n_tips + 2
# Phylum nodes = n_tips + 3 to n_tips + 2 + n_phyla
# Tips = 1 to n_tips

tip_order <- tax_hierarchy[['phylotype_id']]
root <- n_tips + 1
domain_node <- n_tips + 2
phylum_nodes <- (n_tips + 3):(n_tips + 2 + n_phyla)
names(phylum_nodes) <- phyla

# Build edge matrix
edges <- matrix(NA, nrow = 1 + n_phyla + n_tips, ncol = 2)
edge_len <- numeric(nrow(edges))
edge_count <- 0

# Root -> domain
edge_count <- edge_count + 1
edges[edge_count, ] <- c(root, domain_node)
edge_len[edge_count] <- 1

# Domain -> each phylum
for (i in seq_along(phyla)) {
  edge_count <- edge_count + 1
  edges[edge_count, ] <- c(domain_node, phylum_nodes[phyla[i]])
  edge_len[edge_count] <- 1
}

# Each phylum -> its phylotypes
for (i in seq_len(n_tips)) {
  p <- tax_hierarchy[['phylum']][i]
  edge_count <- edge_count + 1
  edges[edge_count, ] <- c(phylum_nodes[p], i)
  edge_len[edge_count] <- 1
}

# Create phylo object
tax_tree <- list(
  edge = edges,
  edge.length = edge_len,
  Nnode = 2 + n_phyla,
  tip.label = tip_order,
  node.label = c('Root', 'Bacteria', phyla)
)
class(tax_tree) <- 'phylo'
attr(tax_tree, 'order') <- 'cladewise'

# Verify
cat('Tree tips:', length(tax_tree[['tip.label']]), '\n')
cat('Tree nodes:', tax_tree[['Nnode']], '\n')

# ---- 5. Build continent matrix ----
cat('Preparing continent matrix...\n')
continent_mat <- tax_hierarchy %>%
  select(phylotype_id, continents) %>%
  rowwise() %>%
  mutate(
    Africa = grepl('Africa', continents),
    Asia = grepl('Asia', continents),
    Europe = grepl('Europe', continents),
    `North America` = grepl('North America', continents),
    Oceania = grepl('Oceania', continents),
    `South America` = grepl('South America', continents)
  ) %>%
  ungroup() %>%
  select(-continents)

continent_long <- continent_mat %>%
  pivot_longer(-phylotype_id, names_to = 'continent', values_to = 'present') %>%
  filter(present) %>%
  select(-present)

# ---- 6. Plot ----
cat('Plotting...\n')

# QC annotation: contamination risk score
# High risk: single member, low prevalence, single bioproject
# Low risk: multiple studies, high prevalence, multiple members
tax_hierarchy <- tax_hierarchy %>%
  mutate(
    novelty_label = ifelse(
      novelty_class == 'novel_all_refs',
      'Novel (all refs <50% ID)',
      'Known (not in HOMD)'
    ),
    gg2_cat = case_when(
      is.na(gg2_best_id) | gg2_best_id == '' ~ 'No GG2 hit',
      as.numeric(gg2_best_id) < 70 ~ '<70% (low)',
      as.numeric(gg2_best_id) < 90 ~ '70-90%',
      TRUE ~ '>=90%'
    ),
    n_bioprojects_cat = case_when(
      n_bioprojects <= 1 ~ '1 (single study)',
      n_bioprojects <= 5 ~ '2-5',
      TRUE ~ '>5'
    ),
    singleton = n_members == 1
  )

# Colors
novelty_colors <- c('Novel (all refs <50% ID)' = '#E41A1C', 'Known (not in HOMD)' = '#377EB8')
gg2_colors <- c('No GG2 hit' = '#999999', '<70% (low)' = '#FF9999', '70-90%' = '#FFCC00', '>=90%' = '#66CC66')
bp_colors <- c('1 (single study)' = '#FF9999', '2-5' = '#FFCC00', '>5' = '#66CC66')

# Plot data (label column for tree matching)
plot_data <- tax_hierarchy %>%
  mutate(label = phylotype_id) %>%
  select(label, everything())

continent_long2 <- continent_long %>%
  mutate(label = phylotype_id)

# Create the cladogram (no branch length meaning)
p <- ggtree(tax_tree, layout = 'rectangular', linewidth = 0.3, branch.length = 'none') %<+% plot_data +
  geom_tiplab(aes(label = paste0(phylum, ' | ', label)), size = 2, hjust = -0.05) +
  geom_tippoint(aes(color = novelty_label), size = 2, alpha = 0.85) +
  scale_color_manual(values = novelty_colors, name = 'Novelty Class') +
  xlim(0, 12) +
  theme(
    legend.position = 'right',
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.margin = margin(5, 5, 5, 5)
  )

# --- Annotation rings (using tree data from %<+%) ---
# Ring 1: n_bioprojects
p <- p + new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(x = 1, y = label, fill = n_bioprojects_cat),
    width = 0.6,
    pwidth = 0.15,
    offset = 0.03,
    axis.params = list(title = 'Bioprojects', title.size = 2.5, text.size = 1.5, vjust = 0.5)
  ) +
  scale_fill_manual(values = bp_colors, name = 'Bioprojects', guide = guide_legend(ncol = 1, keyheight = 0.5, keywidth = 0.5))

# Ring 2: GG2 best identity
p <- p + new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(x = 1, y = label, fill = gg2_cat),
    width = 0.6,
    pwidth = 0.15,
    offset = 0.03,
    axis.params = list(title = 'GG2 ID', title.size = 2.5, text.size = 1.5, vjust = 0.5)
  ) +
  scale_fill_manual(values = gg2_colors, name = 'GG2 Identity', guide = guide_legend(ncol = 1, keyheight = 0.5, keywidth = 0.5))

# Ring 3: singleton (contamination risk)
p <- p + new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(x = 1, y = label, fill = singleton),
    width = 0.6,
    pwidth = 0.1,
    offset = 0.03,
    axis.params = list(title = 'Singleton', title.size = 2.5, text.size = 1.5, vjust = 0.5)
  ) +
  scale_fill_manual(values = c('TRUE' = '#FF6666', 'FALSE' = '#66CC99'), name = 'Single\nmember?', guide = guide_legend(ncol = 1, keyheight = 0.5, keywidth = 0.5))

# Ring 4: median rel abundance (bar)
p <- p + new_scale_fill() +
  geom_fruit(
    geom = geom_bar,
    mapping = aes(x = med_rel_abund, y = label, fill = novelty_label),
    stat = 'identity',
    orientation = 'y',
    width = 0.6,
    pwidth = 0.25,
    axis.params = list(title = 'Median\nrel. abund.', title.size = 2.5, text.size = 1.5, vjust = 0.5)
  ) +
  scale_fill_manual(values = novelty_colors, guide = 'none')

# Ring 5: subject prevalence (bar)
p <- p + new_scale_fill() +
  geom_fruit(
    geom = geom_bar,
    mapping = aes(x = prevalence, y = label, fill = novelty_label),
    stat = 'identity',
    orientation = 'y',
    width = 0.6,
    pwidth = 0.25,
    axis.params = list(title = 'Subject\nprevalence', title.size = 2.5, text.size = 1.5, vjust = 0.5)
  ) +
  scale_fill_manual(values = novelty_colors, guide = 'none')

# Ring 6: continent heatmap (long format data needed for one-to-many)
p <- p + new_scale_fill() +
  geom_fruit(
    data = continent_long2,
    geom = geom_tile,
    mapping = aes(x = continent, y = label, fill = continent),
    pwidth = 0.6,
    offset = 0.05,
    axis.params = list(
      title = 'Continent',
      title.size = 2.5,
      text.size = 1.8,
      angle = 30,
      hjust = 1
    )
  ) +
  scale_fill_manual(
    values = c(
      'Africa' = '#FFD700',
      'Asia' = '#FF6347',
      'Europe' = '#4169E1',
      'North America' = '#32CD32',
      'Oceania' = '#FF69B4',
      'South America' = '#9370DB'
    ),
    name = 'Continent',
    guide = guide_legend(ncol = 1, keyheight = 0.5, keywidth = 0.5)
  )

# ---- 7. Save ----
out_pdf <- file.path(BASE, 'figures/fig5_phylotype_tree.pdf')
out_png <- file.path(BASE, 'figures/fig5_phylotype_tree.png')
cat('Saving to', out_pdf, '\n')
ggsave(out_pdf, p, width = 16, height = 12, dpi = 300)
ggsave(out_png, p, width = 16, height = 12, dpi = 300)
cat('Done!\n')
