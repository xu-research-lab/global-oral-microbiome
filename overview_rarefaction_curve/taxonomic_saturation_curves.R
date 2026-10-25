suppressPackageStartupMessages({
  library(Matrix)
  library(hdf5r)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mgcv)
  library(scales)
})

base_sz <- 22

update_geom_defaults("text", list(family = "sans"))
update_geom_defaults("label", list(family = "sans"))

theme_set(theme_minimal(base_size = base_sz, base_family = "sans"))

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
out_dir <- if (length(script_file)) dirname(normalizePath(script_file)) else getwd()

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (!length(hit)) stop("None of these paths exists:\n", paste(paths, collapse = "\n"))
  hit[[1]]
}

biom_file <- first_existing(c(
  file.path(out_dir, "..", "data", "metadata-and-biom-v2-qc", "biom_oral_qc.biom"),
  file.path(out_dir, "..", "input_1021", "biom_oral_qc.biom")
))
taxonomy_file <- first_existing(c(
  file.path(out_dir, "..", "data", "metadata-and-biom-v2-qc", "taxonomy.tsv"),
  file.path(out_dir, "..", "input_1021", "taxonomy.tsv")
))

tsv_out <- file.path(out_dir, "rarefaction_pcofg.tsv")
png_out <- file.path(out_dir, "rarefaction_curve.png")


REPLICATES <- 149
todo <- c(1, 10, 100, 500, 1000,
          seq(from = 2000, to = 10000, by = 2000),
          seq(from = 15000, to = 35000, by = 5000),
          seq(from = 50000, to = 110000, by = 15000),
          125000, 150000)

rarefaction_sampling <- function(taxtable) {
  richness <- as(t(taxtable) > 0, "dgCMatrix")

  observed <- sum(Matrix::rowSums(richness) > 0)
  rarefaction <- data.frame(ncol(richness), observed)
  names(rarefaction) <- c("scount", "observed")

  for (scount in todo) {
    message("SAMPLE SIZE: ", scount)
    for (iter in seq(0, REPLICATES)) {
      sub <- richness[, sample.int(ncol(richness), scount), drop = FALSE]
      observed <- sum(Matrix::rowSums(sub) > 0)
      result <- data.frame(scount, observed)
      names(result) <- c("scount", "observed")
      rarefaction <- rbind(rarefaction, result)
    }
  }

  message("SAMPLE SIZE: ALL")
  observed <- sum(Matrix::rowSums(richness) > 0)
  result <- data.frame(ncol(richness), observed)
  names(result) <- c("scount", "observed")
  rarefaction <- rbind(rarefaction, result)

  return(rarefaction)
}

combine_taxa <- function(dataset, labels) {
  keep <- !is.na(labels)
  labels <- labels[keep]
  dataset <- dataset[, keep, drop = FALSE]
  taxa <- sort(unique(labels))
  map <- sparseMatrix(
    i = match(labels, taxa),
    j = seq_along(labels),
    x = 1,
    dims = c(length(taxa), length(labels))
  )
  final <- as(dataset %*% t(map), "dgCMatrix")
  colnames(final) <- taxa
  rownames(final) <- rownames(dataset)
  return(final)
}

message("Reading BIOM table: ", biom_file)
h5 <- H5File$new(biom_file, mode = "r")
feature_ids <- h5[["observation/ids"]][]
sample_ids <- h5[["sample/ids"]][]
raredata <- sparseMatrix(
  i = h5[["observation/matrix/indices"]][] + 1L,
  p = h5[["observation/matrix/indptr"]][],
  x = h5[["observation/matrix/data"]][],
  dims = c(length(sample_ids), length(feature_ids))
)
h5$close_all()
rownames(raredata) <- sample_ids
colnames(raredata) <- feature_ids
message(sprintf("  %d samples x %d features, %d non-zero values",
                nrow(raredata), ncol(raredata), length(raredata@x)))

message("Reading taxonomy: ", taxonomy_file)
taxonomy <- fread(taxonomy_file, sep = "\t", header = TRUE, quote = "",
                  data.table = FALSE, showProgress = FALSE)

extract_rank <- function(taxon, prefix) {
  pattern <- paste0("(?:^|;)\\s*", prefix, "([^;]*)")
  hits <- regmatches(taxon, regexec(pattern, taxon, perl = TRUE))
  value <- vapply(hits, function(x) {
    if (length(x) == 2L) trimws(x[[2]]) else NA_character_
  }, character(1))
  value[!is.na(value) & !nzchar(value)] <- NA_character_
  value
}

tax_row <- match(colnames(raredata), taxonomy[[1]])
if (anyNA(tax_row)) {
  warning(sum(is.na(tax_row)), " features have no taxonomy entry")
}
taxon_strings <- taxonomy[[2]][tax_row]
rank_labels <- lapply(c("p__", "c__", "o__", "f__", "g__"),
                      function(p) extract_rank(taxon_strings, p))

taxphylum <- combine_taxa(raredata, rank_labels[[1]])
taxclass <- combine_taxa(raredata, rank_labels[[2]])
taxorder <- combine_taxa(raredata, rank_labels[[3]])
taxfamily <- combine_taxa(raredata, rank_labels[[4]])
taxgenus <- combine_taxa(raredata, rank_labels[[5]])

for (nm in c("phylum", "class", "order", "family", "genus")) {
  message(sprintf("  %s: %d unique taxa",
                  nm, ncol(get(paste0("tax", nm)))))
}

dropped <- todo[todo > nrow(raredata)]
todo <- todo[todo <= nrow(raredata)]
if (length(dropped)) {
  message("Skipping sizes larger than the dataset: ", paste(dropped, collapse = ", "))
}
message("Compendium sizes: ", paste(todo, collapse = ", "),
        " + full table (", nrow(raredata), ")")

set.seed(42)
message("phylum")
rarefaction.phylum <- rarefaction_sampling(taxphylum)
rarefaction.phylum$level <- "phylum"
message("class")
rarefaction.class <- rarefaction_sampling(taxclass)
rarefaction.class$level <- "class"
message("order")
rarefaction.order <- rarefaction_sampling(taxorder)
rarefaction.order$level <- "order"
message("family")
rarefaction.family <- rarefaction_sampling(taxfamily)
rarefaction.family$level <- "family"
message("genus")
rarefaction.genus <- rarefaction_sampling(taxgenus)
rarefaction.genus$level <- "genus"

message("compiling")
rarefaction <- rbind(rarefaction.phylum, rarefaction.class,
                     rarefaction.order, rarefaction.family,
                     rarefaction.genus)
rarefaction$level <- factor(rarefaction$level,
                            levels = c("phylum", "class", "order", "family", "genus"))

rm(rarefaction.phylum, rarefaction.class, rarefaction.order,
   rarefaction.family, rarefaction.genus)

rarefaction_summary <- rarefaction %>%
  group_by(level, scount) %>%
  summarise(
    sd = sd(observed, na.rm = TRUE),
    observed = mean(observed),
    n_replicates = n(),
    .groups = "drop"
  )

rank_names <- c(Phylum = "phylum", Class = "class", Order = "order",
                Family = "family", Genus = "genus")
pcofg <- data.frame(Samples = sort(unique(rarefaction_summary$scount)))
for (nm in names(rank_names)) {
  d <- rarefaction_summary[rarefaction_summary$level == rank_names[[nm]], ]
  d <- d[match(pcofg$Samples, d$scount), ]
  pcofg[[nm]] <- d$observed
  pcofg[[paste0(nm, "_sd")]] <- d$sd
}
write.table(pcofg, file = tsv_out, sep = "\t", quote = FALSE, row.names = FALSE)
message("Saved ", tsv_out)

df_rare <- rarefaction_summary %>%
  mutate(Taxonomy = tools::toTitleCase(as.character(level))) %>%
  select(Samples = scount, Taxonomy, Unique_taxa = observed) %>%
  mutate(Taxonomy = factor(Taxonomy,
                           levels = c("Genus", "Family", "Order",
                                      "Class", "Phylum")))

rare_plot <- ggplot(df_rare, aes(x = Samples, y = Unique_taxa, color = Taxonomy)) +
  geom_smooth(
    method    = "gam",
    formula   = y ~ s(x, bs = "cs"),
    se        = FALSE,
    linewidth = 0.7
  ) +
  scale_color_manual(
    values = c(
      "Genus"  = "#5C6BC0",
      "Family" = "#FFB300",
      "Order"  = "#EF5350",
      "Class"  = "#26C6DA",
      "Phylum" = "#66BB6A"
    )
  ) +
  scale_x_continuous(
    breaks = seq(0, 80000, by = 20000),
    limits = c(0, 98000),
    labels = c("0", "20K", "40K", "60K", "80K"),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = c(0, 1000, 2000, 3000, 4000, 5000),
    labels = c("0", "1K", "2K", "3K", "4K", "5K"),
    expand = expansion(mult = c(0.1, 0.06))
  ) +

  annotate("text",
           x = 81000, y = 4900, label = "Genus",
           color = "#5C6BC0", fontface = "bold",
           size = base_sz / 3.8, hjust = 0) +
  annotate("text",
           x = 81000, y = 1400, label = "Family",
           color = "#FFB300", fontface = "bold",
           size = base_sz / 3.8, hjust = 0) +
  annotate("text",
           x = 81000, y = 800, label = "Order",
           color = "#EF5350", fontface = "bold",
           size = base_sz / 3.8, hjust = 0) +
  annotate("text",
           x = 81000, y = 350, label = "Class",
           color = "#26C6DA", fontface = "bold",
           size = base_sz / 3.8, hjust = 0) +
  annotate("text",
           x = 81000, y = 5, label = "Phylum",
           color = "#66BB6A", fontface = "bold",
           size = base_sz / 3.8, hjust = 0) +
  labs(x = "Samples", y = "Observed taxa") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = base_sz - 4, base_family = "sans") +
  theme(
    axis.text.x     = element_text(
      size = base_sz * 0.72,
      face = "bold", color = "#303030"
    ),
    axis.text.y     = element_text(
      size = base_sz * 0.72,
      face = "bold", color = "#303030"
    ),
    axis.title.y    = element_text(
      size = base_sz * 0.85,
      face = "bold", color = "#303030",
      margin = margin(r = 6)
    ),
    axis.title.x    = element_text(
      size = base_sz * 0.85,
      face = "bold", color = "#303030",
      margin = margin(t = 10)
    ),
    legend.position = "none",
    plot.margin     = margin(5, 5, 5, 5, unit = "pt"),
    plot.background = element_blank(),
    panel.background = element_blank()
  )

ggsave(png_out, plot = rare_plot, width = 7, height = 5.5,
       device = "png", dpi = 600, bg = "white")
message("Saved ", png_out)

final <- rarefaction_summary %>%
  filter(scount == max(scount)) %>%
  mutate(level = tools::toTitleCase(as.character(level))) %>%
  arrange(desc(observed)) %>%
  mutate(percent_of_total = 100 * observed / c(
    Phylum = ncol(taxphylum), Class = ncol(taxclass), Order = ncol(taxorder),
    Family = ncol(taxfamily), Genus = ncol(taxgenus)
  )[level])
print(as.data.frame(final[, c("level", "observed", "sd", "percent_of_total")]),
      digits = 4)
