suppressPackageStartupMessages({
  library(data.table)
  library(rhdf5)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  key <- paste0("--", name, "=")
  hit <- grep(paste0("^", key), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", key), "", hit[[1]])
}

paths <- list(
  biom = get_arg("biom", "/beegfs/dongbiao/oral/table_gg2_id99_cleaned_62782_115345_all.biom"),
  oral_metadata = get_arg("oral-metadata", "/beegfs/huangxiaochang/global_oral/metadata_v2_merge/metadata-v2-merge-oral.csv"),
  outdir = get_arg("outdir", "QC"),
  prefix = get_arg("prefix", "oral_otu_sample_frequency")
)

dir.create(paths$outdir, recursive = TRUE, showWarnings = FALSE)

normalize_sample_id <- function(x) {
  x <- trimws(as.character(x))
  x <- basename(x)
  x <- sub("\\.gz$", "", x, ignore.case = TRUE)
  x <- sub("\\.(fastq|fq|fasta|fa)$", "", x, ignore.case = TRUE)
  x <- sub("(_R?[12])$", "", x, ignore.case = TRUE)
  x
}

message("Reading oral metadata...")
metadata <- fread(paths$oral_metadata)
stopifnot("sample_id" %in% names(metadata))
oral_ids <- unique(normalize_sample_id(metadata$sample_id))

message("Reading BIOM sample and OTU IDs...")
sample_ids <- h5read(paths$biom, "/sample/ids")
otu_ids <- h5read(paths$biom, "/observation/ids")
normalized_biom_ids <- normalize_sample_id(sample_ids)

oral_sample_index <- which(normalized_biom_ids %in% oral_ids)
if (length(oral_sample_index) == 0) {
  stop("No oral metadata samples matched BIOM sample IDs after normalization.")
}

unmatched_oral_ids <- setdiff(oral_ids, normalized_biom_ids)
if (length(unmatched_oral_ids) > 0) {
  fwrite(
    data.table(sample_id = unmatched_oral_ids),
    file.path(paths$outdir, paste0(paths$prefix, "_oral_metadata_samples_not_in_biom.csv"))
  )
}

message("Reading sparse matrix pointers and indices...")
indptr <- h5read(paths$biom, "/sample/matrix/indptr")
indices <- h5read(paths$biom, "/sample/matrix/indices")
if (length(indptr) != length(sample_ids) + 1) {
  stop("BIOM sample matrix indptr length does not equal sample_count + 1.")
}

message("Counting OTU presence across oral samples...")
otu_presence_count <- integer(length(otu_ids))
for (j in oral_sample_index) {
  start <- indptr[j] + 1L
  end <- indptr[j + 1L]
  if (end >= start) {
    otu_idx <- indices[start:end] + 1L
    otu_presence_count[otu_idx] <- otu_presence_count[otu_idx] + 1L
  }
}

oral_sample_count <- length(oral_sample_index)
otu_freq <- data.table(
  otu_id = as.character(otu_ids),
  oral_sample_presence_count = otu_presence_count,
  oral_sample_count = oral_sample_count,
  oral_sample_frequency = otu_presence_count / oral_sample_count
)

fwrite(otu_freq, file.path(paths$outdir, paste0(paths$prefix, "_table.csv")))

summary_table <- data.table(
  metric = c(
    "oral_metadata_unique_samples",
    "oral_samples_matched_in_biom",
    "oral_metadata_samples_not_in_biom",
    "total_otus",
    "otus_present_in_at_least_one_oral_sample",
    "otus_absent_from_oral_samples"
  ),
  value = c(
    length(oral_ids),
    oral_sample_count,
    length(unmatched_oral_ids),
    length(otu_ids),
    sum(otu_presence_count > 0),
    sum(otu_presence_count == 0)
  )
)
fwrite(summary_table, file.path(paths$outdir, paste0(paths$prefix, "_summary.csv")))

plot_data <- otu_freq[oral_sample_presence_count > 0]
fwrite(plot_data, file.path(paths$outdir, paste0(paths$prefix, "_present_only_table.csv")))

p <- ggplot(plot_data, aes(oral_sample_frequency)) +
  geom_histogram(bins = 200, fill = "#4C78A8", color = "white", linewidth = 0.1) +
  scale_x_log10(labels = percent_format(accuracy = 0.1)) +
  scale_y_log10(labels = comma) +
  labs(
    title = "OTU prevalence across oral samples",
    subtitle = paste0(
      "Oral samples matched in BIOM: ", comma(oral_sample_count),
      "; OTUs present in oral samples: ", comma(sum(otu_presence_count > 0)),
      " / ", comma(length(otu_ids)), "; zero-frequency OTUs excluded"
    ),
    x = "Fraction of oral samples where OTU is present",
    y = "Number of OTUs (log10)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

png_file <- file.path(paths$outdir, paste0(paths$prefix, "_histogram.png"))
pdf_file <- file.path(paths$outdir, paste0(paths$prefix, "_histogram.pdf"))
ggsave(png_file, p, width = 9, height = 5.5, dpi = 300)
ggsave(pdf_file, p, width = 9, height = 5.5)

message("Done.")
message("Oral metadata unique samples: ", length(oral_ids))
message("Oral samples matched in BIOM: ", oral_sample_count)
message("Unmatched oral metadata samples: ", length(unmatched_oral_ids))
message("OTUs present in >=1 oral sample: ", sum(otu_presence_count > 0), "/", length(otu_ids))
message("Table: ", file.path(paths$outdir, paste0(paths$prefix, "_table.csv")))
message("Summary: ", file.path(paths$outdir, paste0(paths$prefix, "_summary.csv")))
message("PNG: ", png_file)
message("PDF: ", pdf_file)
