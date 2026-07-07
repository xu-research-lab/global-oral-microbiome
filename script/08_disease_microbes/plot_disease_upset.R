#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(UpSetR)
})

BASE <- "/beegfs/huangxiaochang/global_oral/new_analysis/maaslin3"
RESULTS_DIR <- file.path(BASE, "disease_results")
OUT_DIR <- file.path(BASE, "disease_shared_plots")
Q_THRESHOLD <- 0.05

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

disease_labels <- c(
  "Dental_Caries" = "Dental Caries",
  "Head_and_Neck_Cancer" = "H&N Cancer",
  "HIV" = "HIV",
  "OSCC" = "OSCC",
  "Periodontal_Disease" = "Periodontal Dz",
  "RA" = "RA",
  "Breast_Cancer" = "Breast Cancer",
  "CRC" = "CRC",
  "IBD" = "IBD",
  "Lung_Cancer" = "Lung Cancer"
)

clean_genus <- function(x) {
  x <- trimws(as.character(x))
  has_genus <- grepl("g__", x)
  extracted <- regmatches(x, regexpr("g__[^;]*", x))
  x[has_genus] <- sub("^g__", "", extracted)
  x <- sub("_.+$", "", x)
  ifelse(is.na(x) | x == "" | x == "g__", NA_character_, x)
}

files <- list.files(RESULTS_DIR, pattern = "^all_results.tsv$",
                    recursive = TRUE, full.names = TRUE)
if (length(files) == 0) stop("No all_results.tsv files found: ", RESULTS_DIR)

disease_sets <- list()
for (f in sort(files)) {
  disease_dir <- basename(dirname(f))
  disease <- ifelse(disease_dir %in% names(disease_labels),
                    disease_labels[[disease_dir]], gsub("_", " ", disease_dir))
  res <- read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
  genera <- res %>%
    filter(metadata == "Disease_Status",
           !is.na(qval_individual),
           qval_individual <= Q_THRESHOLD) %>%
    mutate(Genus = clean_genus(feature)) %>%
    filter(!is.na(Genus)) %>%
    pull(Genus) %>%
    unique()
  if (length(genera) > 0) disease_sets[[disease]] <- genera
}

if (length(disease_sets) == 0) {
  stop("No significant disease-associated genera found at q <= ", Q_THRESHOLD)
}

upset_order <- names(sort(sapply(disease_sets, length), decreasing = TRUE))
upset_input <- fromList(disease_sets[upset_order])

write.csv(
  data.frame(Disease = upset_order,
             SignificantGenera = sapply(disease_sets[upset_order], length),
             row.names = NULL),
  file.path(OUT_DIR, "upset_set_sizes.csv"),
  row.names = FALSE
)

pdf(file.path(OUT_DIR, "disease_genera_upset.pdf"), width = 13, height = 7.2)
upset(
  upset_input,
  nsets = length(upset_order),
  nintersects = 80,
  order.by = "freq",
  mb.ratio = c(0.55, 0.45),
  matrix.color = "#2C7FB8",
  main.bar.color = "#2C7FB8",
  sets.bar.color = "#41B6C4",
  mainbar.y.label = "Shared genera count",
  sets.x.label = "Significant genera per disease",
  point.size = 3.6,
  line.size = 1.1,
  text.scale = c(1.45, 1.25, 1.1, 1.1, 0.95, 0.95)
)
dev.off()

png(file.path(OUT_DIR, "disease_genera_upset.png"),
    width = 13, height = 7.2, units = "in", res = 300)
upset(
  upset_input,
  nsets = length(upset_order),
  nintersects = 80,
  order.by = "freq",
  mb.ratio = c(0.55, 0.45),
  matrix.color = "#2C7FB8",
  main.bar.color = "#2C7FB8",
  sets.bar.color = "#41B6C4",
  mainbar.y.label = "Shared genera count",
  sets.x.label = "Significant genera per disease",
  point.size = 3.6,
  line.size = 1.1,
  text.scale = c(1.45, 1.25, 1.1, 1.1, 0.95, 0.95)
)
dev.off()

cat("Saved UpSet plot to ", OUT_DIR, "\n", sep = "")
print(data.frame(Disease = upset_order,
                 SignificantGenera = sapply(disease_sets[upset_order], length),
                 row.names = NULL))
