#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
})

BASE <- "/beegfs/huangxiaochang/global_oral/new_analysis/maaslin3"
RESULT_FILE <- file.path(BASE, "site_stratified_results", "Saliva", "all_results.tsv")
OUT_DIR <- file.path(BASE, "site_stratified_saliva_heatmaps")
TAXONOMY <- "/beegfs/huangxiaochang/global_oral/taxonomy.tsv"
Q_THRESHOLD <- 0.1
VARIABLES <- c("Age_Group", "Sex", "Geography", "Region_custom")

refs <- c(
  Age_Group = "Reference: 18-60 years",
  Sex = "Reference: F",
  Geography = "Reference: North America",
  Region_custom = "Reference: V4"
)

titles <- c(
  Age_Group = "Age group",
  Sex = "Sex",
  Geography = "Geography",
  Region_custom = "16S region"
)

value_orders <- list(
  Age_Group = c("0-6 months", "6-36 months", "3-12 years", "12-18 years", "13-18 years", "18-60 years", ">60 years", "＞60 years", "Unknown"),
  Sex = c("F", "M", "Unknown"),
  Geography = c("North America", "Asia", "Europe", "South America", "Oceania", "Africa", "Unknown"),
  Region_custom = c("V4", "V1-V2", "V1-V3", "V2", "V2-V3", "V3", "V3-V4", "V3-V5", "V4-V5", "V5-V7", "V6-V8", "Unknown")
)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

clean_genus <- function(x) {
  x <- trimws(as.character(x))
  has_genus <- grepl("g__", x)
  extracted <- regmatches(x, regexpr("g__[^;]*", x))
  x[has_genus] <- sub("^g__", "", extracted)
  x <- sub("_.+$", "", x)
  ifelse(is.na(x) | x == "" | x == "g__", NA_character_, x)
}

extract_rank <- function(taxon, prefix) {
  parts <- strsplit(as.character(taxon), ";", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  hit <- parts[startsWith(parts, prefix)]
  if (length(hit) == 0 || hit[1] == prefix) return(NA_character_)
  hit[1]
}

build_family_map <- function() {
  if (!file.exists(TAXONOMY)) return(setNames(character(), character()))
  tax <- read.delim(TAXONOMY, stringsAsFactors = FALSE, check.names = FALSE)
  tax_col <- if ("Taxon" %in% colnames(tax)) "Taxon" else colnames(tax)[2]
  fam <- vapply(tax[[tax_col]], extract_rank, character(1), prefix = "f__")
  gen <- vapply(tax[[tax_col]], extract_rank, character(1), prefix = "g__")
  gen_clean <- clean_genus(gen)
  df <- data.frame(Genus = gen_clean, Family = fam, stringsAsFactors = FALSE) %>%
    filter(!is.na(Genus), !is.na(Family), Family != "f__") %>%
    count(Genus, Family, name = "n") %>%
    group_by(Genus) %>%
    arrange(desc(n), Family, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()
  setNames(df$Family, df$Genus)
}

make_axis_labels <- function(genera, family_map) {
  family <- family_map[genera]
  family[is.na(family)] <- "f__"
  genus <- paste0("g__", genera)
  family_width <- max(nchar(family), na.rm = TRUE) + 2
  labels <- paste0(str_pad(family, width = family_width, side = "right"), genus)
  names(labels) <- genera
  labels
}

prepare_variable <- function(results, variable) {
  var_data <- results %>%
    filter(metadata == variable, !is.na(qval_individual), !is.na(coef), !is.na(Genus)) %>%
    group_by(Genus, value) %>%
    arrange(qval_individual, desc(abs(coef)), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      q_plot = ifelse(qval_individual <= 0, 1e-300, qval_individual),
      signed_logq = pmin(-log10(q_plot), 10) * sign(coef),
      sig_label = ifelse(qval_individual < Q_THRESHOLD, ifelse(coef > 0, "+", "-"), ""),
      Direction = ifelse(coef > 0, "+", "-")
    )

  sig_genera <- var_data %>%
    filter(qval_individual < Q_THRESHOLD) %>%
    distinct(Genus) %>%
    pull(Genus)
  if (length(sig_genera) == 0) return(data.frame())

  var_data %>% filter(Genus %in% sig_genera)
}

order_genera <- function(data) {
  data %>%
    mutate(sig_score = ifelse(qval_individual < Q_THRESHOLD, signed_logq, 0)) %>%
    group_by(Genus) %>%
    summarise(
      mean_signed = mean(sig_score, na.rm = TRUE),
      min_q = min(qval_individual, na.rm = TRUE),
      max_abs = max(abs(coef), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(mean_signed, min_q, desc(max_abs), Genus) %>%
    pull(Genus)
}

plot_variable <- function(results, variable, family_map) {
  data <- prepare_variable(results, variable)
  if (nrow(data) == 0) {
    cat("Skip ", variable, ": no significant Saliva associations at q < ", Q_THRESHOLD, "\n", sep = "")
    return(NULL)
  }

  genus_order <- order_genera(data)
  observed_values <- unique(data$value)
  preferred <- value_orders[[variable]]
  value_levels <- c(preferred[preferred %in% observed_values], sort(setdiff(observed_values, preferred)))

  grid <- expand.grid(
    Genus = genus_order,
    value = value_levels,
    stringsAsFactors = FALSE
  ) %>%
    left_join(
      data %>% select(Genus, value, feature, coef, qval_individual, signed_logq, sig_label, Direction, model, N, N_not_zero),
      by = c("Genus", "value")
    ) %>%
    mutate(
      signed_logq = ifelse(is.na(signed_logq), 0, signed_logq),
      sig_label = ifelse(is.na(sig_label), "", sig_label),
      Genus = factor(Genus, levels = rev(genus_order)),
      value = factor(value, levels = value_levels)
    )

  axis_labels <- make_axis_labels(genus_order, family_map)

  write.table(
    data %>% arrange(qval_individual, desc(abs(coef))),
    file.path(OUT_DIR, paste0("saliva_site_", variable, "_significant_source.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    grid,
    file.path(OUT_DIR, paste0("saliva_site_", variable, "_heatmap_data.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    data.frame(Genus = genus_order, AxisLabel = axis_labels[genus_order], stringsAsFactors = FALSE),
    file.path(OUT_DIR, paste0("saliva_site_", variable, "_feature_labels.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  col_width <- ifelse(variable == "Geography", 0.82, ifelse(variable == "Region_custom", 0.58, 0.52))
  width <- max(ifelse(variable == "Geography", 6.8, 4.8), length(value_levels) * col_width + 2.4)
  height <- max(4.5, length(genus_order) * 0.17 + 1.3)

  p <- ggplot(grid, aes(x = value, y = Genus, fill = signed_logq)) +
    geom_tile(color = "grey86", linewidth = 0.14) +
    geom_text(aes(label = sig_label), size = 2.25, fontface = "bold") +
    scale_y_discrete(labels = axis_labels) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, limits = c(-10, 10), name = "signed\n-log10(q)"
    ) +
    labs(
      title = titles[[variable]],
      subtitle = paste0(refs[[variable]], "; q < ", Q_THRESHOLD),
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8.2, face = "bold"),
      axis.text.y = element_text(size = 8.4, family = "mono", hjust = 0),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.key.width = unit(0.75, "cm"),
      legend.key.height = unit(0.25, "cm"),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6.5),
      plot.title = element_text(hjust = 0.5, size = 10.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 7.5),
      plot.margin = margin(6, 6, 6, 6)
    )

  pdf_file <- file.path(OUT_DIR, paste0("saliva_site_", variable, "_heatmap.pdf"))
  png_file <- file.path(OUT_DIR, paste0("saliva_site_", variable, "_heatmap.png"))
  ggsave(pdf_file, p, width = width, height = height, limitsize = FALSE)
  ggsave(png_file, p, width = width, height = height, dpi = 300, limitsize = FALSE)

  cat(sprintf("Saved %s: %d genera, %d levels\n", variable, length(genus_order), length(value_levels)))
  invisible(data.frame(
    Variable = variable,
    SignificantGenera = length(unique(data$Genus[data$qval_individual < Q_THRESHOLD])),
    SignificantAssociations = sum(data$qval_individual < Q_THRESHOLD, na.rm = TRUE),
    Levels = paste(value_levels, collapse = ";"),
    PDF = pdf_file,
    PNG = png_file,
    stringsAsFactors = FALSE
  ))
}

if (!file.exists(RESULT_FILE)) stop("Missing result file: ", RESULT_FILE)
results <- read.delim(RESULT_FILE, stringsAsFactors = FALSE, check.names = FALSE)
results$Genus <- clean_genus(results$feature)
family_map <- build_family_map()

summaries <- lapply(VARIABLES, function(v) plot_variable(results, v, family_map))
summaries <- bind_rows(summaries)
write.csv(summaries, file.path(OUT_DIR, "saliva_site_heatmap_summary.csv"), row.names = FALSE)

all_sig <- results %>%
  filter(metadata %in% VARIABLES, !is.na(qval_individual), qval_individual < Q_THRESHOLD, !is.na(Genus)) %>%
  group_by(metadata, Genus, value) %>%
  arrange(qval_individual, desc(abs(coef)), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(metadata, qval_individual, desc(abs(coef)))
write.table(all_sig, file.path(OUT_DIR, "saliva_site_all_significant_q0.1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Done. Output: ", OUT_DIR, "\n", sep = "")
