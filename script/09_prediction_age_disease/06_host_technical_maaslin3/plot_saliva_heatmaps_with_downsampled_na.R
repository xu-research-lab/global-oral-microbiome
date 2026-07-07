#!/usr/bin/env Rscript
# ============================================================================
# Same as plot_saliva_heatmaps.R, but replaces North America results with the
# stratified downsampled North America results.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

BASE <- "/beegfs/huangxiaochang/global_oral/new_analysis/maaslin3"
RESULTS_DIR <- file.path(BASE, "results")
DS_NA_DIR <- file.path(BASE, "saliva_sex_downsampled_na", "results", "Saliva_North_America")
OUT_DIR <- file.path(BASE, "saliva_heatmaps_downsampled_na")
TAXONOMY <- "/beegfs/huangxiaochang/global_oral/taxonomy.tsv"
Q_THRESHOLD <- 0.1

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

geography_order <- c("North America", "Asia", "Europe", "South America", "Oceania", "Africa")
variables <- c("Age_Group", "Sex", "Region_custom")

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

collect_saliva_results <- function() {
  files <- list.files(RESULTS_DIR, pattern = "^all_results.tsv$", recursive = TRUE, full.names = TRUE)
  files <- files[grepl("/Saliva_", files)]
  if (length(files) == 0) stop("No Saliva all_results.tsv files found under ", RESULTS_DIR)

  rows <- list()
  for (f in sort(files)) {
    group <- basename(dirname(f))
    geography <- sub("^Saliva_", "", group)
    geography <- gsub("_", " ", geography)

    # Replace North America with downsampled results
    if (geography == "North America") {
      f <- file.path(DS_NA_DIR, "all_results.tsv")
      if (!file.exists(f)) stop("Downsampled North America results not found at ", f)
      group <- "Saliva_North_America_Downsampled"
    }

    res <- read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
    res$Group <- group
    res$Site <- "Saliva"
    res$Geography <- geography
    res$FacetX <- geography
    res$Genus <- clean_genus(res$feature)
    rows[[group]] <- res
  }
  bind_rows(rows) %>% filter(!is.na(Genus), Geography != "Unknown")
}

order_features <- function(data) {
  data %>%
    group_by(Genus) %>%
    summarize(
      min_q = min(qval_individual, na.rm = TRUE),
      mean_coef = mean(coef, na.rm = TRUE),
      first_geo = min(match(Geography, geography_order), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      sig = min_q < Q_THRESHOLD,
      abs_coef = abs(mean_coef)
    ) %>%
    arrange(desc(sig), first_geo, min_q, desc(abs_coef), Genus) %>%
    pull(Genus)
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

variable_title <- function(variable) {
  switch(variable,
         "Age_Group" = "Age-Associated",
         "Sex" = "Sex-Associated",
         "Region_custom" = "Region-Associated",
         variable)
}

reference_text <- function(variable) {
  switch(variable,
         "Age_Group" = "Reference: 18-60 years",
         "Sex" = "Reference: Female",
         "Region_custom" = "Reference: V4",
         "Reference level")
}

plot_variable <- function(all_results, variable, family_map,
                         consistent = FALSE,
                         output_suffix = "heatmap",
                         title_suffix = "") {
  data <- all_results %>%
    filter(metadata == variable, !is.na(qval_individual)) %>%
    group_by(Genus, Geography, value) %>%
    arrange(qval_individual, desc(abs(coef)), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      q_plot = ifelse(qval_individual <= 0, 1e-300, qval_individual),
      signed_logq = pmin(-log10(q_plot), 10) * sign(coef),
      sig_label = ifelse(qval_individual < Q_THRESHOLD, ifelse(coef > 0, "+", "-"), ""),
      Direction = ifelse(coef > 0, "+", "-")
    )

  sig_rows <- data %>% filter(qval_individual < Q_THRESHOLD)
  if (nrow(sig_rows) == 0) return(NULL)

  if (!consistent) {
    kept_pairs <- sig_rows %>% distinct(value, Geography)
    kept_genera <- sig_rows %>% distinct(value, Genus)
    data <- data %>%
      inner_join(kept_pairs, by = c("value", "Geography")) %>%
      inner_join(kept_genera, by = c("value", "Genus"))
  } else {
    kept <- sig_rows %>%
      distinct(value, Genus, Geography, Direction) %>%
      count(value, Genus, Direction, name = "SignificantGeographies") %>%
      filter(SignificantGeographies >= 2)
    if (nrow(kept) == 0) return(NULL)
    data <- data %>%
      inner_join(kept, by = c("value", "Genus", "Direction"))
  }

  if (nrow(data) == 0) return(NULL)

  genus_order <- order_features(data)
  value_levels <- unique(data$value[order(data$value)])

  grid_parts <- list()
  for (v in value_levels) {
    facet_data <- data %>% filter(value == v)
    facet_geos <- geography_order[geography_order %in% unique(facet_data$Geography)]
    facet_genera <- genus_order[genus_order %in% unique(facet_data$Genus)]
    if (length(facet_geos) == 0 || length(facet_genera) == 0) next
    grid_parts[[as.character(v)]] <- expand.grid(
      Genus = facet_genera,
      Geography = facet_geos,
      value = v,
      stringsAsFactors = FALSE
    )
  }

  grid <- bind_rows(grid_parts) %>%
    left_join(
      data %>% select(Genus, Geography, value, coef, qval_individual, signed_logq, sig_label),
      by = c("Genus", "Geography", "value")
    ) %>%
    mutate(
      signed_logq = ifelse(is.na(signed_logq), 0, signed_logq),
      sig_label = ifelse(is.na(sig_label), "", sig_label),
      Genus = factor(Genus, levels = rev(genus_order)),
      Geography = factor(Geography, levels = geography_order),
      value = factor(value, levels = value_levels)
    )

  axis_labels <- make_axis_labels(genus_order, family_map)
  write.table(grid, file.path(OUT_DIR, paste0("saliva_", variable, "_", output_suffix, "_data.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(
    data.frame(Genus = genus_order, AxisLabel = axis_labels[genus_order], stringsAsFactors = FALSE),
    file.path(OUT_DIR, paste0("saliva_", variable, "_", output_suffix, "_feature_labels.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  max_geos_per_facet <- grid %>% distinct(value, Geography) %>% count(value) %>% pull(n) %>% max()
  width <- max(6.2, length(value_levels) * max_geos_per_facet * 0.30)
  height <- max(5.2, length(genus_order) * 0.18)

  p <- ggplot(grid, aes(x = Geography, y = Genus, fill = signed_logq)) +
    geom_tile(color = "grey85", linewidth = 0.14) +
    geom_text(aes(label = sig_label), size = 2.2, fontface = "bold") +
    facet_wrap(~ value, nrow = 1, scales = "free_x") +
    scale_y_discrete(labels = axis_labels) +
    scale_x_discrete(drop = TRUE) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#D73027",
      midpoint = 0, limits = c(-10, 10), name = "signed\n-log10(q)"
    ) +
    labs(
      title = paste0("Saliva ", variable_title(variable), title_suffix,
                     " (North America Downsampled)"),
      subtitle = paste0(reference_text(variable),
                        "; q < ", Q_THRESHOLD,
                        "; +/− indicates coefficient direction; NA n=2470"),
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7.5, face = "bold"),
      axis.text.y = element_text(size = 8.2, family = "mono", hjust = 0),
      strip.text = element_text(size = 8.2, face = "bold"),
      panel.grid = element_blank(),
      panel.spacing = unit(0.2, "lines"),
      legend.position = "bottom",
      legend.key.width = unit(0.75, "cm"),
      legend.key.height = unit(0.25, "cm"),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6.5),
      plot.title = element_text(hjust = 0.5, size = 11.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 8.0),
      plot.margin = margin(6, 6, 6, 6)
    )

  pdf_file <- file.path(OUT_DIR, paste0("saliva_", variable, "_", output_suffix, ".pdf"))
  png_file <- file.path(OUT_DIR, paste0("saliva_", variable, "_", output_suffix, ".png"))
  ggsave(pdf_file, p, width = width, height = height, limitsize = FALSE)
  ggsave(png_file, p, width = width, height = height, limitsize = FALSE, dpi = 300)
  cat(sprintf("Saved %s %s: %d genera, %d facets (%.1f x %.1f in)\n",
              variable, output_suffix, length(genus_order), length(value_levels), width, height))
}

cat("Loading Saliva MaAsLin3 stratified results (with downsampled North America)...\n")
all_results <- collect_saliva_results()
family_map <- build_family_map()

summary_df <- all_results %>%
  filter(metadata %in% variables, !is.na(qval_individual), qval_individual < Q_THRESHOLD) %>%
  group_by(metadata, Geography) %>%
  summarize(SignificantGenera = n_distinct(Genus), .groups = "drop") %>%
  arrange(metadata, Geography)
write.csv(summary_df, file.path(OUT_DIR, "saliva_significant_genera_counts.csv"), row.names = FALSE)
print(summary_df)

for (variable in variables) {
  plot_variable(
    all_results, variable, family_map,
    consistent = FALSE,
    output_suffix = "heatmap",
    title_suffix = ""
  )
  plot_variable(
    all_results, variable, family_map,
    consistent = TRUE,
    output_suffix = "consistent_heatmap",
    title_suffix = " consistent"
  )
}

cat("Done. Outputs in ", OUT_DIR, "\n", sep = "")
