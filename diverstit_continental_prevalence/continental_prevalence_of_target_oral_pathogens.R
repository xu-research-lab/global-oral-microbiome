suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

find_script_dir <- function() {
  if (interactive()) return(normalizePath(getwd()))
  args <- commandArgs(trailingOnly = FALSE)
  f <- args[grep("^--file=", args)]
  if (length(f) == 1) return(dirname(normalizePath(sub("^--file=", "", f))))
  normalizePath(getwd())
}
script_dir <- find_script_dir()

PROJECT_ROOT <- normalizePath(file.path(script_dir, ".."))
DATA_DIR     <- file.path(PROJECT_ROOT, "data", "metadata-and-biom-v2-qc")

BASE_DIR  <- PROJECT_ROOT
PANEL_DIR <- Sys.getenv(
  "FIG3_PANEL_DIR",
  unset = file.path(script_dir, "results", "Fig3")
)
INPUT_CACHE_DIR <- file.path(PANEL_DIR, ".input_cache_0922")

INPUT_BIOM <- file.path(DATA_DIR, "biom_oral_qc_dedup.biom")
INPUT_TAXONOMY <- file.path(DATA_DIR, "taxonomy.tsv")

metadata_path <- Sys.getenv("FIG1_METADATA",
                            unset = file.path(DATA_DIR,
                                              "metadata_oral_qc_dedup.tsv"))
abund_file <- Sys.getenv(
  "FIG1_PHYLUM_TABLE",
  unset = file.path(INPUT_CACHE_DIR, "genus-relative-table_0922.tsv")
)

species_file <- Sys.getenv(
  "FIG1_SPECIES_TABLE",
  unset = file.path(INPUT_CACHE_DIR, "species-relative-table_0922.tsv")
)

PREV_THRESH <- 0.0001
TOP_N       <- 7L
FIG_DPI     <- 600L
BASE_SZ     <- 18
OUT_PREFIX  <- "Fig3"

SITE_DISPLAY <- c("Saliva", "Tooth", "Oral Mucosa", "Oral Swab", "Tongue",
                  "Gingiva")
SITE_LOWER   <- tolower(SITE_DISPLAY)

CONTINENTS <- c("North America", "South America", "Europe", "Asia", "Africa",
                "Oceania")

NAME_TO_ABBR <- c("North America" = "NA", "South America" = "SA",
                  "Europe" = "EU", "Asia" = "AS", "Africa" = "AF",
                  "Oceania" = "OC")

CONT_ABBR <- c("AS", "EU", "AF", "OC", "NA", "SA")

CONTINENT_FACET_ABBR <- c("NA", "EU", "AS", "OC", "AF", "SA")

AGE_ORDER <- c("0-6 months", "6-36 months", "3-12 years", "12-18 years",
               "18-60 years", ">60 years")
AGE_SHORT <- c("0-6 mo", "6-36 mo", "3-12 yr", "12-18 yr", "18-60 yr", ">60 yr")

GENUS_TOP7 <- c("Streptococcus", "Veillonella_A", "Prevotella", "Rothia",
                "Granulicatella", "Haemophilus_D", "Neisseria")

LOW_FILL  <- "#FFFFFF"
HIGH_FILL <- "#084594"
NA_FILL   <- "#E8E8E8"

pick_meta_col <- function(md, preferred, fallback = NULL) {
  if (preferred %in% names(md)) return(preferred)
  if (!is.null(fallback) && fallback %in% names(md)) return(fallback)
  preferred
}

short_taxon_label <- function(x) sub("_[0-9]+$", "", x)

to_title <- function(x) {
  vapply(strsplit(x, " "), function(parts) {
    paste0(toupper(substr(parts, 1, 1)), substr(parts, 2, nchar(parts)),
           collapse = " ")
  }, character(1))
}

fmt_pct_label <- function(x) sprintf("%.0f", round(x * 100))

ensure_prepared_inputs <- function() {
  input_files <- c(INPUT_BIOM, metadata_path, INPUT_TAXONOMY)
  missing_inputs <- input_files[!file.exists(input_files)]
  if (length(missing_inputs) > 0) {
    stop("Figure 3 input not found: ", paste(missing_inputs, collapse = ", "))
  }

  newest_input <- max(file.info(input_files)$mtime, na.rm = TRUE)
  stale <- function(path) {
    !file.exists(path) || file.info(path)$mtime < newest_input
  }
  if (!stale(abund_file) && !stale(species_file)) return(invisible(NULL))

  prep_script <- file.path(script_dir, "figure3_prepare_inputs_0922.py")
  if (!file.exists(prep_script)) {
    stop("input conversion script not found: ", prep_script)
  }

  python_bin <- Sys.getenv("FIG3_PYTHON", unset = "")
  if (!nzchar(python_bin)) {
    found <- Sys.which(c("python3", "python"))
    found <- found[nzchar(found)]
    if (length(found) == 0) stop("Python is required to convert the BIOM input")
    python_bin <- unname(found[1])
  }

  cat("-- rebuilding Figure 3 input cache --\n")
  status <- system2(
    python_bin,
    args = c(
      shQuote(prep_script),
      "--biom", shQuote(INPUT_BIOM),
      "--metadata", shQuote(metadata_path),
      "--taxonomy", shQuote(INPUT_TAXONOMY),
      "--out-dir", shQuote(INPUT_CACHE_DIR)
    )
  )
  if (!identical(status, 0L)) stop("Figure 3 input conversion failed")
  invisible(NULL)
}

prepare_data <- function() {
  ensure_prepared_inputs()
  if (!file.exists(metadata_path)) stop("metadata not found: ", metadata_path)
  if (!file.exists(abund_file))    stop("genus table not found: ", abund_file)

  md <- data.table::fread(metadata_path, header = TRUE, sep = "\t")
  col_site <- pick_meta_col(md, "site_cleaned", "site")
  col_ageg <- pick_meta_col(md, "age_cleaned_group", "age_group")
  col_cont <- pick_meta_col(md, "continent_group", "continent")

  rel <- data.table::fread(abund_file, header = TRUE, sep = "\t", skip = 1)
  data.table::setnames(rel, 1, "Taxon")
  sample_cols <- setdiff(names(rel), "Taxon")
  rel_long <- data.table::melt(rel, id.vars = "Taxon",
                               measure.vars = sample_cols,
                               variable.name = "sample_id",
                               value.name = "Abundance")
  rel_long[, sample_id := as.character(sample_id)]

  rel_long[, Genus := Taxon]
  rel_long <- rel_long[!is.na(Genus) & Genus != ""]
  rel_long[, Short := short_taxon_label(Genus)]
  rel_long[, Present := !is.na(Abundance) & Abundance >= PREV_THRESH]
  rm(rel)

  total_samples <- data.table::uniqueN(rel_long$sample_id)
  rank_dt <- rel_long[Present == TRUE,
                      .(Detected = data.table::uniqueN(sample_id)),
                      by = Short]
  rank_dt[, Prevalence := Detected / total_samples]
  rank_dt <- rank_dt[order(-Prevalence, Short)]
  top9 <- head(rank_dt$Short, TOP_N)

  meta_sub <- md[, .(sample_id,
                     site = get(col_site),
                     cont = get(col_cont),
                     age  = get(col_ageg))]
  meta_sub[, sample_id := as.character(sample_id)]
  meta_sub <- unique(meta_sub, by = "sample_id")

  clean_site_age_cont <- function(d) {
    d[, site := tolower(trimws(as.character(site)))]
    d[, cont := trimws(as.character(cont))]
    d[, age  := gsub("＞", ">", as.character(age))]
    bad_site <- c("other", "others", "other (specify)", "others (specify)")
    d <- d[!is.na(site) & site != "" & !(site %in% bad_site)]
    d
  }
  meta_clean <- clean_site_age_cont(copy(meta_sub))
  present_long <- merge(rel_long[Present == TRUE, .(sample_id, Short)],
                        meta_sub, by = "sample_id", all.x = FALSE)
  present_long <- clean_site_age_cont(present_long)
  present_long[, site_display := factor(to_title(site), levels = SITE_DISPLAY)]

  N_by_site_cont <- meta_clean[!is.na(cont) & cont != "",
                               .(N = .N), by = .(site, cont)]
  N_by_site_age  <- meta_clean[!is.na(age) & age != "",
                               .(N = .N), by = .(site, age)]

  count_present <- function(pd, strata_col) {
    cnt <- pd[, .(N_present = data.table::uniqueN(sample_id)),
              by = c("site", strata_col, "Short")]
    data.table::setnames(cnt, strata_col, "Stratum")
    cnt
  }

  pd_a_obs <- present_long[cont %in% CONTINENTS]
  pd_a_obs[, Contin := NAME_TO_ABBR[cont]]
  cnt_a <- count_present(pd_a_obs, "Contin")
  Nc <- N_by_site_cont[cont %in% CONTINENTS]
  Nc[, Contin := NAME_TO_ABBR[cont]]
  grid_a <- data.table::CJ(site = SITE_LOWER, Stratum = CONT_ABBR, Short = top9)
  pd_a <- merge(grid_a, cnt_a, by = c("site", "Stratum", "Short"), all.x = TRUE)
  pd_a <- merge(pd_a, Nc[, .(site, Stratum = Contin, N)],
                by = c("site", "Stratum"), all.x = TRUE)
  pd_a[, Prevalence := fifelse(N > 0, fifelse(is.na(N_present), 0, N_present) / N,
                              NA_real_)]
  pd_a[, site_display := factor(to_title(site), levels = SITE_DISPLAY)]

  pd_b_obs <- present_long[age %in% AGE_ORDER]
  pd_b_obs[, Ageg := factor(age, levels = AGE_ORDER, labels = AGE_SHORT)]
  cnt_b <- count_present(pd_b_obs, "Ageg")
  Nb <- N_by_site_age[age %in% AGE_ORDER]
  Nb[, Ageg := factor(age, levels = AGE_ORDER, labels = AGE_SHORT)]
  grid_b <- data.table::CJ(site = SITE_LOWER, Stratum = AGE_SHORT, Short = top9)
  pd_b <- merge(grid_b, cnt_b, by = c("site", "Stratum", "Short"), all.x = TRUE)
  pd_b <- merge(pd_b, Nb[, .(site, Stratum = Ageg, N)],
                by = c("site", "Stratum"), all.x = TRUE)
  pd_b[, Prevalence := fifelse(N > 0, fifelse(is.na(N_present), 0, N_present) / N,
                              NA_real_)]
  pd_b[, site_display := factor(to_title(site), levels = SITE_DISPLAY)]

  N_by_continent <- meta_clean[cont %in% CONTINENTS, .(Sample_N = .N),
                               by = .(Continent = cont)]
  continent_count_obs <- present_long[cont %in% CONTINENTS,
                                      .(Count = data.table::uniqueN(sample_id)),
                                      by = .(Continent = cont, Short)]
  continent_count_obs[, Rank := data.table::frank(-Count, ties.method = "first"),
                      by = Continent]
  pd_continent_count <- merge(
    continent_count_obs[, .(Continent, Short, Count, Rank)],
    N_by_continent, by = "Continent", all.x = TRUE
  )
  pd_continent_count[, Prevalence := Count / Sample_N]
  pd_continent_count[, ContLabel := factor(NAME_TO_ABBR[Continent],
                                           levels = CONTINENT_FACET_ABBR)]

  N_by_site <- meta_clean[, .(Sample_N = .N), by = .(Site = site)]
  site_count_obs <- present_long[, .(Count = data.table::uniqueN(sample_id)),
                                 by = .(Site = site, Short)]
  site_count_obs[, Rank := data.table::frank(-Count, ties.method = "first"),
                 by = Site]
  pd_site_count <- merge(
    site_count_obs[, .(Site, Short, Count, Rank)],
    N_by_site, by = "Site", all.x = TRUE
  )
  pd_site_count[, Prevalence := Count / Sample_N]
  pd_site_count[, SiteLabel := factor(to_title(Site), levels = SITE_DISPLAY)]

  prev_a_obs <- pd_a[!is.na(Prevalence)]
  prev_b_obs <- pd_b[!is.na(Prevalence)]

  list(rank_dt = rank_dt, top9 = top9, total_samples = total_samples,
       pd_a = pd_a, pd_b = pd_b,
       pd_continent_count = pd_continent_count,
       pd_site_count = pd_site_count,
       prev_a_obs = prev_a_obs, prev_b_obs = prev_b_obs,
       N_by_site_cont = N_by_site_cont, N_by_site_age = N_by_site_age,
       meta_clean = meta_clean)
}

prepare_species_data <- function(meta_clean) {
  if (!file.exists(species_file)) {
    cat("\n-- species-relative table not found:", species_file, "--\n")
    cat("   (run species_relative_prep.py first) -> species figure skipped --\n")
    return(NULL)
  }

  rel <- data.table::fread(species_file, header = TRUE, sep = "\t", skip = 1)
  data.table::setnames(rel, 1, "Species")
  sample_cols <- setdiff(names(rel), "Species")
  rel_long <- data.table::melt(rel, id.vars = "Species",
                               measure.vars = sample_cols,
                               variable.name = "sample_id",
                               value.name = "Abundance")
  rel_long[, sample_id := as.character(sample_id)]
  rel_long[, Present := !is.na(Abundance) & Abundance >= PREV_THRESH]
  rm(rel)

  total_samples <- data.table::uniqueN(rel_long$sample_id)
  rank_s <- rel_long[Present == TRUE,
                     .(Detected = data.table::uniqueN(sample_id)),
                     by = Species]
  rank_s[, Prevalence := Detected / total_samples]
  rank_s <- rank_s[order(-Prevalence, Species)]

  present_long <- merge(rel_long[Present == TRUE, .(sample_id, Species)],
                        meta_clean, by = "sample_id", all.x = FALSE)
  present_long[, Contin := NAME_TO_ABBR[cont]]
  cnt <- present_long[cont %in% CONTINENTS,
                      .(N_present = data.table::uniqueN(sample_id)),
                      by = .(Contin, Species)]
  Nc <- meta_clean[cont %in% CONTINENTS, .(N = .N), by = cont]
  Nc[, Contin := NAME_TO_ABBR[cont]]

  grid_s <- data.table::CJ(Stratum = CONT_ABBR, Species = rank_s$Species)
  pd <- merge(grid_s, cnt[, .(Stratum = Contin, Species, N_present)],
              by = c("Stratum", "Species"), all.x = TRUE)
  pd <- merge(pd, Nc[, .(Stratum = Contin, N)],
              by = "Stratum", all.x = TRUE)
  pd[, Prevalence := fifelse(N > 0, fifelse(is.na(N_present), 0, N_present) / N,
                            NA_real_)]
  pd[, Short := Species]

  list(pd = pd, rank_s = rank_s, total_samples = total_samples)
}

write_out <- function(plot, stem, width, height, dpi = FIG_DPI) {
  png_path <- file.path(PANEL_DIR, paste0(OUT_PREFIX, "_", stem, "_0922.png"))
  grDevices::png(png_path, width = width, height = height,
                 units = "in", res = dpi, bg = "white")
  print(plot)
  grDevices::dev.off()
  cat("  wrote", basename(png_path), "\n")

  jpeg_path <- file.path(PANEL_DIR, paste0(OUT_PREFIX, "_", stem, "_0922.jpg"))
  grDevices::jpeg(jpeg_path, width = width, height = height,
                  units = "in", res = dpi, quality = 95)
  print(plot)
  grDevices::dev.off()
  cat("  wrote", basename(jpeg_path), "\n")
}

build_heatmap <- function(pd, xorder, rot_x = FALSE, base = BASE_SZ,
                          phyla = GENUS_TOP7, facet_sites = TRUE,
                          italic_y = FALSE, label_size = 3.8,
                          axis_mult = 1.0, square_cells = FALSE,
                          strip_wrap_width = 8, wrap_y = FALSE) {
  pd <- copy(pd)
  pd[, Stratum := factor(as.character(Stratum), levels = xorder)]
  pd[, Phylum_f := factor(Short, levels = rev(phyla))]
  p <- ggplot(pd, aes(x = Stratum, y = Phylum_f, fill = Prevalence)) +
    geom_tile(color = "white", linewidth = 0.55) +
    geom_text(aes(label = ifelse(!is.na(Prevalence),
                                 paste0(fmt_pct_label(Prevalence), "%"),
                                 ""),
                  color = ifelse(!is.na(Prevalence) & Prevalence >= 0.5,
                                 "white", "#333333")),
              size = label_size, show.legend = FALSE) +
    scale_color_identity() +
    scale_y_discrete(labels = if (wrap_y) function(x) ifelse(nchar(x) > 30L, gsub(" ", "\n", x), x) else waiver()) +
    scale_fill_gradient(low = LOW_FILL, high = HIGH_FILL, limits = c(0, 1),
                        breaks = c(0, 0.25, 0.5, 0.75, 1),
                        labels = scales::label_percent(accuracy = 1),
                        na.value = NA_FILL, name = "Prevalence") +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = base) +
    theme(
      axis.text         = element_text(color = "#333333"),
      axis.text.x       = element_text(size = rel(1.08 * axis_mult),
                                       face = "bold",
                                       angle = if (rot_x) 40 else 0,
                                       hjust = if (rot_x) 1 else 0.5,
                                       vjust = if (rot_x) 1 else 0.5),
      axis.text.y       = element_text(size = rel(1.18 * axis_mult),
                                       face = if (italic_y) "bold.italic" else "bold"),
      axis.ticks        = element_blank(),
      axis.line         = element_blank(),
      panel.border      = element_rect(color = "#BEBEBE", fill = NA,
                                       linewidth = 0.5),
      panel.spacing     = unit(0.6, "lines"),
      panel.grid        = element_blank(),
      strip.text        = element_text(size = rel(1.3), face = "bold",
                                       color = "#333333",
                                       lineheight = 0.95,
                                       margin = margin(b = 6, t = 3)),
      strip.background  = element_blank(),
      legend.position   = "right",
      legend.title      = element_text(size = rel(1.05), face = "bold",
                                       color = "#333333"),
      legend.text       = element_text(size = rel(1.0), color = "#333333"),
      legend.key.height = unit(1.0, "cm"),
      legend.key.width  = unit(0.32, "cm"),
      plot.margin       = margin(12, 14, 10, 10),
      plot.background   = element_rect(fill = "white", color = NA)
    )
  if (facet_sites) p <- p + facet_wrap(~ site_display, ncol = 3,
                                       labeller = ggplot2::labeller(
                                         site_display = ggplot2::label_wrap_gen(strip_wrap_width)))
  if (square_cells) p <- p + coord_equal()
  p
}

main <- function() {
  if (!dir.exists(PANEL_DIR)) dir.create(PANEL_DIR, recursive = TRUE)
  cat("Figure 3 panel C renderer (species x continent)\n")
  cat("  script dir       :", script_dir, "\n")
  cat("  dedup BIOM       :", INPUT_BIOM, "\n")
  cat("  dedup metadata   :", metadata_path, "\n")
  cat("  taxonomy         :", INPUT_TAXONOMY, "\n")
  cat("  species table    :", species_file, "\n")
  cat("  PANEL_DIR        :", PANEL_DIR, "\n")
  cat("  detection cutoff : RA >= ", PREV_THRESH, " (",
      PREV_THRESH * 100, "%)\n", sep = "")

  dat <- prepare_data()

  cat("\n-- Panel C: species prevalence by continent --\n")
  sp <- prepare_species_data(dat$meta_clean)
  if (!is.null(sp)) {
    sp$rank_s[, Rank := seq_len(.N)]
    cat("  n =", sp$total_samples,
        "samples; species ranked by global prevalence:\n")
    print(sp$rank_s[, .(Rank, Species, GlobalPrevalence = Prevalence,
                        Detected)], digits = 4)
    cat("  continent denominators:\n")
    print(sp$pd[Species == sp$rank_s$Species[1],
                .(Continent = as.character(Stratum), Sample_N = N)],
          row.names = FALSE)
    p_c <- build_heatmap(sp$pd, CONT_ABBR, rot_x = FALSE, base = 16,
                         phyla = sp$rank_s$Species, facet_sites = FALSE,
                         italic_y = TRUE, label_size = 3.8, axis_mult = 1.0,
                         square_cells = TRUE, wrap_y = TRUE)
    p_c <- p_c + theme(plot.margin = margin(3, 14, 3, 10))
    write_out(p_c, "panelC_species_continent_prevalence", 7.0, 3.7)
    data.table::fwrite(
      sp$pd[, .(Continent = as.character(Stratum),
                Species = as.character(Species), Prevalence, Sample_N = N)],
      file.path(
        PANEL_DIR,
        "Fig3_panelC_species_continent_prevalence_stats_0922.csv"
      )
    )
  }

  cat("\nPanel C written to:", PANEL_DIR, "\n")
}

if (!identical(Sys.getenv("FIG3_NO_MAIN"), "1")) main()
