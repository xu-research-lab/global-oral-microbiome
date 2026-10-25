suppressPackageStartupMessages({
  library(png)
  library(grid)
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
BASE_DIR     <- PROJECT_ROOT

PANEL_DIR <- file.path(script_dir, "results", "Fig1")
BASE_DIR    <- Sys.getenv("FIG1_BASE_DIR",  unset = BASE_DIR)
PANEL_DIR   <- Sys.getenv("FIG1_PANEL_DIR", unset = PANEL_DIR)
AUTO_INSTALL <- as.logical(Sys.getenv("FIG1_AUTO_INSTALL", unset = "FALSE"))


FIG1_TEXT_SCALE  <- as.numeric(Sys.getenv("FIG1_TEXT_SCALE",  unset = "1.18"))
FIG1_PANEL_SCALE <- as.numeric(Sys.getenv("FIG1_PANEL_SCALE", unset = "0.88"))
fig1_base_sz <- function(env_var, default) {
  configured <- Sys.getenv(env_var, unset = "")
  base <- if (nzchar(configured)) as.numeric(configured) else default
  base * FIG1_TEXT_SCALE
}
fig1_panel_in <- function(x) x * FIG1_PANEL_SCALE

fig1_site_levels <- c(
  "Saliva", "Tooth", "Oral Mucosa", "Oral Swab",
  "Tongue", "Gingiva", "Others"
)
fig1_site_colors <- c(
  "Saliva"      = "#0072B2",
  "Tooth"       = "#D55E00",
  "Oral Mucosa" = "#009E73",
  "Oral Swab"   = "#CC79A7",
  "Tongue"      = "#E69F00",
  "Gingiva"     = "#6A3D9A",
  "Others"      = "#56B4E9"
)

metadata_path <- Sys.getenv(
  "FIG1_METADATA",
  unset = file.path(DATA_DIR, "metadata_oral_qc.tsv")
)

panel_png <- c(
  C      = "global_sampling_world_map_final_cropped_260928.png"
)

save_plot_png_svg <- function(plot, png_path, width, height, bg = "white") {
  ggplot2::ggsave(png_path, plot = plot, width = width, height = height,
                  device = "png", dpi = 600, bg = bg)
  invisible(png_path)
}

load_pkgs <- function(pkgs, auto_install = FALSE) {
  pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(pkgs) == 0) return(invisible(NULL))
  if (auto_install) {
    install.packages(pkgs, repos = "http://cran.us.r-project.org")
    pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(pkgs) == 0) return(invisible(NULL))
  }
  paste("missing package(s):", paste(pkgs, collapse = ", "))
}
load_libs <- function(pkgs, auto_install = FALSE) {
  ok <- load_pkgs(pkgs, auto_install)
  if (!is.null(ok)) return(ok)
  suppressPackageStartupMessages(lapply(pkgs, library, character.only = TRUE))
  invisible(NULL)
}

pick_meta_col <- function(md, preferred, fallback = NULL) {
  if (preferred %in% names(md)) return(preferred)
  if (!is.null(fallback) && fallback %in% names(md)) return(fallback)
  preferred
}

render_panel_C_world_map <- function(out_dir) {
  if (!is.null(err <- load_libs(c("tidyverse", "sf", "countrycode", "rnaturalearth",
                                  "rnaturalearthdata", "ggrepel", "scales",
                                  "patchwork", "ggtext"), AUTO_INSTALL))) return(err)
  if (!file.exists(metadata_path))
    return(sprintf("metadata table not found: %s", metadata_path))

  sf_use_s2(FALSE)

  base_sz  <- fig1_base_sz("FIG1_C_BASE_SZ", 28)
  text_sz  <- base_sz / .pt

  theme_map <- theme_void(base_size = base_sz, base_family = "sans") +
    theme(
      legend.title   = element_text(size = base_sz * 0.9, face = "bold", color = "#303030"),
      legend.text    = element_text(size = base_sz * 0.8, color = "#303030"),
      legend.position = "top",
      legend.box      = "horizontal",
      legend.spacing.x = unit(3, "mm"),
      legend.key      = element_blank(),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin      = margin(5, 5, 5, 5, "pt")
    )

  continent_order <- c("NAM", "SA", "EU", "AS", "AF", "OC")

  md_raw <- read_tsv(metadata_path, show_col_types = FALSE)
  col_site    <- pick_meta_col(md_raw, "site_cleaned",    "site")
  col_country <- pick_meta_col(md_raw, "country_cleaned", "country")
  col_age     <- pick_meta_col(md_raw, "age_cleaned",     "age")
  col_bmi     <- pick_meta_col(md_raw, "bmi_cleaned",     "bmi")
  col_sex     <- pick_meta_col(md_raw, "sex_cleaned",     "sex")

  metadata <- tibble::tibble(
    BioProject      = md_raw[["bioproject_id"]],
    subject_id      = md_raw[[pick_meta_col(md_raw, "subject_id_cleaned", "subject_id")]],
    Country_Cleaned = md_raw[[col_country]],
    Site_Cleaned    = md_raw[[col_site]]
  ) %>%
    mutate(
      Site_Cleaned = str_to_title(Site_Cleaned),
      Site_Cleaned = case_when(Site_Cleaned == "Other" ~ "Others", TRUE ~ Site_Cleaned)
    ) %>%
    filter(!(Site_Cleaned %in% c("Nasal", "Tonsils")), !is.na(Country_Cleaned),
           Country_Cleaned != "") %>%
    mutate(iso3c = countrycode(Country_Cleaned, origin = "country.name", destination = "iso3c"))

  south_america_iso3c <- c("ARG", "BOL", "BRA", "CHL", "COL", "ECU", "GUY", "PRY",
                           "PER", "SUR", "URY", "VEN", "FLK", "GUF")

  meta_clean <- metadata %>%
    mutate(
      Temp_Continent = countrycode(iso3c, "iso3c", "continent"),
      Continent = case_when(
        Temp_Continent == "Americas" & (iso3c %in% south_america_iso3c) ~ "SA",
        Temp_Continent == "Americas" & !(iso3c %in% south_america_iso3c) ~ "NAM",
        Temp_Continent == "Europe"   ~ "EU",
        Temp_Continent == "Asia"     ~ "AS",
        Temp_Continent == "Africa"   ~ "AF",
        Temp_Continent == "Oceania"  ~ "OC",
        iso3c == "ATA"               ~ "Antarctica",
        TRUE                         ~ Temp_Continent
      ),
      age_num = suppressWarnings(as.numeric(md_raw[[col_age]])),
      bmi_num = suppressWarnings(as.numeric(md_raw[[col_bmi]])),
      sex_label = case_when(
        tolower(md_raw[[col_sex]]) %in% c("f", "female") ~ "F",
        tolower(md_raw[[col_sex]]) %in% c("m", "male")   ~ "M",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!Continent %in% c("Other", "Antarctica")) %>%
    mutate(Continent = factor(Continent, levels = continent_order))

  # Count subjects once per (BioProject, subject_id) across the whole compendium,
  # assigning each subject to the continent of its first sample, so the six
  # continent cards sum to the global unique-subject total.
  subject_counts <- meta_clean %>%
    distinct(BioProject, subject_id, .keep_all = TRUE) %>%
    count(Continent, name = "n_subj")

  continent_summary <- meta_clean %>%
    group_by(Continent) %>%
    summarise(
      n_countries = n_distinct(Country_Cleaned),
      n_bp   = n_distinct(BioProject),
      n_samp = n(),
      .groups = "drop"
    ) %>%
    left_join(subject_counts, by = "Continent") %>%
    mutate(n_subj = coalesce(n_subj, 0L)) %>%
    mutate(
      Card_Title = case_when(
        Continent == "NAM" ~ "NA",
        Continent == "SA"  ~ "SA",
        Continent == "EU"  ~ "EU",
        Continent == "AS"  ~ "AS",
        Continent == "AF"  ~ "AF",
        Continent == "OC"  ~ "OC"
      ),
      card_label = paste0(
        "<span style='font-size: ", round(base_sz * 0.55), "pt; color:#1F2933;'><b>", Card_Title, "</b></span><br>",
        "<span style='font-size: ", round(base_sz * 0.35), "pt; color:#475569;'>",
        "Countries = ", n_countries, "<br>",
        "Bioprojects = ", n_bp, "<br>",
        "Subjects = ", scales::comma(n_subj), "<br>",
        "Samples = ", scales::comma(n_samp),
        "</span>"
      )
    )

  continent_coords <- data.frame(
    Continent = c("NAM", "SA", "EU", "AS", "AF", "OC"),
    lon_pie   = c(-163, -131,  -66,  132,   -40,   54),
    lat_pie   = c(  15,  -15,   30,   19,   -35,  -32)
  )
  continent_summary <- continent_summary %>%
    left_join(continent_coords, by = "Continent") %>%
    mutate(lon_lab = lon_pie + 9.5)

  cat("Total unique (BioProject, subject_id) subjects:",
      sum(continent_summary$n_subj), "\n")

  country_summary <- meta_clean %>%
    group_by(Country_Cleaned, Continent) %>%
    summarise(
      n_samples = n(),
      n_studies = n_distinct(BioProject),
      .groups = "drop"
    ) %>%
    mutate(iso3c = countrycode(Country_Cleaned, "country.name", "iso3c"))

  world_sf <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf") %>%
    filter(name != "Antarctica") %>%
    mutate(iso3c = iso_a3) %>%
    mutate(
      Temp_Continent = countrycode(iso3c, "iso3c", "continent"),
      Continent = case_when(
        Temp_Continent == "Americas" & (iso3c %in% south_america_iso3c) ~ "SA",
        Temp_Continent == "Americas" & !(iso3c %in% south_america_iso3c) ~ "NAM",
        Temp_Continent == "Europe"   ~ "EU",
        Temp_Continent == "Asia"     ~ "AS",
        Temp_Continent == "Africa"   ~ "AF",
        Temp_Continent == "Oceania"  ~ "OC",
        name == "Somaliland"         ~ "AF",
        name == "N. Cyprus"          ~ "AS",
        name == "Kosovo"             ~ "EU",
        TRUE                         ~ NA_character_
      )
    ) %>%
    select(-Temp_Continent)

  set.seed(42)
  usa_lower48 <- sf::st_as_sfc(
    sf::st_bbox(c(xmin = -125, ymin = 24.4, xmax = -66.9, ymax = 49.4),
                crs = sf::st_crs(world_sf))
  )
  bp_points <- lapply(seq_len(nrow(country_summary)), function(i) {
    cs <- country_summary[i, ]
    idx <- which(world_sf$iso3c == cs$iso3c)
    if (length(idx) == 0) return(NULL)
    geom_i <- world_sf$geometry[idx]
    if (identical(cs$iso3c, "USA")) {
      geom_i <- sf::st_intersection(geom_i, usa_lower48)
      saved_seed <- get(".Random.seed", envir = globalenv())
      set.seed(1)
      pts <- tryCatch(st_sample(geom_i, size = cs$n_studies), error = function(e) NULL)
      assign(".Random.seed", saved_seed, envir = globalenv())
    } else {
      pts <- tryCatch(st_sample(geom_i, size = cs$n_studies), error = function(e) NULL)
    }

    if (is.null(pts) || length(pts) == 0) {
      if (cs$iso3c == "BRB") {
        return(data.frame(
          iso3c = "BRB", lon = rep(-59.5, cs$n_studies),
          lat = rep(13.2, cs$n_studies), n_studies = cs$n_studies
        ))
      }
      return(NULL)
    }
    coords <- st_coordinates(pts)
    data.frame(
      iso3c = cs$iso3c,
      lon   = coords[, 1],
      lat   = coords[, 2],
      n_studies = cs$n_studies
    )
  }) %>%
    bind_rows()

  site_continent_data <- meta_clean %>%
    filter(!is.na(Site_Cleaned)) %>%
    count(Continent, Site_Cleaned) %>%
    group_by(Continent) %>%
    arrange(Site_Cleaned, .by_group = TRUE) %>%
    mutate(
      prop = n / sum(n),
      ymax = cumsum(prop),
      ymin = cumsum(prop) - prop
    ) %>%
    ungroup() %>%
    mutate(Site_Cleaned = factor(Site_Cleaned, levels = fig1_site_levels))

  pie_grobs <- lapply(continent_order, function(cont) {
    df <- site_continent_data %>% filter(Continent == cont)
    if (nrow(df) == 0) return(grid::nullGrob())
    p <- ggplot(df, aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax,
                        fill = Site_Cleaned)) +
      geom_rect(color = "white", linewidth = 0.3) +
      coord_polar("y", start = 0) +
      xlim(0, 1) +
      scale_fill_manual(values = fig1_site_colors, limits = fig1_site_levels) +
      theme_void() +
      theme(legend.position = "none", plot.background = element_blank(),
            panel.background = element_blank())
    ggplotGrob(p)
  })

  r_pie <- 7.5
  pie_layers <- lapply(seq_along(continent_order), function(i) {
    cont <- continent_order[i]
    coords <- continent_coords %>% filter(Continent == cont)
    annotation_custom(
      grob = pie_grobs[[i]],
      xmin = coords$lon_pie - r_pie, xmax = coords$lon_pie + r_pie,
      ymin = coords$lat_pie - r_pie, ymax = coords$lat_pie + r_pie
    )
  })

  p_map_base <- ggplot() +
    geom_sf(data = world_sf, fill = "#F1F3F5", color = "#CBD5E1", linewidth = 0.3) +
    geom_point(data = bp_points, aes(x = lon, y = lat, size = n_studies),
               shape = 21, fill = "#0F766E", color = "white", stroke = 0.3,
               alpha = 0.65) +
    scale_size_continuous(
      trans  = "log10",
      name   = "Number of samples",
      range  = c(1.3, 6.8),
      limits = c(1, 150),
      breaks = c(1, 10, 100),
      labels = c("1", "10", "100"),
      guide  = guide_legend(override.aes = list(color = "white", fill = "#0F766E",
                                                shape = 21, stroke = 0.45,
                                                alpha = 0.95))
    ) +
    pie_layers +
    geom_richtext(
      data = continent_summary,
      aes(x = lon_lab, y = lat_pie, label = card_label, hjust = 0),
      vjust = 0.5,
      fill = "white", color = "#475569",
      label.color = "#E2E8F0", label.size = 0.2,
      label.padding = unit(0.4, "lines"),
      label.r = unit(2, "mm"),
      lineheight = 1.3,
      alpha = 0.97
    ) +
    coord_sf(xlim = c(-180, 180), ylim = c(-55, 75), expand = FALSE, clip = "off") +
    theme_map

  p_site_legend <- ggplot(data.frame(Site = factor(fig1_site_levels,
                                                   levels = fig1_site_levels)),
                          aes(x = 1, y = 1, fill = Site)) +
    geom_point(shape = 21, size = 0, color = NA, stroke = 0) +
    scale_fill_manual(
      name = "Sampling Site",
      values = fig1_site_colors,
      limits = fig1_site_levels,
      guide = guide_legend(override.aes = list(size = 5, color = "white",
                                               stroke = 0.5, alpha = 1))
    ) +
    theme_void() +
    theme(
      legend.position = "left",
      legend.title = element_text(size = base_sz * 0.55, face = "bold", color = "#303030"),
      legend.text = element_text(size = base_sz * 0.55, color = "#303030"),
      legend.key.size = unit(5, "mm"),
      legend.background = element_blank(),
      plot.background = element_blank(),
      panel.background = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )

  p_final_map <- p_map_base +
    inset_element(
      p_site_legend,
      left = 0.01, bottom = 0.05, right = 0.20, top = 0.35,
      align_to = "full", clip = FALSE
    ) +
    plot_annotation(
      theme = theme(
        plot.background  = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = NA, color = NA)
      )
    )

  map_out <- file.path(out_dir, unname(panel_png["C"]))
  if (!grepl("\\.png$", map_out, ignore.case = TRUE)) {
    map_out <- paste0(map_out, ".png")
  }
  grDevices::png(map_out, width = fig1_panel_in(18), height = fig1_panel_in(7.2), units = "in",
                 res = 600, bg = "white")
  print(p_final_map)
  grDevices::dev.off()
  if (!file.exists(map_out)) {
    stop("world map output was not written: ", map_out)
  }
  invisible(NULL)
}

cat("Script dir :", script_dir, "\n")
cat("BASE_DIR   :", BASE_DIR, "\n")
cat("PANEL_DIR  :", PANEL_DIR, "\n")
if (!dir.exists(PANEL_DIR)) dir.create(PANEL_DIR, recursive = TRUE)

cat("\n== rendering panel d - global sampling map ... ")
ans <- tryCatch(suppressMessages(suppressWarnings(render_panel_C_world_map(PANEL_DIR))),
                error = function(e) paste("error:", conditionMessage(e)))
cat(if (is.null(ans)) "OK\n" else paste("SKIPPED:", ans, "\n", sep = " "))
cat("\nDone. Panel written to:", PANEL_DIR, "\n")
