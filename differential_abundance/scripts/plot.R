sys.source(file.path(script_dir, "plot_common.R"), envir = environment())
if (analysis == "sex") {
  sys.source(file.path(script_dir, "plot_sex.R"), envir = environment())
} else {
  sys.source(file.path(script_dir, "selection.R"), envir = environment())
  read_plot_table <- function(metric) {
    fread(file.path(
      figures_dir,
      paste0("figureB_", metric, "_data_2method.tsv")
    ), data.table = FALSE)
  }
  heat <- read_plot_table("heatmap")
  bubble <- read_plot_table("bubble")
  prev <- read_plot_table("prevalence")
  rows <- heat %>%
    distinct(feature, y_plot) %>%
    arrange(y_plot)
  stopifnot(
    !anyDuplicated(rows$feature), !anyDuplicated(rows$y_plot),
    all(bubble$y_plot == rows$y_plot[match(bubble$feature, rows$feature)]),
    all(prev$y_plot == rows$y_plot[match(prev$feature, rows$feature)]),
    isTRUE(all.equal(prev$prevalence, 100 * prev$n_present / prev$n_total))
  )
  if (analysis == "age") {
    sys.source(file.path(script_dir, "plot_age.R"), envir = environment())
    display_rows <- diff(range(rows$y_plot)) + 1
    save_final(
      render_age(heat, bubble, prev), 7.7 * 3 / 4,
      max(5.3333, 0.13 * display_rows + 1.5)
    )
  } else {
    sys.source(file.path(script_dir, "plot_ova.R"), envir = environment())
    save_final(render_ova(heat, bubble, prev), 120 / 25.4, 205 / 25.4)
  }
}
