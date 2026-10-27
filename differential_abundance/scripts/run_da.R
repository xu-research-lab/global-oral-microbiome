#!/usr/bin/env Rscript
# Unified entry point for the four retained differential-abundance figures.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || any(args %in% c("--help", "-h"))) {
  cat("Usage: Rscript differential_abundance/scripts/run_da.R <age|continent|site|sex>\n",
    "  --stage=plot|prepare|consensus|fit|finalize|all  (default: plot)\n",
    "  --method=all|linda|maaslin2|ancombc2   (comma-separated; fit only)\n",
    "  --only=INDEX_OR_COMPARISON             (fit only)\n",
    "  --output-dir=PATH                      (final plot output directory)\n",
    "Environment: GLOBAL_ORAL_ROOT, DA_RESULTS_ROOT\n",
    sep = ""
  )
  quit(save = "no", status = 0)
}
analysis <- args[1]
if (!analysis %in% c("age", "continent", "site", "sex")) stop("Unknown analysis: ", analysis)
options <- args[-1]
allowed <- c("--stage=", "--method=", "--only=", "--output-dir=")
for (a in options) if (!any(startsWith(a, allowed))) stop("Unknown option: ", a)
arg_value <- function(prefix, default = "") {
  hit <- options[startsWith(options, prefix)]
  if (length(hit)) substring(tail(hit, 1), nchar(prefix) + 1L) else default
}
stage <- arg_value("--stage=", "plot")
if (!stage %in% c("plot", "prepare", "consensus", "fit", "finalize", "all")) stop("Unknown stage: ", stage)
script_file <- sub("^--file=", "", commandArgs()[startsWith(commandArgs(), "--file=")][1])
project_root <- Sys.getenv("GLOBAL_ORAL_ROOT", "")
if (!nzchar(project_root)) project_root <- file.path(dirname(normalizePath(script_file)), "../..")
project_root <- normalizePath(project_root, mustWork = TRUE)
source(file.path(project_root, "differential_abundance/scripts/config.R"))
env <- da_profile(analysis, project_root)
env$method_arg <- arg_value("--method=", "all")
env$only_arg <- arg_value("--only=")
env$plot_output_dir <- arg_value("--output-dir=", env$figures_dir)
selected_stages <- switch(stage,
  all = c("fit", "consensus", "prepare", "plot"),
  finalize = c("consensus", "prepare", "plot"),
  stage
)
for (s in selected_stages) {
  file <- switch(s,
    fit = "fit.R",
    consensus = "consensus.R",
    prepare = "prepare_plot_data.R",
    plot = "plot.R"
  )
  message(analysis, ": ", s)
  sys.source(file.path(env$script_dir, file), envir = env)
}
