#!/usr/bin/env Rscript
# plot_jaccard_facets_paired.R — paired-run variant of plot_jaccard_facets.R
# (the ORIGINAL script is kept byte-identical; only OUT_DIR differs and
# reads/writes the paired-pipeline results):
#   in : result_paired/jaccard/jaccard_<dis>_vs_HC.csv
#   out: result_paired/jaccard/jaccard_facets.png (+ .pdf)

suppressMessages(library(ggplot2))

ROOT <- normalizePath(file.path(dirname(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE))), ".."))
OUT_DIR <- file.path(ROOT, "result_paired", "jaccard")

a <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
    k <- which(a == flag)
    if (length(k) && length(a) >= k + 1) return(a[k + 1])
    if (length(k)) return("")
    default
}
dis_arg <- get_arg("--dis", "PD,DC,OSCC")
out_arg <- get_arg("--out", "_paired")

DIS <- strsplit(dis_arg, ",", fixed = TRUE)[[1]]
DIS_LAB <- c(PD = "Periodontal disease", DC = "Dental caries",
             OSCC = "OSCC")

long <- do.call(rbind, lapply(DIS, function(dis) {
    m <- read.csv(file.path(OUT_DIR, sprintf("jaccard_%s_vs_HC.csv", dis)),
                  row.names = 1, check.names = FALSE)
    data.frame(
        dis     = factor(dis, levels = DIS),
        dis_mod = factor(rep(sub("^[A-Z]+-", "", colnames(m)),
                             each = nrow(m)),
                         levels = sub("^[A-Z]+-", "", colnames(m))),
        hc      = factor(rep(rownames(m), times = ncol(m)),
                         levels = rev(rownames(m))),
        J       = as.vector(as.matrix(m))
    )
}))

p <- ggplot(long, aes(dis_mod, hc, fill = J)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", J)), size = 2.0,
              color = "grey15") +
    facet_grid(. ~ dis, scales = "free_x", space = "free_x",
               labeller = labeller(dis = DIS_LAB)) +
    scale_fill_gradient(low = "#f7fbff", high = "#08306b",
                        name = "Jaccard", limits = c(0, 1),
                        breaks = seq(0, 1, 0.25),
                        guide = guide_colourbar(barwidth = 0.3,
                                                barheight = 2.1,
                                                frame.colour = "grey40",
                                                ticks.colour = "grey40")) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(
        strip.background = element_rect(fill = "grey85", color = NA),
        strip.text       = element_text(face = "bold", size = 9,
                                        color = "grey15",
                                        margin = margin(2, 0, 2, 0)),
        panel.grid       = element_blank(),
        panel.border     = element_rect(color = "grey40", linewidth = 0.4),
        axis.text.x      = element_text(size = 7.5),
        axis.text.y      = element_text(size = 7.5),
        axis.ticks       = element_line(linewidth = 0.3),
        legend.text      = element_text(size = 7),
        legend.title     = element_text(size = 8),
        legend.key.spacing.y = grid::unit(0.1, "cm"),
        plot.margin      = margin(2, 4, 2, 2)
    )

n_hc <- length(unique(long$hc))
h_png <- max(2.07, 0.30 + 0.17 * n_hc)
w_png <- switch(length(DIS), 2.8, 5.0, 7.0)
base <- if (nzchar(out_arg)) paste0("jaccard_facets", out_arg) else "jaccard_facets"
ggsave(file.path(OUT_DIR, paste0(base, ".png")), p,
       width = w_png, height = h_png, dpi = 300)
ggsave(file.path(OUT_DIR, paste0(base, ".pdf")), p,
       width = w_png, height = h_png, device = cairo_pdf)
cat("done ->", file.path(OUT_DIR, paste0(base, ".png")), "\n")