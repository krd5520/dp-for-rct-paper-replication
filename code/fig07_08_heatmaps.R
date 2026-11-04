# fig07_08_heatmaps.R
# Figure 7 (fig:liberia-ciover-subred, fig:liberia-ciover-full) and
# Figure 8 (fig:liberia-absdiff-subred, fig:liberia-absdiff-full).
#
# Reviewer 3:
#  * "Why is there a row for Confidential in Figures 7 and 8?"
#    -> methods come from as_method(), which has no "Confidential" level, so
#       the empty row is gone.
#  * "In Fig 7 the scale goes from bad to good. In Fig 8 with the same
#    colours, the scale goes from good to bad. It would be better to use
#    different colours for the different measures and apply the same light to
#    dark scale direction."
#    -> CI overlap uses blues, absolute difference uses oranges, and in BOTH
#       darker means closer to the confidential result (better utility).
#  * The three covariate sets are facet rows of one plot, replacing the two
#    stacked images with \phantom padding in section-real-world.tex. Every
#    legend bin is shown, even if no cell falls in it (the current Fig 7
#    legend silently omits 0.125-0.25).
#  * Absolute difference is binned by decades (0.001, 0.01, ...) instead of
#    the unrounded log-spaced cut points (0.007, 0.041, 0.222, ...).
#
# Input: liberia_utility.csv (covset, treatment, response, method,
#        ci_overlap, absdiff).

COVSET_LEVELS <- c("Full Baseline Covariates", "Subset of Covariates",
                   "Reduced Subset of Covariates")
TREAT_LEVELS  <- c("Cash Only", "Therapy Only", "Both Therapy and Cash")
RESP_LEVELS   <- c("fam_asb_lt", "stealnb_ltav", "carryweapon_ltav",
                   "asbhostilstd_ltav")

lib <- read_input("liberia_utility.csv")
lib <- lib[lib$method != "Confidential", ]  # in case the input still has it
lib$method    <- as_method(lib$method)
lib$covset    <- factor(lib$covset,    levels = intersect(COVSET_LEVELS, lib$covset))
lib$treatment <- factor(lib$treatment, levels = TREAT_LEVELS)
lib$response  <- factor(lib$response,  levels = RESP_LEVELS)

# CI overlap: 8 equal bins on [0, 1]; high overlap = dark.
ov_breaks <- seq(0, 1, by = 0.125)
lib$ov_bin <- cut(lib$ci_overlap, ov_breaks, include.lowest = TRUE)

# Absolute difference: decade bins; SMALL difference = dark.
ad_breaks <- c(0, 0.001, 0.01, 0.1, 1, 10, Inf)
ad_labels <- c("≤ 0.001", "0.001–0.01", "0.01–0.1", "0.1–1",
               "1–10", "> 10")
lib$ad_bin <- cut(lib$absdiff, ad_breaks, labels = ad_labels, include.lowest = TRUE)

plot_heat <- function(d, fillvar, palette, direction, legend_title) {
  ggplot(d, aes(response, method, fill = .data[[fillvar]])) +
    geom_tile(colour = "white", linewidth = 0.4) +
    facet_grid(covset ~ treatment) +
    scale_y_discrete(limits = rev) +
    scale_fill_brewer(legend_title, palette = palette, direction = direction,
                      drop = FALSE) +
    labs(x = "Response variable", y = NULL) +
    theme_paper(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid = element_blank(),
          strip.text.y = element_text(angle = 0, hjust = 0),
          legend.position = "right")
}

p7 <- plot_heat(lib, "ov_bin", "Blues", 1,
                "CI overlap\n(darker = closer to\nconfidential)")
p8 <- plot_heat(lib, "ad_bin", "Oranges", -1,
                "Absolute difference\n(darker = closer to\nconfidential)")

save_fig(p7, "cioverlap_heatmap_allcovsets", width = 9, height = 6.5)
save_fig(p8, "absdiff_heatmap_allcovsets",   width = 9, height = 6.5)
