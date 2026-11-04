# fig04_cioverlap.R
# Figure 4 (fig:sim2-cioverlap-hist): distribution of CI overlap, Models 1-9.
#
# Reviewer 3:
#  * "The x-axis scale overlaps and is confusing. It looks as if it goes from
#    0.000 to 1.000!" -> breaks at 0, 0.5, 1 only, labelled "0", "0.5", "1",
#    with extra panel spacing.
#  * "The colour choice is poor ... overlapping histograms are [not] a good
#    choice here. Coloured step functions might work, faceting might work."
#    -> three alternatives, none of which uses overlapping filled bars:
#       A) coloured step-function (outline) histograms, one per method;
#       B) empirical CDFs (also step functions), which show the point masses
#          of GenM-1/GenM-2 at 0 and 0.5 as vertical jumps;
#       C) faceting: one row per method, one column per model, so no
#          histograms overlap at all.
#    The text on p.27 ("CI overlaps for GenModel are often near 0 or near 0.5";
#    "MV Histogram has a large spread") is visible in each.
#
# Input: sim2_models.csv (model, n_cont, n_bin, method, rep, ci_overlap, ...).

BIN_WIDTH <- 0.1

sim2 <- read_input("sim2_models.csv")
ov <- sim2[sim2$model %in% paste("Model", 1:9), ]
ov$method <- as_method(ov$method)
ov$cont <- factor(paste(ov$n_cont, "continuous"),
                  levels = paste(sort(unique(ov$n_cont)), "continuous"))
ov$bin  <- factor(paste(ov$n_bin, "binary"),
                  levels = paste(sort(unique(ov$n_bin)), "binary"))

x_overlap <- scale_x_continuous("CI overlap", limits = c(0, 1),
                                breaks = c(0, 0.5, 1),
                                labels = c("0", "0.5", "1"),
                                expand = expansion(add = 0.03))
col_method <- scale_colour_manual("Method", values = METHOD_COLORS)

# --- A) step-function histograms --------------------------------------------
# Histogram counts are computed here and drawn with geom_step over the bin
# edges, so each method is a single unfilled outline. Bins are closed on the
# left, [a, b), with the last bin closed on both sides so overlap = 1 counts.
step_hist <- function(d) {
  edges <- seq(0, 1, by = BIN_WIDTH)
  out <- lapply(split(d, list(d$method, d$cont, d$bin), drop = TRUE), function(g) {
    n <- as.vector(table(cut(g$ci_overlap, edges, right = FALSE,
                             include.lowest = TRUE)))
    data.frame(method = g$method[1], cont = g$cont[1], bin = g$bin[1],
               x = edges, count = c(n, n[length(n)]))
  })
  do.call(rbind, out)
}

step_hist_model <- function(d) {
  edges <- seq(0, 1, by = BIN_WIDTH)
  out <- lapply(split(d, list(d$model, d$method), drop = TRUE), function(g) {
    n <- as.vector(table(cut(g$ci_overlap, edges, right = FALSE,
                             include.lowest = TRUE)))
    data.frame(method = g$method[1], model = g$model[1],
               x = edges, count = c(n, n[length(n)]))
  })
  do.call(rbind, out)
}



pA <- ggplot(step_hist(ov), aes(x, count, colour = method)) +
  geom_step(direction = "hv", linewidth = 0.8) +
  facet_grid(bin ~ cont) +
  x_overlap + col_method +
  labs(y = "Number of repetitions") +
  theme_paper()
save_fig(pA, "sim2_cioverlap_step", width = 7, height = 6)

pA2 <- ggplot(step_hist_model(ov), aes(x, count, colour = model)) +
  geom_step(direction = "hv", linewidth = 0.8) +
  facet_wrap( ~ method, scales="free_y") +
  x_overlap + #col_method +
  labs(y = "Number of repetitions") +
  theme_paper()
save_fig(pA2, "sim2_cioverlap_step_model_facet_method", width = 7, height = 6)

# --- B) empirical CDFs ------------------------------------------------------
pB <- ggplot(ov, aes(ci_overlap, colour = method)) +
  stat_ecdf(linewidth = 0.8, pad = TRUE) +
  facet_grid(bin ~ cont) +
  x_overlap + col_method +
  scale_y_continuous("Share of repetitions with overlap ≤ x",
                     breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  theme_paper()
save_fig(pB, "sim2_cioverlap_ecdf", width = 8, height = 6)

# --- C) faceted histograms, one method per row ------------------------------
lab <- sprintf("%s\n%d cont., %d bin.", ov$model, ov$n_cont, ov$n_bin)
ov$model_lab <- factor(lab, levels = unique(lab[order(match(ov$model, MODEL_LEVELS))]))
pC <- ggplot(ov, aes(ci_overlap, fill = method)) +
  geom_histogram(binwidth = BIN_WIDTH, boundary = 0, closed = "left",
                 colour = "white", linewidth = 0.2) +
  facet_grid(method ~ model_lab) +
  scale_x_continuous("CI overlap", limits = c(0, 1), breaks = c(0, 0.5, 1),
                     labels = c("0", "0.5", "1")) +
  scale_fill_manual(values = METHOD_COLORS, guide = "none") +
  labs(y = "Number of repetitions") +
  theme_paper(base_size = 9) +
  theme(strip.text.y = element_text(angle = 0, hjust = 0),
        panel.spacing.x = unit(0.6, "lines"))
save_fig(pC, "sim2_cioverlap_facet", width = 11, height = 4.5)


pC2 <- ggplot(ov, aes(ci_overlap, fill = method)) +
  geom_histogram(binwidth = BIN_WIDTH, boundary = 0, closed = "left",
                 colour = "white", linewidth = 0.2) +
  facet_grid(method ~ model_lab, scales="free_y") +
  scale_x_continuous("CI overlap", limits = c(0, 1), breaks = c(0, 0.5, 1),
                     labels = c("0", "0.5", "1")) +
  scale_fill_manual(values = METHOD_COLORS, guide = "none") +
  labs(y = "Number of repetitions") +
  theme_paper(base_size = 9) +
  theme(strip.text.y = element_text(angle = 0, hjust = 0),
        panel.spacing.x = unit(0.6, "lines"))
save_fig(pC2, "sim2_cioverlap_facet_freey", width = 11, height = 4.5)

