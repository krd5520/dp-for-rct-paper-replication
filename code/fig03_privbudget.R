# fig03_privbudget.R
# Figure 3 (fig:sim-privbudget-boxplots) and its appendix companion.
#
# Reviewer 3:
#  * "Why do the x-axis scales in Figure 3 include 5?"
#    -> epsilon is drawn as a discrete axis with exactly the budgets that were
#       run (0.5, 1, 2, 4), so no budget that was not simulated is shown.
#  * "Fig 18 is not a version of Fig 3" (caption cross-reference)
#    -> we now produce the true with-outliers version of Fig 3
#       (sim1_treff_BudgetX_AllMethods_KeepOutliers.png), which is what the
#       Fig 3 caption promises. It replaces the reference to Fig 18.
#
# Input: sim1_privbudget.csv (method, epsilon, rep, treff, conf_treff,
#        true_treff). Epsilon values not wanted in the plot (e.g. 5000) are
#        dropped via EPS_SHOW.

EPS_SHOW <- c(0.5, 1, 2, 4)

sim1 <- read_input("sim1_privbudget.csv")
sim1 <- sim1[sim1$epsilon %in% EPS_SHOW, ]
sim1$method  <- as_method(sim1$method)
sim1$eps_lab <- factor(format(sim1$epsilon, drop0trailing = TRUE),
                       levels = format(sort(unique(sim1$epsilon)),
                                       drop0trailing = TRUE))

refs <- data.frame(
  ref   = c("Confidential estimate", "True treatment effect"),
  value = c(unique(sim1$conf_treff), unique(sim1$true_treff))
)
stopifnot(nrow(refs) == 2)  # one confidential and one true value

plot_privbudget <- function(d, keep_outliers) {
  p <- ggplot(d, aes(x = eps_lab, y = treff)) +
    geom_hline(data = refs,
               aes(yintercept = value, colour = ref, linetype = ref),
               linewidth = 0.7) +
    geom_boxplot(width = 0.55, outlier.shape = if (keep_outliers) 19 else NA,
                 outlier.size = 1.2) +
    facet_wrap(~ method, nrow = 1) +
    scale_colour_manual(NULL, values = c("Confidential estimate" = COL_CONFIDENTIAL,
                                         "True treatment effect" = COL_TRUE)) +
    scale_linetype_manual(NULL, values = c("Confidential estimate" = "dashed",
                                           "True treatment effect" = "solid")) +
    labs(x = expression("Privacy budget (" * epsilon * ")"),
         y = "Sanitized treatment effect") +
    theme_paper()
  if (!keep_outliers) {
    # Zoom to the whiskers rather than deleting points, so the boxes still
    # summarise all 20 repetitions.
    w <- unlist(tapply(d$treff, list(d$method, d$eps_lab),
                       function(v) range(boxplot.stats(v)$stats)))
    p <- p + coord_cartesian(ylim = range(c(w, refs$value)))
  }
  p
}

save_fig(plot_privbudget(sim1, keep_outliers = FALSE),
         "sim1_treff_BudgetX_AllMethods_RemoveOutliers", width = 7, height = 3.6)
save_fig(plot_privbudget(sim1, keep_outliers = TRUE),
         "sim1_treff_BudgetX_AllMethods_KeepOutliers", width = 7, height = 3.6)
