# 01_make_mock_data.R
# Writes MOCK inputs with the same structure as the real ones (see README.md)
# so the figure scripts can be run and reviewed without the simulation and
# Liberia outputs. The numbers are made up; never use these figures in the
# paper.

set.seed(20261018)
mockdir <- file.path(CODEDIR, "data-mock")
dir.create(mockdir, showWarnings = FALSE, recursive = TRUE)

reps <- 1:20

# --- Simulation 1: privacy budget (Fig 3 and its with-outliers version) ----
eps <- c(0.5, 1, 2, 4)
spread <- c("MV Histogram" = 0.18, "Hybrid" = 0.09, "GenModel" = 0.20)
sim1 <- do.call(rbind, lapply(names(spread), function(m)
  do.call(rbind, lapply(eps, function(e) {
    tr <- rnorm(length(reps), 5.2, spread[[m]] / sqrt(e))
    if (m == "GenModel") {  # GenModel is prone to large outliers
      k <- sample(reps, 2)
      tr[k] <- tr[k] + sample(c(-1, 1), 2, TRUE) * runif(2, 1, 4)
    }
    data.frame(method = m, epsilon = e, rep = reps, treff = tr)
  }))))
sim1$conf_treff <- 5.2
sim1$true_treff <- 5.0
write.csv(sim1, file.path(mockdir, "sim1_privbudget.csv"), row.names = FALSE)

# --- Simulation 2: number/type of covariates (Figs 4, 5, 6, 18) -------------
models <- data.frame(
  model  = c(paste("Model", 1:9), "Unbounded", "Finite Set"),
  n_cont = c(2, 5, 10, 2, 5, 10, 2, 5, 10, NA, NA),
  n_bin  = c(2, 2, 2, 5, 5, 5, 10, 10, 10, NA, NA)
)
models$conf_treff <- round(rnorm(nrow(models), 5.1, 0.05), 3)
sd2 <- c("MV Histogram" = 0.2, "Hybrid" = 0.1, "GenM-1" = 0.8, "GenM-2" = 0.35)
sim2 <- do.call(rbind, lapply(seq_len(nrow(models)), function(i)
  do.call(rbind, lapply(names(sd2), function(m) {
    tr <- rnorm(length(reps), models$conf_treff[i], sd2[[m]])
    ov <- switch(m,
      "MV Histogram" = runif(length(reps), 0, 1),
      "Hybrid"       = pmin(1, rbeta(length(reps), 6, 1.5)),
      # GenModel CI overlap piles up at exactly 0 and near 0.5
      sample(c(0, 0.5, 0.6), length(reps), TRUE, c(0.45, 0.45, 0.1)))
    cbind(models[i, ], method = m, row.names = NULL, rep = reps, treff = tr, ci_overlap = ov)
  }))))
write.csv(sim2, file.path(mockdir, "sim2_models.csv"), row.names = FALSE)

# --- Liberia application (Figs 7 and 8) -------------------------------------
lib <- expand.grid(
  covset    = c("Full Baseline Covariates", "Subset of Covariates",
                "Reduced Subset of Covariates"),
  treatment = c("Cash Only", "Therapy Only", "Both Therapy and Cash"),
  response  = c("fam_asb_lt", "stealnb_ltav", "carryweapon_ltav",
                "asbhostilstd_ltav"),
  method    = c("MV Histogram", "Hybrid", "GenM-1", "GenM-2", "GenM-9"),
  stringsAsFactors = FALSE
)
lib$ci_overlap <- round(runif(nrow(lib)), 3)
lib$absdiff    <- signif(10^runif(nrow(lib), -3.5, 1.5), 3)
write.csv(lib, file.path(mockdir, "liberia_utility.csv"), row.names = FALSE)

message("Mock data written to ", mockdir)
