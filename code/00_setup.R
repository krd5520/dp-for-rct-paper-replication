# 00_setup.R
# Shared configuration for the figures revised in response to Reviewer 3
# (HDSR-2023-078R1). Sourced by every fig*.R script via main.R.
#
# Dependencies: ggplot2 (>= 3.4, for `linewidth`) and scales only.

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

# ---------------------------------------------------------------------------
# Paths. DATADIR holds the input CSVs described in README.md; by default the
# mock data written by 01_make_mock_data.R. Point it at the real simulation /
# Liberia outputs by setting the environment variable R3_DATADIR.
# ---------------------------------------------------------------------------
CODEDIR <- if (exists("CODEDIR")) CODEDIR else "code"
DATADIR <- Sys.getenv("R3_DATADIR", file.path(CODEDIR, "data-mock"))
OUTDIR  <- Sys.getenv("R3_OUTDIR",  file.path(CODEDIR, "output"))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

read_input <- function(name) {
  f <- file.path(DATADIR, name)
  if (!file.exists(f)) stop("Missing input file: ", f, " (see code/README.md)")
  read.csv(f, stringsAsFactors = FALSE)
}

save_fig <- function(p, name, width, height, dpi = 300) {
  f <- file.path(OUTDIR, paste0(name, ".png"))
  ggsave(f, p, width = width, height = height, dpi = dpi, bg = "white")
  message("Wrote ", f)
  invisible(f)
}

# ---------------------------------------------------------------------------
# Methods: one fixed order and one fixed colour per method, used in every
# figure. Okabe-Ito palette (colour-blind safe). "Confidential" is a
# reference, never a method, so it is not in this list (R3, Figs 7/8).
# ---------------------------------------------------------------------------
# "GenModel" is the generic label used when epsilon varies along an axis
# (Fig 3); "GenM-<eps>" when the budget is fixed (\dpmb{} in notation.tex).
METHOD_LEVELS <- c("MV Histogram", "Hybrid", "GenModel",
                   "GenM-1", "GenM-2", "GenM-9")
METHOD_COLORS <- c(
  "MV Histogram" = "#009E73",  # bluish green
  "Hybrid"       = "#E69F00",  # orange
  "GenModel"     = "#0072B2",  # blue
  "GenM-1"       = "#0072B2",  # blue
  "GenM-2"       = "#CC79A7",  # reddish purple
  "GenM-9"       = "#56B4E9"   # sky blue
)

# Reference values are drawn in the same colours everywhere.
COL_CONFIDENTIAL <- "black"  # contrasts with every method colour
COL_TRUE         <- "grey60"

MODEL_LEVELS <- c(paste("Model", 1:9), "Unbounded", "Finite Set")

as_method <- function(x) {
  bad <- setdiff(unique(x), METHOD_LEVELS)
  if (length(bad)) stop("Unknown method label(s): ", paste(bad, collapse = ", "))
  factor(x, levels = intersect(METHOD_LEVELS, unique(x)))
}

theme_paper <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text       = element_text(face = "bold"),
      legend.position  = "bottom",
      panel.spacing    = unit(1, "lines")
    )
}
