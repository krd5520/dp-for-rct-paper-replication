# main.R
# Regenerate all figures revised in response to Reviewer 3.
# Run from the repository root:   Rscript code/main.R
# With real inputs:               R3_DATADIR=/path/to/csvs Rscript code/main.R

CODEDIR <- "code"
source(file.path(CODEDIR, "00_setup.R"))

if (DATADIR == file.path(CODEDIR, "data-mock")) {
  message("R3_DATADIR not set: using MOCK data. Do not put these figures in the paper.")
  source(file.path(CODEDIR, "01_make_mock_data.R"))
}

for (f in c("fig03_privbudget.R", "fig04_cioverlap.R",
            "fig05_06_treff_models.R", "fig07_08_heatmaps.R")) {
  message("== ", f)
  source(file.path(CODEDIR, f), local = new.env())
}
