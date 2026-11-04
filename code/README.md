# Figures revised for Reviewer 3 (HDSR-2023-078R1)

R code that redraws the figures criticised by Reviewer 3 (the graphics
reviewer) in the 2026-08-08 decision letter. It needs only `ggplot2` (>= 3.4)
and `scales`.

```bash
Rscript code/main.R                              # mock data -> code/output/
R3_DATADIR=/path/to/csvs Rscript code/main.R     # real data -> code/output/
```

The simulation and Liberia outputs are **not** in this repository. Without
`R3_DATADIR`, `main.R` runs `01_make_mock_data.R` and draws the figures from
made-up numbers, so the design can be reviewed. Do not use those figures in
the paper. PNGs are written to `code/output/` (override with `R3_OUTDIR`) and
keep the file names used by the paper's `\includegraphics`, so once they are
generated from real data they can be copied over `figures-main/` and
`tables-and-figures-appendix/`.

## Input files (`R3_DATADIR`)

All files are long-format CSV, one row per estimate. Method labels must be
one of `MV Histogram`, `Hybrid`, `GenModel`, `GenM-1`, `GenM-2`, `GenM-9`
(see `METHOD_LEVELS` in `00_setup.R`). Rows labelled `Confidential` are
dropped from the heatmaps.

| File | Columns | Used by |
|---|---|---|
| `sim1_privbudget.csv` | `method`, `epsilon`, `rep`, `treff`, `conf_treff`, `true_treff` | Fig 3 + appendix |
| `sim2_models.csv` | `model` (`Model 1`…`Model 9`, `Unbounded`, `Finite Set`), `n_cont`, `n_bin`, `method`, `rep`, `treff`, `conf_treff`, `ci_overlap` | Figs 4, 5, 6, 18 |
| `liberia_utility.csv` | `covset`, `treatment`, `response`, `method`, `ci_overlap`, `absdiff` | Figs 7, 8 |

## What changed, comment by comment

| Reviewer 3 comment | Script | Change |
|---|---|---|
| Fig 3: the caption calls Fig 18 a version of Fig 3, but it isn't | `fig03_privbudget.R` | Also writes `sim1_treff_BudgetX_AllMethods_KeepOutliers.png`, the actual with-outliers version of Fig 3. **The caption needs a text edit (below).** |
| Fig 3: why does the x axis include 5? | `fig03_privbudget.R` | ε is a discrete axis showing only the budgets that were run (0.5, 1, 2, 4). |
| Fig 4: overlapping x-axis labels (reads as "0.000 to 1.000") | `fig04_cioverlap.R` | Ticks only at 0, 0.5 and 1, with wider panel spacing. |
| Fig 4: poor colours, overlapping histograms; suggests step functions or faceting | `fig04_cioverlap.R` | Three alternatives: `sim2_cioverlap_step` (coloured step outlines), `sim2_cioverlap_ecdf` (ECDFs), `sim2_cioverlap_facet` (method × model grid). Colour-blind-safe Okabe–Ito palette; each method has the same colour in every figure. |
| Figs 5/6: triangles for point estimates are distracting | `fig05_06_treff_models.R` | The confidential estimate is a thick black bar, wider than the box. `R3_CONF_MARK=point` gives a filled circle instead. Boxes are no longer coloured by model (redundant with the x axis), and the model labels no longer overlap. |
| Fig 18 is Fig 6 plus Fig 5's third panel | `fig05_06_treff_models.R` | Drawn by the same function as Figs 5 and 6. Either recaption it (below) or drop it. |
| Figs 7/8: why is there a Confidential row? | `fig07_08_heatmaps.R` | Removed. |
| Figs 7/8: the scales run in opposite directions with the same colours | `fig07_08_heatmaps.R` | Blues for CI overlap and oranges for absolute difference. In **both**, darker means closer to the confidential result, and the legend titles say so. Absolute difference is binned by decades. All three covariate sets are one faceted figure, which replaces the two stacked images. |

## Text edits still needed (not R)

1. **Fig 3 caption** (`section-evaluation.tex`, `fig:sim-privbudget-boxplots`):
   point "plots including outliers" at a new appendix figure that includes
   `sim1_treff_BudgetX_AllMethods_KeepOutliers.png`, not at
   `fig:apdx-sim-privbudget`.
2. **Fig 18 caption** (`section-appendix-plots-and-tables.tex`): if the figure
   is kept, describe it as GenM-1 and GenM-2 (Figs 5 and 6) on a common scale.
3. **Fig 9 sentence** (`section-real-world.tex`, paragraph before
   `fig:liberia-some-diag-plt`). The reviewer flags this sentence as unclear:
   "When a model is used to generate synthetic data and is fitted to the
   resulting synthetic data, the underlying assumptions of regression model
   will always meet by the synthetic data ones whether or not they were in the
   original data." A suggested rewrite: "Synthetic data generated from a
   regression model satisfy that model's assumptions by construction. Fitting
   the same model to such data therefore yields well-behaved diagnostics even
   when the assumptions are violated in the original data."
4. **Figs 7/8 LaTeX**: replace the two `\includegraphics` per figure (and the
   `\phantom` padding) with the single `*_heatmap_allcovsets.png`.
