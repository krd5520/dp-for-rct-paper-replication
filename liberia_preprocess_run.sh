#!/bin/bash
#
# Bash script to run all the paper analyses
# Alternatively:
# - open the "dp-for-rct.Rproj" in Rstudio, and "knit" the README.Rmd

cd program

R CMD BATCH liberia_replicate_original_table2b.R
R CMD BATCH liberia_subset_process_data.R
R CMB BATCH liberia_variable_tables.R