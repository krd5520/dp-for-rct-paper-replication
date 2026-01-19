#!/bin/bash
#
# Bash script to run all the paper analyses
# Alternatively:
# - open the "dp-for-rct.Rproj" in Rstudio, and "knit" the README.Rmd

cd program

R CMD BATCH simulations.R
R CMD BATCH simulations_tables_and_figures.R