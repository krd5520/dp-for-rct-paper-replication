Replication packages for 

Webb, K.R., Mukherjee, S., Mustafi, A.,  Slavković, A., and Vilhuber, L. (2026), "Assessing Utility of Differential Privacy for RCTs" https://arxiv.org/html/2309.14581v2

You can replicate the analysis done but used two .sh files. However, these are specific to Windows and may require slightly different commands for other OS. Alternate instructions on how to reproduce the analysis with the R code directly will be provided in detail below.

You will need the following R libraries:
  * rporjroot
  * haven
  * dplyr
  * labelled
  * sandwich
  * knitr
  * kableExtra

To replicate the simulations from the paper in Section 4:Evaluating the Mechanisms, run the `run_simulations.sh file`. Alternatively, you can run `simulation.R` and then `simulations_tables_figures.R`.

To replicate the analysis that is performed on the Blattman et al. paper, download the replication package from: 

Blattman, Christopher, Jamison, Julian C., and Sheridan, Margaret. Replication data for: Reducing Crime and Violence: Experimental Evidence from Cognitive Behavioral Therapy in Liberia. Nashville, TN: American Economic Association [publisher], 2025. Ann Arbor, MI: Inter-university Consortium for Political and Social Research [distributor], 2025-12-21. https://doi.org/10.3886/E113056V2

You will need to place the `STYL_Final.dta` file from Blattman et al. in a folder `data/liberia_replication_data`. You should create this folder in the main folder 'dp-for-rct-paper-replication`. Then run `liberia_preprocess_run.sh`. 
Alternatively, in `programs` folder, `liberia_subset_process_data.R` will take the `STYL.dta` data file in, clean the data, preserve only the variables we need, and save the data in it a specified output location.
You can specify the location of the `STYL_Final.dta` file and where the .RDS output files are saved by changing the following found in lines 5-7 of the R file.

outputpath=paste0(basepath,"/output")
input.datapath=paste0(basepath,"/data/liberia_replication_data/STYL_Final.dta")
out.datapath=paste0(basepath,"/data")

To replicate the original analysis by Blattman et al. to get their Table 2b using R, run `liberia_replicate_original_table2b.R` from the programs folder. 

To make the tables and figures from Webb et al. run `


