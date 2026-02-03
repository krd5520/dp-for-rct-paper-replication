---
title: README for "Assessing Utility of Differential Privacy for RCTs' (Kaitlyn R. Webb, Soumya Mukherjee, Aratrika Mustafi, Aleksandra Slavković, Lars Vilhuber)
contributors:
  - Kaitlyn Webb
version: 1.2.rc1
---


## Overview

The code in this replication package constructs the analysis of proposed privacy-preserving mechanisms to create synthetic data for randomized control trial (RCT) replication packages. It includes analysis on simulated data and a utility comparison by replicating the analysis of Blattman et. al (2017).

To get started with the repository run the following commands to clone the project locally:

```bash
git clone https://github.com/krd5520/dp-for-rct-paper-replication.git
```


## Directories Overview
Here is a brief visual overview of the project repository:
```
dp-for-rct-paper-replication/
├── program/ # Contains main files to run simulations and replicate Blattman et al.
│   ├── functions/ #directory contains .R files with the functions called by other .R files in programs/ directory
│   ├── liberia_compare_methods_covariate_sets.R #compares utility of our proposed mechanisms on the Blattman et al. (2025) replication data for three sets of covariates
│   ├── liberia_replicate_original_table2b.R #replicates Blattman et al. (2017) table 2b to verify the analysis process
│   ├── liberia_subset_process_data.R #reads in Blattman et al. (2025) data and cleans it
│   ├── liberia_variable_tables.R # produces latex (.tex) tables of the variables used in analysis by Blattman et al. (2016; 2025)
│   ├── simulations.R #runs simulations
│   └── simulation_table_and_figures.R #creates various plots and tables (as .tex Latex files)
├── DPrct.Rproj			# R project to contain and run the various code within	
├── global-libaries.R # file to install all required libaries
├── liberia_preprocess_run.sh # Shell script which can clean the Blattman et al. (2025) data, replicate the analysis done by Blattman et al. (2016), and produce the variable tables.
├── LICENSE.txt # license for the code and project
├── run_simulations.sh # Shell scrpt to run the simulations and make the figures and tables
└── README.md
```

To run the code you will need to create `data/` and `output/` directories and a sub-directory `data/liberia_replication_data/`.


## Data Availability and Provenance Statements


### Statement about Rights

- [X] I certify that the author(s) of the manuscript have legitimate access to and permission to use the data used in this manuscript. 
- [X] I certify that the author(s) of the manuscript have documented permission to redistribute/publish the data contained within this replication package. Appropriate permission are documented in the LICENSE.txt file file.


### License for Data

Simulation data are licensed under Creative Commons/CC-BY-NC license. See LICENSE.txt for details.


### Summary of Availability

- [X] All data **are** publicly available.
- [ ] Some data **cannot be made** publicly available.
- [ ] **No data can be made** publicly available.



#### Summary of Data Availability


| Data.Name  | Data.Files | Location | Provided | Citation |
| -- | -- | -- | -- | -- | 
| “Replication data for: Reducing Crime and Violence” | User download | data/liberia_replication_data | FALSE | Blattman et al. (2025) |


### Details on Replication data for: Reducing Crime and Violence

The paper uses data from "Replication data for: Reducing Crime and Violence: Experimental Evidence from Cognitive Behavioral Therapy in Liberia" (Blattman et al., 2025). Data is subject to a redistribution restriction, but can be freely downloaded from https://doi.org/10.3886/E113056V2. Select 'Download this project' at the bottom of the page and access an ICPSR account by creating one or through your institution.
We have not redistributed the files, but used them in accordance with its Modified BSD License and Creative Commons Attribution 4.0 International Public License. 

Within th zipped download folder there is:
- `20150503_data/`
  - `dtafiles/STYL_Final.dta`: Save this file as `STYL_FINAL` in directory `data/liberia_replication_data`. 
  - `dofiles/AER_1.do`: Includes analysis details that are references as comments in `program/liberia_replicate_original_table2b.R` to justify variable selection. The file can be viewed through Stata or Notepad text reader.
  - `dofiles/AER_2.do`: Includes analysis details that are references as comments in `program/liberia_replicate_original_table2b.R` to justify variable selection. The file can be viewed through Stata or Notepad text reader
- `LICENSE.txt`: Details the packages redistribution restrictions and its Modified BSD License and Creative Commons Attribution 4.0 International Public License.


We accessed and downloaded the data January 2026 (corresponding to the updates Blattman et al. (2025) made to the package in December 2025). 

Datafile: `data/liberia_replication_data/STYL_Final.dta` (not provided)


## Computational requirements

The required R packages can be installed by running `global-libraries.R`. 

There are several packages that will be installed from the CRAN. Additionally, a package from GitHub (repository `krd5520/DPrct`) is used for the proposed mechanisms. It will be installed using the following R code: 

``` R
remotes::install_github("krd5520/DPrct@2479df9fafbc622b7c270e9910f688df0ca77a6f") #https://github.com/krd5520/DPrct

```


### Software Requirements


- [X] The replication package contains one or more programs to install all dependencies and set up the necessary directory structure. 

- R (code was last run with version 4.2.3)
  - `dplyr`, version 1.1.4
  - `devtools`, version 2.4.6
  - `rprojroot`,  version 2.1.1
  - `ggplot2`, version 4.0.1
  - `knitr`, version 1.51
  - `kableExtra`, version 1.4.0
  - `sandwich`, version 3.1-1
  - `markdown`, version 2.0
  - `haven`, version 2.5.5
  - `labelled`, version 2.16.0
  - `remotes`, version 2.5.0
  - the program `global-libraries.R`  will install all dependencies locally, and should be run once.

Optional portions of the code use bash scripting, which may require Linux.


### Controlled Randomness

> INSTRUCTIONS: Some estimation code uses random numbers, almost always provided by pseudorandom number generators (PRNGs). For reproducibility purposes, these should be provided with a deterministic seed, so that the sequence of numbers provided is the same for the original author and any replicators. While this is not always possible, it is a requirement by many journals' policies. The seed should be set once, and not use a time-stamp. If using parallel processing, special care needs to be taken. If using multiple programs in sequence, care must be taken on how to call these programs, ideally from a main program, so that the sequence is not altered. 
> INSTRUCTIONS: If no PRNG is used, check the appropriate box.
> INSTRUCTIONS: If despite attempts to control for randomness, the results are not fully reproducible, please provide a detailed explanation of why, and ideally what kind of instability in numbers a replicator should expect. 

- [X] Random seed is set within each .R file in `programs/` directory
- [ ] The analysis relies on random number generation, but setting a seed is not possible
- [ ] No Pseudo random generator is used in the analysis described here.

### Memory, Runtime, Storage Requirements


#### Summary time to reproduce

Approximate time needed to reproduce the analyses on a standard (CURRENT YEAR) desktop machine:

- [ ] <10 minutes
- [ ] 10-60 minutes
- [ ] 1-2 hours
- [ ] 2-8 hours
- [X] 8-24 hours
- [ ] 1-3 days
- [ ] 3-14 days
- [ ] > 14 days

#### Summary of required storage space

Approximate storage space needed:

- [ ] < 25 MBytes
- [ ] 25 MB - 250 MB
- [ ] 250 MB - 2 GB
- [X] 2 GB - 25 GB
- [ ] 25 GB - 250 GB
- [ ] > 250 GB

- [ ] Not feasible to run on a desktop machine, as described below.

#### Computational Details

The code was last run on a **4-core Intel-based laptop with Microsoft Windows 11 Home version 10.0.26100Build26100 with 84GB of free space**. 

Detail on computation time can be found in Webb et al. (2026)

## Description of programs/code

- Programs in `programs/liberia_compare_methods_covariate_sets.R` compare the utility of our proposed mechanisms on the Blattman et al. (2025) replication data for three sets of covariates.
- Programs in `programs/liberia_replicate_original_table2b.R` replicate Blattman et al. (2017) table 2b to verify the analysis process. Script in `liberia_preprocess_run.sh` will run various cleaning code and table generation code on Blattman et al. (2025) data.
- Programs in `programs/liberia_subset_process_data.R` read in Blattman et al. (2025) data and cleans it. Script in `liberia_preprocess_run.sh` will run various cleaning code and table generation code on Blattman et al. (2025) data.
- Programs in `programs/liberia_variable_tables.R` produce latex (.tex) tables of the variables used in analysis by Blattman et al. (2016; 2025). Script in `liberia_preprocess_run.sh` will run various cleaning code and table generation code on Blattman et al. (2025) data.
- Programs in `programs/simulations.R` run simulations comparing utility across privacy budgets and utility across covariate models. Script in `run_simulations.sh` will run all simulation code. 
- Programs in `programs/simulation_table_and_figures.R` create various plots and tables (as .tex Latex files) from analysis on the simulation data. Script in `run_simulations.sh` will run all simulation code.

### License for Code

The code is licensed under a Modified BSD License. See LICENSE.txt file for details.

## Instructions to Replicators

- Open `DPrct.Rproj` in Rstudio or similar program
- Edit variables in `user_defined_variables.R` to define paths and outputs
- Run `global-libararies.R` once on a new system to set up the working environment. 
- Download the data files referenced above. Each should be stored in directory `data/liberia_replication_data/`, in the format that you download them in. Extract `STLY_Final.dta` from zipped download and place in `data/liberia_replication_data/`
- Run `run_simulations.R` or:
  - Run `programs/simulations.R`
  - Then run `programs/simulation_table_and_figures.R`
- Run `liberia_preprocess_run.sh`or:
  - Run `programs/liberia_replicate_original_table2b.R`
  - Then run `programs/liberia_subset_process_data.R`
  - Finally, run `programs/liberia_variable_tables.R`
- Run `programs/liberia_compare_methods_coveriate_sets.R`

### Details on various programs

#### Running Shell Scripts (*Optional*)
To simply rerun the simulation code, you can run the following in a Linux terminal.

```bash
cd dp-for-rct-paper-replication

chmod +x run_simulations.sh
./run_simulations.sh

```

Similarly for cleaning the Blattman et al. (2025) data and generating basic tables.

```bash
cd dp-for-rct-paper-replication

chmod +x liberia_preprocess_run.sh
./run_simulations.sh

```


## List of tables and programs


> INSTRUCTIONS: Your programs should clearly identify the tables and figures as they appear in the manuscript, by number. Sometimes, this may be obvious, e.g. a program called "`table1.do`" generates a file called `table1.png`. Sometimes, mnemonics are used, and a mapping is necessary. In all circumstances, provide a list of tables and figures, identifying the program (and possibly the line number) where a figure is created.
>
> NOTE: If the public repository is incomplete, because not all data can be provided, as described in the data section, then the list of tables should clearly indicate which tables, figures, and in-text numbers can be reproduced with the public material provided.

The provided code reproduces:

- [X] All numbers provided in text in the paper
- [X] All tables and figures in the paper
- [ ] Selected tables and figures in the paper, as explained and justified below.


| Figure/Table #    | Program                  | Line Number | Output file                      | Note                            |
|-------------------|--------------------------|-------------|----------------------------------|---------------------------------|
| Table 1           | 02_analysis/table1.do    |             | summarystats.csv                 ||
| Table 2           | 02_analysis/table2and3.do| 15          | table2.csv                       ||
| Table 3           | 02_analysis/table2and3.do| 145         | table3.csv                       ||
| Figure 1          | n.a. (no data)           |             |                                  | Source: Herodus (2011)          |
| Figure 2          | 02_analysis/fig2.do      |             | figure2.png                      ||
| Figure 3          | 02_analysis/fig3.do      |             | figure-robustness.png            | Requires confidential data      |

## References

Webb, Kaitlyn R., Mukherjee, Soumya, Mustafi, Aratrika,  Slavković, Aleksandra, and Vilhuber, Lars. 2026. "Assessing Utility of Differential Privacy for RCTs." https://arxiv.org/html/2309.14581v2

Blattman, Christopher, Jamison, Julian C., and Sheridan, Margaret. 2025. "Replication data for: Reducing Crime and Violence: Experimental Evidence from Cognitive Behavioral Therapy in Liberia."" Nashville, TN: American Economic Association [publisher]. Ann Arbor, MI: Inter-university Consortium for Political and Social Research [distributor], 2025-12-21. https://doi.org/10.3886/E113056V2

Blattman, Christopher, Jamison, Julian C., and Sheridan, Margaret. 2017. "Reducing Crime and Violence: Experimental Evidence from Cognitive Behavioral Therapy in Liberia." American Economic Review 107 (4): 1165–1206.
---

## Acknowledgements

Some content on this page was copied from [Hindawi](https://www.hindawi.com/research.data/#statement.templates). Other content was adapted  from [Fort (2016)](https://doi.org/10.1093/restud/rdw057), Supplementary data, with the author's permission.
