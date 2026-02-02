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

To simply rerun the simulation code, you can run the following in a Linux terminal.
```bash
cd dp-for-rct-paper-replication

chmod +x run_simulations.sh
./run_simulations.sh

```

Alternatively, you can open `DPrct.Rproj` in Rstudio or similar program and run the R files in the following order: 
  1. global_libraries.R
  2. 

## Directories Overview
Here is a brief visual overview of the project repository:
```
dp-for-rct-paper-replication/
├── compare_methods/
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
│   ├── DataDiag/					# Contains intermediate datasets and diagnostic files
│   ├── NAICS6_Pyfunctions/			# Contains python helper libraries
│   ├── GUIcode/			# Work in-progress GUI creation code to walk a user through creating the config file
│   ├── generateMicrodata.py		# The main script
│   ├── main.py		              # The main script
│   └── config.yaml					# Contains configurable parameters and model selections
└── README.md
```


## Data Availability and Provenance Statements


### Statement about Rights

- [X] I certify that the author(s) of the manuscript have legitimate access to and permission to use the data used in this manuscript. 
- [X] I certify that the author(s) of the manuscript have documented permission to redistribute/publish the data contained within this replication package. Appropriate permission are documented in the LICENSE.txt file file.


### License for Code

The code, scripts, and programs are licensed under a Modified BSD License. All other objects are licensed under Creative Commons/CC-BY-NC license. See LICENSE.txt for details.


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
  -`dtafiles/STYL_Final.dta`: Save this file as `STYL_FINAL` in directory `data/liberia_replication_data`. 
  - `dofiles/AER_1.do`: Includes analysis details that are references as comments in `program/liberia_replicate_original_table2b.R` to justify variable selection. The file can be viewed through Stata or Notepad text reader.
  - `dofiles/AER_2.do`: Includes analysis details that are references as comments in `program/liberia_replicate_original_table2b.R` to justify variable selection. The file can be viewed through Stata or Notepad text reader
- `LICENSE.txt`: Details the packages redistribution restrictions and its Modified BSD License and Creative Commons Attribution 4.0 International Public License.


We accessed and downloaded the data January 2026 (corresponding to the updates Blattman et al. (2025) made to the package in December 2025). 

Blattman, Christopher, Jamison, Julian C., and Sheridan, Margaret. Replication data for: Reducing Crime and Violence: Experimental Evidence from Cognitive Behavioral Therapy in Liberia. Nashville, TN: American Economic Association [publisher], 2025. Ann Arbor, MI: Inter-university Consortium for Political and Social Research [distributor], 2025-12-21. https://doi.org/10.3886/E113056V2

Datafile: `data/liberia_replication_data/STYL_Final.dta` (not provided)


## Computational requirements

The required R packages can be installed by running `global-libraries.R`. 
The CRAN packages required are as following with the version we used for our analysis:
 
  
Additionally, a package from GitHub (repository `krd5520/DPrct`) is used for the proposed mechanisms. It can be installed using the following R code. 

``` R
remotes::install_github("krd5520/DPrct@2479df9fafbc622b7c270e9910f688df0ca77a6f") #https://github.com/krd5520/DPrct

```


### Software Requirements

> INSTRUCTIONS: List all of the software requirements, up to and including any operating system requirements, for the entire set of code. It is suggested to distribute most dependencies together with the replication package if allowed, in particular if sourced from unversioned code repositories, Github repos, and personal webpages. In all cases, list the version *you* used. All packages should be listed in human-readable form in this README, but should also be included in a setup or install script.

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

> INSTRUCTIONS: Memory and compute-time requirements may also be relevant or even critical. Some example text follows. It may be useful to break this out by Table/Figure/section of processing. For instance, some estimation routines might run for weeks, but data prep and creating figures might only take a few minutes. You should also describe how much storage is required in addition to the space visible in the typical repository, for instance, because data will be unzipped, data downloaded, or temporary files written.

#### Summary time to reproduce

Approximate time needed to reproduce the analyses on a standard (CURRENT YEAR) desktop machine:

- [ ] <10 minutes
- [ ] 10-60 minutes
- [ ] 1-2 hours
- [ ] 2-8 hours
- [ ] 8-24 hours
- [ ] 1-3 days
- [ ] 3-14 days
- [ ] > 14 days

#### Summary of required storage space

Approximate storage space needed:

- [ ] < 25 MBytes
- [ ] 25 MB - 250 MB
- [ ] 250 MB - 2 GB
- [ ] 2 GB - 25 GB
- [ ] 25 GB - 250 GB
- [ ] > 250 GB

- [ ] Not feasible to run on a desktop machine, as described below.

#### Computational Details

The code was last run on a **4-core Intel-based laptop with MacOS version 10.14.4 with 200GB of free space**. 

Portions of the code were last run on a **32-core Intel server with 1024 GB of RAM, 12 TB of fast local storage**. Computation took **734 hours**. 

Portions of the code were last run on a **12-node AWS R3 cluster, consuming 20,000 core-hours, with 2TB of attached storage**.  

> INSTRUCTIONS: Identifiying hardware and OS can be obtained through a variety of ways:
> Some of these details can be found as follows:
>
> - (Windows) by right-clicking on "This PC" in File Explorer and choosing "Properties"
> - (Mac) Apple-menu > "About this Mac"
> - (Linux) see code in [linux-system-info.sh](https://github.com/AEADataEditor/replication-template/blob/master/tools/linux-system-info.sh)`


## Description of programs/code

> INSTRUCTIONS: Give a high-level overview of the program files and their purpose. Remove redundant/ obsolete files from the Replication archive.

- Programs in `programs/01_dataprep` will extract and reformat all datasets referenced above. The file `programs/01_dataprep/main.do` will run them all.
- Programs in `programs/02_analysis` generate all tables and figures in the main body of the article. The program `programs/02_analysis/main.do` will run them all. Each program called from `main.do` identifies the table or figure it creates (e.g., `05_table5.do`).  Output files are called appropriate names (`table5.tex`, `figure12.png`) and should be easy to correlate with the manuscript.
- Programs in `programs/03_appendix` will generate all tables and figures  in the online appendix. The program `programs/03_appendix/main-appendix.do` will run them all. 
- Ado files have been stored in `programs/ado` and the `main.do` files set the ADO directories appropriately. 
- The program `programs/00_setup.do` will populate the `programs/ado` directory with updated ado packages, but for purposes of exact reproduction, this is not needed. The file `programs/00_setup.log` identifies the versions as they were last updated.
- The program `programs/config.do` contains parameters used by all programs, including a random seed. Note that the random seed is set once for each of the two sequences (in `02_analysis` and `03_appendix`). If running in any order other than the one outlined below, your results may differ.

### (Optional, but recommended) License for Code

> INSTRUCTIONS: Most journal repositories provide for a default license, but do not impose a specific license. Authors should actively select a license. This should be provided in a LICENSE.txt file, separately from the README, possibly combined with the license for any data provided. Some code may be subject to inherited license requirements, i.e., the original code author may allow for redistribution only if the code is licensed under specific rules - authors should check with their sources. For instance, some code authors require that their article describing the econometrics of the package be cited. Licensing can be complex. Some non-legal guidance may be found [here](https://social-science-data-editors.github.io/guidance/Licensing_guidance.html).

The code is licensed under a MIT/BSD/GPL [choose one!] license. See LICENSE.txt file for details.

## Instructions to Replicators

> INSTRUCTIONS: The first two sections ensure that the data and software necessary to conduct the replication have been collected. This section then describes a human-readable instruction to conduct the replication. This may be simple, or may involve many complicated steps. It should be a simple list, no excess prose. Strict linear sequence. If more than 4-5 manual steps, please wrap a main program/Makefile around them, in logical sequences. Examples follow.

- Edit `programs/config.do` to adjust the default path
- Run `programs/00_setup.do` once on a new system to set up the working environment. 
- Download the data files referenced above. Each should be stored in the prepared subdirectories of `data/`, in the format that you download them in. Do not unzip. Scripts are provided in each directory to download the public-use files. Confidential data files requested as part of your FSRDC project will appear in the `/data` folder. No further action is needed on the replicator's part.
- Run `programs/01_main.do` to run all steps in sequence.

### Details on various programs

- `programs/00_setup.do`: will create all output directories, install needed ado packages. 
   - If wishing to update the ado packages used by this archive, change the parameter `update_ado` to `yes`. However, this is not needed to successfully reproduce the manuscript tables. 
- `programs/01_dataprep`:  
   - These programs were last run at various times in 2018. 
   - Order does not matter, all programs can be run in parallel, if needed. 
   - A `programs/01_dataprep/main.do` will run them all in sequence, which should take about 2 hours.
- `programs/02_analysis/main.do`.
   - If running programs individually, note that ORDER IS IMPORTANT. 
   - The programs were last run top to bottom on July 4, 2019.
- `programs/03_appendix/main-appendix.do`. The programs were last run top to bottom on July 4, 2019.
- Figure 1: The figure can be reproduced using the data provided in the folder “2_data/data_map”, and ArcGIS Desktop (Version 10.7.1) by following these (manual) instructions:
  - Create a new map document in ArcGIS ArcMap, browse to the folder
“2_data/data_map” in the “Catalog”, with files  "provinceborders.shp", "lakes.shp", and "cities.shp". 
  - Drop the files listed above onto the new map, creating three separate layers. Order them with "lakes" in the top layer and "cities" in the bottom layer.
  - Right-click on the cities file, in properties choose the variable "health"... (more details)

## List of tables and programs


> INSTRUCTIONS: Your programs should clearly identify the tables and figures as they appear in the manuscript, by number. Sometimes, this may be obvious, e.g. a program called "`table1.do`" generates a file called `table1.png`. Sometimes, mnemonics are used, and a mapping is necessary. In all circumstances, provide a list of tables and figures, identifying the program (and possibly the line number) where a figure is created.
>
> NOTE: If the public repository is incomplete, because not all data can be provided, as described in the data section, then the list of tables should clearly indicate which tables, figures, and in-text numbers can be reproduced with the public material provided.

The provided code reproduces:

- [ ] All numbers provided in text in the paper
- [ ] All tables and figures in the paper
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

> INSTRUCTIONS: As in any scientific manuscript, you should have proper references. For instance, in this sample README, we cited "Ruggles et al, 2019" and "DESE, 2019" in a Data Availability Statement. The reference should thus be listed here, in the style of your journal:

Steven Ruggles, Steven M. Manson, Tracy A. Kugler, David A. Haynes II, David C. Van Riper, and Maryia Bakhtsiyarava. 2018. "IPUMS Terra: Integrated Data on Population and Environment: Version 2 [dataset]." Minneapolis, MN: *Minnesota Population Center, IPUMS*. https://doi.org/10.18128/D090.V2

Department of Elementary and Secondary Education (DESE), 2019. "Student outcomes database [dataset]" *Massachusetts Department of Elementary and Secondary Education (DESE)*. Accessed January 15, 2019.

U.S. Bureau of Economic Analysis (BEA). 2016. “Table 30: "Economic Profile by County, 1969-2016.” (accessed Sept 1, 2017).

Inglehart, R., C. Haerpfer, A. Moreno, C. Welzel, K. Kizilova, J. Diez-Medrano, M. Lagos, P. Norris, E. Ponarin & B. Puranen et al. (eds.). 2014. World Values Survey: Round Six - Country-Pooled Datafile Version: http://www.worldvaluessurvey.org/WVSDocumentationWV6.jsp. Madrid: JD Systems Institute.

---

## Acknowledgements

Some content on this page was copied from [Hindawi](https://www.hindawi.com/research.data/#statement.templates). Other content was adapted  from [Fort (2016)](https://doi.org/10.1093/restud/rdw057), Supplementary data, with the author's permission.
