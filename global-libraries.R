####################################
# global libraries used everywhere #
####################################
mran.date <- "2026-02-02"

get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

if (get_os()=="linux") {
## Rstudio Package Manager
  #if (getOption("repos")["CRAN"]=="@CRAN@") {
  options(repos = c(REPO_NAME = paste0("https://packagemanager.posit.co/cran/__linux__/focal/",mran.date)))
  #} else {
  #message("Repo for CRAN already set")
  #}
} else {
## MRAN
  options(repos=paste0("https://packagemanager.posit.co/cran/",mran.date))
}

getOption("repos")["CRAN"]

pkgTest <- function(x)
{
	if (!require(x,character.only = TRUE))
	{
		install.packages(x,dep=TRUE)
		if(!require(x,character.only = TRUE)) stop("Package not found")
	}
	return("OK")
}

global.libraries <- c("remotes", #version 2.5.0
                      "dplyr", #version 1.1.4
                      "devtools", #version 2.4.6
                      "rprojroot", # version 2.1.1
                      "ggplot2", #version 4.0.1
                      "RColorBrewer", #version 1.1-3
                      "cowplot",
                      "grafify",
                      "knitr", #version 1.51
                      "kableExtra", #version 1.4.0
                      "sandwich", #version 3.1-1
                      "markdown", #version 2.0
                      "haven", #version 2.5.5
                      "labelled" #version 2.16.0
)

results <- sapply(as.list(global.libraries), pkgTest)
cbind(global.libraries,results)

# installing an additional package, version used in analysis
remotes::install_github("krd5520/DPrct@f5dc482f0acb7e609c4667f82b505124db627ffd")#332e40827d7400cdfaf968fb821104d526bae448")#332e40827d7400cdfaf968fb821104d526bae448")#660e72c62ddcca17e6cdda29712bc8b34c745651")#f5dc482f0acb7e609c4667f82b505124db627ffd")#1a5f8eccdeb16bdfd138617778df68b5f2d7f6ce")#f34bb688f4dc7379e4feabdd04b18b6f7cf4e189")#f826a1969dec478029e7a62177d7c8d7c8c499b5")#2479df9fafbc622b7c270e9910f688df0ca77a6f") #https://github.com/krd5520/DPrct
#remotes::install_github("krd5520/DPrct@3c4fb176508891f5fe8cf1a1469c0edf7df537ec")
#remotes::install_github("krd5520/DPrct@2a6b3d35d59154fe9cf36d4c7129a21bd2f7926c")



# if (!require("grateful",character.only = TRUE)) {
#   # install.packages("remotes")
#   remotes::install_github("Pakillo/grateful")
#   if(!require("grateful",character.only = TRUE)) stop("Package not found")
# }
