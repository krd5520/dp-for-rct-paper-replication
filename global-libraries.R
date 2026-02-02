####################################
# global libraries used everywhere #
####################################
mran.date <- "2024-12-17"

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

global.libraries <- c("dplyr", #version 1.1.4
                      "devtools", #version 2.4.6
                      "rprojroot", # version 2.1.1
                      "ggplot2", #version 4.0.1
                      "RColorBrewer", #version 1.1-3
                      "knitr", #version 1.51
                      "kableExtra", #version 1.4.0
                      "sandwich", #version 3.1-1
                      "markdown", #version 2.0
                      "haven", #version 2.5.5
                      "labelled" #version 2.16.0
)

results <- sapply(as.list(global.libraries), pkgTest)
cbind(libraries,results)

# installing an additional package, version used in analysis
remotes::install_github("krd5520/DPrct/tree/660e72c62ddcca17e6cdda29712bc8b34c745651") #https://github.com/krd5520/DPrct


# if (!require("grateful",character.only = TRUE)) {
#   # install.packages("remotes")
#   remotes::install_github("Pakillo/grateful")
#   if(!require("grateful",character.only = TRUE)) stop("Package not found")
# }
