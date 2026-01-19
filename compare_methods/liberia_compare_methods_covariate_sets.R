devtools::install_github("krd5520/DPrct")
#renv::install("krd5520/DPrct",verbose=F,lock=T)

basepath = rprojroot::find_rstudio_root_file()
main.start=proc.time()



source(paste0(basepath,"/functions/directory_setup.R"))
source(paste0(basepath,"/functions/main_liberia_compare_methods.R"))

####### Testing the multivariate histogram method on the Liberia Data
#file path and file location parameters
basepath = rprojroot::find_rstudio_root_file()
file.suffix="0123"
outputpath=paste0(basepath,"/output/liberia")
datapath=paste0(basepath,"/data/liberia_subset_nolabels.Rda")

#####################
use.allobs=FALSE #use the data set that keeps possible error in odd observation
use.obs1=TRUE #use the data set that sets odd observation to value 1
## otherwise use the data set that removes the odd observation
###########

###### Privacy global paramaters ####3
bins.param=2/3 #bins for multivariate histograms
synthdata.budget=c(1,0) #overall budget epsilon, delta
use.continuous.noise=TRUE #whether uniform noise should be added to midpoint of continuous variables

#### parameters for DP Model- Based
# privacy budget variables #
mv.prop.budget=0.5
## within each response, the residual variance gets 10% of the budget
win.y.resvar=0.1
##      the treatment coefficients get 30% of the budget (i.e. 10% for each)
win.y.trcoef=0.3
##
# other params #
bound.means=100
bound.sds=c(2^(-15),2^(15))
n.iters=5000
range.alpha=0.05 #for range function
ci.confidence=0.95
########

each.rseed=2
############ END PARAMETERS ##################
dir.set=directory_setup(basepath,outputpath,list("liberia/liberia_checkpoints"))
out=main_liberia_compare_methods(datapath=datapath,file.suffix=file.suffix,
                                 outputpath=outputpath,
                                 bins.param=bins.param,
                                 synthdata.budget=synthdat.budget,
                                 use.continuous.noise=use.continuous.noise,
                                 mv.prop.budget=mv.prop.budget,win.y.resvar=win.y.resvar,
                                 win.y.trcoef=win.y.trcoef,
                                 bound.means=bound.means,
                                 bound.sds=bound.sds,
                                 n.iters=n.iters,
                                 range.alpha=range.alpha,
                                 ci.confidence=ci.confidence,
                                 each.rseed=each.rseed,
                                 use.allobs=use.allobs,
                                 use.obs1=use.obs1)
