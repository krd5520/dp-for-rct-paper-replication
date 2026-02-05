#devtools::install_github("krd5520/DPrct")

#renv::install("krd5520/DPrct",verbose=F,lock=T)
main.start=proc.time()

####### Testing the multivariate histogram method on the Liberia Data
#file path and file location parameters
basepath = rprojroot::find_rstudio_root_file()
file.suffix="v0321_all_eps2_binssixth_nostandardize"
outputpath=paste0(basepath,"/output")
in.datapath=paste0(basepath,"/data/liberia_subset_nolabels.Rda")
out.datapath=paste0(basepath,"/data")
override.checkpoints=T
#####################
use.allobs=FALSE #use the data set that keeps possible error in odd observation
use.obs1=TRUE #use the data set that sets odd observation to value 1
## otherwise use the data set that removes the odd observation
###########

###### Privacy global paramaters ####3
bins.param=1/6 #bins for multivariate histograms
synthdata.budget=c(2,0) #overall budget epsilon, delta
use.continuous.noise=TRUE #whether uniform noise should be added to midpoint of continuous variables
use.genmod.full=T
use.mvhist.full=T
which.sets=c("Reduced","Subset","Full")

#### parameters for DP Model- Based
# privacy budget variables #
mv.prop.budget=0.5
## within each response, the residual variance gets 20% of the budget
win.y.resvar=0.1
##      the treatment coefficients get 30% of the budget (i.e. 10% for each)
win.y.trcoef=0.3
##
# other params #
bound.means=50
bound.sds=c(2^(-15),2^(15))
n.iters=5000
range.alpha=0.05 #for range function
ci.confidence=0.95

#regression diagnostic plot properties
reg.assumption.plots=T
reg.plot.props=list("pt.sz"=0.3,"ln.sz"=0.7,"bs.sz"=9,"ncol"=2,"width"=8,"height"=7)

#response.vars.all = c('fam_asb_lt',
#                      paste0(c('stealnb','carryweapon','asbhostilstd'),"_ltav"))

#### Other Budgets for Full DP
### otherbudgetsDP list with each element representing a budget. Each element should have mv.epsilon,mv.delta,per.y.eps,per.y.del
otherbudgetsDP=list("DP-Mb Fix MV Budget Double Overall"=c(synthdata.budget,1,0),
                    "DP-Mb Fix MV Budget to Each Y Budget"=c(synthdata.budget,2,0))

########
each.rseed=1
standardize.covar=F
all.standardize.covar=F

## Bounds parameters
load(paste0(outputpath,"/liberia/liberia_bounds.Rda"))
continuous.vars=names(bounds.list)

if(all.standardize.covar==T){
  standardize.vars=var.df$variable[!grepl("Binary|Standard",var.df$type)]
  standardize.vars=continuous.vars[continuous.vars%in%standardize.vars]
  not.std=lapply(var.df$variable[grepl("Standard",var.df$type)],function(x)c(var.df$min[var.df$variable==x],var.df$max[var.df$variable==x]))
  names(not.std)=var.df$variable[grepl("Standard",var.df$type)]
  yes.std=rep(list(c(-n.stdev,n.stdev)),length(standardize.vars))
  names(yes.std)=standardize.vars
  continuous.limits=c(yes.std,not.std)
}else if(standardize.covar==T){
standardize.vars=var.df$variable[var.df$addstandardize==T]
continuous.limits=bounds.std.list
}else{
  standardize.vars=NULL
  continuous.limits=bounds.list
}
n.std.dev=n.stdev


time.outputs.file=paste0(outputpath,"/liberia/liberia_diagnostic_",file.suffix,".txt")
param.print=c(paste0("Add continuous noise: ",use.continuous.noise,"; Standardize Some Covariates: ",
                     standardize.covar,";Use GenModel on Full Covariates: ",use.genmod.full,"; each.rseed: ",each.rseed),
              paste0("Continuous Limits: ",paste(names(continuous.limits),"=",continuous.limits)),
              paste0("Standardized Deviations Limit=+/-",n.std.dev),
              paste0("General Parameters: bin concentration=",bins.param," mean bounds=",bound.means,"; s.e. bounds=",paste0(round(bound.sds,4),collapse=", "),
                     "; range alpha=",range.alpha,"; Number of Proxies B=",n.iters),
              paste0("Base Privacy Budget: Epsilon=",synthdata.budget[1]," Delta=",synthdata.budget[2]),
              paste("Additional Budgets for GenModel:\n",paste(names(otherbudgetsDP),"=",otherbudgetsDP)),
              paste0("GenMod Privacy Allocations:"," MV Histogram proportion=",mv.prop.budget,
                     "; Residual Variance proportion=",win.y.resvar,
                     "; Treatment Effect proportion=",win.y.trcoef))


## set up timeout file
if(file.exists(time.outputs.file)==FALSE){
  print('making new file')
  writeLines(c(paste0(rep("#",20),collapse=""),file.suffix,paste(c("user:","system:","elapsed:"),proc.time(),collapse=", "),paste0(rep("#",20),collapse=""),param.print),time.outputs.file)
}else{
  CON=file(time.outputs.file,"a")
  writeLines(c(" "," ",paste0(rep("#",20),collapse=""),file.suffix,paste(c("user:","system:","elapsed:"),proc.time(),collapse=", "),paste0(rep("#",20),collapse=""),param.print),CON)
  close(CON)
}
############ END PARAMETERS ##################
source(paste0(basepath,"/program/functions/directory_setup.R"))
source(paste0(basepath,"/program/functions/main_liberia_compare_methods.R"))
source(paste0(basepath,"/program/functions/liberia_compare_methods_onecovset.R"))
source(paste0(basepath,"/program/functions/reg_assumptions_function.R"))

dir.set=directory_setup(basepath,outputpath,list("liberia/liberia_checkpoints",
                                                 "figures",paste0("figures/liberia_",file.suffix)))

out=main_liberia_compare_methods(in.datapath=in.datapath,file.suffix=file.suffix,
                                 outputpath=outputpath,out.datapath=out.datapath,
                                 bins.param=bins.param,
                                 synthdata.budget=synthdata.budget,
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
                                 use.obs1=use.obs1,
                                 use.genmod.full=use.genmod.full,
                                 otherbudgetsDP=otherbudgetsDP,
                                 override.checkpoints = override.checkpoints,
                                 reg.plot.props=reg.plot.props,
                                 reg.assumption.plots=reg.assumption.plots,
                                 n.std.dev=n.std.dev,continuous.limits=continuous.limits,
                                 continuous.vars=continuous.vars,
                                 standardize.vars=standardize.vars,
                                 diagnostic.file=time.outputs.file,which.sets=which.sets,var.df=var.df,cat.vars=cat.vars)

comp.time=(proc.time()-main.start)[[3]]
hours=comp.time%/%(60*60)
minutes=(comp.time%%(60*60))%/%60
seconds=(comp.time%%(60*60))%%60
print(paste("Computation Time:",hours,"hours,",minutes,"minutes and",seconds,"seconds."))
