### Simulations Code
#options(verbose=F)
#install.packages("devtools","tidyr","ggplot2")
devtools::install_github("krd5520/DPrct")

library(tidyr)
library(dplyr)


####### Testing the multivariate histogram method on the Liberia Data
basepath = rprojroot::find_rstudio_root_file()
main.start=proc.time()

source(paste0(basepath,"/program/functions/generate_simulated_dataset_function.R"))
source(paste0(basepath,"/program/functions/simple_simulation_functions.R"))
source(paste0(basepath,"/program/functions/reg_assumptions_function.R"))
source(paste0(basepath,"/program/functions/sim_compare_method_function.R"))
source(paste0(basepath,"/program/functions/main_simulation_modelcompare.R"))
source(paste0(basepath,"/program/functions/main_simulation_budgetcompare.R"))
source(paste0(basepath,"/program/functions/directory_setup.R"))
source(paste0(basepath,"/program/functions/main_simulation_nonoise.R"))
source(paste0(basepath,"/program/functions/sim_nonoise_dpmb.R"))


################## FILE PARAMETERS ##############################
output.folder=paste0(basepath,"/output")
file.suffix="0320_loweffect" #"v0225_ncat" #added to the files names for all generated files
#suffixes for saving the checkpoint in simulation 1,
#    and saving the check and diagnostic plots in simulation 2
sim1.conf.suffix="conf_sim_budgets"
sim2.conf.suffix="conf_sim_models"
save.plot.png=T

dir.set=directory_setup(basepath,output.folder,list("simulations","figures","tables"))
time.outputs.file=paste0(output.folder,"/simulations/time_out_",file.suffix,".txt")

## In Simulation 2, diagnostic plots are saved (QQ Norm Plot, and Residual vs. Fit)
# pt.sz is size of points, ln.sz is size of line, bs.sz is base size for text
# width1 and height1 are for the size of the png saved for the standard models
# additional models are saved using width2 and height2
conf.plot=list("pt.sz1"=0.3,"ln.sz1"=0.7,"bs.sz1"=8,
               "pt.sz2"=0.4,"ln.sz2"=0.9,"bs.sz2"=9,
               "ncol1"=3,"ncol2"=2,
               "width1"=8.5,"height1"=5,
               "width2"=5,"height2"=2)


################### SET OVERALL PARAMETERS #####################
nsims=20 #number of repetitions
nobs=1000 #number of observations
ntreat=2 #number of treatments (including control)
#treat.effect=5 #treatment effect
treat.effect=-0.25
intercept=1 #model intercept
#resid.sd=2 #residual standard deviation
resid.sd=0.87
se.normal=T #use central limit theory to treat the collection of se values from proxies as normal sample.

#names for treatment and response variables
response.vars="y"
treat.vars="t1"

cat.probs=list(c(1/2,1/2))
#geometric sequence generated coefficients are defined by the start and ratio
#parameters. If start=c and ratio=r, the coefficients are c, c*(r^1), c*(r^2),....
cat.gseq=c(3,7/11)#(start,ratio)
cont.gseq=c(0.99,2/3)#(start,cont)
finite.bounds=c(1,1.2)
cont.funcs=list(stats::rnorm,runif.round)
cont.params=list(list("n"=nobs,"mean"=0,"sd"=2),
                 list("n"=nobs,"min"=finite.bounds[1],"max"=finite.bounds[2],"d"=2))
n.std.dev=5
standardize.y=F
standardize.x=T
std.bounds=c(-n.std.dev,n.std.dev)

#### Privacy Mechanism Parameters ######
### other global parameters
bins.param=1/6 #bins for multivariate histograms

use.continuous.noise=F #whether uniform noise should be added to midpoint of continuous variables

#### parameters for Full DP Model- Based
# privacy budget variables #
mv.prop.budget=0.5
## within each response, the residual variance gets 10% of the budget
win.y.resvar=0.1
##      the treatment mod.coefs get 30% of the budget (i.e. 10% for each)
win.y.trcoef=0.1
##
# other params #
bound.means=50
bound.sds=c(2^(-15),2^(15))
n.iters.fulldp=5000
range.alpha=0.05

### Analysis Parameters
stderr.func=function(mod,std.type="HC1"){
  sqrt(diag(sandwich::vcovHC(mod,type=std.type)))
}
ci.confidence=0.95
############


############ Compare Privacy Budget Parameters #################
sim1.simdata.rseed=2
sim1.san.rseed=2
#overall budgets epsilon, delta
sim1.synthdata.budget=list("MassiveEps"=c(5000,0),
                      "PureDP_EpsHalf"=c(1/2,0),
                      #"ApproxDP_EpsHalf"=c(1/2,1/nobs),
                      "PureDP_Eps1"=c(1,0),
                      #"ApproxDP_Eps1"=c(1,1/nobs),
                      "PureDP_Eps2"=c(2,0),
                      #"PureDP_Eps3"=c(3,0),
                      "PureDP_Eps4"=c(4,0),
                      #"PureDP_Eps5"=c(5,0),
                      "PureDP_Eps6"=c(6,0))
sim1.ncat=4
sim1.ncont=4
response.limits1=c(-10,25)



########### Compare Regression Model Parameters ###################
#combinations of number of categorical and continuous covariates
sim2.simdata.rseed=2
sim2.san.rseed=1
combos.NCov=data.frame("NCat"=rep(c(2,4,6),each=3),
                       "NCont"=rep(c(2,4,6),3),
                       "Name"=paste("Model",seq(1,9)))
sim2.eps=1
sim2.delta=0
dpmodel.higher.budget=TRUE
response.limits2=c(-10,25)

### Define Additional Regression Models
mx.ncont=max(combos.NCov$NCont)
add.ncat=2
add.ncont=mx.ncont%/%2
cat.formula=paste0("x",seq(mx.ncont+1,mx.ncont+add.ncat),collapse="+")
base.formula=paste0(response.vars,"~",paste0(treat.vars,collapse="+"))
#additional formulas one with only the unbounded continuous and another with only the finite continuous
add.formulas=NULL#c(
  #paste0(base.formula,"+",paste0("x",seq(1,mx.ncont-1,2),collapse="+"),"+",cat.formula))#,
  #paste0(base.formula,"+",paste0("x",seq(2,mx.ncont,2),collapse="+"),"+",cat.formula))
#names(add.formulas)=c("Unbounded")#,"Finite Set")


num.cores=1#parallel::detectCores()-1
##cl=parallel::makeCluster(num.cores)
##doParallel::registerDoParallel(cl)

param.print=c(paste0("Repetitions: ",nsims,"; Observations:",nobs,"; Add continuous noise: ",use.continuous.noise,"; Standardize y: ",standardize.y),
              paste0("Continuous Limits: Finite=",paste0(finite.bounds,collapse=", "),"; Standardized Deviations Limit=+/-",n.std.dev,
                     "; y=",paste0(response.limits1,collapse=", "),"(sim1) ",paste0(response.limits2,collapse=", "),"(sim2)"),
              paste0("General Parameters: bin concentration=",bins.param," mean bounds=",bound.means,"; s.e. bounds=",paste0(round(std.bounds,4),collapse=", "),
                     "; range alpha=",range.alpha,"; Number of Proxies B=",n.iters.fulldp),
              paste0("GenMod Privacy Allocations:"," MV Histogram proportion=",mv.prop.budget,"; Residual Variance proportion=",win.y.resvar,"; Treatment Effect proportion=",win.y.trcoef))


## set up timeout file
if(file.exists(time.outputs.file)==FALSE){
  print('making new file')
  writeLines(c(paste0(rep("#",20),collapse=""),file.suffix,paste(c("user:","system:","elapsed:"),proc.time(),collapse=", "),paste0(rep("#",20),collapse=""),param.print),time.outputs.file)
}else{
  CON=file(time.outputs.file,"a")
  writeLines(c(" "," ",paste0(rep("#",20),collapse=""),file.suffix,paste(c("user:","system:","elapsed:"),proc.time(),collapse=", "),paste0(rep("#",20),collapse=""),param.print),CON)
  close(CON)
}


######## RUN Simulations ###############

###### SIM 1: Budget ######

#get continuous limits, continuous variables, and standardize variables based on inputs
#cont.limits1=rep(list(std.bounds,finite.bounds),sim1.ncont%/%2)


cont.vars.all1=c(paste0("x",seq(1,sim1.ncont)),response.vars)

if(standardize.x==TRUE){
  std.vars.all1=c(paste0("x",seq(1,sim1.ncont)))
  cont.limits1=rep(list(std.bounds),sim1.ncont)
}else{
  std.vars.all1=c(paste0("x",seq(1,sim1.ncont,2)))
  cont.limits1=rep(list(std.bounds,finite.bounds),sim1.ncont%/%2)
  if(sim1.ncont%%2==1){
   cont.limits1=c(cont.limits1,list(std.bounds))
  }
}
if(standardize.y==FALSE){
  cont.limits1=c(cont.limits1,list(response.limits1))
}else{
  cont.limits1=c(cont.limits1,rep(list(std.bounds),length(response.vars)))
  std.vars.all1=c(std.vars.all1,response.vars)
}
names(cont.limits1)=cont.vars.all1

sim1.param.message=c(" "," ",paste(paste0(rep("#",5),collapse=""),"Simulation 1: Budget (confidential data seed=",sim1.simdata.rseed,", sanitize data seeds=",sim1.san.rseed,")",paste0(rep("#",5),collapse="")),
  paste0("Sim 1 Model:\n ",paste0(response.vars,"~",intercept," + ",paste0(treat.effect,"*T",seq(1,ntreat),collapse=" + ")," + ",
                              paste0(round(cont.gseq[1]*(cont.gseq[2]^seq(1,sim1.ncont)),2),"*a",seq(1,sim1.ncont),collapse=" + ")," + ",
                              paste0(round(cat.gseq[1]*(cat.gseq[2]^seq(1,sim1.ncat)),2),"*b",seq(1,sim1.ncat),collapse=" + ")," + ",
                              " N(0,",round(resid.sd^2,2),")")),
  "where 't' variables are treatments, 'a' variables are continuous, and 'b' variables are binary.\n",
  "Budgets: ",
  paste0(names(sim1.synthdata.budget),":",sapply(sim1.synthdata.budget,function(x)paste0("(epsilon, delta)=(",x[1],", ",x[2],")"))))

CON=file(time.outputs.file,"a")
writeLines(sim1.param.message,CON)
close(CON)

sink(type="message")
sim1=main_simulation_budgetcompare(nsims=nsims,nobs=nobs,response=response.vars,treat.vars=treat.vars,
                                   ncat=sim1.ncat,ncont=sim1.ncont,
                                       intercept=intercept,treat.effect=treat.effect,
                                       cat.gseq=cat.gseq,cont.gseq=cont.gseq,resid.sd=resid.sd,
                                       cat.probs=cat.probs,cont.funcs=cont.funcs,cont.params=cont.params,
                                       synthdata.budget=sim1.synthdata.budget,
                                       sim.data.rseed=sim1.simdata.rseed,
                                       conf.suffix=sim1.conf.suffix,file.suffix=file.suffix,
                                       output.folder=output.folder,
                                       bins.param=bins.param,use.continuous.noise=use.continuous.noise,
                                       mv.prop.budget=mv.prop.budget,win.y.resvar=win.y.resvar,
                                       win.y.trcoef=win.y.trcoef,bound.means=bound.means,bound.sds=bound.sds,
                                       range.alpha=range.alpha,n.iters.fulldp=n.iters.fulldp,
                                       ci.confidence=ci.confidence,stderr.func=stderr.func,
                                       san.rseed=sim1.san.rseed,
                                   se.normal=se.normal,n.std.dev=n.std.dev,continuous.limits=cont.limits1,
                                   continuous.vars=cont.vars.all1,standardize.vars=std.vars.all1,diagnostic.file=time.outputs.file)

sink()
comp.time=(proc.time()-main.start)[[3]]
hours=comp.time%/%(60*60)
minutes=(comp.time%%(60*60))%/%60
seconds=(comp.time%%(60*60))%%60
sim1.time.message=paste("Sim 1 Computation Time:",hours,"hours,",minutes,"minutes and",seconds,"seconds.")
print(sim1.time.message)
CON=file(time.outputs.file,"a")
writeLines(sim1.time.message,CON)
close(CON)

save(sim1,file=paste0(output.folder,"/simulations/simulation1_",file.suffix,".Rda"))

## Note: If going to standardize the response. Need to do that before fitting the confidential model as well.

####### SIM 2: Models ####
source(paste0(basepath,"/program/functions/main_simulation_modelcompare.R"))
sim2.start=proc.time()

#get continuous limits, continuous variables, and standardize variables based on inputs


cont.vars.all2=c(paste0("x",seq(1,sim1.ncont)),response.vars)

if(standardize.x==TRUE){
  std.vars.all2=c(paste0("x",seq(1,mx.ncont)))
  cont.limits2=rep(list(std.bounds),mx.ncont)
}else{
  std.vars.all2=c(paste0("x",seq(1,mx.ncont,2)))
  cont.limits2=rep(list(std.bounds,finite.bounds),mx.ncont%/%2)
  if(mx.ncont%%2==1){
    cont.limits2=c(cont.limits2,list(std.bounds))
  }
}
if(standardize.y==FALSE){
  cont.limits2=c(cont.limits2,list(response.limits2))
}else{
  cont.limits2=c(cont.limits2,rep(list(std.bounds),length(response.vars)))
  std.vars.all2=c(std.vars.all2,response.vars)
}
names(cont.limits2)=cont.vars.all2

sim2.param.message=c(" ",paste(paste0(rep("#",5),collapse=""),"Simulation 2: Models (confidential data seed=",sim2.simdata.rseed,", sanitize data seeds=",sim2.san.rseed,")",paste0(rep("#",5),collapse="")),
                     paste0("Sim 2 Generating Model:\n ",paste0(response.vars,"~",intercept," + ",paste0(treat.effect,"*t",seq(1,ntreat),collapse=" + "),"\n + ",
                                                              paste0(round(cont.gseq[1]*(cont.gseq[2]^seq(1,mx.ncont)),2),"*a",seq(1,mx.ncont),collapse=" + "),"\n + ",
                                                              paste0(round(cat.gseq[1]*(cat.gseq[2]^seq(1,max(combos.NCov$NCat))),2),"*b",seq(1,max(combos.NCov$NCat)),collapse=" + ")," + ",
                                                              " N(0,",round(resid.sd^2,2),")")),
                     "where 't' variables are treatments, 'a' variables are continuous and 'b' variables are binary.\n",
                     paste0("Privacy Budgets: Base (epsilon,delta)=(",sim2.eps,", ",sim2.delta,")",
                            ifelse(dpmodel.higher.budget==TRUE,paste0("Higher GenMod Budget (epsilon,delta)=(",sim2.eps/mv.prop.budget,", ",sim2.delta/mv.prop.budget,")")," ")),
                     "Models Fitted: ",
                     apply(combos.NCov,1,function(rw)paste0(rw[3],": m continuous=",rw[2]," m binary=",rw[1])),
                     "Additional Models:",
                     paste0(names(add.formulas),": ",add.formulas)
                     )

CON=file(time.outputs.file,"a")
writeLines(sim2.param.message,CON)
close(CON)

sim2=main_simulation_modelcompare(combos.NCov=combos.NCov,nsims=nsims,
                             nobs=nobs,response=response.vars,treat.vars=treat.vars,
                             add.formulas=add.formulas,add.names=names(add.formulas),
                             add.ncat=add.ncat,add.ncont=add.ncont,
                             intercept=intercept,treat.effect=treat.effect,
                             cat.gseq=cat.gseq,cont.gseq=cont.gseq,resid.sd=resid.sd,
                             cat.probs=cat.probs,cont.funcs=cont.funcs,cont.params=cont.params,
                             priv.eps=sim2.eps,priv.delta=sim2.delta,
                             dpmodel.higher.budget=dpmodel.higher.budget,
                             sim.data.rseed=sim2.simdata.rseed,
                             conf.suffix=sim2.conf.suffix,file.suffix=file.suffix,
                             output.folder=output.folder,conf.plot=conf.plot,
                             bins.param=bins.param,use.continuous.noise=use.continuous.noise,
                             mv.prop.budget=mv.prop.budget,win.y.resvar=win.y.resvar,
                             win.y.trcoef=win.y.trcoef,bound.means=bound.means,bound.sds=bound.sds,
                             range.alpha=range.alpha,n.iters.fulldp=n.iters.fulldp,
                             ci.confidence=ci.confidence,stderr.func=stderr.func,
                             san.rseed=sim2.san.rseed,se.normal=se.normal,continuous.limits=cont.limits2,
                             continuous.vars=cont.vars.all2,standardize.vars=std.vars.all2,diagnostic.file=time.outputs.file)#save.plot.png=save.plot.png)
sim2_nonoise=NULL
# source(paste0(basepath,"/program/functions/generate_simulated_dataset_function.R"))
# source(paste0(basepath,"/program/functions/simple_simulation_functions.R"))
# source(paste0(basepath,"/program/functions/reg_assumptions_function.R"))
# source(paste0(basepath,"/program/functions/sim_compare_method_function.R"))
# source(paste0(basepath,"/program/functions/main_simulation_modelcompare.R"))
# source(paste0(basepath,"/program/functions/main_simulation_budgetcompare.R"))
# source(paste0(basepath,"/program/functions/directory_setup.R"))
# source(paste0(basepath,"/program/functions/main_simulation_nonoise.R"))
# source(paste0(basepath,"/program/functions/sim_nonoise_dpmb.R"))
# sim2_nonoise=main_simulation_nonoise(combos.NCov=combos.NCov,nsims=nsims,
#                                   nobs=nobs,
#                                   response=response.vars,
#                                   treat.vars=treat.vars,
#                                   add.formulas=add.formulas,
#                                   add.names=names(add.formulas),
#                                   add.ncat=add.ncat,
#                                   add.ncont=add.ncont,
#                                   intercept=intercept,
#                                   treat.effect=treat.effect,
#                                   cat.gseq=cat.gseq,
#                                   cont.gseq=cont.gseq,
#                                   resid.sd=resid.sd,
#                                   cat.probs=cat.probs,
#                                   cont.funcs=cont.funcs,
#                                   cont.params=cont.params,
#                                   sim.data.rseed=sim2.simdata.rseed,
#                                   conf.suffix=sim2.conf.suffix,
#                                   file.suffix=paste0("nonoise_",file.suffix),
#                                   output.folder=output.folder,
#                                   conf.plot=conf.plot,
#                                   bins.param=bins.param,
#                                   use.continuous.noise=use.continuous.noise,
#                                   n.iters.fulldp=n.iters.fulldp,
#                                   ci.confidence=ci.confidence,
#                                   stderr.func=stderr.func,
#                                   san.rseed=sim2.san.rseed,se.normal=se.normal)#save.plot.png=save.plot.png)
#
# save(sim2_nonoise,file=paste0(output.folder,"/simulations/simulations_nonoise",file.suffix,".Rda"))

comp.time=(proc.time()-sim2.start)[[3]]
hours=comp.time%/%(60*60)
minutes=(comp.time%%(60*60))%/%60
seconds=(comp.time%%(60*60))%%60

sim2.time.message=paste("Sim 2: Computation Time:",hours,"hours,",minutes,"minutes and",seconds,"seconds.")
print(sim2.time.message)
CON=file(time.outputs.file,"a")
writeLines(sim2.time.message,CON)
close(CON)

save(sim2,file=paste0(output.folder,"/simulations/simulation2_",file.suffix,".Rda"))



comp.time=(proc.time()-main.start)[[3]]
hours=comp.time%/%(60*60)
minutes=(comp.time%%(60*60))%/%60
seconds=(comp.time%%(60*60))%%60

total.time.message=paste("TOTAL: Computation Time:",hours,"hours,",minutes,"minutes and",seconds,"seconds.")
print(total.time.message)
CON=file(time.outputs.file,"a")
writeLines(total.time.message,CON)
close(CON)
save(sim1,sim2,sim2_nonoise,file=paste0(output.folder,"/simulations/simulations_",file.suffix,".Rda"))
