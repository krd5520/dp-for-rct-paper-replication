#inputs should include synthdata.budget, bins.params, use.continuous.noise,
#   treat.vars,cont.vars.all, mv.prop.budget,win.y.resvar,win.y.trcoef,bound.means,
#   bound.sds,n.iters.fulldp, range.alpha,mod.names, ci.confidence,rseed,treat.effect,
#   stderr.func,se.normal, dpmodel.higher.budget, priv.eps,priv.delta, combos.NCov

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
file.suffix="v0220_test" #added to the files names for all generated files
#suffixes for saving the checkpoint in simulation 1,
#    and saving the check and diagnostic plots in simulation 2
sim1.conf.suffix="conf_sim_budgets"
sim2.conf.suffix="conf_sim_models"
save.plot.png=T

dir.set=directory_setup(basepath,output.folder,
                        list("simulations","figures","tables",
                             paste0("simulations/sims_",file.suffix),paste0("simulations/sims_",file.suffix,"/sim2")))


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
nsims=5 #number of repetitions
nobs=1000 #number of observations
ntreat=2 #number of treatments (including control)
treat.effect=5 #treatment effect
intercept=1 #model intercept
resid.sd=2 #residual standard deviation
se.normal=T #use central limit theory to treat the collection of se values from proxies as normal sample.

#names for treatment and response variables
response.vars="y"
treat.vars="t1"

cat.probs=list(c(1/2,1/2))
#geometric sequence generated coefficients are defined by the start and ratio
#parameters. If start=c and ratio=r, the coefficients are c, c*(r^1), c*(r^2),....
cat.gseq=c(3,7/11)#(start,ratio)
cont.gseq=c(0.99,2/3)#(start,cont)
cont.funcs=list(stats::rnorm,runif.round)
cont.params=list(list("n"=nobs,"mean"=0,"sd"=2),
                 list("n"=nobs,"min"=0,"max"=0.2,"d"=2))


#### Privacy Mechanism Parameters ######
### other global paramaters
bins.param=2/3 #bins for multivariate histograms

use.continuous.noise=T #whether uniform noise should be added to midpoint of continuous variables

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
n.iters.fulldp=6000
range.alpha=0.05

### Analysis Parameters
stderr.func=function(mod,std.type="HC1"){
  sqrt(diag(sandwich::vcovHC(mod,type=std.type)))
}
ci.confidence=0.95
############


########### Compare Regression Model Parameters ###################
#combinations of number of categorical and continuous covariates
sim.data.rseed=1
san.rseed=1
combos.NCov=data.frame("NCat"=rep(c(4,10,20),each=3),
                       "NCont"=rep(c(2,5,10),3),
                       "Name"=paste("Model",seq(1,9)))
priv.eps=1
priv.delta=0
dpmodel.higher.budget=TRUE

### Define Additional Regression Models
mx.ncont=max(combos.NCov$NCont)
add.ncat=4
add.ncont=mx.ncont%/%2
cat.formula=paste0("x",seq(mx.ncont+1,mx.ncont+add.ncat),collapse="+")
base.formula=paste0(response.vars,"~",paste0(treat.vars,collapse="+"))
#additional formulas one with only the unbounded continuous and another with only the finite continuous
add.formulas=c(
  paste0(base.formula,"+",paste0("x",seq(1,mx.ncont-1,2),collapse="+"),"+",cat.formula),
  paste0(base.formula,"+",paste0("x",seq(2,mx.ncont,2),collapse="+"),"+",cat.formula))
names(add.formulas)=c("Unbounded","Finite Set")

ntreat=length(treat.effect)
reg.formulas.sim2=sim2_models_list(response=response,treat.vars=treat.vars,combos.NCov=combos.NCov,
                                   add.formulas=add.formulas,add.names=add.names)
mod.names=names(reg.formulas.sim2)

if(is.null(add.formulas)==FALSE){
  add.NCov=data.frame("NCat"=add.ncat,"NCont"=add.ncont,"Name"=add.names)
  combos.NCov=dplyr::bind_rows(combos.NCov,add.NCov)
}
colnames(combos.NCov)=c("NCat","NCont","Model")


response="y"
treat.vars="t1",

dpmodel.higher.budget=TRUE,

conf.plot=list("pt.sz1"=1,"ln.sz1"=1.5,"bs.sz1"=12,
               "pt.sz2"=1,"ln.sz2"=1.5,"bs.sz2"=15,
               "ncol1"=3,"ncol2"=2,
               "width1"=4,"height1"=2.5,
               "width2"=6,"height1"=3.5),

se.normal=T

ncat=max(combos.NCov$NCat)
ncont=max(combos.NCov$NCont)
mod.coefs=sim_coefs(ncat=ncat,ncat.groups=sapply(cat.probs,length),ncont=ncont,
                    cat.coef.start=cat.gseq[1],cat.coef.ratio=cat.gseq[2],
                    cont.coef.start=cont.gseq[1],cont.coef.ratio=cont.gseq[2],
                    treat.effect=treat.effect,intercept=intercept)
model.params=list("model.coefs"=mod.coefs,"residual.err.sd"=resid.sd)
cont.vars.all=c("y",paste0("x",seq(1,ncont)))
synthdata.budget=rep(list("Budget1"=c(priv.eps,priv.delta)),length(reg.formulas.sim2)) #overall budgets epsilon, delta


sim.dfs=generate_simulated_dataset(num.cat=ncat,num.cont=ncont,num.treat=ntreat,nobs=nobs,
                                   mod.coefs=mod.coefs,residual.sd=resid.sd,
                                   cat.probs=cat.probs,cont.funcs=cont.funcs,
                                   cont.params=cont.params,
                                   rseed=sim.data.rseed)




conf.summaries=lapply(seq(1,length(reg.formulas.sim2)),
                      function(idx)summary(glm(reg.formulas.sim2[idx],data=sim.dfs))$coefficients)

conf.plots= lapply(seq(1,length(reg.formulas.sim2)-length(add.formulas)),
                   function(idx)reg_assumptions(reg.formulas.sim2[idx],
                                                data=sim.dfs,
                                                response.var=response,
                                                pt.sz=conf.plot$pt.sz1,
                                                ln.sz=conf.plot$ln.sz1,
                                                bs.sz=conf.plot$bs.sz1,
                                                mod.name=mod.names[idx]))
save.plot.file=paste0(output.folder,"/figures/simulations_diagnostic_",conf.suffix,"_",file.suffix,".png")
ggplot2::ggsave(filename=save.plot.file,bg="white",
                width=conf.plot$width1,height=conf.plot$height1,
                plot=cowplot::plot_grid(plotlist=conf.plots,ncol=conf.plot$ncol1))
#dev.off()
if(is.null(add.formulas)==FALSE){
  add.plots=lapply(seq(1+length(reg.formulas.sim2)-length(add.formulas),length(reg.formulas.sim2)),
                   function(idx)reg_assumptions(reg.formulas.sim2[idx],
                                                data=sim.dfs,
                                                response.var=response,
                                                pt.sz=conf.plot$pt.sz2,
                                                ln.sz=conf.plot$ln.sz2,
                                                bs.sz=conf.plot$bs.sz2,
                                                mod.name=mod.names[idx]))
  save.plot.name.xtra=paste0(output.folder,"/figures/simulations_diagnostic_",conf.suffix,"_",file.suffix,"_xtra.png")
  ggplot2::ggsave(filename=save.plot.name.xtra,bg="white"
                  ,width=conf.plot$width2,height=conf.plot$height2,
                  plot=cowplot::plot_grid(plotlist=add.plots,ncol=conf.plot$ncol2))
  #dev.off()
}else{
  add.plots=NULL
}

save(sim.dfs,
     reg.formulas.sim2,synthdata.budget, bins.params, use.continuous.noise,
     treat.vars,cont.vars.all, mv.prop.budget,win.y.resvar,win.y.trcoef,bound.means,
     bound.sds,n.iters.fulldp, range.alpha,mod.names, ci.confidence,rseed,treat.effect,
     stderr.func,se.normal, dpmodel.higher.budget, priv.eps,priv.delta, combos.NCov,
     file=paste0(output.folder,"/simulations/sims_",file.suffix,"/sim2_inputs_",conf.suffix,"_",file.suffix,".Rda"))




