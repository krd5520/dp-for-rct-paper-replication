devtools::install_github("krd5520/DPrct",force=T)

####### Testing the multivariate histogram method on the Liberia Data
basepath = rprojroot::find_rstudio_root_file()
main.start=proc.time()

source(paste0(basepath,"/functions/generate_simulated_dataset_function.R"))
#function to generate a geometric sequence a_n=start(ratio^n) for n=0,...,end
geomSeq <- function(start,ratio,end){
  end=end-1
  start*ratio**(0:end)
}

#########
#### Simulation Data Parameters #####
nsims=20
nobs=1000 #number of observations
ntreat=2 #number of treatments (including control)
treat.effect=5 #treatment effect
intercept=1 #model intercept

#combinations of number of categorical and continuous covariates
combos.num.covariate=data.frame("NCat"=c(rep(4,3),rep(10,3),rep(20,3)),
                               #"NCont"=rep(c(4,10,20),3),
                               "NCont"=rep(c(2,5,10),3))


#categorical covariates
ncat=max(combos.num.covariate$NCat)
cat.probs=rep(list(c(1/2,1/2)),ncat) #for generating the data
cat.coef=geomSeq(3,7/11,length(unlist(cat.probs))-ncat) #model coefficient

ncont=max(combos.num.covariate$NCont)
#function to generate and round n observations from a uniform distribution to d number of digits
runif.round=function(n,min,max,d=2){
  round(runif(n=n,min=min,max=max),digits=d)
}
#half (rounded up) of the continuous covariates are Normal(0,1),
# the other half are Uniform(0,1) rounded to 2 digits
#
half.ncont=(ncont%/%2)
cont.funcs=rep(list(stats::rnorm,runif.round),
               times=c(half.ncont+(1*(ncont%%2)),half.ncont))
cont.params=rep(list(list("n"=nobs,"mean"=0,"sd"=2),
                     list("n"=nobs,"min"=0,"max"=0.2,"d"=2)),
                times=c(half.ncont+(1*(ncont%%2)),half.ncont))
cont.coef=geomSeq(0.99,2/3,ncont) #mod.coefs

#model mod.coefs and residual error standard deviation
mod.coefs=c(intercept,treat.effect,cont.coef,cat.coef)
resid.sd=2

########## Privacy Mechanism Parameters #############
### other global paramaters
bins.param=2/3 #bins for multivariate histograms

use.continuous.noise=TRUE #whether uniform noise should be added to midpoint of continuous variables
ci.confidence=0.95

#### parameters for Full DP Model- Based
# privacy budget variables #
mv.prop.budget=0.2
## within each response, the residual variance gets 10% of the budget
win.y.resvar=0.1
##      the treatment mod.coefs get 30% of the budget (i.e. 10% for each)
win.y.trcoef=0.3
##
# other params #
bound.means=100
bound.sds=c(2^(-15),2^(15))
n.iters.fulldp=5000
range.alpha=0.05
############





base.formula="y~t1"

reg.formulas.sim2=sapply(seq(1,nrow(combos.num.covariate)),
                                function(idx)paste0(base.formula,"+",paste0("x",seq(1,combos.num.covariate$NCont[idx]),collapse="+"),
                                                    "+",
                                                    paste0("x",seq(ncont+1,ncont+combos.num.covariate$NCat[idx]),collapse="+")))
names(reg.formulas.sim2)=paste0("Model",seq(1,length(combos.num.covariate$NCont)))#"Cont",combos.num.covariate$NCont,"_Bin",combos.num.covariate$NCat)
reg.formulas.allonetypecont=c(paste0(base.formula,"+",
                                     paste0("x",seq(1,ncont,2),collapse="+"),
                                     "+",paste0("x",seq(ncont+1,ncont+4),collapse="+")),
                              paste0(base.formula,"+",
                                     paste0("x",seq(2,ncont,2),collapse="+"),
                                     "+",paste0("x",seq(ncont+1,ncont+4),collapse="+")))
names(reg.formulas.allonetypecont)=c("UnBounded_Cont10_Bin4","FiniteSet_Cont10_Bin4")
reg.formulas.sim2=c(reg.formulas.sim2,reg.formulas.allonetypecont)
#reg.formulas.sim2=#paste(paste0("y",seq(1,length(reg.formulas.sim2))),reg.formulas.sim2)
response.vars="y"#paste0("y",seq(1,length(reg.formulas.sim2)))
mod.names=c(paste("Model",seq(1,length(combos.num.covariate$NCont))),"Unbounded","FiniteSet")

sim.dfs=generate_simulated_dataset(num.cat=ncat,num.cont=ncont,num.treat=ntreat,nobs=nobs,
                                   mod.coefs=mod.coefs,residual.sd=resid.sd,
                                   cat.probs=cat.probs,cont.funcs=cont.funcs,
                                   cont.params=cont.params,#reg.mods = reg.formulas.sim2,
                                   rseed=1)

indic.factor=sapply(colnames(sim.dfs),function(x)is.numeric(sim.dfs[,x]))
cont.vars.all=colnames(sim.dfs)[indic.factor]
synthdata.budget=rep(list("PureDP_Eps1"=c(1/2,0)),length(reg.formulas.sim2)) #overall budgets epsilon, delta



#true.coef.vals=c(Treatment effect 1, intercept)
true.coef.vals=c(treat.effect)
stderr.func=function(mod,std.type="HC1"){
  sqrt(diag(sandwich::vcovHC(mod,type=std.type)))
}

conf.summaries=lapply(seq(1,length(reg.formulas.sim2)),function(idx)summary(glm(reg.formulas.sim2[idx],data=sim.dfs))$coefficients)

source(paste0(basepath,"/functions/reg_assumptions_function.R"))
tempplot= lapply(seq(1,length(reg.formulas.sim2)-2),
                function(idx)reg_assumptions(reg.formulas.sim2[idx],
                                             data=sim.dfs,
                                             response.var=response.vars[idx],
                                             pt.sz=4,
                                             ln.sz=2.5,
                                             bs.sz=33,
                                             mod.name=mod.names[idx]))
cowplot::plot_grid(plotlist=tempplot,ncol=2)
png(paste0(basepath,"/output/figures/diagnosticFitPlots_simulations_numCovariates.png"),width=4000,height=2500)
print(cowplot::plot_grid(plotlist=tempplot,ncol=3))
dev.off()
#save(out.df,comp.times,file=paste0(basepath,"/inst/tables/ITT_compare_methods_data.Rda"))

source(paste0(basepath,"/functions/sim_compare_method_function.R"))

sim2.out=parallel::mclapply(seq(1,nsims),function(idx)sim_compare_methods(simdata=sim.dfs,
                                                             reg.models=reg.formulas.sim2,
                                                             synthdata.budget=synthdata.budget,
                                                             bins.param=bins.param,
                                                             use.continuous.noise=use.continuous.noise,
                                                             treat.vars="t1",
                                                             cont.as.cont=cont.vars.all,
                                                             mv.prop.budget=mv.prop.budget,
                                                             win.y.resvar=win.y.resvar,
                                                             win.y.trcoef=win.y.trcoef,
                                                             bound.means=bound.means,
                                                             bound.sds=bound.sds,
                                                             n.iters.fulldp=n.iters.fulldp,
                                                             range.alpha=range.alpha,
                                                             mod.names=mod.names,
                                                             ci.confidence=ci.confidence,
                                                             rseed=idx,
                                                             sim=paste0("Sim",idx),
                                                             include.full.model.based=TRUE,
                                                             true.coef.vals=true.coef.vals,
                                                             se.func=stderr.func))
model.params=list("model.coefs"=mod.coefs,"residual.err.sd"=resid.sd)

#NOTE: RERUN.
save(sim2.out,sim.dfs,model.params, conf.summaries,
     file=paste0(basepath,"/output/compare_covariate_number_simulations_v1205.Rda"))


