devtools::install_github("krd5520/DPrct")

####### Testing the multivariate histogram method on the Liberia Data
basepath = rprojroot::find_rstudio_root_file()
main.start=proc.time()

source(paste0(basepath,"/functions/generate_simulated_dataset_function.R"))
#function to generate a geometric sequence a_n=start(ratio^n) for n=0,...,end
geomSeq <- function(start,ratio,end){
  end=end-1
  start*ratio**(0:end)
}

nsims=1
nobs=1000 #number of observations
ntreat=2 #number of treatments
treat.effect=5 #treatment effect
intercept=1 #model intercept

#categorical covariates
ncat=4
cat.probs=rep(list(c(1/2,1/2)),ncat) #for generating the data
cat.coef=geomSeq(3,7/11,ncat) #model coefficient

ncont=4
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
synthdata.budget=list("MassiveEps"=c(5000,0),
                      #"PureDP_EpsHalf"=c(1/2,0),
                      #"ApproxDP_EpsHalf"=c(1/2,1/nobs),
                      "PureDP_Eps1"=c(1,0),
                      #"ApproxDP_Eps1"=c(1,1/nobs),
                      #"PureDP_Eps2half"=c(2.5,0),
                      "PureDP_Eps2"=c(5,0)) #overall budgets epsilon, delta

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
n.iters.fulldp=10000
range.alpha=0.05
############


sim.dfs=generate_simulated_dataset(num.cat=ncat,num.cont=ncont,num.treat=ntreat,nobs=nobs,
                                   mod.coefs=mod.coefs,residual.sd=resid.sd,
                                   cat.probs=cat.probs,cont.funcs=cont.funcs,cont.params=cont.params,
                                   rseed=1)
indic.factor=sapply(colnames(sim.dfs),function(x)is.numeric(sim.dfs[,x]))
cont.vars.all=colnames(sim.dfs)[indic.factor]


reg.formulas.sim1=rep(paste0("y~t1+",paste0("x",seq(1,ncat+ncont),collapse="+")),times=length(synthdata.budget))
names(reg.formulas.sim1)=names(synthdata.budget)
response.vars=rep("y",times=length(reg.formulas.sim1))

synthdata.budget.eps=sapply(synthdata.budget,function(x)x[1])
synthdata.budget.del=sapply(synthdata.budget,function(x)x[2])

mod.names=names(synthdata.budget)

#true.coef.vals=c(Treatment effect 1, intercept)
true.coef.vals=c(treat.effect)
stderr.func=function(mod,std.type="HC1"){
  sqrt(diag(sandwich::vcovHC(mod,type=std.type)))
}

#source(paste0(basepath,"/functions/reg_assumptions_function.R"))
#png(paste0(basepath,"/output/figures/diagnosticFitPlots_Simulation_privbudget.png"),width=1500,height=900)
#reg_assumptions(reg.formulas.sim1[1],data=sim.dfs,pt.sz=5,ln.sz=2)
#dev.off()


#save(out.df,comp.times,file=paste0(basepath,"/inst/tables/ITT_compare_methods_data.Rda"))

source(paste0(basepath,"/functions/sim_compare_method_function.R"))

sim1.out=#parallel::mclapply(
  lapply(
    seq(1,nsims),
    function(idx)
      sim_compare_methods(simdata=sim.dfs,
                          reg.models=reg.formulas.sim1,
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
ITT1=dplyr::bind_rows(lapply(sim1.out,"[[",1)) #combine all repetitions
View(ITT1)
#NOTE: RERUN
#save(sim1.out,sim.dfs,model.params,file=paste0(basepath,"/output/compare_privbudget_simulations_v1206.Rda"))


load(paste0(basepath,"/output/compare_privbudget_simulations_v1206.Rda"))

#get table of Average Treatment Effect, Average Overlap, Average Absolute Diff, Max Std. Err
ITT1=dplyr::bind_rows(lapply(sim1.out,"[[",1)) #combine all repetitions
ITT1$groupid=paste0(ITT1$Model,"_",ITT1$Method,"_",ITT1$Treatment)

ITT1$groupid=paste0(ITT1$Model,"_",ITT1$Method,"_",ITT1$Treatment)
summary.tab=ITT1%>%
  #filter(Model!="MassiveEps")%>%
  group_by(Model,Method,Treatment,groupid)%>%
  summarise("Avg.ITT"=mean(ITT),
            "Avg.CI.overlap"=mean(CI.overlap,na.rm=T),
            "Avg.Abs.Err"=mean(Abs.Err,na.rm=T),
            "Avg.Rel.Err"=mean(Rel.Err,na.rm=T),
            "Max.StdErr"=max(StdErr) )
summary.tab$Epsilon=1
summary.tab$Epsilon[grepl("EpsHalf",summary.tab$groupid)]=0.5
summary.tab$Epsilon[grepl("Eps2",summary.tab$groupid)]=2
summary.tab$Epsilon[grepl("Eps2half",summary.tab$groupid)]=2.5
summary.tab$Epsilon[grepl("MassiveEps",summary.tab$groupid)]=5000
summary.tab$Epsilon[grepl("Confidential",summary.tab$groupid)]=0
summary.tab$Epsilon[grepl("Eps5",summary.tab$groupid)]=5
summary.tab$Delta=ifelse(grepl("Approx",summary.tab$groupid),"1/1000","0")#c(rep("$\\frac{1}{1000}$",12)," "," ",rep("0",24))
summary.tab$Algorithm="DP Model-based"
summary.tab$Algorithm[grepl("Confidential",summary.tab$groupid)]="Confidential"
summary.tab$Algorithm[grepl("MV",summary.tab$groupid)]="MV Histogram"
summary.tab$Algorithm[grepl("Hybrid",summary.tab$groupid)]="Hybrid"

summary.tab.treff=summary.tab[summary.tab$Treatment!="control",]%>%
  #filter(Epsilon!=2000)%>%
  arrange(desc(Algorithm))%>%arrange(desc(Delta))%>%arrange(desc(Epsilon))
#copy and paste the output of the following into latex to get table
# OR use R markdown with latex capabilities and chunk output="asis"
knitr::kable(summary.tab.treff[c(nrow(summary.tab.treff),seq(1,nrow(summary.tab.treff)-1)),c(10,11,12,5,6,8)],
             col.names = c("Epsilon","Delta","Algorithm",
                           "Mean $\\treff$ Estimate","Mean $\\cioverlap$",
                           "Max Std. Err."),
             digits=2,format="latex",booktabs=T,
             linesep=c("\\addlinespace",rep(c("","","\\addlinespace"),7)),
             escape=F)


load(paste0(basepath,"/output/compare_covariate_number_simulations_v1129.Rda"))
ITT2=dplyr::bind_rows(lapply(sim2.out,"[[",1))
ITT2$groupid=paste0(ITT2$Model,"_",ITT2$Method,"_",ITT2$Treatment)
ITT2$NCont=as.numeric(sapply(strsplit(ITT2$Model,"_"),function(x)gsub("Cont","",x[1])))
ITT2$NCat=as.numeric(sapply(strsplit(ITT2$Model,"_"),function(x)gsub("Bin","",x[2])))
ITT2$NCovariates=ITT2$NCat+ITT2$NCont
summary.tab=ITT2%>%group_by(Method,Treatment,Model,NCont,NCat,NCovariates)%>%
  summarise( "Avg.ITT"=mean(ITT),
             "Avg.CI.overlap"=mean(CI.overlap,na.rm=T),
             #"Avg.Abs.Err"=mean(Abs.Err,na.rm=T),
             "Avg.Rel.Err"=mean(Rel.Err,na.rm=T),
             "Max.StdErr"=max(StdErr)
  )
#reorder Method factor levels
summary.tab$Method=factor(summary.tab$Method,levels=c("Confidential","MV Histogram","Hybrid","DP Model-based"))
#summary.tab$ModelLab=factor(summary.tab$Model,levels=c(paste("Model",seq(1,9)),"FiniteSet","Unbounded"))
summary.tab.catcont=summary.tab%>%#filter(!grepl("Finite|Bounded",Model))%>%
  arrange(Method)%>%arrange(NCat)%>%arrange(NCont)
rownames(summary.tab.catcont)=NULL
knitr::kable(summary.tab.catcont[summary.tab.catcont$Treatment!="control",c(4,5,1,7,8,9)],
             col.names = c("Cont.","Cat.","Algorithm","Mean $\\treff$ Estimate","Mean $\\cioverlap$",
                           #"Mean Abs. Diff. ($\\treff=5$)",
                           "Max Std. Err."),
             digits=2,format="latex",booktabs=T,linesep=c("","","","\\addlinespace"),escape=F)

