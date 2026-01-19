devtools::install_github("krd5520/DPrct")

#renv::install("krd5520/DPrct",verbose=F,lock=T)
main.start=proc.time()

####### Testing the multivariate histogram method on the Liberia Data
#file path and file location parameters
basepath = rprojroot::find_rstudio_root_file()
file.suffix="v0212"
outputpath=paste0(basepath,"/output")
in.datapath=paste0(basepath,"/data/liberia_subset_nolabels.Rda")
out.datapath=paste0(basepath,"/data")
override.checkpoints=T

source(paste0(basepath,"/program/functions/generate_simulated_dataset_function.R"))
source(paste0(basepath,"/program/functions/simple_simulation_functions.R"))
source(paste0(basepath,"/program/functions/reg_assumptions_function.R"))
source(paste0(basepath,"/program/functions/sim_compare_method_function.R"))
source(paste0(basepath,"/program/functions/main_simulation_modelcompare.R"))
source(paste0(basepath,"/program/functions/main_simulation_budgetcompare.R"))
source(paste0(basepath,"/program/functions/directory_setup.R"))
source(paste0(basepath,"/program/functions/main_simulation_nonoise.R"))
source(paste0(basepath,"/program/functions/sim_nonoise_dpmb.R"))

alt_multivariate_histogram<-function(data,continuous.vars=NULL,
                                 num.bin=NULL,bin.param=NA,
                                 which.cont.out=FALSE,return.time=TRUE){
  ### Check Inputs ###
  stopifnot(base::is.data.frame(data)) #check data input

  ### Multivariate Number of Bins Vector ####
  ## each entry is the number of bins for the corresponding column in data

  #get number of continuous variables
  if(base::is.null(continuous.vars)){ #if continuous.vars=NULL
    num.continuous=0 #number continuous variables=0, no need num.bin or bin.param
    #mx.n.combos=prod(sapply(colnames(data),function(x)length(table(data[,x],useNA = "ifany"))))
    #return(base::data.frame(base::table(data)))
  }else{ #if continuous.vars is supplied
    cont.data=data[,continuous.vars,drop=F] #subset data to be continuous variables
    if(base::ncol(cont.data)==0){
      warning(paste0("No continuous variables detected. Dimension of cont.data is...",dim(cont.data),"... and length of continuous.vars is...",length(continuous.vars)))

    }

    #number of continuous variables
    num.continuous=ncol(cont.data)
  }

  if(num.continuous>0){ #if there are continuous variables
    if(base::sum(base::sapply(cont.data,base::is.numeric))!=num.continuous){ #check columns are numeric
      stop(
        base::paste("Continuous variables selected are not numeric values. Columns Selected:",
                    base::paste(base::colnames(cont.data),collapse=", ")))
    } #end if continuous variables are not numeric

    #if there are continuous variables, either num.bin or bin.param must be numeric
    stopifnot(base::is.numeric(num.bin)|base::is.numeric(bin.param))
    if(base::is.numeric(num.bin)==T){ #num.bin must be a single value or a value for each column
      stopifnot(length(num.bin)==1|base::length(num.bin)==num.continuous)

      if(base::length(num.bin)==1){ #if num.bin is single value
        num.bin=base::rep(num.bin,num.continuous)
      } #end if num.bin single value

    }else{ #if num.bin not supplied bin.param must be a single value between 0 and 1
      stopifnot(bin.param>0&bin.param<1)
      num.bin=base::rep(base::ceiling((base::nrow(cont.data)^bin.param)),num.continuous)
    }
    cont.gr.bin=sapply(seq(1,num.continuous),function(i)length(unique(cont.data[,i]))>num.bin[i])
    if(sum(!cont.gr.bin)>0){
      print("in sum>0")
      cont.data=cont.data[,cont.gr.bin,drop=F]
      num.continuous=ncol(cont.data)
    }
    which.cont=base::colnames(data)%in%base::colnames(cont.data) #logical if continuous variable
    print(sum(which.cont))
    temp=base::sapply(base::seq(1,num.continuous),
                      function(i)DPrct::continuous_bins(cont.data[,i],num.bin[i]))
    print(dim(temp))
    print(nrow(data))
    print(dim(cont.data))
    data[,which.cont]=base::sapply(base::seq(1,num.continuous),
                                   function(i)DPrct::continuous_bins(cont.data[,i],num.bin[i]))
  }# end if there are continuous variables
  mv.histogram=dplyr::count(data,data[seq(1,nrow(data)),],name="Freq")

  #if which.cont.out=TRUE return mv.histogram and which.cont,else only mv.histogram
  if(which.cont.out==TRUE){
    return(base::list(mv.histogram,which.cont))
  }else{
    return(mv.histogram)
  }
}
nsims=20 #number of repetitions
nobs=1000 #number of observations
ntreat=2 #number of treatments (including control)
treat.effect=5 #treatment effect
intercept=1 #model intercept
resid.sd=2 #residual standard deviation

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

use.continuous.noise=TRUE #whether uniform noise should be added to midpoint of continuous variables

#### parameters for Full DP Model- Based
# privacy budget variables #
mv.prop.budget=0.5
## within each response, the residual variance gets 10% of the budget
win.y.resvar=0.1
##      the treatment mod.coefs get 30% of the budget (i.e. 10% for each)
win.y.trcoef=0.1
##
# other params #
bound.means=100
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
sim1.simdata.rseed=1
sim1.san.rseed=1
#overall budgets epsilon, delta
sim1.synthdata.budget=list("MassiveEps"=c(5000,0),
                           "PureDP_EpsHalf"=c(1/2,0),
                           "ApproxDP_EpsHalf"=c(1/2,1/nobs),
                           "PureDP_Eps1"=c(1,0),
                           "ApproxDP_Eps1"=c(1,1/nobs),
                           "PureDP_Eps2"=c(2,0),
                           "PureDP_Eps4"=c(4,0))
sim1.ncat=4
sim1.ncont=4



########### Compare Regression Model Parameters ###################
#combinations of number of categorical and continuous covariates
sim2.simdata.rseed=1
sim2.san.rseed=1
combos.NCov=data.frame("NCat"=rep(c(4,10,20),each=3),
                       "NCont"=rep(c(2,5,10),3),
                       "Name"=paste("Model",seq(1,9)))
sim2.eps=1
sim2.delta=0
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



treat.effect=5
ntreat=length(treat.effect)
reg.formulas.sim2=sim2_models_list(response=response.vars,treat.vars=treat.vars,combos.NCov=combos.NCov,
                                   add.formulas=add.formulas,add.names=names(add.formulas))
mod.names=names(reg.formulas.sim2)

ncat=max(combos.NCov$NCat)
ncont=max(combos.NCov$NCont)
mod.coefs=sim_coefs(ncat=ncat,ncat.groups=sapply(cat.probs,length),ncont=ncont,
                    cat.coef.start=cat.gseq[1],cat.coef.ratio=cat.gseq[2],
                    cont.coef.start=cont.gseq[1],cont.coef.ratio=cont.gseq[2],
                    treat.effect=treat.effect,intercept=intercept)
model.params=list("model.coefs"=mod.coefs,"residual.err.sd"=resid.sd)
cont.vars.all=c("y",paste0("x",seq(1,ncont)))
#synthdata.budget=rep(list("Budget1"=c(priv.eps,priv.delta)),length(reg.formulas.sim2)) #overall budgets epsilon, delta

treat.effect=5
sim.dfs=generate_simulated_dataset(num.cat=ncat,num.cont=ncont,num.treat=ntreat,nobs=nobs,
                                   mod.coefs=mod.coefs,residual.sd=resid.sd,
                                   cat.probs=cat.probs,cont.funcs=cont.funcs,
                                   cont.params=cont.params,
                                   rseed=1)
tempmv=alt_multivariate_histogram(sim.dfs,c("y",paste0("x",seq(1,10,2))),num.bin = 667)
temp2=DPrct::synthdata_perturb_mvhist(data=sim.dfs[,colnames(sim.dfs)!="y"],1,0,continuous.vars=paste0("x",seq(1,10,2)),
                                     num.bin=NULL,bin.param=2/3,
                                     add.cont.variation = F,
                                     treatment.colname="t1",with.treatment = F,within.blocks = F,
                                     which.cont.out=FALSE,return.time=FALSE,rseed=1)


head(temp2)

conf.model=glm(reg.formulas.sim2[1],data=sim.dfs)
source(paste0(basepath,"/program/functions/generate_simulated_dataset_function.R"))
source(paste0(basepath,"/program/functions/simple_simulation_functions.R"))
source(paste0(basepath,"/program/functions/reg_assumptions_function.R"))
source(paste0(basepath,"/program/functions/sim_compare_method_function.R"))
source(paste0(basepath,"/program/functions/main_simulation_modelcompare.R"))
source(paste0(basepath,"/program/functions/main_simulation_budgetcompare.R"))
source(paste0(basepath,"/program/functions/directory_setup.R"))
source(paste0(basepath,"/program/functions/main_simulation_nonoise.R"))
source(paste0(basepath,"/program/functions/sim_nonoise_dpmb.R"))
temp=nonoise_dpmb(formula=reg.formulas.sim2[9],
                  confidential.data=sim.dfs,
                  synth.data=temp2,
                  continuous.vars=paste0("x",seq(1,10,2)),
                  num.bin=NULL,
                  bin.param=NA,
                  num.iters=5000,
                  ci.alpha=0.05,
                  assign.type="simple",
                  blocks=NULL,
                  within.blocks=FALSE,
                  treatment.var=NULL,
                  rseed=NA,
                  clusters=NULL,
                  return.time=FALSE,
                  return.confidential.table=FALSE,
                  return.san.summary=TRUE,
                  use.san.residerror=TRUE)
head(temp)




sim_compare_methods(simdata=sim.dfs,
                    reg.models=reg.formulas.sim2,
                    synthdata.budget=synthdata.budget.double,
                    bins.param=bins.param,
                    use.continuous.noise=use.continuous.noise,
                    treat.vars=treat.vars,
                    cont.as.cont=cont.vars.all,
                    mv.prop.budget=mv.prop.budget,
                    win.y.resvar=win.y.resvar,
                    win.y.trcoef=win.y.trcoef,
                    bound.means=bound.means,
                    bound.sds=bound.sds,
                    n.iters.fulldp=n.iters.fulldp,
                    range.alpha=range.alpha,
                    mod.names=mod.names.double,
                    ci.confidence=ci.confidence,
                    rseed=san.rseed+idx,
                    sim=paste0("Sim",idx),
                    include.full.model.based=TRUE,
                    only.full.model.based=TRUE,
                    true.coef.vals=c(treat.effect),
                    se.func=stderr.func,
                    se.normal=se.normal))

sim_compare_methods(simdata=sim.dfs,
                    reg.models=reg.formulas.sim1,
                    synthdata.budget=synthdata.budget,
                    bins.param=bins.param,
                    use.continuous.noise=use.continuous.noise,
                    treat.vars=treat.vars,
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
                    rseed=san.rseed+idx,
                    sim=paste0("Sim",idx),
                    include.full.model.based=TRUE,
                    true.coef.vals=c(treat.effect),
                    se.func=stderr.func,se.normal=se.normal)


#renv::install("krd5520/DPrct",verbose=F,lock=T)
main.start=proc.time()

####### Testing the multivariate histogram method on the Liberia Data
#file path and file location parameters
basepath = rprojroot::find_rstudio_root_file()
file.suffix="v0212_test"
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
bins.param=2/3 #bins for multivariate histograms
synthdata.budget=c(1,0) #overall budget epsilon, delta
use.continuous.noise=FALSE #whether uniform noise should be added to midpoint of continuous variables

#### parameters for DP Model- Based
# privacy budget variables #
mv.prop.budget=0.5
## within each response, the residual variance gets 20% of the budget
win.y.resvar=0.1
##      the treatment coefficients get 30% of the budget (i.e. 10% for each)
win.y.trcoef=0.3
##
# other params #
bound.means=75
bound.sds=c(2^(-15),2^(15))
n.iters=5000
range.alpha=0.05 #for range function
ci.confidence=0.95

#regression diagnostic plot properties
reg.assumption.plots=T
reg.plot.props=list("pt.sz"=0.3,"ln.sz"=0.7,"bs.sz"=9,"ncol"=2,"width"=8,"height"=7)


#### Other Budgets for Full DP
### otherbudgetsDP list with each element representing a budget. Each element should have mv.epsilon,mv.delta,per.y.eps,per.y.del
otherbudgetsDP=list("DP-Mb Fix MV Budget Double Overall"=c(synthdata.budget,synthdata.budget/8),
                    "DP-Mb Fix MV Budget to Each Y Budget"=c(synthdata.budget,synthdata.budget),
                    "DP-Mb Fix Overall Budget"=c(synthdata.budget/2,synthdata.budget/16))

########

each.rseed=2
############ END PARAMETERS ##################
source(paste0(basepath,"/program/functions/directory_setup.R"))
source(paste0(basepath,"/program/functions/main_liberia_compare_methods.R"))
source(paste0(basepath,"/program/functions/liberia_compare_methods_onecovset.R"))
source(paste0(basepath,"/program/functions/reg_assumptions_function.R"))
dir.set=directory_setup(basepath,outputpath,list("liberia/liberia_checkpoints",
                                                 "figures",paste0("figures/liberia_",file.suffix)))



red.sub.covariates=c('age_b', 'asbhostil_b', 'drugssellever_b', 'drinkboozeself_b',
                     'druggrassself_b', 'harddrugsever_b', 'steals_b')

sub.covariates=c(red.sub.covariates, 'wealth_indexstd_b',
                 'school_b', 'cognitive_score_b')

baseline.full=c('age_b', 'livepartner_b', 'mpartners_b', 'hhunder15_b', 'famseeoften_b',
                'muslim_b', 'school_b', 'schoolbasin_b', 'literacy_b', 'mathscore_b', 'health_resc_b',
                'disabled_b', 'depression_b', 'distress_b', 'rel_commanders_b', 'faction_b', 'warexper_b',
                'profitsump99avg7d_b', 'wealth_indexstd_b', 'homeless_b', 'slphungry7dx_b', 'savstockp99_b',
                'loan50_b', 'loan300_b', 'illicit7da_zero_b', 'agricul7da_zero_b', 'nonagwage7da_zero_b',
                'allbiz7da_zero_b', 'nonaghigh7da_zero_b', 'agriculeveramt_b', 'nonagbizeveramt_b',
                'nonaghigheveramt_b', 'drugssellever_b', 'drinkboozeself_b', 'druggrassself_b',
                'grassdailyuser_b', 'harddrugsever_b', 'harddrugsdailyuser_b', 'steals_b',
                'stealnb_nonviol_b', 'stealnb_felony_b', 'disputes_all_b', 'asbhostil_b',
                'conscientious_b', 'neurotic_b', 'grit_b', 'rewardresp_b', 'locuscontr_b', 'impulsive_b',
                'selfesteem_b', 'patient_game_real_b', 'inconsistent_game_resc_b', 'risk_game_resc_b',
                'timedecl_b', 'riskdecl_b', 'cognitive_score_b', 'ef_score_b')
########

#list all possible response variables for round==5
response.vars.all = c('fam_asb_lt',
                      paste0(c('drugssellever','stealnb','disputes_all_z','carryweapon',
                               'arrested','asbhostilstd','domabuse_z'),"_ltav"))
treat.vars=c("cashassonly","tpassonly","tpcashass") #treatment variables
block.vars=c("tp_strata_alt","cg_strata") #block variables
mod.names=response.vars.all
# mod.names=c("Anti-social Behaviors, z-score",
#             "Usually Sells Drugs",
#             "# of thefts/robberies in past 2 weeks",
#             "Disputes and fights in past 2 weeks, z-score",
#             "Carries a weapon on body",
#             "Arrested in past 2 weeks",
#             "Aggressive behaviors, z-score",
#             "Verbal/physical abuse of partner, z-score")
treat.names=c("Therapy Only","Cash Only","Both")


## read in data
load(in.datapath)
col.class=attr(liberia.sub,"column.classifications")
cat.vars=col.class$cat.vars
cont.as.cont=col.class$cont.as.cont
cont.as.cat=col.class$cont.as.cat
#list2env(attr(liberia.sub,"column.classifications"),envir=enviroment())#globalenv())
cont.vars.all=c(cont.as.cont,cont.as.cat,response.vars.all)
cont.vars.all=cont.vars.all[!duplicated(cont.vars.all)]

if(use.allobs==TRUE){
  liberia.sub=liberia.sub.allobs
}else if(use.obs1==TRUE){
  liberia.sub=liberia.sub.obs1
}
liberia.sub[,colnames(liberia.sub)!="treatment"]=apply(liberia.sub[,colnames(liberia.sub)!="treatment"],2,
                                                       function(x)as.numeric(as.character(x)))
liberia.sub[,block.vars]=apply(liberia.sub[,block.vars],2,
                               function(x)as.factor(as.character(x)))

#########
mv.bins=nrow(liberia.sub)^bins.param

#save(out.df,comp.times,file=paste0(basepath,"/inst/tables/ITT_compare_methods_data.Rda"))

print("Starting Covariate Set Comparisons")
temp=liberia_compare_methods_onecovset(liberia.sub=liberia.sub,
                                             synthdata.budget=synthdata.budget,
                                             mv.bins=mv.bins,use.continuous.noise=use.continuous.noise,
                                             covariate.vars=red.sub.covariates,
                                             block.vars=block.vars,
                                       treat.vars=treat.vars,
                                       response.vars.all=response.vars.all,
                                             cont.as.cont=cont.as.cont[cont.as.cont %in% red.sub.covariates],
                                             mv.prop.budget=mv.prop.budget,win.y.resvar=win.y.resvar,win.y.trcoef=win.y.trcoef,
                                             bound.means=bound.means,bound.sds=bound.sds,n.iters=n.iters,range.alpha=range.alpha,
                                             mod.names=mod.names,treat.names=treat.names,
                                             mult.test.correct=c("treatments","responses"),bonferroni.npvals=1,
                                             ci.confidence=ci.confidence,
                                             rseed=each.rseed,
                                             sim="Reduced",
                                             include.dp.model=TRUE,
                                             otherbudgetsDP=NULL,fname=red.file,
                                             reg.assumption.plots=reg.assumption.plots,
                                             reg.plot.path=paste0(outputpath,"/figures/liberia_",file.suffix),
                                             reg.plot.prefix="liberia_diagnostics_full",
                                             reg.plot.props=reg.plot.props,se.normal=T)

temp=liberia_compare_methods_onecovset(liberia.sub,
                                           synthdata.budget=c(1,0),mv.bins=2/3,
                                           use.continuous.noise=T,
                                           covariate.vars,block.vars,treat.vars,response.vars.all,
                                           cont.as.cont,
                                           mv.prop.budget=NULL,win.y.resvar=NULL,win.y.trcoef=NULL,
                                           bound.means=NULL,bound.sds=NULL,n.iters=NA,range.alpha=0.05,
                                           mod.names,treat.names,
                                           mult.test.correct=c("treatments","responses"),bonferroni.npvals=NA,
                                           ci.confidence=0.95,
                                           rseed=NA,
                                           sim=NA,
                                           include.dp.model=TRUE,
                                           otherbudgetsDP=NULL,otherbudgetsDP.names=NULL,fname="liberia.Rda",
                                           reg.assumption.plots=F,
                                           reg.plot.path=NULL,
                                           reg.plot.prefix="liberia_diagnostics",
                                           reg.plot.props=NULL,
                                           se.normal=T)
