
## Function to compare simulated datasets.
## TO DO: Document inputs and outputs
## INPUTS:
### simdata: simulated data frame with all necessary variables for the reg.models inputted
### reg.models: vector of strings of the regression models to fit
### synthdata.budget: List of vectors of (epsilon, delta) values for each reg.models entry
### bins.params: bin parameter for discritizing continuous variables
### use.continuous.noise: logical for if uniform noise should be added to the bin medians after sampling perturbed histogram
### treat.vars: string with treatment variable column name
### cont.as.cont: vector of column names for continuous values to be discretized
##### INPUTS for DP Model-based Algorithm (see dp_synthdata)
### mv.prop.budget: proportion of privacy budget to be allocated to MV Histogram step
### win.y.resvar, win.y.trcoef: proportion of privacy budget per response variable to be
###             allocated to the residual variance and treatment coefficient respectively.
### bound.means, bound.sds, range.alpha: parameters for dp_confidence_intervals()
### n.iters.fulldp: number of proxy response variables/models fit to generate a sample of coefficients and residuals
### mod.names: vector of strings for the model names corresponding to each reg.models
### treat.names: vector of strings (or just a string) with the treatment name
### bonferroni.npvals: number of comparisons to adjust the p-values with bonferroni correction
### ci.confidence: confidence level for CI of treatment effect from sanitized/confident datasets
### rseed: if not NA, then it sets the random seed at beginning of the function application
### sim: to be added as a column to the output (useful for when combining several column outputs from lapply or mclappy function)
### include.full.model.based: logical for if the DP Model-based should be included in comparisons
### true.coef.vals: value to compare estimated treatment effects to in absolute and relative difference metrics
### se.func: function to get the standard errors for treatment effects from sanitized/confidential datasets
### n.obs: number of observations starting at observation 1 to be used from the simulated dataset
##
##OUTPUTS:
### treff.table: table of with rows for each model, treatment level, and Method combination and columns:
###             (1) ITT estimates and control means, (2) Std. Errors, (3) (Adjusted) Pvalues,
###             (4) Treatment name, (5) Model name (6) CI overlap, (7) Absolute Difference, (8) Method
### comp.times.summary: table of computation times for each method
### comp.times.tab: More detailed table of computation times
### list.of.synthdatas: list of sanitized datasets for each method
### simdata: original simulated dataset from input
sim_compare_methods=function(simdata,
                             reg.models,
                             synthdata.budget,bins.param,use.continuous.noise,
                             treat.vars,
                             cont.as.cont,
                             mv.prop.budget=NA,win.y.resvar=NA,win.y.trcoef=NA,
                             bound.means=NULL,bound.sds=NULL,n.iters.fulldp=NA,range.alpha=0.05,
                             mod.names,treat.names,
                             bonferroni.npvals=NA,
                             ci.confidence=0.95,
                             rseed=NA,
                             sim=NA,
                             include.full.model.based=TRUE,
                             only.full.model.based=FALSE,
                             true.coef.vals=NULL,
                             se.func=NULL,
                             n.obs=NULL,#se.normal=T,
                             num.cores=1,
                             #n.std.dev=15,continuous.limits=NULL,standardize.col=NULL,
                             diagnostic.file=NULL){

  start.time=proc.time() #start time
  if(is.na(rseed)==FALSE){ #set seed if needed
    set.seed(rseed)
  }
  if(include.full.model.based==TRUE){ #check inputs for full.model.based
    stopifnot((is.na(mv.prop.budget)==FALSE)&&(mv.prop.budget>0 & mv.prop.budget<1))
    stopifnot((is.na(win.y.resvar)==FALSE)&(is.na(win.y.trcoef)==FALSE))
    stopifnot((win.y.resvar>0)&(win.y.trcoef>0)&((win.y.resvar+win.y.trcoef)<1))
    stopifnot(is.null(bound.means)==FALSE&is.null(bound.sds)==FALSE&is.na(n.iters.fulldp)==FALSE)

  }

  all.mod.vars=unique(unlist(lapply(reg.models,function(x)all.vars(as.formula(x)))))
  #subset data
  simdata=simdata[,colnames(simdata)%in%c(all.mod.vars,"treatment","control"),drop=F]
  simdata[,treat.vars]=apply(simdata[,treat.vars,drop=F],2,function(x)as.numeric(as.character(x)))
  #extract epsilons and deltas
  synthdata.budget.eps=sapply(synthdata.budget,function(x)x[1])
  synthdata.budget.del=sapply(synthdata.budget,function(x)x[2])
  se.normal=TRUE

  ### Constructing regression model formulas for each model
  n.regs=length(reg.models)


  fit.conf.models.start=proc.time() #start time for confidential data
  comp.times.summary=list("Set Up Time"=(fit.conf.models.start-start.time)[[3]])

  #function to get treatment effects, confidence intervals, absolute, error, and relative error
  get_conf_treff<-function(idx,data=simdata,
                           formulas=reg.models,
                           tr.vars=treat.vars,
                           mod.nms=mod.names,
                           ci.alpha=1-ci.confidence,
                           true.vals=true.coef.vals,
                           stderr.func=se.func,
                           nobs=n.obs){
    start.time.internal=proc.time()
    if(is.null(nobs)==FALSE){
      data=data[seq(1,nobs[idx]),]
    }
    mod.vars=base::all.vars(as.formula(formulas[idx]))
    mod.fit=stats::glm(formulas[idx],family="gaussian",data=data)
    mod.summary=summary(mod.fit)$coefficients
    mod.coefs=rownames(mod.summary)

    sub.summary=mod.summary[mod.coefs%in%tr.vars,,drop=F] #regression summary info for treat.vars
    #if no bonferroni.npvals supplied, use number of coefficients for treat.vars
    estimates=base::unname(sub.summary[,1])

    control.responses=unlist(data[data$control==1,mod.vars[1]])

    if(is.null(stderr.func)==TRUE){ #if no stderr.func provided, use model summary std.err and pvalue
      se.num=sub.summary[,2]
      pval=base::unname(sub.summary[,4])
    }else{ #if stderr.func provided,  get std.err and pvalue from it
      se.num=stderr.func(mod.fit)
      se.num=se.num[names(se.num)%in%c(treat.vars)]
      z=abs(estimates)/se.num
      pval=pnorm(z,lower.tail=FALSE)
    }
    #std. Err for control mean
    control.se=base::sqrt(stats::var(control.responses,na.rm=T)/sum(!is.na(control.responses)))

    #warning("before tab.out")
    #initalize table
    tab.out=data.frame("ITT"=c(estimates,mean(control.responses,na.rm=T)),
                       "StdErr"=c(se.num,control.se),
                       "Pvalue"=c(pval,NA),
                       "Treatment"=c(rownames(sub.summary),"control"),
                       "Model"=rep(mod.nms[idx],length(pval)+1),
                       "CI.overlap"=rep(NA,length(pval)+1),
                       "Abs.Err"=rep(NA,length(pval)+1),
                       "Rel.Err"=rep(NA,length(pval)+1))
    if(is.null(true.vals)==FALSE){ #if true value provided, add Abs.Err and Rel Err
      true.vals=c(true.vals,mean(control.responses))
      tab.out$Abs.Err[seq(1,length(true.vals))]=abs(tab.out$ITT[seq(1,length(true.vals))]-true.vals)
      tab.out$Rel.Err[seq(1,length(true.vals))]=tab.out$Abs.Err[seq(1,length(true.vals))]/true.vals
    }
    #get upper and lower confidence interval bounds
    zalph2=qnorm(min(ci.alpha/2,(1-ci.alpha)/2),lower.tail = F)
    conf.ci.up=tab.out$ITT+(zalph2*tab.out$StdErr)
    conf.ci.lw=tab.out$ITT-(zalph2*tab.out$StdErr)
    #return table and computation time
    return(list(tab.out,(proc.time()-start.time.internal)[[3]]))
  }


  if(length(unique(reg.models))>1){ #if the regression formulas are difference
    #get confidential data results for each regression formula
    conf.outs=parallel::mclapply(seq(1,n.regs),get_conf_treff,mc.cores = num.cores)
    #warning("after conf.out")
    conf.treffs=dplyr::bind_rows(lapply(conf.outs,function(x)data.frame(x[[1]])))
    conf.comp.times=t(data.frame(sapply(conf.outs,"[[",2)))
  }else{ #if only 1 unique regression formula, get the results for confidential data with that formula
    conf.outs=get_conf_treff(1)
    conf.treffs=data.frame(conf.outs[[1]])
    conf.comp.times=cbind(data.frame("c1"=conf.outs[[2]]),data.frame(matrix(NA,nrow=1,ncol=n.regs-1)))
  }
  #warning(paste("after conf.treffs. dim of conf.comp.times=",paste0(dim(conf.comp.times),collapse=", ")))
  colnames(conf.comp.times)=mod.names
  row.names(conf.comp.times)="fit.conf.model"
  conf.treffs$Method="Confidential"
  #warning("after post-process conf.treffs")


  get.synths.start=proc.time() #start time for all sanitizing algorithms
  comp.times.summary=list("Fit Confidential Model(s)"=(get.synths.start-fit.conf.models.start)[[3]])
  print("Starting Synthetic Data Methods")

  simdata=simdata[,all.mod.vars]
  if(only.full.model.based==FALSE){
  ############# Multivariate Histogram Fully Synthetic epsilon=1 ############################
  #Internal Function: to get multivariate histogram synthetic data for each model
  get_synth_mvhist=function(idx,data=simdata,
                            formulas=reg.models,
                            priv.budget.eps=synthdata.budget.eps,
                            priv.budget.del=synthdata.budget.del,
                            cont.vars=cont.as.cont,
                            bparam=bins.param,
                            use.cont.noise=use.continuous.noise,
                            tr.vars=treat.vars,
                            mod.nms=mod.names,
                            compare.table=conf.treffs,
                            ci.alpha=1-ci.confidence,
                            true.vals=true.coef.vals,
                            stderr.func=se.func,
                            nobs=n.obs,
                            #cont.lims=continuous.limits,std.vars=standardize.col,
                            diag.file=diagnostic.file
  ){
    mod.vars=base::all.vars(as.formula(formulas[idx]))

    #first entry of mod.vars is "~"
    cont.vars=cont.vars[cont.vars%in%mod.vars[seq(2,length(mod.vars))]]
    cont.vars=cont.vars[!(cont.vars%in%c(tr.vars))]
    #std.vars=std.vars[std.vars%in%mod.vars[seq(2,length(mod.vars))]]
    #std.vars=std.vars[!(std.vars%in%c(tr.vars))]
    if(length(cont.vars)==0){ #if no continuous variables, vector of F
      cont.vars=rep(F,ncol(data[,colnames(data)%in%mod.vars,drop=F]))

    }
    if(is.null(nobs)==FALSE){ #subset data if needed
      data=data[seq(1,nobs[idx]),]
    }
    # if(is.null(cont.lims)==FALSE){
    # if(is.null(names(cont.lims))==FALSE){
    #   if(ncol(data[,cont.vars,drop=F])>0){
    #     cont.lims=cont.lims[names(cont.lims)%in% colnames(data[,cont.vars,drop=F])]
    #   }
    # }else if(ncol(data[,cont.vars,drop=F])>0){
    #   cont.lims=cont.lims[[seq(1,ncol(data[,cont.vars,drop=F])>0)]]
    # }
    # }

    #print(paste0("std.vars in get synthmvhist ",idx," are",paste0(std.vars,collapse=", ")))

    mv.out=DPrct::synthdata_perturb_mvhist(data=data[,colnames(data)%in%mod.vars,drop=F],
                                           epsilon=priv.budget.eps[idx],
                                           delta=priv.budget.del[idx],
                                           continuous.vars=cont.vars,
                                           bin.param = bparam,
                                           add.cont.variation=use.cont.noise,
                                           treatment.colname=tr.vars,
                                           assign.type="complete",
                                           return.time=TRUE,
                                           with.treatment=FALSE, #treatment is kept in the mv histogram
                                           within.blocks=FALSE,
                                           blocks=NULL,
                                           block.sizes=NULL)#,
                                           #continuous.limits=cont.lims,
                                           #standardize.cont=std.vars,
                                           #diagnostic.file=diag.file,quietly=T
                                           #)


    mv.synth.comp.time=attr(mv.out,"comp.time") #computation time

    #add control and treatment columns
    mv.out$treatment="control"
    for(trv in treat.vars){
      mv.out$treatment[mv.out[,trv]==1]=trv
    }
    mv.out$control=as.factor(ifelse(mv.out$treatment=="control",1,0))

    #add model names if needed
    if(is.null(mod.nms)==TRUE){
      mod.nms=paste0("Model",seq(1,idx))
    }


    mv.out=dplyr::mutate_if(mv.out,is.character,as.factor)

    # Fit model on the sanitized data
    start.fit.time=proc.time()

    mod.fit=stats::glm(unname(formulas[idx]),family="gaussian",
                       data=mv.out[,!(colnames(mv.out)%in%c("control","treatment"))])
    mod.summary=summary(mod.fit)$coefficients
    mod.coefs=rownames(mod.summary)

    sub.summary=mod.summary[mod.coefs%in%tr.vars,,drop=F] #regression summary info for treat.vars
    #if no bonferroni.npvals supplied, use number of coefficients for treat.vars
    estimates=base::unname(sub.summary[,1])
    control.responses=unlist(mv.out[mv.out$control==1,mod.vars[1]])

    if(is.null(stderr.func)==TRUE){
      se.num=sub.summary[,2]
      pval=base::unname(sub.summary[,4])
    }else{
      se.num=stderr.func(mod.fit)
      se.num=se.num[names(se.num)%in%c(treat.vars)]
      z=abs(estimates)/se.num
      pval=pnorm(z,lower.tail=FALSE)
    }
    control.se=base::sqrt(stats::var(control.responses,na.rm=T)/sum(!is.na(control.responses)))


    tab.out=data.frame("ITT"=c(estimates,mean(control.responses,na.rm=T)),
                       "StdErr"=c(se.num,control.se),
                       "Pvalue"=c(pval,NA),
                       "Treatment"=c(rownames(sub.summary),"control"),
                       "Model"=rep(mod.nms[idx],length(pval)+1),
                       "CI.overlap"=rep(NA,length(pval)+1),
                       "Abs.Err"=rep(NA,length(pval)+1),
                       "Rel.Err"=rep(NA,length(pval)+1))
    if(is.null(compare.table)==FALSE){
      if(nrow(compare.table)!=nrow(tab.out)){
        compare.table=compare.table[compare.table$Model==mod.nms[idx],]
      }
      zalph2=qnorm(min(ci.alpha/2,(1-ci.alpha)/2),lower.tail = F)
      compare.up=compare.table$ITT+(zalph2*compare.table$StdErr)
      compare.lw=compare.table$ITT-(zalph2*compare.table$StdErr)
      hybrid.up=tab.out$ITT+(zalph2*tab.out$StdErr)
      hybrid.lw=tab.out$ITT-(zalph2*tab.out$StdErr)
      tab.out$CI.overlap=DPrct::CI_overlap(compare.lw,compare.up,hybrid.lw,hybrid.up)
    }
    if(is.null(true.vals)==FALSE){
      true.vals=c(true.vals,mean(data[data$control==1,mod.vars[1]],rm.na=T))
      tab.out$Abs.Err=abs(tab.out$ITT-true.vals)
      tab.out$Rel.Err=tab.out$Abs.Err/true.vals
    }

    return(list(tab.out,c(mv.synth.comp.time[c(1,3,4)],"fit.san.model.utility.measures"=(proc.time()-start.fit.time)[[3]]),mv.out))
  }

  mv.outs=parallel::mclapply(seq(1,n.regs),get_synth_mvhist,mc.cores = num.cores)
  mv.comp.times=sapply(mv.outs,"[[",2)
  colnames(mv.comp.times)=mod.names
  mv.treffs=dplyr::bind_rows(lapply(mv.outs,function(x)data.frame(x[[1]])))
  mv.treffs$Method="MV Histogram"


  comp.times.summary=c(comp.times.summary,
                       list("All MV Histogram"=(proc.time()-get.synths.start)[[3]]))

  print("Done with Perturb MV Histogram")
  ################ Hybrid DP epsilon=1 #######################
  all.y.hybrid.start=proc.time()
  get_synth_hybriddp=function(idx,
                              synth.data=NULL,
                              formulas=reg.models,
                              confidential.data=simdata,
                              tr.vars=treat.vars,
                              priv.budget.eps=synthdata.budget.eps,
                              priv.budget.del=synthdata.budget.del,
                              continuous.vars=cont.as.cont,
                              bparam=bins.param,
                              add.cont.noise=use.continuous.noise,
                              mod.nms=mod.names,
                              compare.table=conf.treffs,
                              ci.alpha=1-ci.confidence,
                              true.vals=true.coef.vals,
                              stderr.func=se.func,
                              nobs=n.obs,
                              #cont.lims=continuous.limits,
                              #std.vars=standardize.col,
                              diag.file=diagnostic.file

  ){
    if(is.null(nobs)==FALSE){
      confidential.data=confidential.data[seq(1,nobs[idx]),]
    }
    mod.vars=base::all.vars(as.formula(formulas[idx]))
    pred.vars=mod.vars[seq(2,length(mod.vars))]
    pred.vars=pred.vars[!(pred.vars%in%c(tr.vars))]
    continuous.vars=continuous.vars[continuous.vars%in%pred.vars]
    #std.vars=std.vars[std.vars%in%pred.vars]
    if(length(continuous.vars)==0){ #if no continuous variables, vector of F
      continuous.vars=rep(F,ncol(confidential.data[,colnames(confidential.data)%in%mod.vars,drop=F]))
    }
    covar.data=confidential.data[,colnames(confidential.data)%in%pred.vars]

    # if(is.null(cont.lims)==FALSE){
    #   if(is.null(names(cont.lims))==FALSE){
    #     if(ncol(covar.data[,continuous.vars,drop=F])>0){
    #       cont.lims=cont.lims[names(cont.lims)%in% colnames(covar.data[,continuous.vars,drop=F])]
    #     }
    #   }else if(ncol(covar.data[,continuous.vars,drop=F])>0){
    #     cont.lims=cont.lims[seq(1,ncol(covar.data[,continuous.vars,drop=F])>0)]
    #   }
    # }

    #print(paste0("std.vars in get synthhybrid ",idx," are",paste0(std.vars,collapse=", ")))
    mv.covariates=DPrct::synthdata_perturb_mvhist(data=covar.data,
                                                  epsilon=priv.budget.eps[idx],
                                                  delta=priv.budget.del[idx],
                                                  continuous.vars=continuous.vars,
                                                  num.bin=NULL,
                                                  bin.param = bparam,
                                                  add.cont.variation=add.cont.noise,
                                                  treatment.colname="treatment",
                                                  assign.type="complete",
                                                  return.time=TRUE,
                                                  with.treatment=TRUE,
                                                  within.blocks=FALSE,
                                                  blocks=NULL,
                                                  block.sizes=NULL,
                                                  conditions=c(tr.vars,"control"))#,
                                                  #continuous.limits=cont.lims,
                                                  #standardize.cont=std.vars,diagnostic.file=diag.file,quietly=T)
    for(tvar in c(tr.vars,"control")){
      mv.covariates[,tvar]=ifelse(mv.covariates$treatment==tvar,1,0)
    }
    #liberia.sub[,treat.vars]=apply(liberia.sub[,treat.vars],2,as.factor)
    mv.covariates=mv.covariates[,!(colnames(mv.covariates)%in%c("control"))]
    mv.covariate.comp.time=attr(mv.covariates,"comp.time")

    temp.out=DPrct::hybrid_synth(formula=formulas[idx],
                                 confidential.data=confidential.data[,colnames(confidential.data)%in%mod.vars,drop=F],
                                 synth.data=mv.covariates,
                                 continuous.vars=continuous.vars,
                                 num.bin=NULL,
                                 bin.param=bparam,
                                 assign.type="complete",
                                 blocks=NULL,
                                 within.blocks=FALSE,
                                 treatment.var=tr.vars,
                                 rseed=NA,
                                 returntypes=c("synth.data","comp.time"))

    hybrid.out=temp.out[[1]]
    hybrid.out$treatment="control"
    for(trv in tr.vars){
      hybrid.out$treatment[hybrid.out[,trv]==1]=trv
      hybrid.out[,trv]=as.numeric(as.character(hybrid.out[,trv]))
    }
    hybrid.out$control=as.factor(ifelse(hybrid.out$treatment=="control",1,0))
    start.fit.time=proc.time()
    # Fit SynthModel
    if(is.null(mod.nms)==TRUE){
      mod.nms=paste0("Model",seq(1,idx))
    }
    mod.fit=stats::glm(formulas[idx],family="gaussian",data=hybrid.out)
    mod.summary=summary(mod.fit)$coefficients
    mod.coefs=rownames(mod.summary)

    sub.summary=mod.summary[mod.coefs%in%tr.vars,,drop=F] #regression summary info for treat.vars
    # #if no bonferroni.npvals supplied, use number of coefficients for treat.vars
    estimates=base::unname(sub.summary[,1])
    control.responses=unlist(hybrid.out[hybrid.out[,"control"]==1,mod.vars[1]])

    if(is.null(stderr.func)==TRUE){
      se.num=sub.summary[,2]
      pval=base::unname(sub.summary[,4])
    }else{
      se.num=stderr.func(mod.fit)
      se.num=se.num[names(se.num)%in%c(treat.vars)]
      z=abs(estimates)/se.num
      pval=pnorm(z,lower.tail=FALSE)
    }
    control.se=base::sqrt(stats::var(control.responses,na.rm=T)/sum(!is.na(control.responses)))


    tab.out=data.frame("ITT"=c(estimates,mean(control.responses,na.rm=T)),
                       "StdErr"=c(se.num,control.se),
                       "Pvalue"=c(pval,NA),
                       "Treatment"=c(rownames(sub.summary),"control"),
                       "Model"=rep(mod.nms[idx],length(pval)+1),
                       "CI.overlap"=rep(NA,length(pval)+1),
                       "Abs.Err"=rep(NA,length(pval)+1),
                       "Rel.Err"=rep(NA,length(pval)+1))
    if(is.null(compare.table)==FALSE){
      if(nrow(compare.table)!=nrow(tab.out)){
        compare.table=compare.table[compare.table$Model==mod.nms[idx],]
      }
      zalph2=qnorm(min(ci.alpha/2,(1-ci.alpha)/2),lower.tail = F)
      compare.up=compare.table$ITT+(zalph2*compare.table$StdErr)
      compare.lw=compare.table$ITT-(zalph2*compare.table$StdErr)
      hybrid.up=tab.out$ITT+(zalph2*tab.out$StdErr)
      hybrid.lw=tab.out$ITT-(zalph2*tab.out$StdErr)
      tab.out$CI.overlap=DPrct::CI_overlap(compare.lw,compare.up,hybrid.lw,hybrid.up)
    }
    if(is.null(true.vals)==FALSE){
      true.vals=c(true.vals,mean(confidential.data[confidential.data$control==1,mod.vars[1]],rm.na=T))
      tab.out$Abs.Err=abs(tab.out$ITT-true.vals)
      tab.out$Rel.Err=tab.out$Abs.Err/true.vals
    }

    tempcomp=temp.out[[2]]
    comp.times=c("mv.hist.san.predictors"=unname(mv.covariate.comp.time[4]),
                 "fit.conf.model"=tempcomp$fit.model.time,
                 "san.response.time"=tempcomp$san.response.time,
                 "total.synth"=tempcomp$total_time,
                 "fit.san.model.utility.measures"=(proc.time()-start.fit.time)[[3]])
    return(list(tab.out,comp.times,hybrid.out))
  }

  hybrid.outs=parallel::mclapply(seq(1,n.regs),get_synth_hybriddp,mc.cores = num.cores)
  hybrid.comp.times=sapply(hybrid.outs,"[[",2)
  colnames(hybrid.comp.times)=mod.names
  hybrid.treffs=dplyr::bind_rows(lapply(hybrid.outs,function(x)data.frame(x[[1]])))
  hybrid.treffs$Method="Hybrid"


  comp.times.summary=c(comp.times.summary,
                       list("Hybrid"=(proc.time()-all.y.hybrid.start)[[3]]))

  print("Done with Hybrid Method")
  }

  ####################

  if(include.full.model.based==TRUE){
    ################ My Full DP Synthetic epsilon=1 #######################
    synthdata.start=proc.time()



    #get synthetic

    #internal function to get synthetic data set, fit model.
    #Returns table of ITT for treatments, with ci overlap, relative error, abs error; computation times, and the dataset
    get_synth_dpmodel=function(idx,
                               synth.data=NULL,
                               formulas=reg.models,
                               confidential.data=simdata,
                               tr.vars=treat.vars,
                               mv.portion.budget=mv.prop.budget,
                               per.response.eps=synthdata.budget.eps,
                               per.response.del=synthdata.budget.del,
                               bd.sd=bound.sds,
                               bd.mean=bound.means,
                               nits=n.iters.fulldp,
                               resvar.prop=win.y.resvar,
                               trcoef.prop=win.y.trcoef,
                               budget.prop.list=win.y.budget.props,
                               alph=range.alpha,
                               continuous.vars=cont.as.cont,
                               bparam=bins.param,
                               add.cont.noise=use.continuous.noise,
                               mod.nms=mod.names,
                               compare.table=conf.treffs,
                               ci.alpha=1-ci.confidence,
                               true.vals=true.coef.vals,
                               stderr.func=se.func,
                               nobs=n.obs,
                               #cont.lims=continuous.limits,
                               #std.vars=standardize.col,
                               diag.file=diagnostic.file

    ){
      if(is.null(nobs)==FALSE){
        confidential.data=confidential.data[seq(1,nobs[idx]),]
      }
      mod.vars=base::all.vars(as.formula(formulas[idx]))
      pred.vars=mod.vars[seq(2,length(mod.vars))]
      per.y.prop=(1-mv.prop.budget)
      win.y.x.coefs=(1-resvar.prop-trcoef.prop)
      ############
      n.x.coefs=ncol(
        stats::model.matrix(~.,
                            data=confidential.data[,colnames(confidential.data)%in%pred.vars]))-length(tr.vars)
      n.treat.coefs=length(tr.vars)
      win.y.budget.props=c(rep(trcoef.prop/n.treat.coefs,n.treat.coefs),
                           rep(win.y.x.coefs/n.x.coefs,n.x.coefs),resvar.prop)
      pred.vars=mod.vars[seq(2,length(mod.vars))]
      pred.vars=pred.vars[!(pred.vars%in%c(tr.vars))]
      continuous.vars=continuous.vars[continuous.vars%in%pred.vars]
      #std.vars=std.vars[std.vars%in%pred.vars]
      if(length(continuous.vars)==0){ #if no continuous variables, vector of F
        continuous.vars=rep(F,ncol(confidential.data[,colnames(confidential.data)%in%mod.vars,drop=F]))
      }
      covar.data=confidential.data[,colnames(confidential.data)%in%pred.vars]


      # if(is.null(cont.lims)==FALSE){
      #   if(is.null(names(cont.lims))==FALSE){
      #     if(ncol(covar.data[,continuous.vars,drop=F])>0){
      #       cont.lims=cont.lims[names(cont.lims)%in% colnames(covar.data[,continuous.vars,drop=F])]
      #     }
      #   }else if(ncol(covar.data[,continuous.vars,drop=F])>0){
      #     cont.lims=cont.lims[seq(1,ncol(covar.data[,continuous.vars,drop=F])>0)]
      #   }
      # }

      #print(paste0("std.vars in get synthdpmodel ",idx," are",paste0(std.vars,collapse=", ")))
      mv.covariates=DPrct::synthdata_perturb_mvhist(data=covar.data,
                                                    epsilon=per.response.eps[idx]*mv.portion.budget,
                                                    delta=per.response.del[idx]*mv.portion.budget,
                                                    continuous.vars=continuous.vars,
                                                    num.bin=NULL,
                                                    bin.param = bparam,
                                                    add.cont.variation=add.cont.noise,
                                                    treatment.colname="treatment",
                                                    assign.type="complete",
                                                    return.time=TRUE,
                                                    with.treatment=TRUE,
                                                    within.blocks=FALSE,
                                                    blocks=NULL,
                                                    block.sizes=NULL,
                                                    conditions=c(tr.vars,"control"))#,
                                                    #continuous.limits=cont.lims,
                                                    #standardize.cont=std.vars,diagnostic.file=diag.file,quietly=T)
      for(tvar in c(tr.vars,"control")){
        mv.covariates[,tvar]=ifelse(mv.covariates$treatment==tvar,1,0)
      }
      #liberia.sub[,treat.vars]=apply(liberia.sub[,treat.vars],2,as.factor)
      mv.covariates=mv.covariates[,!(colnames(mv.covariates)%in%c("control"))]
      mv.covariate.comp.time2=attr(mv.covariates,"comp.time")
      # if(length(bd.mean)<n.x.coefs){
      #   if(is.list(bd.mean)==FALSE){
      #     bd.mean=as.list(bd.mean)
      #   }
      #   bd.mean=rep(bd.mean,(n.x.coefs%/%length(bd.mean)+1))
      # }
      #warning("in sim_compare_method before dp_synthdata")
      temp.out=DPrct::dp_synthdata(formula=formulas[idx],
                                   confidential.data=confidential.data[,colnames(confidential.data)%in%mod.vars,drop=F],
                                   synth.data=mv.covariates,
                                   synth.epsilon=NULL,#per.response.eps[idx]*mv.portion.budget,
                                   synth.delta=NULL,#per.response.del[idx]*mv.portion.budget,
                                   continuous.vars=continuous.vars,
                                   num.bin=NULL,
                                   bin.param=bparam,
                                   assign.type="simple",
                                   blocks=NULL,
                                   within.blocks=FALSE,
                                   epsilon.list=as.list(per.response.eps[idx]*(1-mv.portion.budget)*budget.prop.list),
                                   delta.list=as.list(per.response.del[idx]*(1-mv.portion.budget)*budget.prop.list),
                                   bd.sd.list=bd.sd,
                                   bd.mean.list=bd.mean,
                                   alphas.list=alph,
                                   num.iters=nits,
                                   treatment.var=tr.vars,
                                   rseed=NA,
                                   return.time=TRUE,
                                   return.confidential.table=FALSE,
                                   return.san.summary=FALSE,
                                   use.san.residerror=TRUE,
                                   return.treatment.pvalue=FALSE,
                                   return.treatment.CI=FALSE,
                                   treat.ppart.list=NULL)

      dpmod.out=temp.out[[1]]
      dpmod.out$treatment="control"
      for(trv in tr.vars){
        dpmod.out$treatment[dpmod.out[,trv]==1]=trv
        dpmod.out[,trv]=as.numeric(as.character(dpmod.out[,trv]))
      }
      dpmod.out$control=as.factor(ifelse(dpmod.out$treatment=="control",1,0))

      # Fit SynthModel
      if(is.null(mod.nms)==TRUE){
        mod.nms=paste0("Model",seq(1,idx))
      }

      start.fit.time=proc.time()
      mod.fit=stats::glm(formulas[idx],family="gaussian",data=dpmod.out)
      mod.summary=summary(mod.fit)$coefficients
      mod.coefs=rownames(mod.summary)

      sub.summary=mod.summary[mod.coefs%in%tr.vars,,drop=F] #regression summary info for treat.vars
      # #if no bonferroni.npvals supplied, use number of coefficients for treat.vars
      estimates=base::unname(sub.summary[,1])
      control.responses=unlist(dpmod.out[dpmod.out[,"control"]==1,mod.vars[1]])


      if(is.null(stderr.func)==TRUE){
        se.num=sub.summary[,2]
        pval=base::unname(sub.summary[,4])
      }else{
        se.num=stderr.func(mod.fit)
        se.num=se.num[names(se.num)%in%c(treat.vars)]
        z=abs(estimates)/se.num
        pval=pnorm(z,lower.tail=FALSE)
      }
      control.se=base::sqrt(stats::var(control.responses,na.rm=T)/sum(!is.na(control.responses)))



      tab.out=data.frame("ITT"=c(estimates,mean(control.responses,na.rm=T)),
                         "StdErr"=c(se.num,control.se),
                         "Pvalue"=c(pval,NA),
                         "Treatment"=c(rownames(sub.summary),"control"),
                         "Model"=rep(mod.nms[idx],length(pval)+1),
                         "CI.overlap"=rep(NA,length(pval)+1),
                         "Abs.Err"=rep(NA,length(pval)+1),
                         "Rel.Err"=rep(NA,length(pval)+1))
      if(is.null(compare.table)==FALSE){
        if(nrow(compare.table)!=nrow(tab.out)){
          compare.table=compare.table[compare.table$Model==mod.nms[idx],]
        }
        zalph2=qnorm(min(ci.alpha/2,(1-ci.alpha)/2),lower.tail = F)
        compare.up=compare.table$ITT+(zalph2*compare.table$StdErr)
        compare.lw=compare.table$ITT-(zalph2*compare.table$StdErr)
        dpmod.up=tab.out$ITT+(zalph2*tab.out$StdErr)
        dpmod.lw=tab.out$ITT-(zalph2*tab.out$StdErr)
        tab.out$CI.overlap=DPrct::CI_overlap(compare.lw,compare.up,dpmod.lw,dpmod.up)
      }
      if(is.null(true.vals)==FALSE){
        true.vals=c(true.vals,mean(confidential.data[confidential.data$control==1,mod.vars[1]],rm.na=T))
        tab.out$Abs.Err=abs(tab.out$ITT-true.vals)
        tab.out$Rel.Err=tab.out$Abs.Err/true.vals
      }


      tempcomp=temp.out[[2]]
      comp.times=c("mv.hist.san.predictors"=unname(mv.covariate.comp.time2[4]),
                   "fit.conf.model"=tempcomp$conf.model.fit,
                   "iter.proxy.fit"=tempcomp$iter.proxy.fit,
                   "san.params"=tempcomp$san.summary.time,
                   "san.response.time"=tempcomp$san.response.time,
                   "total.synth"=tempcomp$total_time,
                   "fit.san.model.utility.measures"=(proc.time()-start.fit.time)[[3]])

      return(list(tab.out,comp.times,dpmod.out))
    }

    fulldp.outs=parallel::mclapply(seq(1,n.regs),get_synth_dpmodel,mc.cores = num.cores)
    fulldp.comp.times=sapply(fulldp.outs,"[[",2)
    colnames(fulldp.comp.times)=mod.names
    fulldp.treffs=dplyr::bind_rows(lapply(fulldp.outs,function(x)data.frame(x[[1]])))
    fulldp.treffs$Method="DP Model-based"


    comp.times.summary=c(comp.times.summary,
                         list("All DP Model-based"=(proc.time()-get.synths.start)[[3]]))

    print("Done with Full DP Model-Based")
  }

  ##############
  comp.times.summary=c(comp.times.summary,
                       list("All Synthetic Data Methods"=(proc.time()-get.synths.start)[[3]]))

  table.start=proc.time()
  if(only.full.model.based==FALSE){
  treff.table=rbind(conf.treffs,mv.treffs,hybrid.treffs,fulldp.treffs)
  comp.times.tab=rbind(conf.comp.times,mv.comp.times,hybrid.comp.times,fulldp.comp.times)
  list.of.synthdatas=list("MV Histogram"=lapply(mv.outs,"[[",3),
                          "Hybrid"=lapply(hybrid.outs,"[[",3),
                          "DP Model-based"=lapply(fulldp.outs,"[[",3))

  }else{
    treff.table=rbind(conf.treffs,fulldp.treffs)
    comp.times.tab=rbind(conf.comp.times,fulldp.comp.times)
    list.of.synthdatas=list("DP Model-based"=lapply(fulldp.outs,"[[",3))
  }
  if(is.na(sim)==FALSE){
    treff.table$Sim=sim
    print(paste("Completed Simulation:",sim))
  }

  end.time=proc.time()
  comp.times.summary=c(comp.times.summary,list("Make Table"=(end.time-table.start)[[3]],
                                               "Total Time"=(end.time-start.time)[[3]])
  )
  attr(treff.table,"all.comp.times")=list("Summary"=comp.times.summary,
                                          "Per.Model"=comp.times.tab)

  return(list(treff.table,comp.times.summary,comp.times.tab,list.of.synthdatas,simdata))
}
