#no noise dpmb


###### Adapt sim_compare_method function
nonoise_sim_compare_methods=function(simdata,
                             reg.models,
                             bins.param,use.continuous.noise=T,
                             treat.vars,
                             cont.as.cont,
                             n.iters.fulldp=NA,
                             mod.names,treat.names,
                             bonferroni.npvals=NA,
                             ci.confidence=0.95,
                             rseed=NA,
                             sim=NA,
                             include.full.model.based=TRUE,
                             only.full.model.based=FALSE,
                             true.coef.vals=NULL,
                             se.func=NULL,
                             n.obs=NULL,
                             se.normal){

  start.time=proc.time() #start time
  if(is.na(rseed)==FALSE){ #set seed if needed
    set.seed(rseed)
  }
  if(include.full.model.based==TRUE){ #check inputs for full.model.based
    stopifnot(is.na(n.iters.fulldp)==FALSE)
  }

  all.mod.vars=unique(unlist(lapply(reg.models,function(x)all.vars(as.formula(x)))))
  #subset data
  simdata=simdata[,colnames(simdata)%in%c(all.mod.vars,"treatment","control"),drop=F]
  simdata[,treat.vars]=apply(simdata[,treat.vars,drop=F],2,function(x)as.numeric(as.character(x)))

  ### Constructing regression model formulas for each model
  n.regs=length(reg.models)

  fit.conf.models.start=proc.time() #start time for confidential data
  comp.times.summary=list("Set Up Time"=(fit.conf.models.start-start.time)[[3]])

  #function to get treatment effects, confidence intervals, absolute, error, and relative error

  if(length(unique(reg.models))>1){ #if the regression formulas are different
    #get confidential data results for each regression formula
    conf.outs=parallel::mclapply(seq(1,n.regs),function(i)get_conf_treff(i,data=simdata,
                                 formulas=reg.models,
                                 tr.vars=treat.vars,
                                 mod.nms=mod.names,
                                 ci.alpha=1-ci.confidence,
                                 true.vals=true.coef.vals,
                                 stderr.func=se.func,
                                 nobs=n.obs))
    #warning("after conf.out")
    conf.treffs=dplyr::bind_rows(lapply(conf.outs,function(x)data.frame(x[[1]])))
    conf.comp.times=t(data.frame(sapply(conf.outs,"[[",2)))
  }else{ #if only 1 unique regression formula, get the results for confidential data with that formula
    conf.outs=get_conf_treff(1,data=simdata,
                             formulas=reg.models,
                             tr.vars=treat.vars,
                             mod.nms=mod.names,
                             ci.alpha=1-ci.confidence,
                             true.vals=true.coef.vals,
                             stderr.func=se.func,
                             nobs=n.obs)
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
  #print("Starting Synthetic Data Methods")

  if(only.full.model.based==FALSE){
    ############# Multivariate Histogram Fully Synthetic############################
    #Internal Function: to get multivariate histogram synthetic data for each model
    get_synth_mvhist=function(idx,data=simdata,
                              formulas=reg.models,
                              cont.vars=cont.as.cont,
                              bparam=bins.param,
                              use.cont.noise=use.continuous.noise,
                              tr.vars=treat.vars,
                              mod.nms=mod.names,
                              compare.table=conf.treffs,
                              ci.alpha=1-ci.confidence,
                              true.vals=true.coef.vals,
                              stderr.func=se.func,
                              nobs=n.obs
    ){
      mod.vars=base::all.vars(as.formula(formulas[idx]))
      #predictor.vars=mod.vars[-1]
      #first entry of mod.vars is "~"
      cont.vars=cont.vars[cont.vars%in%mod.vars]
      cont.vars=cont.vars[!(cont.vars%in%c(tr.vars,"control"))]
      if(length(cont.vars)==0){ #if no continuous variables, vector of F
        cont.vars=rep(F,ncol(data[,colnames(data)%in%mod.vars,drop=F]))
      }
      if(is.null(nobs)==FALSE){ #subset data if needed
        data=data[seq(1,nobs[idx]),]
      }
      mv.out=DPrct::synthdata_perturb_mvhist(data=data[,colnames(data)%in%mod.vars,drop=F],
                                             epsilon=1,
                                             delta=0,
                                             continuous.vars=cont.vars,
                                             bin.param = bparam,
                                             add.cont.variation=use.cont.noise,
                                             treatment.colname=tr.vars,
                                             assign.type="complete",
                                             return.time=TRUE,
                                             with.treatment=FALSE, #treatment is kept in the mv histogram
                                             within.blocks=FALSE,
                                             blocks=NULL,
                                             block.sizes=NULL,
                                             perturb=F)

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
      # Fit model on the sanitized data
      start.fit.time=proc.time()
      mod.fit=stats::glm(formulas[idx],family="gaussian",data=data)
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
        mv.up=tab.out$ITT+(zalph2*tab.out$StdErr)
        mv.lw=tab.out$ITT-(zalph2*tab.out$StdErr)
        tab.out$CI.overlap=DPrct::CI_overlap(compare.lw,compare.up,mv.lw,mv.up)
      }
      if(is.null(true.vals)==FALSE){
        true.vals=c(true.vals,mean(data[data$control==1,mod.vars[1]],rm.na=T))
        tab.out$Abs.Err=abs(tab.out$ITT-true.vals)
        tab.out$Rel.Err=tab.out$Abs.Err/true.vals
      }

      return(list(tab.out,c(mv.synth.comp.time[c(1,3,4)],
                            "fit.san.model.utility.measures"=(proc.time()-start.fit.time)[[3]]),
                  mv.out))
    }

    mv.outs=parallel::mclapply(seq(1,n.regs),get_synth_mvhist)
    mv.comp.times=sapply(mv.outs,"[[",2)
    colnames(mv.comp.times)=mod.names
    mv.treffs=dplyr::bind_rows(lapply(mv.outs,function(x)data.frame(x[[1]])))
    mv.treffs$Method="MV Histogram"


    comp.times.summary=c(comp.times.summary,
                         list("All MV Histogram"=(proc.time()-get.synths.start)[[3]]))

    #print("Done with Perturb MV Histogram")
    ################ Hybrid DP #######################
    all.y.hybrid.start=proc.time()
    get_synth_hybriddp=function(idx,
                                synth.data=NULL,
                                formulas=reg.models,
                                confidential.data=simdata,
                                tr.vars=treat.vars,
                                continuous.vars=cont.as.cont,
                                bparam=bins.param,
                                add.cont.noise=use.continuous.noise,
                                mod.nms=mod.names,
                                compare.table=conf.treffs,
                                ci.alpha=1-ci.confidence,
                                true.vals=true.coef.vals,
                                stderr.func=se.func,
                                nobs=n.obs

    ){
      if(is.null(nobs)==FALSE){
        confidential.data=confidential.data[seq(1,nobs[idx]),]
      }
      mod.vars=base::all.vars(as.formula(formulas[idx]))
      pred.vars=mod.vars[seq(2,length(mod.vars))]
      pred.vars=pred.vars[!(pred.vars%in%c(tr.vars,"control"))]
      continuous.vars=continuous.vars[continuous.vars%in%pred.vars]
      if(length(continuous.vars)==0){
        continuous.vars=rep(F,ncol(confidential.data[,colnames(confidential.data)%in%mod.vars,drop=F]))
      }
      covar.data=confidential.data[,colnames(confidential.data)%in%pred.vars]
      mv.covariates=DPrct::synthdata_perturb_mvhist(data=covar.data,
                                                    epsilon=1,
                                                    delta=0,
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
                                                    conditions=c(tr.vars,"control"),
                                                    perturb=F)
      for(tvar in c(tr.vars,"control")){
        mv.covariates[,tvar]=ifelse(mv.covariates$treatment==tvar,1,0)
      }
      #liberia.sub[,treat.vars]=apply(liberia.sub[,treat.vars],2,as.factor)
      #mv.covariates=mv.covariates[,!(colnames(mv.covariates)%in%c("control"))]
      mv.covariate.comp.time=attr(mv.covariates,"comp.time")

      start.synth.y=proc.time()
      confidential.data.mod=confidential.data[,colnames(confidential.data)%in%mod.vars,drop=F]
      conf.model<-stats::glm(formula=formulas[idx],family="gaussian",data=confidential.data.mod)
      synth.model.data=synth.data[,colnames(synth.data)%in%mod.vars]
      sim.response=DPrct::simulate_response_glm(conf.model,newdata=mv.covariates[,colnames(mv.covariates)%in%mod.vars[-1],drop=F])

      hybrid.out=cbind(mv.covariates,sim.response)
      colnames(hybrid.out)=c(colnames(hybrid.out)[seq(1,ncol(hybrid.out)-1)],mod.vars[1])
      tempcomp=(proc.time()-start.synth.y)[[3]]

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

      comp.times=c("mv.hist.san.predictors"=unname(mv.covariate.comp.time[4]),
                   "total.synth"=tempcomp)
      return(list(tab.out,comp.times,hybrid.out))
    }

    hybrid.outs=parallel::mclapply(seq(1,n.regs),get_synth_hybriddp)
    hybrid.comp.times=sapply(hybrid.outs,"[[",2)
    colnames(hybrid.comp.times)=mod.names
    hybrid.treffs=dplyr::bind_rows(lapply(hybrid.outs,function(x)data.frame(x[[1]])))
    hybrid.treffs$Method="Hybrid"


    comp.times.summary=c(comp.times.summary,
                         list("Hybrid"=(proc.time()-all.y.hybrid.start)[[3]]))

    #print("Done with Hybrid Method")
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
                               nits=n.iters.fulldp,
                               continuous.vars=cont.as.cont,
                               bparam=bins.param,
                               add.cont.noise=use.continuous.noise,
                               mod.nms=mod.names,
                               compare.table=conf.treffs,
                               ci.alpha=1-ci.confidence,
                               true.vals=true.coef.vals,
                               stderr.func=se.func,
                               nobs=n.obs

    ){
      if(is.null(nobs)==FALSE){
        confidential.data=confidential.data[seq(1,nobs[idx]),]
      }
      mod.vars=base::all.vars(as.formula(formulas[idx]))
      pred.vars=mod.vars[seq(2,length(mod.vars))]

      ############
      n.x.coefs=ncol(
        stats::model.matrix(~.,
                            data=confidential.data[,colnames(confidential.data)%in%pred.vars]))-length(tr.vars)
      n.treat.coefs=length(tr.vars)
      pred.vars=pred.vars[!(pred.vars%in%c(tr.vars))]
      continuous.vars=continuous.vars[continuous.vars%in%pred.vars]
      if(length(continuous.vars)==0){
        continuous.vars=rep(F,ncol(confidential.data[,colnames(confidential.data)%in%mod.vars,drop=F]))
      }
      covar.data=confidential.data[,colnames(confidential.data)%in%pred.vars]
      mv.covariates=DPrct::synthdata_perturb_mvhist(data=covar.data,
                                                    epsilon=1,
                                                    delta=0,
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
                                                    conditions=c(tr.vars,"control"),
                                                    perturb=F)
      for(tvar in c(tr.vars,"control")){
        mv.covariates[,tvar]=ifelse(mv.covariates$treatment==tvar,1,0)
      }

      mv.covariates=mv.covariates[,!(colnames(mv.covariates)%in%c("control"))]
      mv.covariate.comp.time2=attr(mv.covariates,"comp.time")

      temp.out=nonoise_dpmb(formula=formulas[idx],
                                   confidential.data=confidential.data[,colnames(confidential.data)%in%mod.vars,drop=F],
                                   synth.data=mv.covariates,
                                   continuous.vars=continuous.vars,
                                   num.bin=NULL,
                                   bin.param=bparam,
                            ci.alpha=ci.alpha,
                                   assign.type="simple",
                                   blocks=NULL,
                                   within.blocks=FALSE,
                            clusters=NULL,
                                   num.iters=nits,
                                   treatment.var=tr.vars,
                                   rseed=NA,
                                   return.time=TRUE,
                                   return.confidential.table=FALSE,
                                   return.san.summary=FALSE,
                                   use.san.residerror=TRUE)

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
      #warning(paste("tempcomp is",paste(tempcomp,collapse=", ")))
      comp.times=c("mv.hist.san.predictors"=unname(mv.covariate.comp.time2[4]),
                   "fit.conf.model"=tempcomp$conf.model.fit,
                   "iter.proxy.fit"=tempcomp$iter.proxy.fit,
                   "san.params"=tempcomp$san.summary.time,
                   "san.response.time"=tempcomp$san.response.time,
                   "total.synth"=tempcomp$total_time,
                   "fit.san.model.utility.measures"=(proc.time()-start.fit.time)[[3]])

      return(list(tab.out,comp.times,dpmod.out))
    }

    fulldp.outs=parallel::mclapply(seq(1,n.regs),get_synth_dpmodel)
    fulldp.comp.times=sapply(fulldp.outs,"[[",2)
    colnames(fulldp.comp.times)=mod.names
    fulldp.treffs=dplyr::bind_rows(lapply(fulldp.outs,function(x)data.frame(x[[1]])))
    fulldp.treffs$Method="DP Model-based"


    comp.times.summary=c(comp.times.summary,
                         list("All DP Model-based"=(proc.time()-get.synths.start)[[3]]))

    #print("Done with Full DP Model-Based")
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

############


main_simulation_nonoise<-function(combos.NCov,nobs,nsims,
                                       response="y",treat.vars="t1",
                                       add.formulas=NULL,add.names=NULL,add.ncat=NULL,add.ncont=NULL,
                                       intercept,treat.effect,cat.gseq,cont.gseq,resid.sd,
                                       cat.probs,cont.funcs,cont.params,
                                       sim.data.rseed=1,
                                       conf.suffix="sim_models",file.suffix="",
                                       output.folder="~",
                                       conf.plot=list("pt.sz1"=1,"ln.sz1"=1.5,"bs.sz1"=12,
                                                      "pt.sz2"=1,"ln.sz2"=1.5,"bs.sz2"=15,
                                                      "ncol1"=3,"ncol2"=2,
                                                      "width1"=4,"height1"=2.5,
                                                      "width2"=6,"height1"=3.5),
                                       bins.param,use.continuous.noise=TRUE,
                                       n.iters.fulldp,
                                       ci.confidence=0.95,stderr.func=NULL,
                                       san.rseed=1

){
  main.start=proc.time()
  ntreat=length(treat.effect)
  reg.formulas.sim2=sim2_models_list(response=response,treat.vars=treat.vars,combos.NCov=combos.NCov,
                                     add.formulas=add.formulas,add.names=add.names)
  mod.names=names(reg.formulas.sim2)

  ncat=max(combos.NCov$NCat)
  ncont=max(combos.NCov$NCont)
  mod.coefs=sim_coefs(ncat=ncat,ncat.groups=sapply(cat.probs,length),ncont=ncont,
                      cat.coef.start=cat.gseq[1],cat.coef.ratio=cat.gseq[2],
                      cont.coef.start=cont.gseq[1],cont.coef.ratio=cont.gseq[2],
                      treat.effect=treat.effect,intercept=intercept)
  model.params=list("model.coefs"=mod.coefs,"residual.err.sd"=resid.sd)
  cont.vars.all=c("y",paste0("x",seq(1,ncont)))



  gen.sim.start=proc.time()
  times.list=list("SetUp"=(gen.sim.start-main.start)[[3]])
  sim.dfs=generate_simulated_dataset(num.cat=ncat,num.cont=ncont,num.treat=ntreat,nobs=nobs,
                                     mod.coefs=mod.coefs,residual.sd=resid.sd,
                                     cat.probs=cat.probs,cont.funcs=cont.funcs,
                                     cont.params=cont.params,
                                     rseed=sim.data.rseed)


  conf.analysis.start=proc.time()
  times.list=c(times.list,list("Simulate Data"=(conf.analysis.start-gen.sim.start)[[3]]))


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

  conf.analysis.end=proc.time()
  times.list=c(times.list,list("conf.analysis"=(conf.analysis.end-conf.analysis.start)[[3]]))
  save(conf.plots,add.plots,conf.summaries,sim.dfs,
       reg.formulas.sim2,model.params,times.list,
       file=paste0(output.folder,"/simulations/checkpoint_",conf.suffix,"_",file.suffix,".Rda"))


  #print("done with conf and simdata")
  start.simulations=proc.time()
  sim2.out=parallel::mclapply(seq(1,nsims),function(idx)nonoise_sim_compare_methods(simdata=sim.dfs,
                                                                            reg.models=reg.formulas.sim2,
                                                                            bins.param=bins.param,
                                                                            use.continuous.noise=use.continuous.noise,
                                                                            treat.vars=treat.vars,
                                                                            cont.as.cont=cont.vars.all,
                                                                          n.iters.fulldp=n.iters.fulldp,
                                                                            mod.names=mod.names,
                                                                            ci.confidence=ci.confidence,
                                                                            rseed=san.rseed+idx,
                                                                            sim=paste0("Sim",idx),
                                                                            include.full.model.based=TRUE,
                                                                            true.coef.vals=c(treat.effect),
                                                                            se.func=stderr.func))



  end.simulations=proc.time()
  times.list=c(times.list,list("simulations"=(end.simulations-start.simulations)[[3]]))

  ITT2=dplyr::bind_rows(lapply(sim2.out,"[[",1))

  if(is.null(add.formulas)==FALSE){
    add.NCov=data.frame("NCat"=add.ncat,"NCont"=add.ncont,"Name"=add.names)
    combos.NCov=dplyr::bind_rows(combos.NCov,add.NCov)
  }
  colnames(combos.NCov)=c("NCat","NCont","Model")
  ITT2=dplyr::left_join(ITT2,combos.NCov,by="Model")
  ITT2$NCovariates=ITT2$NCat+ITT2$NCont
  summary.tab=ITT2%>%group_by(Method,Treatment,Model,NCont,NCat,NCovariates)%>%
    summarise( "Avg.ITT"=mean(ITT),
               "Avg.CI.overlap"=mean(CI.overlap,na.rm=T),
               "Avg.Abs.Err"=mean(Abs.Err,na.rm=T),
               "Avg.Rel.Err"=mean(Rel.Err,na.rm=T),
               "Max.StdErr"=max(StdErr)
    )
  #reorder Method factor levels
  summary.tab$Method=factor(summary.tab$Method,levels=c("Confidential","MV Histogram","Hybrid","DP Model-based"))
  summary.tab$ModelLab=factor(summary.tab$Model,levels=combos.NCov$Model)
  summary.tab.catcont=summary.tab%>%#filter(!grepl("Finite|Bounded",Model))%>%
    arrange(Method)%>%arrange(NCat)%>%arrange(NCont)
  rownames(summary.tab.catcont)=NULL


  algo.param2=list("bin.param"=bins.param,
                   "DPMb.niters"=n.iters.fulldp)

  main.end=proc.time()
  times.list=c(times.list,list("summarize"=(main.end-end.simulations)[[3]],
                               "total"=(main.end-main.start)[[3]]))
  #NOTE: RERUN.
  save.file.name=paste0(output.folder,"/simulations/results_models_",file.suffix,".Rda")
    save(sim2.out,ITT2,summary.tab,summary.tab.catcont,
         combos.NCov,times.list,algo.param2,
         file=save.file.name)
  return(list(ITT2,times.list,summary.tab,algo.param2))
}


