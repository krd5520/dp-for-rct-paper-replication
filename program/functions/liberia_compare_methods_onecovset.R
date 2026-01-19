### otherbudgetsDP list with each element representing a budget. Each element should have mv.epsilon,mv.delta,per.y.eps,per.y.del
liberia_compare_methods_onecovset=function(liberia.sub,
                                 synthdata.budget,mv.bins,use.continuous.noise,
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
                                 se.normal=T,
                                 nstddev=5,cont.limits=NULL,
                                 std.vars=NULL,diag.file=NULL,use.mvhist=T,checkpointfname="checkpoint.Rda"){

  start.time=proc.time()
  if(is.na(rseed)==FALSE){
    set.seed(rseed)
  }

  ### Processing inputs...
  n.response=length(response.vars.all)

  if((reg.assumption.plots==TRUE)&(is.null(reg.plot.path))){
    stop("No file path for regression diagnostic plots provided in reg.plot.path input.")
  }
  if((reg.assumption.plots==TRUE)&(is.null(reg.plot.props)==TRUE)){
    warning("No regression diagnostic plot properties supplied in reg.plot.props. Using defaults: pt.sz=4,ln.sz=2.5,bs.sz=33, ncol=3,width=3500,height=200")
    reg.plot.props=list("pt.sz"=4,"ln.sz"=2.5,"bs.sz"=33,"ncol"=3,"width"=3500,"height"=2000)
  }

  ## budget and methd parameter inputs....
  if(include.dp.model==TRUE){
    stopifnot((is.null(win.y.resvar)==FALSE)&(is.null(win.y.trcoef)==FALSE))
    stopifnot((win.y.resvar>0)&(win.y.trcoef>0)&((win.y.resvar+win.y.trcoef)<1))
    stopifnot(is.null(bound.means)==FALSE&is.null(bound.sds)==FALSE&is.na(n.iters)==FALSE)
  }


  if(length(synthdata.budget)==1){
    synthdata.budget.eps=synthdata.budget
    synthdata.budget.del=0
  }else{
    synthdata.budget.eps=synthdata.budget[1]
    synthdata.budget.del=synthdata.budget[2]
  }
  #print(paste(synthdata.budget.eps,synthdata.budget.del))

  if(is.null(otherbudgetsDP)==FALSE){
    other.mv.eps=sapply(otherbudgetsDP,"[[",1)
    other.mv.del=sapply(otherbudgetsDP,"[[",2)
    other.per.y.del=sapply(otherbudgetsDP,"[[",4)
    other.per.y.eps=sapply(otherbudgetsDP,"[[",3)
    other.eps=(other.mv.eps+(other.per.y.eps*n.response))
    other.del=(other.mv.del+(other.per.y.del*n.response))
    other.mv.prop=other.mv.eps/(other.eps)
    n.otherb=length(otherbudgetsDP)
    if(is.null(otherbudgetsDP.names)==TRUE){
      if(is.null(names(otherbudgetsDP))==TRUE){
        otherbudgetsDP.names=names(otherbudgetsDP)
        }else{
          possible.names=paste0("DP-Mb_Eps",
                               other.eps,
                               "_Del",other.del)
          dup.names=duplicated(possible.names)
          if(sum(dup.names)>0){
            dup.names=possible.names[dup.names]
            for(nm in dup.names){
              dup.idx=which(possible.names==nm)
              possible.names[dup.idx]=paste0(possible.names[dup.idx],"_MVProp",other.mv.prop[dup.idx])
            }
          }
            otherbudgetsDP.names=possible.names
        }
    }
    all.DP.Mb.budgets=data.frame("Name"=otherbudgetsDP.names,"Epsilon"=other.eps,
                                 "Delta"=other.del,"MV.Prop"=other.mv.prop)
    if(include.dp.model==TRUE){
      all.DP.Mb.budgets=
        rbind(data.frame("Name"="DP-Mb Same Privacy","Epsilon"=synthdata.budget.eps,
                         "Delta"=synthdata.budget.del,"MV.Prop"=mv.prop.budget),
              all.DP.Mb.budgets)

    }
  }else{
    n.otherb=0
    all.DP.b.budgets=NULL
  }

  liberia.sub=liberia.sub[,colnames(liberia.sub)%in%c(treat.vars,block.vars,response.vars.all,covariate.vars)]

  std.vars=colnames(liberia.sub)[colnames(liberia.sub)%in%std.vars]
  cont.as.cont=colnames(liberia.sub)[colnames(liberia.sub)%in%cont.as.cont]
  cont.limits=cont.limits[colnames(liberia.sub)[colnames(liberia.sub)%in%names(cont.limits)]]

  ### Constructing regression model formulas for each response variable
  families=rep("gaussian",length(response.vars.all))

  covariate.formula=paste(c(covariate.vars,block.vars),collapse="+")
  treat.expand.formula=paste(treat.vars,collapse="+")

  reg.formulas=paste(response.vars.all,
                     paste(treat.expand.formula,covariate.formula,sep="+"),
                     sep="~")


  #variables
  covariate.data=liberia.sub[,colnames(liberia.sub)%in%c(covariate.vars,block.vars)]
  #warning(paste0("There are ",ncol(covariate.data)," columns in covariate data.
  #                with ",length(block.vars),"block vars and ",length(covariate.vars),".
  #                 There are ",ncol(liberia.sub[,!(colnames(liberia.sub)%in%c("treatment","control"))])," columns in
  #                 liberia.sub without treatment or control. There are ",length(treat.vars)," treatment variables."))

  if(sum(colnames(liberia.sub)%in%std.vars)==0){
    std.vars=NULL
  }else{
  std.vars=colnames(liberia.sub)[colnames(liberia.sub)%in%std.vars]
  }
  if(sum(colnames(liberia.sub)%in%cont.as.cont)==0){
    cont.as.cont=NULL
    cont.limits=NULL
  }else{
  cont.as.cont=colnames(liberia.sub)[colnames(liberia.sub)%in%cont.as.cont]
  cont.limits=cont.limits[names(cont.limits)%in%cont.as.cont]
  cont.limits=cont.limits[colnames(liberia.sub)[colnames(liberia.sub)%in%names(cont.limits)]]
  }
  #### intialize time keeping ###
  comp.times.perturb.hist=data.frame("SyntheticMethod"=c("MV Histogram","Hybrid","DP-Mb Same Privacy"),
                                     "nHistVars"=c(ncol(liberia.sub[,!(colnames(liberia.sub)%in%c("treatment","control"))]),
                                                   ncol(covariate.data),ncol(covariate.data)),
                                     "SynthTreatmentAssigned"=c(F,T,T),
                                     "FullDP"=c(T,F,T),
                                     "MakeMVHist.Time"=rep(NA,3),
                                     "SanHistProps.Time"=rep(NA,3),
                                     "SampleHist.Time"=rep(NA,3),
                                     "Total.Time"=rep(NA,3))
  comp.times.by.response=data.frame("SyntheticMethod"=c(rep("Hybrid",3),rep("DP-Mb Same Privacy",5)),
                                    "nPredictorVars"=rep(ncol(covariate.data)+length(block.vars)+length(treat.vars),8),
                                    "Timed.Step"=c("total.time","san.response","fit.conf.model",
                                                   "total.time","fit.conf.model","iter.proxy.fit","san.summary","san.response"))
  y.times.df=as.data.frame(matrix(rep(NA,n.response*8),ncol=n.response))
  colnames(y.times.df)=response.vars.all
  comp.times.by.response=cbind(comp.times.by.response,y.times.df)

  if(include.dp.model==FALSE){
    comp.times.perturb.hist=comp.times.perturb.hist[1:3,]
    comp.times.by.response=comp.times.by.response[1:3,]
  }
  if(is.null(otherbudgetsDP)==FALSE){
    comp.times.perturb.hist.other=
      data.frame("SyntheticMethod"=unlist(otherbudgetsDP.names),
                        "nHistVars"=rep(ncol(covariate.data),n.otherb),
                        "SynthTreatmentAssigned"=rep(T,n.otherb),
                        "FullDP"=rep(T,n.otherb),
                        "MakeMVHist.Time"=rep(NA,n.otherb),
                        "SanHistProps.Time"=rep(NA,n.otherb),
                        "SampleHist.Time"=rep(NA,n.otherb),
                        "Total.Time"=rep(NA,n.otherb))
    comp.times.by.response.other=
            data.frame("SyntheticMethod"=rep(unlist(otherbudgetsDP.names),each=5),
                       "nPredictorVars"=rep(ncol(covariate.data)+length(block.vars)+length(treat.vars),5*n.otherb),
                       "Timed.Step"=rep(c("total.time","fit.conf.model","iter.proxy.fit","san.summary","san.response"),n.otherb))
    y.times.df.other=as.data.frame(matrix(rep(NA,n.response*5*n.otherb),ncol=n.response))
    colnames(y.times.df.other)=response.vars.all
    comp.times.by.response.other=cbind(comp.times.by.response.other,y.times.df.other)
    comp.times.by.response=rbind(comp.times.by.response,comp.times.by.response.other)
    comp.times.perturb.hist=rbind(comp.times.perturb.hist,comp.times.perturb.hist.other)
  }
  if(is.na(sim)==FALSE){
    comp.times.perturb.hist$Sim=rep(sim,nrow(comp.times.perturb.hist))
    comp.times.by.response$Sim=rep(sim,nrow(comp.times.by.response))
  }


  get.synths.start=proc.time()
  comp.times.summary=list("Set Up Time"=(get.synths.start-start.time)[[3]])


  print("Starting Synthetic Data Methods")

  ############# Multivariate Histogram Fully Synthetic epsilon=2 ############################
  #get multivariate histogram synthetic data
  if(use.mvhist==T){
  mv.liberia.sub=liberia.sub

  mv.liberia.sub$treatment="control"
  for(tvar in treat.vars){
    mv.liberia.sub$treatment[unlist(mv.liberia.sub[,tvar])==1]=tvar
  }
  #mv.liberia.sub$treatment=as.factor(mv.liberia.sub$treatment)

  mv.out=DPrct::synthdata_perturb_mvhist(data=mv.liberia.sub[,!(colnames(mv.liberia.sub)%in%c(treat.vars,"control"))],
                                         epsilon=synthdata.budget.eps,
                                         delta=synthdata.budget.del,
                                         continuous.vars=cont.as.cont[cont.as.cont%in%colnames(covariate.data)],
                                         cat.vars=c("treatment",block.vars),
                                         num.bin=mv.bins,
                                         add.cont.variation=use.continuous.noise,
                                         treatment.colname=treat.vars,
                                         assign.type="block",
                                         return.time=TRUE,
                                         with.treatment=FALSE, #treatment is kept in the mv histogram
                                         within.blocks=FALSE,
                                         blocks=block.vars,
                                         block.sizes=NULL,
                                         rseed=NA,continuous.limits=cont.limits,
                                         standardize.cont=std.vars)


  for(tvar in c("control",treat.vars)){
    mv.out[,tvar]=ifelse(mv.out$treatment==tvar,1,0)
  }

  mv.synth.comp.time=attr(mv.out,"comp.time") #computation time
  comp.times.perturb.hist[1,5:8]=mv.synth.comp.time
  comp.times.summary=c(comp.times.summary,
                       list("Perturb MV Histogram"=(proc.time()-get.synths.start)[[3]]))

  print("Done with Perturb MV Histogram")
  }
  #stop()
  ################ Hybrid DP epsilon=1 #######################
  all.y.hybrid.start=proc.time()

  if(sum(colnames(covariate.data)%in%std.vars)==0){
    std.vars=NULL
  }else{
    std.vars=colnames(covariate.data)[colnames(covariate.data)%in%std.vars]
  }
  if(sum(colnames(covariate.data)%in%cont.as.cont)==0){
    cont.as.cont=NULL
    cont.limits=NULL
  }else{
    cont.as.cont=colnames(covariate.data)[colnames(covariate.data)%in%cont.as.cont]
    cont.limits=cont.limits[names(cont.limits)%in%cont.as.cont]
    cont.limits=cont.limits[colnames(covariate.data)[colnames(covariate.data)%in%names(cont.limits)]]
  }

  mv.covariates=DPrct::synthdata_perturb_mvhist(data=covariate.data[,!(colnames(covariate.data)%in%c("treatment","control",treat.vars))],
                                                epsilon=synthdata.budget.eps,
                                                delta=synthdata.budget.del,
                                                continuous.vars=cont.as.cont,
                                                num.bin=mv.bins,
                                                add.cont.variation=use.continuous.noise,
                                                treatment.colname="treatment",
                                                assign.type="block",
                                                return.time=TRUE,
                                                with.treatment=TRUE,
                                                within.blocks=FALSE,
                                                blocks=block.vars,
                                                block.sizes=NULL,
                                                conditions=c(treat.vars,"control"),
                                                rseed=NA,
                                                factorial=TRUE,
                                                continuous.limits=cont.limits,
                                                standardize.cont=std.vars)
  #mv.covariates$treatment=ifelse("control"==1,"control",NA) ###NOTE IF SOMETHING BREAKS IT MIGHT BE THIS CHANGE!
  #print(sum(grepl(paste0("control|treatment|",paste0(treat.vars,collapse="|")),colnames(mv.covariates))))
  #print(head(mv.covariates[,colnames(mv.covariates)%in% c("control",treat.vars,"treatment")]))

  mv.covariates$treatment="control"
  for(tr.var in c("control",treat.vars)){
  mv.covariates$treatment[unlist(mv.covariates[,tr.var])==1]=tr.var
  }
  #print(head(mv.covariates[,colnames(mv.covariates)%in% c("control",treat.vars,"treatment")]))

  #liberia.sub[,treat.vars]=apply(liberia.sub[,treat.vars],2,as.factor)
  #mv.covariates=mv.covariates[,colnames(mv.covariates)!=c("control")]
  mv.covariate.comp.time=attr(mv.covariates,"comp.time")




  print("Done with MV Histogram for Hybrid")


  get_synth_response_hybriddp=function(idx,synth.data=mv.covariates[,colnames(mv.covariates)!="treatment"],
                                       formulas=reg.formulas,
                                       y.vars=response.vars.all,
                                       fams=families,
                                       confidential.data=liberia.sub,
                                       tr.vars=treat.vars){
    temp.out=DPrct::hybrid_synth(formulas[idx],
                                 family=fams[idx],
                                 confidential.data=confidential.data,
                                 epsilon=NA,
                                 delta=0,
                                 treatment.var=tr.vars,
                                 synth.data=synth.data,
                                 continuous.vars=NULL,
                                 num.bin=NULL,
                                 bin.param=NA,
                                 add.cont.variation=FALSE,
                                 assign.type="simple",
                                 blocks=NULL,
                                 within.blocks=FALSE,
                                 clusters=NULL,
                                 returntypes=c("synth.data","comp.time"),
                                 rseed=NA)
    response.var=y.vars[idx]
    temp.synth=temp.out[[1]]
    list(temp.synth[,response.var],temp.out[[2]])
  }

  synth.response=parallel::mclapply(seq(1,n.response),get_synth_response_hybriddp)
  hybrid.dp.comp.times=sapply(synth.response,"[[",2)
  colnames(hybrid.dp.comp.times)=paste0(response.vars.all,".comp.time")
  synth.ys=dplyr::bind_cols(lapply(synth.response,"[[",1))
  hybrid.synthdata=cbind(synth.ys,mv.covariates)
  #print(sum(grepl(paste0("control|treatment|",paste0(treat.vars,collapse="|")),colnames(hybrid.synthdata))))
  #print(head(hybrid.synthdata[,colnames(hybrid.synthdata)%in% c("control",treat.vars,"treatment")]))

  #hybrid.synthdata$control=ifelse(hybrid.synthdata$treatment=="control",1,0)
  all.y.hybrid.comp.time=(proc.time()-all.y.hybrid.start)[[3]]

  comp.times.perturb.hist[2,5:8]=mv.covariate.comp.time
  comp.times.summary=c(comp.times.summary,
                       list("Hybrid"=all.y.hybrid.comp.time))
  for(nc in seq(1,n.response)){
    comp.times.by.response[1:3,3+nc]=unlist(hybrid.dp.comp.times[,nc])
  }

  hybrid.dp.comp.times=list("total.hybrid.time"=all.y.hybrid.comp.time,
                            "hybrid.response.times"=hybrid.dp.comp.times,
                            "mv.T.X.time"=mv.covariate.comp.time)

  print("Done with Not-DP Hybrid Method")

  ####################

  ################ My DP Synthetic epsilon=1 #######################
  dp_model_based_full=function(liberia.sub,
                               covariate.data,
                               reg.formulas,families,n.response,
                               mv.priv.budget,per.response.priv.budget,
                               cont.as.cont,block.vars,treat.vars,response.vars.all,
                               mv.bins, use.continuous.noise,
                               win.y.resvar,win.y.trcoef,
                               bound.sds,bound.means,n.iters,range.alpha,
                               method.name="DP-Mb Same Privacy",use.se.normal=se.normal,
                               continuous.limits=cont.limits,
                               standardize.cont=std.vars){
    synthdata.start=proc.time()





    ################### Get Parameters ##############################
    ##each of the responses get equal proportion of the the overall budget after the sythetic covariate cost
    ## within each response, the residual variance gets 10% of the budget
    ##      the other coefficients (and the intercept) get the remaining portion of budget within each response
    win.y.x.coefs=(1-win.y.resvar-win.y.trcoef)
    ############

    ##############################
    n.x.coefs=length(covariate.vars)+1+length(c(t(unique(liberia.sub[,block.vars[1]])),
                                                t(unique(liberia.sub[,block.vars[2]]))))
    n.treat.coefs=length(treat.vars)
    win.y.budget.props=c(rep(win.y.trcoef/n.treat.coefs,n.treat.coefs),
                         rep(win.y.x.coefs/n.x.coefs,n.x.coefs),win.y.resvar)


      per.y.eps=per.response.priv.budget[1]
      per.y.del=per.response.priv.budget[2]


    #get synthetic covariate and treatment data
    #covariate.data from Hybrid DP section
    mv.covariates2=DPrct::synthdata_perturb_mvhist(data=covariate.data[,!(colnames(covariate.data)%in%c(treat.vars,"control","treatment"))],
                                                   epsilon=mv.priv.budget[1],
                                                   delta=mv.priv.budget[2],
                                                   continuous.vars=cont.as.cont,
                                                   num.bin=mv.bins,
                                                   add.cont.variation=use.continuous.noise,
                                                   treatment.colname="treatment",
                                                   assign.type="block",
                                                   return.time=TRUE,
                                                   rseed=NA,
                                                   with.treatment=TRUE,
                                                   within.blocks=FALSE,
                                                   blocks=block.vars,
                                                   block.sizes=NULL,
                                                   conditions=c(treat.vars,"control"),
                                                   factorial=TRUE,
                                                   continuous.limits=continuous.limits,
                                                   standardize.cont=standardize.cont)
    # for(tr.var in c(treat.vars,"control")){
    #   mv.covariates2[,tr.var]=ifelse(mv.covariates2$treatment==tr.var,1,0)
    # }
    #mv.covariates2$treatment=ifelse("control"==1,"control",NA) ###CHANGED THIS

    mv.covariates2$treatment="control"
    for(tr.var in c("control",treat.vars)){
      mv.covariates2$treatment[unlist(mv.covariates2[,tr.var])==1]=tr.var
    }
  #print(head(mv.covariates2[,colnames(mv.covariates2)%in% c("control",treat.vars,"treatment")]))


    #mv.covariates2=mv.covariates2[,!(colnames(mv.covariates2)%in%c("control"))]
    mv.covariate.comp.time2=attr(mv.covariates2,"comp.time")



    print("Done with MV Hist of Full DP")
    get_synth_response_synthdatadp=function(idx,
                                            synth.data=mv.covariates2,
                                            formulas=reg.formulas,
                                            y.vars=response.vars.all,
                                            fams=families,
                                            confidential.data=liberia.sub,
                                            tr.vars=treat.vars,
                                            per.response.eps=per.y.eps,
                                            per.response.del=per.y.del,
                                            bd.sd=bound.sds,
                                            bd.mean=bound.means,
                                            nits=n.iters,
                                            budget.prop.list=win.y.budget.props,
                                            alph=range.alpha,se.norm=use.se.normal){
      temp.out=DPrct::dp_synthdata(formula=formulas[idx],
                                   confidential.data=confidential.data,
                                   synth.data=synth.data,
                                   synth.epsilon=NULL,
                                   synth.delta=NULL,
                                   continuous.vars=NULL,
                                   num.bin=NULL,
                                   bin.param=NA,
                                   assign.type="simple",
                                   blocks=NULL,
                                   clusters=NULL,
                                   within.blocks=FALSE,
                                   epsilon.list=as.list(per.response.eps*budget.prop.list),
                                   delta.list=as.list(per.response.del*budget.prop.list),
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
                                   treat.ppart.list=NULL,
                                   se.normal=T)

      response.var=y.vars[idx]
      temp.synth=temp.out[[1]]
      list(temp.synth[,response.var],temp.out[[2]])
    }

    synth.response=parallel::mclapply(seq(1,n.response),get_synth_response_synthdatadp)
    #synth.response=lapply(seq(1,n.response),get_synth_response_synthdatadp)
    #warning("outside synth.response mclappy")
    full.dp.synth=dplyr::bind_cols(mv.covariates2,dplyr::bind_cols(sapply(synth.response,"[[",1)))
    # for(tr.var in c(treat.vars,"control")){
    #   mv.covariates2[,tr.var]=ifelse(mv.covariates2$treatment==tr.var,1,0)
    # }
    #full.dp.synth$treatment=ifelse(full.dp.synth$tpassonly==1,"tpassonly","control")
    #full.dp.synth$treatment[full.dp.synth$cashassonly==1]="cashassonly"
    #full.dp.synth$treatment[full.dp.synth$tpcashass==1]="tpcashass"
    #full.dp.synth$control=ifelse(full.dp.synth$treatment=="control",1,0)
    by.y.comp.times=sapply(synth.response,"[[",2)
    full.dp.comp.times=list("total.dp.time"=(proc.time()-synthdata.start)[[3]],
                            "dp.response.time"=by.y.comp.times,
                            "mv.T.X.time"=mv.covariate.comp.time2)
    comp.times.perturb.hist.full=mv.covariate.comp.time2#[3,5:8]=mv.covariate.comp.time2
    comp.times.summary.full=list((proc.time()-synthdata.start)[[3]])
    names(comp.times.summary.full)=method.name
    comp.times.y.full=unlist(by.y.comp.times[,1])
    if(n.response>1){
      for(nc in seq(2,n.response)){
        comp.times.y.full=cbind(comp.times.y.full,unlist(by.y.comp.times[,nc]))#[4:8,3+nc]=unlist(by.y.comp.times[,nc])
      }
    }
    if(is.null(per.y.del)==TRUE){
      per.y.del=0
    }
    full.budget=list("epsilon"=mv.priv.budget[1]+sum(per.y.eps),"delta"=mv.priv.budget[2]+sum(per.y.del))
    return(list(full.dp.synth,
                comp.times.summary.full,
                comp.times.perturb.hist.full,
                comp.times.y.full,
                full.budget,
                synth.response))
  }

  all.DP.Mb.df.list=NULL
  all.DP.Mb.comp.times.summary=NULL
  DP.MB.df.names=NULL


  if(include.dp.model==TRUE){
    mv.priv.budget=unlist(c(synthdata.budget.eps,synthdata.budget.del)*mv.prop.budget)
    per.response.priv.budget=c(synthdata.budget.eps,synthdata.budget.del)*(1-mv.prop.budget)/n.response
      fulldpout=dp_model_based_full(liberia.sub=liberia.sub[,!(colnames(liberia.sub)%in%c("control","treatment"))],
                                    covariate.data=covariate.data,
                                    reg.formulas=reg.formulas,families=families,
                                    n.response=n.response,
                                    mv.priv.budget=mv.priv.budget,
                                    per.response.priv.budget=per.response.priv.budget,
                                    cont.as.cont=cont.as.cont,block.vars=block.vars,
                                    treat.vars=treat.vars,response.vars.all=response.vars.all,
                                    mv.bins=mv.bins, use.continuous.noise=use.continuous.noise,
                                    win.y.resvar=win.y.resvar,win.y.trcoef=win.y.trcoef,
                                    bound.sds=bound.sds,bound.means=bound.means,
                                    n.iters=n.iters,range.alpha=range.alpha,
                                    method.name="DP-Mb",use.se.normal=se.normal,continuous.limits=cont.limits,
                                    standardize.cont=std.vars)
      #print("Done with MV Hist of Full DP before post process")
      all.DP.Mb.df.list=c(all.DP.Mb.df.list,list("DP-Mb Same Privacy"=fulldpout[[1]]))
      all.DP.Mb.comp.times.summary=c(all.DP.Mb.comp.times.summary,fulldpout[[2]])
      comp.times.perturb.hist[comp.times.perturb.hist$SyntheticMethod=="DP-Mb Same Privacy",5:8]=fulldpout[[3]]
      comp.times.by.response[comp.times.by.response$SyntheticMethod=="DP-Mb Same Privacy",
                             seq(4,3+n.response)]=fulldpout[[4]]
      DP.MB.df.names=c(DP.MB.df.names,"DP-Mb Same Privacy")
    print("Done with DP Model-Based")
  }else{
    print("Skipped DP Model-Based")
  }
  if(is.null(otherbudgetsDP)==FALSE){
    for(i in seq(1,n.otherb)){
      name.i=otherbudgetsDP.names[i]
      print(paste("name:",name.i,". budget mv:",other.mv.eps[i],other.mv.del[i],
                  ", per y:",other.per.y.eps[i],other.per.y.del[i]))
      #warning(paste0("in otherbudgetsDP ",name.i))
        fulldpout=dp_model_based_full(liberia.sub=liberia.sub,covariate.data=covariate.data,
                                    reg.formulas=reg.formulas,families=families,
                                    n.response=n.response,
                                    mv.priv.budget=c(other.mv.eps[i],other.mv.del[i]),
                                    per.response.priv.budget=c(other.per.y.eps[i],other.per.y.del[i]),
                                    cont.as.cont=cont.as.cont,block.vars=block.vars,
                                    treat.vars=treat.vars,response.vars.all=response.vars.all,
                                    mv.bins=mv.bins, use.continuous.noise=use.continuous.noise,
                                    win.y.resvar=win.y.resvar,win.y.trcoef=win.y.trcoef,
                                    bound.sds=bound.sds,bound.means=bound.means,
                                    n.iters=n.iters,range.alpha=range.alpha,
                                    method.name=name.i,
                                    use.se.normal = se.normal,
                                    continuous.limits=cont.limits,
                                    standardize.cont=std.vars)
        print("Extra DP-Mb out before pre processing")
      all.DP.Mb.df.list=c(all.DP.Mb.df.list,list(fulldpout[[1]]))
      DP.MB.df.names=c(DP.MB.df.names,name.i)
      all.DP.Mb.comp.times.summary=c(all.DP.Mb.comp.times.summary,fulldpout[[2]])
      comp.times.perturb.hist[comp.times.perturb.hist$SyntheticMethod==name.i,5:8]=fulldpout[[3]]
      comp.times.by.response[comp.times.by.response$SyntheticMethod==name.i,
                             seq(4,3+n.response)]=fulldpout[[4]]
    }
    names(all.DP.Mb.df.list)=DP.MB.df.names
    print("Done with Other Budget DP Model-Based")
  }else{
    print("No Other Budget DP Model-Based")
  }


  ##############
  comp.times.summary=c(comp.times.summary,all.DP.Mb.comp.times.summary,
                       list("All Synthetic Data Methods"=(proc.time()-get.synths.start)[[3]]))

  table.start=proc.time()
  sterr.type="HC1"
  stderr.func=function(mod,std.type=sterr.type){
    sqrt(diag(sandwich::vcovHC(mod,type=std.type)))
  }



  table.data.list=list("Confidential"=liberia.sub,
                    "MV Histogram"=mv.out,
                    "Hybrid"=hybrid.synthdata)

  if(is.null(all.DP.Mb.df.list)==FALSE){
    table.data.list=c(table.data.list,all.DP.Mb.df.list)
  }

  data.list=list()
  listname=names(table.data.list)
  i=1
  for(df in table.data.list){
    print(listname[i])
    i=i+1
    if(sum(colnames(df)=="treatment")==0){
      df$treatment="control"
      for(tvar in treat.vars){
       df$treatment[unlist(df[,tvar])==1]=tvar
      }
    }else if(sum(is.na(df$treatment))>0){
      df$treatment="control"
      for(tvar in treat.vars){
        df$treatment[unlist(df[,tvar])==1]=tvar
      }
    }
    df$control=ifelse(df$treatment=="control",1,0)
    print(paste("treatment column has:",paste(unique(df$treatment),collapse=", "),". There are ",sum(grepl(paste0("control|",paste0(treat.vars,collapse="|")),colnames(df)))," columns of control or treatment levels."))

    print(head(df[,colnames(df)%in%c("treatment","control",treat.vars)]))
    df[,block.vars]=lapply(df[,block.vars],as.factor)
    data.list=c(data.list,list(df))
  }
  names(data.list)=listname
  save(data.list,comp.times.summary,comp.times.by.response,file=checkpointfname)
  print("Before ITTcompare")


# ))
#

  out.df=DPrct::compare_ITTdata(data.list,
                                reg.formulas,
                                response.vars=response.vars.all,
                                families=rep("gaussian",length(reg.formulas)),
                                treat.vars=treat.vars,
                                control.var="control",
                                bonferroni.npvals=bonferroni.npvals,
                                mult.test.correct=mult.test.correct,
                                add.pval.stars=TRUE,
                                model.names=mod.names,
                                stderr.func=stderr.func,
                                treat.names=treat.names,
                                include.ci.overlap = TRUE,
                                ci.overlap.ref.name = "Confidential",
                                ci.confidence=ci.confidence,
                                pivot.xtreat=TRUE)
  if(is.na(sim)==FALSE){
    out.df$Sim=sim
  }

  end.time=proc.time()
  comp.times.summary=c(comp.times.summary,list("Make Table"=(end.time-table.start)[[3]],
                                               "Total Time"=(end.time-start.time)[[3]])
  )
  attr(out.df,"all.comp.times")=list("Summary"=comp.times.summary,
                                     "Perturb.Hist.Steps"=comp.times.perturb.hist,
                                     "By.Response.Steps"=comp.times.by.response)
  out.list=list(out.df,comp.times.summary,comp.times.perturb.hist,comp.times.by.response)
  save(out.list,file=fname)
  #return(list(table.data.list,comp.times.summary,comp.times.perturb.hist,comp.times.by.response))
#}

  ### Plots and Summaries

  if(reg.assumption.plots==TRUE){
  for(i in seq(1,length(table.data.list))){
    plt.name=tolower(unlist(gsub("[^A-z0-9]","",names(table.data.list)[i])))
    plt.df=table.data.list[[i]]
    plt.df=plt.df[,colnames(plt.df)%in%c(response.vars.all,covariate.vars,block.vars,treat.vars)]

    san.plots= lapply(seq(1,n.response),
                       function(idx)reg_assumptions(reg.formulas[idx],
                                                    data=plt.df,
                                                    response.var=response.vars.all[idx],
                                                    pt.sz=reg.plot.props$pt.sz,
                                                    ln.sz=reg.plot.props$ln.sz,
                                                    bs.sz=reg.plot.props$bs.sz,
                                                    long.plot.title=F,
                                                    mod.name=response.vars.all[idx]))
    save.plot.name=paste0(reg.plot.path,"/",reg.plot.prefix,"_",plt.name,".png")
    ggplot2::ggsave(filename=save.plot.name,
        width=reg.plot.props$width,height=reg.plot.props$height,
        plot=cowplot::plot_grid(plotlist=san.plots,ncol=reg.plot.props$ncol),
        bg="white")
    #dev.off()

    save(san.plots,table.data.list,response.vars.all,reg.formulas,
         file=paste0(reg.plot.path,"/",reg.plot.prefix,"_",plt.name,".Rda"))
    print(paste("Done Plotting Regression Diagnostic Plots for",plt.name))

  }
  }

  return(out.list)
}
