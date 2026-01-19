##internal function adapt dp_iter_hybrid
nonoise_iter_hybrid=function(conf.model,
                     confidential.data,
                     num.iters,
                     predictor.vars,
                     response.var){
  yb=stats::simulate(conf.model,nsim=num.iters,type="response") #simulate proxy responses

  #########
  #FUNCTION:attach proxy response to predictor data, fit model, return coefficients
  get_betas=function(coly,
                     conf.predictors=confidential.data[,predictor.vars],
                     conf.mod=conf.model,
                     y.var=response.var
  ){
    iter.data=cbind(conf.predictors,coly)
    colnames(iter.data)=c(colnames(conf.predictors),y.var)
    modb=stats::glm(formula=conf.mod$formula,family = conf.mod$family,data=iter.data)
    return(modb$coefficients)
  }
  #######

  betas=t(sapply(1:num.iters,function(x)get_betas(yb[,x])))
  #use get_betas for each iteration.
  return(betas[,!is.na(conf.model$coefficients)])
}



nonoise_iter_hybrid_mse=function(conf.model,
                        confidential.data,
                        synth.data,
                        num.iters,
                        model.vars,
                        response.var){

  pred.data=synth.data[,model.vars[-1]]

  cat.cols=sapply(seq(1,ncol(pred.data)),function(x)is.character(pred.data[,x]))
  pred.data[,cat.cols]=lapply(pred.data[,cat.cols],factor)



  #yb=DPrct::simulate_response_glm(conf.model,newdata=pred.data,nsim=num.iters)

  yb=DPrct::simulate_response_glm(conf.model,newdata=pred.data,nsim=num.iters,
                              predictor.formula = as.formula(paste0("~",paste0(all.vars(as.formula(conf.model$formula))[-1],collapse="+"))))


  predictor.formula=stats::as.formula(
    paste0("~",as.character(
      as.formula(conf.model$formula))[3]))
  pred.newdata=data.frame(synth.data[,model.vars[-1]])
  colnames(pred.newdata)=model.vars[-1]
  modMat=stats::model.matrix(predictor.formula, pred.newdata)

  #get estimated model coefficients
  mod.coefs=stats::coefficients(conf.model)
  has.coef.name=names(mod.coefs)[!is.na(mod.coefs)] #only columns where coefficient!=NA
  num.coefs=length(has.coef.name)

  modMat=modMat[,colnames(modMat)%in%has.coef.name]
  #########
  #FUNCTION:attach proxy response to predictor data, fit model, return coefficients
  get_betas_and_residuals=function(coly,
                                   synth.predictors=synth.data[,model.vars[-1]],
                                   conf.mod=conf.model,
                                   y.var=response.var,
                                   return.invxtx=FALSE){

    #bind response to predictors
    iter.data=cbind(synth.predictors,coly)
    colnames(iter.data)=c(colnames(synth.predictors),y.var)
    #fit model
    modb=stats::glm(formula=conf.mod$formula,family = conf.mod$family,data=iter.data)
    if(as.numeric(return.invxtx)==1){
      inv.xtx=stats::vcov(modb)*summary(modb)$df.residuals/sum(modb$residuals^2)
    }else{
      inv.xtx=NULL
    }
    #return coefficient estimates and residuals
    return(list("betas"=modb$coefficients,"residuals"=modb$residuals,"inv.xtx"=inv.xtx))
  }
  #######

  #parallelize apply to each iteration
  iter.out=parallel::mclapply(1:num.iters,
                              function(x)get_betas_and_residuals(yb[,x],return.invxtx = x))
  #combine coefficients into matrix
  #column for each coefficient (including intercept), row for each iteration
  betas=t(sapply(1:num.iters,function(idx)iter.out[[idx]]$betas))

  #stack residuals into a vector
  residuals=unlist(lapply(1:num.iters,function(idx)iter.out[[idx]]$residuals))

  #get sanitized estimate of residual standard deviation
  san.mse=mean((residuals)^2)
  nr=nrow(confidential.data) #number of observations

  if(sum(colSums(modMat==0)==nr)>0){
    warning("some columns in modMat have only 0 values")
  }

  #covariance matrix of coefficients is sigma^2
  inv.xtx=iter.out[[1]]$inv.xtx

  cov.mat=(san.mse/(nr-num.coefs))*nr*(inv.xtx)#t(modMat)%*%modMat/nr)

  return(list("iter.betas"=betas[,colnames(betas)%in%has.coef.name],"cov.mat"=cov.mat,"san.mse"=san.mse))
}


nonoise_dpmb=function(formula,
                      confidential.data,
                      synth.data,
                      continuous.vars=NULL,
                      num.bin=NULL,
                      bin.param=NA,
                      num.iters=1e4,
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
                      use.san.residerror=TRUE){
  start.time=proc.time()
  time.values=NULL
  #print("Inside dp_synthdata")

  if(is.na(rseed)==FALSE){
    set.seed(rseed)
  }

  psf=DPrct::parse_formula(formula=formula,confidential.data = confidential.data,
                    treatment.var = treatment.var)
  formula=psf[[1]]
  model.vars=psf[[2]]
  response.var=psf[[3]]
  synth.vars=psf[[4]]
  treatment.var=psf[[5]]

  synth.vars=c(model.vars[-1],blocks,clusters)
  synth.vars=synth.vars[!duplicated(synth.vars)]


  confidential.data=confidential.data[,colnames(confidential.data)%in%c(response.var,synth.vars)]

  model.fit.start=proc.time()
  #### Step 0: In paper
  #fit the confidential data to the model and get estimated coefficient and standard error
  conf.model<-stats::glm(formula=formula,family="gaussian",data=confidential.data)
  ####
  if(return.time==TRUE){
    time.values=c(time.values,
                  list("conf.model.fit"=(proc.time()-model.fit.start)[[3]]))
  }


  # if synthetic data is supplied, check it is data.frame w/ all the model.vars
  num.coefs=sum(!is.na(conf.model$coefficients))
  #get model matrix
  predictor.formula=stats::as.formula(paste0("~",as.character(formula)[3]))
  pred.newdata=as.data.frame(synth.data[,model.vars[-1]])
  colnames(pred.newdata)=model.vars[-1]
  cat.cols=sapply(seq(1,ncol(pred.newdata)),function(x)is.character(pred.newdata[,x]))
  pred.newdata[,cat.cols]=lapply(pred.newdata[,cat.cols],as.factor)
  modMat=stats::model.matrix(predictor.formula, pred.newdata)




  iter.proxy.fit.start=proc.time()
  if(use.san.residerror==TRUE){

    list2env(nonoise_iter_hybrid_mse(conf.model=conf.model,
                            confidential.data=confidential.data,
                            synth.data=pred.newdata,
                            num.iters=num.iters,
                            model.vars=model.vars,
                            response.var=response.var),envir=globalenv())


    coef.names=colnames(iter.betas)
  }else{#don't use residerror
    iter.betas=nonoise_iter_hybrid(conf.model=conf.model,
                           confidential.data=confidential.data,
                           num.iters=num.iters,
                           predictor.vars=model.vars[-1],
                           response.var=model.vars[1])
    cov.mat=NULL
    coef.names=colnames(iter.betas)
  }

  if(return.time==TRUE){
    iter.proxy.fit.stop=proc.time()
    time.values=c(time.values,
                  list("iter.proxy.fit.time"=(iter.proxy.fit.stop-iter.proxy.fit.start)[[3]]))
  }


  num.coefs=ncol(iter.betas)
  if(is.null(cov.mat)==TRUE){ #if not using sanitized residuals

    cov.mat=stats::vcov(conf.model)
  }else{
    san.mse=(residuals(conf.model))^2/nrow(confidential.data)
  }
  san.summary=matrix(nrow=num.coefs,ncol=7)
  #get estimate and standard error
  san.summary[,1]=rowMeans(t(iter.betas))
  san.summary[,2]=sapply(seq(1,num.coefs),function(i)base::sqrt(stats::var(iter.betas[,i],na.rm=T)))
  san.summary[,6]=abs(san.summary[,1])/san.summary[,2]
  san.summary[,7]=pnorm(san.summary[,5],lower.tail=FALSE)
  zalph2=qnorm(min(ci.alpha/2,(1-ci.alpha)/2),lower.tail = F)
  san.summary[,4]=san.summary[,1]+(zalph2*san.summary[,2])
  san.summary[,3]=san.summary[,1]-(zalph2*san.summary[,2])
  san.summary[,5]=1-min(ci.alpha/2,(1-ci.alpha)/2)
  rownames(san.summary)=colnames(iter.betas)
  colnames(san.summary)=c("san.coef","san.sd","san.CI.lower","san.CI.upper",
                          "confidence.level","san.z.value","san.pval")

  if(return.time==TRUE){
    san.summ.stop=proc.time()
    time.values=c(time.values,
                  list("san.summary.time"=(san.summ.stop-iter.proxy.fit.stop)[[3]]))
  }

  san.coefs=san.summary[,1]
  has.coef.name=names(san.coefs)[!is.na(san.coefs)] #only columns where coefficient!=NA
  san.coef.std.err=san.summary[,2]

  #warning("before model.matrix")

  #model predictor matrix
  modMat=stats::model.matrix(predictor.formula, pred.newdata)
  modMat.cnames=colnames(modMat)
  not.in.modMat=has.coef.name[!(has.coef.name %in% modMat.cnames)]
  if(length(not.in.modMat)>0){ #if some variables don't appear
    zeros.mat=matrix(0,ncol=length(not.in.modMat),nrow=nrow(modMat))
    modMat=cbind(modMat,zeros.mat)
    colnames(modMat)=c(modMat.cnames,not.in.modMat)
  }
  modMat=modMat[,has.coef.name]

  ############################## PROBLEM HERE B/C MISSING A LEVEL OF BLOCK.
  #coef.rv=mvtnorm::rmvnorm(1, san.coefs[names(san.coefs)%in%has.coef.name],
  #                         cov.mat[rownames(cov.mat)%in%has.coef.name,colnames(cov.mat)%in%has.coef.name])

  sim.res=mvtnorm::rmvnorm(1, rep(0,nrow(modMat)), (san.mse*diag(nrow=nrow(modMat))))
  sim.response=as.numeric(c(modMat%*%matrix(san.coefs[names(san.coefs)%in% has.coef.name],ncol=1))+c(sim.res))

  if(return.time==TRUE){
    san.y.stop=proc.time()
    time.values=c(time.values,
                  list("san.response.time"=(san.y.stop-san.summ.stop)[[3]]))
  }

  #add synthetic response and reorder columns to match confidential data
  curr.synth.colnames=colnames(synth.data)
  conf.colnames=colnames(confidential.data)
  synth.data$response=sim.response
  colnames(synth.data)=c(curr.synth.colnames,response.var)
  synth.data=synth.data[,conf.colnames[conf.colnames %in% colnames(synth.data)]]
  output=list("synth.data"=synth.data)

  if(return.confidential.table==TRUE){
    output=c(output,list("confidential.summary"=summary(conf.model)))
  }
  if(return.san.summary==TRUE){
    san.summary=san.summary[,apply(san.summary,2,function(col)!all(is.na(col)))]
    output=c(output,list("san.summary"=san.summary))
  }
  if(return.time==TRUE){
    stop.time=proc.time()
    output=c(output,
             list("comp.time"=c(list("total.dp_synthdata.time"=(stop.time-start.time)[[3]]),
                                time.values)))
  }
  return(output)
}
