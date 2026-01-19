### Small functions for simulations

sim_coefs<-function(ncat,ncat.groups=2,ncont,
                    cat.coef.start,cat.coef.ratio,
                    cont.coef.start,cont.coef.ratio,
                    treat.effect,intercept){
#function to generate a geometric sequence a_n=start(ratio^n) for n=0,...,end
geomSeq <- function(start,ratio,end){
  end=end-1
  start*ratio**(0:end)
}
if(ncat>0){
  stopifnot(length(ncat.groups)==1|length(ncat.groups)==ncat)
    num.cat.coefs=sum(ncat*ncat.groups)-ncat
    cat.coef=geomSeq(cat.coef.start,cat.coef.ratio,num.cat.coefs) #model coefficient
}else{
  cat.coef=NULL
}
if(ncont>0){
  cont.coef=geomSeq(cont.coef.start,cont.coef.ratio,ncont) #mod.coefs
}else{
  cont.coef=NULL
}
mod.coefs=c(intercept,treat.effect,cont.coef,cat.coef)
return(mod.coefs)
}

sim2_models_list<-function(response,treat.vars,combos.NCov,
                           add.formulas=NULL,add.names=NULL){
  ncat=max(combos.NCov$NCat)
  ncont=max(combos.NCov$NCont)
  if(length(treat.vars)>1){
    base.formula=paste0(response,"~",paste0(treat.vars,collapse="+"))
  }else{
    base.formula=paste0(response,"~",treat.vars)
  }

  reg.formulas.sim2=sapply(seq(1,nrow(combos.NCov)),
                           function(idx)paste0(base.formula,"+",
                                               paste0("x",seq(1,combos.NCov$NCont[idx]),collapse="+"),
                                               "+",
                                               paste0("x",seq(ncont+1,ncont+combos.NCov$NCat[idx]),collapse="+")))
  names(reg.formulas.sim2)=combos.NCov$Name
  if(is.null(add.formulas)==FALSE){
    reg.formulas.sim2=c(reg.formulas.sim2,add.formulas)
    if(is.null(add.names)==TRUE){
      if(is.null(names(add.formulas))==FALSE){
        names(reg.formulas.sim2)=c(combos.NCov$Name,names(add.formulas))
      }else{
      names(reg.formulas.sim2)=c(combos.NCov$Name,paste0("Extra Model ",seq(1,length(add.formulas))))
      }
    }else{
      names(reg.formulas.sim2)=c(combos.NCov$Name,add.names)
    }
  }

  return(reg.formulas.sim2)
}

#function to generate and round n observations from a uniform distribution to d number of digits
runif.round=function(n,min,max,d){
  round(runif(n=n,min=min,max=max),digits=d)
}

#function to get summary values from confidential data
get_conf_treff<-function(idx,data,
                         formulas=reg.models,
                         tr.vars="t1",
                         mod.nms,
                         ci.alpha=0.05,
                         true.vals,
                         stderr.func=NULL,
                         nobs,mse.out=T){
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
  if(mse.out==T){
    mse=(residuals(mod.fit))^2/nobs
  }else{
    mse=NULL
  }
  return(list(tab.out,(proc.time()-start.time.internal)[[3]],mse))
}
