## Generate a simulated dataset
## INPUTS:
### num.cat, num.cont are the number of continuous and categorical variables. Defaults to 0.
### num.treat is the number of treatment and control groups
### nobs is the number of observations to be generated
### mod.coefs are the coefficients for the model to be used in generating the response variable.
### residual.sd is the residual standard deviation of the model to be used in generating the response variables
### cat.counts is a list of vectors where each list entry corresponds to a categorical variable, the vector
###     is the number of observations in each category of the variable. If the list entry is NULL or
###     cat.counts==NULL,then cat.props will be used. If there is only one element and cat.probs==NULL,
###     then these counts will be used for all categorical variables.
### cat.probs is a list of the vectors of probabilities where each list entry corresponds to a categorical
###     variable, the vector is the probability a given observation is in each category of the variable.
###     The entries will be used sequentially to fill in gaps of cat.counts if cat.counts has NULL entries.
###     If there is only one element and cat.counts==NULL, this probabilities will be used for all variables.
### cont.funcs is a list of functions (defaults to stats::rnorm). Each list entry corresponds to the function
###     to generate the continuous variable. Additional arguments for the functions can be included in cont.params
### cont.params is a list of lists. Each first level list corresponds to a continuous variable. The sublist is
###     any additional arguments to be used in the corresponding cont.funcs
### rseed if rseed!=NA, then it is a value to set the random seed. Otherwise, no seed will be set.
generate_simulated_dataset<-function(num.cat=0,num.cont=0,num.treat=2,nobs,mod.coefs,residual.sd,
                                     cat.counts=NULL,cat.probs=NULL,cont.funcs=NULL,cont.params=NULL,reg.mods=NULL,rseed=NA,...){
  if(is.na(rseed)==FALSE){
    set.seed(rseed)
  }
  stopifnot(num.cat>=0&num.cont>=0&num.treat>=0&nobs>=1&residual.sd>0)
  if(num.cat+num.cont==0){
    message("No continuous or covariate columns generated. Only the treatment and response.")
  }
  if(num.treat==1){
    warning("num.treat should include the control group. num.treat=1 is adjusted to be num.treat=2, to include a control group.")
    num.treat=2
  }

  #generate continuous columns
  if(num.cont>0){
    if(is.null(cont.funcs)==TRUE){ #default to normal(0,1) generation
      stopifnot(length(cont.funcs)==length(cont.params))
      cont.funcs=rep(list(stats::rnorm),num.cont)
      cont.params=rep(list(list("n"=nobs)),num.cont)
      nfuncs=num.cont
      warning("No functions provided for continuous covariates. Using normal distribution.")
    }else{
      if(is.list(cont.funcs)==FALSE){
        cont.funcs=list(cont.funcs)
      }
      if(is.list(cont.params)==FALSE){
        cont.params=list(cont.params)
      }
      stopifnot(length(cont.funcs)==length(cont.params))
      nfuncs=length(cont.funcs)
      if(nfuncs<num.cont){
        all.rep=num.cont%/%nfuncs
        num.remain=num.cont%%nfuncs
        rep.times=rep(all.rep,nfuncs)
        if(num.remain>0){
          rep.times[seq(1,num.remain)]=all.rep+1
        }
        cont.funcs=rep(cont.funcs,times=rep.times)
        cont.params=rep(cont.params,times=rep.times)
      }
    }

    cont.mat=sapply(seq(1,num.cont),
                    function(i)base::do.call(cont.funcs[[i]],cont.params[[i]]))
    if((num.cont>2)&(nfuncs<num.cont)){
      divisor=num.cont%/%nfuncs
      remain=num.cont%%nfuncs
      indices=seq(1,num.cont)
      if(remain>0){
      #get sets of indices which have all function types.
      #example case: 3 functions, 8 continuous covariates
      # if A1 A2 A3 B1 B2 B3 C1 C2 are columns of 3 function A, 3 function B, and 2 function C
      # then full.sets.indics=list(c(1,4,7),c(2,5,8))
      full.sets.indices=lapply(seq(1,divisor),
                      function(x)c(indices[seq(x,x+(divisor*remain),by=1+(divisor))],
                                   indices[seq(x+(divisor*remain)+remain,num.cont,by=divisor)]))
      #remain indices which do not have a full set of functions
      #example case: continues from above; remain.sets.indices=list(c(3,6))
      #then reordered cont.mat is A1 B1 C1 A2 B2 C2 A3 B3
      if(remain>1){
        remain.sets.indices=indices[seq(divisor+1,(divisor+1)*remain,by=divisor+1)]
      }else if(remain==1){
        remain.sets.indices=indices[divisor+1]
      }else{
        remain.sets.indices=NULL
      }
      reorder.indices=c(unlist(full.sets.indices),remain.sets.indices)
      }else{
        reorder.indices=unlist(lapply(seq(1,divisor),
                                 function(x)indices[seq(x,num.cont,by=divisor)]))
        }
      cont.mat=cont.mat[,reorder.indices]
      colnames(cont.mat)=paste0("x",seq(1,num.cont))
    }
  }else{
    cont.mat=NULL
  }
  # generate treatments
  treat.col=data.frame("treatment"=randomizr::complete_ra(N=nobs,num_arms = num.treat,conditions=c(paste0("treatment",seq(1,num.treat-1)),"control")))
  treat.mat=stats::model.matrix(~.,data=treat.col)
  treat.mat=data.frame(treat.mat[,seq(2,ncol(treat.mat)),drop=F])


  #generate categorical columns
  cat.mat=NULL
  if(num.cat>0){ #if there are categorical columns
    if(is.null(cat.counts)==FALSE){  #### if cat.counts is not null, we use it to define counts
      stopifnot(is.list(cat.counts)==T)
      if(!((sapply(cat.counts,sum)==nobs)|(sapply(cat.counts,sum)==0))){
        stop("cat.counts entries must be NULL or sum to the number of observations")
      }
      if(length(cat.counts)==1){
        stopifnot(sum(cat.counts[[1]])==nobs)
        #if cat.probs is a list, and together cat.counts and cat.probs account for all num.cat variables
        if((is.list(cat.probs)==TRUE)&(length(cat.probs)+1==num.cat)){
          cat.count=c(cat.count,lapply(seq(1,length(cat.probs)),function(x)stats::rmultinom(1,nobs,cat.probs[[x]])))
        }else if((is.list(cat.probs)==FALSE)&(sum(cat.probs)==1)&(num.cat==2)){
          #if cat.probs is a vector of probabilities summing to 1 and there are only 2 categorical variables....
          cat.counts=c(cat.counts,stats::rmultinom(1,nobs,cat.probs[[x]]))
        }else if(is.null(cat.probs)==TRUE){
          #if cat.counts is a list of length 1, and no cat.probs is provided.
          #Default to all categorical variables have the same counts
          cat.counts=rep(list(cat.counts),num.cat)
        }
      }else{#list length>1
        null.count.indic=sapply(cat.counts,is.null)
        null.prop.indic=sapply(cat.probs,is.null)
        if(sum(!null.count.indic)!=num.cat){ #not enough categorical variable counts for num.cat
          if((is.null(cat.probs)==TRUE)|((length(null.prop.indic)>0)&&(sum(!null.count.indic)+sum(!null.prop.indic)!=num.cat))){
            stop("Number of cat.counts and cat.probs do not match the number of categorical variables, num.cat.")
          }else{ #there are enough cat.counts and cat.probs to make num.cat columns
            cat.counts=c(cat.counts,rep(list(NULL),num.cat-length(cat.counts)))
            null.count.indic=sapply(cat.counts,is.null)
            null.count.index=which(null.count.indic)
            notnull.prop.index=which(!null.prop.indic)
            cat.count.probs=lapply(notnull.prop.index,function(x)stats::rmultinom(1,nobs,cat.probs[[x]]))
            for(i in seq(1,length(null.count.index))){
              cat.count[[null.count.index[i]]]=cat.count.probs[[i]]
            }
          } #end there are enough cat.counts and cat.probs
        } #end length cat.count != num.cat
      }
    }else{ #if cat.counts is NULL
      stopifnot(is.null(cat.probs)==FALSE)
      if((is.list(cat.probs)==TRUE&length(cat.probs)==1)){ #if only one cat.probs supplied, use it for all variables
        stopifnot(sum(cat.probs[[1]])==1)
        cat.counts=stats::rmultinom(num.cat,nobs,cat.probs[[1]])
      }else if(is.list(cat.probs)==FALSE){
        stopifnot(sum(cat.probs)==1)
        cat.counts=stats::rmultinom(num.cat,nobs,cat.probs)
      }else{
        stopifnot(length(cat.probs)==num.cat)
        cat.counts=sapply(seq(1,num.cat),function(x)stats::rmultinom(1,nobs,cat.probs[[x]]))
      }
    }
      #for each categorical column
      for(catidx in seq(1,num.cat)){
        cprefix=LETTERS[catidx%%26] #get a letter (letters repeat after 26 columns)
        #get group numbers
        if(is.list(cat.counts)==TRUE){
          groups=sample(rep(seq(1,length(cat.counts[[catidx]])),times=cat.counts[[catidx]]))
        }else{
          groups=sample(rep(seq(1,length(cat.counts[,catidx])),times=cat.counts[,catidx]))
        }
        cat.mat=cbind(cat.mat,paste0(cprefix,groups))
      }
      colnames(cat.mat)=paste0("x",seq(num.cont+1,num.cat+num.cont))
  }

  sim.data=data.frame(dplyr::bind_cols(treat.mat,cont.mat,cat.mat))
  colnames(sim.data)=c(
                       paste0("t",seq(1,ncol(treat.mat))),
                       paste0("x",seq(1,num.cat+num.cont)))
  modMat=stats::model.matrix(formula(paste0("~",paste0(colnames(sim.data),collapse="+"))),
                             data=sim.data)
    #regression matrix form
  if(is.null(reg.mods)==TRUE){
    fitted.y=modMat%*%mod.coefs #fitted response
    res.err=stats::rnorm(nobs,mean=0,sd=residual.sd) #random residual
    #combine response, treatment indicators, continuous and categorical columns
    sim.data$y=unlist(as.numeric(fitted.y+res.err))
    sim.data$treatment=treat.col
  }else{
    if(length(reg.mods)==1){
    pred.vars=base::all.vars(stats::as.formula(reg.mods))[3]
    pred.var.indx=which(colnames(sim.data)%in%pred.vars)
    fitted.y=modMat[,pred.var.indx]%*%mod.coefs[pred.var.indx]
    res.err=stats::rnorm(nobs,mean=0,sd=residual.sd) #random residual
    sim.data$y=unlist(as.numeric(fitted.y+res.err))
    }else{
      for(i in seq(1,length(reg.mods))){
        pred.vars=base::all.vars(stats::as.formula(reg.mods[i]))[3]
        pred.var.indx=which(colnames(sim.data)%in%pred.vars)
        fitted.y=modMat[,pred.var.indx,drop=F]%*%mod.coefs[pred.var.indx]
        res.err=stats::rnorm(nobs,mean=0,sd=residual.sd) #random residual
        sim.data[,paste0("y",i)]=unlist(as.numeric(fitted.y+res.err))
      }
    }
  }
    sim.data$control=ifelse(rowSums(treat.mat)==0,1,0)
    sim.data$treatment=as.factor(unlist(treat.col))
    levels(sim.data$treatment)=c("control",paste0("t",seq(1,num.treat-1)))#nlevels(sim.data$treatment)-1))))
    return(sim.data)
}



