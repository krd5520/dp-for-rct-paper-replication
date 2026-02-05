main_simulation_budgetcompare<-function(nobs,ncat,ncont,nsims,
                         response="y",treat.vars="t1",
                         intercept,treat.effect,cat.gseq,cont.gseq,resid.sd,
                         cat.probs,cont.funcs,cont.params,
                         synthdata.budget,
                         sim.data.rseed=1,
                         conf.suffix="sim_budgets",file.suffix="",
                         output.folder="~",
                         bins.param,use.continuous.noise=TRUE,
                         mv.prop.budget,win.y.resvar,win.y.trcoef,
                         bound.means,bound.sds,range.alpha,
                         n.iters.fulldp,
                         ci.confidence,stderr.func=NULL,
                         san.rseed=1,se.normal=F,num.cores=1,
                         #n.std.dev=15,continuous.limits=NULL,
                         continuous.vars=NULL,
                         #standardize.vars=NULL,
                         #diagnostic.file=NULL,
                         save.csv.data.dir=NULL){
  main.start=proc.time()
  # ntreat=length(treat.vars)+1
  # if("control" %in% treat.vars){
  #   ntreat=length(treat.vars)+1
  # }else{
  #   treat.vars=c(treat.vars,"control")
    ntreat=length(treat.vars)
  #}

  reg.formulas.sim1=rep(paste0(response,"~",paste0(treat.vars,collapse="+"),"+",
                               paste0("x",seq(1,ncat+ncont),collapse="+")),
                        times=length(synthdata.budget))
  names(reg.formulas.sim1)=names(synthdata.budget)
  response.vars=rep(response,times=length(reg.formulas.sim1))

  synthdata.budget.eps=sapply(synthdata.budget,function(x)x[1])
  synthdata.budget.del=sapply(synthdata.budget,function(x)x[2])

  mod.names=names(synthdata.budget)

  print("before mod.coefs")

  mod.coefs=sim_coefs(ncat=ncat,ncat.groups=sapply(cat.probs,length),ncont=ncont,
                      cat.coef.start=cat.gseq[1],cat.coef.ratio=cat.gseq[2],
                      cont.coef.start=cont.gseq[1],cont.coef.ratio=cont.gseq[2],
                      treat.effect=treat.effect,intercept=intercept)
  model.params=list("model.coefs"=mod.coefs,"residual.err.sd"=resid.sd)

  if(is.null(continuous.vars)==TRUE){
    cont.vars.all=c(paste0("x",seq(1,ncont)),"y")
  }else{
    cont.vars.all=continuous.vars
  }

  # if(is.null(standardize.vars)==TRUE){
  #   std.vars.all=c(paste0("x",seq(1,ncont)))
  # }else{
  #   std.vars.all=standardize.vars
  #
  # }
  #
  # if(is.null(continuous.limits)==FALSE){
  #   if(is.null(names(continuous.limits))==TRUE){
  #     names(continuous.limits)=continuous.vars
  #   }else{
  #     continuous.limits=continuous.limits[continuous.vars]
  #   }
  # }


  gen.sim.start=proc.time()
  times.list=list("SetUp"=(gen.sim.start-main.start)[[3]])
  print("before generate_simulated_dataset")#bef)
  sim.dfs=generate_simulated_dataset(num.cat=ncat,num.cont=ncont,num.treat=ntreat,nobs=nobs,
                                     mod.coefs=mod.coefs,residual.sd=resid.sd,
                                     cat.probs=cat.probs,cont.funcs=cont.funcs,
                                     cont.params=cont.params,
                                     rseed=sim.data.rseed)
  print(head(sim.dfs))



  sim.dfs=sim.dfs[,!(colnames(sim.dfs)%in%c("treatment","control"))]
  # if(is.null(diagnostic.file)==FALSE){
  #   out.text=sapply(seq(1,length(continuous.limits)),
  #                   function(i)paste0(names(continuous.limits)[i],
  #                                     ": Limits=(",paste0(continuous.limits[[i]],collapse=","),
  #                                     ") have true range=(",round(min(sim.dfs[,names(continuous.limits)[i]]),3),
  #                                     ", ",round(max(sim.dfs[,names(continuous.limits)[i]]),3)))
  #   CON=file(diagnostic.file,"a")
  #   writeLines(out.text,CON)
  #   close(CON)
  # }



  # if(response%in%std.vars.all){
  #   sim.dfs[,response]=(sim.dfs[,response]-mean(sim.dfs[,response],na.rm=T))/(base::sqrt(stats::var(sim.dfs[,response],na.rm=T)))
  #   std.vars.all=std.vars.all[std.vars.all!=response]
  # }



  times.list=c(times.list,list("Simulate Data"=(proc.time()-gen.sim.start)[[3]]))
  save(sim.dfs,synthdata.budget,reg.formulas.sim1,model.params,times.list,
       file=paste0(output.folder,"/simulations/checkpoints/checkpoint_",conf.suffix,"_",file.suffix,".Rda"))
  if(is.null(save.csv.data.dir)==F){
    #print(colnames(sim.dfs))
    write.csv(sim.dfs,file=paste0(save.csv.data.dir,"/generated_",conf.suffix,"_",file.suffix,".csv"))
  }
  print(cat(paste0("in main_simualtion_budget_compare continuous variables are:", paste0(cont.vars.all,collapse=", "))))

  #print(head(sim.dfs))
  print(cat("Simulated Data generated and Confidential Model Fit."))
  start.simulations=proc.time()
  sim1.out=parallel::mclapply(
      seq(1,nsims),
      function(idx)
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
                            se.func=stderr.func#,#se.normal=se.normal,
                            #continuous.limits=continuous.limits,standardize.col=std.vars.all,
                            #diagnostic.file=diagnostic.file
                            ),mc.cores=num.cores)

  end.simulations=proc.time()
  times.list=c(times.list,list("simulations"=(end.simulations-start.simulations)[[3]]))

  #print("sims mclapply finished")

  save(sim1.out,file=paste0(output.folder,"/simulations/checkpoint_budgets_",file.suffix,".Rda"))
  #get table of Average Treatment Effect, Average Overlap, Average Absolute Diff, Max Std. Err
  ITT1=dplyr::bind_rows(lapply(sim1.out,"[[",1)) #combine all repetitions
  priv.budget.df=data.frame("Model"=names(synthdata.budget),
                            "Epsilon"=synthdata.budget.eps,
                            "Delta"=synthdata.budget.del)
  ITT1=dplyr::left_join(ITT1,priv.budget.df,by="Model")
  summary.tab=ITT1%>%
    group_by(Model,Method,Treatment,Epsilon,Delta)%>%
    summarise("Avg.ITT"=mean(ITT,na.rm=T),
              "Median.ITT"=median(ITT,na.rm = T),
              "Avg.CI.overlap"=mean(CI.overlap,na.rm=T),
              "Avg.Abs.Err"=mean(Abs.Err,na.rm=T),
              "Median.Abs.Err"=median(Abs.Err,na.rm=T),
              "Avg.Rel.Err"=mean(Rel.Err,na.rm=T),
              "Max.StdErr"=max(StdErr))


  summary.tab$Delta=ifelse(summary.tab$Delta==0,"0","1/1000")
  summary.tab$Algorithm="DP Model-based"
  summary.tab$Algorithm[grepl("Confidential",summary.tab$Method)]="Confidential"
  summary.tab$Algorithm[grepl("MV",summary.tab$Method)]="MV Histogram"
  summary.tab$Algorithm[grepl("Hybrid",summary.tab$Method)]="Hybrid"

  summary.tab.treff=summary.tab[summary.tab$Treatment!="control",]%>%
    arrange(desc(Algorithm))%>%
    arrange(desc(Delta))%>%
    arrange(desc(Epsilon))
  rownames(summary.tab.treff)=NULL


  algo.param1=list("bin.param"=bins.param,
                    "total.budget"=synthdata.budget,
                    "DPMb.mvhist.budget"=mv.prop.budget,
                    "DPMb.response.budget"=(1-mv.prop.budget),
                    "DPMb.y.resvar"=(1-mv.prop.budget)*win.y.resvar,
                    "DPMb.y.trcoef"=(1-mv.prop.budget)*win.y.trcoef,
                    "DPMb.y.covcoef"=(1-mv.prop.budget)*(1-win.y.trcoef-win.y.resvar),
                    "DPMb.boundmeans"=bound.means,
                    "DPMb.boundsds"=bound.sds,
                    "DPMb.rangealpha"=range.alpha,
                    "DPMb.niters"=n.iters.fulldp)

  main.end=proc.time()
  times.list=c(times.list,list("summarize"=(main.end-end.simulations)[[3]],
                               "total"=(main.end-main.start)[[3]]))
  #NOTE: RERUN.
  save.file.name=paste0(output.folder,"/simulations/results_budgets_",file.suffix,".Rda")

  #NOTE: RERUN
  save(sim1.out,summary.tab,summary.tab.treff,times.list,ITT1,algo.param1,sim.dfs,
       file=save.file.name)
  return(list(sim1.out,summary.tab.treff,times.list,algo.param1,sim.dfs))
}
