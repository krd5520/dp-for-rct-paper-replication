main_simulation_modelcompare<-function(combos.NCov,nobs,nsims,
                                       response="y",treat.vars="t1",
                                       add.formulas=NULL,add.names=NULL,add.ncat=NULL,add.ncont=NULL,
                                       intercept,treat.effect,cat.gseq,cont.gseq,resid.sd,
                                       cat.probs,cont.funcs,cont.params,
                                       priv.eps=1,priv.delta=0,dpmodel.higher.budget=TRUE,
                                       sim.data.rseed=1,
                                       conf.suffix="sim_models",file.suffix="",
                                       output.folder="~",
                                       conf.plot=list("pt.sz1"=1,"ln.sz1"=1.5,"bs.sz1"=12,
                                                      "pt.sz2"=1,"ln.sz2"=1.5,"bs.sz2"=15,
                                                      "ncol1"=3,"ncol2"=2,
                                                      "width1"=4,"height1"=2.5,
                                                      "width2"=6,"height1"=3.5),
                                       bins.param,use.continuous.noise=TRUE,
                                       mv.prop.budget,win.y.resvar,win.y.trcoef,
                                       bound.means,bound.sds,range.alpha,
                                       n.iters.fulldp,
                                       ci.confidence,stderr.func=NULL,
                                       san.rseed=1,
                                       #se.normal=F,
                                       num.cores=1,continuous.vars=NULL,
                                       save.csv.data.dir=NULL#,standardize.vars=NULL,continuous.limits=NULL,diagnostic.file=NULL

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

  synthdata.budget=rep(list("Budget1"=c(priv.eps,priv.delta)),length(reg.formulas.sim2)) #overall budgets epsilon, delta


  gen.sim.start=proc.time()
  times.list=list("SetUp"=(gen.sim.start-main.start)[[3]])
  sim.dfs=generate_simulated_dataset(num.cat=ncat,num.cont=ncont,num.treat=ntreat,nobs=nobs,
                                     mod.coefs=mod.coefs,residual.sd=resid.sd,
                                     cat.probs=cat.probs,cont.funcs=cont.funcs,
                                     cont.params=cont.params,
                                     rseed=sim.data.rseed)


  sim.dfs=sim.dfs[,!(colnames(sim.dfs)%in%c("control","treatment"))]
  if(is.null(continuous.vars)==TRUE){
    cont.vars.all=c(paste0("x",seq(1,ncont)),"y")
  }else{
    cont.vars.all=continuous.vars
  }
  if(is.null(save.csv.data.dir)==F){
    write.csv(sim.dfs,file=paste0(save.csv.data.dir,"/generated_",conf.suffix,"_",file.suffix,".csv"))
  }
    #print(colnames(sim.dfs))
  #print(cont.vars.all1)

  # if(is.null(standardize.vars)==TRUE){
  #   std.vars.all=c(paste0("x",seq(1,ncont)))
  # }else{
  #   std.vars.all=standardize.vars
  # }
  # if(response%in%std.vars.all){
  #   sim.dfs[,response]=(sim.dfs[,response]-stats::mean(sim.dfs[,response],na.rm=T))/(base::sqrt(stats::var(sim.dfs[,response],na.rm=T)))
  #   std.vars.all=std.vars.all[std.vars.all!=response]
  # }
  #
  # if(is.null(names(continuous.limits))==TRUE){
  # names(continuous.limits)=continuous.vars
  # }
  # continuous.limits=continuous.limits[continuous.vars]


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

  #print("done with conf data")
  conf.analysis.end=proc.time()
  times.list=c(times.list,list("conf.analysis"=(conf.analysis.end-conf.analysis.start)[[3]]))
  save(conf.plots,add.plots,conf.summaries,sim.dfs,
       reg.formulas.sim2,model.params,times.list,
       file=paste0(output.folder,"/simulations/checkpoint_",conf.suffix,"_",file.suffix,".Rda"))


  start.simulations=proc.time()
  sim2.out=parallel::mclapply(seq(1,nsims),function(idx)sim_compare_methods(simdata=sim.dfs,
                                                                            reg.models=reg.formulas.sim2,
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
                                                                            se.func=stderr.func))#,
                                                                            #se.normal=se.normal,continuous.limits=continuous.limits,standardize.col=std.vars.all))

  if(dpmodel.higher.budget==TRUE){
    print("in higher.budget==TRUE")
    synthdata.budget.double=rep(list(c(priv.eps/mv.prop.budget,priv.delta/mv.prop.budget)),length(reg.formulas.sim2)) #overall budgets epsilon, delta
    mod.names.double=paste0("HigherBudget_",names(reg.formulas.sim2))
    sim2.out.onlyfull=parallel::mclapply(seq(1,nsims),function(idx)sim_compare_methods(simdata=sim.dfs,
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
                                                                                       se.func=stderr.func))#,
                                                                                       #se.normal=se.normal,continuous.limits=continuous.limits,standardize.col=std.vars.all))

  }

  end.simulations=proc.time()
  times.list=c(times.list,list("simulations"=(end.simulations-start.simulations)[[3]]))

  save(sim2.out,sim2.out.onlyfull,file=paste0(output.folder,"/simulations/checkpoint_models_mclapplyfin_",file.suffix,".Rda"))
  #get table of Average Treatment Effect, Average Overlap, Average Absolute Diff, Max Std. Err
  ITT2=dplyr::bind_rows(lapply(sim2.out,"[[",1)) #combine all repetitions
  ITT2.onlyfull=dplyr::bind_rows(lapply(sim2.out.onlyfull,"[[",1))
  ITT2=dplyr::bind_rows(ITT2,ITT2.onlyfull)
  if(is.null(add.formulas)==FALSE){
    add.NCov=data.frame("NCat"=add.ncat,"NCont"=add.ncont,"Name"=add.names)
    combos.NCov=dplyr::bind_rows(combos.NCov,add.NCov)
  }
  colnames(combos.NCov)=c("NCat","NCont","Model")
  ITT2=dplyr::left_join(ITT2,combos.NCov,by="Model")
  ITT2$NCovariates=ITT2$NCat+ITT2$NCont


  summary.tab=ITT2%>%group_by(Method,Treatment,Model,NCont,NCat,NCovariates)%>%
    summarise( "Avg.ITT"=mean(ITT),
               "Median.ITT"=median(ITT,na.rm=T),
               "Avg.CI.overlap"=mean(CI.overlap,na.rm=T),
               "Avg.Abs.Err"=mean(Abs.Err,na.rm=T),
               "Median.Abs.Err"=median(Abs.Err,na.rm=T),
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
                   "total.budget"=synthdata.budget,
                   "DPMb.mvhist.budget"=mv.prop.budget,
                   "DPMb.response.budget"=mv.prop.budget,
                   "DPMb.y.resvar"=win.y.resvar,
                   "DPMb.y.trcoef"=win.y.trcoef,
                   "DPMb.y.covcoef"=(1-mv.prop.budget)*(1-win.y.trcoef-win.y.resvar),
                   "DPMb.boundmeans"=bound.means,
                   "DPMb.boundsds"=bound.sds,
                   "DPMb.rangealpha"=range.alpha,
                   "DPMb.niters"=n.iters.fulldp)

  main.end=proc.time()
  times.list=c(times.list,list("summarize"=(main.end-end.simulations)[[3]],
                               "total"=(main.end-main.start)[[3]]))
  #NOTE: RERUN.
  save.file.name=paste0(output.folder,"/simulations/results_models_",file.suffix,".Rda")
  if(dpmodel.higher.budget==TRUE){
    synthdata.budget=c(synthdata.budget,synthdata.budget.double)
    save(sim2.out,sim2.out.onlyfull,ITT2,summary.tab,summary.tab.catcont,
         combos.NCov,synthdata.budget,times.list,
         file=save.file.name)

  }else{
    save(sim2.out,ITT2,summary.tab,summary.tab.catcont,
         combos.NCov,synthdata.budget,times.list,algo.param2,
         file=save.file.name)
  }
  return(list(ITT2,times.list,summary.tab,algo.param2))
}


# main_simulation_modelcompare<-function(combos.NCov,nobs,nsims,
#                                   response="y",treat.vars="t1",
#                                   add.formulas=NULL,add.names=NULL,add.ncat=NULL,add.ncont=NULL,
#                                   intercept,treat.effect,cat.gseq,cont.gseq,resid.sd,
#                                   cat.probs,cont.funcs,cont.params,
#                                   priv.eps=1,priv.delta=0,dpmodel.higher.budget=TRUE,
#                                   sim.data.rseed=1,
#                                   conf.suffix="sim_models",file.suffix="",
#                                   output.folder="~",
#                                   conf.plot=list("pt.sz1"=1,"ln.sz1"=1.5,"bs.sz1"=12,
#                                                  "pt.sz2"=1,"ln.sz2"=1.5,"bs.sz2"=15,
#                                                  "ncol1"=3,"ncol2"=2,
#                                                  "width1"=4,"height1"=2.5,
#                                                  "width2"=6,"height1"=3.5),
#                                   bins.param,use.continuous.noise=TRUE,
#                                   mv.prop.budget,win.y.resvar,win.y.trcoef,
#                                   bound.means,bound.sds,range.alpha,
#                                   n.iters.fulldp,
#                                   ci.confidence,stderr.func=NULL,
#                                   san.rseed=1,
#                                   se.normal=F
#
# ){
#   main.start=proc.time()
#   ntreat=length(treat.effect)
#   reg.formulas.sim2=sim2_models_list(response=response,treat.vars=treat.vars,combos.NCov=combos.NCov,
#                              add.formulas=add.formulas,add.names=add.names)
#   mod.names=names(reg.formulas.sim2)
#
#   ncat=max(combos.NCov$NCat)
#   ncont=max(combos.NCov$NCont)
#   mod.coefs=sim_coefs(ncat=ncat,ncat.groups=sapply(cat.probs,length),ncont=ncont,
#                                 cat.coef.start=cat.gseq[1],cat.coef.ratio=cat.gseq[2],
#                                 cont.coef.start=cont.gseq[1],cont.coef.ratio=cont.gseq[2],
#                                 treat.effect=treat.effect,intercept=intercept)
#   model.params=list("model.coefs"=mod.coefs,"residual.err.sd"=resid.sd)
#   cont.vars.all=c("y",paste0("x",seq(1,ncont)))
#   synthdata.budget=rep(list("Budget1"=c(priv.eps,priv.delta)),length(reg.formulas.sim2)) #overall budgets epsilon, delta
#
#
#   gen.sim.start=proc.time()
#   times.list=list("SetUp"=(gen.sim.start-main.start)[[3]])
#   sim.dfs=generate_simulated_dataset(num.cat=ncat,num.cont=ncont,num.treat=ntreat,nobs=nobs,
#                                      mod.coefs=mod.coefs,residual.sd=resid.sd,
#                                      cat.probs=cat.probs,cont.funcs=cont.funcs,
#                                      cont.params=cont.params,
#                                      rseed=sim.data.rseed)
#
#
#   conf.analysis.start=proc.time()
#   times.list=c(times.list,list("Simulate Data"=(conf.analysis.start-gen.sim.start)[[3]]))
#
#
#   conf.summaries=lapply(seq(1,length(reg.formulas.sim2)),
#                         function(idx)summary(glm(reg.formulas.sim2[idx],data=sim.dfs))$coefficients)
#
#   conf.plots= lapply(seq(1,length(reg.formulas.sim2)-length(add.formulas)),
#                    function(idx)reg_assumptions(reg.formulas.sim2[idx],
#                                                 data=sim.dfs,
#                                                 response.var=response,
#                                                 pt.sz=conf.plot$pt.sz1,
#                                                 ln.sz=conf.plot$ln.sz1,
#                                                 bs.sz=conf.plot$bs.sz1,
#                                                 mod.name=mod.names[idx]))
#   save.plot.file=paste0(output.folder,"/figures/simulations_diagnostic_",conf.suffix,"_",file.suffix,".png")
#   ggplot2::ggsave(filename=save.plot.file,bg="white",
#       width=conf.plot$width1,height=conf.plot$height1,
#       plot=cowplot::plot_grid(plotlist=conf.plots,ncol=conf.plot$ncol1))
#   #dev.off()
#   if(is.null(add.formulas)==FALSE){
#     add.plots=lapply(seq(1+length(reg.formulas.sim2)-length(add.formulas),length(reg.formulas.sim2)),
#                     function(idx)reg_assumptions(reg.formulas.sim2[idx],
#                                                  data=sim.dfs,
#                                                  response.var=response,
#                                                  pt.sz=conf.plot$pt.sz2,
#                                                  ln.sz=conf.plot$ln.sz2,
#                                                  bs.sz=conf.plot$bs.sz2,
#                                                  mod.name=mod.names[idx]))
#     save.plot.name.xtra=paste0(output.folder,"/figures/simulations_diagnostic_",conf.suffix,"_",file.suffix,"_xtra.png")
#     ggplot2::ggsave(filename=save.plot.name.xtra,bg="white"
#         ,width=conf.plot$width2,height=conf.plot$height2,
#         plot=cowplot::plot_grid(plotlist=add.plots,ncol=conf.plot$ncol2))
#     #dev.off()
#   }else{
#     add.plots=NULL
#   }
#
#   print("done with conf data")
#   conf.analysis.end=proc.time()
#   times.list=c(times.list,list("conf.analysis"=(conf.analysis.end-conf.analysis.start)[[3]]))
#   save(conf.plots,add.plots,conf.summaries,sim.dfs,
#        reg.formulas.sim2,model.params,times.list,
#        file=paste0(output.folder,"/simulations/checkpoint_",conf.suffix,"_",file.suffix,".Rda"))
#
#
#   start.simulations=proc.time()
#   sim2.out=parallel::mclapply(seq(1,nsims),function(idx)sim_compare_methods(simdata=sim.dfs,
#                                                                             reg.models=reg.formulas.sim2,
#                                                                             synthdata.budget=synthdata.budget,
#                                                                             bins.param=bins.param,
#                                                                             use.continuous.noise=use.continuous.noise,
#                                                                             treat.vars=treat.vars,
#                                                                             cont.as.cont=cont.vars.all,
#                                                                             mv.prop.budget=mv.prop.budget,
#                                                                             win.y.resvar=win.y.resvar,
#                                                                             win.y.trcoef=win.y.trcoef,
#                                                                             bound.means=bound.means,
#                                                                             bound.sds=bound.sds,
#                                                                             n.iters.fulldp=n.iters.fulldp,
#                                                                             range.alpha=range.alpha,
#                                                                             mod.names=mod.names,
#                                                                             ci.confidence=ci.confidence,
#                                                                             rseed=san.rseed+idx,
#                                                                             sim=paste0("Sim",idx),
#                                                                             include.full.model.based=TRUE,
#                                                                             true.coef.vals=c(treat.effect),
#                                                                             se.func=stderr.func,
#                                                                             se.normal=se.normal))
#
#   if(dpmodel.higher.budget==TRUE){
#     print("in higher.budget==TRUE")
#   synthdata.budget.double=rep(list(c(priv.eps/mv.prop.budget,priv.delta/mv.prop.budget)),length(reg.formulas.sim2)) #overall budgets epsilon, delta
#   mod.names.double=paste0("HigherBudget_",names(reg.formulas.sim2))
#   sim2.out.onlyfull=parallel::mclapply(seq(1,nsims),function(idx)sim_compare_methods(simdata=sim.dfs,
#                                                                                      reg.models=reg.formulas.sim2,
#                                                                                      synthdata.budget=synthdata.budget.double,
#                                                                                      bins.param=bins.param,
#                                                                                      use.continuous.noise=use.continuous.noise,
#                                                                                      treat.vars=treat.vars,
#                                                                                      cont.as.cont=cont.vars.all,
#                                                                                      mv.prop.budget=mv.prop.budget,
#                                                                                      win.y.resvar=win.y.resvar,
#                                                                                      win.y.trcoef=win.y.trcoef,
#                                                                                      bound.means=bound.means,
#                                                                                      bound.sds=bound.sds,
#                                                                                      n.iters.fulldp=n.iters.fulldp,
#                                                                                      range.alpha=range.alpha,
#                                                                                      mod.names=mod.names.double,
#                                                                                      ci.confidence=ci.confidence,
#                                                                                      rseed=san.rseed+idx,
#                                                                                      sim=paste0("Sim",idx),
#                                                                                      include.full.model.based=TRUE,
#                                                                                      only.full.model.based=TRUE,
#                                                                                      true.coef.vals=c(treat.effect),
#                                                                                      se.func=stderr.func,
#                                                                                      se.normal=se.normal))
#
#   }
#
#   end.simulations=proc.time()
#   times.list=c(times.list,list("simulations"=(end.simulations-start.simulations)[[3]]))
#
#   ITT2=dplyr::bind_rows(lapply(sim2.out,"[[",1))
#   ITT2.onlyfull=dplyr::bind_rows(lapply(sim2.out.onlyfull,"[[",1))
#   ITT2=dplyr::bind_rows(ITT2,ITT2.onlyfull)
#   if(is.null(add.formulas)==FALSE){
#     add.NCov=data.frame("NCat"=add.ncat,"NCont"=add.ncont,"Name"=add.names)
#     combos.NCov=dplyr::bind_rows(combos.NCov,add.NCov)
#   }
#   colnames(combos.NCov)=c("NCat","NCont","Model")
#   ITT2=dplyr::left_join(ITT2,combos.NCov,by="Model")
#   ITT2$NCovariates=ITT2$NCat+ITT2$NCont
#   summary.tab=ITT2%>%group_by(Method,Treatment,Model,NCont,NCat,NCovariates)%>%
#     summarise( "Avg.ITT"=mean(ITT),
#                "Avg.CI.overlap"=mean(CI.overlap,na.rm=T),
#                "Avg.Abs.Err"=mean(Abs.Err,na.rm=T),
#                "Avg.Rel.Err"=mean(Rel.Err,na.rm=T),
#                "Max.StdErr"=max(StdErr)
#     )
#   #reorder Method factor levels
#   summary.tab$Method=factor(summary.tab$Method,levels=c("Confidential","MV Histogram","Hybrid","DP Model-based"))
#   summary.tab$ModelLab=factor(summary.tab$Model,levels=combos.NCov$Model)
#   summary.tab.catcont=summary.tab%>%#filter(!grepl("Finite|Bounded",Model))%>%
#     arrange(Method)%>%arrange(NCat)%>%arrange(NCont)
#   rownames(summary.tab.catcont)=NULL
#
#
#   algo.param2=list("bin.param"=bins.param,
#                   "total.budget"=synthdata.budget,
#                   "DPMb.mvhist.budget"=mv.prop.budget,
#                   "DPMb.response.budget"=mv.prop.budget,
#                   "DPMb.y.resvar"=win.y.resvar,
#                   "DPMb.y.trcoef"=win.y.trcoef,
#                   "DPMb.y.covcoef"=(1-mv.prop.budget)*(1-win.y.trcoef-win.y.resvar),
#                   "DPMb.boundmeans"=bound.means,
#                   "DPMb.boundsds"=bound.sds,
#                   "DPMb.rangealpha"=range.alpha,
#                   "DPMb.niters"=n.iters.fulldp)
#
#   main.end=proc.time()
#   times.list=c(times.list,list("summarize"=(main.end-end.simulations)[[3]],
#                                "total"=(main.end-main.start)[[3]]))
#   #NOTE: RERUN.
#   save.file.name=paste0(output.folder,"/simulations/results_models_",file.suffix,".Rda")
#   if(dpmodel.higher.budget==TRUE){
#     synthdata.budget=c(synthdata.budget,synthdata.budget.double)
#     save(sim2.out,sim2.out.onlyfull,ITT2,summary.tab,summary.tab.catcont,
#          combos.NCov,synthdata.budget,times.list,
#          file=save.file.name)
#
#   }else{
#     save(sim2.out,ITT2,summary.tab,summary.tab.catcont,
#          combos.NCov,synthdata.budget,times.list,algo.param2,
#          file=save.file.name)
#   }
#   return(list(ITT2,times.list,summary.tab,algo.param2))
# }
