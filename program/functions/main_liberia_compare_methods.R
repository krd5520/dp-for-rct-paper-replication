main_liberia_compare_methods<-function(in.datapath,file.suffix,outputpath,out.datapath,
                                       bins.param,
                                                 synthdata.budget,use.continuous.noise,
                                                 mv.prop.budget,win.y.resvar,win.y.trcoef,
                                                 bound.means,bound.sds,n.iters,range.alpha=0.05,
                                                 ci.confidence=0.95,each.rseed=NULL,
                                       use.allobs=FALSE,use.obs1=FALSE,
                                       use.genmod.full=FALSE,
                                       otherbudgetsDP=NULL,otherbudgetsDP.names=NULL,
                                       reg.assumption.plots=F,
                                       reg.plot.props=NULL,
                                       override.checkpoints=F,
                                       n.std.dev=5,continuous.limits=NULL,continuous.vars=NULL,
                                       standardize.vars=NULL,diagnostic.file=NULL,
                                       which.sets=c("Reduced","Subset","Full"),
                                       var.df=NULL,cat.vars=NULL,response.vars=NULL){
  start.main=proc.time()

  ## Considering 3 sets of covariates
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
  if(is.null(response.vars)==TRUE){
  response.vars.all = c('fam_asb_lt',
                        paste0(c('drugssellever','stealnb','disputes_all_z','carryweapon',
                                 'arrested','asbhostilstd','domabuse_z'),"_ltav"))
  }else{
    response.vars.all=response.vars
  }
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
  if(is.null(continuous.vars)==T){
    col.class=attr(liberia.sub,"column.classifications")
    if(is.null(cat.vars)==TRUE){
    cat.vars=col.class$cat.vars
    }
    cont.as.cont=col.class$cont.as.cont
    cont.as.cat=col.class$cont.as.cat
    #list2env(attr(liberia.sub,"column.classifications"),envir=enviroment())#globalenv())
    cont.vars.all=c(cont.as.cont,cont.as.cat,response.vars.all)
    cont.vars.all=cont.vars.all[!duplicated(cont.vars.all)]
  }else{
    cont.as.cont=continuous.vars
  }

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

  #if(is.null(cat.vars)==FALSE){
  #  liberia.sub[,cat.vars]=apply(liberia.sub[,cat.vars],2,
  #                                 function(x)as.factor(as.character(x)))
  #}

  if(is.null(var.df)==FALSE){
    index.vars.df=var.df[var.df$type=="Index",]
    cont.as.cat=index.vars.df$variable[1+index.vars.df$max-index.vars.df$min<=mv.bins]
    need.cat.levels=index.vars.df$variable[1+index.vars.df$max-index.vars.df$min>index.vars.df$nVals]
    if(sum(need.cat.levels%in%cont.as.cat)>0){
      cont.as.cat=cont.as.cat[!(cont.as.cat%in%need.cat.levels)]
    }
    cont.as.cont=cont.as.cont[!(cont.as.cont%in%c(cat.vars,cont.as.cat))]
  }

  #save(out.df,comp.times,file=paste0(basepath,"/inst/tables/ITT_compare_methods_data.Rda"))



  #liberia.sub[,!(colnames(liberia.sub)%in%c("treatment","control"))]

  if(is.null(continuous.limits)==FALSE){
    if(is.null(names(continuous.limits))==TRUE){
      names(continuous.limits)=continuous.vars
    }else{
      continuous.limits=continuous.limits[continuous.vars]
    }
  }
  continuous.limits=continuous.limits[cont.as.cont]

  print("Starting Covariate Set Comparisons")

  full.comparison=NULL
  summary.times=NULL
  perturb.hist.times=NULL
  by.response.times=NULL
  if("Full"%in%which.sets){
  full.out.file=paste0(outputpath,"/liberia/liberia_checkpoints/liberia_full_baseline_",file.suffix,".Rda")
  if((file.exists(full.out.file)==TRUE)&(override.checkpoints==FALSE)){
    load(full.out.file)
    full.out=out.list
    #add if file exists
  }else{
    start.full=proc.time()

  full.out=liberia_compare_methods_onecovset(liberia.sub=liberia.sub,
                                   synthdata.budget=synthdata.budget,
                                   mv.bins=mv.bins,use.continuous.noise=use.continuous.noise,
                                   covariate.vars=sub.covariates,
                                   block.vars=block.vars,treat.vars=treat.vars,response.vars.all=response.vars.all,
                                   cont.as.cont=cont.as.cont,
                                   mv.prop.budget=mv.prop.budget,win.y.resvar=win.y.resvar,win.y.trcoef=win.y.trcoef,
                                   bound.means=bound.means,bound.sds=bound.sds,n.iters=n.iters,range.alpha=range.alpha,
                                   mod.names=mod.names,treat.names=treat.names,
                                   mult.test.correct=c("treatments","responses"),bonferroni.npvals=1,
                                   rseed=each.rseed,
                                   sim="Full Baseline",
                                   include.dp.model=use.genmod.full,
                                   otherbudgetsDP=otherbudgetsDP,
                                   ci.confidence=ci.confidence,fname=full.out.file,
                                   reg.assumption.plots=reg.assumption.plots,
                                   reg.plot.path=paste0(outputpath,"/figures/liberia_",file.suffix),
                                   reg.plot.prefix="liberia_diagnostics_full",
                                   reg.plot.props=reg.plot.props,se.normal=T,

                                   nstddev=n.std.dev,cont.limits=continuous.limits,
                                   std.vars=standardize.vars,diag.file=diagnostic.file,checkpointfname=paste0(outputpath,"/liberia/liberia_checkpoints/full_raw_data_",file.suffix,".Rda"))
  CON=file(diagnostic.file,"a")
  writeLines(c(" "," ",paste0(rep("#",20),collapse=""),paste0("Full Baseline: ",length(baseline.full)," covariates. Use GenModel? ",use.genmod.full),
               paste("Time Elapsed:",(start.full-proc.time())[[3]]),paste0(rep("#",20),collapse="")),CON)
  close(CON)

  #save(full.out,file=full.out.file)
  }
  full.comparison=full.out[[1]]
  summary.times=unlist(full.out[[2]])
  perturb.hist.times=full.out[[3]]
  by.response.times=full.out[[4]]
  print("Done with Full Covariates")
  print(".")
  }else{
    print("Skipped Full Covariates")
    print(".")
  }

  if("Subset"%in%which.sets){
  sub.out.file=paste0(outputpath,"/liberia/liberia_checkpoints/liberia_subset_",file.suffix,".Rda")
  if((file.exists(sub.out.file)==TRUE)&(override.checkpoints==FALSE)){
    load(sub.out.file)
    sub.out=out.list
    #add if file exists
  }else{
    start.sub=proc.time()
  sub.out=liberia_compare_methods_onecovset(liberia.sub=liberia.sub,
                                  synthdata.budget=synthdata.budget,
                                  mv.bins=mv.bins,use.continuous.noise=use.continuous.noise,
                                  covariate.vars=sub.covariates,
                                  block.vars=block.vars,treat.vars=treat.vars,response.vars.all=response.vars.all,
                                  cont.as.cont=cont.as.cont,
                                  mv.prop.budget=mv.prop.budget,win.y.resvar=win.y.resvar,win.y.trcoef=win.y.trcoef,
                                  bound.means=bound.means,bound.sds=bound.sds,n.iters=n.iters,range.alpha=range.alpha,
                                  mod.names=mod.names,treat.names=treat.names,
                                  mult.test.correct=c("treatments","responses"),bonferroni.npvals=1,
                                  rseed=each.rseed,
                                  sim="Subset",
                                  include.dp.model=TRUE,
                                  otherbudgetsDP=otherbudgetsDP,
                                  ci.confidence=ci.confidence,fname=sub.out.file,
                                  reg.assumption.plots=reg.assumption.plots,
                                  reg.plot.path=paste0(outputpath,"/figures/liberia_",file.suffix),
                                  reg.plot.prefix="liberia_diagnostics_subset",
                                  reg.plot.props=reg.plot.props,se.normal=T,

                                  nstddev=n.std.dev,cont.limits=continuous.limits,
                                  std.vars=standardize.vars,diag.file=diagnostic.file,checkpointfname=paste0(outputpath,"/liberia/liberia_checkpoints/subset_raw_data_",file.suffix,".Rda"))
  CON=file(diagnostic.file,"a")
  writeLines(c(" "," ",paste0(rep("#",20),collapse=""),paste0("Subset: ",length(sub.covariates)," covariates."),
               paste("Time Elapsed:",(start.sub-proc.time())[[3]]),paste0(rep("#",20),collapse="")),CON)
  close(CON)
  #save(sub.out,file=sub.out.file)
  }
  full.comparison=dplyr::bind_rows(full.comparison,
                                   sub.out[[1]])
  summary.times=dplyr::bind_rows(summary.times,
                                 unlist(sub.out[[2]]))
  perturb.hist.times=dplyr::bind_rows(perturb.hist.times,
                                      sub.out[[3]])
  by.response.times=dplyr::bind_rows(by.response.times,
                                     sub.out[[4]])
  print("Done with Subset Covariates")
  print(".")
  }else{
    print("Skipped Subset Covariates")
    print(".")
  }

  if("Reduced"%in%which.sets){
  red.out.file=paste0(outputpath,"/liberia/liberia_checkpoints/liberia_reduced_subset_",file.suffix,".Rda")
  if((file.exists(red.out.file)==TRUE)&(override.checkpoints==FALSE)){
    load(red.out.file)
    red.out=out.list
    #add if file exists
  }else{
    start.red=proc.time()
  red.out=liberia_compare_methods_onecovset(liberia.sub=liberia.sub,
                                  synthdata.budget=synthdata.budget,
                                  mv.bins=mv.bins,use.continuous.noise=use.continuous.noise,
                                  covariate.vars=red.sub.covariates,
                                  block.vars=block.vars,treat.vars=treat.vars,response.vars.all=response.vars.all,
                                  cont.as.cont=cont.as.cont,
                                  mv.prop.budget=mv.prop.budget,win.y.resvar=win.y.resvar,win.y.trcoef=win.y.trcoef,
                                  bound.means=bound.means,bound.sds=bound.sds,n.iters=n.iters,range.alpha=range.alpha,
                                  mod.names=mod.names,treat.names=treat.names,
                                  mult.test.correct=c("treatments","responses"),bonferroni.npvals=1,
                                  ci.confidence=ci.confidence,
                                  rseed=each.rseed,
                                  otherbudgetsDP=otherbudgetsDP,
                                  sim="Reduced Subset",
                                  include.dp.model=TRUE,
                                  fname=red.out.file,
                                  reg.assumption.plots=reg.assumption.plots,
                                  reg.plot.path=paste0(outputpath,"/figures/liberia_",file.suffix),
                                  reg.plot.prefix="liberia_diagnostics_reduced",
                                  reg.plot.props=reg.plot.props,se.normal=T,
                                  nstddev=n.std.dev,cont.limits=continuous.limits,
                                  std.vars=standardize.vars,diag.file=diagnostic.file,checkpointfname=paste0(outputpath,"/liberia/liberia_checkpoints/reduced_raw_data_",file.suffix,".Rda"))
  CON=file(diagnostic.file,"a")
  writeLines(c(" "," ",paste0(rep("#",20),collapse=""),paste0("Reduced Subset: ",length(red.sub.covariates)," covariates."),
               paste("Time Elapsed:",(start.red-proc.time())[[3]]),paste0(rep("#",20),collapse="")),CON)
  close(CON)

  #save(red.out,file=red.out.file)
  }
  full.comparison=dplyr::bind_rows(full.comparison,red.out[[1]])
  summary.times=dplyr::bind_rows(summary.times,unlist(red.out[[2]]))
  perturb.hist.times=dplyr::bind_rows(perturb.hist.times,red.out[[3]])
  by.response.times=dplyr::bind_rows(by.response.times,red.out[[4]])
  print("Done with Reduced Subset Covariates")
  print(".")
  }else{
    print("Skipped Reduced Subset Covariates")
    print(".")
  }



  tempCovariateSet=c("Full Baseline","Subset","Reduced Subset")
  tempnPredictors=c(length(block.vars)+length(treat.vars)+length(baseline.full),
                              length(block.vars)+length(treat.vars)+length(sub.covariates),
                              length(block.vars)+length(treat.vars)+length(red.sub.covariates))
  which.sets.nm=c("Full","Subset","Reduced")
  which.sets.idx=c(1,2,3)
  if(length(which.sets)==1){
    which.sets.idx=which(which.sets.nm==which.sets)
  }else if(length(which.sets==2)){
    which.sets.idx=which.sets.idx[which.sets.nm%in%which.sets]
  }
  summary.times$CovariateSet=tempCovariateSet[which.sets.idx]
  summary.times$nPredictors=tempnPredictors[which.sets.idx]

  CON=file(diagnostic.file,"a")
  writeLines(c(" "," ",paste0(rep("#",20),collapse=""),
               paste("Time Elapsed:",(start.main-proc.time())[[3]]),paste0(rep("#",20),collapse="")),CON)
  close(CON)

  saveRDS(full.comparison,file=paste0(out.datapath,"/liberia_compare_methods_data_",file.suffix,".Rds"))
  save(full.comparison,summary.times,perturb.hist.times,by.response.times,
       #full.comparison.2budgets,summary.times.2budgets,perturb.hist.times.2budgets,by.response.times.2budgets,
       file=paste0(outputpath,"/liberia/liberia_ITT_compare_methods_data_",file.suffix,".Rda"))
}
