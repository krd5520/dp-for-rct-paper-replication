replicate_table2b_and_diagnostic_plots=function(
    liberia.sub,
           covariate.vars=c('age_b', 'livepartner_b', 'mpartners_b', 'hhunder15_b', 'famseeoften_b',
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
                            'timedecl_b', 'riskdecl_b', 'cognitive_score_b', 'ef_score_b'),
           block.vars=c("tp_strata_alt","cg_strata"),
           treat.vars=c("cashassonly","tpassonly","tpcashass"),
           response.vars=c('fam_asb_lt','drugssellever_ltav','stealnb_ltav','disputes_all_z_ltav','carryweapon_ltav',
                           'arrested_ltav','asbhostilstd_ltav','domabuse_z_ltav'),
           output_dir=paste0(rprojroot::find_rstudio_root_file(),"/output"),
           figures_output_path=paste0(rprojroot::find_rstudio_root_file(),"/output/figures/liberia_replicate2b.png"),
           table_output_path=paste0(rprojroot::find_rstudio_root_file(),"/output/tables/liberia_replicate2b.tex"),
           table_figures_rda_path=paste0(rprojroot::find_rstudio_root_file(),"/output/liberia_replication/liberia_replicate2b_table_and_figures_variables.Rda"),
           code_source=rprojroot::find_rstudio_root_file(),
           plot.params=list("point.sz"=4,
                            "line.sz"=2.5,
                            "base.sz"=33,
                            "combined_width"=3000,
                            "combined_height"=2500,
                            "ncol_plotgrid"=2),
           quietly=FALSE){
    main.start=proc.time()

  source(paste0(code_source,"/program/functions/reg_assumptions_function.R"))
  covariate.formula=paste(c(covariate.vars,block.vars),collapse="+")

  #AER_1.do sets block variables to factors
  liberia.sub[,block.vars]=
    apply(liberia.sub[,block.vars],2,
          function(x)as.factor(as.character(x)))

  reg.formula=paste(response.vars,
                    paste(c(treat.vars,covariate.vars,block.vars),collapse="+"),sep="~")

  #in AER_1.do file the ",r" indicates that robust (sandwich estimator) of variance

  ## Uncomment the large chunk below to brute force investigation into which robust estimator is used.
  # sterr.func=function(mod){
  #   sqrt(diag(sandwich::vcovHC(mod)))
  # }
  #
  # attempt_sandwich=try(DPrct::ITTtable(data=liberia.sub,
  #                     families=rep("gaussian",length(reg.formula)),
  #                     reg.models=reg.formula,
  #                     treat.vars=treat.vars,
  #                     control.var="control",
  #                     stderr.func=sterr.func,
  #                     add.pval.stars = FALSE,
  #                     mult.test.correct=NULL))
  # if(is.na(attempt_sandwich)==T){
  #   cat("BJS (2025) uses a robust estimator of variance. The default sandwhich estimator in R does not work though. We will investigate other robust estimators." )
  # }
  # #The robust sandwich method does not work. However there are other types of it.
  #
  # sterr.type="HC0"
  # sterr.func=function(mod,std.type=sterr.type){
  #   sqrt(diag(sandwich::vcovHC(mod,type=std.type)))
  # }
  #
  # type.list=c("HC0","HC1","HC2","HC3","HC4","HC4m","HC5","const")
  # for(sterr.type in type.list){
  #   tab=DPrct::ITTtable(data=liberia.sub,
  #                       families=rep("gaussian",length(reg.formula)),
  #                       reg.models=reg.formula,
  #                       treat.vars=treat.vars,
  #                       control.var = "control",
  #                       stderr.func = sterr.func,
  #                       add.pval.stars = FALSE,
  #                       mult.test.correct=NULL)
  #   if((sum(is.nan(tab$StdErr_cashassonly))==0)&&
  #      (sum(abs(round(tab$StdErr_cashassonly,3)-
  #               c(0.097,0.030,0.388,0.090,0.035,0.025,0.107,0.113)))<10^(-6))){
  #     cat(paste0("ITT table values with sandwhich estimator for robust standard errors type: ",sterr.type))
  #     tempdf=tab[,!grepl("Pvalue",colnames(tab))]
  #     print(tempdf[,colnames(tempdf)!="nObs"])
  #   }
  #
  # }

  #### Anaylsis is done with additive block indicators, _ltav response variables,
  ## and robust HC type="HC1" errors.
  ## (Upon further investigation this corresponds to the default robust errors in Stata)
  sterr.type="HC1"
  sterr.func=function(mod,std.type=sterr.type){
    sqrt(diag(sandwich::vcovHC(mod,type=std.type)))
  }
  liberia.rep.tab=DPrct::ITTtable(data=liberia.sub,
                                  families=rep("gaussian",length(reg.formula)),
                                  reg.models=reg.formula,
                                  treat.vars=treat.vars,
                                  control.var = "control",
                                  stderr.func = sterr.func,
                                  add.pval.stars = TRUE,
                                  mult.test.correct=NULL)


  #ITT result tables

  #adjust pvalues like the analysis of Table 2b by BJS
  pvalcol.indic=grepl("Pvalue",colnames(liberia.rep.tab))
  liberia.rep.tab[liberia.rep.tab$Response=="fam_asb_lt",pvalcol.indic]=
    liberia.rep.tab[liberia.rep.tab$Response=="fam_asb_lt",pvalcol.indic]*9
  liberia.rep.tab[liberia.rep.tab$Response!="fam_asb_lt",pvalcol.indic]=
    liberia.rep.tab[liberia.rep.tab$Response!="fam_asb_lt",pvalcol.indic]*21
  if(quietly==F){
    cat("\nWe use Bonferroni adjustment for multiple testing correction in our ITT tables, but BJS(2017) uses a different correction method. Our p-values won't match exactly.")
  }
  for(colidx in which(pvalcol.indic)){
    liberia.rep.tab[,colidx]=pmin(rep(1,sum(pvalcol.indic)),t(unlist(liberia.rep.tab[,colidx])))
  }

  mod.names=c("Anti-social Behaviors, z-score",
              "Usually Sells Drugs",
              "# of thefts/robberies in past 2 weeks",
              "Disputes and fights in past 2 weeks, z-score",
              "Carries a weapon on body",
              "Arrested in past 2 weeks",
              "Aggressive behaviors, z-score",
              "Verbal/physical abuse of partner, z-score")
  liberia.rep.tab$Outcome=mod.names
  liberia.rep.tab=liberia.rep.tab[,c("Outcome","ControlMean",
                                     paste0(c("ITT","StdErr","Pvalue"),"_tpassonly"),
                                     paste0(c("ITT","StdErr","Pvalue"),"_cashassonly"),
                                     paste0(c("ITT","StdErr","Pvalue"),"_tpcashass"))]

  rep2btab=knitr::kable(liberia.rep.tab,digits=3,
                        col.names=c("Outcome","Control mean",rep(c("ITT","Std. Err","Adj.p-value"),3)),
                        format="latex",booktabs=TRUE)
  rep2btab=kableExtra::add_header_above(rep2btab,header=c(" "=2,"Therapy Only"=3,"Cash Only"=3,"Both"=3))

  dir.create(dirname(table_output_path),showWarnings = FALSE,recursive=TRUE)
  writeLines(rep2btab,table_output_path)
  if(quietly==F){
    cat(paste0("\nReplication of Table2b in Blattman et al. (2017) is saved as a LaTex file: ",table_output_path))
  }

  ### Regression Diagnostic Plots
  list.diagnostic.liberia.plots= lapply(seq(1,length(response.vars)),
                                        function(idx)reg_assumptions(data=liberia.sub,
                                                                     reg.formula[idx],
                                                                     response.var=response.vars[idx],
                                                                     pt.sz=plot.params$point.sz,
                                                                     ln.sz=plot.params$line.sz,
                                                                     bs.sz=plot.params$base.sz))
  dir.create(dirname(figures_output_path),showWarnings = FALSE,recursive=TRUE)
  png(figures_output_path,
      width=plot.params$combined_width,height=plot.params$combined_height)
  print(cowplot::plot_grid(plotlist=list.diagnostic.liberia.plots,
                           ncol=plot.params$ncol_plotgrid))
  dev.off()

  dir.create(dirname(table_figures_rda_path),showWarnings = FALSE,recursive=TRUE)
  save(liberia.rep.tab,list.diagnostic.liberia.plots,
       file=table_figures_rda_path)

  time.liberia.rep2b.orig=(proc.time()-main.start)[[3]]
  if(quietly==F){
    cat(paste0("\nSaving grid of diagnostic plots for each response variable as png :",figures_output_path,
        "\nSaving dataframe of ITT table as 'liberia.rep.tab' and list of plots as 'list.diagnostic.liberia.plots' as an Rda file.",
        "\nLocation: ",table_figures_rda_path))
  }
  return(cowplot::plot_grid(plotlist=list.diagnostic.liberia.plots,
                            ncol=plot.params$ncol_plotgrid))
  }

#liberia.sub=extract_styl_final("C:/Users/Kaitlyn/Downloads/113056-V2.zip")

#outplot=replicate_table2b_and_diagnostic_plots(liberia.sub)

