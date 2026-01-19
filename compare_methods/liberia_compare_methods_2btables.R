####### Testing the multivariate histogram method on the Liberia Data
basepath = rprojroot::find_rstudio_root_file()
main.start=proc.time()
###

load(paste0(basepath,"/output/ITT_compare_methods_data_v0128.Rda"))

mod.names=c("Anti-social Behaviors, z-score",
            "Usually Sells Drugs",
            "# of thefts/robberies in past 2 weeks",
            "Disputes and fights in past 2 weeks, z-score",
            "Carries a weapon on body",
            "Arrested in past 2 weeks",
            "Aggressive behaviors, z-score",
            "Verbal/physical abuse of partner, z-score")
response.vars.all = c('fam_asb_lt',
                      paste0(c('drugssellever','stealnb','disputes_all_z','carryweapon',
                               'arrested','asbhostilstd','domabuse_z'),"_ltav"))

full.comparison$ModelName=factor(full.comparison$ModelName,levels=mod.names)
#levels(full.comparison$ModelName)=response.vars.all
head(full.comparison)
adjp.in.func=(colnames(full.comparison)[grepl("Pvalue",colnames(full.comparison))])[1]=="AdjPvalue"
adj.out=TRUE
if(adjp.in.func==TRUE){
  bonf.adj=NULL
  pvalvar.name="AdjPvalue"
  pvallab.name="Adj. p-value"
}else{
  bonf.adj=1
  pvalvar.name="Pvalue"
  pvallab.name="p-value"
  if(adj.out==TRUE){
    adj.fam=9
    adj.components=3*7
    pvallab.name="Adj. p-value"
  }
}


transform.mv.hist=dplyr::mutate_if(perturb.hist.times,grepl("Time",colnames(perturb.hist.times)),function(x)x*100)
#transform.mv.hist$Algorithm=ifelse(transform.mv.hist$SyntheticMethod=="Hybrid","Hybrid",
#                                   ifelse(transform.mv.hist$SyntheticMethod=="Full DP Model-Based","DP Model-based","MV Histogram"))
transform.mv.hist$Sim[transform.mv.hist$Sim=="Full Baseline"]="Full"
#comp time tables
mvhist_comptime=knitr::kable(transform.mv.hist[c(1,2,seq(4,9)),c(1,9,2,5,6,7)],format='latex',booktabs=T,
                             col.names=c("Method","Covariate Set","Number of Variable in Histogram",
                                        "Make Histogram","Sanitize Proportions","Sample Sanitized Histogram"))
mvhist_comptime

by.response.times$Sim=factor(by.response.times$Sim,levels=c("Reduced Subset","Subset","Full Baseline"))
by.response.times=dplyr::arrange(by.response.times,Sim)
byy.df=by.response.times[,c(3:11)]
byy.df$Step=ifelse(by.response.times$Timed.Step=="fit.conf.model","Fit model with $\\data$",
                        ifelse(by.response.times$Timed.Step=="total.time","Total",
                               ifelse(by.response.times$Timed.Step=="san.summary","Sanitize model parameters",
                                      ifelse(by.response.times$Timed.Step=="iter.proxy.fit","Generate and fit proxy","Generate $\\sanby$"))))
byy.df=byy.df[,c("Step",colnames(by.response.times)[4:11])]
transform.byy.df=byy.df#dplyr::mutate_at(byy.df,colnames(by.response.times)[4:11],function(x)x*100)

byresponse_comptime=knitr::kable(transform.byy.df,format='latex',booktabs=T,escape=F)

starti=1
for(cov.name in c("Reduced Subset (12 predictors)","Subset (15 predictors)","Full (60 predictors)")){
  byresponse_comptime=kableExtra::pack_rows(byresponse_comptime,cov.name,
                                            starti,starti+ifelse(cov.name=="Full (60 predictors)",2,7),
                                            hline_before = T,hline_after = T)
  starti=starti+ifelse(cov.name=="Full (60 predictors)",3,8)
}

starti=1
for(alg.name in c(rep(c("Hybrid (Algorithm \\ref{algo:hybrid})","DP Model-based (Algorithm \\ref{algo:fullmodelbased})"),2),"Hybrid (Algorithm \\ref{algo:hybrid})")){
  nrowspack=ifelse(alg.name=="Hybrid (Algorithm \\ref{algo:hybrid})",3,5)
  byresponse_comptime=kableExtra::pack_rows(byresponse_comptime,alg.name,
                                            starti,starti+nrowspack-1,
                                            escape=F,bold=F,italic=T,indent=F)
  starti=starti+nrowspack
}
byresponse_comptime

sub.summary.times=summary.times[,c(8,9,1,2,3,7,5)]
knitr::kable(sub.summary.times, format="latex",booktabs=T,digits=2,
             col.names = c("Covariate Set","Number of Predictors","Set Up Time","MV Histogram","Hybrid","DP Model-based","Make Table"))

save(mvhist_comptime,byresponse_comptime,
     file=paste0(basepath,"/output/tables/computation_times_covsets_methods.Rda"))

#ITT result tables
full.comparison$ModelName=as.character(full.comparison$ModelName)
#adjust pvalues a
if(adj.out==TRUE&adjp.in.func==FALSE){
  pvalcol.indic=grepl("Pvalue",colnames(full.comparison))
  full.comparison[full.comparison$ModelName=="fam_asb_lt",pvalcol.indic]=
    full.comparison[full.comparison$ModelName=="fam_asb_lt",pvalcol.indic]*adj.fam
  full.comparison[full.comparison$ModelName!="fam_asb_lt",pvalcol.indic]=
    full.comparison[full.comparison$ModelName!="fam_asb_lt",pvalcol.indic]*adj.components
  for(colidx in which(pvalcol.indic)){
    full.comparison[,colidx]=pmin(rep(1,nrow(full.comparison)),t(unlist(full.comparison[,colidx])))
  }
}


sub.compare=full.comparison[full.comparison$data %in% c("MV Histogram", "Hybrid","DP Model-based"),]
cnames=grepl("Model|nObs|ITT|StdErr|Pvalue|overlap|Sim|data",colnames(full.comparison))
sub.compare.all=full.comparison[,cnames]
sub.compare.all=sub.compare.all[,c("ModelName","nObs","ITT_Control",
                                   paste(c("ITT","StdErr",pvalvar.name,"ci.overlap"),"tpassonly",sep="_"),
                                   paste(c("ITT","StdErr",pvalvar.name,"ci.overlap"),"cashassonly",sep="_"),
                                   paste(c("ITT","StdErr",pvalvar.name,"ci.overlap"),"tpcashass",sep="_"),
                                   "Sim","data")]
#sub.compare.all$ModelName=paste0(sub.compare.all$ModelName,ifelse(sub.compare$nObs<946,"\\dagger",""))
sub.compare.san=sub.compare.all[sub.compare.all$data!="Confidential",]
sub.compare.conf=sub.compare.all[sub.compare.all$data=="Confidential",]
tab.df.conf=sub.compare.conf[,c(1,3,4,5,6,8,9,10,12,13,14)]

#original data table
tab.conf=knitr::kable(tab.df.conf,format="latex",digits=2,booktabs=T,
                      col.names = c("Response","Control Mean",rep(c("ITT","Std. Error",pvallab.name),3)))
tab.conf=kableExtra::add_header_above(tab.conf,
                                      header=c(" "=2,
                                               "Therapy Only"=3,"Cash Only"=3,"Both"=3),align="c")
starti=1
for(cov.name in c("Full (60 predictors)","Subset (15 predictors)","Reduced Subset (12 predictors)")){
  tab.conf=kableExtra::pack_rows(tab.conf,cov.name,starti,starti+7)
  starti=starti+8
}
tab.conf


tab.list=NULL
for(cov.name in c("Full Baseline","Subset","Reduced Subset")){
  tab.df.full=sub.compare.san[(sub.compare.san$Sim==cov.name),#&(sub.compare.san$data!="DP Model-Based"),
                         c("ModelName","ITT_Control",
                           paste(c("ITT","StdErr",pvalvar.name,"ci.overlap"),"tpassonly",sep="_"),
                           paste(c("ITT","StdErr",pvalvar.name,"ci.overlap"),"cashassonly",sep="_"),
                           paste(c("ITT","StdErr",pvalvar.name,"ci.overlap"),"tpcashass",sep="_"),"data","Sim")]
  print(tab.df.full$data)
  print(tab.df.full$Sim)
  tab.df=tab.df.full[,seq(1,14)]
  rownames(tab.df)=NULL
  tab.temp=knitr::kable(tab.df[,seq(1,14)],format="latex",digits=2,booktabs=T,
                        col.names = c("Response","Control Mean",rep(c("ITT","Std. Error",pvallab.name,"CI overlap"),3)))
  tab.temp=kableExtra::add_header_above(tab.temp,
                                        header=c(" "=3,
                                                 "Therapy Only"=4,"Cash Only"=4,"Both"=4),align="c")
  method.rows=c("MV Histogram","Hybrid")
  if(cov.name!="Full Baseline"){
    method.rows=c(method.rows,"DP Model-based")
  }
  starti=1
  for(data.name in method.rows){
    print(data.name)
    tab.temp=kableExtra::pack_rows(tab.temp,data.name,starti,starti+7)
    starti=starti+8
  }
  templist=list(tab.temp)
  names(templist)=c(cov.name)
  tab.list=c(tab.list,templist)
}

tab.list[[1]]
tab.list[[2]]
tab.list[[3]]
bycovset.tablist=tab.list

save(tab.conf,bycovset.tablist,
     file=paste0(basepath,"/output/tables/ITTtables_covsets_methods.Rda"))



