#### function to create summary table of computation times and table of by response computation times
# by.response.times and summary.times are the dataframes saved from main_liberia_compare_methods
# response.vars is a list of response
liberia_comptimes=function(by.response.times,summary.times,response.vars,
                           by.response.fname,summary.fname,outputpath,
                           secs.100ths=FALSE,change.DPMb.names=NULL){

  ###### by response comp time summary table #########
  lvl.sims=unique(by.response.times$Sim) #Sim Levels (covariate set)
  reduce.name.indic=grepl("Reduce",lvl.sims) #reduced covariate subset indicator
  #order covariate set levels as Reduced, Subset, Full, followed by any other potential sim levels
  lvl.sims=c(lvl.sims[reduce.name.indic],lvl.sims[(!reduce.name.indic)&(grepl("Subset",lvl.sims))],
             lvl.sims[grepl("Full",lvl.sims)],lvl.sims[!grepl("Reduce|Subset|Full",lvl.sims)])
  by.response.times$Sim=factor(as.character(by.response.times$Sim),levels=lvl.sims)

  #change DP Model-based names based on input
  if(is.null(change.DPMb.names)==FALSE){
    for(old.nm in names(change.DPMb.names)){
      by.response.times$SyntheticMethod[by.response.times$SyntheticMethod==old.nm]=change.DPMb.names[[old.nm]]
    }
    order.methods=c("MV Histogram","Hybrid",unname(sapply(change.DPMb.names,"[[",1)))
  }else{
    method.names=unique(by.response.times$SyntheticMethod)
    order.methods=c("MV Histogram","Hybrid",method.names[grepl("DP",method.names,ignore.case = T)])
  }
  #order the methods as MV Histogram, Hybrid, and then DP Model Based (in input order of change.DPMb.names)
  by.response.times$SyntheticMethod=factor(as.character(by.response.times$SyntheticMethod),
                                           levels=order.methods)
  by.response.times=dplyr::arrange(by.response.times,SyntheticMethod)
  by.response.times=dplyr::arrange(by.response.times,Sim)

  #select columns
  byy.df=by.response.times[,c("SyntheticMethod","Sim",colnames(by.response.times)[colnames(by.response.times)%in%response.vars])]

  #change step names
  byy.df$Step=ifelse(by.response.times$Timed.Step=="fit.conf.model","Fit model with $\\data$",
                     ifelse(by.response.times$Timed.Step=="total.time","Total",
                            ifelse(by.response.times$Timed.Step=="san.summary","Sanitize model parameters",
                                   ifelse(by.response.times$Timed.Step=="iter.proxy.fit","Generate and fit proxy",
                                          "Generate $\\sanby$"))))
  byy.df=byy.df[byy.df$Step!="Fit model with $\\data$",]
  #make step the first column
  byy.tab.cols=colnames(by.response.times)[colnames(by.response.times)%in%response.vars]
  byy.df=byy.df[,c("Step","SyntheticMethod",byy.tab.cols,"Sim")]
  if(secs.100ths==TRUE){
    transform.byy.df=dplyr::mutate_at(byy.df,byy.tab.cols,function(x)x*100)
  }else{
    transform.byy.df=byy.df
  }
  tab.df.temp=transform.byy.df[,c("Step",byy.tab.cols)]
  newformat.colnames=paste0("\\rotatebox{83}{",gsub("_","\\_",response.vars,fixed=T),"}")
  colnames(tab.df.temp)=c("Step",newformat.colnames)
  #computation time table by response

  row.names(tab.df.temp)=NULL
  byresponse_comptime=knitr::kable(tab.df.temp,booktabs=T,
                                   format='latex',escape=F)

  #covset.df=data.frame(table(as.character(transform.byy.df$Sim)))
  starti=1
  for(cov.name in lvl.sims){ #group rows by covariate set
    sim.freq=sum(transform.byy.df$Sim==cov.name)#covset.df$Freq[covset.df$Var1==cov.name]
    byresponse_comptime=kableExtra::pack_rows(byresponse_comptime,cov.name,
                                              starti,sim.freq-1,
                                              hline_before = T,hline_after = T)
    starti=starti+sim.freq
  }

  starti=1
  #reduced
  for(alg.name in order.methods[-1]){ #group rows by algorithm
    nrowspack=sum(transform.byy.df$SyntheticMethod==alg.name)/ifelse(grepl("DP",alg.name,ignore.case = T)==TRUE,2,3)
    byresponse_comptime=kableExtra::pack_rows(byresponse_comptime,alg.name,
                                              starti,starti+nrowspack-1,
                                              escape=F,bold=F,italic=T,indent=F)
    starti=starti+nrowspack
  }
  #subset
  for(alg.name in order.methods[-1]){ #group rows by algorithm
    nrowspack=sum(transform.byy.df$SyntheticMethod==alg.name)/ifelse(grepl("DP",alg.name,ignore.case = T)==TRUE,2,3)
    byresponse_comptime=kableExtra::pack_rows(byresponse_comptime,alg.name,
                                              starti,starti+nrowspack-1,
                                              escape=F,bold=F,italic=T,indent=F)
    starti=starti+nrowspack
  }
  #reduced
  for(alg.name in order.methods[2]){ #group rows by algorithm
    nrowspack=sum(transform.byy.df$SyntheticMethod==alg.name)/ifelse(grepl("DP",alg.name,ignore.case=T)==TRUE,2,3)
    byresponse_comptime=kableExtra::pack_rows(byresponse_comptime,alg.name,
                                              starti,starti+nrowspack-1,
                                              escape=F,bold=F,italic=T,indent=F)
    starti=starti+nrowspack
  }

  by.response.file=paste0(outputpath,"/",by.response.fname,".tex")
  writeLines(c("\\begin{table}","\\centering"),by.response.file)
  CON=file(by.response.file,"a")
  writeLines(byresponse_comptime,CON)
  writeLines(c("\\caption{Computation time per response variable for each steps of the model-based algorithms.}",
                     paste0("\\label{tab:liberia-comptimes-byresponse}"),"\\end{table}"),CON)
  close(CON)
  #writeLines(byresponse_comptime,by.response.file)

  ######### comp time summary table #################

  #select relevant columns in order Covariate Set, number of predictors, set up time,
  #     column for each algorithm & budget combo, make table time
  summary.method.colnames=c("Perturb MV Histogram","Hybrid",colnames(summary.times)[grepl("DP",colnames(summary.times),ignore.case=T)])

  col.select=c("CovariateSet","nPredictors","Set Up Time",summary.method.colnames,"Make Table")
  sub.summary.times=summary.times[,col.select]
  rownames(sub.summary.times)=NULL
  summary.tab.out=knitr::kable(sub.summary.times, format="latex",escape=F,digits=2,booktabs=T,
                               col.names = c("Covariate Set","Number of Predictors","Set Up Time",order.methods,"Make Table"))

  summary.file=paste0(outputpath,"/",summary.fname,".tex")
  writeLines(c("\\begin{table}","\\centering"),summary.file)
  CON=file(summary.file,"a")
  writeLines(summary.tab.out,CON)
  writeLines(c("\\caption{Summary of computation times per algorithm and covariate set.}",
                     paste0("\\label{tab:liberia-comptimes-summary}"),"\\end{table}"),CON)
  close(CON)
  #writeLines(summary.tab.out,summary.file)

  save(byresponse_comptime,summary.tab.out,
       file=paste0(outputpath,"/computation_times_covsets_methods.Rda"))
  return(list(byresponse_comptime,summary.tab.out))
}


liberia_table2b_covsets_bymethod=function(full.comparison,n.response=8,
                                          tab.digits=2,
                                          pvallab.name="Adj. p-value",
                                          pvalvar.name="AdjPvalue",
                                          include.control.mean=T,
                                          conf.tex.cols=NULL,#c(0.45,0.07,rep(c(0.07,0.075,0.07),3)),
                                          alg.tex.cols=NULL,#c(0.21,0.072,rep(c(0.061,0.055,0.055,0.055),3)),
                                          outputpath,conf.2b.fname,san.tex.prefix,out.fname,packrows=T){


  lvl.sims=levels(full.comparison$Sim)
  order.methods=levels(full.comparison$data)

  levels(full.comparison$ModelName)=gsub("_","\\_",levels(full.comparison$ModelName),fixed=T)
  #subset columns
  cnames=grepl("Model|nObs|ITT|StdErr|Pvalue|overlap|Sim|data",colnames(full.comparison))
  sub.compare.all=full.comparison[,cnames]
  if(include.control.mean==TRUE){
    sub.compare.all=sub.compare.all[,c("ModelName","ITT_Control",
                                       paste(c("ITT","StdErr",pvalvar.name,"ci.overlap"),"tpassonly",sep="_"),
                                       paste(c("ITT","StdErr",pvalvar.name,"ci.overlap"),"cashassonly",sep="_"),
                                       paste(c("ITT","StdErr",pvalvar.name,"ci.overlap"),"tpcashass",sep="_"),"data")]
    conf.col.names=c("Response","Control Mean",rep(c("ITT","Std. Error",pvallab.name),3))
    alg.col.names= c("Response","Control Mean",rep(c("ITT","Std. Error",pvallab.name,"CI overlap"),3))
    header.n=2
  }else{
    sub.compare.all=sub.compare.all[,c("ModelName",
                                       paste(c("ITT","StdErr",pvalvar.name,"ci.overlap"),"tpassonly",sep="_"),
                                       paste(c("ITT","StdErr",pvalvar.name,"ci.overlap"),"cashassonly",sep="_"),
                                       paste(c("ITT","StdErr",pvalvar.name,"ci.overlap"),"tpcashass",sep="_"),"data")]
    conf.col.names=c("Response",rep(c("ITT","Std. Error",pvallab.name),3))
    alg.col.names= c("Response",rep(c("ITT","Std. Error",pvallab.name,"CI overlap"),3))
    header.n=1
  }

conf.file=paste0(outputpath,"/",conf.2b.fname,".tex")
#confidential table remove nobs, ci overlap columns, Sim and data columns
tab.df.conf=sub.compare.all[sub.compare.all$data=="Confidential",
                            grepl("Model|ITT|StdErr|Pvalue",colnames(sub.compare.all))]
rownames(tab.df.conf)=NULL
#original data table
tab.conf=knitr::kable(tab.df.conf,format="latex",digits=2,escape = F,booktabs=T,
                      col.names = conf.col.names)#c("Response","Control Mean",rep(c("ITT","Std. Error",pvallab.name),3)))
tab.conf=kableExtra::add_header_above(tab.conf,
                                      header=c(" "=header.n,
                                               "Therapy Only"=3,"Cash Only"=3,"Both"=3),align="c")

if(packrows==TRUE){
starti=1
for(cov.name in lvl.sims){ #group covariate sets
  tab.conf=kableExtra::pack_rows(tab.conf,cov.name,starti,starti+n.response-1)
  starti=starti+n.response #8 response variables
}
}


if(is.null(conf.tex.cols)==FALSE){
  for(coli in seq(1,length(conf.tex.cols))){
    tex.colspec=paste0(">{\\\\raggedleft\\\\arraybackslash}p{",conf.tex.cols[coli],"\\\\textwidth}")
    tab.conf=kableExtra::column_spec(tab.conf,column=coli,
                                     latex_valign = "p",
                                     width=paste0("{",conf.tex.cols[coli],"\\\\textwidth}"),
                                     latex_column_spec=tex.colspec)
  }
}

writeLines(c("\\begin{table}","\\centering","\\small"),conf.file)
CON=file(conf.file,"a")
writeLines(tab.conf,CON)
writeLines(c("\\caption{Original data results across the three covariate sets.}\\label{tab:apdx-itt-orig}","\\end{table}"),CON)
close(CON)

print("Done with Confidential data ITT table")


tab.list=NULL
for(alg.name in order.methods){ #for each algorithm
  alg.fname=tolower(unlist(gsub("[^A-z0-9]","",alg.name)))
  alg.fname=gsub("\\","",alg.fname,fixed=T)
  alg.file=paste0(outputpath,"/",san.tex.prefix,"_",alg.fname,".tex")

  tab.df=sub.compare.all[(sub.compare.all$data==alg.name),#&(sub.compare.san$data!="DP Model-Based"),
                         grepl("Model|ITT|StdErr|Pvalue|ci.overlap",colnames(sub.compare.all))]

  rownames(tab.df)=NULL
  tab.temp=knitr::kable(tab.df,format="latex",digits=tab.digits,booktabs=T,escape=F,
                        col.names = alg.col.names)#c("Response","Control Mean",rep(c("ITT","Std. Error",pvallab.name,"CI overlap"),3)))
  tab.temp=kableExtra::add_header_above(tab.temp,
                                        header=c(" "=header.n,
                                                 "Therapy Only"=4,"Cash Only"=4,"Both"=4),align="c")
  if(packrows==TRUE){
    print("inpackrows")
  #group by covariate set
  starti=1
  if(grepl("DP",alg.name,ignore.case=T)==T){
    covset.names=lvl.sims[!grepl("Full",lvl.sims)]
  }else{
    covset.names=lvl.sims
  }
  for(cov.name in covset.names){
    tab.temp=kableExtra::pack_rows(tab.temp,cov.name,starti,starti+n.response-1)
    starti=starti+n.response
  }
  }

  if(is.null(alg.tex.cols)==FALSE){
    for(coli in seq(1,length(alg.tex.cols))){
      tex.colspec=paste0(">{\\\\raggedleft\\\\arraybackslash}p{",alg.tex.cols[coli],"\\\\textwidth}")
      tab.temp=kableExtra::column_spec(tab.temp,column=coli,
                                       latex_valign = "p",
                                       width=paste0("{",alg.tex.cols[coli],"\\\\textwidth}"),
                                       latex_column_spec=tex.colspec)
    }
  }


  writeLines(c("\\begin{table}","\\centering","\\small"),alg.file)
  CON=file(alg.file,"a")
  writeLines(tab.temp,CON)
  writeLines(c(paste("\\caption{Synthetic data from the",alg.name,
                     "method results across the three covariate sets.}",
                     paste0("\\label{tab:apdx-itt-",alg.fname,"}"))),CON)
  close(CON)

starti=1
print(paste("Finished ITT table for", alg.name, "data."))


templist=list(tab.temp)
#names(templist)=c(cov.name)
tab.list=c(tab.list,templist)
}

bymethod.tablist=tab.list
final.file=paste0(outputpath,"/",out.fname,".Rda")

save(tab.conf,bymethod.tablist,
     file=paste0(outputpath,"/",out.fname,".Rda"))

return(list(tab.conf,bymethod.tablist))
}


liberia_table2b_oneresponse=function(full.comparison,response.name="fam_asb_lt",
                                            tab.digits=2,include.pval=F,
                                     tex.cols=NULL,
                                     pvallab.name="Adj. p-value",pvalvar.name="AdjPvalue",
                                          outputpath,tex.prefix,out.fname,packrows=T,includeconf=F){

  y.fname=tolower(unlist(gsub("[^A-z0-9]","",response.name)))
  itt.file=paste0(outputpath,"/",tex.prefix,"_",y.fname,".tex")
  #order by covariate set
  lvl.sims=levels(full.comparison$Sim)
  order.methods=levels(full.comparison$data)

  full.comparison=full.comparison[(full.comparison$ModelName==response.name),]#

  if(includeconf==F){
    full.comparison=full.comparison[full.comparison$data!="Confidential",]
  }

  n.methods=length(order.methods)+ifelse(includeconf==T,0,-1)
  full.comparison=dplyr::arrange(full.comparison,data)
  full.comparison=dplyr::arrange(full.comparison,Sim)


  #subset columns
  cnames=grepl("ITT|StdErr|Pvalue|overlap|diff|Sim|data",colnames(full.comparison))
  sub.compare.all=full.comparison[,cnames]
  sub.compare.all=sub.compare.all[,c("data",
                                     paste(c("ITT","StdErr",pvalvar.name,"ci.overlap","abs.diff"),"tpassonly",sep="_"),
                                     paste(c("ITT","StdErr",pvalvar.name,"ci.overlap","abs.diff"),"cashassonly",sep="_"),
                                     paste(c("ITT","StdErr",pvalvar.name,"ci.overlap","abs.diff"),"tpcashass",sep="_")
                        )]
  if(include.pval==F){
    sub.compare.all=sub.compare.all[,!grepl("Pvalue",colnames(sub.compare.all))]
    tab.colnames=c("Method",rep(c("ITT","Std. Error","CI overlap","Abs. Diff."),3))
    header.n=4
  }else{
    tab.colnames=c("Method",rep(c("ITT","Std. Error",pvallab.name,"CI overlap","Abs. Diff"),3))
    header.n=5
  }

  if(includeconf==FALSE){
  #confidential table remove nobs, and ci overlap columns
  tab.df=sub.compare.all[sub.compare.all$data!="Confidential",]#colnames(sub.compare.all)!="Sim"]
  }else{
  tab.df=sub.compare.all
}

  rownames(tab.df)=NULL
  #original data table
  tab.itt=knitr::kable(tab.df,format="latex",digits=2,escape=F,booktabs=T,
                        col.names = c(tab.colnames))

  tab.itt=kableExtra::add_header_above(tab.itt,
                                        header=c(" "=1,
                                                 "Therapy Only"=header.n,"Cash Only"=header.n,"Both"=header.n
                                                 ),align="c")

  if(packrows==TRUE){
  starti=1
  for(cov.name in lvl.sims){ #group covariate sets
    pack.n=ifelse(grepl("Full",cov.name)==TRUE,2,n.methods)
    tab.itt=kableExtra::pack_rows(tab.itt,cov.name,starti,starti+pack.n-1)
    starti=starti+pack.n
  }
  }

  if(is.null(tex.cols)==FALSE){
  for(coli in seq(1,length(tex.cols))){
    tex.colspec=paste0(">{\\\\raggedleft\\\\arraybackslash}p{",tex.cols[coli],"\\\\textwidth}")
    if(coli%in%c(1,1+header.n,1+(2*header.n))){
    tex.colspec=paste0(tex.colspec,"|")
    }
    tab.itt=kableExtra::column_spec(tab.itt,column=coli,
                                     latex_valign = "p",border_right = ifelse(coli%in%c(1,1+header.n,1+(2*header.n)),TRUE,FALSE),
                                     width=paste0(tex.cols[coli],"\\\\textwidth"))#,
                                     #latex_column_spec=tex.colspec)
  }
  }


  writeLines(c("\\begin{table}","\\centering","\\small"),itt.file)
  CON=file(itt.file,"a")
  writeLines(tab.itt,CON)
  writeLines(c(
    paste("\\caption{Utility measure for synthetic datasets from the various methods across the three covariate sets for the",
          response.name,"response variable.}",paste0("\\label{tab:apdx-itt-",response.name,"}")),
    "\\end{table}"),CON)
  close(CON)


  print(paste("Finished Utility table for", response.name, "response variable."))


  save(tab.itt,
       file=paste0(outputpath,"/",out.fname,"_",y.fname,".Rda"))
  return(tab.itt)
}

pretable_process=function(full.comparison,n.response=8,adj.out=T,
                        change.DPMb.names=NULL,
                        change.mod.names=NULL){
  lvl.sims=unique(full.comparison$Sim) #Sim Levels (covariate set)
  reduce.name.indic=grepl("Reduce",lvl.sims) #reduced covariate subset indicator
  #order covariate set levels as Reduced, Subset, Full, followed by any other potential sim levels
  lvl.sims=c(lvl.sims[reduce.name.indic],lvl.sims[(!reduce.name.indic)&(grepl("Subset",lvl.sims))],
             lvl.sims[grepl("Full",lvl.sims)],lvl.sims[!grepl("Reduce|Subset|Full",lvl.sims)])
  full.comparison$Sim=factor(as.character(full.comparison$Sim),levels=lvl.sims)

   #change model names
  if(is.null(change.mod.names)==FALSE){
    full.comparison$ModelName=factor(as.character(full.comparison$ModelName),levels=names(change.mod.names))
    levels(full.comparison$ModelName)=unname(sapply(change.mod.names,"[[",1))
  }else{
    full.comparison$ModelName=factor(full.comparison$ModelName,levels=unique(full.comparison$ModelName))
  }

  #change DP Model-based names
  if(is.null(change.DPMb.names)==FALSE){
    for(old.nm in names(change.DPMb.names)){
      full.comparison$data[full.comparison$data==old.nm]=change.DPMb.names[[old.nm]]
    }
    order.methods=c("MV Histogram","Hybrid",sapply(change.DPMb.names,"[[",1))
  }else{
    method.names=unique(full.comparison$data)
    order.methods=c("MV Histogram","Hybrid",method.names[grepl("DP",method.names,ignore.case = T)])
  }
  full.comparison$data=factor(as.character(full.comparison$data),levels=c("Confidential",order.methods))
  full.comparison=dplyr::arrange(full.comparison,ModelName)
  full.comparison=dplyr::arrange(full.comparison,data)
  full.comparison=dplyr::arrange(full.comparison,Sim)

  #adjust p-values
  adjp.in.func=(colnames(full.comparison)[grepl("Pvalue",colnames(full.comparison))])[1]=="AdjPvalue"
  if(adjp.in.func==TRUE){
    bonf.adj=NULL
    pvalvar.name="AdjPvalue"
    pvallab.name="Adj. p-value"
  }else{
    bonf.adj=1
    pvalvar.name="Pvalue"
    pvallab.name="p-value"
    if(adj.out==TRUE){
      adj.fam=3 #for each treatment
      adj.components=3*(n.response-1) #for each treatment and response variable (other than fam_asb_lt)
      pvallab.name="Adj. p-value"
    }
  }


  #ITT result tables
  full.comparison$ModelName=as.character(full.comparison$ModelName)
  #adjust pvalues
  if(adj.out==TRUE&adjp.in.func==FALSE){
    print("Adjusting Pvalue for Multiple comparisons...")
    pvalcol.indic=grepl("Pvalue",colnames(full.comparison))
    full.comparison[full.comparison$ModelName=="fam_asb_lt",pvalcol.indic]=
      full.comparison[full.comparison$ModelName=="fam_asb_lt",pvalcol.indic]*adj.fam
    full.comparison[full.comparison$ModelName!="fam_asb_lt",pvalcol.indic]=
      full.comparison[full.comparison$ModelName!="fam_asb_lt",pvalcol.indic]*adj.components
    for(colidx in which(pvalcol.indic)){
      full.comparison[,colidx]=pmin(rep(1,nrow(full.comparison)),t(unlist(full.comparison[,colidx])))
    }
    pvalvar.name="AdjPvalue"
    colnames(full.comparison)=gsub("^Pvalue","AdjPvalue",colnames(full.comparison))
  }


  full.comparison$ModelName=factor(full.comparison$ModelName,levels=unique(full.comparison$ModelName))
  #add relative and absolute difference
  itt.indic=grepl("ITT",colnames(full.comparison))
  n.itt=sum(itt.indic)

  diff.mat=matrix(rep(NA,nrow(full.comparison)*(n.itt*2)),ncol=n.itt*2)
  colnames(diff.mat)=c(gsub("ITT","abs.diff",colnames(full.comparison[,itt.indic])),
                      gsub("ITT","rel.diff",colnames(full.comparison[,itt.indic])))


  for(cov.set in lvl.sims){ #for each covariate set, get rel and abs diff
    print(paste("Adding Relative and Absolute Difference for", cov.set," covariate set."))
    cov.indic=full.comparison$Sim==cov.set #indicator of covariate set
    conf.cov.indic=((full.comparison$Sim==cov.set)&(full.comparison$data=="Confidential")&(!is.na(full.comparison$ITT_Control)))


    #subset data into confidential and sanitized data
    conf.itt=full.comparison[conf.cov.indic,itt.indic,drop=F]
    san.itt=full.comparison[(full.comparison$data!="Confidential")&cov.indic,itt.indic,drop=F]

    conf.itt=conf.itt[!is.na(conf.itt[,1]),]

    #repeat confidential data to make same number of rows as sanitized data
    multiplicity=nrow(san.itt)/nrow(conf.itt) #number of sanitized methods

    conf.itt=dplyr::bind_rows(base::replicate(multiplicity,conf.itt,simplify=F))
    #absolute difference
    absdiff=abs(san.itt-conf.itt)
    colnames(absdiff)=gsub("ITT","abs.diff",colnames(san.itt))

    #relative difference
    reldiff=absdiff/conf.itt
    colnames(reldiff)=gsub("ITT","rel.diff",colnames(san.itt))

    for(i in seq(1,n.itt)){
      diff.mat[((full.comparison$data!="Confidential")&cov.indic),i]=absdiff[,i]
      diff.mat[((full.comparison$data!="Confidential")&cov.indic),i+n.itt]=reldiff[,i]
    }
  } #end for loop of covariate sets
  full.comparison=cbind(full.comparison,diff.mat)
  return(full.comparison)
}

abs_diff=function(conf.vals,san.vals){
  abs(conf.vals-san.vals)
}

rel_diff=function(conf.vals,san.vals){
  abs_diff(conf.vals,san.vals)/conf.vals
}


