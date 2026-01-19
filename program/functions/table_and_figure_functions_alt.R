heatmap_plot=function(data,value.col,x.col,y.col,scale.colors=NULL,low.val.red=TRUE,
                      plt.title=NULL,value.name=NULL,x.name=NULL,y.name=NULL,facet.col=NULL,
                      base.size=50,reverse.y=F,reverse.x=F,cut.bounds=NULL,log.transform.value=FALSE){
  if(is.null(scale.colors)==TRUE){
    scale.colors=RColorBrewer::brewer.pal(9,"YlOrRd")
    if(low.val.red==FALSE){
      scale.colors=base::rev(scale.colors)
    }
  }
  if(is.null(value.name)==TRUE){
    value.name=colnames(data[,value.col,drop=F])
  }
  if(is.null(x.name)==TRUE){
    x.name=colnames(data[,x.col,drop=F])
  }
  if(is.null(y.name)==TRUE){
    y.name=colnames(data[,y.col,drop=F])
  }
  if(is.null(plt.title)==TRUE){
    plt.title=stringi::stri_trans_totitle(paste("Comparing",value.name, "Across",
                                                colnames(data[,x.col,drop=F]),"and",
                                                colnames(data[,value.col,drop=F])))
  }

  dval=data[,value.col,drop=T]
  if(is.null(cut.bounds)==TRUE){
    if(log.transform.value==TRUE){
      tval=round(log(dval),3)
      data[,value.col,drop=F]=tval
      cut.breaks=c(-Inf,seq(min(tval,na.rm=T),max(tval,na.rm=T),length.out=length(scale.colors)-2),Inf)
        cut.labs=e_breaks_labeller(cut.breaks)
        discrete.values=cut(tval,breaks=cut.breaks,labels=cut.labs)
    }else{
      discrete.values=cut(data[,value.col,drop=T],breaks=length(scale.colors),include.lowest=T)
      }
    }else{
    cut.breaks=seq(min(cut.bounds),max(cut.bounds),length.out=length(scale.colors))
    discrete.values=cut(data[,value.col,drop=T],breaks=cut.breaks,include.lowest=T)
  }
  plot.df=data.frame("xcol"=unlist(data[,x.col,drop=T]),
                     "ycol"=unlist(data[,y.col,drop=T]),
                     "valuecol"=discrete.values)
  if(is.null(facet.col)==FALSE){
    plot.df$facetcol=data[,facet.col,drop=T]
  }

  plot.df=cbind(plot.df,data)
  plot.out=ggplot2::ggplot(plot.df,ggplot2::aes(x=xcol,y=ycol,fill=valuecol))+
    ggplot2::geom_tile()+
    ggplot2::scale_fill_manual(values=scale.colors)+
    ggplot2::labs(x=x.name,y=y.name,fill=value.name)+
    ggplot2::ggtitle(plt.title)+
    ggplot2::theme_minimal(base_size=base.size)
  if(reverse.y==TRUE){
    plot.out=plot.out+ggplot2::scale_y_discrete(limits=rev(levels(plot.df$ycol)))
  }
  if(reverse.x==TRUE){
    plot.out=plot.out+ggplot2::scale_x_discrete(limits=rev(levels(plot.df$xcol)))
  }
  if(is.null(facet.col)==FALSE){
    facet.vals=data[,facet.col,drop=F]
    plot.out=plot.out+ggplot2::facet_grid(~facetcol)
  }
  return(plot.out)
}

heatmap_x_covset=function(data,value.col,x.col,y.col,scale.colors=NULL,low.val.red=TRUE,
                                     value.name=NULL,x.name=NULL,y.name=NULL,facet.col=NULL,
                                     base.size=50,reverse.y=F,reverse.x=F,cut.bounds=NULL,
                           extra.theme=NULL,extra.facet=NULL,
                          log.transform.value=FALSE,
                           treat.vars=c("tpassonly","cashassonly","tpcashass","Control"),
                           treat.names=c("Therapy Only","Cash Only","Both Therapy and Cash","Control Mean"),
                          plot.file.prefix,outputpath,save.individual=F,
                          combo.whs=NULL,individual.whs=NULL,...){
if(is.null(value.name)==TRUE){
  value.name=value.col
}
  if(is.null(combo.whs)==TRUE){
    combo.whs=c(NA,NA,1)
  }
  if(is.null(individual.whs)==TRUE){
    individual.whs=c(NA,NA,1)
  }
  cov.sets=unique(data$Sim)

  df.coln=colnames(data)
  sub.df=data[data$data!="Confidential",
                      grepl(paste(x.col,y.col,value.col,"Sim",sep="|"),df.coln)]
  df.long=tidyr::pivot_longer(sub.df,cols=df.coln[grepl(value.col,df.coln)],
                                 names_to = "Treatment",names_prefix = paste0(value.col,"_"),
                                 values_to = "value")

  if(is.null(treat.names)==FALSE){
  for(i in seq(1,length(treat.vars))){
    df.long$Treatment[df.long$Treatment==treat.vars[i]]=treat.names[i]
  }
    treat.vars=treat.names
  }
  df.long$Treatment=factor(df.long$Treatment,levels=treat.vars)




  df.long.plt=df.long[!grepl("Control",df.long$Treatment),]
  if(sum(grepl("Reduce",cov.sets))>0){
  plot.reduce=heatmap_plot(df.long.plt[grepl("Reduced",df.long.plt$Sim),],
                                   value.col="value",x.col=x.col,y.col=y.col,
                                   reverse.y=reverse.y,reverse.x=reverse.x,
                                   scale.colors=scale.colors,
                                   cut.bounds=cut.bounds,
                                   plt.title=paste("Reduced Subset of Covariates"),
                           log.transform.value=log.transform.value,
                                   value.name=value.name,x.name=x.name,y.name=y.name,base.size=base.size)
  }else{
    plot.reduce=NULL
  }
  if(sum(cov.sets=="Subset")>0){
  plot.subset=heatmap_plot(df.long.plt[df.long.plt$Sim=="Subset",],
                           value.col="value",x.col=x.col,y.col=y.col,
                           reverse.y=reverse.y,reverse.x=reverse.x,
                           scale.colors=scale.colors,
                           cut.bounds=cut.bounds,
                           log.transform.value=log.transform.value,
                           plt.title=paste("Subset of Covariates"),
                           value.name=value.name,x.name=x.name,y.name=y.name,base.size=base.size)
  }else{
    plot.subset=NULL
  }
  if(sum(grepl("Full",cov.sets))>0){
  plot.full=heatmap_plot(df.long.plt[grepl("Full",df.long.plt$Sim),],
                          value.col="value",x.col=x.col,y.col=y.col,
                          reverse.y=reverse.y,reverse.x=reverse.x,
                          scale.colors=scale.colors,
                          cut.bounds=cut.bounds,
                         log.transform.value=log.transform.value,
                          plt.title=paste("Full Baseline Covariates"),
                          value.name=value.name,x.name=x.name,y.name=y.name,base.size=base.size)
  }else{
    plot.full=NULL
  }
  if(is.null(extra.theme)==FALSE){
    plot.reduce=plot.reduce+extra.theme
    plot.subset=plot.subset+extra.theme
    plot.full=plot.full+extra.theme
    }
  if(is.null(extra.facet)==FALSE){
    plot.reduce=plot.reduce+extra.facet
    plot.subset=plot.subset+extra.facet
    plot.full=plot.full+extra.facet
  }

  if(length(cov.sets)==2){
  plt.legend=cowplot::get_legend(plot.subset)

  subplots=cowplot::plot_grid(
    cowplot::plot_grid(plot.reduce+ggplot2::theme(legend.position="none"),plot.subset+ggplot2::theme(legend.position="none"),ncol=1),
    plt.legend,ncol=2,rel_widths=c(1,0.2))

  combo.file=paste0(outputpath,"/",plot.file.prefix,"_heatmap_combo_subsets.png")
  ggplot2::ggsave(combo.file,plot = subplots,
                  bg="white",
                  scale=combo.whs[3],
                  width=combo.whs[1],
                  height=combo.whs[2],...)
  }
  if(length(cov.sets)==3){
    plt.legend=cowplot::get_legend(plot.full)
    subplots=cowplot::plot_grid(
      cowplot::plot_grid(plot.reduce+ggplot2::theme(legend.position="none"),plot.subset+ggplot2::theme(legend.position="none"),ncol=1),
      plt.legend,ncol=2,rel_widths=c(1,0.2))

    combo.file=paste0(outputpath,"/",plot.file.prefix,"_heatmap_combo_subsets.png")
    ggplot2::ggsave(combo.file,plot = subplots,
                    bg="white",
                    scale=combo.whs[3],
                    width=combo.whs[1],
                    height=combo.whs[2],...)
  full.file=paste0(outputpath,"/",plot.file.prefix,"_heatmap_full.png")
  ggplot2::ggsave(full.file,plot = plot.full,
                  bg="white",
                  scale=individual.whs[3],
                  width=individual.whs[1],
                  height=individual.whs[2],
                  ...)

  }
  if(length(cov.sets)==1){
    save.individual=T
  }

  if(save.individual==TRUE){
    if(sum(grepl("Reduce",cov.sets))>0){
    red.file=paste0(outputpath,"/",plot.file.prefix,"_heatmap_reduced.png")
    ggplot2::ggsave(red.file,plot = plot.reduce,
                    bg="white",scale=individual.whs[3],
                    width=individual.whs[1],
                    height=individual.whs[2],
                    ...)
    }
    if(sum(cov.sets=="Subset")>0){
    sub.file=paste0(outputpath,"/",plot.file.prefix,"_heatmap_subset.png")
    ggplot2::ggsave(sub.file,plot = plot.subset,bg="white",
                    scale=individual.whs[3],
                    width=individual.whs[1],
                    height=individual.whs[2],
                    ...)
    }
    if(sum(grepl("Full",cov.sets))>0){
      full.file=paste0(outputpath,"/",plot.file.prefix,"_heatmap_full.png")
      ggplot2::ggsave(full.file,plot = plot.full,
                      bg="white",
                      scale=individual.whs[3],
                      width=individual.whs[1],
                      height=individual.whs[2],
                      ...)
    }
  }

  return(df.long)
}



itt_pval_table=function(data,alpha=0.1,not.signif.color="grey",not.signif.color.name="gray",
                        pvalvar.name="AdjPvalue",pvallab.name="Adj. P-value",tab.digits=2,
                        tex.cols=NULL,
                        outputpath,tex.prefix,out.fname){


  lvl.sims=levels(data$Sim)
  method.names=levels(data$data)
  n.methods=length(method.names)
  mod.names=levels(data$ModelName)

  data$Method=data$data
  data=dplyr::arrange(data,Method)
  data=dplyr::arrange(data,ModelName)
  data=dplyr::arrange(data,Sim)


  #subset to needed columns, exclude control mean statistics
  cnames=grepl("ITT|Pvalue|Model|Method|Sim",colnames(data))
  itt.df=data[,cnames]
  itt.df=itt.df[,!grepl("Control",colnames(itt.df),ignore.case=T)]
  treat.vars=c("tpassonly","cashassonly","tpcashass")
  for(trvar in treat.vars){
    pvals=itt.df[,paste0(c("ITT",pvalvar.name),"_",trvar)]
    colnames(pvals)=c("ITT","Pvalue")
    pvals$Pvalue[is.na(pvals$Pvalue)]=0
    itt.df[,paste0("ITT_",trvar)]= kableExtra::cell_spec(round(itt.df[,paste0("ITT_",trvar)],digits=tab.digits),
                                                    "latex",
                       color = dplyr::case_when(
                         pvals$Pvalue >alpha ~ not.signif.color,
                         .default = "black"
                       ))
  }


  print("after cell_spec")
  covset.tab.list=NULL
  for(covset in lvl.sims){
    print(covset)
    covset.fname=tolower(unlist(gsub("[^A-z0-9]","",covset)))
    temp.file=paste0(outputpath,"/",tex.prefix,"_",covset.fname,".tex")

    sub.df=itt.df[itt.df$Sim==covset,c("Method",paste0("ITT_",treat.vars))]


    pack.n=length(unique(sub.df$Method))#ifelse(grepl("Full",covset),3,n.methods)
    response.col=rep(" ",pack.n*length(mod.names))
    first.y.row=paste0("\\multirow{",pack.n,"}{*}{\\rotatebox{90}{",gsub("_","\\_",mod.names,fixed=T),"}}")
    print(paste("length first.y.row=",length(first.y.row),""))
    response.col[seq(1,length(response.col),pack.n)]=first.y.row
    sub.df$Response=response.col
    sub.df=sub.df[,c("Response","Method",paste0("ITT_",treat.vars))]
    row.names(sub.df)=NULL

    tab.temp= knitr::kable(sub.df,format="latex",
                           booktabs=T,digits=tab.digits,linesep=c("",""),escape = F,
                           row.names = NA,
                           col.names = c("Response","Method","Therapy Only","Cash Only","Both"))


    print("after tab.temp")

    # starti=1
    # for(y.name in mod.names){
    #
    #
    #   tab.temp=kableExtra::pack_rows(tab.temp,y.name,starti,starti+pack.n-1)
    # starti=starti+pack.n
    # }

    if(is.null(tex.cols)==FALSE){
      for(coli in seq(1,length(tex.cols))){
        tex.colspec=paste0(">{\\\\raggedleft\\\\arraybackslash}p{",tex.cols[coli],"\\\\textwidth}")
        tab.temp=kableExtra::column_spec(tab.temp,column=coli,
                                        latex_valign = "p",
                                        width=paste0("{",tex.cols[coli],"\\\\textwidth}"),
                                        latex_column_spec=tex.colspec)
      }
    }



    writeLines(c("\\begin{table}","\\centering"),temp.file)
    CON=file(temp.file,"a")
    writeLines(tab.temp,CON)
    writeLines(c(
      paste("\\caption{The estimated treatment effect for each treatment is shown across the various datasets for each response and the",
      covset,"covariates. Values in \\textcolor{",not.signif.color,"}{",as.character(not.signif.color.name),
      "} are not significant at the $\\alpha=",alpha,"$ level.}"),
      paste0("\\label{tab:itt-pval-",covset.fname,"}"),
      "\\end{table}"),CON)
    close(CON)

    tab.temp=list(tab.temp)
    names(tab.temp)=covset

    covset.tab.list=c(covset.tab.list,tab.temp)

  }

  save(covset.tab.list,itt.df,
       file=paste0(outputpath,"/",out.fname,".Rda"))

  return(covset.tab.list)
}

e_breaks_labeller=function(cut.breaks,include.lowest=F,dig.lab=3){
  lower=round(exp(cut.breaks),dig.lab)[seq(1,length(cut.breaks)-1)]
  upper=lower[-1]
  labs=paste0("(",lower,", ",upper,"]")
  labs=c(paste0("<=",lower[2]),labs[seq(2,length(labs)-1)],paste0(">",upper[length(upper)]))
  labs
}

one_response_treateffects=function(data,reg.model,response=NULL,
                                   treat.vars,treat.names=NULL,scale.colors=NULL,line.sz=2,
                                   plt.title=NULL,base.size=11,
                                   outputpath=NULL,plot.file.prefix="liberia_",
                                   plot.save.size=c(3.5,1.5,1),cov.set="Reduced Subset",only.high.genm=F){

  if(is.null(scale.colors)==T){
    scale.colors=RColorBrewer::brewer.pal(length(treat.vars),"Dark2")
  }

  sub.data=data[,grepl("Model|data|Sim|ITT|StdErr|lower|upper",colnames(data),ignore.case = F)]
  sub.data=sub.data[,grepl(paste0("Model|data|Sim|",paste0(treat.vars,collapse="|")),colnames(sub.data))]
  if(only.high.genm==TRUE){
    keepdata=unique(sub.data$dataPlot)
    keepdata=c(keepdata[1:3],keepdata[length(keepdata)])
    sub.data=sub.data[sub.data$dataPlot%in%keepdata,]
  }
  long.df=tidyr::pivot_longer(sub.data,names_to=c("measure","treatment"),names_sep="_",
                              cols=colnames(sub.data)[grepl("ITT|StdErr|lower|upper",colnames(sub.data),ignore.case = F)],
                              values_to="value")
  pltdata=tidyr::pivot_wider(long.df,values_from = "value",names_from = "measure")
  pltdata$treatment=factor(pltdata$treatment,levels=c("tpassonly","cashassonly","tpcashass"))
  levels(pltdata$treatment)=c("Therapy","Cash","Both")

  if(is.null(response)){
    response=unique(pltdata$ModelName)
  }
  if(is.null(treat.names)==FALSE){
    pltdata$treatment=factor(pltdata$treatment,levels=treat.vars)
    levels(pltdata$treatment)=treat.names
  }

  if(length(response)==1){
  if(is.null(plt.title)==TRUE){
    plt.title=response
  }

  pltdata=pltdata[(pltdata$ModelName==response)&(pltdata$Sim==cov.set),]
  plt=ggplot2::ggplot(data=pltdata,ggplot2::aes(x=dataPlot,y=ITT,color=treatment,ymin=ci.lower,ymax=ci.upper))+
    ggplot2::geom_point()+
    ggplot2::geom_linerange(linewidth=line.sz)+
    ggplot2::scale_color_manual(values=scale.colors)+
    ggplot2::labs(x="Method",y="Est. Effect",fill="Treatment")+
    ggplot2::ggtitle(plt.title)+
    ggplot2::theme_minimal(base_size=base.size)+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=60,hjust=1))+
    ggplot2::facet_wrap(~pltdata$treatment,scales="free")

  plt.file=paste0(outputpath,plot.file.prefix,"_",response,"_",cov.set,"_ci_plot.png")
  ggplot2::ggsave(plt.file,plot = plt,bg="white",
                  scale=plot.save.size[3],
                  width=plot.save.size[1],
                  height=plot.save.size[2])
  }else{
    for(yvar in response){
      if(is.null(plt.title)==TRUE){
        yplt.title=yvar
      }else{
        yplt.title=paste(plt.title,yvar)
      }

      pltdata.y=pltdata[(pltdata$ModelName==yvar)&(pltdata$Sim==cov.set),]
      plt=ggplot2::ggplot(data=pltdata.y,ggplot2::aes(x=dataPlot,y=ITT,ymin=ci.lower,ymax=ci.upper))+
        ggplot2::geom_point()+
        ggplot2::geom_linerange(linewidth=line.sz)+
        ggplot2::scale_color_manual(values=scale.colors)+
        ggplot2::labs(x="Method",y="Est. Effect")+
        ggplot2::ggtitle(yplt.title)+
        ggplot2::theme_minimal(base_size=base.size)+
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle=60,hjust=1))+
        ggplot2::facet_wrap(~pltdata.y$treatment,scales="free")

      plt.file=paste0(outputpath,plot.file.prefix,"_",yvar,"_",cov.set,"_ci_plot.png")
      ggplot2::ggsave(plt.file,plot = plt,bg="white",
                      scale=plot.save.size[3],
                      width=plot.save.size[1],
                      height=plot.save.size[2])
    }
  }
  return(pltdata)
}

# itt_controlmeans_table=function(data,metrics=NULL,metric.colnames=NULL,
#                         pvalvar.name="AdjPvalue",pvallab.name="Adj. P-value",tab.digits=2,
#                         tex.cols=NULL,
#                         outputpath,tex.prefix,out.fname){
#
#
#   lvl.sims=levels(data$Sim)
#   method.names=levels(data$data)
#   n.methods=length(method.names)
#   mod.names=levels(data$ModelName)
#
#   data$Method=data$data
#   data=dplyr::arrange(data,Method)
#   data=dplyr::arrange(data,ModelName)
#   data=dplyr::arrange(data,Sim)
#
#
#   #subset to needed columns, exclude control mean statistics
#   cnames=grepl("Control|Model|Method|Sim",colnames(data))
#   itt.df=data[,cnames]
#
#   if(is.null(metrics)==FALSE){
#     itt.df=itt.df[,grepl(paste0(metrics,"|Model|Method|Sim"),colnames(itt.df))]
#     print(dim(itt.df))
#   }
#
#   if(is.null(metric.colnames)==TRUE){
#     metric.colnames=colnames(itt.df)[grepl(paste0(metrics),colnames(itt.df))]
#     print(metric.colnames)
#   }
#
#   covset.tab.list=NULL
#   for(covset in lvl.sims){
#     covset.fname=tolower(unlist(gsub("[^A-z0-9]","",covset)))
#     temp.file=paste0(outputpath,"/",tex.prefix,"_",covset.fname,".tex")
#
#     sub.df=itt.df[itt.df$Sim==covset,]
#     pack.n=ifelse(grepl("Full",covset),3,n.methods)
#     response.col=rep(" ",pack.n*length(mod.names))
#     first.y.row=paste0("\\multirow{",pack.n,"}{*}{\\rotatebox{90}{",mod.names,"}")
#     response.col[seq(1,length(response.col),pack.n)]=first.y.row
#     sub.df$Response=response.col
#     sub.df=sub.df[,c("Response","Method",colnames(sub.df)[!grepl("Method|Response",colnames(sub.df))])]
#     row.names(sub.df)=NULL
#
#     tab.temp= knitr::kable(sub.df,format="latex",
#                            booktabs=T,digits=tab.digits,linesep=c("",""),escape = F,
#                            row.names = NA,
#                            col.names = c("Response","Method",metric.colnames))
#
#
#
#     # starti=1
#     # for(y.name in mod.names){
#     #
#     #
#     #   tab.temp=kableExtra::pack_rows(tab.temp,y.name,starti,starti+pack.n-1)
#     # starti=starti+pack.n
#     # }
#
#     if(is.null(tex.cols)==FALSE){
#       for(coli in seq(1,length(tex.cols))){
#         tex.colspec=paste0(">{\\\\raggedleft\\\\arraybackslash}p{",tex.cols[coli],"\\\\textwidth}")
#         tab.temp=kableExtra::column_spec(tab.temp,column=coli,
#                                          latex_valign = "p",
#                                          width=paste0("{",tex.cols[coli],"\\\\textwidth}"),
#                                          latex_column_spec=tex.colspec)
#       }
#     }
#
#
#
#     writeLines(c("\\begin{table}","\\centering"),temp.file)
#     CON=file(temp.file,"a")
#     writeLines(tab.temp,CON)
#     writeLines(c(
#       paste("\\caption{The mean of the response for the control group is shown across the various datasets for each response and the",
#             covset,"covariates.}\\label{tab:itt-controlmean-",covset.fname,"}"),
#       "\\end{table}"),CON)
#     close(CON)
#
#     tab.temp=list(tab.temp)
#     names(tab.temp)=covset
#
#     covset.tab.list=c(covset.tab.list,tab.temp)
#
#   }
#
#   save(covset.tab.list,itt.df,
#        file=paste0(outputpath,"/",out.fname,".Rda"))
#
#   return(covset.tab.list)
# }
#


