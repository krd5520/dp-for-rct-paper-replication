
heatmap_plot=function(data,value.col,x.col,y.col,scale.colors=NULL,low.val.red=TRUE,
                      plt.title=NULL,value.name=NULL,x.name=NULL,y.name=NULL,base.size=50){
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
  discrete.values=cut(data[,value.col,drop=T],breaks=length(scale.colors)-1)
  plot.df=data.frame("x"=unlist(data[,x.col,drop=T]),"y"=unlist(data[,y.col,drop=T]),"value"=discrete.values)
  plot.out=ggplot2::ggplot(data,ggplot2::aes(x=x,y=y,fill=value))+
    ggplot2::geom_tile()+
    ggplot2::scale_fill_manual(values=scale.colors)+
    ggplot2::labs(x=x.col,y=y.col,fill=value.col)+
    ggplot2::ggtitle(plt.title)+
    ggplot2::theme_minimal(base_size=base.size)
  return(plot.out)
}

itt_pval_table=function(data,alpha=0.05,adj.out=TRUE,
                        significance.colors=c("white","grey"),mod.names=NULL){
  #if(is.null(mod.names)==FALSE){
  #data$ModelName=factor(data$ModelName,levels=mod.names)
  #}else{
  #  data$ModelName=as.factor(data$Model)
  #}
  #adjp.in.func=(colnames(data)[grepl("Pvalue",colnames(data))])[1]=="AdjPvalue"
  #if(adjp.in.func==TRUE){
  #  bonf.adj=NULL
  #  pvalvar.name="AdjPvalue"
  # pvallab.name="Adj. p-value"
  #}else{
  #  bonf.adj=1
  #  pvalvar.name="Pvalue"
  #  pvallab.name="p-value"
  #  if(adj.out==TRUE){
   #   adj.fam=9
  #    adj.components=3*7
  #    pvallab.name="Adj. p-value"
  #  }
  #}
#
#   if(adj.out==TRUE&adjp.in.func==FALSE){
#     pvalcol.indic=grepl("Pvalue",colnames(data))
#     data[data$ModelName=="fam_asb_lt",pvalcol.indic]=
#       data[data$ModelName=="fam_asb_lt",pvalcol.indic]*adj.fam
#     data[data$ModelName!="fam_asb_lt",pvalcol.indic]=
#       data[data$ModelName!="fam_asb_lt",pvalcol.indic]*adj.components
#     for(colidx in which(pvalcol.indic)){
#       data[,colidx]=pmin(rep(1,nrow(data)),t(unlist(data[,colidx])))
#     }
#   }


  cnames=grepl("ITT|Pvalue|Model",colnames(data))
  itt.df=data[,cnames]
  itt.df=itt.df[,!grepl("control",colnames(itt.df))]
  itt.df=tidyr::pivot_longer(itt.df,
                             cols=colnames(ittdf)[grepl("ITT|Pvalue")],
                             names_to = "Metric_Treatment")
  splt.cols=base::strsplit(itt.df$Metric_Treatment,"_")
  itt.df$metric=sapply(splt.cols,"[[",1)
  itt.df$treatment=sapply(splt.cols,"[[",2)
  pval.df=itt.df[grepl("Pvalue",itt.df$metric),]
  conf.significant=unlist(pval.df[pval.df$Method=="Confidential",pvalvar.name,drop=T])<=alpha
  itt.df=itt.df[!grepl("Pvalue",itt.df$metric),]
  itt.df=tidyr::pivot_wider(id_cols=c(ModelName,Sim),names_from =c(treatment,Method),values_from=value)
  TF.significant=rep(conf.significant,each=12)

  col.order=paste
  #sub.compare.all$ModelName=paste0(sub.compare.all$ModelName,ifelse(sub.compare$nObs<946,"\\dagger",""))
  sub.compare.san=sub.compare.all[sub.compare.all$data!="Confidential",]
  sub.compare.conf=sub.compare.all[sub.compare.all$data=="Confidential",]
  tab.df.conf=sub.compare.conf[,c(1,3,4,5,6,8,9,10,12,13,14)]

}
  itt.tab=kableExtra::row_spec(itt.tab,1,color="black",
                                  background=ifelse(df$reduced.subset==TRUE,"redsubsetcol",ifelse(df$subset==TRUE,"subsetcol","white")))
