library(ggplot2)
library(dplyr) #to use %>% command
####### Testing the multivariate histogram method on the Liberia Data
basepath = rprojroot::find_rstudio_root_file()

outputpath=paste0(basepath,"/output")
file.suffix="v0304"

input.plots=list("base.size"=13,"var.plot.size"=c(2300,1500),
                 "ln.sz"=1,"B.plot.size"=c(3500,2700),
                 "err.plot.size"=c(3500,1000))

source(paste0(basepath,"/program/functions/directory_setup.R"))
dir.set=directory_setup(basepath,outputpath,list("figures",paste0("figures/asymptotics_",file.suffix)))

B=5000
resid.var.rel.width=0.25
low.resid.var=c(0.005,0.05,1,1.5)
low.n.obs=c(10,20,30,40)
low.seed=1

high.resid.var=c(5,10,50,100)
high.n.obs=c(15,30,45,60)
high.seed=1

B.vals=c(100,5000,10000)
B.resid.var=c(0.005,1,50)
B.n.obs=c(10,20,30,40)
B.seed=1
B.rel.height=0.1

err.B.vals=c(100,5000,7000,10000)
err.resid.var=c(0.005,1,50)
err.n.obs=seq(5,200,5)
err.replicates=100
err.seed=1

plot.fname.pref=paste0(outputpath,"/figures/asymptotics_",file.suffix,"/")


sim_sample_mse=function(resid.sd=sqrt(resid.var),nobs=n.obs){
  resid=rnorm(nobs,0,resid.sd)
  return(mean(resid^2))
}

sample_mse_data=function(n.vals,B,resid.var){
  sample.mse=dplyr::bind_rows(
    lapply(n.vals,
           function(x)data.frame("n"=rep(x,2*B),
                                 "B"=rep(B,2*B),
                                 "var"=rep(resid.var,2*B),
                                 "dist"=c(rep("MSE sample",B),rep("Normal",B)),
                                 "mse"=c(replicate(B,sim_sample_mse(resid.sd = sqrt(resid.var),nobs=x)),
                                         rnorm(B,resid.var,sqrt(2*(resid.var^2)/x))))))
  sample.mse$ndist=paste(as.character(sample.mse$n),sample.mse$dist)
  sample.mse$ndist=factor(sample.mse$ndist,levels=unlist(
    lapply(as.character(n.vals),function(n)c(paste(n,"MSE sample"),paste(n,"Normal")))))
  return(sample.mse)
}

err_sample_mse=function(n.vals,B,resid.var,sim=1){
  sample.mse=dplyr::bind_rows(
    lapply(n.vals,
           function(x)data.frame("n"=x,
                                 "B"=B,
                                 "var"=resid.var,
                                 "mse"=mean(replicate(B,sim_sample_mse(resid.sd = sqrt(resid.var),nobs=x))),
                                 "sim"=sim)))
  sample.mse$err=sample.mse$mse-sample.mse$var
  return(sample.mse)
}

plot_densities=function(i,df.list=data.list,leg.vect=NULL,ln.sz=10,base.size=11){
  if(is.null(leg.vect)==TRUE){
    leg.vect=rep("right",i)
  }
  #df=df.list[[i]]
  #df.norm=df[df$dist=="Normal",]
  #df.mse=df[df$dist!="Normal",]
  ggplot(data=df.list[[i]],aes(x=mse,color=as.factor(n),linetype=as.factor(dist)))+
    geom_density(fill=NA,show.legend = F)+
    stat_density(geom="line",position="identity",size=0)+
    ggplot2::theme_minimal(base_size=base.size)+
    labs(color="n",linetype="Distribution")+
    scale_color_brewer(palette="Dark2")+
    ggtitle(names(df.list)[i])+
    theme(legend.position = leg.vect[i],axis.title.x = element_blank())+
    guides(color=guide_legend(override.aes = list(size=10)),
           linetype=guide_legend(override.aes = list(size=10)))
}

plot_sample_means=function(err.df,leg.vect=NULL,ln.sz=10,base.size=11,ndigit=3,use.facet=T,use.points=T,use.se=F){
  if(is.null(leg.vect)==TRUE){
    leg.vect="right"
  }
  if(is.character(err.df$var)==FALSE){
    err.df$varLab=paste0("Var=",round(err.df$var,ndigit))
  }else{
    err.df$varLab=err.df$var
  }

  if(is.factor(err.df$B)==FALSE){
    err.df$BLab=paste0("B=",err.df$B)
    err.df$BLab=factor(err.df$BLab,levels=paste0("B=",unique(err.df$B)))
    err.df$B=factor(err.df$B,levels=unique(err.df$B))
  }else{
    err.df$BLab=err.df$B
  }
  plt= ggplot(data=err.df,aes(x=n,y=abs(err),color=B,fill=B))
  if(use.points==T){
    plt=plt+
      geom_point(shape=20)
  }
  plt=plt+
    geom_smooth(se=use.se)+
    ggplot2::theme_minimal(base_size=base.size)+
    labs(x="n",color="B",y="Absolute Error")+
    scale_color_brewer(palette="Dark2")+
    scale_fill_brewer(palette="Set2")+
    ggtitle("Sample Mean Estimator of Variance using MSE")+
    theme(legend.position = leg.vect)
  if(use.facet==TRUE){
    plt=plt+
      facet_grid(varLab~BLab,scale="free")
  }else{
    plt=plt+facet_wrap((~varLab),scale="free")
  }
  plt
}

#https://github.com/wilkelab/cowplot/issues/202
troubleshoot_get_legend <- function(plot, legend = NULL) {

  gt <- ggplotGrob(plot)

  pattern <- "guide-box"
  if (!is.null(legend)) {
    pattern <- paste0(pattern, "-", legend)
  }

  indices <- grep(pattern, gt$layout$name)

  not_empty <- !vapply(
    gt$grobs[indices],
    inherits, what = "zeroGrob",
    FUN.VALUE = logical(1)
  )
  indices <- indices[not_empty]

  if (length(indices) > 0) {
    return(gt$grobs[[indices[1]]])
  }
  return(NULL)
}


get_legend_special=function(df.list,position="right",ln.sz=10,base.size=11,cpal="Dark2"){
  if(position=="bottom"){
    print("here")
    plt1=ggplot(data=df.list[[1]],aes(y=as.numeric(mse),x=as.numeric(mse),linetype=dist))+
      geom_line(size=ln.sz)+
      ggplot2::theme_minimal(base_size=base.size)+
      labs(linetype="Distribution")+
      theme(legend.position = "bottom")
      plt2=ggplot(data=df.list[[1]],aes(y=as.numeric(mse),x=as.numeric(mse),color=as.factor(n)))+
        geom_line(size=ln.sz)+
        ggplot2::theme_minimal(base_size=base.size)+
        labs(color="n")+
        scale_color_brewer(palette=cpal)+
        theme(legend.position = "bottom")
      legend1 = troubleshoot_get_legend(plt1)
      legend2 = troubleshoot_get_legend(plt2)
      leg=cowplot::plot_grid(legend1,legend2,nrow=1)

  }else{
  plt=ggplot(data=df.list[[1]],aes(x=as.numeric(mse),y=as.numeric(mse),color=as.factor(n),linetype=as.factor(dist)))+
    geom_line(size=ln.sz)+
    ggplot2::theme_minimal(base_size=base.size)+
    labs(color="n",linetype="Distribution")+
    scale_color_brewer(palette=cpal)+
    theme(legend.position = position)
  leg=cowplot::get_legend(plt)
  }
  return(leg)
}
## low variance plots
set.seed(low.seed)
low.data.list=lapply(low.resid.var,function(var)sample_mse_data(n.vals=low.n.obs,B=B,resid.var=var))
names(low.data.list)=paste0("Residual variance=",low.resid.var)

ndata=length(low.data.list)
leglow=get_legend_special(df.list=low.data.list,position="right",
                                       ln.sz=input.plots[["ln.sz"]],
                                       base.size = input.plots[["base.size"]])

legend.vector=rep("none",ndata)
cowplot::plot_grid(
  cowplot::plot_grid(
    plotlist=lapply(seq(1,ndata),function(idx)plot_densities(i=idx,
                                                             df.list=low.data.list,
                                                             leg.vect=legend.vector,
                                                             ln.sz=input.plots[["ln.sz"]],
                                                             base.size = input.plots[["base.size"]])),
    ncol=2),
  leglow, ncol=2,rel_widths = c(1,resid.var.rel.width))

ggsave(filename=paste0(plot.fname.pref,"mse_normal_lowvar.png"),
  width=((input.plots[["var.plot.size"]])[1]),
  height=((input.plots[["var.plot.size"]])[2]),
  units="px",bg="white")


## high variance plots

high.data.list=lapply(high.resid.var,function(var)sample_mse_data(n.vals=high.n.obs,B=B,resid.var=var))
names(high.data.list)=paste0("Residual variance=",high.resid.var)

ndata=length(high.data.list)
leghigh=get_legend_special(df.list=high.data.list,position="right",
                       ln.sz=input.plots[["ln.sz"]],
                       base.size = input.plots[["base.size"]])


legend.vector=rep("none",ndata)
cowplot::plot_grid(
  cowplot::plot_grid(
    plotlist=lapply(seq(1,ndata),function(idx)plot_densities(i=idx,
                                                             df.list=high.data.list,
                                                             leg.vect=legend.vector,
                                                             ln.sz=input.plots[["ln.sz"]],
                                                             base.size = input.plots[["base.size"]])),
    ncol=2),
  leghigh, ncol=2,rel_widths = c(1,resid.var.rel.width))

ggsave(filename=paste0(plot.fname.pref,"mse_normal_highvar.png"),
       width=((input.plots[["var.plot.size"]])[1]),
       height=((input.plots[["var.plot.size"]])[2]),
       units="px",bg="white")


## B variance plots
set.seed(B.seed)

B.data.list=list()
for(B in B.vals){
  temp.data.list=lapply(B.resid.var,function(var)sample_mse_data(n.vals=B.n.obs,B=B,resid.var=var))
  names(temp.data.list)=paste0("Residual variance=",B.resid.var,", B=",B)
  B.data.list=c(B.data.list,temp.data.list)

}

ndata=length(B.data.list)
legB=get_legend_special(df.list=B.data.list,position="bottom",
                       ln.sz=input.plots[["ln.sz"]],
                       base.size = input.plots[["base.size"]]*2)



legend.vector=rep("none",ndata)
cowplot::plot_grid(
  cowplot::plot_grid(
    plotlist=lapply(seq(1,ndata),function(idx)plot_densities(i=idx,
                                                             df.list=B.data.list,
                                                             leg.vect=legend.vector,
                                                             ln.sz=input.plots[["ln.sz"]],
                                                             base.size = input.plots[["base.size"]])),
    ncol=3),
  legB, ncol=1,rel_heights = c(1,B.rel.height))

ggsave(filename=paste0(plot.fname.pref,"mse_normal_B.png"),
       width=((input.plots[["B.plot.size"]])[1]),
       height=((input.plots[["B.plot.size"]])[2]),
       units="px",bg="white")


## Get Sample Mean Error


mse_func_for_in_replicate=function(B.v=B,var.v=var,n.v=err.n.obs,sim=1){
  temp.data.list=lapply(B.v,function(B)err_sample_mse(n.vals=n.v,B=B,resid.var=var.v,sim=sim))
  temp.data=dplyr::bind_rows(temp.data.list)
  temp.data
}

set.seed(err.seed)

err.data.list=data.frame()
for(var in err.resid.var){
  temp.data.list=lapply(seq(1,err.replicates),function(i)mse_func_for_in_replicate(B.v=err.B.vals,var.v=var,n.v=err.n.obs,sim=i))
  err.data.list=c(err.data.list,temp.data.list)
}

err.df=dplyr::bind_rows(err.data.list)

plot_sample_means(err.df[err.df$B>150,],leg.vect=NULL,ln.sz=10,base.size=11,ndigit=3,use.facet=F,use.points = F,use.se=T)

ggsave(filename=paste0(plot.fname.pref,"err_mse_sampmean.png"),
       width=((input.plots[["err.plot.size"]])[1]),
       height=((input.plots[["err.plot.size"]])[2]),
       units="px",bg="white")

plot_sample_means(err.df[err.df$B<150,],leg.vect=NULL,ln.sz=10,base.size=11,ndigit=3,use.points = F,use.se=T)

ggsave(filename=paste0(plot.fname.pref,"err_mse_sampmean_Blow.png"),
       width=((input.plots[["err.plot.size"]])[1]),
       height=((input.plots[["err.plot.size"]])[2]),
       units="px",bg="white")



save(high.data.list,low.data.list,
     high.resid.var,high.n.obs,low.resid.var,low.n.obs,
     B.data.list,B.vals,B.resid.var,B.n.obs,
     err.df,err.B.vals,err.resid.var,err.n.obs,err.replicates,
     file=paste0(plot.fname.pref,"mse_normal_data.Rda"))




# colnames(tempdf)="error"
# tempdf$gen=rownames(tempdf)
# tempdf$n=gsub("[^0-9]*","",tempdf$gen)
# tempdf$normal=grepl("Norm",tempdf$gen)
# tempdf$truevar=var
#     tempdf$B=B
#     mean.data.list=c(mean.data.list,list(tempdf))
#   }
# }
#
# data.list=list()
# mean.data.list=list()
# plot.list=list()
# for(B in B.vals){
# for(var in resid.var){
# set.seed(1)
#
#
# sample.mse=dplyr::bind_rows(
#   lapply(n.obs,
#          function(x)data.frame("n"=c(rep(as.character(x),B),
#                                      rep(paste0("Norm",as.character(x)),B)),
#                                "mse"=c(replicate(B,sim_sample_mse(resid.sd = sqrt(var),nobs=x)),
#                                        rnorm(B,var,sqrt(2*(var^2)/x))))))
# sample.mse$n=factor(sample.mse$n,levels=unlist(
#   lapply(as.character(n.obs),function(n)c(n,paste0("Norm",n)))))
# past.names=names(data.list)
# data.list=c(data.list,list(sample.mse))
# names(data.list)=c(past.names,paste0("Var=",var,", B=",B))
#
# tempdf=data.frame(tapply(sample.mse$mse,sample.mse$n,function(x)mean(x)-var))
# colnames(tempdf)="error"
# tempdf$gen=rownames(tempdf)
# tempdf$n=gsub("[^0-9]*","",tempdf$gen)
# tempdf$normal=grepl("Norm",tempdf$gen)
# tempdf$truevar=var
# tempdf$B=B
# mean.data.list=c(mean.data.list,list(tempdf))
# }
# }
#
# err.df=dplyr::bind_rows(mean.data.list)
# plot.list=lapply(seq(1,length(data.list)),function(i)ggplot(data=data.list[[i]],aes(x=mse,color=n))+
#                    geom_density(fill=NA)+
#                    theme_bw()+
#                    scale_color_brewer(palette="Paired")+
#                    ggtitle(names(data.list)[i])+
#                    theme(legend.position = ifelse(i==length(data.list),"bottom","none")))
#
# cowplot::plot_grid(plotlist=plot.list,ncol=3,rel_heights =  = c(rep(1,length(B.vals)-1),1.2))
# rownames(err.df)=NULL
#
#
# library(ggplot2)
# library(cowplot)
# library(dplyr)
# n.obs=seq(10,100,10)
# B.vals=c(50,seq(100,5000,100))
# resid.var=c(0.05)
#
# data.list=list()
# mean.data.list=list()
# plot.list=list()
# for(B in B.vals){
#   for(var in resid.var){
#     set.seed(1)
#     sim_sample_mse=function(resid.sd=sqrt(resid.var),nobs=n.obs){
#       resid=rnorm(nobs,0,resid.sd)
#       return(mean(resid^2))
#     }
#     sample.mse=dplyr::bind_rows(
#       lapply(n.obs,
#              function(x)data.frame("n"=c(rep(as.character(x),B),
#                                          rep(paste0("Norm",as.character(x)),B)),
#                                    "mse"=c(replicate(B,sim_sample_mse(resid.sd = sqrt(var),nobs=x)),
#                                            rnorm(B,var,sqrt(2*(var^2)/x))),
#                                    "B"=rep(B,2*B),
#                                    "var"=rep(var,2*B))))
#     sample.mse$n=factor(sample.mse$n,levels=unlist(
#       lapply(as.character(n.obs),function(n)c(n,paste0("Norm",n)))))
#     past.names=names(data.list)
#     data.list=c(data.list,list(sample.mse))
#     names(data.list)=c(past.names,paste0("Var=",var,", B=",B))
#
#     tempdf=data.frame(tapply(sample.mse$mse,sample.mse$n,function(x)mean(x)-var))
#     colnames(tempdf)="error"
#     tempdf$gen=rownames(tempdf)
#     tempdf$n=gsub("[^0-9]*","",tempdf$gen)
#     tempdf$normal=grepl("Norm",tempdf$gen)
#     tempdf$truevar=var
#     tempdf$B=B
#     mean.data.list=c(mean.data.list,list(tempdf))
#   }
# }
# full.df=dplyr::bind_rows(data.list)
# full.df$normal=ifelse(grepl("Norm",full.df$n)==T,"NormalApprox","MSEsample")
# full.df$diff.from.true=full.df$mse-full.df$var
# full.df$n=gsub("[^0-9]","",full.df$n)
# full.df=full.df%>%select(-var)%>%pivot_wider(id_cols=c(B,n),names_from = normal,values_from = c(mse,diff.from.true))#%>%mutate(diff.from.normapprox=)
# head(full.df)
#
# err.df=dplyr::bind_rows(mean.data.list)
# err.df$NormalLab=ifelse(err.df$normal==TRUE,"Gaussian Appoximation","MSE Samples")
# rownames(err.df)=NULL
# ggplot(err.df,aes(y=error,x=B,color=as.character(n)))+
#   geom_point()+
#   geom_smooth(se=F)+
#   theme_bw()+
#   scale_color_brewer(palette="Paired")+
#   labs(y="Samp.Mean- True Var",x="Num. Proxies",color="n observation")+
#   facet_wrap(~NormalLab,scales="free")
#
# head(err.df)
# err.df.wide=err.df%>%select(n,B,error,truevar,NormalLab)%>%pivot_wider(values_from = error,names_from = NormalLab)
# colnames(err.df.wide)=c(colnames(err.df.wide)[seq(1,length(colnames(err.df.wide))-2)],"MSEsampErr","NormApproxErr")
# err.df.wide$diff.approx.mse=err.df.wide$MSEsampErr-err.df.wide$NormApproxErr
# head(err.df.wide)
# ggplot(err.df.wide,aes(y=diff.approx.mse,x=B,color=as.character(n)))+
#   geom_point()+
#   geom_smooth(se=F)+
#   theme_bw()+
#   scale_color_brewer(palette="Paired")+
#   labs(y="MSE.Samp.Mean- Norm.Approx.Samp.Mean",x="Num. Proxies",color="n observation")
#
#
# n.obs=seq(10,150,2)
# B.vals=c(10,50, 100,500,1000,5000)#seq(100,5000,100))
# resid.var=c(0.05)
#
# data.list=list()
# mean.data.list=list()
# plot.list=list()
# for(B in B.vals){
#   for(var in resid.var){
#     set.seed(1)
#     sim_sample_mse=function(resid.sd=sqrt(resid.var),nobs=n.obs){
#       resid=rnorm(nobs,0,resid.sd)
#       return(mean(resid^2))
#     }
#     sample.mse=dplyr::bind_rows(
#       lapply(n.obs,
#              function(x)data.frame("n"=c(rep(as.character(x),B),
#                                          rep(paste0("Norm",as.character(x)),B)),
#                                    "mse"=c(replicate(B,sim_sample_mse(resid.sd = sqrt(var),nobs=x)),
#                                            rnorm(B,var,sqrt(2*(var^2)/x))),
#                                    "B"=rep(B,2*B),
#                                    "var"=rep(var,2*B))))
#     sample.mse$n=factor(sample.mse$n,levels=unlist(
#       lapply(as.character(n.obs),function(n)c(n,paste0("Norm",n)))))
#     past.names=names(data.list)
#     data.list=c(data.list,list(sample.mse))
#     names(data.list)=c(past.names,paste0("Var=",var,", B=",B))
#
#     tempdf=data.frame(tapply(sample.mse$mse,sample.mse$n,function(x)mean(x)-var))
#     colnames(tempdf)="error"
#     tempdf$gen=rownames(tempdf)
#     tempdf$n=gsub("[^0-9]*","",tempdf$gen)
#     tempdf$normal=grepl("Norm",tempdf$gen)
#     tempdf$truevar=var
#     tempdf$B=B
#     mean.data.list=c(mean.data.list,list(tempdf))
#   }
# }
# full.df=dplyr::bind_rows(data.list)
# full.df$normal=ifelse(grepl("Norm",full.df$n)==T,"NormalApprox","MSEsample")
# full.df$diff.from.true=full.df$mse-full.df$var
# full.df$n=gsub("[^0-9]","",full.df$n)
# full.df=full.df%>%select(-var)%>%pivot_wider(id_cols=c(B,n),names_from = normal,values_from = c(mse,diff.from.true))#%>%mutate(diff.from.normapprox=)
# head(full.df)
#
# err.df=dplyr::bind_rows(mean.data.list)
# err.df$NormalLab=ifelse(err.df$normal==TRUE,"Gaussian Appoximation","MSE Samples")
# rownames(err.df)=NULL
# ggplot(err.df,aes(y=error,x=as.numeric(n),color=as.character(B)))+
#   geom_point()+
#   geom_smooth(se=F)+
#   theme_bw()+
#   scale_color_brewer(palette="Paired")+
#   labs(y="Samp.Mean- True Var",x="n observations",color="Num. Proxies")+
#   facet_wrap(B~NormalLab,scales="free")
#
# ggplot(err.df,aes(y=error,x=as.numeric(n),color=NormalLab))+
#   geom_point()+
#   geom_smooth(se=F)+
#   theme_bw()+
#   scale_color_brewer(palette="Paired")+
#   labs(y="Samp.Mean- True Var",x="n observations",color=" ")+
#   facet_wrap(~B,scales="free")
#
# head(err.df)
# err.df.wide=err.df%>%select(n,B,error,truevar,NormalLab)%>%pivot_wider(values_from = error,names_from = NormalLab)
# colnames(err.df.wide)=c(colnames(err.df.wide)[seq(1,length(colnames(err.df.wide))-2)],"MSEsampErr","NormApproxErr")
# err.df.wide$diff.approx.mse=err.df.wide$MSEsampErr-err.df.wide$NormApproxErr
# head(err.df.wide)
# ggplot(err.df.wide,aes(y=diff.approx.mse,x=as.numeric(n),color=as.character(B)))+
#   geom_point()+
#   geom_smooth(se=F)+
#   theme_bw()+
#   scale_color_brewer(palette="Paired")+
#   labs(y="MSE.Samp.Mean- Norm.Approx.Samp.Mean",x="n observations",color="Num.Proxies")#+
#   #facet_wrap(~B,scales="free")
