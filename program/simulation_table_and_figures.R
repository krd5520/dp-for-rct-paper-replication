devtools::install_github("krd5520/DPrct")
library(dplyr) #to use %>% command
####### Testing the multivariate histogram method on the Liberia Data
basepath = rprojroot::find_rstudio_root_file()
main.start=proc.time()

outputpath=paste0(basepath,"/output")
file.suffix="0320_loweffect"
truetreat=-0.25
#file.suffix="0320_5effect"
truetreat=5
sim2.conf.suffix="conf_sim_models"
save.plot.png=T

sim1.plots=list("y"="treff","x"="Budget",
                "facet"="wrap","methods"="all","outliers"="keep",
                "base.size"=13,"plot.size"=c(2000,12000))
DPMbplot.name=c("DP-Mb(1)"="GenM-1",
                "DP-Mb(2)"="GenM-2")
DPMbtab.name=c("DP-Mb(1)"="\\dpmb{1}",
               "DP-Mb(2)"="\\dpmb{2}")
longDPMb="GenModel"

plot.file.identifier=c("sim1_treff_BudgetX_AllMethods_KeepOutliers",
                       "sim2_CIoverlap_ModelsX_AllMethods_KeepOutliers",
                       "sim2_treff_ModelX_NotDPMb_KeepOutliers",
                       "sim2_treff_ModelX_OnlyDPMb_RemoveOutliers",
                       "sim2_treff_ModelX_OnlyDPMb_KeepOutliers",
                       "sim2_CIoverlap_MethodsX",
                       "sim2_treff_MethodsX_KeepOutliers")
table.digits=2
sim1.base.size=13
sim2.base.size=13
plot.sizes=list("row1col3"=c(2000,1000), #sim1 plot
                "row1col2"=c(1500,1000), #ModelX OnlyDPmb and  ModelX NoDPmb
                "row2col2"=c(2000,2000), #ModelX All Methods
                "row3col3"=c(1500,1500),
                "row1col1"=c(1000,1000))
plot.dims=c("row1col3","row2col2","row1col2","row1col2","row1col2","row3col3","row3col3")



plot.file.names=paste0(outputpath,"/figures/simplots_",file.suffix,"/",plot.file.identifier,".png")
plot.file.names.nonoise=paste0(outputpath,"/figures/simplots_",file.suffix,"/nonoise",plot.file.identifier,".png")

plot.fname.pref=paste0(outputpath,"/figures/simplots_",file.suffix,"/")
plot.file.names.nonoise=paste0(outputpath,"/figures/simplots_",file.suffix,"/nonoise",plot.file.identifier,".png")

source(paste0(basepath,"/program/functions/directory_setup.R"))
dir.set=directory_setup(basepath,outputpath,list("tables","figures","simulations",paste0("figures/simplots_",file.suffix)))


#get ylims
get_ylims=function(value,group,quantiles.val=c(0.1,0.9),center=NA){
  quants=tapply(value,group,function(x)quantile(x,quantiles.val))
  if(is.na(center)==F){
  ylim.diff=ceiling(max(abs(center-unlist(quants)))/0.25)*0.25
  min.y=center-ylim.diff
  max.y=center+ylim.diff
  }else{
    min.y=floor(min(abs(unlist(quants)))*0.25)/0.25
    max.y=ceiling(max(abs(unlist(quants)))/0.25)*0.25
  }
  return(c(min.y,max.y))
}
###############################
############### Sim 1 ##########
#################################
################### Tables ###########################
#load simulation data for privacy budgets
load(paste0(outputpath,"/simulations/results_budgets_",file.suffix,".Rda"))

load(paste0(outputpath,"/simulations/simulation1_",file.suffix,".Rda"))
#load(paste0(outputpath,"/simulations/sims_",file.suffix,"/results_budgets_",file.suffix,".Rda"))

summary.tab.treff$Algorithm[grepl("DP",summary.tab.treff$Algorithm,ignore.case = T)]="\\longdpmb"
col.order=c("Epsilon","Delta","Algorithm","Avg.ITT",#"Median.ITT",
            "Avg.CI.overlap",#"Avg.Abs.Err",
            "Median.Abs.Err")#,"Max.StdErr")
conf.row=grepl("Confidential",summary.tab.treff$Method)
summary.tab.treff.table=rbind(summary.tab.treff[conf.row,col.order],summary.tab.treff[!conf.row,col.order])
#copy and paste the output of the following into latex to get table
# OR use R markdown with latex capabilities and chunk output="asis"
sim1tab=knitr::kable(summary.tab.treff.table,
             col.names = c("\\epsilon","\\delta","Algorithm",
                           "Mean $\\treff$ Estimate",#"Median $\\treff$ Estimate",
                           "Mean $\\cioverlap$", #"Mean Abs. Diff.",
                           "Median Abs. Diff."),
                          # "Max Std. Err."),
             digits=table.digits,format="latex",booktabs=T,
             linesep=c("\\addlinespace",rep(c("","","\\addlinespace"),7)),
             escape=F)
sim1tab
writeLines(sim1tab,paste0(outputpath,"/tables/sim1_budget_",file.suffix,".tex"))


summary.tab.treff.table2=tidyr::pivot_wider(summary.tab.treff[!conf.row,col.order],names_from = Algorithm,values_from = c(Avg.ITT,Avg.CI.overlap,Median.Abs.Err))
summary.tab.treff.table2$Epsilon=as.character(summary.tab.treff.table2$Epsilon)
summary.tab.treff.table2$Delta[summary.tab.treff.table2$Delta!="0"]="0.001"
new.col.order=c("Epsilon","Delta",
                paste0(c("Avg.ITT_","Avg.CI.overlap_","Median.Abs.Err_"),"MV Histogram"),
                paste0(c("Avg.ITT_","Avg.CI.overlap_","Median.Abs.Err_"),"Hybrid"),
                paste0(c("Avg.ITT_","Avg.CI.overlap_","Median.Abs.Err_"),"\\longdpmb"))
#copy and paste the output of the following into latex to get table
# OR use R markdown with latex capabilities and chunk output="asis"
sim1tab2=knitr::kable(summary.tab.treff.table2[,new.col.order],
                     col.names = c("$\\epsilon$","$\\delta$",
                                   rep(c("Mean $\\widetilde{\\treff}$",
                                         "Mean $\\cioverlap$","Median $\\absdiff$"),3)),
                     # "Max Std. Err."),
                     digits=table.digits,format="latex",booktabs=T,
                     #linesep=rep(c("","","\\addlinespace"),7),
                     escape=F)
sim1tab2=kableExtra::add_header_above(sim1tab2,c("","","MV Histogram"=3,"Hybrid"=3,"\\longdpmb"=3))
sim1tab2
writeLines(sim1tab2,paste0(outputpath,"/tables/sim1_budget_wide_",file.suffix,".tex"))

treat.itt=ITT1%>%filter(Treatment!="control"&!grepl("Approx",Model))
treat.itt$MethMod=paste0(treat.itt$Method,"_",treat.itt$Model)
treat.itt.notdpmb=treat.itt%>%filter(grepl("Confidential|Hybrid|Histogram",Method))
treat.itt.onlydpmb=treat.itt%>%filter(!grepl("Confidential|Hybrid|Histogram",Method))

temp.notdpmb=dplyr::bind_rows(tapply(treat.itt.notdpmb$StdErr,treat.itt.notdpmb$MethMod,summary))
min(temp.notdpmb$Min.)
max(temp.notdpmb$Max.)
temp.onlydpmb=dplyr::bind_rows(tapply(treat.itt.onlydpmb$StdErr,treat.itt.onlydpmb$MethMod,summary))
min(temp.onlydpmb$Min.)
max(temp.onlydpmb$Max.)
tapply(treat.itt.onlydpmb$StdErr,treat.itt.onlydpmb$MethMod,summary)
################### Plots ###########################
plot.idx=1
sim1.plots=list("y"="treff","x"="Budget",
                "facet"="wrap","methods"="all","outliers"="keep",
                "base.size"=13,"plot.size"=c(2000,12000))


plotdata.itt=ITT1%>%filter(Treatment!="control")%>%
  filter(Method!="Confidential")%>%
  filter(grepl("Pure",Model))
plotdata.itt$Method=factor(plotdata.itt$Method,
                           levels=c("Confidential","MV Histogram","Hybrid","DP Model-based"),
                           labels=c("Confidential","MV Histogram","Hybrid",longDPMb))
dta.nrow=nrow(plotdata.itt)
min.x=floor(min(plotdata.itt$Epsilon)-0.5)
max.x=ceiling(max(plotdata.itt$Epsilon)+0.5)
sim1plot=ggplot2::ggplot(data=plotdata.itt,
                ggplot2::aes(y=ITT,x=Epsilon,group=Epsilon))+
  ggplot2::geom_boxplot()+
  ggplot2::geom_line(
    ggplot2::aes(x=seq(min.x,max.x,length.out=dta.nrow),y=rep(truetreat,dta.nrow),
                 colour="True Treatment Effect"))+
  ggplot2::geom_line(
    ggplot2::aes(x=seq(min.x,max.x,length.out=dta.nrow),
                 y=rep(unlist(ITT1$ITT[(ITT1$Treatment!="control")&(ITT1$Method=="Confidential")]),dta.nrow/length(unique(ITT1$Sim))),
                 colour="Confidential Estimated Treatment Effect"),linetype="dashed")+
  ggplot2::labs(x=expression("Privacy Budget"~(epsilon)),
                y="Treatment Effect",color=" ")+
  ggplot2::theme_minimal(base_size=sim1.base.size)+
  ggplot2::theme(legend.position = "bottom")+
  ggplot2::facet_wrap(~plotdata.itt$Method)
sim1plot
#if(save.plot.png==T){
ggplot2::ggsave(filename=paste0(plot.fname.pref,"sim1_treff_BudgetX_NotMVHist_KeepOutliers.png"),
                width=((plot.sizes[["row1col3"]])[1]),
                height=((plot.sizes[["row1col3"]])[2]),
                units="px",bg="white")

for(m in unique(plotdata.itt$Method)){
  tempdf=plotdata.itt[plotdata.itt$Method==m,]
  ntemp=nrow(tempdf)
  sim1plot=ggplot2::ggplot(data=tempdf,
                           ggplot2::aes(y=ITT,x=Epsilon,group=Epsilon))+
    ggplot2::geom_boxplot()+
    ggplot2::geom_line(
      ggplot2::aes(x=seq(min.x,max.x,length.out=ntemp),y=rep(truetreat,ntemp),
                   colour="True Treatment Effect"))+
    ggplot2::geom_line(
      ggplot2::aes(x=seq(min.x,max.x,length.out=ntemp),
                   y=rep(unlist(ITT1$ITT[(ITT1$Treatment!="control")&(ITT1$Method=="Confidential")]),ntemp/length(unique(ITT1$Sim))),
                   colour="Confidential Estimated Treatment Effect"),linetype="dashed")+
    ggplot2::labs(x=expression("Privacy Budget"~(epsilon)),
                  y="Treatment Effect",color=" ")+
    ggplot2::theme_minimal(base_size=sim1.base.size)+
    ggplot2::theme(legend.position = "bottom")+
    ggplot2::ggtitle(m)
  sim1plot
  ggplot2::ggsave(filename=paste0(plot.fname.pref,"sim1_treff_BudgetX_",m,"_KeepOutliers.png"),
                  width=((plot.sizes[["row1col1"]])[1]),
                  height=((plot.sizes[["row1col1"]])[2]),
                  units="px",bg="white")


}
ylims.plt1=get_ylims(plotdata.itt$ITT,plotdata.itt$Method)
sim1plot.noout=ggplot2::ggplot(data=plotdata.itt,
                         ggplot2::aes(y=ITT,x=Epsilon,group=Epsilon))+
  ggplot2::geom_boxplot(outlier.shape=NA)+
  ggplot2::geom_line(
    ggplot2::aes(x=seq(min.x,max.x,length.out=dta.nrow),y=rep(truetreat,dta.nrow),
                 colour="True Effect"))+
  ggplot2::geom_line(
    ggplot2::aes(x=seq(min.x,max.x,length.out=dta.nrow),
                 y=rep(unlist(ITT1$ITT[(ITT1$Treatment!="control")&(ITT1$Method=="Confidential")]),dta.nrow/length(unique(ITT1$Sim))),
                 colour="Confidential Estimate"),linetype="dashed")+
  ggplot2::labs(x=expression("Privacy Budget"~(epsilon)),
                y="Treatment Effect",color=" ")+
  ggplot2::theme_minimal(base_size=sim1.base.size)+
  ggplot2::theme(legend.position = "bottom")+
  ggplot2::scale_y_continuous(limits = ylims.plt1)+
  ggplot2::facet_wrap(~plotdata.itt$Method)
sim1plot.noout
#if(save.plot.png==T){
ggplot2::ggsave(filename=paste0(plot.fname.pref,"sim1_treff_BudgetX_AllMethods_RemoveOutliers.png"),
                width=((plot.sizes[["row1col3"]])[1]),
                height=((plot.sizes[["row1col3"]])[2]),
                units="px",bg="white")
#}
plot.idx=plot.idx+1

nomv.pltdata=plotdata.itt[plotdata.itt$Method!="MV Histogram",]
ylims.plt1=get_ylims(nomv.pltdata$ITT,nomv.pltdata$Method,center=truetreat)
dta.nrow=nrow(nomv.pltdata)
sim1plot.noout.nomv=ggplot2::ggplot(data=nomv.pltdata,
                               ggplot2::aes(y=ITT,x=Epsilon,group=Epsilon))+
  ggplot2::geom_boxplot(outlier.shape=NA)+
  ggplot2::geom_line(
    ggplot2::aes(x=seq(min.x,max.x,length.out=dta.nrow),y=rep(truetreat,dta.nrow),
                 colour="True Effect"))+
  ggplot2::geom_line(
    ggplot2::aes(x=seq(min.x,max.x,length.out=dta.nrow),
                 y=rep(unlist(ITT1$ITT[(ITT1$Treatment!="control")&(ITT1$Method=="Confidential")]),dta.nrow/length(unique(ITT1$Sim))),
                 colour="Confidential Estimate"),linetype="dashed")+
  ggplot2::labs(x=expression("Privacy Budget"~(epsilon)),
                y="Treatment Effect",color=" ")+
  ggplot2::theme_minimal(base_size=sim1.base.size)+
  ggplot2::theme(legend.position = "bottom")+
  ggplot2::scale_y_continuous(limits = ylims.plt1)+
  ggplot2::facet_wrap(~nomv.pltdata$Method)
sim1plot.noout.nomv
#if(save.plot.png==T){
ggplot2::ggsave(filename=paste0(plot.fname.pref,"sim1_treff_BudgetX_NoMVHist_RemoveOutliers.png"),
                width=((plot.sizes[["row1col2"]])[1]),
                height=((plot.sizes[["row1col2"]])[2]),
                units="px",bg="white")
#}
plot.idx=plot.idx+1

onlymv.pltdata=plotdata.itt[plotdata.itt$Method=="MV Histogram",]
ylims.plt1=get_ylims(onlymv.pltdata$ITT,onlymv.pltdata$Method,center=truetreat)
dta.nrow=nrow(onlymv.pltdata)
sim1plot.noout.onlymv=ggplot2::ggplot(data=onlymv.pltdata,
                                    ggplot2::aes(y=ITT,x=Epsilon,group=Epsilon))+
  ggplot2::geom_boxplot(outlier.shape=NA)+
  ggplot2::geom_line(
    ggplot2::aes(x=seq(min.x,max.x,length.out=dta.nrow),y=rep(truetreat,dta.nrow),
                 colour="True Effect"))+
  ggplot2::geom_line(
    ggplot2::aes(x=seq(min.x,max.x,length.out=dta.nrow),
                 y=rep(unlist(ITT1$ITT[(ITT1$Treatment!="control")&(ITT1$Method=="Confidential")]),dta.nrow/length(unique(ITT1$Sim))),
                 colour="Confidential Estimate"),linetype="dashed")+
  ggplot2::labs(x=expression("Privacy Budget"~(epsilon)),
                y="Treatment Effect",color=" ")+
  ggplot2::theme_minimal(base_size=sim1.base.size)+
  ggplot2::theme(legend.position = "none")+
  ggplot2::scale_y_continuous(limits = c(ylims.plt1[1],5.5))+
  ggplot2::ggtitle("MV Histogram")
#  ggplot2::facet_wrap(~nomv.pltdata$Method)
sim1plot.noout.onlymv
#if(save.plot.png==T){
ggplot2::ggsave(filename=paste0(plot.fname.pref,"sim1_treff_BudgetX_OnlyMVHist_RemoveOutliers.png"),
                width=((plot.sizes[["row1col1"]])[1]),
                height=((plot.sizes[["row1col1"]])[2]),
                units="px",bg="white")
#}
plot.idx=plot.idx+1


ylims.plt1=get_ylims(nomv.pltdata$Abs.Err,nomv.pltdata$Method,center=0)
abserr.noout.nomv=ggplot2::ggplot(data=nomv.pltdata,
                                    ggplot2::aes(y=Abs.Err,x=Epsilon,group=Epsilon))+
  ggplot2::geom_boxplot(outlier.shape=NA)+
  ggplot2::geom_line(
    ggplot2::aes(x=seq(min.x,max.x,length.out=nrow(nomv.pltdata)),y=rep(0.2,nrow(nomv.pltdata)),
                 colour="Confidential Effect Error"))+
  ggplot2::labs(x=expression("Privacy Budget"~(epsilon)),
                y="Absolute Error",color=" ")+
  ggplot2::theme_minimal(base_size=sim1.base.size)+
  ggplot2::theme(legend.position = "bottom")+
  ggplot2::scale_y_continuous(limits = c(0,ylims.plt1[2]))+
  ggplot2::facet_wrap(~nomv.pltdata$Method)
abserr.noout.nomv
#if(save.plot.png==T){
ggplot2::ggsave(filename=paste0(plot.fname.pref,"sim1_AbsErr_BudgetX_NoMVHist_RemoveOutliers.png"),
                width=((plot.sizes[["row1col2"]])[1]),
                height=((plot.sizes[["row1col2"]])[2]),
                units="px",bg="white")
#}
plot.idx=plot.idx+1

ylims.plt1=get_ylims(onlymv.pltdata$Abs.Err,onlymv.pltdata$Method,center=0)
abserr.noout.onlymv=ggplot2::ggplot(data=onlymv.pltdata,
                                      ggplot2::aes(y=Abs.Err,x=Epsilon,group=Epsilon))+
  ggplot2::geom_boxplot(outlier.shape=NA)+
  ggplot2::geom_line(
    ggplot2::aes(x=seq(min.x,max.x,length.out=nrow(onlymv.pltdata)),y=rep(0.2,nrow(onlymv.pltdata)),
                 colour="Confidential Effect Error"))+
  ggplot2::labs(x=expression("Privacy Budget"~(epsilon)),
                y="Treatment Effect",color=" ")+
  ggplot2::labs(x=expression("Privacy Budget"~(epsilon)),
                y="Absolute Error",color=" ")+
  ggplot2::theme_minimal(base_size=sim1.base.size)+
  ggplot2::theme(legend.position = "none")+
  ggplot2::scale_y_continuous(limits = c(0,ylims.plt1[2]))+
  ggplot2::ggtitle("MV Histogram")
#  ggplot2::facet_wrap(~nomv.pltdata$Method)
abserr.noout.onlymv
#if(save.plot.png==T){
ggplot2::ggsave(filename=paste0(plot.fname.pref,"sim1_AbsErr_BudgetX_OnlyMVHist_RemoveOutliers.png"),
                width=((plot.sizes[["row1col1"]])[1]),
                height=((plot.sizes[["row1col1"]])[2]),
                units="px",bg="white")
#}
plot.idx=plot.idx+1


abserr.keepout=ggplot2::ggplot(data=plotdata.itt,
                                  ggplot2::aes(y=Abs.Err,x=Epsilon,group=Epsilon))+
  ggplot2::geom_boxplot(outlier.shape=20)+
  ggplot2::geom_line(
    ggplot2::aes(x=seq(min.x,max.x,length.out=nrow(plotdata.itt)),y=rep(0.2,nrow(plotdata.itt)),
                 colour="Confidential Effect Error"))+
  ggplot2::theme_minimal(base_size=sim1.base.size)+
  ggplot2::labs(x=expression("Privacy Budget"~(epsilon)),
                y="Absolute Error",color=" ")+
  ggplot2::theme(legend.position = "bottom")+
  #ggplot2::scale_y_continuous(limits = ylims.plt1)+
  ggplot2::facet_wrap(~plotdata.itt$Method)
abserr.keepout
#if(save.plot.png==T){
ggplot2::ggsave(filename=paste0(plot.fname.pref,"sim1_AbsErr_BudgetX_AllMethods_KeepOutliers.png"),
                width=((plot.sizes[["row1col3"]])[1]),
                height=((plot.sizes[["row1col3"]])[2]),
                units="px",bg="white")
#}
plot.idx=plot.idx+1


# #ylims.plt1=get_ylims(nomv.pltdata$Abs.Err,nomv.pltdata$Method,center=5.20)
# ciover.noout=ggplot2::ggplot(data=plotdata.itt,
#                                   ggplot2::aes(y=CI.overlap,x=Epsilon,group=Epsilon))+
#   ggplot2::geom_boxplot(outlier.shape=NA)+
#   ggplot2::labs(x=expression("Privacy Budget"~(epsilon)),
#                 y="CI Overlap",color=" ")+
#   ggplot2::theme_minimal(base_size=sim1.base.size)+
#   ggplot2::theme(legend.position = "bottom")+
#   ggplot2::scale_y_continuous(limits = c(0,1))+
#   ggplot2::facet_wrap(~nomv.pltdata$Method)
# ciover.noout
# #if(save.plot.png==T){
# ggplot2::ggsave(filename=paste0(plot.fname.pref,"sim1_CIoverlap_BudgetX_AllMethods_RemoveOutliers.png"),
#                 width=((plot.sizes[["row1col3"]])[1]),
#                 height=((plot.sizes[["row1col3"]])[2]),
#                 units="px",bg="white")
# #}
# plot.idx=plot.idx+1



ciover.keepout=ggplot2::ggplot(data=plotdata.itt,
                               ggplot2::aes(y=CI.overlap,x=Epsilon,group=Epsilon))+
  ggplot2::geom_boxplot(outlier.shape=20)+
  ggplot2::labs(x=expression("Privacy Budget"~(epsilon)),
                y="CI Overlap",color=" ")+
  ggplot2::theme_minimal(base_size=sim1.base.size)+
  ggplot2::theme(legend.position = "bottom")+
  #ggplot2::scale_y_continuous(limits = ylims.plt1)+
  ggplot2::facet_wrap(~plotdata.itt$Method)
ciover.keepout
#if(save.plot.png==T){
ggplot2::ggsave(filename=paste0(plot.fname.pref,"sim1_CIoverlap_BudgetX_AllMethods_KeepOutliers.png"),
                width=((plot.sizes[["row1col3"]])[1]),
                height=((plot.sizes[["row1col3"]])[2]),
                units="px",bg="white")
#}


save(sim1tab,sim1plot,plotdata.itt,summary.tab.treff.table,
     sim1.base.size,min.x,max.x,dta.nrow,
     file=paste0(outputpath,"/simulations/sim1_tables_and_figures_",file.suffix,".Rda"))

###############################
############### Sim 2 ##########
 #################################
#load simulation data from models
load(paste0(basepath,"/output/simulations/results_models_",file.suffix,".Rda"))
#data.sim2=(sim2.out[[5]])[[4]]

######## Table #########
sim1.summary.cond1=summary.tab.catcont$Treatment!="control"
sim1.summary.cond2=((summary.tab.catcont$Method!="Confidential")|
  ((summary.tab.catcont$Method=="Confidential")&grepl("Higher",summary.tab.catcont$Model)))
col.order1=c("Model","NCont","NCat","Method","Avg.ITT",
             #"Median.ITT",
             "Avg.CI.overlap",#"Avg.Abs.Err",
             "Median.Abs.Err",
             #"Max.StdErr",
             "ModelLab")
summary.treat.catcont=summary.tab.catcont[sim1.summary.cond1&sim1.summary.cond2,col.order1]
combos.NCov$Name=paste0("HigherBudget_",combos.NCov$Model)
for(md.nm in combos.NCov$Name){
  summary.treat.catcont$NCat[summary.treat.catcont$Model==md.nm]=combos.NCov$NCat[combos.NCov$Name==md.nm]
  summary.treat.catcont$NCont[summary.treat.catcont$Model==md.nm]=combos.NCov$NCont[combos.NCov$Name==md.nm]
  summary.treat.catcont$ModelLab[summary.treat.catcont$Model==md.nm]=combos.NCov$Model[combos.NCov$Name==md.nm]
}
summary.treat.catcont$Algorithm=as.character(summary.treat.catcont$Method)
higher.indic=(grepl("HigherBudget",summary.treat.catcont$Model)&(summary.treat.catcont$Method=="DP Model-based"))
summary.treat.catcont$Algorithm[summary.treat.catcont$Method=="DP Model-based"]="DP-Mb (1)"
summary.treat.catcont$Algorithm[higher.indic]="DP-Mb (2)"
  #paste("Higher Budget",summary.treat.catcont$Algorithm[higher.indic])
#summary.treat.catcont$ModelLab=factor(summary.treat.catcont$ModelLab,levels=mod.order)
table(summary.treat.catcont$Algorithm)
summary.treat.catcont$Algorithm=factor(summary.treat.catcont$Algorithm,levels=c("Confidential","MV Histogram","Hybrid","DP-Mb (1)","DP-Mb (2)"),
                                       labels=c("Confidential","MV Histogram","Hybrid","\\dpmb{1}","\\dpmb{2}"))

summary.tab.catcont.table=summary.treat.catcont%>%arrange(Algorithm)%>%arrange(ModelLab)
col.order=c("ModelLab","NCont","NCat","Algorithm","Avg.ITT",#"Median.ITT",
            "Avg.CI.overlap",#"Avg.Abs.Err",
            "Median.Abs.Err")#,"Max.StdErr")
summary.tab.catcont.table=summary.tab.catcont.table[!is.na(summary.tab.catcont.table$Algorithm),col.order]
rownames(summary.tab.catcont.table)=NULL
sim2table=knitr::kable(summary.tab.catcont.table,
             col.names = c("Model","Cont.","Cat.","Algorithm","Mean $\\treff$ Estimate",#"Median $\\treff$ Estimate",
                           "Mean $\\cioverlap$",
                           #"Mean Abs. Diff. ($\\treff=5$)",
                           "Median Abs. Diff. ($\\treff=5$)"),
                           #"Max Std. Err."),
             digits=table.digits,format="latex",booktabs=T,
             linesep=c("","","","","\\addlinespace"),escape=F)
sim2table
writeLines(sim2table,paste0(basepath,"/output/tables/sim2_models_",file.suffix,".tex"))

treat.itt=ITT2%>%filter(Treatment!="control")
treat.itt$MethMod=paste0(treat.itt$Method,"_",treat.itt$Model)
treat.itt.conf=treat.itt%>%filter(grepl("Confidential",Method))
treat.itt.notdpmb=treat.itt%>%filter(grepl("Hybrid|Histogram",Method))
treat.itt.onlydpmb=treat.itt%>%filter(!grepl("Confidential|Hybrid|Histogram",Method))

temp.conf=dplyr::bind_rows(tapply(treat.itt.conf$StdErr,treat.itt.conf$MethMod,summary))
min(temp.conf$Min.)
max(temp.conf$Max.)
temp.notdpmb=dplyr::bind_rows(tapply(treat.itt.notdpmb$StdErr,treat.itt.notdpmb$MethMod,summary))
min(temp.notdpmb$Min.)
max(temp.notdpmb$Max.)
temp.onlydpmb=dplyr::bind_rows(tapply(treat.itt.onlydpmb$StdErr,treat.itt.onlydpmb$MethMod,summary))
min(temp.onlydpmb$Min.)
max(temp.onlydpmb$Max.)
for(method in unique(treat.itt$Method)){
  temp=treat.itt%>%filter(Method==method)
  temp.se=dplyr::bind_rows(tapply(temp$StdErr,temp$MethMod,summary))
  print(paste(method, "has min", round(min(temp.se$Min.),4),"and max",round(max(temp.se$Max.),4)))
}

tapply(treat.itt$StdErr,treat.itt$MethMod,min)
tapply(treat.itt$StdErr,treat.itt$MethMod,max)

tapply(treat.itt.notdpmb$StdErr,treat.itt.notdpmb$MethMod,summary)
tapply(treat.itt.onlydpmb$StdErr,treat.itt.onlydpmb$MethMod,summary)


############ PLOT ##################
sim2.select.cond1=((ITT2$Method!="Confidential")|
                ((ITT2$Method=="Confidential")&(!grepl("Higher",ITT2$Model))))
ITT2.dedup=ITT2[sim2.select.cond1,]
ITT2.dedup$MethodLab=ITT2.dedup$Method
ITT2.dedup$MethodLab[grepl("Higher",ITT2.dedup$Model)]="DP-Mb (2)"
ITT2.dedup$MethodLab[(!grepl("Higher",ITT2.dedup$Model))&(ITT2.dedup$Method=="DP Model-based")]="DP-Mb (1)"
ITT2.dedup$MethodLab=factor(ITT2.dedup$MethodLab,
                            levels=c("Confidential","MV Histogram","Hybrid",
                                     "DP-Mb (1)","DP-Mb (2)"),
                            labels=c("Confidential","MV Histogram","Hybrid",unname(DPMbplot.name)))
ITT2.dedup$ModelLab=gsub("HigherBudget_","",ITT2.dedup$Model)
ITT2.dedup$ModelLab=factor(ITT2.dedup$ModelLab,levels=c(paste("Model",seq(1,9)),"Unbounded","Finite Set"))

sim2.plotdata=ITT2.dedup%>%
  filter(Treatment!="control")#%>%filter(Method!="Confidential")
not.conf=(sim2.plotdata$Method!="Confidential")
sim2.plotdata$Method=as.character(sim2.plotdata$Method)
for(md.nm in combos.NCov$Name){
  sim2.plotdata$NCat[sim2.plotdata$Model==md.nm]=combos.NCov$NCat[combos.NCov$Name==md.nm]
  sim2.plotdata$NCont[sim2.plotdata$Model==md.nm]=combos.NCov$NCont[combos.NCov$Name==md.nm]
}

sim2.plotlist=NULL
######## Plot 1
cioverlap.by.ncovXmethod.plot=
  ggplot2::ggplot(data=sim2.plotdata[not.conf,],
                  ggplot2::aes(y=CI.overlap,
                               group=ModelLab,
                               x=ModelLab,
                               color=ModelLab))+
  ggplot2::geom_boxplot()+
  ggplot2::labs(x="Model",y="CI Overlap")+
  ggplot2::theme_minimal(base_size=sim2.base.size)+
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=35,vjust=1.05,hjust=0.95),
                 legend.position = "none")+
  grafify::scale_color_grafify("kelly")+
  ggplot2::facet_wrap(~MethodLab)
cioverlap.by.ncovXmethod.plot
#if(save.plot.png==T){
  ggplot2::ggsave(filename=paste0(plot.fname.pref,"sim2_CIoverlap_ModelX_AllMethods_KeepOutliers.png"),
                  width=((plot.sizes[["row2col2"]])[1]),
                  height=((plot.sizes[["row2col2"]])[2]),
                  units="px",bg="white")
#}

sim2.plot.list=c(sim2.plotlist,list("CI Overlap By Method And Model"=cioverlap.by.ncovXmethod.plot))

##### Plot 2 (and variations) ###

only.conf=sim2.plotdata[sim2.plotdata$Method=="Confidential",]
conf.rw=nrow(only.conf)


ylims.v=get_ylims(sim2.plotdata$ITT,sim2.plotdata$MethodLab)


sim2.itt.by.method.model=function(data, only.conf.data, ylims=NULL,
                                  bs.sz=sim2.base.size,outlier.shp=NA,freescale=F){
conf.rw=nrow(only.conf.data)
conf.df.full=NULL
for(MeLab in unique(data$MethodLab)){
  add.method.conf=only.conf.data
  add.method.conf$MethodLab=MeLab
  conf.df.full=c(conf.df.full,list(add.method.conf))
}
conf.df.full=dplyr::bind_rows(conf.df.full)

treff.plot=
  ggplot2::ggplot(data=data,
                  ggplot2::aes(y=ITT,
                               group=ModelLab,
                               x=ModelLab,
                               color=ModelLab))+
  ggplot2::geom_boxplot(outlier.shape = outlier.shp)+
  grafify::scale_color_grafify("kelly")+
  ggplot2::theme_minimal(base_size=sim2.base.size)+
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=35,vjust=1.05,hjust=0.95),
                 legend.position = "none")#margin=ggplot2::margin(c(0,0,0,0))),
                 #axis.title.x = ggplot2::element_text(margin=ggplot2::margin(c(0,0,0,0))),
                 #legend.title=ggplot2::element_blank(),legend.position = "none",
                 #legend.box.margin = ggplot2::margin(c(0,0,0,0)),
                 #legend.margin = ggplot2::margin(c(0,0,0,0)),
                 #legend.spacing.x = ggplot2::unit(0,"mm"),
                 #legend.spacing.y = ggplot2::unit(0,"mm"))#,
                 #plot.margin=ggplot2::margin(c(0,0,0,0)),
                 #panel.spacing.y = ggplot2::unit(0,"mm"))
if(is.null(ylims)==FALSE){
  treff.plot=treff.plot+ggplot2::scale_y_continuous(limits = ylims)
}

  treff.plot=treff.plot+
    ggplot2::geom_point(ggplot2::aes(x=conf.df.full$ModelLab,
                                   y=conf.df.full$ITT,
                                   shape="Confidential Estimator",
                                   color="Confidential"),
                      color="black")+
  ggplot2::labs(x="Model",y="Treatment Effect",shape=" ")+
  ggplot2::scale_shape_manual(values=c(17))+
  ggplot2::guides(color="none",
                  shape=ggplot2::guide_legend(position="bottom"))
  if(freescale==T){
    treff.plot=treff.plot+ggplot2::facet_wrap(~MethodLab,scales = "free")
  }else{
    treff.plot=treff.plot+ggplot2::facet_wrap(~MethodLab)
}
return(treff.plot)
}
treff.bymethodmod.noout=sim2.itt.by.method.model(sim2.plotdata[not.conf,],only.conf.data=only.conf,ylims=ylims.v)
treff.bymethodmod.noout
ggplot2::ggsave(filename=paste0(plot.fname.pref,"sim2_treff_ModelX_AllMethods_RemoveOutliers.png"),
                width=((plot.sizes[["row2col2"]])[1]),
                height=((plot.sizes[["row2col2"]])[2]),
                units="px",bg="white")

treff.bymethodmod.free=sim2.itt.by.method.model(sim2.plotdata[not.conf&(sim2.plotdata$MethodLab!="GenM-1"),],only.conf.data=only.conf,ylims=c(3.7,6.4),outlier.shp = 20)
treff.bymethodmod.free
ggplot2::ggsave(filename=paste0(plot.fname.pref,"sim2_treff_ModelX_ThreeMethods_KeepOutliers.png"),
                width=((plot.sizes[["row1col3"]])[1]),
                height=((plot.sizes[["row1col3"]])[2]),
                units="px",bg="white")
treff.bymethodmod.dpmb1=sim2.itt.by.method.model(sim2.plotdata[not.conf&(sim2.plotdata$MethodLab=="GenM-1"),],only.conf.data=only.conf,outlier.shp = 20)
treff.bymethodmod.dpmb1
ggplot2::ggsave(filename=paste0(plot.fname.pref,"sim2_treff_ModelX_GenM1_KeepOutliers.png"),
                width=((plot.sizes[["row1col1"]])[1]),
                height=((plot.sizes[["row1col1"]])[2]),
                units="px",bg="white")
treff.bymod.noDPmb.wout=sim2.itt.by.method.model(sim2.plotdata[not.conf&(sim2.plotdata$Method!="DP Model-based"),],
                                                 only.conf,outlier.shp=20)
treff.bymod.noDPmb.wout

#if(save.plot.png==T){
  ggplot2::ggsave(filename=paste0(plot.fname.pref,"sim2_treff_ModelsX_NotDPMb_KeepOutliers.png"),
                  width=((plot.sizes[["row1col2"]])[1]),
                  height=((plot.sizes[["row1col2"]])[2]),
                  units="px",bg="white")
#}

treff.bymod.DPmb.noout=sim2.itt.by.method.model(sim2.plotdata[not.conf&(sim2.plotdata$Method=="DP Model-based"),],
                                                only.conf,ylims=ylims.v,outlier.shp=NA)
treff.bymod.DPmb.noout
#if(save.plot.png==T){
  ggplot2::ggsave(filename=paste0(plot.fname.pref,"sim2_treff_ModelsX_OnlyDPMb_RemoveOutliers.png"),
                  width=((plot.sizes[["row1col2"]])[1]),
                  height=((plot.sizes[["row1col2"]])[2]),
                  units="px",bg="white")
#}

treff.bymod.DPmb.wout=sim2.itt.by.method.model(sim2.plotdata[not.conf&(sim2.plotdata$Method=="DP Model-based"),],
                                                 only.conf,outlier.shp=20)
treff.bymod.DPmb.wout
#if(save.plot.png==T){
  ggplot2::ggsave(filename=paste0(plot.fname.pref,"sim2_treff_ModelsX_OnlyDPMb_KeepOutliers.png"),
                  width=((plot.sizes[["row1col2"]])[1]),
                  height=((plot.sizes[["row1col2"]])[2]),
                  units="px",bg="white")
#}




sim2.plot.list=c(sim2.plot.list,
                 list("TreatEffect Boxplot All Method Without Outliers"=treff.bymethodmod.noout,
                   "TreatEffect Boxplot Not DP-Mb WITH Outliers"=treff.bymod.noDPmb.wout,
                      "TreatEffect Boxplot Only DP-Mb Without Outliers"=treff.bymod.DPmb.noout,
                      "TreatEffect Boxplot Only DP-Mb WITH Outliers"=treff.bymod.noDPmb.wout))


##### Plot 3

sim2.plotdata.noadd=sim2.plotdata[grepl("Model",sim2.plotdata$Model),]
not.conf.noadd=(sim2.plotdata.noadd$Method!="Confidential")
sim2.plotdata.noadd.noconf=sim2.plotdata.noadd[not.conf.noadd,]
sim2.plotdata.noadd.noconf$NCatlab=factor(paste("Binary",sim2.plotdata.noadd.noconf$NCat),levels=paste("Binary",c(2,5,10)))
sim2.plotdata.noadd.noconf$NContlab=factor(paste("Continuous",sim2.plotdata.noadd.noconf$NCont),levels=paste("Continuous",c(2,5,10)))
sim2.plotdata.noadd.noconf$MethodLab=factor(sim2.plotdata.noadd.noconf$MethodLab,levels=levels(sim2.plotdata.noadd$MethodLab)[2:5])

cioverlap.by.ncatXncontXmethod.plot=
  ggplot2::ggplot(data=sim2.plotdata.noadd.noconf,
                  ggplot2::aes(y=CI.overlap,group=MethodLab,x=MethodLab,color=MethodLab))+
  ggplot2::geom_boxplot()+
  ggplot2::labs(x="Algorithm",y="CI Overlap",color="Algorithm")+
  ggplot2::theme_bw(base_size=sim2.base.size)+
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=35,vjust=1.05,hjust=0.95),
                 legend.position = "none")+
  #ggplot2::scale_x_discrete(labels=c("MV Histogram","Hybrid","DP Mo"))+
  ggplot2::scale_color_brewer(palette="Dark2")+
  ggplot2::facet_grid(NCatlab~NContlab)
cioverlap.by.ncatXncontXmethod.plot
#if(save.plot.png==T){
  ggplot2::ggsave(filename=paste0(plot.fname.pref,"sim2_CIoverlap_MethodsX.png"),
                  width=((plot.sizes[["row3col3"]])[1]),
                  height=((plot.sizes[["row3col3"]])[2]),
                  units="px",bg="white")
#}

#}
  ylims.v.noadd.noconf=get_ylims(sim2.plotdata.noadd.noconf$ITT,sim2.plotdata.noadd.noconf$MethodLab)



treff.by.ncatXncontXmethod.plot.no.outliers=
  ggplot2::ggplot(data=sim2.plotdata.noadd.noconf,
                  ggplot2::aes(y=ITT,group=MethodLab,x=MethodLab,color=MethodLab))+
  ggplot2::geom_boxplot(outlier.shape = NA)+
  ggplot2::labs(x="Algorithm",y="Treatment Effect",color="Algorithm")+
  ggplot2::theme_bw(base_size=sim2.base.size)+
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=35,vjust=1.05,hjust=0.95),
                 legend.position = "none")+
  ggplot2::scale_y_continuous(limits = ylims.v.noadd.noconf)+
  #ggplot2::scale_x_discrete(labels=c("MV Histogram","Hybrid","DP Mo"))+
  ggplot2::scale_color_brewer(palette="Dark2")+
  ggplot2::facet_grid(NCatlab~NContlab)
treff.by.ncatXncontXmethod.plot.no.outliers

#if(save.plot.png==T){
ggplot2::ggsave(filename=paste0(plot.fname.pref,"sim2_treff_MethodsX_RemoveOutliers.png"),
                width=((plot.sizes[["row3col3"]])[1]),
                height=((plot.sizes[["row3col3"]])[2]),
                units="px",bg="white")
#}

treff.by.ncatXncontXmethod.plot.with.outliers=
  ggplot2::ggplot(data=sim2.plotdata.noadd.noconf,
                  ggplot2::aes(y=ITT,group=MethodLab,x=MethodLab,color=MethodLab))+
  ggplot2::geom_boxplot(outlier.shape = 20)+
  ggplot2::labs(x="Algorithm",y="Treatment Effect",color="Algorithm")+
  ggplot2::theme_bw(base_size=sim2.base.size)+
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=35,vjust=1.05,hjust=0.95),
                 legend.position = "none")+
  #ggplot2::scale_y_continuous(limits = ylims.v.noadd.noconf)+
  #ggplot2::scale_x_discrete(labels=c("MV Histogram","Hybrid","DP Mo"))+
  ggplot2::scale_color_brewer(palette="Dark2")+
  ggplot2::facet_grid(NCatlab~NContlab)
treff.by.ncatXncontXmethod.plot.with.outliers

#if(save.plot.png==T){
ggplot2::ggsave(filename=paste0(plot.fname.pref,"sim2_treff_MethodsX_KeepOutliers.png"),
                width=((plot.sizes[["row3col3"]])[1]),
                height=((plot.sizes[["row3col3"]])[2]),
                units="px",bg="white")

absdiff.by.ncatXncontXmethod.plot.no.outliers=
  ggplot2::ggplot(data=sim2.plotdata.noadd.noconf,
                  ggplot2::aes(y=Abs.Err,group=MethodLab,x=MethodLab,color=MethodLab))+
  ggplot2::geom_boxplot(outlier.shape = NA)+
  ggplot2::labs(x="Algorithm",y="Absolute Difference",color="Algorithm")+
  ggplot2::theme_bw(base_size=sim2.base.size)+
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=35,vjust=1.05,hjust=0.95),
                 legend.position = "none")+
  ggplot2::scale_y_continuous(limits = c(0,2.3))+#,#ylims.v.noadd.noconf)+
  #ggplot2::scale_x_discrete(labels=c("MV Histogram","Hybrid","DP Mo"))+
  ggplot2::scale_color_brewer(palette="Dark2")+
  ggplot2::facet_grid(NCatlab~NContlab)
absdiff.by.ncatXncontXmethod.plot.no.outliers

ggplot2::ggsave(filename=paste0(plot.fname.pref,"sim2_absdiff_MethodsX_RemoveOutliers.png"),
                width=((plot.sizes[["row3col3"]])[1]),
                height=((plot.sizes[["row3col3"]])[2]),
                units="px",bg="white")

absdiff.by.ncatXncontXmethod.plot.with.outliers=
  ggplot2::ggplot(data=sim2.plotdata.noadd.noconf,
                  ggplot2::aes(y=Abs.Err,group=MethodLab,x=MethodLab,color=MethodLab))+
  ggplot2::geom_boxplot(outlier.shape = 20)+
  ggplot2::labs(x="Algorithm",y="Absolute Difference",color="Algorithm")+
  ggplot2::theme_bw(base_size=sim2.base.size)+
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=35,vjust=1.05,hjust=0.95),
                 legend.position = "none")+
  #ggplot2::scale_y_continuous(limits = c(0,2.3))+#,#ylims.v.noadd.noconf)+
  #ggplot2::scale_x_discrete(labels=c("MV Histogram","Hybrid","DP Mo"))+
  ggplot2::scale_color_brewer(palette="Dark2")+
  ggplot2::facet_grid(NCatlab~NContlab)
absdiff.by.ncatXncontXmethod.plot.with.outliers

ggplot2::ggsave(filename=paste0(plot.fname.pref,"sim2_absdiff_MethodsX_KeepOutliers.png"),
                width=((plot.sizes[["row3col3"]])[1]),
                height=((plot.sizes[["row3col3"]])[2]),
                units="px",bg="white")
#}

cioverlap.dist.plot=
  ggplot2::ggplot(data=sim2.plotdata.noadd.noconf,
                  ggplot2::aes(x=CI.overlap,color=MethodLab,fill=MethodLab))+
  ggplot2::geom_histogram(bins=15,position="identity",alpha=0.2)+
  ggplot2::labs(x="CI Overlap",color="Algorithm",y="Frequency",fill="Algorithm")+
  ggplot2::theme_bw(base_size=sim2.base.size)+
  #ggplot2::theme(axis.text.x=ggplot2::element_text(angle=35,vjust=1.05,hjust=0.95),
  #               legend.position = "none")+
  ggplot2::theme(legend.position = "bottom")+
  #ggplot2::scale_x_discrete(labels=c("MV Histogram","Hybrid","DP Mo"))+
  ggplot2::scale_color_brewer(palette="Dark2")+
  ggplot2::scale_fill_brewer(palette="Set2")+
  ggplot2::facet_grid(NCatlab~NContlab,scales = "free_y")
cioverlap.dist.plot
ggplot2::ggsave(filename=paste0(plot.fname.pref,"sim2_cioverlap_hist.png"),
                width=((plot.sizes[["row3col3"]])[1]),
                height=((plot.sizes[["row3col3"]])[2]),
                units="px",bg="white")

sim2.plot.list=c(sim2.plot.list,list("CI Overlap by Number of Variables and Methods"=cioverlap.by.ncatXncontXmethod.plot,
                                    "TreatEffect by Number of Variables and Methods Without Outliers"=treff.by.ncatXncontXmethod.plot.no.outliers,
                 "TreatEffect by Number of Variables and Methods WITH Outliers"=treff.by.ncatXncontXmethod.plot.with.outliers,
                 "AbsDiff by Number of Variables and Methods Without Outliers"=absdiff.by.ncatXncontXmethod.plot.no.outliers,
                 "AbsDiff by Number of Variables and Methods WITH Outliers"=treff.by.ncatXncontXmethod.plot.with.outliers))

save(sim2.plot.list,
     sim2.plotdata.noadd.noconf,sim2.plotdata,only.conf,
     not.conf,ylims.v,ylims.v.noadd.noconf,
     sim2table,sim2.base.size,
     file=paste0(outputpath,"/simulations/sim2_tables_and_figures_",file.suffix,".Rda"))


###############################
############### Sim 2 No Noise ##########
#################################
#load simulation data from models
load(paste0(basepath,"/output/simulations/simulations_nonoise",file.suffix,".Rda"))

t1=sim2_nonoise[[1]]
t3=sim2_nonoise[[3]]
summary.tab.catcont.nn=sim2_nonoise[[3]]

######## Table #########
sim1.summary.cond1=summary.tab.catcont.nn$Treatment!="control"
# sim1.summary.cond2=((summary.tab.catcont$Method!="Confidential")|
#                       ((summary.tab.catcont$Method=="Confidential")&grepl("Higher",summary.tab.catcont$Model)))
col.order1=c("Model","NCont","NCat","Method","Avg.ITT","Avg.CI.overlap","Avg.Abs.Err","Max.StdErr","ModelLab")
summary.treat.catcont.nn=summary.tab.catcont.nn[sim1.summary.cond1,col.order1]

for(md.nm in combos.NCov$Name){
  summary.treat.catcont.nn$NCat[summary.treat.catcont.nn$Model==md.nm]=combos.NCov$NCat[combos.NCov$Name==md.nm]
  summary.treat.catcont.nn$NCont[summary.treat.catcont.nn$Model==md.nm]=combos.NCov$NCont[combos.NCov$Name==md.nm]
  summary.treat.catcont.nn$ModelLab[summary.treat.catcont.nn$Model==md.nm]=combos.NCov$Name[combos.NCov$Name==md.nm]
}
summary.treat.catcont.nn$Algorithm=as.character(summary.treat.catcont.nn$Method)
#higher.indic=(grepl("HigherBudget",summary.treat.catcont$Model)&(summary.treat.catcont$Method=="DP Model-based"))
summary.treat.catcont.nn$Algorithm[summary.treat.catcont.nn$Method=="DP Model-based"]="DP-Mb"
#summary.treat.catcont$Algorithm[higher.indic]="DP-Mb (2)"
#paste("Higher Budget",summary.treat.catcont$Algorithm[higher.indic])
#summary.treat.catcont$ModelLab=factor(summary.treat.catcont$ModelLab,levels=mod.order)
table(summary.treat.catcont.nn$Algorithm)
summary.treat.catcont.nn$Algorithm=factor(summary.treat.catcont.nn$Algorithm,levels=c("Confidential","MV Histogram","Hybrid","DP-Mb"),
                                       labels=c("Confidential","MV Histogram","Hybrid","\\longdpmb"))

summary.tab.catcont.table=summary.treat.catcont.nn%>%arrange(Algorithm)%>%arrange(ModelLab)
col.order=c("ModelLab","NCont","NCat","Algorithm","Avg.ITT","Avg.CI.overlap","Avg.Abs.Err","Max.StdErr")
sim2table=knitr::kable(summary.tab.catcont.table[!is.na(summary.tab.catcont.table$Algorithm),col.order],
                       col.names = c("Model","Cont.","Cat.","Algorithm","Mean $\\treff$ Estimate","Mean $\\cioverlap$",
                                     "Mean Abs. Diff. ($\\treff=5$)",
                                     "Max Std. Err."),
                       digits=table.digits,format="latex",booktabs=T,
                       linesep=c("","","","\\addlinespace"),escape=F)
sim2table
writeLines(sim2table,paste0(basepath,"/output/tables/sim2_models_nonoise_",file.suffix,".tex"))

############ PLOT ##################

# sim2.select.cond1=((ITT2$Method!="Confidential")|
#                      ((ITT2$Method=="Confidential")&(!grepl("Higher",ITT2$Model))))
ITT2.dedup=sim2_nonoise[[1]]#ITT2#[sim2.select.cond1,]
ITT2.dedup$MethodLab=ITT2.dedup$Method
#ITT2.dedup$MethodLab[grepl("Higher",ITT2.dedup$Model)]="DP-Mb (2)"
#ITT2.dedup$MethodLab[(!grepl("Higher",ITT2.dedup$Model))&(ITT2.dedup$Method=="DP Model-based")]="DP-Mb"
ITT2.dedup$MethodLab=factor(ITT2.dedup$MethodLab,
                            levels=c("Confidential","MV Histogram","Hybrid",
                                     "DP Model-based"),
                            labels=c("Confidential","MV Histogram","Hybrid",longDPMb))
#ITT2.dedup$ModelLab=gsub("HigherBudget_","",ITT2.dedup$Model)
ITT2.dedup$ModelLab=factor(ITT2.dedup$Model,levels=c(paste("Model",seq(1,9)),"Unbounded","Finite Set"))

sim2.plotdata=ITT2.dedup%>%
  filter(Treatment!="control")#%>%filter(Method!="Confidential")
not.conf=(sim2.plotdata$Method!="Confidential")
sim2.plotdata$Method=as.character(sim2.plotdata$Method)
for(md.nm in combos.NCov$Name){
  sim2.plotdata$NCat[sim2.plotdata$Model==md.nm]=combos.NCov$NCat[combos.NCov$Name==md.nm]
  sim2.plotdata$NCont[sim2.plotdata$Model==md.nm]=combos.NCov$NCont[combos.NCov$Name==md.nm]
}

sim2.plotlist=NULL
######## Plot 1

cioverlap.by.ncovXmethod.plot=
  ggplot2::ggplot(data=sim2.plotdata[not.conf,],
                  ggplot2::aes(y=CI.overlap,
                               group=ModelLab,
                               x=ModelLab,
                               color=ModelLab))+
  ggplot2::geom_boxplot()+
  ggplot2::labs(x="Model",y="CI Overlap")+
  ggplot2::theme_minimal(base_size=sim2.base.size)+
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=35,vjust=1.05,hjust=0.95),
                 legend.position = "none")+
  grafify::scale_color_grafify("kelly")+
  ggplot2::facet_wrap(~MethodLab)
cioverlap.by.ncovXmethod.plot
#if(save.plot.png==T){
ggplot2::ggsave(filename=paste0(plot.fname.pref,"nonoise_sim2_CIoverlap_ModelX_KeepOutliers.png"),
                 width=((plot.sizes[["row1col3"]])[1]),
                 height=((plot.sizes[["row1col3"]])[2]),
                 units="px",bg="white")
#}
plot.idx=plot.idx+1
sim2.plot.list=c(sim2.plotlist,list("NoNoise CI Overlap By Method And Model"=cioverlap.by.ncovXmethod.plot))

##### Plot 2 (and variations) ###

only.conf=sim2.plotdata[sim2.plotdata$Method=="Confidential",]
conf.rw=nrow(only.conf)

ylims.v.nonoise=get_ylims(sim2.plotdata$ITT,sim2.plotdata$MethodLab)
treff.bymod.wout=sim2.itt.by.method.model(sim2.plotdata[not.conf,],only.conf.data=only.conf,freescale = T,outlier.shp=20)
treff.bymod.wout

ggplot2::ggsave(filename=paste0(plot.fname.pref,"nonoise_sim2_treff_ModelX_AllMethods_KeepOutliers.png"),
                width=((plot.sizes[["row1col3"]])[1]),
                height=((plot.sizes[["row1col3"]])[2]),
                units="px",bg="white")




sim2.plot.list=c(sim2.plot.list,
                 list("NoNoise TreatEffect Boxplot WITH Outliers"=treff.bymod.wout))


##### Plot 3

sim2.plotdata.noadd=sim2.plotdata[grepl("Model",sim2.plotdata$Model),]
not.conf.noadd=(sim2.plotdata.noadd$Method!="Confidential")
sim2.plotdata.noadd.noconf=sim2.plotdata.noadd[not.conf.noadd,]
sim2.plotdata.noadd.noconf$NCatlab=factor(paste("Binary",sim2.plotdata.noadd.noconf$NCat),levels=paste("Binary",c(4,10,20)))
sim2.plotdata.noadd.noconf$NContlab=factor(paste("Continuous",sim2.plotdata.noadd.noconf$NCont),levels=paste("Continuous",c(2,5,10)))
sim2.plotdata.noadd.noconf$MethodLab=factor(sim2.plotdata.noadd.noconf$MethodLab,levels=levels(sim2.plotdata.noadd$MethodLab)[2:5])

cioverlap.by.ncatXncontXmethod.plot=
  ggplot2::ggplot(data=sim2.plotdata.noadd.noconf,
                  ggplot2::aes(y=CI.overlap,group=MethodLab,x=MethodLab,color=MethodLab))+
  ggplot2::geom_boxplot()+
  ggplot2::labs(x="Algorithm",y="CI Overlap",color="Algorithm")+
  ggplot2::theme_bw(base_size=sim2.base.size)+
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=35,vjust=1.05,hjust=0.95),
                 legend.position = "none")+
  #ggplot2::scale_x_discrete(labels=c("MV Histogram","Hybrid","DP Mo"))+
  ggplot2::scale_color_brewer(palette="Dark2")+
  ggplot2::facet_grid(NCatlab~NContlab)
cioverlap.by.ncatXncontXmethod.plot
#if(save.plot.png==T){
ggplot2::ggsave(filename=paste0(plot.fname.pref,"nonoise_sim2_CIOverlap_MethodsX_GridMod.png"),
                width=((plot.sizes[["row3col3"]])[1]),
                height=((plot.sizes[["row3col3"]])[2]),
                units="px",bg="white")
#}

treff.by.ncatXncontXmethod.plot=
  ggplot2::ggplot(data=sim2.plotdata.noadd.noconf,
                  ggplot2::aes(y=ITT,group=MethodLab,x=MethodLab,color=MethodLab))+
  ggplot2::geom_boxplot(outlier.shape = NA)+
  ggplot2::labs(x="Algorithm",y="Treatment Effect",color="Algorithm")+
  ggplot2::theme_bw(base_size=sim2.base.size)+
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=35,vjust=1.05,hjust=0.95),
                 legend.position = "none")+
  #ggplot2::scale_y_continuous(limits = ylim.v.nonoise)+
  #ggplot2::scale_x_discrete(labels=c("MV Histogram","Hybrid","DP Mo"))+
  ggplot2::scale_color_brewer(palette="Dark2")+
  ggplot2::facet_grid(NCatlab~NContlab)
treff.by.ncatXncontXmethod.plot

#if(save.plot.png==T){
ggplot2::ggsave(filename=paste0(plot.fname.pref,"nonoise_sim2_treff_MethodsX_GridMod.png"),
                width=((plot.sizes[["row3col3"]])[1]),
                height=((plot.sizes[["row3col3"]])[2]),
                units="px",bg="white")
#}
plot.idx=plot.idx+1
sim2.plot.list=c(sim2.plot.list,list("NoNoise CI Overlap by Number of Variables and Methods"=cioverlap.by.ncatXncontXmethod.plot,
                                     "NoNoise TreatEffect by Number of Variables and Methods Without Outliers"=treff.by.ncatXncontXmethod.plot))

save(sim2.plot.list,
     sim2.plotdata.noadd.noconf,sim2.plotdata,only.conf,
     not.conf,min.y,max.y,
     sim2table,sim2.base.size,
     file=paste0(outputpath,"/simulations/sim2_tables_and_figures_nonoise",file.suffix,".Rda"))


#### On this version, something went wrong with the diagonostic plot save, so I re-generated them
source(paste0(basepath,"/program/functions/reg_assumptions_function.R"))
load(paste0(basepath,"/output/simulations/checkpoint_",sim2.conf.suffix,"_",file.suffix,".Rda"))
conf.plot=list("pt.sz1"=4,"ln.sz1"=2.5,"bs.sz1"=33,
               "pt.sz2"=2,"ln.sz2"=1,"bs.sz2"=14,
               "ncol1"=3,"ncol2"=2,
               "width1"=3500,"height1"=2000,
               "width2"=1000,"height2"=300)
mod.names=names(reg.formulas.sim2)
add.plots=lapply(seq(1+length(reg.formulas.sim2)-2,length(reg.formulas.sim2)),
                 function(idx)reg_assumptions(reg.formulas.sim2[idx],
                                              data=sim.dfs,
                                              response.var="y",
                                              pt.sz=conf.plot$pt.sz2,
                                              ln.sz=conf.plot$ln.sz2,
                                              bs.sz=conf.plot$bs.sz2,
                                              mod.name=mod.names[idx]))
save.plot.name.xtra=paste0(basepath,"/output/figures/diagnostic_",sim2.conf.suffix,"_",file.suffix,"_xtra.png")
png(filename=save.plot.name.xtra
    ,width=conf.plot$width2,height=conf.plot$height2)
print(cowplot::plot_grid(plotlist=add.plots,ncol=conf.plot$ncol2))
dev.off()

conf.plots=lapply(seq(1,length(reg.formulas.sim2)-2),
                 function(idx)reg_assumptions(reg.formulas.sim2[idx],
                                              data=sim.dfs,
                                              response.var="y",
                                              pt.sz=conf.plot$pt.sz1,
                                              ln.sz=conf.plot$ln.sz1,
                                              bs.sz=conf.plot$bs.sz1,
                                              mod.name=mod.names[idx]))

save.plot.name=paste0(basepath,"/output/figures/diagnostic_",sim2.conf.suffix,"_",file.suffix,".png")
png(filename=save.plot.name
    ,width=conf.plot$width1,height=conf.plot$height1)
print(cowplot::plot_grid(plotlist=conf.plots,ncol=conf.plot$ncol1))
dev.off()
#
#
head(sim1)
