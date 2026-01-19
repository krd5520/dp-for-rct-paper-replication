devtools::install_github("krd5520/DPrct")
library(dplyr) #to use %>% command
####### Testing the multivariate histogram method on the Liberia Data
basepath = rprojroot::find_rstudio_root_file()
main.start=proc.time()

outputpath=paste0(basepath,"/output")
file.suffix="v1214_2"
table.digits=2
sim1.base.size=13
sim2.base.size=13
sim1.plotsize=list("width"=5000,"height"=5000)
sim2.plotsize=list("width"=5000,"height"=5000)

###############################
############### Sim 2 ##########
#################################
################### Tables ###########################
#load simulation data for privacy budgets
load(paste0(outputpath,"/simulations/results_budget_",file.suffix,".Rda"))


colnames(summary.tab.treff)
col.order=c("Epsilon","Delta","Algorithm","Avg.ITT","Avg.CI.overlap","Avg.Abs.Err","Max.StdErr")
conf.row=grepl("Confidential",summary.tab.treff$Method)
summary.tab.treff.table=rbind(summary.tab.treff[conf.row,col.order],summary.tab.treff[!conf.row,col.order])
#copy and paste the output of the following into latex to get table
# OR use R markdown with latex capabilities and chunk output="asis"
sim1tab=knitr::kable(summary.tab.treff.table,
             col.names = c("\\epsilon","\\delta","Algorithm",
                           "Mean $\\treff$ Estimate","Mean $\\cioverlap$", "Mean Abs. Diff.",
                           "Max Std. Err."),
             digits=table.digits,format="latex",booktabs=T,
             linesep=c("\\addlinespace",rep(c("","","\\addlinespace"),7)),
             escape=F)
sim1tab
writeLines(sim1tab,file=paste0(outputpath,"/tables/sim1_budget_",file.suffix,".tex"))


################### Plots ###########################
plotdata.itt=ITT1%>%filter(Treatment!="control")%>%
  filter(Method!="Confidential")%>%
  filter(Delta="0")
dta.nrow=nrow(plotdata.itt)
min.x=floor(min(plotdata.itt$Epsilon)-0.5)
max.x=ceiling(max(plotdata.ittt$Epsilon)+0.5)
sim1plot=ggplot2::ggplot(data=plotdata.itt,
                ggplot2::aes(y=ITT,x=Epsilon,group=Epsilon))+
  ggplot2::geom_boxplot()+
  ggplot2::geom_line(
    ggplot2::aes(x=seq(min.x,max.x,length.out=dta.nrow),y=rep(5,dta.nrow),
                 colour="True Treatment Effect"))+
  ggplot2::geom_line(
    ggplot2::aes(x=seq(min.x,max.x,length.out=dta.nrow),
                 y=rep(unlist(ITT1$ITT[(ITT1$Treatment!="control")&(ITT1$Method=="Confidential")]),dta.nrow/20),
                 colour="Confidential Estimated Treatment Effect"),linetype="dashed")+
  ggplot2::labs(x=expression("Privacy Budget"~(epsilon)),
                y="Treatment Effect",color=" ")+
  ggplot2::theme_minimal(base_size=sim1.base.size)+
  ggplot2::theme(legend.position = "bottom")+
  ggplot2::facet_wrap(~as.factor(plotdata.itt$Method))
sim1plot
png(filename=paste0(outputpath,"/figures/sim1_treateffectDist_",file.suffix,".png"),
    width=sim1.plotsize$width,height=sim1.plotsize$height)
print(sim1plot)
dev.off()


save(sim1tab,sim1plot,plotdata.itt,summary.tab.treff.table,
     filename=paste0(outputpath,"/simulations/sim1_tables_and_figures_",file.suffix,".Rda"))

###############################
############### Sim 2 ##########
 #################################
#load simulation data from models
load(paste0(basepath,"/output/simulations/results_models_",file.suffix,".Rda"))


######## Table #########
col.order=c("NCont","NCat","Algorithm","Avg.ITT","Avg.CI.overlap","Avg.Abs.Err","Max.StdErr")
conf.row=grepl("Confidential",summary.tab.catcont$Method)
summary.tab.catcont.table=rbind(summary.tab.catcont[conf.row,col.order],summary.tab.catcont[!conf.row,col.order])

sim2table=knitr::kable(summary.tab.catcont.table,
             col.names = c("Cont.","Cat.","Algorithm","Mean $\\treff$ Estimate","Mean $\\cioverlap$",
                           "Mean Abs. Diff. ($\\treff=5$)",
                           "Max Std. Err."),
             digits=table.digits,format="latex",booktabs=T,
             linesep=c("","","","\\addlinespace"),escape=F)
sim2table
writeLines(sim2tab,file=paste0(basepath,"/output/tables/sim2_models_",file.suffix,".tex"))

############ PLOT ##################
plotdata.cioverlap=ITT2%>%
  filter(Treatment!="control")%>%
  filter(Method!="Confidential")
plotdata.cioverlap$Method=as.character(plotdata.cioverlap$Method)

plotdata.cioverlap.noadd=plotdata.cioverlap%>%filter(grepl("Model",Model))
cioverlap.by.ncovXmethod.plot=
  ggplot2::ggplot(data=plotdata.cioverlap,
                ggplot2::aes(y=CI.overlap,
                             group=Model,
                             x=Model,
                             color=as.factor(Model)))+
  ggplot2::geom_boxplot()+
  ggplot2::labs(x="Model",y="CI Overlap")+
  ggplot2::theme_minimal(base_size=sim2.base.size)+
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=35,vjust=1.05,hjust=0.95),
                 legend.position = "none")+
  ggplot2::scale_color_brewer(palette="YlGnBu")+
 ggplot2::facet_wrap(~as.factor(plotdata.cioverlap$Method))
cioverlap.by.ncovXmethod.plot


cioverlap.by.ncatXncontXmethod.plot=
  ggplot2::ggplot(data=plotdata.cioverlap.noadd,
                ggplot2::aes(y=CI.overlap,group=Method,x=as.factor(desc(Method)),color=Method))+
  ggplot2::geom_boxplot()+
  ggplot2::ggtitle("CI Overlap for Sanitizing Algorithms across Models")+
  ggplot2::labs(x="Algorithm",y="CI Overlap",color="Algorithm")+
  ggplot2::theme_minimal(base_size=sim2.base.size)+
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=35,vjust=1.05,hjust=0.95),
                 legend.position = "none")+
  ggplot2::scale_x_discrete(labels=c("MV Histogram","Hybrid","DP Model-based"))+
  ggplot2::scale_color_brewer(palette="Dark2")+
  ggplot2::facet_grid(plotdata.cioverlap$Ncatlab~plotdata.cioverlap$Ncontlab)
cioverlap.by.ncatXncontXmethod.plot

png(paste0(outputpath,"/figures/sim2_CIoverlap_by_ncatcont_",file.suffix,".png"),width=sim2.plotsize$width,height=sim2.plotsize$height)
print(cioverlap.by.ncatXncontXmethod.plot)
dev.off()

plotdata.itt2=ITT2%>%
  filter(Treatment!="control")%>%
  filter(Method!="Confidential")
plotdata.itts2$Method=as.character(plotdata.itt2$Method)
sim2.conf.df=plotdata.itt2=ITT2%>%
  filter(Treatment!="control")%>%
  filter(Method=="Confidential")%>%select("Model","ITT")
colnames(sim2.conf.df)=c("Model","Conf.ITT")
plotdata.itt2=left_join(plotdata.itt2,sim2.conf.df,by=Model)

plt.nrow=nrow(plotdata.itt2)
itt.by.modelXmethod.plot=
  ggplot2::ggplot(data=plotdata.itt2,
                  ggplot2::aes(y=ITT,group=Method,x=desc(Method)))+
  ggplot2::geom_boxplot()+
  ggplot2::ggtitle("Estimated Treatment Effect for Sanitizing Algorithms across Models")+
  ggplot2::labs(x="Algorithm",y="Estimated Treatment Effect",color=" ")+
  ggplot2::geom_line(
    ggplot2::aes(x=desc(Method),y=rep(5,plt.nrow),
               colour="True Treatment Effect"))+
  ggplot2::geom_line(
    ggplot2::aes(x=desc(Method),y=Conf.ITT,
                 colour="Confidential Estimated Treatment Effect"),linetype="dashed")+
  ggplot2::theme_minimal(base_size=sim2.base.size)+
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=35,vjust=1.05,hjust=0.95),
                 legend.position = "none")+
  ggplot2::scale_x_discrete(labels=c("MV Histogram","Hybrid","DP Model-based"))+
  ggplot2::scale_color_brewer(palette="Dark2")+
  ggplot2::facet_grid(~plotdata.itt2$Model)
itt.by.modelXmethod.plot


png(filename=paste0(basepath,"/output/figures/sim2_treateffectDist_",file.suffix,".png"),
    width=sim2.plotsize$width,height=sim2.plotsize$height)
print(itt.by.modelXmethod.plot)
dev.off()

save(sim2tab,summary.tab.catcont.table,
     cioverlap.by.ncovXmethod.plot,cioverlap.by.ncatXncontXmethod.plot,itt.by.modelXmethod.plot,
     plotdata.itt2,plotdata.cioverlap,
     filename=paste0(outputpath,"/simulations/sim2_tables_and_figures_",file.suffix,".Rda"))

