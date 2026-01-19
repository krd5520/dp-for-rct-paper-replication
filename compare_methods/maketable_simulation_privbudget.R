devtools::install_github("krd5520/DPrct")
library(dplyr) #to use %>% command
####### Testing the multivariate histogram method on the Liberia Data
basepath = rprojroot::find_rstudio_root_file()
main.start=proc.time()

################### Tables ###########################
#load simulation data for privacy budgets
load(paste0(basepath,"/output/compare_privbudget_simulations_v1212.Rda"))

#get table of Average Treatment Effect, Average Overlap, Average Absolute Diff, Max Std. Err
ITT1=dplyr::bind_rows(lapply(sim1.out,"[[",1)) #combine all repetitions
ITT1$groupid=paste0(ITT1$Model,"_",ITT1$Method,"_",ITT1$Treatment)

ITT1$groupid=paste0(ITT1$Model,"_",ITT1$Method,"_",ITT1$Treatment)
summary.tab=ITT1%>%
  #filter(Model!="MassiveEps")%>%
  group_by(Model,Method,Treatment,groupid)%>%
  summarise("Avg.ITT"=mean(ITT),
            "Avg.CI.overlap"=mean(CI.overlap,na.rm=T),
            "Avg.Abs.Err"=mean(Abs.Err,na.rm=T),
            #"Avg.Rel.Err"=mean(Rel.Err,na.rm=T),
            "Max.StdErr"=max(StdErr) )
summary.tab$Epsilon=1
summary.tab$Epsilon[grepl("EpsHalf",summary.tab$groupid)]=0.5
summary.tab$Epsilon[grepl("Eps2",summary.tab$groupid)]=5 #because on 12/12 I misnamed this privacy budget
summary.tab$Epsilon[grepl("Eps2half",summary.tab$groupid)]=2.5
summary.tab$Epsilon[grepl("MassiveEps",summary.tab$groupid)]=5000
summary.tab$Epsilon[grepl("Confidential",summary.tab$groupid)]=0
summary.tab$Epsilon[grepl("Eps5",summary.tab$groupid)]=5
summary.tab$Delta=ifelse(grepl("Approx",summary.tab$groupid),"1/1000","0")#c(rep("$\\frac{1}{1000}$",12)," "," ",rep("0",24))
summary.tab$Algorithm="DP Model-based"
summary.tab$Algorithm[grepl("Confidential",summary.tab$groupid)]="Confidential"
summary.tab$Algorithm[grepl("MV",summary.tab$groupid)]="MV Histogram"
summary.tab$Algorithm[grepl("Hybrid",summary.tab$groupid)]="Hybrid"

summary.tab.treff=summary.tab[summary.tab$Treatment!="control",]%>%
  #filter(Epsilon!=2000)%>%
  arrange(desc(Algorithm))%>%arrange(desc(Delta))%>%arrange(desc(Epsilon))
colnames(summary.tab.treff)
#copy and paste the output of the following into latex to get table
# OR use R markdown with latex capabilities and chunk output="asis"
knitr::kable(summary.tab.treff[c(nrow(summary.tab.treff),seq(1,nrow(summary.tab.treff)-1)),c(9,10,11,5,6,8)],
             col.names = c("\\epsilon","\\delta","Algorithm",
                           "Mean $\\treff$ Estimate","Mean $\\cioverlap$", #"Mean Abs. Diff.",
                           "Max Std. Err."),
             digits=2,format="latex",booktabs=T,
             linesep=c("\\addlinespace",rep(c("","","\\addlinespace"),7)),
             escape=F)

summary.tab.control=summary.tab[grepl("control",summary.tab$groupid),]%>%
  arrange(desc(Algorithm))%>%arrange(desc(Delta))%>%arrange(desc(Epsilon))
#copy and paste the output of the following into latex to get table
# OR use R markdown with latex capabilities and chunk output="asis"
knitr::kable(summary.tab.control[c(19,seq(1,18)),c(6,7,8,2,3,4,5)],
             col.names = c("Epsilon","Delta","Algorithm","Mean $\\treff$ Estimate","Mean $\\cioverlap$",
                           "Mean Abs. Diff. ($\\treff=5$)",
                           "Max Std. Err."),
             digits=3,format="latex",booktabs=T,linesep=c("","","","\\addlinespace"),escape=F)



# #load simulation data from models
# load(paste0(basepath,"/output/compare_covariate_number_simulations_v1205.Rda"))
# ITT2=dplyr::bind_rows(lapply(sim2.out,"[[",1))
# ITT2$groupid=paste0(ITT2$Model,"_",ITT2$Method,"_",ITT2$Treatment)
# #ITT2$NCont=as.numeric(sapply(strsplit(ITT2$Model,"_"),function(x)gsub("Cont","",x[1])))
# #ITT2$NCat=as.numeric(sapply(strsplit(ITT2$Model,"_"),function(x)gsub("Bin","",x[2])))
# ITT2$NCovariates=ITT2$NCat+ITT2$NCont
# summary.tab=ITT2%>%group_by(Method,Treatment,Model)%>%
#   summarise( "Avg.ITT"=mean(ITT),
#     "Avg.CI.overlap"=mean(CI.overlap,na.rm=T),
#             "Avg.Abs.Err"=mean(Abs.Err,na.rm=T),
#            # "Avg.Rel.Err"=mean(Rel.Err,na.rm=T),
#             "Max.StdErr"=max(StdErr)
#            )
# #reorder Method factor levels
# summary.tab$Method=factor(summary.tab$Method,levels=c("Confidential","MV Histogram","Hybrid","DP Model-based"))
# summary.tab$ModelLab=factor(summary.tab$Model,levels=c(paste("Model",seq(1,9)),"FiniteSet","Unbounded"))
# summary.tab.models=summary.tab%>%#filter(!grepl("Finite|Bounded",Model))%>%
#   arrange(Method)%>%arrange(ModelLab)
# knitr::kable(summary.tab.models[summary.tab.models$Treatment!="control",c(8,1,4,5,7)],
#              col.names = c("Model","Algorithm","Mean $\\treff$ Estimate","Mean $\\cioverlap$",
#                                                                     #"Mean Abs. Diff. ($\\treff=5$)",
#                                                                     "Max Std. Err."),
#              digits=2,format="latex",booktabs=T,linesep=c("","","","\\addlinespace"),escape=F)
#

#load simulation data from models
load(paste0(basepath,"/output/compare_covariate_number_simulations_v1129.Rda"))
ITT2=dplyr::bind_rows(lapply(sim2.out,"[[",1))
ITT2$groupid=paste0(ITT2$Model,"_",ITT2$Method,"_",ITT2$Treatment)
ITT2$NCont=as.numeric(sapply(strsplit(ITT2$Model,"_"),function(x)gsub("Cont","",x[1])))
ITT2$NCat=as.numeric(sapply(strsplit(ITT2$Model,"_"),function(x)gsub("Bin","",x[2])))
ITT2$NCovariates=ITT2$NCat+ITT2$NCont
summary.tab=ITT2%>%group_by(Method,Treatment,Model,NCont,NCat,NCovariates)%>%
  summarise( "Avg.ITT"=mean(ITT),
             "Avg.CI.overlap"=mean(CI.overlap,na.rm=T),
             "Avg.Abs.Err"=mean(Abs.Err,na.rm=T),
             # "Avg.Rel.Err"=mean(Rel.Err,na.rm=T),
             "Max.StdErr"=max(StdErr)
  )
#reorder Method factor levels
summary.tab$Method=factor(summary.tab$Method,levels=c("Confidential","MV Histogram","Hybrid","DP Model-based"))
#summary.tab$ModelLab=factor(summary.tab$Model,levels=c(paste("Model",seq(1,9)),"FiniteSet","Unbounded"))
summary.tab.catcont=summary.tab%>%#filter(!grepl("Finite|Bounded",Model))%>%
  arrange(Method)%>%arrange(NCat)%>%arrange(NCont)
rownames(summary.tab.catcont)=NULL
knitr::kable(summary.tab.catcont[summary.tab.catcont$Treatment!="control",c(4,5,1,7,8,9,10)],
             col.names = c("Cont.","Cat.","Algorithm","Mean $\\treff$ Estimate","Mean $\\cioverlap$",
                           "Mean Abs. Diff. ($\\treff=5$)",
                           "Max Std. Err."),
             digits=2,format="latex",booktabs=T,linesep=c("","","","\\addlinespace"),escape=F)

############ PLOTS ##################
plotdata.cioverlap=ITT2%>%filter(Treatment!="control")%>%filter(Method!="Confidential")
plotdata.cioverlap$Method=as.character(plotdata.cioverlap$Method)
# ggplot2::ggplot(data=plotdata.cioverlap,
#                 ggplot2::aes(x=NCovariates,y=CI.overlap,color=Method))+
#   ggplot2::geom_boxplot()+
#   #ggplot2::geom_point()+
#   #ggplot2::geom_smooth(se=F)+
#   ggplot2::scale_color_brewer(palette="Dark2")#+
#  # ggplot2::facet_wrap(~ITT2$Method,scales="free")
ggplot2::ggplot(data=plotdata.cioverlap,
                ggplot2::aes(y=CI.overlap,group=NCovariates,x=NCovariates,color=as.factor(NCovariates)))+
  ggplot2::geom_boxplot()+
  #ggplot2::geom_point()+
  #ggplot2::geom_smooth(se=F)+
  ggplot2::scale_color_brewer(palette="YlGnBu")+
 ggplot2::facet_wrap(~as.factor(plotdata.cioverlap$Method))


plotdata.cioverlap$Ncatlab=factor(paste(plotdata.cioverlap$NCat,"binary covariates"),levels=paste(c(4,10,20),"binary covariates"))
#levels(plotdata.cioverlap$Ncatlab)=c("4 binary covariates","10 binary covariates", "20 binary covariates")
plotdata.cioverlap$Ncontlab=factor(paste(plotdata.cioverlap$NCont,"continuous covariates"),levels=paste(c(2,5,10),"continuous covariates"))
#levels(plotdata.cioverlap$Ncontlab)=c("2 continuous covariates","5 continuous covariates","10 continuous covariates")
#levels(plotdata.cioverlap$Method)=c("MV Histogram","Hybrid","DP Model-based")

ncat.by.ncont.by.method.plot=
  ggplot2::ggplot(data=plotdata.cioverlap,
                ggplot2::aes(y=CI.overlap,group=Method,x=as.factor(desc(Method)),color=Method))+
  ggplot2::geom_boxplot()+
  #ggplot2::geom_point()+
  #ggplot2::geom_smooth(se=F)+
  ggplot2::ggtitle("CI Overlap for Sanitizing Algorithms across Models")+
  ggplot2::labs(x="Algorithm",y="CI Overlap",color="Algorithm")+
  ggplot2::theme_bw()+
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=35,vjust=1.05,hjust=0.95),
                 legend.position = "none")+
  ggplot2::scale_x_discrete(labels=c("MV Histogram","Hybrid","DP Model-based"))+
  ggplot2::scale_color_brewer(palette="Dark2")+
  ggplot2::facet_grid(plotdata.cioverlap$Ncatlab~plotdata.cioverlap$Ncontlab)
ncat.by.ncont.by.method.plot

png(paste0(basepath,"/output/figures/ModelsCIoverlapPlot.png"),width=5000,height=5000)
ncat.by.ncont.by.method.plot
dev.off()

plotdata.itt=ITT1%>%filter(Treatment!="control")%>%filter(Method!="Confidential")%>%
  filter(grepl("Pure",Model))
plotdata.itt$Epsilon=1
plotdata.itt$Epsilon[grepl("EpsHalf",plotdata.itt$Model)]=0.5
plotdata.itt$Epsilon[grepl("Eps2",plotdata.itt$Model)]=5 #because on 12/12 I misnamed this privacy budget
plotdata.itt$Epsilon[grepl("Eps2half",plotdata.itt$Model)]=2.5
plotdata.itt$Epsilon[grepl("Eps5",plotdata.itt$Model)]=5
#plotdata.itt$Method=as.character(plotdata.itt$Method)
# ggplot2::ggplot(data=plotdata.cioverlap,
#                 ggplot2::aes(x=NCovariates,y=CI.overlap,color=Method))+
#   ggplot2::geom_boxplot()+
#   #ggplot2::geom_point()+
#   #ggplot2::geom_smooth(se=F)+
#   ggplot2::scale_color_brewer(palette="Dark2")#+
#  # ggplot2::facet_wrap(~ITT2$Method,scales="free")
ggplot2::ggplot(data=plotdata.itt,
                ggplot2::aes(y=ITT,x=Epsilon,group=Epsilon))+#,x=NCovariates,color=as.factor(NCovariates)))+
  ggplot2::geom_boxplot()+
  #ggplot2::geom_abline(intercept=5,slope=0,color="black")+
  ggplot2::geom_line(ggplot2::aes(x=seq(0,6,length.out=240),y=rep(5,240),colour="True Treatment Effect"))+
  ggplot2::geom_line(ggplot2::aes(x=seq(0,6,length.out=240),y=rep(unlist(ITT1$ITT[(ITT1$Treatment!="control")&(ITT1$Method=="Confidential")]),12),colour="Confidential Estimated Treatment Effect"),linetype="dashed")+
  ggplot2::labs(x="Privacy Budget",y="Treatment Effect",color=" ")+
    ggplot2::theme_minimal(base_size=13)+
  ggplot2::theme(legend.position = "bottom")+
 # ggplot2::scale_colour_brewer(palette="Dark2")+
  #ggplot2::annotate("text",x=2.5,y=4.8,label="True Treatment Effect",color="black")+
  #ggplot2::annotate("text",x=2,y=5.5,label="Confidential Estimated Treatment Effect",color="grey37")
  ggplot2::facet_wrap(~as.factor(plotdata.itt$Method))


plotdata.cioverlap$Ncatlab=factor(paste(plotdata.cioverlap$NCat,"binary covariates"),levels=paste(c(4,10,20),"binary covariates"))
#levels(plotdata.cioverlap$Ncatlab)=c("4 binary covariates","10 binary covariates", "20 binary covariates")
plotdata.cioverlap$Ncontlab=factor(paste(plotdata.cioverlap$NCont,"continuous covariates"),levels=paste(c(2,5,10),"continuous covariates"))
#levels(plotdata.cioverlap$Ncontlab)=c("2 continuous covariates","5 continuous covariates","10 continuous covariates")
#levels(plotdata.cioverlap$Method)=c("MV Histogram","Hybrid","DP Model-based")

ncat.by.ncont.by.method.plot=
  ggplot2::ggplot(data=plotdata.cioverlap,
                  ggplot2::aes(y=ITT,group=Method,x=as.factor(desc(Method)),color=Method))+
  ggplot2::geom_boxplot()+
  #ggplot2::geom_point()+
  #ggplot2::geom_smooth(se=F)+
  ggplot2::ggtitle("Estimated Treatment Effect for Sanitizing Algorithms across Models")+
  ggplot2::labs(x="Algorithm",y="Estimated Treatment Effect",color="Algorithm")+
  ggplot2::theme_bw()+
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=35,vjust=1.05,hjust=0.95),
                 legend.position = "none")+
  ggplot2::scale_x_discrete(labels=c("MV Histogram","Hybrid","DP Model-based"))+
  ggplot2::scale_color_brewer(palette="Dark2")+
  ggplot2::facet_grid(plotdata.cioverlap$Ncatlab~plotdata.cioverlap$Ncontlab)
ncat.by.ncont.by.method.plot

png(paste0(basepath,"/output/figures/ModelsCIoverlapPlot.png"),width=5000,height=5000)
ncat.by.ncont.by.method.plot
dev.off()

time.df=dplyr::bind_rows(lapply(sim2.out,function(x)as.data.frame(attr(x[[1]],"all.comp.times")$Per.Model)))
time.df$measured=rownames(time.df)
rownames(time.df)=NULL
time.df$Sim=as.factor(rep(seq(1,20),each=14))
time.df$Method=rep(c("Confidential",rep("MV Histogram",4),rep("Hybrid",4),rep("DP Model-based",5)),20)
time.df=time.df[((time.df$Method=="MV Histogram")&(grepl("synthdata|utility",time.df$measured)))|(time.df$Method!="MV Histogram"),]
time.df.summ=time.df[time.df$Method=="Confidential",c(seq(1,9),11,12)]
for(algo in c("MV Histogram","Hybrid","DP Model-based")){
  print(algo)
  sub.df=time.df[time.df$Method==algo,c(1:10)]
  print(sub.df)
  sub.df=time.df%>%group_by(Sim)%>%summarise_if(is.numeric,sum)
  sub.df$Method=algo

  time.df.summ=dplyr::bind_rows(time.df.summ,sub.df)
}
T.time.df=tidyr::pivot_longer(time.df.summ,cols=colnames(time.df)[1:9],names_to="Model",values_to = "Time")
T.time.df$NCont=rep(c(2,5,10),240)
T.time.df$NCat=rep(rep(c(4,10,20),each=3),80)
T.time.df$NCovariates=T.time.df$NCat+T.time.df$NCont

ggplot2::ggplot(data=T.time.df%>%filter(Method!="Confidential"),ggplot2::aes(x=NCovariates,y=Time,color=Method,group=Method))+
  ggplot2::geom_point()+
  ggplot2::geom_smooth(se=F)+
  ggplot2::facet_wrap(~T.time.df$Method[T.time.df$Method!="Confidential"])
