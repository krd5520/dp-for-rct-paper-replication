#devtools::install_github("krd5520/DPrct")

####### Processing  Liberia Data
basepath = rprojroot::find_rstudio_root_file()
main.start=proc.time()


#file.suffix="v0321_all_eps2_binssixth_nostandardize"
file.suffix="v0321_all_eps2_binssixth_allstandardize"
outputpath=paste0(basepath,"/output")
in.datapath=paste0(outputpath,"/liberia/liberia_ITT_compare_methods_data_",file.suffix,".Rda")
out.datapath=paste0(basepath,"/data")

load(in.datapath)
#pretable_process parameters
#change the names of DP model-based method
# change.DPMb.names.tables=list("DP-Mb Same Privacy"="\\dpmb{3}",
#                        "DP-Mb_Eps11_Del0"="\\dpmb{7}",
#                        "DP-Mb_Eps19_Del0"="\\dpmb{11}")
# change.method.names.plots=list("Confidential"="Confidential",
#                                "MV Histogram"="MV Histogram",
#                                "Hybrid"="Hybrid",
#   "\\dpmb{3}"="GenM-3",
#                              "\\dpmb{7}"="GenM-7", #(3,.375)
#                              "\\dpmb{11}"="GenM-11")

#pretable_process parameters
# #change the names of DP model-based method
# change.DPMb.names.tables=list("DP-Mb Same Privacy"="\\dpmb{2}",
#                               "DP-Mb_Eps10_Del0"="\\dpmb{6}",
#                               "DP-Mb_Eps18_Del0"="\\dpmb{10}")
# change.method.names.plots=list("Confidential"="Confidential",
#                                "MV Histogram"="MV Histogram",
#                                "Hybrid"="Hybrid",
#                                "\\dpmb{2}"="GenM-2",
#                                "\\dpmb{6}"="GenM-6", #(3,.375)
#                                "\\dpmb{10}"="GenM-10")

#change the names of DP model-based method
change.DPMb.names.tables=list("DP-Mb Same Privacy"="\\dpmb{2}",
                              "DP-Mb_Eps10_Del0"="\\dpmb{10}",
                              "DP-Mb_Eps18_Del0"="\\dpmb{18}")
change.method.names.plots=list("Confidential"="Confidential",
                               "MV Histogram"="MV Histogram",
                               "Hybrid"="Hybrid",
                               "\\dpmb{2}"="GenM-2",
                               "\\dpmb{10}"="GenM-10", #(3,.375)
                               "\\dpmb{18}"="GenM-18")
adj.out=TRUE #adjust p-values

#liberia_table2b_covset_bymethod paramater
include.control.mean=F #include control mean in itt tables like 2b?
#column widths for latex tables in
itt.conf.tex.cols=c(0.23,rep(c(0.07,0.07,0.07),3)) #confidential data itt 2b table
itt.alg.tex.cols=c(0.21,rep(c(0.072,0.055,0.055,0.068),3)) #each synthetic data itt 2b table
itt.famasb.tex.cols=c(0.17,rep(c(0.061,0.055,0.055,0.055),3))

#liberia_table2b_oneresponse inputs
itt.pval.tex.cols=c(0.1,0.25,0.2,0.2,0.2)
alpha=0.05
insignificant.color="gray"
insignificant.color.name="gray"

#heat map figure inputs
cioverlap.cols=RColorBrewer::brewer.pal(9,"Blues")
absdiff.cols=RColorBrewer::brewer.pal(9,"Purples")
standard.theme=ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 27, vjust = 1.15, hjust=1),
                              title=ggplot2::element_text(size=6),
                              axis.title = ggplot2::element_blank())
standard.facet=ggplot2::facet_wrap(~Treatment)
combo.whs=c(6,4.5,1) #combination plot width height scale for png saved
individual.whs=c(6,2,1) #individual plot width height scale for png saved
heat.map.bs.sz=8

color.treff=RColorBrewer::brewer.pal(3,"Dark2")

diag.plt.whs=c(4.25,1.75,1) #individual plot width, height, scale for png save on diagnostic plots

#variables etc.
response.vars.all = unique(full.comparison$ModelName)
main.text.response.vars=c('fam_asb_lt','stealnb_ltav','carryweapon_ltav','asbhostilstd_ltav')

treat.vars=c("cashassonly","tpassonly","tpcashass","Control")
treat.names=c("Cash Only","Therapy Only","Both","Control Mean")



source(paste0(basepath,"/program/functions/table_and_figure_functions_alt.R"))
source(paste0(basepath,"/program/functions/liberia_tables_comptime_2b.R"))
source(paste0(basepath,"/program/functions/directory_setup.R"))

dir.set=directory_setup(basepath,outputpath,
                        list("tables",paste0("tables/liberia_",file.suffix),
                             paste0("tables/liberia_",file.suffix,"/main"),
                             "figures",paste0("figures/liberia_",file.suffix),
                             paste0("figures/liberia_",file.suffix,"/main"),
                             paste0("figures/liberia_",file.suffix,"/diagnostic_main")))

tables=list(liberia_comptimes(by.response.times=by.response.times,summary.times=summary.times,
                           response.vars=response.vars.all,
                           by.response.fname="liberia_comptime_byresponse",
                           summary.fname="liberia_summary_time",
                           outputpath=paste0(outputpath,"/tables/liberia_",file.suffix),
                           secs.100ths=FALSE,
                           change.DPMb.names=change.DPMb.names.tables))

comparison.df=pretable_process(full.comparison=full.comparison,adj.out=T,
                        change.DPMb.names=change.DPMb.names.tables,
                        change.mod.names=NULL)
#comparison.df.remove.obs=comparison.df

if(length(unique(full.comparison$Sim))>1){
tables=c(tables,list(liberia_table2b_covsets_bymethod(full.comparison=comparison.df,n.response=length(unique(comparison.df$ModelName)),
                        tab.digits=2,include.control.mean=include.control.mean,
                        outputpath=paste0(outputpath,"/tables/liberia_",file.suffix),
                        conf.2b.fname="liberia_allcovsets_confidential",
                        san.tex.prefix="liberia_allcovsets_sanitized",
                        out.fname="liberia_allcovsets_allmethods",
                        conf.tex.cols=itt.conf.tex.cols,alg.tex.cols = itt.alg.tex.cols,packrows=T)))
}
print("here")
comparison.df$dataPlot=factor(as.character(comparison.df$data),levels=names(change.method.names.plots))
levels(comparison.df$dataPlot)=unname(sapply(change.method.names.plots,"[[",1))

for(yname in main.text.response.vars){
tables=c(tables,list(liberia_table2b_oneresponse(full.comparison=comparison.df,
                            response.name=yname,
                                 tab.digits=2,include.pval=F,
                            tex.cols=itt.famasb.tex.cols,
                                 outputpath=paste0(outputpath,"/tables/liberia_",file.suffix,"/main"),
                                 tex.prefix="liberia_itt",
                                 out.fname="liberia_itt",includeconf=T,packrows=T)))
temp=one_response_treateffects(data=comparison.df,response=yname,
                               treat.vars=treat.vars[1:3],treat.names=treat.names[1:3],
                               scale.colors=color.treff,
                               line.sz=1,
                               plt.title=NULL,base.size=8,
                               outputpath=paste0(outputpath,"/figures/liberia_",file.suffix,"/main/"),plot.file.prefix="liberia_",
                               plot.save.size=individual.whs,cov.set="Reduced Subset",only.high.genm=T)
temp=one_response_treateffects(data=comparison.df,response=yname,
                               treat.vars=treat.vars[1:3],treat.names=treat.names[1:3],
                               scale.colors=color.treff,
                               line.sz=1,
                               plt.title=NULL,base.size=8,
                               outputpath=paste0(outputpath,"/figures/liberia_",file.suffix,"/main/"),plot.file.prefix="liberia_",
                               plot.save.size=individual.whs,cov.set="Subset",only.high.genm=T)

}

tables=c(tables,list(itt_pval_table(data=comparison.df,alpha=alpha,
                                    not.signif.color = insignificant.color,
                                    not.signif.color.name=insignificant.color.name,
                   tex.cols=itt.pval.tex.cols,
                   outputpath=paste0(outputpath,"/tables/liberia_",file.suffix),
                   tex.prefix="liberia_itt_pval",
                   out.fname="liberia_itt_pval")))




ci.plot.list=heatmap_x_covset(comparison.df,
                            value.col="ci.overlap",x.col="dataPlot",y.col="ModelName",
                            reverse.y=T,reverse.x=F,
                            scale.colors=cioverlap.cols,cut.bounds=c(0,1),
                            value.name="CI Overlap",x.name="Method",y.name="Response",base.size=heat.map.bs.sz,
                            extra.theme=standard.theme,extra.facet=standard.facet,
                            outputpath=paste0(outputpath,"/figures/liberia_",file.suffix),
                            plot.file.prefix="cioverlap",
                            combo.whs=combo.whs,
                            individual.whs=individual.whs,
                            treat.vars=treat.vars,
                            treat.names = treat.names)


adiff.plot.list=heatmap_x_covset(comparison.df,
                                 value.col="abs.diff",x.col="dataPlot",y.col="ModelName",
                                 reverse.y=T,reverse.x=F,
                                 scale.colors=absdiff.cols,
                                 value.name="Absolute Difference",x.name="Method",y.name="Response",base.size=heat.map.bs.sz,
                                 extra.theme=standard.theme,
                                 extra.facet=standard.facet,
                                 outputpath=paste0(outputpath,"/figures/liberia_",file.suffix),
                                 plot.file.prefix="absdiff",
                                 log.transform.value=TRUE,
                                 combo.whs=combo.whs,
                                 individual.whs=individual.whs,
                                 treat.vars=treat.vars,
                                 treat.names = treat.names)


print("Before Main Text Variable Stuff")
### MAIN TEXT Variables
main.text.df=comparison.df[comparison.df$ModelName%in%main.text.response.vars,]
main.text.df$ModelName=factor(as.character(main.text.df$ModelName),levels=main.text.response.vars)
rownames(main.text.df)=NULL



tables=c(tables,list(itt_pval_table(data=main.text.df,alpha=alpha,
                                    not.signif.color = insignificant.color,
                                    not.signif.color.name=insignificant.color.name,
                                    tex.cols=itt.pval.tex.cols,
                                    outputpath=paste0(outputpath,"/tables/liberia_",file.suffix,"/main"),
                                    tex.prefix="liberia_itt_pval",
                                    out.fname="liberia_itt_pval")))


standard.theme=ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 27, vjust = 1.15, hjust=1),
                              title=ggplot2::element_text(size=6),
                              axis.title = ggplot2::element_blank())
standard.facet=ggplot2::facet_wrap(~Treatment)

ci.plot.list=heatmap_x_covset(main.text.df,
                              value.col="ci.overlap",x.col="ModelName",y.col="dataPlot",
                              reverse.y=T,reverse.x=F,
                              scale.colors=cioverlap.cols,cut.bounds=c(0,1),
                              value.name="CI Overlap",x.name="Response",y.name="Method",base.size=heat.map.bs.sz,
                              extra.theme=standard.theme,extra.facet=standard.facet,
                              outputpath=paste0(outputpath,"/figures/liberia_",file.suffix,"/main"),
                              plot.file.prefix="cioverlap",
                              combo.whs=combo.whs,
                              individual.whs=individual.whs,
                              treat.vars=treat.vars,
                              treat.names = treat.names)

adiff.plot.list=heatmap_x_covset(main.text.df,
                                 value.col="abs.diff",x.col="ModelName",y.col="dataPlot",
                                 reverse.y=T,reverse.x=F,
                                 scale.colors=absdiff.cols,
                                 value.name="Absolute Difference",x.name="Response",y.name="Method",base.size=heat.map.bs.sz,
                                 extra.theme=standard.theme,extra.facet=standard.facet,
                                 outputpath=paste0(outputpath,"/figures/liberia_",file.suffix,"/main"),
                                 plot.file.prefix="absdiff",
                                 log.transform.value=TRUE,
                                 combo.whs=combo.whs,
                                 individual.whs=individual.whs,
                                 treat.vars=treat.vars,
                                 treat.names = treat.names)


diag.response.vars=main.text.response.vars[2:4]
main.idx=which(response.vars.all %in% diag.response.vars)
method.names=c("Confidential","MV Histogram","Hybrid",names(change.DPMb.names.tables)[1])
method.fname=c(tolower(unlist(gsub("[^A-z0-9]","",method.names))))
print(method.fname)

method.names=unname(change.method.names.plots[1:4])
for(m in seq(1,length(method.fname))){
  load(paste0(outputpath,"/figures/liberia_",file.suffix,"/liberia_diagnostics_reduced_",method.fname[m],".Rda"))
  for(i in seq(1,length(main.idx))){
    temp.title=cowplot::ggdraw()+cowplot::draw_label(method.names[m])
    temp.plot=cowplot::plot_grid(temp.title,san.plots[[main.idx[i]]],ncol=1,rel_heights = c(0.1,1))
    individual.fname=paste0(outputpath,"/figures/liberia_",file.suffix,"/diagnostic_main/liberia_diag_select_red_",method.fname[m],"_",main.text.response.vars[i+1],".png")
    ggplot2::ggsave(filename=individual.fname,plot=temp.plot,bg="white",scale=diag.plt.whs[3],width=diag.plt.whs[1],height=diag.plt.whs[2])
  }
}

method.names=names(change.DPMb.names.tables)
method.fname=c(tolower(unlist(gsub("[^A-z0-9]","",method.names))))
method.names=unname(change.method.names.plots[seq(4,3+length(method.fname))])
for(m in seq(1,length(method.fname))){
  load(paste0(outputpath,"/figures/liberia_",file.suffix,"/liberia_diagnostics_reduced_",method.fname[m],".Rda"))
    temp.title=cowplot::ggdraw()+cowplot::draw_label(method.names[m],hjust=0)
    temp.plot=cowplot::plot_grid(cowplot::plot_grid(temp.title,NULL,NULL,ncol=3,rel_widths = c(1,1,1)),san.plots[[main.idx[3]]],ncol=1,rel_heights = c(0.1,1))
    individual.fname=paste0(outputpath,"/figures/liberia_",file.suffix,"/diagnostic_main/liberia_diag_select_red_",method.fname[m],"_",main.text.response.vars[4],".png")
    ggplot2::ggsave(filename=individual.fname,plot=temp.plot,bg="white",scale=diag.plt.whs[3],width=diag.plt.whs[1],height=diag.plt.whs[2])
  }

