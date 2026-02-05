##############################################
####### USER DEFINED GLOBAL VARIABLES ########
##############################################

#Note: BJS (2017) refers to Blattman et al. (2017) and BJS (2025) refers to Blattman et al. (2025)

####### File locations ####################
code_path = rprojroot::find_rstudio_root_file() #path to main dp-for-rct-paper-replication folder
zipped_BJS_download="C:/Users/Kaitlyn/Downloads/113056-V2.zip" #path to download BJS (2025) zipped materials

#location to save subset of BJS (2025) data which only includes rows and columns for Table 2b of BJS (2017).
BJS_data_path=paste0(code_path,"/data/LiberiaRound5.Rda")

##############################################
############# Simulation Variables ###########
##############################################

simulations.file.suffix="02022026"
simulations.generated.data.directory=paste0(code_path,"/data/simulated_data")
simulations.output.directory=paste0(code_path,"/output")

##### Simulation Plot Parameters #####

# In Simulation 2, diagnostic plots are saved (QQ Norm Plot, and Residual vs. Fit)
# pt.sz is size of points, ln.sz is size of line, bs.sz is base size for text
# width1 and height1 are for the size of the png saved for the standard models
# additional models are saved using width2 and height2
simulation.plot.parameters=list("pt.sz1"=0.3,"ln.sz1"=0.7,"bs.sz1"=8,
               "pt.sz2"=0.4,"ln.sz2"=0.9,"bs.sz2"=9,
               "ncol1"=3,"ncol2"=2,
               "width1"=8.5,"height1"=5,
               "width2"=5,"height2"=2)


#####################################################
######## Liberia Study Replication Variables #########
#####################################################

