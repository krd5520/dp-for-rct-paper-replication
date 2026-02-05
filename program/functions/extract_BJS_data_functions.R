
## Given zip_path which is the file path to the zipped Blattman et al. data and an output directory to save the relevant data
# extract just the STYL_Final.dta file and place it in the specified output directory.
extract_styl_final<-function(zip_path,
                             output_dir=paste0(rprojroot::find_rstudio_root_file(),"/data/liberia_replication_data"),
                             output_subset_data_filename="LiberiaRound5.Rda",
                             quietly=F){
  main.start=proc.time()
  if(!file.exists(zip_path)){ #check zip file exists
    stop('Zip file not found')
  }

  dir.create(output_dir,showWarnings = F,recursive=T) #create output directory

  in_zip=unzip(zip_path,list=TRUE)$Name
  styl=grep("STYL_Final\\.dta$",in_zip,value=T)

  if(length(styl)==0){
    stop(paste0("\nSTYL_Final.dta not found inside ",zip_path,".\n Contents: ",paste(in_zip,collapse=", ")))
  }

  unzip(zip_path, files=styl[1],exdir=output_dir)

  final_path<-file.path(output_dir,"STYL_Final.dta")

  extracted_path<-file.path(output_dir,styl[1])

  #if final_path not equal to extracted path, then there are directories that STYL_Final.dta is in within output_dir
  if(final_path!=extracted_path){
    file.rename(extracted_path,final_path) #move STYL_Final.dta to output_dir level

    #get a list of the sub-directories of output_dir that STYL_Final.dta was originally in
    extracted_dirs_list=list.dirs(dirname(extracted_path),full.names = T, recursive = T)
    final_dirs_list=list.dirs(dirname(final_path),full.names = T, recursive = T)
    to_remove_dirs=setdiff(final_dirs_list,extracted_dirs_list)#R.utils::getRelativePath(extracted_path,output_dir)

    for(remove_dir in rev(to_remove_dirs)){ #remove each of the empty subdirectories
      if(remove_dir!=output_dir){
        unlink(remove_dir,recursive=T)
      }
    }

  }

  liberia.sub=which_variables(final_path,file.path(output_dir,output_subset_data_filename),quietly=quietly)
  file.remove(final_path)
  if(quietly==F){
    timeextract=proc.time()-main.start
    cat(paste0("\nExtract and subset BJS data time: ",round(timeextract[[3]],2)," seconds."),"")
  }
  return(liberia.sub) #don't return anything
}


which_variables=function(stata_data_path,subset_data_path, quietly=F){
  #Replicating Table 2b in Original Paper using AER.do file as reference
  # The AER.do and the data STYL_Final.dta file can be found in
  # basepath/data/liberia_replication_data/AER.do and
  # basepath/data/liberia_replication_data/STYL_Final.dta
  # respectively or from the original source:
  ## Blattman, Christopher, Jamison, Julian C., and Sheridan, Margaret.
  ### Replication data for: Reducing Crime and Violence: Experimental Evidence from
  ### Cognitive Behavioral Therapy in Liberia. Nashville, TN: American Economic
  ### Association [publisher], 2017. Ann Arbor, MI: Inter-university Consortium for
  ### Political and Social Research [distributor], 2019-10-12. https://doi.org/10.3886/E113056V1

  #'base' variable in AER_1.do
  baseline.full=c('age_b', 'livepartner_b', 'mpartners_b', 'hhunder15_b', 'famseeoften_b',
                  'muslim_b', 'school_b', 'schoolbasin_b', 'literacy_b', 'mathscore_b', 'health_resc_b',
                  'disabled_b', 'depression_b', 'distress_b', 'rel_commanders_b', 'faction_b', 'warexper_b',
                  'profitsump99avg7d_b', 'wealth_indexstd_b', 'homeless_b', 'slphungry7dx_b', 'savstockp99_b',
                  'loan50_b', 'loan300_b', 'illicit7da_zero_b', 'agricul7da_zero_b', 'nonagwage7da_zero_b',
                  'allbiz7da_zero_b', 'nonaghigh7da_zero_b', 'agriculeveramt_b', 'nonagbizeveramt_b',
                  'nonaghigheveramt_b', 'drugssellever_b', 'drinkboozeself_b', 'druggrassself_b',
                  'grassdailyuser_b', 'harddrugsever_b', 'harddrugsdailyuser_b', 'steals_b',
                  'stealnb_nonviol_b', 'stealnb_felony_b', 'disputes_all_b', 'asbhostil_b',
                  'conscientious_b', 'neurotic_b', 'grit_b', 'rewardresp_b', 'locuscontr_b', 'impulsive_b',
                  'selfesteem_b', 'patient_game_real_b', 'inconsistent_game_resc_b', 'risk_game_resc_b',
                  'timedecl_b', 'riskdecl_b', 'cognitive_score_b', 'ef_score_b')


  #'dvs_t2' in AER_1.do
  response.stem=c('fam_asb','drugssellever','stealnb','disputes_all_z','carryweapon',
                  'arrested','asbhostilstd','domabuse_z')

  #'strata' in AER_1.do
  block.vars=c("tp_strata_alt","cg_strata") #block variables

  #focusing on round==5 (term=lt in AER.do file)
  #read in the liberia data
  liberia.full=haven::read_stata(stata_data_path)
  #Data from:
  ##Blattman, Christopher, Jamison, Julian C., and Sheridan, Margaret. Replication
  ##data for: Reducing Crime and Violence: Experimental Evidence from Cognitive Behavioral
  ##Therapy in Liberia. Nashville, TN: American Economic Association [publisher], 2017.
  ##Ann Arbor, MI: Inter-university Consortium for Political and Social Research [distributor],
  ##2019-10-12. https://doi.org/10.3886/E113056V1
  liberia.full=labelled::remove_attributes(liberia.full,"format.stata")


  y.sub.data=liberia.full[,grepl(paste0(response.stem,collapse="|"),colnames(liberia.full))]
  # Instructions.pdf say "_st" and "_b" are short term and pretreatment values
  if(quietly==F){
    varnames=names(sapply(y.sub.data[,grepl("_lt",colnames(y.sub.data))],function(x)attr(x,"label"),simplify=TRUE))
    varlabels=unname(sapply(y.sub.data[,grepl("_lt",colnames(y.sub.data))],function(x)attr(x,"label"),simplify=TRUE))
    cat("\nBased on the variable names and labels (seen below) Table2b in Blattman et al. uses the suffix '_lt' on 'fam_asb' response and '_ltav' on the other response variables.")
    print(data.frame("Variable Names"=varnames,"Labels"=varlabels))
  }
  #fam_asb_lt is the correct choice. The other repsonse variables have "_ltav"

  #list all possible response variables for round==5
  response.vars.allltav = c('fam_asb_lt',paste0(response.stem[-1],"_ltav"))
  treat.vars=c("cashassonly","tpassonly","tpcashass") #treatment variables

  col.vars=c(response.vars.allltav,
             treat.vars,block.vars,baseline.full,c("control"))

  #Section about ITT tables in AER_1.do, we see the regression model for long term restricts the data to unfound_wave==0 and round==5
  row.conditions=((liberia.full$round==5)&(liberia.full$unfound_wave==0))
  #subset data
  liberia.sub=liberia.full[row.conditions,col.vars[!duplicated(col.vars)]]

  if(quietly==FALSE){
    cat(paste0("\nRestricting data to round 5 and unfound_wave 0 rows (",nrow(liberia.sub),"/",nrow(liberia.full),").",
                 "\nKeeping only ",length(colnames(liberia.sub)),"/",length(colnames(liberia.full)),") columns which were used in Table 2b of Blattman et al. (2017).",
                 "\nAdd a column 'treatment' which combined the three treatment variable indicators."))
  }

    #combine treatment variables into 1 column
    # (this is useful for the dp_synthdata and hybrid_dp function)
    liberia.sub$treatment=ifelse(liberia.sub$cashassonly==1,"cashassonly","tpassonly")
    liberia.sub$treatment[liberia.sub$tpcashass==1]="tpcashass"
    liberia.sub$treatment[liberia.sub$control==1]="control"

    if(quietly==FALSE){
                cat(paste0("\nSaving subsetted data to ",subset_data_path))
  }

  save(liberia.sub,file=subset_data_path)#save smaller data
  return(liberia.sub)

}

#x=extract_styl_final("C:/Users/Kaitlyn/Downloads/113056-V2.zip")
#which_variables("C:/Users/Kaitlyn/Documents/GitHub/dp-for-rct-paper/data/liberia_replication_data/STYL_Final.dta")

