####### Processing  Liberia Data
basepath = rprojroot::find_rstudio_root_file()
main.start=proc.time()
outputpath=paste0(basepath,"/output")
source(paste0(basepath,"/program/functions/directory_setup.R"))

dir.set=directory_setup(basepath,outputpath,list("tables"))

########## Global for Model Variables ############
load(paste0(basepath,"/data/liberia_round5_labels.Rda"))

#A reduced subset of the original set
red.sub.covariates=c('age_b', 'asbhostil_b', 'drugssellever_b', 'drinkboozeself_b',
                     'druggrassself_b', 'harddrugsever_b', 'steals_b')
#A subset of the original set of covariates
sub.covariates=c(red.sub.covariates, 'wealth_indexstd_b','school_b', 'cognitive_score_b')

response.all = c('fam_asb_lt',
                      paste0(c('drugssellever','stealnb','disputes_all_z','carryweapon',
                               'arrested','asbhostilstd','domabuse_z'),"_ltav"))
treat.vars=c("cashassonly","tpassonly","tpcashass") #treatment variables
block.vars=c("tp_strata_alt","cg_strata") #block variables

n.stdev=5

adjust_label_pre=function(labstr){
  labstr=gsub("[^!]*:[a-z]{2}[0-9]*[^a-zA-Z0-9]","",labstr)
  labstr=gsub("[^!]*: [a-z]{2}[0-9]*[^a-zA-Z0-9]","",labstr)
  labstr=gsub("[A-Za-z]*: [0-9]*[^a-zA-Z0-9]","",labstr)
  if(length(labstr)>1){
    out=sapply(labstr,function(x)stringr::str_to_sentence(stringr::str_to_lower(base::trimws(x))))
  }else{
    out=stringr::str_to_sentence(stringr::str_to_lower(base::trimws(labstr)))
  }
  out
}

adjust_label_post=function(labstr){
  labstr=gsub("2 weeks (el) or 4 weeks (bl)","2 or 4 weeks",labstr)
  labstr=gsub(" ( [^!]*)","",labstr)
  labstr=gsub("[^!]*:","",labstr)
  labstr=gsub("Indicator for ","",labstr)
  labstr=gsub("Hours per week","Weekly hours",labstr,ignore.case = T)
  labstr=gsub("\\([^!]*","",labstr)
  labstr=gsub("add risk averse","",labstr)
  labstr=gsub("-- [^!]*","",labstr)
  if(length(labstr)>1){
    out=sapply(labstr,function(x)base::trimws(stringr::str_to_sentence(stringr::str_to_lower(x))))
  }else{
    out=base::trimws(stringr::str_to_sentence(stringr::str_to_lower(labstr)))
  }
  out
}

adjust_index_label=function(df.range){
  labstr=df.range$label

  #handle case when no numbers given
  hasindexlab=grepl("Index",labstr,ignore.case = T)
  df.range$min[hasindexlab]=df.range$obs.min[hasindexlab]
  df.range$max[hasindexlab]=df.range$obs.max[hasindexlab]

  #handle case when (#-#) is the range expressed in variable label
  hasindexrange=grepl("[0-9]-[0-9]*",labstr)
  index.range=sapply(labstr[hasindexrange],function(x)regmatches(x,regexpr("[0-9]*-[0-9]*",x)))#gsub(") [^!]*","",labstr[hasindexrange])
  df.range$min[hasindexrange]=as.numeric(sapply(strsplit(index.range,"-"),"[[",1))
  df.range$max[hasindexrange]=as.numeric(sapply(strsplit(index.range,"-"),"[[",2))
  df.range$label[hasindexrange]=gsub(" ([0-9]*-[0-9]*)[^!]*","",df.range$label[hasindexrange])

  #handle case when # to # is the range expressed in variable label
  hasindexrangeto=grepl("[0-9] to [0-9]*",labstr)
  index.rangeto=sapply(labstr[hasindexrangeto],function(x)regmatches(x,regexpr("[0-9]* to [0-9]*",x)))#gsub(") [^!]*","",labstr[hasindexrange])
  df.range$min[hasindexrangeto]=as.numeric(sapply(strsplit(index.rangeto," to "),"[[",1))
  df.range$max[hasindexrangeto]=as.numeric(sapply(strsplit(index.rangeto," to "),"[[",2))
  df.range$label[hasindexrangeto]=gsub(" ([0-9]* to [0-9]*)[^!]*","",df.range$label[hasindexrangeto])

  # #handle case when (#-#) is the range expressed in variable label
  # hasindexrangedash=grepl("[0-9] -[0-9]*",labstr)
  # index.rangedash=sapply(labstr[hasindexrangedash],function(x)regmatches(x,regexpr("[0-9]* -[0-9]*",x)))#gsub(") [^!]*","",labstr[hasindexrange])
  # print(head(index.rangedash))
  # df.range$min[hasindexrangedash]=as.numeric(sapply(strsplit(index.rangedash,"-"),"[[",1))
  # df.range$max[hasindexrangedash]=as.numeric(sapply(strsplit(index.rangedash,"-"),"[[",2))
  # df.range$label[hasindexrangedash]=gsub(" ([0-9]* -[0-9]*)[^!]*","",df.range$label[hasindexrangedash])

  rownames(df.range)=NULL
  return(df.range)

}

adjust_binary_label=function(df.range,liberia.sub){
  indicator.desc=grepl("Indicator",df.range$label,ignore.case = T)
  lowcats=(df.range$nVals<3)
  df.range$label[indicator.desc]=gsub("Indicator for","",df.range$label[indicator.desc],ignore.case=T)
  df.range$type[lowcats|indicator.desc]="Binary"
  df.range$min[lowcats|indicator.desc]=0
  df.range$max[lowcats|indicator.desc]=1
  return(df.range)
}

manual_adjust=function(df.range,nstddev,liberia.sub){
  y0=c("stealnb_ltav","Average count of stealing events in past 2 weeks","Average Discrete","0","70") #5 stealing a day for past 2 weeks
  y1=c("carryweapon_ltav","Average indicator of carrying a weapon","Average Binary","0","1")
  y2=c("arrested_ltav","Average indicator if arrested in past 2 weeks","Average Binary","0","1")
  man0=c("age_b","Age","Discrete","18","35") #from paper (turns out its wrong)
  man1=c("mpartners_b","Count of partners","Discrete","0","5")  #made up this bounds
  man2=c("hhunder15_b","Count of children<15 in household","Discrete","0","25") #made up this bound
  man3=c("famseeoften_b","Sees family often","Categorical","1","4") #technically alsp have 97 and 98 for "don't know" and "refuse to answer" but no observations have these, so removing.
  man4=c("school_b","Years of school","Discrete","0","18") #up to masters' level
  man5=c("rel_commanders_b","Relations to commanders index","Index","0","4")
  man6=c("warexper_b","War experience index","Index","0","12")
  man7=c("slphungry7dx_b","Days slept hungry","Discrete","0","7")
  man8=c("illicit7da_zero_b","Hours in illicit activities","Continuous","0","168") #assuming this is also per week?
  man9=c("agricul7da_zero_b","Hours per week in agricultue","Continuous","0","168")
  man10=c("nonagwage7da_zero_b","Hours per week in low-skill wage labor","Continuous","0","168")
  man11=c("allbiz7da_zero_b","Hours per week in low-skill business","Continuous","0","168")
  man12=c("nonaghigh7da_zero_b","Hours per week in high-skill work","Continuous","0","168")
  man13=c("agriculeveramt_b","Years in agriculture","Discrete","0","35")
  man14=c("nonagbizeveramt_b","Years in non-agriculture business","Discrete","0","35")
  man15=c("nonaghigheveramt_b","Years in high-skill work","Discrete","0","35")
  man16=c("stealnb_nonviol_b","Count of nonviolent stealing events","Discrete","0","70") #5 stealing incidents a day for 2 weeks?
  man17=c("stealnb_felony_b","Count of felony stealing events","Discrete","0","70")

  #In the paper dispute_all_b is listed as an index 0-9, but it takes values 0-30 in the data
  #after some investigation it appears that disputes_all_b is a sum of number of fights and the indicator of if fined.
  #The documentation suggests it was meant to be a sum of indicators for 8 types of fights and the fine indicator.
  # Thus it should be treated as a discrete and lower bounded variable (upper bound is 2.5 fights per day in last 2 weeks?)
  man18=c("disputes_all_b","Count of disputes and fights in past 2 weeks","Discrete","0","35")
  man19=c("riskdecl_b","Declared risk appetite","Standardized",as.character(-nstddev),as.character(nstddev))
  man20=c("timedecl_b","Declared patience","Standardized",as.character(-nstddev),as.character(nstddev))

  missing=colnames(liberia.sub)[!(colnames(liberia.sub)%in%df.range$variable)]
  if(length(missing)>0){
    message(paste("Some variables have no labels from the dataset: ",paste(missing,collapse=", ")))
  }
  #Liberia national income average in 2010 (World Inequality Database) converted to 2010 USD (BLS.gov CPI Inflation Calculator)
  xtra0=c("profitsump99avg7d_b","Weekly cash earnings USD","Top-Coded Continuous","0","1000")
  #Liberia national wealth per person in 2010 (World Inequality Database) converted to 2010 USD  (BLS.gov CPI Inflation Calculator)
  xtra1=c("savstockp99_b","Stock saving USD","Top-Coded Continuous","0","2200")
  # conscientious_b is labeled differently in paper than in data
  xtra2=c("conscientious_b","Conscientious index","Index","0","24")
  # grit_b is labeled differently in paper than in data
  xtra3=c("grit_b","Grit index","Index","0","21")
  # rewardresp_b is labeled differently in paper than in data
  xtra4=c("rewardresp_b","Reward responsiveness index","Index","0","24")
  # impulsive_b is labeled differently in paper than in data
  xtra5=c("impulsive_b","Impulsiveness index","Index","0","21")


  manual.mat=matrix(c(y0,y1,y2,man0,man1,man2,man3,man4,man5,man6,man7,man8,man9,man10,
                      man11,man12,man13,man14,man15,man16,man17,man18,man19,man20,
                      xtra0,xtra1,xtra2,xtra3,xtra4,xtra5),
                    ncol=5,byrow=T)
  colnames(manual.mat)=c("variable","label","type","min","max")
  manual.mat=as.data.frame(manual.mat)
  manual.mat$min=as.numeric(manual.mat$min)
  manual.mat$max=as.numeric(manual.mat$max)
  order.vars=df.range$variable[(df.range$variable %in% manual.mat$variable)]
  rw.reorder=sapply(order.vars,function(x)which(manual.mat$variable==x))
  manual.mat=manual.mat[unlist(rw.reorder),]
  df.range.rm=df.range[(df.range$variable %in% manual.mat$variable),]
  df.range[df.range$variable %in%manual.mat$variable,c("variable","label","type","min","max")]=manual.mat
  #df.range[df.range$variable %in%manual.mat$variable,c("variable","label","type","min","max")]=
  #  lapply(df.range[df.range$variable %in%manual.mat$variable,c("variable","label","type","min","max")],
  #         function(x)manual.mat[manual.mat$variable==x[1],])
  #df.range=dplyr::bind_rows(df.range[!(df.range$variable %in% manual.mat$variable),],manual.mat)
  df.range.rm=df.range.rm[is.na(df.range.rm$min)==FALSE,]
  vars.rm=sapply(df.range.rm$variable,function(x)
    ((df.range.rm$min[df.range.rm$variable==x]!=manual.mat$min[manual.mat$variable==x])|(df.range.rm$max[df.range.rm$variable==x]!=manual.mat$max[manual.mat$variable==x])|(df.range.rm$type[df.range.rm$variable==x]!=manual.mat$type[manual.mat$variable==x])))
  if(sum(vars.rm)>0){
    message(paste("Manual variable differs from automatic generated information:",paste(names(vars.rm)[vars.rm],collapse=", ")))

    df.range.rm=dplyr::left_join(df.range.rm[df.range.rm$variable %in%c(names(vars.rm)[vars.rm]),c("variable","type","min","max","obs.min","obs.max")],
                                 manual.mat[manual.mat$variable%in%c(names(vars.rm)[vars.rm]),c("variable","type","min","max")],
                                 by="variable",suffix=c("auto","manual"))
    if(nrow(df.range.rm)==1){
      for(coli in seq(1,3)){
        if(df.range.rm[1,coli+1]!=df.range.rm[1,coli+6]){
          message(paste(gsub("auto","",colnames(df.range.rm)[coli+1]),"is auto-generated as",df.range.rm[1,coli+1],", but manually overridden to be",df.range.rm[1,coli+6]))
        }
      }
    }else{
      for(rw in seq(1,nrow(df.range.rm))){
        for(coli in seq(1,3)){
          if(df.range.rm[rw,coli+1]!=df.range.rm[rw,coli+6]){
            message(paste("For", df.range.rm$variable[rw],"the",gsub("auto","",colnames(df.range.rm)[coli+1]),"is auto-generated as",df.range.rm[rw,coli+1],", but manually overridden to be",df.range.rm[rw,coli+6]))
          }
        }
      }
    }
  }

  return(df.range)
}

variable_info=function(labstr,liberia.sub,sub.covariates,red.sub.covariates,response.vars,nstddev=n.stdev,block.vars,treat.vars){
  liberia.sub=labelled::remove_labels(liberia.sub)
  liberia.sub=labelled::remove_val_labels(liberia.sub)
  labstr=c(labstr,c("profitsump99avg7d_b"="Unknown","savstockp99_b"="Unknown"))
    df.range=data.frame("variable"=names(labstr),"label"=adjust_label_pre(labstr),
                        "type"=NA,"min"=NA,"max"=NA,
                        "obs.min"=NA,"obs.max"=NA,
                        "subset"=(names(labstr) %in%sub.covariates),
                        "reduced.subset"=(names(labstr) %in% red.sub.covariates),
                        "model.role"="covariate")
    df.range$nVals=sapply(df.range$variable,function(x)length(unique(unlist(liberia.sub[,x]))))
    df.range$obs.min=sapply(df.range$variable,function(x)min(unlist(liberia.sub[,x]),na.rm=T))
    df.range$obs.max=sapply(df.range$variable,function(x)max(unlist(liberia.sub[,x]),na.rm=T))
    df.range$nNA=sapply(df.range$variable,function(x)sum(is.na(unlist(liberia.sub[,x]))))

    #role in model (i.e. response, covariate, block)
    df.range$model.role[(df.range$variable %in%response.vars)]="response"
    df.range$model.role[(df.range$variable %in%c(treat.vars,"control"))]="treatment"
    df.range$model.role[(df.range$variable %in%block.vars)]="block"

    #for block variables,
    df.range$type[(df.range$variable %in%block.vars)]="Categorical"
    df.range$min[(df.range$variable %in%block.vars)]=df.range$obs.min[(df.range$variable %in%block.vars)]
    df.range$max[(df.range$variable %in%block.vars)]=df.range$obs.max[(df.range$variable %in%block.vars)]

    hasindexrange=grepl("[0-9]-[0-9]*",df.range$label)|grepl("[0-9] to [0-9]*",df.range$label)
    df.range$type[((hasindexrange)|(grepl("index",df.range$label,ignore.case=T)))]="Index"

    #print("before index variables")
    #print(head(df.range))

    #for index variables
    df.range=adjust_index_label(df.range=df.range)

    #for binary variables
    df.range=adjust_binary_label(df.range=df.range,liberia.sub=liberia.sub)

    #print("after binary variables")
    #print(head(df.range))

    #for standardized variables
    standardize=((grepl("Standardize",df.range$label,ignore.case = T))|(grepl("_score",df.range$variable))|(grepl("std_",df.range$variable))|(grepl("_z_",df.range$variable))|(grepl("fam_asb_lt",df.range$variable)))
    df.range$type[standardize]="Standardized"
    df.range$min[standardize]=-nstddev
    df.range$max[standardize]=nstddev

    df.range$type[df.range$model.role=="response"]=paste("Average",df.range$type[df.range$model.role=="response"])
    df.range$label=adjust_label_post(df.range$label)

    df.range=manual_adjust(df.range = df.range,nstddev = nstddev,liberia.sub=liberia.sub)

    return(df.range)
}



outside_range=function(df.range,liberia.sub){
  liberia.sub=labelled::remove_labels(liberia.sub)
  liberia.sub=labelled::remove_val_labels(liberia.sub)
  below.min=rep(F,nrow(df.range))
  above.max=rep(F,nrow(df.range))
  below.min[!is.na(df.range$min)]=(df.range$obs.min[!is.na(df.range$min)]<df.range$min[!is.na(df.range$min)])
  above.max[!is.na(df.range$max)]=(df.range$obs.max[!is.na(df.range$max)]>df.range$max[!is.na(df.range$max)])
  out.range.var=df.range$variable[below.min|above.max]
  print(out.range.var)
  df.range$pOut=0
  get_pOut=function(var,lib.df=liberia.sub,df.r=df.range){
    notna=(is.na(unlist(lib.df[,var]))==FALSE)
    sumout=sum(lib.df[notna,var]<df.r$min[df.r$variable==var])
    if(is.na(df.r$max[df.r$variable==var])==FALSE){
      sumout=sumout+sum(lib.df[notna,var]>df.r$max[df.r$variable==var])
    }
    sumout/sum(notna)
  }
  df.range$pOut[df.range$variable%in%out.range.var]=sapply(df.range$variable[df.range$variable%in%out.range.var],get_pOut)
  return(df.range)
}




all.labs=unlist(lapply(liberia.sub,function(x)attr(x,"label")))
tvar=variable_info(labstr=all.labs,
              liberia.sub=liberia.sub,sub.covariates=sub.covariates,
              red.sub.covariates = red.sub.covariates,
              response.vars = response.all,block.vars=block.vars,treat.vars=treat.vars)

var.df=outside_range(tvar,liberia.sub)
var.df[var.df$pOut>0,]
#despite there only being 168 hours in a week some respondents have reported higher hours/week in low-skill business or labor
#despite reporting they are focused on ages 18-35, ages 11 to 44 were sampled
var.df$min[tvar$variable=="age_b"]=5 #DOL supplied statistic about 17.4% of 5-17 year olds work (https://www.dol.gov/sites/dolgov/files/ILAB/Liberia%20CL%20Reference%20Cards%20.pdf)
var.df$max[tvar$variable=="age_b"]=65
var.df=outside_range(var.df,liberia.sub)
View(var.df[(var.df$pOut>0)&(!grepl("Standardize",var.df$type)),])

cat.vars.df=var.df[var.df$type=="Binary"|var.df$type=="Categorical",]
index.vars.df=var.df[var.df$type=="Index",]
need.cat.levels=index.vars.df$variable[index.vars.df$max+1>index.vars.df$nVals]
cat.vars=cat.vars.df$variable
# 5 variables do not take all possible values: "conscientious_b","neurotic_b","rewardresp_b","locuscontr_b","selfesteem_b"
var.df$addstandardize=ifelse(var.df$type=="Continuous",T,
                             ifelse(var.df$type=="Top-Coded Continuous",T,
                                    ifelse((var.df$type=="Discrete")&(var.df$max>20),T,F)))

View(var.df)


bounds.list=lapply(var.df$variable[(!(var.df$model.role%in%c("block","treatment")))&(!(var.df$variable%in%cat.vars))],function(x)c(var.df$min[var.df$variable==x],var.df$max[var.df$variable==x]))
names(bounds.list)=var.df$variable[(!(var.df$model.role%in%c("block","treatment")))&(!(var.df$variable%in%cat.vars))]
bounds.std.list=bounds.list
bounds.std.list[var.df$variable[var.df$addstandardize==T]]=rep(list(c(-n.stdev,n.stdev)),sum(var.df$addstandardize==T))
save(var.df,bounds.list,bounds.std.list,n.stdev,cat.vars,file=paste0(outputpath,"/liberia/liberia_bounds.Rda"))

################## Investigate disputes situation
liberia.full=haven::read_stata(paste0(basepath,"/data/liberia_replication_data/STYL_Final.dta"))

#These nine columns are summed to create disputes column in STYL_Constuction
find.cols=c(paste0("fight",c("smlneigh","bigneigh","phys","weapon","smlleader","bigleader","smlpolice","bigpolice"),"_b"),"fine_b","disputes_all_b")
fights.df=liberia.full[,find.cols]
fights.df$sumfights=rowSums(fights.df[,seq(1,9)],na.rm=T)
tail(fights.df[,c(10,11)])
head(fights.df[fights.df$sumfights>6,])
fights.df$diff=fights.df$disputes_all_b-fights.df$sumfights
dplyr::arrange(fights.df[,c(10,11,12)],diff)
lapply(fights.df,function(x)attr(x,"format.stata"))
apply(fights.df,2,function(x)paste("min",min(x,na.rm=T),"max",max(x,na.rm=T),"mean",mean(x,na.rm=T)))


###########################

response.and.known.vars=c(treat.vars,block.vars,response.all,"control","treatment")
#covariate column labels
x.clabs=unlist(lapply(liberia.sub[,!(colnames(liberia.sub)%in%response.and.known.vars)],
                    function(x)attr(x,"label")))

#data with column name, column label, indicator for subset, indicator for reduced subset,
# number unique values, number of NA
df.covars=data.frame("variable"=names(x.clabs),
                   "label"=x.clabs,
                   "subset"=(names(x.clabs) %in%sub.covariates),
                   "reduced.subset"=(names(x.clabs) %in% red.sub.covariates))
df.covars$nVals=unlist(sapply(liberia.sub[,colnames(liberia.sub)%in%df.covars$variable],
                            function(x)length(unique(x))))
df.covars$nNA=unlist(sapply(liberia.sub[,colnames(liberia.sub)%in%df.covars$variable],
                          function(x)sum(is.na(x))))
#add stars for variables missing observations
df.covars$variable[df.covars$nNA>0]=paste0(df.covars$variable[df.covars$nNA>0],"*")
df.covars$variable[df.covars$nNA>5]=paste0(df.covars$variable[df.covars$nNA>5],"*")

select.cols=c("variable","label","subset","reduced.subset")
#binary variable data
df.vars.bin=df.covars[df.covars$nVals<3,select.cols]
df.vars.bin$label=gsub("Indicator for","",df.vars.bin$label,ignore.case=T)
#not binary variable data
df.vars.cont=df.covars[df.covars$nVals>2,c(select.cols,"nVals")]
#columns with finite values (i.e. index, categorical, discrete bounded variables)
## index indicates it has finite or discrete values, sum is sum of an index,
## hungry is for "how many days have you been hungry in last week",
## members is about number of members in the household
df.index=df.vars.cont[grepl("index|members|sum|hungry",
                            df.vars.cont$label,ignore.case=T)&df.vars.cont$nVals<51,]
df.index=df.index[df.index$variable!="disputes_all_b",] #note sure what to count disputes_all_b as?


#all the other variables
df.vars.other=df.vars.cont[!(df.vars.cont$variable%in%df.index$variable),]
missing=colnames(liberia.sub)[!(colnames(liberia.sub)%in%c(response.and.known.vars,names(x.clabs)))]
if(length(missing)>0){
  miss.vars=data.frame("variable"=missing,
                       "label"=rep("No data label",length(missing)),
                       "subset"=(missing %in%sub.covariates),
                       "reduced.subset"=(missing %in% red.sub.covariates))
  df.vars.other=dplyr::bind_rows(df.vars.other,miss.vars)
}
datas.covars= list("Bin"=df.vars.bin,"Index"=df.index,"Other"=df.vars.other)

covar.tab.list=NULL
cov.red.count=NULL
cov.sub.count=NULL
cov.full.count=NULL
for(j in c(1,2,3)){ #for each variable data frame
  df.name=names(datas.covars)[j]
  df=datas.covars[[j]]
  df=dplyr::arrange(df,variable)
  df=dplyr::arrange(df,desc(subset)) #ordered so reduced subset, subset, and then full
  df=dplyr::arrange(df,desc(reduced.subset))
  row.names(df)=NULL
  cov.full.count=c(cov.full.count,nrow(df))
  cov.sub.count=c(cov.sub.count,sum(df$subset))
  cov.red.count=c(cov.red.count,sum(df$reduced.subset))

  cols.idx=c(1,2)
  cols.nm=c("Variable Name","Meaning")

  var.tab=knitr::kable(df[,cols.idx],format="latex",booktabs=T,
                       col.names=cols.nm)
  #group rows by reduced subset, subset, and full
  n.ingroup=c(sum(df$reduced.subset), #in reduced only
              sum(df$subset)-sum(df$reduced.subset), #in reduced and subset
              nrow(df)-sum(df$subset)) #in only full
  group.name=c("In All Levels of Covariate Sets",
               "In Subset and Full Set of Covariates",
               "Only in Full Set of Covariates")
  starti=1
  for(i in seq(1,3)){ #double check this
    if(n.ingroup[i]-1>0){
    var.tab=kableExtra::pack_rows(var.tab,group.name[i],starti,starti+n.ingroup[i]-1,bold=F,italic=T)
    starti=starti+n.ingroup[i]
    }
  }
  #add color based on set. redsubsetcol and subsetcol defined in latex document
  var.tab=kableExtra::column_spec(var.tab,1,color="black",
                                  background=ifelse(df$reduced.subset==TRUE,"redsubsetcol",ifelse(df$subset==TRUE,"subsetcol","white")))
  covar.tab.list=c(covar.tab.list,list(var.tab))
}


covar.tab.list[[1]]
covar.tab.list[[2]]
covar.tab.list[[3]]
cov.tab.fname=c("binary","count_index","cont")
cov.tab.desc=c("binary indicator covariates","covariates that integers because they are counts or bounded index scores","continuous covariates")
for(i in seq(1,3)){
temp.file=paste0(outputpath,"/tables/liberia_covariate_tables_",cov.tab.fname[i],".tex")

caption.txt=paste("\\caption{There are", cov.full.count[i],cov.tab.desc[i],"in the full baselines set.")
if(cov.sub.count[i]>0){
  caption.txt=paste(caption.txt,"The subset has",cov.sub.count[i],"(in \\textcolor{subsetcol}{light blue}")
if(cov.red.count[i]>0){
  caption.txt=paste(caption.txt," and \\textcolor{redsubsetcol}{dark blue}) and the reduced subset has",cov.red.count[i],"(in \\textcolor{redsubsetcol}{dark blue}).}")
}else{
  caption.txt=paste(caption.txt,").}")
}
}else{
  caption.txt=paste(caption.txt,"}")
}

writeLines(c("\\begin{table}","\\small","\\centering"),temp.file)
CON=file(temp.file,"a")
writeLines(covar.tab.list[[i]], CON)
writeLines(c(
  caption.txt,
  paste0("\\label{tab:cov-",cov.tab.fname[i],"}"),
  "\\end{table}"),CON)
close(CON)
}



#repeat process with response variables
y.clabs=unlist(sapply(liberia.sub[,response.vars.all],
                    function(x)attr(x,"label")))
df.y.vars=data.frame("variable"=names(y.clabs),
                     "label"=y.clabs,
                     "subset"=(names(y.clabs) %in%sub.covariates),
                     "reduced.subset"=(names(y.clabs) %in% red.sub.covariates))
df.y.vars$nVals=unlist(sapply(liberia.sub[,colnames(liberia.sub)%in%df.y.vars$variable],function(x)length(unique(x))))
df.y.vars$nNA=unlist(sapply(liberia.sub[,colnames(liberia.sub)%in%df.y.vars$variable],function(x)sum(is.na(x))))
df.covars$nNA=unlist(sapply(liberia.sub[,colnames(liberia.sub)%in%df.covars$variable],
                             function(x)sum(is.na(x))))
#add stars for variables missing observations
df.y.vars$variable[df.y.vars$nNA>0]=paste0(df.y.vars$variable[df.y.vars$nNA>0],"*")
df.y.vars$variable[df.y.vars$nNA>5]=paste0(df.y.vars$variable[df.y.vars$nNA>5],"*")

rownames(df.y.vars)=NULL
cols.keep=c("variable","label","nVals")
cols.nm=c("Variable Name","Description","Unique Values")
y.var.tab=knitr::kable(df.y.vars[,cols.keep],format="latex",booktabs=T,
                     col.names=cols.nm)
y.var.tab
temp.file=paste0(outputpath,"/tables/liberia_responseVariable_tables.tex")

caption.txt="\\caption{BJS analyzed eight outcome variables. Each is the long term average across 12-13 months. There are four ($\\phantom{.}^*$) and 217 ($\\phantom{.}^{**}$) missing observations from some variables.}\\label{tab:response-vars}"

writeLines(c("\\begin{table}","\\centering"),temp.file)
CON=file(temp.file,"a")
writeLines(y.var.tab, CON)
writeLines(c(
  caption.txt,
  "\\end{table}"),CON)
close(CON)

save(df.y.vars,df.covars,y.var.tab,covar.tab.list,
     file=paste0(outputpath,"/liberia/liberia_variable_set_tables.Rda"))
summary(liberia.sub[,colnames(liberia.sub)%in%df.vars.other$variable])
x.clabs[grepl("dispute",x.clabs)]


