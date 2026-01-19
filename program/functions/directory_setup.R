directory_setup<-function(basepath=NULL,outputpath=NULL,out.subfolders=list("figures","simulations","tables")){
  if(is.null(basepath)==TRUE){
    basepath = rprojroot::find_rstudio_root_file()
  }
  if(is.null(outputpath)==TRUE){
    outputpath=paste0(basepath,"/output")
  }
  subf.none=is.null(out.subfolders)
  if(subf.none==FALSE){
  subpaths=lapply(out.subfolders,function(sub)paste0(outputpath,"/",sub))
  }else{
    subpaths=NULL
  }
  if(!dir.exists(outputpath)){
    dir.create(outputpath)
    if(subf.none==FALSE){
      for(subp in subpaths){
        dir.create(subp,recursive=T)
      }
    }
  }else{
    if(subf.none==FALSE){
    for(subp in subpaths){
      if(!dir.exists(subp)){
        dir.create(subp,recursive=T)
      }
    }
    }
  }
  outlist=c(list(basepath,outputpath),subpaths)
  if(subf.none==FALSE){
    names(outlist)=c("basepath","outputpath",sapply(out.subfolders,function(x)paste0(x,"path")))
  }else{
    names(outlist)=c("basepath","outputpath")
  }
  return(outlist)
}
