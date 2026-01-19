### Function to Check Regression Assumptions ###
reg_assumptions=function(reg.model,data,family="gaussian",response.var=NULL,mod.name=NULL,
                         pt.sz=7,ln.sz=3,bs.sz=50,long.plot.title=T,...){

  mod.fit=stats::glm(reg.model,data=data,family=family,... )
  if(is.null(response.var)==TRUE){
    response.var=names(attr(mod.fit$terms,"dataClasses"))[1]
  }
  if(is.null(mod.name)==TRUE){
    mod.name=response.var
  }
  if(long.plot.title==T){
    fit.lab=paste("Fitted",response.var)
    fit.res.title=paste0(mod.name,": Residual vs. Fitted")
    qq.title=paste0(mod.name,": Normal Q-Q Plot")
  }else{
    fit.lab=paste("Fitted y")
    fit.res.title=mod.name
    qq.title=mod.name
  }
  data.temp=data.frame("fitted"=mod.fit$fitted,"residuals"=mod.fit$residuals,"zeros"=0)
  plot.fit.resid=ggplot2::ggplot(data.temp,ggplot2::aes(x=fitted,y=residuals))+
    ggplot2::geom_point(size=pt.sz)+
    ggplot2::geom_smooth(se=F,col="red",linewidth=ln.sz,method = 'loess',formula = 'y ~ x')+
    ggplot2::geom_line(ggplot2::aes(x=fitted,y=zeros),col="black",linetype=2,linewidth=ln.sz)+
    ggplot2::labs(x=fit.lab,y="Residuals")+
    ggplot2::ggtitle(fit.res.title)+
    ggplot2::theme_minimal(base_size=bs.sz)

  data.temp$norm=stats::rnorm(nrow(data.temp))
  plot.qqnorm=ggplot2::ggplot(data=data.temp,ggplot2::aes(sample=residuals))+
    ggplot2::geom_qq_line(color="red",linewidth=ln.sz)+
    ggplot2::geom_qq(size=pt.sz,color="black")+
    ggplot2::labs(x="Theoretical Quantiles",y="Residuals Quantiles")+
    ggplot2::ggtitle(qq.title)+
    ggplot2::theme_minimal(base_size=bs.sz)


  #variance.test=car::ncvTest.glm(mod.fit)
  return(cowplot::plot_grid(plot.fit.resid,plot.qqnorm,ncol=2))
}

