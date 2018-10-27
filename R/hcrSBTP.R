hcrSBTP<-function(yrs,
                  control,
                  ref,
                  target,
                  cpue,
                  catch,
                  wt=0.5){
  #SBT  
  coef       =1+control["k1"]
  flag       =c(cpue<ref)
  coef[,flag]=1-control["k2",flag]
  bit        =target%*%((cpue%/%ref)%^%(coef))
  res        =catch*(1-wt)+bit*(wt)

  dmns=dimnames(catch)
  dmns$year=yrs
  dmns$iter=dimnames(cpue)$iter

  res=FLQuant(rep(c(res),each=length(yrs)),dimnames=dmns)

  return(res)}