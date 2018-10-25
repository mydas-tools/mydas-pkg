hcrSBTP<-function(yrs,
                  control=FLPar(c(k1=0.25,k2=0.25)),
                  ref,
                  target,
                  cpue,
                  catch){
  
  coef       =-control["k1"]
  flag       =cpue>ref
  coef[,flag]=control["k2",flag]
  
  bit     =target%*%(cpue%/%ref)%^%coef
  res     =(catch+bit)

  dmns=dimnames(catch)
  dmns$year=yrs
  dmns$iter=dimnames(cpue)$iter
  
  res=FLQuant(rep(c(res),each=length(yrs)),dimnames=dmns)
  
  return(res)}