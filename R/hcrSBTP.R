hcrSBTP<-function(yrs,
                  control=c(k1=0.25,k2=0.25),
                  catch,
                  cpue,
                  ref,
                  target){
  
  flag    =cpue<ref
  bit     =target*(cpue/ref)^(1+ifelse(flag,-control[1],control[2]))
  res     =(catch+bit)/2
  
  dmns=dimnames(catch)
  dmns$year=yrs
  dmns$iter=dimnames(cpue)$iter
  
  res=FLQuant(rep(c(res),each=length(yrs)),dimnames=dmns)
  
  return(res)}