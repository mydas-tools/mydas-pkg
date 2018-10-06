hcrSBT1<-function(yrs,
                  control=c(k1=2.0,k2=3.0,gamma=1),
                  index,
                  catch,...){
  
  lambda=as.FLQuant(ddply(as.data.frame(index%/%apply(index,6,mean)), 
                          .(iter), with,  data.frame(data=coefficients(lm(data~year))[2])))
  
  flag  =lambda<0
  lambda=abs(lambda)
  res   =1+ifelse(flag,-control["k1"],control["k2"])*exp(log(lambda)*ifelse(flag,control["gamma"],1))
  res   =res%*%apply(catch,6,mean)
  
  dmns     =dimnames(catch)
  dmns$year=yrs
  dmns$iter=dimnames(index)$iter
  
  res      =FLQuant(rep(c(res),each=length(yrs)),dimnames=dmns)
  
  return(res)}

hcrSBT2=function(yrs,
                  control=c(k1=0.25,k2=0.25),
                  catch,
                  cpue,
                  ref,
                  target){
  
  flag    =cpue<ref
  bit     =target*(cpue/ref)*(1+ifelse(flag,-control[1],control[2]))
  res     =(catch+bit)/2
  
  dmns=dimnames(catch)
  dmns$year=yrs
  dmns$iter=dimnames(cpue)$iter
  
  res=FLQuant(rep(c(res),each=length(yrs)),dimnames=dmns)
  
  return(res)}

hcrSBT2orig=function(adult,  juve,
                     yrAdult,yrJuve,
                     refJuve=-(1:5),
                     tac,tarCatch,
                     k1=0.25,k2=0.75,lag=1,interval=3){
  
  adultIdx=adult[,ac(dims(adult)$maxyear)]
  adultRef=aaply(adult[,ac(yrAdult)],3:6,mean)
  flag    =adultIdx<adultRef
  cBit    =tarCatch*(adultIdx/adultRef)*(1+ifelse(flag,-k1,k1))
  
  juveIdx =aaply(juve[,ac(dims(juve)$maxyear+refJuve)],3:6,mean)
  juveRef =aaply(juve[,ac(yrJuve) ],3:6,mean)
  flag    =juveIdx<juveRef
  rBit    =(juveIdx/juveRef)*(1+ifelse(flag,-k2,k2))
  
  # cat('ref Juve:',   as.integer(mean(refJuve)),
  #     '\t Juve:',    as.integer(mean(juve)),
  #     '\t ratio:',   mean(juve/refJuve),
  #     '\t rBit:',    mean(rBit),'\n')
  
  res =0.5*(tac+cBit*rBit)
  
  #   cat('TAC:',        as.integer(mean(tac)),
  #       '\t ratio:',   as.integer((mean(adult/refAdult))),
  #       '\t delta:',   as.integer((mean(cBit))),
  #       '\t New TAC:', as.integer(mean(res)),
  #       '\t rBit:',    mean(rBit),'\n')
  
  dmns=dimnames(tac)
  dmns$year=as.character(as.integer(dmns$year)+lag+seq(interval)-1)
  dmns$iter=dimnames(adult)$iter
  
  res=FLQuant(rep(c(res),each=interval),dimnames=dmns)
  
  return(res)}

hcrSBT2orig=function(adult,juve,
                 yrAdult,yrJuve,
                 refJuve=-(1:5),
                 tac,tarCatch,k1=0.25,er=0.75,lag=1,interval=3){
  
  adultIdx=adult[,ac(dims(adult)$maxyear)]
  adultRef=aaply(adult[,ac(yrAdult)],3:6,mean)
  flag    =adultIdx<adultRef
  cBit    =tarCatch*(adultIdx/adultRef)*(1+ifelse(flag,-k1,k1))
  
  juveIdx =aaply(juve[,ac(dims(juve)$maxyear+refJuve)],3:6,mean)
  juveRef =aaply(juve[,ac(yrJuve) ],3:6,mean)
  flag    =juveIdx<juveRef
  rBit    =(juveIdx/juveRef)*(1+ifelse(flag,k2,-k2))
  
  # cat('ref Juve:',   as.integer(mean(refJuve)),
  #     '\t Juve:',    as.integer(mean(juve)),
  #     '\t ratio:',   mean(juve/refJuve),
  #     '\t rBit:',    mean(rBit),'\n')
  
  res =0.5*(tac+cBit*rBit)
  
  #   cat('TAC:',        as.integer(mean(tac)),
  #       '\t ratio:',   as.integer((mean(adult/refAdult))),
  #       '\t delta:',   as.integer((mean(cBit))),
  #       '\t New TAC:', as.integer(mean(res)),
  #       '\t rBit:',    mean(rBit),'\n')
  
  dmns=dimnames(tac)
  dmns$year=as.character(as.integer(dmns$year)+lag+seq(interval)-1)
  dmns$iter=dimnames(adult)$iter
  
  res=FLQuant(rep(c(res),each=interval),dimnames=dmns)
  
  return(res)}
