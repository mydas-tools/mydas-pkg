#Safety, is stock above PRI
#Time to recovery to green quandrant and once  recovered does it stay there?
#Discounted yield
#AAV in yield

avFn<-function(object){
  
  o1 =object[-length(object)]
  o2 =object[-1]
  
  return(abs((o2-o1)/o1))}

green<-function(flag,year){
  flag=flag>=1
  year=as.numeric(year)
  
  if (sum(flag)==0){
    y=max(year)+1
    p=0
    n=0
  }else{
    y   =year[flag][1]
    flag=flag[year>y]
    n   =length(flag[flag])
    p   =mean(flag)}
  
  return(data.frame(year=y,n=n,p=p))}

udFn<-function(object){
  
  o1 =object[-1]
  o2 =object[-length(object)]
  
  return(sum(o1<o2))/(length(object)-1)}

dRate=function(x,r,wtAv=FALSE) {
  if (wtAv) 
     return( sum(x/(1+r)^(0:(length(x)-1)))/sum(1/(1+r)^(0:(length(x)-1))))
 return(sum(x/(1+r)^(0:(length(x)-1))))
 }

dRate2=function(x,r) {
  wt =1/((1+r)^(0:(length(x)-1)))
  res=sum(x*wt)/sum(wt)
  res}

trend=function(x) {
  #x=x%/%apply(x,6,mean)
  #as.FLQuant(ddply(as.data.frame(x), .(unit,season,area,iter), with,
  #          data.frame(data=lm(data~year)$coefficients[2])))
  
  x=x/mean
  data.frame(data=lm(x~seq(length(x))$coefficients[2]))
  }

randomness=function(obj) {
  #as.FLQuant(ddply(as.data.frame(obj), .(unit,season,area,iter), with,
  #                 data.frame(data=difference.sign.test(data)[["p.value"]])))
  difference.sign.test(obj)[["p.value"]]}

smryStat<-function(dat,dr=0){
  with(dat, data.frame(safety  =min(rec_hat/virgin_rec,na.rm=T),
                       blim    =mean(ssb/virgin_ssb,na.rm=T),
                       ftar    =min(fbar/msy_harvest,na.rm=T),
                       btar    =min(ssb/msy_ssb,na.rm=T),
                       kobe    =green((ssb/msy_ssb>1)&(fbar/msy_harvest<1),year)[1:2],
                       yield   =dRate(catch/msy_yield,dr)/length(catch),
                       yieldAav=mean(avFn(pmax(catch,msy_yield*0.1))),
                       yieldRnd=randomness(catch),
                       ssbRnd  =randomness(ssb),
                       fbRnd   =randomness(fbar)))}
