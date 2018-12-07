source('~/Desktop/sea++/mydas/pkg/R/genTime.R')
gTime=c(round(mydas:::genTime(eq)))
genTime<-function(object,f=0)
  

source('~/Desktop/sea++/mydas/pkg/R/oemLn.R')
alk=mydas:::setALK(lh)  
setALK<-function(par,age=0:40,cv=0.2,lmax=1.2){
  
  
lfd=mydas:::lenSample(catch.n(om)[,20:65],alk,nsample=500)
lenSample<-function(object,alk,nsample){
  
source('~/Desktop/sea++/mydas/pkg/R/omOut.R')
ts   =mydas:::omSmry(om,eq,lh)
omSmry<-function(x,y="missing",z="missing"){
  
source('~/Desktop/flr/mpb/R/biodyn-sra.R')

source('~/Desktop/sea++/mydas/pkg/R/mlz.R')
res=mydas:::mlz(mnLen[,ac(40:60)],prior)
mlz<-function(object,params,ss=500){
  
source('~/Desktop/sea++/mydas/pkg/R/lbspr.R')
lb=mydas:::lbspr(lfd,prior)  
lbspr<-function(object,params){
  
source('~/Desktop/sea++/mydas/pkg/R/hcr.R')
source('~/Desktop/sea++/mydas/pkg/R/mseXSA.R')
source('~/Desktop/sea++/mydas/pkg/R/hcrSBTD.R')
source('~/Desktop/sea++/mydas/pkg/R/mseMPB.R')

prior=FLife:::priors(lh)
priors<-function(object,eq=lhEql(lhPar(object))){
  
mp=mpb:::setMP(as(window(om,end=55),"biodyn"),
               r =median(prior["r"],na.rm=T),
               k =median(prior["v"],na.rm=T),
               b0=0.8,
               p =median(mpb:::p(prior["bmsy"]/prior["v"]),na.rm=TRUE))
setMP<-function(om,r,k,p,b0){
  
  
p =median(mpb:::p(prior["bmsy"]/prior["v"]),na.rm=TRUE))
p<-function(from){
  
srDev=FLife:::rlnoise(nits,FLQuant(0,dimnames=list(year=1:100)),0.3,b=0.0)
setMethod("rnoise", signature(n='numeric', len="missing"),
          function(n=n,len=len,sd=0.3,b=0,burn=0,trunc=0,what=c("year","cohort","age")) {
            noiseFn(n,sd,b,burn,trunc)})

setMethod("rlnoise", signature(n='numeric', len="FLQuant"),
          function(n=n,len=len,sd=0.3,b=0,burn=0,trunc=0,what=c("year","cohort","age")) {
            exp(rnoise(n,len,sd,b,burn,trunc,what))})
