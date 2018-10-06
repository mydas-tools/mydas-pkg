setALK<-function(par,age=0:40,cv=0.2,lmax=1.2){
  
  alk=mdply(data.frame(age=age), function(age,par,cv,lmax){
    res=adply(par,2,function(x,age,cv,lmax){
        ln =seq(1:ceiling(x["linf"]*lmax))
        len=vonB(age,x)
        p  =c(pnorm(1,len,len*cv),
              dnorm(ln[-c(1,length(ln))],len,len*cv),
              pnorm(ln[length(ln)],len,len*cv,lower.tail=FALSE))
  
        data.frame(len=ln,p=p/sum(p))},
      age=age,cv=cv,lmax=lmax)},par=lh,cv=cv,lmax=lmax)
  
  alk=cast(alk,age~len~iter,fun=sum,value="p")
  alk=FLPar(array(alk,dim=unlist(llply(dimnames(alk),length)),dimnames=dimnames(alk)))
  
  alk}

lenSample<-function(object,alk,nsample){
  
  res=mdply(expand.grid(iter=seq(dim(object)[6]),year=seq(dim(object)[2])),
         function(iter,year){ 
            lfd=object[,year,,,,iter]%*%alk[,,iter,drop=T]
            data.frame(length=dimnames(lfd)[[2]],data=apply(rmultinom(nsample,1,prob=lfd),1,sum))})

  res=as.FLQuant(res)
  dimnames(res)[2:6]=dimnames(object)[2:6]
  res}

if (FALSE){
  library(FLCore)
  library(FLife)
  library(plyr)
  library(dplyr)

  load("/home/laurence/Desktop/sea++/mydas/tasks/task4/data/turbot.RData")
  alk=setALK(FLPar(lh[,1]))
  lfd=lenSample(stock.n(om)[,95:100,,,,1:2],alk,nsample=5000)
    
  ggplot(melt(lfd))+
    geom_histogram(aes(Var.3,weight=value),binwidth=1)+
    facet_grid(year~iter)+
    xlab("Length (cm)")+ylab("Frequency")+
    scale_x_continuous(limits=c(0,45))  
  }
