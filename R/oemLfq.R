utils::globalVariables(c("dnorm","age","unit","season","area"))

### Create lengths-at-age 
oemLfq<-function(om,lh,
                 n=catch.n(om),ln=vonB(ages(catch.n(om)),lh),sd =ln*0.2,bin=0:ceiling(max(ln)*1.10)+0.5){
  
  sim=function(ln,sd,n,bin) {data.frame(length=bin,data=dnorm(bin,ln,sd)*n)}
  
  lfq=ddply(model.frame(FLQuants(ln=ln,sd=sd,n=n)),.(age,year,unit,season,area,iter), 
            with, sim(ln,sd,n,bin), bin=bin)
  
  ### sum up by length 
  ddply(lfq,.(length,year,unit,season,area,iter), 
        with, data.frame(freq=sum(data)))}
