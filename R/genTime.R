genTime<-function(object,f=0){
  
     object@fbar=object@fbar[,1]%=%f
     n=stock.n(object)[,1]
     
     apply(stock.wt(object)%*%n%*%ages(stock.wt(object)),6,sum)%/%
     apply(stock.wt(object)%*%n,6,sum)}

#growth=FLife:::vonB
aMsy<-function(object){
  
  object@fbar=object@fbar[,1]%=%c(refpts(object)["msy","harvest"])
  n=stock.n(object)[,1]
  
  apply(stock.wt(object)%*%n%*%ages(stock.wt(object)),6,sum)%/%
    apply(stock.wt(object)%*%n,6,sum)}
