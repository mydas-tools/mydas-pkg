genTime<-function(object){
  
     fbar(object)=fbar(object)[,1]%=%0
     n=stock.n(object)[,1]
     
     apply(stock.wt(object)%*%n%*%ages(stock.wt(object)),6,sum)%/%
     apply(stock.wt(object)%*%n,6,sum)}
