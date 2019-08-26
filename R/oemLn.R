utils::globalVariables(c("adply","dnorm","lh","llply","pnorm","rmultinom"))

#' alk
#' 
#' @title alk
#' 
#' @description 
#' @author Laurence Kell, Sea++
#' 
#' @name alk
#' @param object \code{FLPar} with `linf` or  \code{numeric}
#' @param ... any additional arguments
# #' @param age \code{numeric} 0:40 by default
# #' @param cv  \code{numeric} 0.2 by default
# #' @param lmax  \code{numeric} maximum size relative to 'linf' 1.2 by default
#'  
#' @aliases alk alk-method alk,FLPar-method
#' 
#' @docType methods
#' @export alk
#'  
#' @rdname alk
#' @seealso lenSample
#' 
#' @examples
#' \dontrun{
#' ak=alk(FLPar(linf=50))
#' }
setGeneric('alk', function(object,...) standardGeneric('alk'))

setMethod('alk', signature(object="FLPar"),
          function(object,age=0:40,cv=0.2,lmax=1.2,...)
            
            setALKFn(object,age,cv,lmax))

setALKFn<-function(par,age=0:40,cv=0.2,lmax=1.2){
  
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

#' lenSample
#' 
#' @title lenSample 
#' 
#' @description 
#' @author Laurence Kell, Sea++  
#' 
#' @name lenSample
#' 
#' @param object \code{FLStock} blah,blah,blah,...
#' @param alk \code{FLPar} with ALK
#' @param nsample \code{numeric} sample size, default is 500
#' @param ... any additional arguments
#' 
#' @docType methods
#' 
#' @export lenSample
#' @rdname lenSample
#' @seealso setALK
#' 
#' @aliases lenSample  
#'          lenSample-method  
#'          lenSample,FLQuant,FLPar,missing-method
#'          lenSample,FLQuant,FLPar,numeric-method
#' 
#' @examples
#' \dontrun{
#' lfd=lenSample(catch.n(ple4)[,ac(2000:2005)],alk,nsample=100)
#' }
setGeneric('lenSample', function(object,alk,nsample,...) standardGeneric('lenSample'))

setMethod('lenSample', signature(object="FLQuant",alk="FLPar",nsample="missing"),
          function(object,alk,nsample=500,...)
            
            lenSampleFn(object,alk,nsample))
setMethod('lenSample', signature(object="FLQuant",alk="FLPar",nsample="numeric"),
          function(object,alk,nsample,...)
            
            lenSampleFn(object,alk,nsample))

lenSampleFn<-function(object,alk,nsample){
  
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
  alk=alk(FLPar(lh[,1]))
  lfd=lenSample(stock.n(om)[,95:100,,,,1:2],alk,nsample=5000)
    
  ggplot(melt(lfd))+
    geom_histogram(aes(Var.3,weight=value),binwidth=1)+
    facet_grid(year~iter)+
    xlab("Length (cm)")+ylab("Frequency")+
    scale_x_continuous(limits=c(0,45))  
  }
