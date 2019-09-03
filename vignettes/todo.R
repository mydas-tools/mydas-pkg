#source('~/Desktop/flr/mpb/R/biodyn-sra.R')
#source('~/Desktop/sea++/mydas/pkg/R/mlz.R')

res=mydas:::mlz(mnLen[,ac(40:60)],prior)
mlz<-function(object,params,ss=500){
#' setALK
#' 
#' @title setALK 
#' 
#' @description 
#' @author Laurence Kell, Sea++
#'  
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @export 
#' @docType methods
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' alk=mydas:::setALK(lh)
#' }
  
source('~/Desktop/sea++/mydas/pkg/R/lbspr.R')
lb=mydas:::lbspr(lfd,prior)  
lbspr<-function(object,params){
#' setALK
#' 
#' @title setALK 
#' 
#' @description 
#' @author Laurence Kell, Sea++
#'  
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @export
#' @docType methods
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' alk=mydas:::setALK(lh)
#' }
  
source('~/Desktop/sea++/mydas/pkg/R/hcr.R')
source('~/Desktop/sea++/mydas/pkg/R/mseXSA.R')
source('~/Desktop/sea++/mydas/pkg/R/hcrSBTD.R')
source('~/Desktop/sea++/mydas/pkg/R/mseMPB.R')

prior=FLife:::priors(lh)
priors<-function(object,eq=lhEql(lhPar(object))){
#' setALK
#' 
#' @title setALK 
#' 
#' @description 
#' @author Laurence Kell, Sea++
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @export
#' @docType methods
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' alk=mydas:::setALK(lh)
#' }
  
mp=mpb:::setMP(as(window(om,end=55),"biodyn"),
               r =median(prior["r"],na.rm=T),
               k =median(prior["v"],na.rm=T),
               b0=0.8,
               p =median(mpb:::p(prior["bmsy"]/prior["v"]),na.rm=TRUE))
setMP<-function(om,r,k,p,b0){
#' setALK
#' 
#' @title setALK 
#' 
#' @description 
#' @author Laurence Kell, Sea++
#' 
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @export
#' @docType methods
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' alk=mydas:::setALK(lh)
#' }
  
  
  
srDev=FLife:::rlnoise(nits,FLQuant(0,dimnames=list(year=1:100)),0.3,b=0.0)
setMethod("rnoise", signature(n='numeric', len="missing"),
          function(n=n,len=len,sd=0.3,b=0,burn=0,trunc=0,what=c("year","cohort","age")) {
            noiseFn(n,sd,b,burn,trunc)})

setMethod("rlnoise", signature(n='numeric', len="FLQuant"),
          function(n=n,len=len,sd=0.3,b=0,burn=0,trunc=0,what=c("year","cohort","age")) {
            exp(rnoise(n,len,sd,b,burn,trunc,what))})
#' setALK
#' 
#' @title setALK 
#' 
#' @description 
#' @author Laurence Kell, Sea++
#' @name 
#' @param 
#' 
#' @aliases
#' 
#' @export
#' @docType methods
#' 
#' @rdname 
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' alk=mydas:::setALK(lh)
#' }
