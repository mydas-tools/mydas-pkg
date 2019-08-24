#' hcr
#' 
#' @title hcr 
#' 
#' @description 
#' @author Laurence Kell, Sea++
#'  
#' @name hcr
#' 
#' @export hcr
#' @docType methods
#' 
#' @rdname hcr
#' 
#' @examples
#' \dontrun{
#' data(pl4)
#' }
#' 
hcrConstantCatch<-function(yrs,catch,...){
  res=FLQuant(c(apply(catch,6,mean)), 
              dimnames=list(year=yrs,iter=dimnames(catch)$iter)) 
  res}