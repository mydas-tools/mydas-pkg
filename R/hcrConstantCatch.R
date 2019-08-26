#' hcrConstantCatch
#' 
#' @title hcrConstantCatch 
#' 
#' @description 
#' @author Laurence Kell, Sea++ 
#' 
#' @param yrs \code{numeric} target years to take catch over
#' @param catch an \code{FLQuant} with catch to average over by iter
#' @param ... any additional arguments
#'  
#' @name hcrConstantCatch
#' 
#' @export hcrConstantCatch
#' @docType methods
#' 
#' @rdname hcrConstantCatch
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