utils::globalVariables(c("ages"))
#' gt
#' 
#' @title gt 
#' 
#' @description Generation time: a measure of the distance between generations rather than the time
#' taken for a population to renew itself. This may be the a mother-daughter distance, i.e. the 
#' average age of mothers at birth of their daughters, or may not take sex into account at all.
#' 
#' @author Laurence Kell, Sea++
#'  
#' @import FLCore
#' @import methods
#' 
#' @name gt
#' @param object \code{FLBRP}
#' @param ... any additional arguments
#' 
#' @aliases gt gt-method gt,FLBRP-method
#' 
#' @export gt
#' @docType methods
#' 
#' @rdname gt
#' 
#' @examples
#' \dontrun{
#' data(ple4)
#' gt(FLBRP(ple4))
#' }
setGeneric('gt', function(object,...) standardGeneric('gt'))

setMethod('gt', signature(object="FLBRP"),
          function(object,f=0,...)
            
            genTimeFn(object,f))

genTimeFn<-function(object,f=0){
  
     object@fbar=object@fbar[,1]%=%f
     n=stock.n(object)[,1]
     
     apply((stock.wt(object)%*%mat(object)%*%n)%*%ages(stock.wt(object)),6,sum)%/%
     apply(stock.wt(object)%*%mat(object)%*%n,6,sum)}

aMsy<-function(object){
  
  object@fbar=object@fbar[,1]%=%c(refpts(object)["msy","harvest"])
  n=stock.n(object)[,1]
  
  apply(stock.wt(object)%*%mat(object)%*%n%*%ages(stock.wt(object)),6,sum)%/%
  apply(stock.wt(object)%*%mat(object)%*%n,6,sum)}
