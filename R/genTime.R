utils::globalVariables(c("ages"))
#' gt
#' 
#' @title gt 
#' 
#' @description Calculates the generation time
#' @author Laurence Kell, Sea++
#'  
#' @import FLCore
#' @import methods
#' 
#' @name gt
#' @param object \code{FLBRP}
#' @param ... other arguments
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
#' data(pl4)
#' alk=gt(FLBRP(ple4))
#' }
setGeneric('gt', function(object,...) standardGeneric('gt'))

setMethod('gt', signature(object="FLBRP"),
          function(object,f=0,...)
            
            genTimeFn(object,f))

genTimeFn<-function(object,f=0){
  
     object@fbar=object@fbar[,1]%=%f
     n=stock.n(object)[,1]
     
     apply(stock.wt(object)%*%n%*%ages(stock.wt(object)),6,sum)%/%
     apply(stock.wt(object)%*%n,6,sum)}

aMsy<-function(object){
  
  object@fbar=object@fbar[,1]%=%c(refpts(object)["msy","harvest"])
  n=stock.n(object)[,1]
  
  apply(stock.wt(object)%*%n%*%ages(stock.wt(object)),6,sum)%/%
    apply(stock.wt(object)%*%n,6,sum)}
