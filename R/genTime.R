#' gt
#' 
#' @title gt 
#' 
#' @description 
#' @author Laurence Kell, Sea++
#'  
#' @name gt
#' @param object \code{FLBRP}
#' 
#' @aliases
#' 
#' @export gt
#' @docType method
#' 
#' @rdname gt
#' @seealso 
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
