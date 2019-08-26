utils::globalVariables(c("mdply","cast","melt","vonB","LBSPRfit","difference.sign.test","lm"))

#' lbspr
#' 
#' @title lbspr 
#' 
#' @description A wrapper for LBSPR  see https://cran.r-project.org/web/packages/LBSPR/vignettes/LBSPR.html
#'  This alloe fitting to empirical length data to provide an estimate of the spawning potential 
#'  ratio (SPR) for use in MSE.
#' 
#' @author Laurence Kell, Sea++
#'  
#' @param object \code{data.frame} with length frequency sample
#' @param params \code{FLPar} with values for linf, a50, ato95 and mk
#' @param ... any additional arguments
#' 
#' @aliases lbspr 
#'          lbspr-method 
#'          lbspr,data.frame,FLPar-method
#'          
#' @docType methods
#' 
#' @rdname lbspr
#' 
#' @export lbspr
#' 
#' @examples
#' \dontrun{
#' lfd=lenSample(catch.n(om)[,20:65],ak,nsample=500)
#' lb=lbspr(lfd,prior) 
#' }
setGeneric('lbspr', function(object,params,...) standardGeneric('lbspr'))

lbsprFn<-function(len,params,species="",units="cm"){
  
  pars        =new("LB_pars")
  pars@Linf   =c(params["linf"]) 
  pars@L50    =vonB(c(params["a50"]),params) 
  pars@L95    =vonB(c(params["a50"])+c(params["ato95"]),params)
  pars@MK     =c(params["mk"])
  pars@Species=species
  pars@L_units=units
  
  #labs=dimnames(len)[[1]]
  #brks=cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ),
  #           upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))
  #mid=aaply(brks,1,mean)
  
  LBlen       =new("LB_lengths")
  LBlen@LMids =as.numeric(dimnames(len)[[1]])
  LBlen@LData =len
  LBlen@Years =as.numeric(dimnames(len)[[2]])
  LBlen@NYears=dim(len)[2] 
  
  res=LBSPRfit(pars,LBlen,verbose=FALSE)
  
  res@Ests}


setMethod('lbspr', signature(object="data.frame",params='FLPar'),
          function(object,params,...){
  
  nits=max(dim(object)[6],dim(params)[2])
  
  if (!(dim(object)[6]%in%c(1,nits)|(dim(params)[2]%in%c(1,nits))))
    stop("iters should be equal to n or 1")
  
  res=mdply(data.frame(iter=seq(nits)), function(iter)
    lbsprFn(iter(object,iter)[drop=T],iter(params,iter)))
  res=data.frame(year=dimnames(object)$year,res)
  
  rtn=FLPar(cast(melt(res,id=c("year","iter")),variable~year~iter),units="NA")
  
  rtn})
