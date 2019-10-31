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

    #zz=file(tempfile(),open="wt");#sink(zz);sink(zz,type="message")
    res=LBSPRfit(pars,LBlen,verbose=FALSE)
    #sink()
    
    res@Ests}
  
setMethod('lbspr', signature(object="FLQuant",params='FLPar'),
    function(object,params,...){
            
      nits=max(dim(object)[6],dim(params)[2])
            
      if (!(dim(object)[6]%in%c(1,nits)|(dim(params)[2]%in%c(1,nits))))
          stop("iters should be equal to n or 1")
       
      yrs=dimnames(object)$year
      
      res=mdply(data.frame(iter=seq(nits)), function(iter)
          lbsprFn(iter(object,iter)[drop=T],iter(params,iter)))
            res=data.frame(year=dimnames(object)$year,res)

      res$year=as.numeric(as.character(res$year))  
      
      rtn=FLQuants(dlply(melt(res,id=c("year","iter")), .(variable), 
            with, as.FLQuant(data.frame(year=as.numeric(as.character(year)),
                                        iter=as.numeric(as.character(iter)),
                                        data=value),unit="NA")))

      names(rtn)=c("sl50","sl95","fm","spr")
      
      FLQuants(rtn)})

setMethod('lbspr', signature(object="data.frame",params='FLPar'),
          function(object,params,...){
            
    flq=as.FLQuant(transmute(object,
                             len =as.numeric(as.character(len)),
                             year=as.numeric(as.character(year)),
                             iter=as.numeric(as.character(iter)),
                             data=data))
            
    flq=dat[order(as.numeric(dimnames(flq)[[1]]))]
    flq[is.na(flq)]=0
            
    lbspr(flq,params)})

setGeneric('fwd.lbspr', function(object,fm,...) standardGeneric('fwd.lbspr'))

lbsprFwdFn<-function(params,fm,binwidth=1,species="",units="cm"){
  
  if (!("ato95"%in%dimnames(params)$params))
    stop("ato95 needs to be in params")
  if (!("mk"%in%dimnames(params)$params))
    stop("mk needs to be in params")
  if (!any(c("a1","sel3")%in%dimnames(params)$params))
    stop("a1 or sel3 needs to be in params")
  if (!any(c("a50","l50")%in%dimnames(params)$params))
    stop("a1 or sel3 needs to be in params")
  
  # New object
  warn=options()$warn
  options(warn=-5)
  pars=new("LB_pars")
  options(warn=warn)
  
  # Set parameters
  pars      =new("LB_pars")
  pars@Linf =c(params["linf"])
  if ("a50"%in%dimnames(params)$params & !("l50"%in%dimnames(params)$params))
    pars@SL50 =c(params["l50"])
  else
    pars@SL50 =vonB(c(params["sel1"]),params)
  pars@L50       =vonB(c(params["a50"]),params)
  pars@L95       =pars@L50+vonB(c(params["ato95"]),params)
  pars@MK        =c(params["mk"])
  if ("a1"%in%dimnames(params)$params & !("sel1"%in%dimnames(params)$params))
    pars@SL50      =vonB(c(params["a1"]),params)
  else
    pars@SL50      =vonB(c(params["sel1"]),params)
  pars@SL95      =pars@L95
  pars@BinWidth  =binwidth
  pars@Species   =species
  pars@L_units   =units
  
  
  pars@BinMax=pars@Linf*1.5
  pars@BinMin=  0
  
  pars@FM=fm # set value for FM
  res    =LBSPRsim(pars)
  res    =LBSPRsim(pars)
  
  res=FLPar(apply(res@pLPop[,2:5]*res@pLPop[,1],2,sum))
  units(res)=units
  res}

setMethod('fwd.lbspr', signature(object='FLPar'),
          function(object,fm,...){
            
          if (dim(object)[2]==1)
              FLPar(lbsprFwdFn(object,fm=fm))
            else
              FLPar(t(aaply(object,2,function(x) lbsprFwdFn(FLPar(x),fm=fm))))
          })






