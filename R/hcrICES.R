#' hcrICES
#' 
#' @title hcrICES 
#' 
#' @description Harvest Control Rule, calculates Total Allowable Catch (TAC) based on a hockey stock harvest control rule.
#' @author Laurence Kell, Sea++
#'  
#' @usage hcrICES(object,eql,sr_deviances,params, 
#'          start=max(dimnames(object)$year)-10, end=start+10, interval=1,
#'          err=NULL,bndTac=c(0,Inf),...)
#'          
#' @param object an object of class \code{FLStock} 
#' @param eql \code{FLBRP} with a stock recruitment relationship used for projection
#' @param sr_deviances \code{FLQuant} recuitment deviates on the log scale, i.e. multiplicative
#' @param params \code{FLPar} HCR parameters, specifying blim, btrig, ftar and fmin
#' @param start \code{numeric} first year for simulation
#' @param end \code{numeric}   last year for simukation
#' @param interval  \code{numeric} time step, 1 year by default
#' @param err \code{FLQuant} assessment error on SSB for year used for HCR
#' @param bndTac \code{numeric} bounds on TAC, by default these are turned off, for 20 percent constraint set to c(0.8,1.2)
#' @param ... any additional arguments
#' 
# #' @aliases hcrICES hcrICES-method hcrICES,FLStock,FLBRP-method
#' 
#' @export hcrICES
#' @docType methods
#' 
#' @rdname hcrICES
#' 
#' @return returns a \code{list} with \code{FLStock} and \code{FLPar} objects for the stock and HCR 
#'
#' @import FLCore 
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' data(pl4)
#' }
#' 
#' 
setGeneric('hcrICES', function(object,eql,...) standardGeneric('hcrICES'))

setMethod('hcrICES', signature(object="FLStock",eql='FLBRP'),
          function(object,eql,sr_deviances,params,
                   start   =max(dimnames(object)$year)-10,
                   end     =start+10,
                   interval=1,
                   err     =NULL,
                   bndTac  =c(0,Inf),...){
        
            if ("perfect"%in%names(list(...)))
              perfect=list(...)$perfect
            else 
              perfect=FALSE
            
            bias<-function(object,eql,err,iYr,lag=1){
              stock.n(object)=stock.n(object)%*%err[,ac(iYr)]
              
              if (lag>0)
                object=fwd(object,catch=catch(object)[,ac(iYr-(lag:0))],sr=eql)
              window(object,end=iYr+1)}
            
            chk=NULL
            ## Loop round years
            cat('\n==')
            for (iYr in seq(start,end-1,interval)){
              cat(iYr,", ",sep="")
              
              stkYrs=iYr
              refYrs=iYr
              hcrYrs=iYr+1
              
              refpts(eql)[]=1
              refpts(eql)=propagate(refpts(eql),dim(object)[6])
              
              if (!is.null(err))
                mp=bias(object,eql,err,iYr,lag=hcrYrs-stkYrs-1)
              else
                mp=object
              
              res=hcrFn(mp,eql,params,
                        stkYrs,
                        refYrs,
                        hcrYrs,
                        perfect=perfect,
                        sr_deviances=sr_deviances)
              
              ##TAC Bounds
              flag=c(ssb(mp)[,ac(stkYrs)]<c(FLPar(params)["btrig"])[1])
              rtn=res[[1]]
              rtn=qmax(rtn,catch(object)[,ac(refYrs)]*bndTac[1])
              rtn=qmin(rtn,catch(object)[,ac(refYrs)]*bndTac[2])
              
              ##If less than BLim then no TAC constraints
              if (any(flag))
                rtn[,,,,,flag]=res[[1]][,,,,,flag]
              
              object=fwd(object,catch=rtn,sr=eql,residuals=sr_deviances)
              
              chk=rbind(chk,data.frame(iter=factor(seq(dim(rtn)[6]),
                                               levels=seq(seq(dim(rtn)[6])),
                                               labels=seq(dim(rtn)[6])),
                                       res[[2]],catch=c(rtn[,1])))
            }
            cat("==\n")
            
            list(object,chk)})

hcrFn<-function(object,eql,params,
                stkYrs=max(as.numeric(dimnames(stock(object))$year))-2,
                refYrs=max(as.numeric(dimnames(catch(object))$year))-1,
                hcrYrs=max(as.numeric(dimnames(stock(object))$year)),
                maxF  =2,
                nyrs  =3,
                sr_deviances=NULL,
                perfect=FALSE,
                ...) {
  
  ## HCR
  setTarget<-function(params,stk,hcrYrs){
    
    if (length(dim(params))==2){
      flag=c(1,2)[1+as.numeric(c(stk)>c(params["btrig","lower"]))]
    }else{
      params=propagate(FLPar(params),dim(stk)[6])
      flag=rep(1,dim(stk)[6])}
    
    a=FLPar(a=array((params['ftar',flag]-params['fmin',flag])/(params['btrig',flag]-params['blim',flag]),c(1,dim(stk)[6])))
    b=FLPar(b=array( params['ftar',flag]-a*params['btrig',flag],c(1,dim(stk)[6])))    
    rtn=(stk%*%a)  
    rtn=FLCore::sweep(rtn,2:6,b,'+')
    
    fmin=FLQuant(c(params['fmin',flag]),dimnames=list(iter=seq(dim(stk)[6])))
    ftar=FLQuant(c(params['ftar',flag]),dimnames=list(iter=seq(dim(stk)[6])))
    
    for (i in seq(dims(object)$iter)){
      FLCore::iter(rtn,i)[]=max(FLCore::iter(rtn,i),FLCore::iter(fmin,i))
      FLCore::iter(rtn,i)[]=min(FLCore::iter(rtn,i),FLCore::iter(ftar,i))} 
    
    rtn=window(rtn,end=max(hcrYrs))
    if (hcrYrs>(stkYrs+1))
      rtn[,ac((stkYrs+1):(hcrYrs-1))]=fbar(object)[,ac((stkYrs+1):(hcrYrs-1))] #####bug
    for (i in hcrYrs)
      rtn[,ac(i)]=rtn[,dimnames(rtn)$year[1]]
    
    chk=data.frame(hcrYrs=hcrYrs,
                   ssb   =c(stk),
                   f     =c(rtn[,ac(hcrYrs)]))
    list(rtn[,ac(hcrYrs)],chk)}
  
  if (!perfect)
    status=FLCore::apply(ssb(object)[,ac(stkYrs)],6,mean)
  else  
    status=FLCore::apply(ssb(object)[,ac(min(hcrYrs))],6,mean)
  
  res=setTarget(params,status,hcrYrs)
  rtn=res[[1]]
  
  hvt=rtn
  
  ## TACs for target F
  object=window(object, end=max(as.numeric(hcrYrs)))
  
  ## short term projection setting of future vectors
  if (!perfect){
    ## set up vectors for in year projection and hcr
    for (i in slotNames(FLStock())[-c(1,4,7,10,18:20)])
      slot(object,i)[,ac(min(hcrYrs))]=apply(slot(object,i)[,ac(stkYrs-(nyrs)-1)],c(1,6),mean)
  
    if (length(hcrYrs)>1)
      for (i in slotNames(FLStock())[-c(1,4,7,10,18:20)])
        for (j in hcrYrs[-1])
          slot(object,i)[,ac(j)]=slot(object,i)[,ac(hcrYrs[1])]

    if (stkYrs<(hcrYrs-1))
       object=fwd(object,fbar=fbar(object)[,ac(max(stkYrs+1):min(as.numeric(hcrYrs)-1))],sr=eql)
    }
  
  hvt[is.na(hvt)]=6.6666
  
  if (!perfect)
    rtn=catch(fwd(object, fbar=hvt,sr=eql))[,ac(hcrYrs)]
  else  
    rtn=catch(fwd(object, fbar=hvt,sr=eql,residuals=sr_deviances))[,ac(hcrYrs)]
  rtn[]=rep(c(apply(rtn,c(3:6),mean)),each=dim(rtn)[2])
  
  return(list(rtn,res[[2]]))}