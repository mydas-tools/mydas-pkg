#' @title hcrICES
#'
#' @description
#' Harvest Control Rule, calculates Total Allowable Catch (TAC) based on a hockey stock harvest control rule.
#'
#' @param object an object of class \code{FLStock} 
#' eql \code{FLBRP} with a stock recruitment relationship used for projection
#' sr_deviances  \code{FLQuant} recuitment deviates on the log scale, i.e. multiplicative
#' params \code{} HCR parameters, specifying blim, btrig, ftar and fmin
#' start  \code{} first year for simulation
#' end  \code{}   last year for simukation  
#' interval  \code{} time step, 1 year by default
#' err  \code{} assessment error on SSB for year used for HCR
#' bndTac  \code{} bounds on TAC, by default these are turned off, for 20% constraint set to c(0.8,1.2)
#' 
#' @return returns a \code{list} with \code{FLStock} and \code{FLPar} objects for the stock and HCR 
#'
#' @export
#' @rdname hcrICES
#'
#' @examples
#' \dontrun{
#' }

setGeneric('hcrICES', function(object,eql,...) standardGeneric('hcrICES'))

setMethod('hcrICES', signature(object="FLStock",eql='FLBRP'),
    function(object,eql,sr_deviances,params,
             start   =max(dimnames(object)$year)-10,
             end     =start+10,
             interval=1,
             err=NULL,
             bndTac=c(0,Inf)){
  
  bias<-function(object,eql,err,iYr,lag=1){
    stock.n(object)=stock.n(object)%*%err[,ac(iYr)]
  
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
    refpts(eql)=propagate(refpts(eql),100)
    
    if (!is.null(err))
      mp=bias(object,eql,err,iYr,lag=hcrYrs-stkYrs-1)
    else
      mp=object
   
    res=hcrFn(mp,eql,params,
              stkYrs,
              refYrs,
              hcrYrs)

    ##TAC Bounds
    flag=c(ssb(mp)[,ac(stkYrs)]<c(FLPar(params)["btrig"])[1])
    rtn=res[[1]]
    rtn=qmax(rtn,catch(object)[,ac(refYrs)]*bndTac[1])
    rtn=qmin(rtn,catch(object)[,ac(refYrs)]*bndTac[2])
 
    ##If less than BLim then no TAC constraints
    if (any(flag))
      rtn[,,,,,flag]=res[[1]][,,,,,flag]
    
    object=fwd(object,catch=rtn,sr=eql,residuals=sr_deviances)
    
    chk=rbind(chk,res[[2]])
  }
  cat("==\n")
  
  list(object,chk)})

hcrFn<-function(object,eql,params,
                stkYrs=max(as.numeric(dimnames(stock(object))$year))-2,
                refYrs=max(as.numeric(dimnames(catch(object))$year))-1,
                hcrYrs=max(as.numeric(dimnames(stock(object))$year)),
                maxF  =2,
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

  status=FLCore::apply(ssb(object)[,ac(stkYrs)],6,mean)
  
  res=setTarget(params,status,hcrYrs)
  rtn=res[[1]]

  hvt=rtn
 
  ## TACs for target F
  object=window(object, end=max(as.numeric(hcrYrs)))
  object[,ac(max(as.numeric(hcrYrs)))]=object[,ac(max(as.numeric(hcrYrs))-1)]

  object=fwd(object,fbar=fbar(object)[,ac(min(as.numeric(hcrYrs)-1))],sr=eql)

  hvt[is.na(hvt)]=6.6666

  rtn=catch(fwd(object, fbar=hvt,sr=eql))[,ac(hcrYrs)]
  rtn[]=rep(c(apply(rtn,c(3:6),mean)),each=dim(rtn)[2])

  return(list(rtn,res[[2]]))}