#' @title hcr
#'
#' @description
#' Harvest Control Rule, calculates F, or Total Allowable Catch (TAC) based on a hockey stock harvest control rule.
#'
#' @param object an object of class \code{biodyn} or
#' @param ... other parameters, i.e.
#' refs \code{FLPar or FLBRP} object with reference points, can be missing if refpts are part of \code{object}
#' params \code{FLPar} object with hockey stick HCR parameters, see hcrParam
#' yr numeric vector with years used to values in HCR
#' byr numeric vector with years used for bounds
#' hyr numeric vector with years to use in projection
#' tac \code{logical} should return value be TAC rather than F?
#' bndF \code{vector} with bounds (i.e.min and max values) on iter-annual variability on  F
#' bndTac \code{vector} with bounds (i.e. min and max values) on iter-annual variability on TAC
#' stkYrs numeric vector with years for calculating stock status, max(as.numeric(dimnames(stock(object))$year)),
#' refYrs numeric vector with years for calculating reference points, max(as.numeric(dimnames(catch(object))$year)),
#' hcrYrs numeric vector with years for setting management, max(as.numeric(dimnames(stock(object))$year)),
#' tacMn \code{logical} should the TACs be the average of the values returned?
#' maxF  =2 numeric vector specifying maximum relative F
#' @aliases hcr,biodyn-method
#'
#' @return \code{FLPar} object with value(s) for F or TAC if tac==TRUE
#'
#' @export
#' @rdname hcr
#'
#' @examples
#' \dontrun{
#' bd   =sim()
#'
#' bd=window(bd,end=29)
#' for (i in seq(29,49,1))
#' bd=fwd(bd,harvest=hcr(bd,yr=i,yr=i+1)$hvt)
#' }

setGeneric('hcr', function(object,refs,...) standardGeneric('hcr'))

setMethod('hcr', signature(object="FLStock",refs='FLBRP'),
          function(object,refs,
          params=hcrParam(ftar =0.70*fmsy(refs),
                          btrig=0.80*bmsy(refs),
                          fmin =0.01*fmsy(refs),
                          blim =0.40*bmsy(refs)),
          stkYrs=max(as.numeric(dimnames(stock(object))$year)),
          refYrs=max(as.numeric(dimnames(catch(object))$year)),
          hcrYrs=max(as.numeric(dimnames(stock(object))$year)),
          tac   =TRUE,
          tacMn =TRUE,
          bndF  =NULL, #c(1,Inf),
          bndTac=NULL,
          maxF  =2,
          ...)
  hcrFn(object,refs,
                params,stkYrs,refYrs,hcrYrs,tac,bndF,bndTac,maxF,...))

setMethod('hcr', signature(object="biodyn",refs='FLPar'),
  function(object,refs=hcrParam(ftar =0.70*mpb:::fmsy(refs),
                                btrig=0.80*mpb:::bmsy(refs),
                                fmin =0.01*mpb:::fmsy(refs),
                                blim =0.40*mpb:::bmsy(refs)),
           params=refs,
           stkYrs=max(as.numeric(dimnames(stock(object))$year)),
           refYrs=max(as.numeric(dimnames(catch(object))$year)),
           hcrYrs=max(as.numeric(dimnames(stock(object))$year)),
           tac   =TRUE,
           bndF  =NULL, #c(1,Inf),
           bndTac=NULL,
           maxF  =2,
           ...)
  hcrFn(object,refs,
        params,stkYrs,refYrs,hcrYrs,tac,bndF,bndTac,maxF,...))

hcrFn<-function(object,refs, 
               params=hcrParam(ftar =0.70*refpts(object)['fmsy'],
                               btrig=0.80*refpts(object)['bmsy'],
                               fmin =0.01*refpts(object)['fmsy'],
                               blim =0.40*refpts(object)['bmsy']),
               stkYrs=max(as.numeric(dimnames(stock(object))$year)),
               refYrs=max(as.numeric(dimnames(catch(object))$year)),
               hcrYrs=max(as.numeric(dimnames(stock(object))$year)),
               tac   =TRUE,
               bndF  =NULL, #c(1,Inf),
               bndTac=NULL, 
               maxF  =2,
               ...) {
  
  ## HCR
  dimnames(params)$params=tolower(dimnames(params)$params)
  params=as(params,'FLQuant')  
  params["blim"]=c(qmax(params["blim"],1))
  
  #if (blim>=btrig) stop('btrig must be greater than blim')
  a=(params['ftar']-params['fmin'])/(params['btrig']-params['blim'])
  b=params['ftar']-a*params['btrig']

  ## Calc F
  # bug 
  #val=(SSB%*%a) %+% b
  # bug stock for biomass
  stk=FLCore::apply(stock(object)[,ac(stkYrs)],6,mean)

  rtn=(stk%*%a)  
  rtn=FLCore::sweep(rtn,2:6,b,'+')

  fmin=as(params['fmin'],'FLQuant')
  ftar=as(params['ftar'],'FLQuant')
  for (i in seq(dims(object)$iter)){
    FLCore::iter(rtn,i)[]=max(FLCore::iter(rtn,i),FLCore::iter(fmin,i))
    FLCore::iter(rtn,i)[]=min(FLCore::iter(rtn,i),FLCore::iter(ftar,i))} 
  
  rtn=window(rtn,end=max(hcrYrs))
  #dimnames(rtn)$year=min(hcrYrs)  

  rtn=window(rtn,end=max(hcrYrs))
  rtn[,ac(hcrYrs)]=rtn[,dimnames(rtn)$year[1]]

  ### Bounds ##################################################################################
  ## F
  if (!is.null(bndF)){  
    
    ref=FLCore::apply(harvest(object)[,ac(refYrs-1)],6,mean)
    
    rtn[,ac(min(hcrYrs))]=qmax(rtn[,ac(min(hcrYrs))],ref*bndF[1])
    rtn[,ac(min(hcrYrs))]=qmin(rtn[,ac(min(hcrYrs))],ref*bndF[2])
    
    for (i in hcrYrs[-1]){
        if (iaF){
          rtn[,ac(i)]=qmax(rtn[,ac(i)],rtn[,ac(i-1)]*bndF[1])
          rtn[,ac(i)]=qmin(rtn[,ac(i)],rtn[,ac(i-1)]*bndF[2])
        }else{
          rtn[,ac(i)]=rtn[,ac(i-1)]}
        
        if (!is.null(maxF)) rtn=qmin(rtn,maxF)}
    }

  hvt=rtn

  ## TAC
  if (!tac)
    return(hvt)
  else{
    ## TACs for target F

    if ("FLBRP"%in%is(object))
      object=fwdWindow(object, refs, end=max(as.numeric(hcrYrs)))
    else
      object=window(object, end=max(as.numeric(hcrYrs)))
    object[,ac(max(as.numeric(hcrYrs)))]=object[,ac(max(as.numeric(hcrYrs))-1)]

    #object[,ac(max(as.numeric(hcrYrs)))]=object[,ac(max(as.numeric(hcrYrs))-1)]
    
    if ("FLStock"%in%is(object))      
       object=fwd(object,fbar=fbar(object)[,ac(min(as.numeric(hcrYrs)-1))],sr=refs)
     else
       object=fwd(object,harvest=harvest(object)[,ac(min(as.numeric(hcrYrs)-1))])

    hvt[is.na(hvt)]=6.6666
    #try(save(object,hvt,refs,hcrYrs,file="/home/laurence/Desktop/tmp/mseXSA4.RData"))
    if ("FLStock"%in%is(object))      
       rtn=catch(fwd(object, fbar=hvt,sr=refs))[,ac(hcrYrs)]
    else
       rtn=catch(fwd(object, harvest=hvt))[,ac(hcrYrs)]
    rtn[]=rep(c(apply(rtn,c(3:6),mean)),each=dim(rtn)[2])

    ## Bounds
    if (!is.null(bndTac)){
      ## Reference TAC for bounds
      ref=c(apply(catch(object)[,ac(refYrs)],6,mean))
      ref=FLQuant(rep(ref,each=dim(rtn)[2]),dimnames=dimnames(rtn))

      rtn=qmax(rtn,ref*bndTac[1])
      rtn=qmin(rtn,ref*bndTac[2])
    }
  }

  return(rtn)}
