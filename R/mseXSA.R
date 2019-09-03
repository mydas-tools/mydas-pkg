utils::globalVariables(c("FLXSA","FLBRP","brp"))

#' @name mseXSA
#' 
#' @title mseXSA 
#' 
#' @description MSE with XSA as the stock assessment in the management procedure 
#' @author Laurence Kell, Sea++
#'     
#' @param om \code{FLStock} object as the operating model
#' @param eq \code{FLBRP}   blah,blah,blah,...
#' @param mp \code{FLStock}  blah,blah,blah,...
#' @param control \code{FLXSA.control}  blah,blah,blah,...
#' @param rf \code{}  blah,blah,blah,...
#' @param ftar \code{numeric}  default is 1.0
#' @param fmin \code{numeric}  default is 0.05 HCR
#' @param bpa \code{numeric}  default is 0.5 HCR
#' @param sigma \code{numeric}  default is 0.3 HCR
#' @param bndTac \code{numeric}  default is c(0.01,100) Bounds on TAC changes
#' @param start \code{numeric}  default is range(om)["maxyear"]-30
#' @param end \code{numeric}  default is range(om)["maxyear"]-interval
#' @param interval  \code{numeric}  default is 3 years over which to run MSE, doesnt work if interval==1, this is a bug
#' @param sr_deviances \code{FLQuant} Stochasticity, either by default or supplied as args rlnoise(dim(om)[6],FLQuant(0,dimnames=list(year=start:(end+interval))),0.3),
#' @param u_deviances   \code{FLQuant} rlnoise(dim(mp)[6],FLQuant(0,dimnames=dimnames(iter(stock.n(om),1))),0.2),
#' @param maxF  \code{numeric} 1.0 Capacity, i.e. F in OM can not be greater than this
#' 
#' @docType methods
#' 
#' @rdname mseXSA
#'
#' @export mseXSA 
#' @import FLBRP
#' 
#' @examples
#' \dontrun{
#' data(pl4)
#' }
mseXSA<-function(  
   #OM as FLStock and FLBRP
   om,eq,
   
   #MP, this is XSA
   mp,control,
   rf="missing",  
   
   #HCR
   ftar=1.0,fmin=0.05,bpa=0.5,sigma=0.3,
   
   #Bounds on TAC changes
   bndTac=c(0.01,100),
   
   #years over which to run MSE, doesnt work if interval==1, this is a bug
   interval=3,start=range(om)["maxyear"]-30,end=range(om)["maxyear"]-interval,
   
   #Stochasticity, either by default or suppliedas args
   sr_deviances=rlnoise(dim(om)[6],FLQuant(0,dimnames=list(year=start:(end+interval))),0.3),
   u_deviances =rlnoise(dim(mp)[6],FLQuant(0,dimnames=dimnames(iter(stock.n(om),1))),0.2),
   
   #Capacity, i.e. F in OM can not be greater than this
   maxF=1.0){ 

  # Check last year so you dont run to the end then crash :-( #####################
  end=min(end,range(om)["maxyear"]-interval)

  # Check iters are consistent across objects #####################################
  ## The OM could have just 1 iter
  if (dims(om)$iter==1 & dims(sr_deviances)$iter>1) 
    om=propagate(om,dims(sr_deviances)$iter)
  
  ## Make sure number of iterations in OM are consistent
  nits=c(om=dims(om)$iter, eq=dims(params(eq))$iter, rsdl=dims(sr_deviances)$iter)
  if (length(unique(nits))>=2 & !(1 %in% nits)) ("Stop, iters not '1 or n' in om")
  if (nits['om']==1) stock(om)=propagate(stock(om),max(nits))

  # Limit on capacity, add to fwd(om,maxF=maxF) so catches dont go crazy ########## 
  maxF=mean(FLQuant(1,dimnames=dimnames(sr_deviances))%*%
              apply(fbar(window(om,end=start)),6,max)*maxF,na.rm=TRUE)

  # Observation Error, setup before looping through years ######################### 
  ## Biology can be different from OM
  pGrp=range(mp)["plusgroup"]
  smp =setPlusGroup(om,pGrp)
  smp=trim(smp,age=range(mp)["min"]:range(mp)["max"])
  
  # Index of abundance ###########################################################
  ## this could be replaced with different series if function modified
  cpue=window(stock.n(smp),end=start-1)[seq(dim(smp)[1]-1)]
  cpue=cpue%*%u_deviances[dimnames(cpue)$age,dimnames(cpue)$year]

  # MP ###########################################################################
  ## no need to add biological parameters at this stage, as set already 
  ## Could use FLBRP to set an average in later versions
  ## rather get rid of stuff that has to be added by OEM and stock assessment fit
  mp=window(mp,end=start-1)

  # Loop round years #############################################################
  cat('\n==')
  for (iYr in seq(start,end,interval)){
    #iYr=start
    cat(iYr,", ",sep="")

  # MP ########################################################################## 
  ## Observation Error, using data from last year back to the last assessment
  smp =trim(setPlusGroup(om[,ac(rev(iYr-seq(interval)))],pGrp),age=range(mp)["min"]:range(mp)["max"])
  
  # CPUE add error ###################################################################
  cpue=window(cpue,end=iYr-1)
  cpue[,ac(iYr-rev(seq(interval)))]=stock.n(smp)[dimnames(cpue)$age,ac(iYr-rev(seq(interval)))]%*%
                  u_deviances[dimnames(cpue)$age,ac(iYr-rev(seq(interval)))]
  idx=FLIndex(index=cpue)
  range(idx)[c("startf","endf")]=c(0.01,0.1)

  #### Management Procedure
  mp=trim(setPlusGroup(om,pGrp),age=range(mp)["min"]:range(mp)["max"],year=dims(om)$minyear:iYr)

  ## fit
  ## Bug with adding range
  xsa=FLXSA(mp,idx,control=control,diag.flag=FALSE)
  range(xsa)[c("min","max","plusgroup")]=range(mp)[c("min","max","plusgroup")]
  mp=mp+xsa
  stock.n(mp)[is.na(stock.n(mp))]=1
  
  ## Update and fill in biological parameters
      
  sr=as.FLSR(window(mp,end=iYr-3),model="geomean")
  sr=fmle(sr,control=list(trace=FALSE),method="L-BFGS-B")
  eql=FLBRP(mp,fbar=FLQuant(seq(0,4,length.out=3))) 

  params(eql)=params(sr)
  mp=fwdWindow(mp,eql,end=iYr+interval)

  ## in year update
  mp=fwd(mp,catch=catch(om)[,ac(iYr+seq(interval))],sr=eql)#,effort_max=maxF)
  #print(plot(FLStocks(MP=mp,OM=window(om,start=dims(mp)$minyear,end=dims(mp)$maxyear))))

  ## HCR
  if ("FLBRP"%in%is(rf))
    hcrPar=icesAR(rf,ftar=ftar,fmin=fmin,bpa=bpa,sigma=sigma)
  else
    hcrPar=icesAR(eql,ftar=ftar,fmin=fmin,bpa=bpa,sigma=sigma)
  
  tac=FLQuant(100,dimnames=list(year=iYr+seq(interval)))

  #save(mp,rf,hcrPar,iYr,interval,bndTac,file="/home/laurence/tmp/hrc.RData")
  tac=hockeyStick(mp,eql,params=hcrPar,
                  stkYrs=iYr-1,
                  refYrs=iYr-1,
                  hcrYrs=iYr+seq(interval),
                  bndTac=bndTac,
                  tac =TRUE)

  tac[is.na(tac)]=1  

  #### Operating Model update
  #try(try(save(om,eq,tac,sr_deviances,maxF,file="/home/laurence/Desktop/tmp/mseXSA3.RData")))

  om =fwd(om,catch=tac,sr=list(model="bevholt",params=params(eq)),residuals=sr_deviances)
          #effort_max=mean(maxF)) 
  }
  
  cat("==\n")
  
  return(om)}


#' icesAR
#' 
#' @title icesAR 
#' 
#' @description # The ICES Advice Rule for stocks where a wide range of stock sizes has been seen 
#' so the stock recruitment relationship can be estimated
#'
#' if SSB>= MSY Btrigger
#'    F=Fmsy
#' else
#'    F=FMSY*SSB/Btrigger  
#'
#' FMSY could be a proxy for FMSY, i.e. F0.1
#'
#' Where
#'   Btrigger is BPA,  i.e. Blim×exp(-1.645×sigma)
#'   Blim is the break point from segmented regression
#'   
#'   # http://ices.dk/sites/pub/Publication%20Reports/Advice/2017/2017/12.04.03.01_Reference_points_for_category_1_and_2.pdf
#'   
#' @name icesAR
#' 
#' @author Laurence Kell, Sea++
#'     
#' @param x \code{FLBRP} with reference points
#' @param ftar default of 1.0
#' @param fmin default of0.05
#' @param bpa default of0.5
#' @param sigma default of 0.3
#' 
#' @docType methods
#' 
#' @rdname icesAR
#'
#' @export  icesAR 
#' 
#' @examples
#' \dontrun{
#' data(pl4)
#' }
icesAR<-function(x,ftar=1.0,fmin=0.05,bpa=0.5,sigma=0.3,refpt="msy"){
  
  hcrParam(
    ftar =refpts(x)[refpt,"harvest"]*ftar,
    btrig=refpts(x)[refpt,    "ssb"]*exp(-1.645*sigma),
    fmin =refpts(x)[refpt,"harvest"]*fmin,
    blim =refpts(x)[refpt,    "ssb"]*bpa)}

