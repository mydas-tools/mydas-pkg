#' mse
#' 
#' @title mseSBTP 
#' 
#' @description 
#' @author Laurence Kell, Sea++
#'  
#' @name mseSBTP
#' 
#' @export mseSBTP
#' 
#' @param control blah,blah,blah,...
#' 
#' @param om \code{FLStock} object as the operating model
#' @param eq blah,blah,blah,...
#' @param sr_deviates blah,blah,blah,...
#' @param u_deviates blah,blah,blah,...
#' @param interval blah,blah,blah,...
#' @param start \code{numeric}  default is range(om)["maxyear"]-30
#' @param end \code{numeric}  default is range(om)["maxyear"]-interval
#' @param yrU blah,blah,blah,...
#' @param yrCatch blah,blah,blah,...
#' @param refU blah,blah,blah,...
#' @param refCatch blah,blah,blah,...
#' @param maxF blah,blah,blah,...
#' @param ... any additional arguments
#' 
#' @docType methods
#' 
#' @export mseSBTP
#' @rdname mseSBTP
#' 
#' @examples
#' \dontrun{
#' data(pl4)
#' }
#' 
mseSBTP<-function(
  #OM as FLStock and FLBRP
  om,eq,
  
  sr_deviates,
  u_deviates,
  
  #years over which to run MSE, doesnt work if interval==1, this is a bug
  interval=1,start=range(om)["maxyear"]-30,end=range(om)["maxyear"]-interval,
  
  control=c(k1=0.25,k2=0.25),
  
  #how many years to go back from current year
  yrU     =1:2,   #to estimate current index level    
  yrCatch =seq(interval)-1,     #to estimate last catch
  refU    ="missing", #years for reference Index
  refCatch="missing", #years for reference catch
  
  #Capacity, i.e. F in OM can not be greater than this
  maxF=1.5){

  ##So you dont run to the end then crash
  end=min(end,range(om)["maxyear"]-interval)
 
  ## Make sure number of iterations are consistent
  nits=c(om=dims(om)$iter, eq=dims(params(eq))$iter, rsdl=dims(sr_deviates)$iter)
  if (length(unique(nits))>=2 & !(1 %in% nits)) ("Stop, iters not '1 or n' in om")
  if (nits['om']==1) stock(om)=propagate(stock(om),max(nits))
  
  ## Limit on capacity, add to fwd(om) if you want
  maxF=median(FLQuant(1,dimnames=dimnames(sr_deviates))%*%apply(fbar(window(om,end=start)),6,max)*maxF)
  
  ## Observation Error (OEM) setup
  cpue=window(ssb(om),end=start)
  
  cpue=cpue%*%u_deviates[,dimnames(cpue)$year]
  
  ## Loop round years
  cat('==')
  for (iYr in seq(start,end,interval)){
    cat(iYr,", ",sep="")
    
    ## Observation Error, using data from last year back to the last assessment
    ## CPUE
    cpue=window(cpue,end=iYr-1)
    cpue[,ac(iYr-(interval:1))]=ssb(om)[,ac(iYr-(interval:1))]%*%u_deviates[,ac(iYr-(interval:1))]

    #### Management Procedure
    ##Constant catch
    #tac=hcrConstantCatch(iYr+1,catch=catch(om)[,ac(iYr-(2:1))]) 
    tac=hcrSBTP(yrs   =iYr+seq(interval),
                control=control,
                ref     =apply(cpue[,     ac(refU)],       6,mean),
                target  =apply(catch(om)[,ac(refCatch)],   6,mean),
                cpue    =apply(cpue[,     ac(iYr-yrU)],    6,mean),
                catch   =apply(catch(om)[,ac(iYr-yrCatch)],6,mean))

    #### Operating Model update
    om =fwd(om,catch=tac,sr=eq,residual=sr_deviates,effort_max=mean(maxF))
    #print(plot(window(om,end=iYr+interval)))
  }
  cat('==\n')
  
  return(om)}
