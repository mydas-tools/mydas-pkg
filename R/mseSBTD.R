#' @import FLCore 
#' 
#' @export mseSBTD
#' @export hcrSBTD

mseSBTD<-function(
  #OM as FLStock and FLBRP
  om,eq,
  
  #MP
  control="missing",
  
  srDev,
  uDev,
  
  #years over which to run MSE, doesnt work if interval==1, this is a bug
  interval=1,start=range(om)["maxyear"]-30,end=range(om)["maxyear"]-interval,
  
  #Capacity, i.e. F in OM can not be greater than this
  nyrs  =5,
  cpueFn=ssb,
  lag   =1,
  maxF  =10){

  ##So you dont run to the end then crash
  end=min(end,range(om)["maxyear"]-interval)
  
  ## Make sure number of iterations are consistent
  nits=c(om=dims(om)$iter, eq=dims(params(eq))$iter, rsdl=dims(srDev)$iter)
  if (length(unique(nits))>=2 & !(1 %in% nits)) ("Stop, iters not '1 or n' in om")
  if (nits['om']==1) stock(om)=propagate(stock(om),max(nits))

  ## Observation Error (OEM) setup
  cpue=window(cpueFn(om),end=start)
  cpue=cpue%*%uDev[,dimnames(cpue)$year]

  ## Loop round years
  cat('\n==')
  for (iYr in seq(start,end,interval)){
    cat(iYr,", ",sep="")

    ## Observation Error, using data from last year back to the last assessment
    ## CPUE
    cpue=window(cpue,end=iYr-lag)
    uYrs=rev(seq(nyrs))-1+lag
    cpue[,ac(iYr-uYrs)]=(cpueFn(om))[,ac(iYr-uYrs)]%*%uDev[,ac(iYr-uYrs)]

    #### Management Procedure
    ##Constant catch
    #tac=hcrConstantCatch(iYr+1,catch=catch(om)[,ac(iYr-(2:1))]) 
    tac=hcrSBTD(iYr+seq(interval),
                control=control,
                cpue[,ac(iYr-uYrs)],
                catch(om)[,ac(iYr-rev(seq(interval))+1)])
   
    #### Operating Model update
    om=fwd(om,catch=tac,sr=eq,residual=srDev,maxF=maxF)
    }
  cat('==\n')
  
  return(om)}
