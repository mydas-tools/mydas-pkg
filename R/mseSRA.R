
mseSRA<-function(
  #OM
  om,eq,
  
  #MP
  mp,
  #http://ices.dk/sites/pub/Publication%20Reports/Advice/2017/2017/12.04.03.01_Reference_points_for_category_1_and_2.pdf
  ftar=1.0,btrig=0.5,fmin=0.05,blim=0.3,        
  
  #years over which to run MSE
  interval=3,start=range(om)["maxyear"]-30,end=range(om)["maxyear"]-interval,
  
  #Stochasticity
  srDev, #=rlnorm(dim(om)[6],FLQuant(0,dimnames=list(year=start:end)),0.3),

  #Capacity, i.e. F in OM can not be greater than this
  maxF=1.5){ 
  
  ## Get number of iterations in OM
  nits=c(om=dims(om)$iter, eq=dims(params(eq))$iter, rsdl=dims(srDev)$iter)
  if (length(unique(nits))>=2 & !(1 %in% nits)) ("Stop, iters not '1 or n' in om")
  if (nits['om']==1) stock(om)=propagate(stock(om),max(nits))
  
  mp=window(mp,end=start)
  
  ## Cut in capacity
  maxF=median(apply(fbar(window(om,end=start)),6,max)*maxF)
  
  ## Loop round years
  for (iYr in seq(start,end-interval,interval)){
    cat(iYr,", ",sep="")
    
    ##OEM
    mp=window(mp,end=iYr-1)
    #bug in window
    catch(mp)[,ac(rev(iYr-seq(interval+1)))]=catch(om)[,ac(rev(iYr-seq(interval+1)))]
    
    ##MP
    mp=sra(mp,ssb(om)[,ac(iYr)]%/%refpts(eq)["virgin","ssb"])
    mp=window(mp,end=iYr)
    
    #bug in window
    catch(mp)[,ac(rev(iYr-seq(interval+1)))]=catch(om)[,ac(rev(iYr-seq(interval+1)))]
    catch(mp)[,ac(iYr)]=catch(om)[,ac(iYr)]
    
    #try(save(mp,om,file="/home/laurence/Desktop/test1.RData"))
    mp=fwd(mp,catch=catch(mp)[,ac(iYr)])
    
    ## HCR
    par=hcrParam(ftar =ftar*refpts( mp)["fmsy"],
                 btrig=btrig*refpts(mp)["bmsy"],
                 fmin =fmin*refpts( mp)["fmsy"],
                 blim =blim*refpts( mp)["bmsy"])
    
    #try(save(mp,par,file="/home/laurence/Desktop/test2.RData"))
    tac=hcr(mp,refs=par,hcrYrs=iYr+seq(interval),tac=TRUE)
    tac[is.na(tac)]=1
    
    #### Operating Model Projectionfor TAC
    #try(save(om,tac,sr,eq,srDev,maxF,file="/home/laurence/Desktop/test3.RData"))
    om =fwd(om,catch=tac,sr=eq,residuals=srDev,effort_max=maxF)  
    print(plot(as(list("MP"=                     window(mp,end=iYr),
                       "OM"=as(window(om,end=iYr+interval),"biodyn")),"biodyns")))
  }
  
  return(om)}
