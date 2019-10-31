# Length-based (comparing outputs for Fishing mortality estimated by year) 
## 1. Mean length at equilibruim (Beverton and Holt 1975)
## 2. LBSPR: Life-history and length distribution based (Hordyk et al. 2014)
## 3. LIME: Accounting for variable recruitment and fishing mortality in length-based stock assessments for data-limited fisheries (Rudd and Thorson, under review)
## 4. Utilizing B-H invariants and catches as a function of size (Kokkalis et al. 2015)

# Catch-based (comparing outputs for Biomass or B/Bmsy estimated by year):
##  5. DB-SRA: Depletion-based Stock reduction analysis (Dick and MacCall 2011)
##  6. Catch-MSY (Martel and Froese 2013)
##  7. COM-SIR (Vasconcellos and Cochrane 2005) - Bayesian approach
##  8. SS-COM: State-space-Catch-Only-method (Thorson et al. 2013)
##  9. Modified panel regression (Costello et al. 2012)
## 10. Ensemble method: (Rosenberg et al. 2017) 


mseEMP<-function(
  #OM as FLStock and FLBRP
  om,eql,
  
  #MP,
  control=c(k1=3.0,k2=1.5,gamma=1,nyrs=5,lag=1,interval=3),
  
  #years over which to run MSE, doesnt work if interval==1, this is a bug
  start=range(om)["maxyear"]-30,interval=3,end=range(om)["maxyear"]-interval,
  
  #Stochasticity, either by default or suppliedas args
  srDev=rlnoise(dim(om)[6],FLQuant(0,dimnames=list(year=start:end)),0.3), 
  uDev =rlnoise(dim(mp)[6],FLQuant(0,dimnames=dimnames(iter(stock.n(om),1))),0.2),
  
  #Capacity, i.e. F in OM can not be greater than this
  maxF=1.5){ 
  
  ##So you dont run to the end then crash
  end=min(end,range(om)["maxyear"]-interval)
  
  ## Make sure number of iterations are consistent
  nits=c(om=dims(om)$iter, eql=dims(params(eql))$iter, rsdl=dims(srDev)$iter)
  if (length(unique(nits))>=2 & !(1 %in% nits)) ("Stop, iters not '1 or n' in om")
  if (nits['om']==1) stock(om)=propagate(stock(om),max(nits))
  
  ## Limit on capacity, add to fwd(om) if you want
  maxF=FLQuant(1,dimnames=dimnames(srDev))%*%apply(fbar(window(om,end=start)),6,max)*maxF
  
  ## Observation Error (OEM) setup 
  pGrp=range(om)["plusgroup"]
  
  cpue=window(stock.n(om),end=start)[dimnames(uDev)$age]
  cpue=cpue%*%uDev[,dimnames(cpue)$year]
  
  ## Loop round years
  cat('\n==')
  for (iYr in seq(start,end,interval)){
    cat(iYr,", ",sep="")
    
    ## Observation Error, using data from last year back to the last assessment
    ## CPUE
    cpue=window(cpue,end=iYr-1)
    cpue[,ac(iYr-(interval:1))]=stock.n(om)[dimnames(cpue)$age,ac(iYr-(interval:1))]%*%uDev[,ac(iYr-(interval:1))]
    #### Management Procedure
    
    u=window(apply(cpue,c(2,6),mean),end=iYr-1)
    tac=hcrSBT1(u,catch(om)[,ac(iYr-1)])
    
    #### Operating Model update
    om =fwd(om,catch=tac,sr=eql,sr.residuals=srDev,maxF=mean(maxF)) 
  }
  cat('==\n')
  
  return(om)}