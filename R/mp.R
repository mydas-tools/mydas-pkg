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


#https://github.com/merrillrudd/LIME
#https://github.com/AdrianHordyk/LBSPR
#https://github.com/datalimited/datalimited
#https://github.com/datalimited/global-status-estimates

hcrSBT1<-function(cpue,tac,k1=2.0,k2=3.0,gamma=1,nyrs=5,lag=1,interval=3){
  
  dat=as.data.frame(cpue[,ac(-((nyrs-1):0)+dims(cpue)$maxyear)])
  lambda=as.FLQuant(ddply(dat, .(iter), with,  data.frame(data=coefficients(lm(data~year))[2])))
  
  flag  =lambda<0
  lambda=abs(lambda)
  
  #res=1+ifelse(flag,-k1,k2)*lambda^ifelse(flag,gamma,1)
  res=1+ifelse(flag,-k1,k2)*exp(log(lambda)*ifelse(flag,gamma,1))
  
  res=res%*%tac
  
  dmns=dimnames(tac)
  dmns$year=as.integer(dmns$year)+lag+seq(interval)-1
  dmns$iter=dimnames(cpue)$iter
  
  res=FLQuant(rep(c(res),each=interval),dimnames=dmns)
  
  return(res)}
