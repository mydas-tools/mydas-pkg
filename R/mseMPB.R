utils::globalVariables(c("objFn","setParams<-","setControl<-","fit"))

#' #http://ices.dk/sites/pub/Publication%20Reports/Advice/2017/2017/12.04.03.01_Reference_points_for_category_1_and_2.pdf

#' mse
#' 
#' @title mseMPB 
#' 
#' @description 
#' @author Laurence Kell, Sea++
#'  
#' @name mseMPB
#' 
#' @aliases mseMPB mseMPB-method mseMPB,FLStock,FLBRP-method
#' 
#' @param om \code{FLStock} object as the operating model
#' @param eq  blah,blah,blah,...
#' @param mp  blah,blah,blah,...
#' @param ftar 1.0  blah,blah,blah,...
#' @param btrig 0.5  blah,blah,blah,...
#' @param fmin 0.05  blah,blah,blah,...
#' @param blim 0.3  blah,blah,blah,...      
#' @param sr_deviates rlnorm(dim(om)[6],FLQuant(0,dimnames=list(year=start:end)),0.3),
#' @param u_deviates  =rlnorm(dim(om)[6],FLQuant(0,dimnames=list(year=start:end)),0.2),
#' @param interval 3  blah,blah,blah,...
#' @param start \code{numeric}  default is range(om)["maxyear"]-30
#' @param end \code{numeric}  default is range(om)["maxyear"]-interval
#' @param maxF 1.5 
#'
#' @docType methods
#' 
#' @rdname mseMPB
#' 
#' @examples
#' \dontrun{
#' data(pl4)
#' }
#' 
#' @export mseMPB
#' 
mseMPB<-function(
  om,eq,
  
  mp,
  #http://ices.dk/sites/pub/Publication%20Reports/Advice/2017/2017/12.04.03.01_Reference_points_for_category_1_and_2.pdf
  ftar=1.0,btrig=0.5,fmin=0.05,blim=0.3,        
  
  #years over which to run MSE
  interval=3,start=range(om)["maxyear"]-30,end=range(om)["maxyear"]-interval,
  
  #Stochasticity
  sr_deviates, #=rlnorm(dim(om)[6],FLQuant(0,dimnames=list(year=start:end)),0.3),
  u_deviates,  #=rlnorm(dim(om)[6],FLQuant(0,dimnames=list(year=start:end)),0.2),
  
  #Capacity, i.e. F in OM can not be greater than this
  maxF=1.5){ 
  
  ## Get number of iterations in OM
  nits=c(om=dims(om)$iter, eq=dims(params(eq))$iter, rsdl=dims(sr_deviates)$iter)
  if (length(unique(nits))>=2 & !(1 %in% nits)) ("Stop, iters not '1 or n' in om")
  if (nits['om']==1) stock(om)=propagate(stock(om),max(nits))
  
  mp=window(mp,end=start)
  
  ## Cut in capacity
  maxF=median(apply(fbar(window(om,end=start)),6,max)*maxF)
  
  #### Observation Error (OEM) setup 
  cpue=window(catch(om)%/%fbar(om),end=start)
  cpue=cpue%*%u_deviates[,dimnames(cpue)$year]
  
  ## Loop round years
  for (iYr in seq(start,end-interval,interval)){
    cat(iYr,", ",sep="")
    
    ##OEM
    mp=window(mp,end=iYr-1)
    #bug in window
    catch(mp)[,ac(rev(iYr-seq(interval+1)))]=catch(om)[,ac(rev(iYr-seq(interval+1)))]
    
    cpue=window(cpue,end=iYr-1)
    cpue[,ac(iYr-(interval:1))]=catch(om)[, ac(iYr-(interval:1))]%*%
      u_deviates[      ,ac(iYr-(interval:1))]%/%
      fbar( om)[ ,ac(iYr-(interval:1))]
    #### Management Procedure
    mp@indices=FLQuants("1"=apply(window(cpue,start=range(mp)["minyear"]),c(2,6),sum))
    
    ##MP
    mp@control["sigma1","val"]=cv(u_deviates)
    mp=fit(mp)
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
    #try(save(om,tac,sr,eq,sr_deviates,maxF,file="/home/laurence/Desktop/test3.RData"))
    om =fwd(om,catch=tac,sr=eq,residuals=sr_deviates,effort_max=maxF)  
    #print(plot(as(list("MP"=                     window(mp,end=iYr),
    #                   "OM"=as(window(om,end=iYr+interval),"biodyn")),"biodyns")))
  }
  
  return(om)}

mpTrace<-function(mp,par,tac){
  
  rtn1=cbind(model.frame(params(mp)),
             model.frame(refpts(mp)),
             model.frame(par),
             model.frame(objFn(mp)))
  
  rtn2=model.frame(tac[["tac"]])
  rtn3=model.frame(tac[["stk"]])
  }
         
icesBD<-function(bd,fmsy=1.0,btrig=0.5,blim=0.3, fmin=0.05){
  # http://ices.dk/sites/pub/Publication%20Reports/Advice/2017/2017/12.04.03.01_Reference_points_for_category_1_and_2.pdf
 
  hcrParam(
    ftar =fmsy(bd),
    btrig=bmsy(bd)*btrig,
    blim =bmsy(bd)*blim,
    fmin =fmsy(bd)*fmin)}

mseMPB2<-function(
  #OM
  om,eq,
  
  #MP
  mp,
  #http://ices.dk/sites/pub/Publication%20Reports/Advice/2017/2017/12.04.03.01_Reference_points_for_category_1_and_2.pdf
  ftar=1.0,btrig=0.5,fmin=0.05,blim=0.3,        
  
  #years over which to run MSE
  interval=3,start=range(om)["maxyear"]-30,end=range(om)["maxyear"]-interval,
  
  #Stochasticity
  sr_deviates,   #=rlnorm(dim(om)[6],FLQuant(0,dimnames=list(                     year=start:end)), sr_deviates),
  u_deviates,    #=rlnorm(dim(om)[6],FLQuant(0,dimnames=list(                     year=start:end)),  u_deviates),
  sel_deviates,  #=rlnorm(dim(om)[6],FLQuant(0,dimnames=list(age=dimnames(om)$age,year=start:end)),Sel_deviates),
  
  #Capacity, i.e. F in OM can not be greater than this
  maxF=1.5){ 
  
  if ("FLQuant"%in%is(  u_deviates)) u_deviates  =FLQuants(u_deviates)
  if ("FLQuant"%in%is(sel_deviates)) sel_deviates=FLQuants(sel_deviates)
  
  ## Get number of iterations in OM
  nits=c(om=dims(om)$iter, eq=dims(params(eq))$iter, rsdl=dims(sr_deviates)$iter)
  if (length(unique(nits))>=2 & !(1 %in% nits)) ("Stop, iters not '1 or n' in om")
  if (nits['om']==1) stock(om)=propagate(stock(om),max(nits))
  
  mp=window(mp,end=start)
  
  ## Cut in capacity
  maxF=median(apply(fbar(window(om,end=start)),6,max)*maxF)

  #### Observation Error (OEM) setup
  nU=max(length(u_deviates),length(sel_deviates))
  cpue=FLQuants()
  for (iU in seq(nU)){
    yrs=dimnames(window(sel_deviates[[iU]],end=min(start,dims(sel_deviates[[iU]])$maxyear)))$year
    yrs=yrs[yrs%in%dimnames(u_deviates[[iU]])$year]
    yrs=yrs[yrs%in%dimnames(m(om))$year]
    u         =catch.wt(om)[,yrs]%*%catch.n(om)[,yrs]%/%fbar(om)[,yrs]%*%sel_deviates[[iU]][,yrs]
    cpue[[iU]]=apply(u,2:6,sum)%*%u_deviates[[iU]][,yrs]}
    
  ## Loop round years
  cat("==")
  for (iYr in seq(start,end-interval,interval)){
    cat(iYr,", ",sep="")

    ##OEM
    mp=window(mp,end=iYr-1)
    #bug in window
    catch(mp)[,ac(rev(iYr-seq(interval+1)))]=catch(om)[,ac(rev(iYr-seq(interval+1)))]

    for (iU in seq(nU)){
      yrs=ac(iYr-(interval:1))
      yrs=yrs[yrs%in%dimnames(u_deviates[[iU]])$year]
      yrs=yrs[yrs%in%dimnames(sel_deviates[[iU]])$year]
      if (length(yrs)>0){
        cpue[[iU]]=window(cpue[[iU]],end=iYr-1)
        u         =(catch.wt(om)[,yrs]%*%catch.n(om)[,yrs]%/%fbar(om)[,yrs])%*%sel_deviates[[iU]][,yrs]
        cpue[[iU]][,yrs]=apply(u,2:6,sum)%*%u_deviates[[iU]][,yrs]}}

    #### Management Procedure
    params(    mp)=params( mp)[c("r","k","p","b0")]
    control(   mp)=control(mp)[c("r","k","p","b0")]
    setParams( mp)=cpue
    setControl(mp)=params(mp)
    
    ##MP
    mp=fit(mp,cpue)
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
    #try(save(om,tac,sr,eq,sr_deviates,maxF,file="/home/laurence/Desktop/test3.RData"))
    om =fwd(om,catch=tac,sr=eq,residuals=sr_deviates)#,effort_max=maxF)  
    #print(plot(as(list("MP"=                     window(mp,end=iYr),
    #                   "OM"=as(window(om,end=iYr+interval),"biodyn")),"biodyns")))
  }
  cat("==\n")
  
  return(om)}

#### Old and not used anymore
## Added OEM needs checking
## To do add TAC bounds  
mseAlbn<-function(
  #OM
  
  om,eql,

  #OEM
  saa,rng,
  
  #MP
  mp,
  ftar =0.7,btrig=0.8,fmin=0.1,blim=0.4,        
  #years over which to run MSE
  start=range(om)["maxyear"]-30,interval=3,end=range(om)["maxyear"]-interval,
  
  #Stochasticity
  sr_deviates=rlnorm(dim(om)[6],FLQuant(0,dimnames=list(year=start:end)),0.3), 
  u_deviates =rlnorm(dim(om)[6],FLQuant(0,dimnames=dimnames(iter(stock(om),1))),0.3),
  
  #Capacity, i.e. F in OM can not be greater than this
  maxF=1.5){ 
  
  ## Get number of iterations in OM
  nits=c(om=dims(om)$iter, eql=dims(params(eql))$iter, rsdl=dims(sr_deviates)$iter)
  if (length(unique(nits))>=2 & !(1 %in% nits)) ("Stop, iters not '1 or n' in om")
  if (nits['om']==1) stock(om)=propagate(stock(om),max(nits))
  
  #### Observation Error (OEM)
  cpue=FLQuants(mlply(seq(length(saa)), function(i){
    res=window(apply(oem(om,saa[[i]]),2,sum),
              start=rng[i,"minyear"],end=min(yrRng[i,"maxyear"],start))
    res%*%u_deviates[[i]][,dimnames(res)$year]
    }))
  
  ## SA
  params=params(mp)
  mp=as(window(om,end=start-1),"biodyn")
  mp@indices=cpue
  mp=fwd(mp,catch=catch(mp))
  
  setParams(mp)=mp@indices
  params(mp)[dimnames(params)$params]=params
  setControl(mp)=params(mp)
  control(mp)[substr(dimnames(control(mp))$params,1,1)=="q",1]=1
  
  ## Cut in capacity
  maxF=FLQuant(1,dimnames=dimnames(sr_deviates))%*%apply(fbar(window(om,end=start)),6,max)*maxF

  ## Loop round years
  for (iYr in seq(start,end-interval,interval)){
    cat(iYr,", ",sep="")
    
    mp=window(mp,end=iYr-1)
    
    ##OEM
    catch(mp)[,ac(iYr-interval:1)]=catch(om)[,ac(iYr-interval:1)]
    
    cpue=FLQuants(llply(cpue,window,end=iYr-1))
    for (i in seq(length(cpue))){
      res=window(apply(oem(om,saa[[i]]),2,sum),start=iYr-interval,end=iYr-1)
      cpue[[i]][,ac(iYr-interval:1)]=res%*%u_deviates[[i]][,ac(iYr-interval:1)]}
    
    ##MP
    params(mp)["r"]=0.3
    mp=fwd(mp,catch=catch(mp))
    mp@indices=cpue
    setParams(mp)=cpue
    setControl(mp)=params(mp)
    mp=fit(mp)
    
    ## HCR
    par=hcrParam(ftar =ftar*refpts(mp)["fmsy"],
                 btrig=btrig*refpts(mp)["bmsy"],
                 fmin =fmin*refpts(mp)["fmsy"],
                 blim =blim*refpts(mp)["bmsy"])
    tac=hcr(mp,params=par,hcrYrs=iYr+seq(interval),tac=TRUE)
    
    #### Operating Model Projectionfor TAC
    om =fwd(om,catch=tac,sr=eql,sr.residuals=sr_deviates,maxF=mean(maxF))  
    }
  
  return(om)}

