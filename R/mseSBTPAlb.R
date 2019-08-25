utils::globalVariables(c("hcrSBT1","hcrSBT2"))
utils::globalVariables(c("hcrSBT2"))

#' mse
#' 
#' @title mseEMPSBT2Alb 
#' 
#' @description 
#' @author Laurence Kell, Sea++
#'  
#' @name mseEMPSBT2Alb
#' 
#' @export mseEMPSBT2Alb
#' @docType methods
#' 
#' @export mseEMPSBT2Alb
#' @rdname mseEMPSBT2Alb
#' 
#' @examples
#' \dontrun{
#' data(pl4)
#' }
#' 
mseEMPSBT2Alb<-function(
  #OM
  om,eq,
  
  #MP
  control=c(k1=0.25,k2=0.25),refYr="missing",
  
  #years over which to run MSE
  interval=3,start=range(om)["maxyear"]-30,end=range(om)["maxyear"]-interval,
  
  #Stochasticity
  srDev,   #=rlnorm(dim(om)[6],FLQuant(0,dimnames=list(                     year=start:end)), srDev),
  uDev,    #=rlnorm(dim(om)[6],FLQuant(0,dimnames=list(                     year=start:end)),  uDev),
  selDev,  #=rlnorm(dim(om)[6],FLQuant(0,dimnames=list(age=dimnames(om)$age,year=start:end)),SelDev),
  
  #Capacity, i.e. F in OM can not be greater than this
  maxF=1.5){ 
  
  if ("FLQuant"%in%is(  uDev)) uDev  =FLQuants(uDev)
  if ("FLQuant"%in%is(selDev)) selDev=FLQuants(selDev)
  
  ## Get number of iterations in OM
  nits=c(om=dims(om)$iter, eq=dims(params(eq))$iter, rsdl=dims(srDev)$iter)
  if (length(unique(nits))>=2 & !(1 %in% nits)) ("Stop, iters not '1 or n' in om")
  if (nits['om']==1) stock(om)=propagate(stock(om),max(nits))
  
  ## Cut in capacity
  maxF=median(apply(fbar(window(om,end=start)),6,max)*maxF)
  
  #### Observation Error (OEM) setup
  nU=max(length(uDev),length(selDev))
  cpue=FLQuants()
  for (iU in seq(nU)){
    yrs=dimnames(window(selDev[[iU]],end=min(start,dims(selDev[[iU]])$maxyear)))$year
    yrs=yrs[yrs%in%dimnames(uDev[[iU]])$year]
    yrs=yrs[yrs%in%dimnames(m(om))$year]
    u         =catch.wt(om)[,yrs]%*%catch.n(om)[,yrs]%/%fbar(om)[,yrs]%*%selDev[[iU]][,yrs]
    cpue[[iU]]=apply(u,2:6,sum)%*%uDev[[iU]][,yrs]}
  
  ## Loop round years
  for (iYr in seq(start,end-interval,interval)){
    cat(iYr,", ",sep="")
    
    ##OEM
    for (iU in seq(nU)){
      yrs=ac(iYr-(interval:1))
      yrs=yrs[yrs%in%dimnames(uDev[[iU]])$year]
      yrs=yrs[yrs%in%dimnames(selDev[[iU]])$year]
      if (length(yrs)>0){
        cpue[[iU]]=window(cpue[[iU]],end=iYr-1)
        u         =(catch.wt(om)[,yrs]%*%catch.n(om)[,yrs]%/%fbar(om)[,yrs])%*%selDev[[iU]][,yrs]
        cpue[[iU]][,yrs]=apply(u,2:6,sum)%*%uDev[[iU]][,yrs]}}
    
    #### Management Procedure
    u=as.FLQuant(ddply(as.data.frame(cpue),.(year,iter), with, data.frame(data=mean(data))))

    tac=hcrSBT2(yrs    =iYr+seq(interval),
                control=control,
                catch  =apply(catch(om)[,ac(iYr-seq(interval)-1)],6,mean),
                cpue   =apply(u[,        ac(iYr-1:interval)],     6,mean),
                ref    =apply(u[,        ac(refYr[1]+-1:1)],      6,mean),
                target =apply(catch(om)[,ac(refYr[2]+-1:1)],      6,mean))
    
    tac[is.na(tac)]=1
    
    #### Operating Model Projectionfor TAC
    #try(save(om,tac,sr,eq,srDev,maxF,file="/home/laurence/Desktop/test3.RData"))
    om =fwd(om,catch=tac,sr=eq,residuals=srDev,effort_max=maxF)  
    print(plot(window(om,end=iYr+dim(tac)[2])))
  }
  
  cat("\n")
  return(om)}
