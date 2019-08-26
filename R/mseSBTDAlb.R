utils::globalVariables(c("mlply","oem","yrRng","control","control<-","maxF.","oem"))

#' mseEMPSBT1Alb
#' 
#' @title mseEMPSBT1Alb
#' 
#' @description 
#' @author Laurence Kell, Sea++
#'  
#' @name mseEMPSBT1Alb 
#' 
#' @param om \code{FLStock} object as the operating model
#' @param eq blah,blah,blah,...
#' @param control blah,blah,blah,...
#' @param interval blah,blah,blah,...
#' @param start \code{numeric}  default is range(om)["maxyear"]-30
#' @param end \code{numeric}  default is range(om)["maxyear"]-interval
#' @param u_deviates blah,blah,blah,...
#' @param sr_deviates  blah,blah,blah,...
#' @param sel_deviates  blah,blah,blah,...
#' @param maxF blah,blah,blah,...
#' 
#' @export mseEMPSBT1Alb 
#' @docType methods
#' 
#' @export mseEMPSBT1Alb 
#' @rdname mseEMPSBT1Alb 
#' 
#' @examples
#' \dontrun{
#' data(pl4)
#' }
#' 
mseEMPSBT1Alb<-function(
  #OM
  om,eq,
  
  #MP
  control=c(k1=2.5,k2=2.5,gamma=1),
  
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
  for (iYr in seq(start,end-interval,interval)){
    cat(iYr,", ",sep="")
    
    ##OEM
    for (iU in seq(nU)){
      yrs=ac(iYr-(interval:1))
      yrs=yrs[yrs%in%dimnames(u_deviates[[iU]])$year]
      yrs=yrs[yrs%in%dimnames(sel_deviates[[iU]])$year]
      if (length(yrs)>0){
        cpue[[iU]]=window(cpue[[iU]],end=iYr-1)
        u         =(catch.wt(om)[,yrs]%*%catch.n(om)[,yrs]%/%fbar(om)[,yrs])%*%sel_deviates[[iU]][,yrs]
        cpue[[iU]][,yrs]=apply(u,2:6,sum)%*%u_deviates[[iU]][,yrs]}}
    
    #### Management Procedure
    u=as.FLQuant(ddply(as.data.frame(u),.(year,iter), with, data.frame(data=mean(data))))
    
    tac=hcrSBT1(iYr+seq(interval),
                control=control,
                u[,ac(ac(iYr-(3:1)))],
                catch(om)[,ac(iYr-(2:1))])
    
    #### Operating Model Projectionfor TAC
    #try(save(om,tac,sr,eq,sr_deviates,maxF,file="/home/laurence/Desktop/test3.RData"))
    om =fwd(om,catch=tac,sr=eq,residuals=sr_deviates,effort_max=maxF)  
    #print(plot(window(om,end=iYr+dim(tac)[2])))
  }
  
  cat("\n")
  return(om)}
