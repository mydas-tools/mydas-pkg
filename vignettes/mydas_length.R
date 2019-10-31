## ---- ini, echo=FALSE, results='hide', message=FALSE, warnings=FALSE, cache=FALSE,eval=!TRUE----
## library(knitr)
## source("R/ini.R")


## ----knitr_init, echo=FALSE, results="hide", eval=!FALSE-----------------
library(knitr)  
## Global options
opts_chunk$set(cache     =TRUE,
               echo      =FALSE,
               eval      =TRUE,
               prompt    =FALSE,
               comment   =NA,
               message   =FALSE,
               warning   =FALSE,
               tidy      =TRUE,
               fig.height=4,
               fig.width =6,
               fig.path  ="tex/simtest/len-",
               cache.path="cache/simtest/len/")


## ------------------------------------------------------------------------
options(digits=3)

iFig=0


## ---- pkgs---------------------------------------------------------------
# Load  packages
library(ggplot2)
library(plyr)
library(reshape)

library(popbio)

library(FLCore)
library(ggplotFL)
library(FLBRP)
library(FLasher)
library(FLife)
library(mydas)


## ---- om-----------------------------------------------------------------
lh=FLPar(c(linf= 59.1,  k=0.28, t0=-0.4, a=0.01111,b=3.15,a50=4.0, l50=43.25),units="NA")
lh=lhPar(lh)
eq=lhEql(lh)

gTime=c(round(gt(eq)))
fbar(eq)=refpts(eq)["msy","harvest"]%*%FLQuant(c(rep(.1,19),
                                              seq(.1,2,length.out=30)[-30],
                                              seq(2,1.0,length.out=gTime)[-1],
                                              rep(1.0,61)))[,1:105]

om=as(eq,"FLStock")
om=fwd(om,f=fbar(om)[,-1], sr=eq)


## ---- ts-----------------------------------------------------------------
plot(FLQuants(om, 
          "f" =   function(x) fbar(x)%/%refpts(eq)["msy","harvest"], 
          "ssb" = function(x) ssb(x)%/%refpts( eq)["msy","ssb"], 
          "catch"=function(x) landings(x)%/%refpts(eq)["msy","yield"],
          "rec" = function(x) rec(x)%/%refpts( eq)["msy","rec"])) + 
  geom_hline(aes(yintercept=1),col="red",linetype=2)+
  theme_bw() 


## ---- len----------------------------------------------------------------
ts   =omSmry(om,eq,lh)

head(ts)


## ---- len-plot,fig.height=3----------------------------------------------
mnLen=as.FLQuant(with(ts,data.frame(data=cln,year=year,iter=iter)))

plot(mnLen)


## ---- priors-------------------------------------------------------------
prior=priors(lh)
prior


## ---- mlz----------------------------------------------------------------
library(MLZ)

res=mydas:::mlz(mnLen[,ac(40:60)],prior)

res


## ---- lbspr--------------------------------------------------------------
library(LBSPR)


## ---- alk----------------------------------------------------------------
ak=alk(lh)    


## ---- lfd----------------------------------------------------------------
lfd=lenSample(catch.n(om)[,20:65],ak,nsample=500)


## ---- oem----------------------------------------------------------------
ggplot(melt(lfd[,seq(1,45,10)]))+
  geom_histogram(aes(length,weight=value),binwidth=1)+
  facet_grid(year~iter,scale="free")+
  xlab("Length (cm)")+ylab("Frequency")+
  coord_cartesian(xlim=c(0,mean(lh["linf"])))


## ------------------------------------------------------------------------

setGeneric('lbspr', function(object,params,...) standardGeneric('lbspr'))

lbsprFn<-function(len,params,species="",units="cm"){
  
  pars        =new("LB_pars")
  pars@Linf   =c(params["linf"]) 
  pars@L50    =vonB(c(params["a50"]),params) 
  pars@L95    =vonB(c(params["a50"])+c(params["ato95"]),params)
  pars@MK     =c(params["mk"])
  pars@Species=species
  pars@L_units=units
  
  #labs=dimnames(len)[[1]]
  #brks=cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ),
  #           upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))
  #mid=aaply(brks,1,mean)
  
  LBlen       =new("LB_lengths")
  LBlen@LMids =as.numeric(dimnames(len)[[1]])
  LBlen@LData =len
  LBlen@Years =as.numeric(dimnames(len)[[2]])
  LBlen@NYears=dim(len)[2] 
  
  res=LBSPRfit(pars,LBlen,verbose=FALSE)
  
  res@Ests}


setMethod('lbspr', signature(object="data.frame",params='FLPar'),
          function(object,params,...){
  
  nits=max(dim(object)[6],dim(params)[2])
  
  if (!(dim(object)[6]%in%c(1,nits)|(dim(params)[2]%in%c(1,nits))))
    stop("iters should be equal to n or 1")
  
  res=mdply(data.frame(iter=seq(nits)), function(iter)
    lbsprFn(iter(object,iter)[drop=T],iter(params,iter)))
  res=data.frame(year=dimnames(object)$year,res)
  
  
  rtn=FLPar(cast(melt(res,id=c("year","iter")),variable~year~iter),units="NA")
  
  rtn})

setMethod('lbspr', signature(object="FLQuant",params='FLPar'),
          function(object,params,...){
            
            nits=max(dim(object)[6],dim(params)[2])
            
            if (!(dim(object)[6]%in%c(1,nits)|(dim(params)[2]%in%c(1,nits))))
              stop("iters should be equal to n or 1")
            
            res=mdply(data.frame(iter=seq(nits)), function(iter)
              lbsprFn(iter(object,iter)[drop=T],iter(params,iter)))
            res=data.frame(year=dimnames(object)$year,res)
            
            
            rtn=FLPar(cast(melt(res,id=c("year","iter")),variable~year~iter),units="NA")
            
            rtn})


## ---- sa-----------------------------------------------------------------
lb=lbspr(lfd,prior)  


## ---- spr----------------------------------------------------------------
ggplot(melt(sweep(lb["SPR"],c(1,3),lb["SPR","40"],"/")))+
  geom_boxplot(aes(ac(year),value))+
  scale_x_discrete(breaks=seq(20,60,10))+
  ylab("SPR")+xlab("Year")+theme_bw()


## ---- fm-----------------------------------------------------------------
ggplot(melt(sweep(lb["FM"],c(1,3),lb["FM","40"],"/")))+
  geom_boxplot(aes(ac(year),value))+
  scale_x_discrete(breaks=seq(20,60,10))+
  ylab("F")+xlab("Year")+theme_bw()

