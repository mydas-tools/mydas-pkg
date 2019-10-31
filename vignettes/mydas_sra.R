## ---- ini, echo=FALSE, results='hide', message=FALSE, warnings=FALSE, cache=FALSE,eval=!TRUE----
## library(knitr)


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
               fig.height=6,
               fig.width =8,
               fig.path  ="tex/simtest/srar-",
               cache.path="cache/simtest/sra/")


## ------------------------------------------------------------------------
options(digits=3)

iFig=0


## ---- pkgs---------------------------------------------------------------
# Load  packages
library(ggplot2)
library(plyr)
library(reshape)

library(FLCore)
library(ggplotFL)
library(FLBRP)
library(FLasher)
library(FLife)

library(mydas)


## ------------------------------------------------------------------------
lh=FLPar(c(linf= 59.1,  k=0.28, t0=-0.4, a=0.01111,b=3.15,a50=4.0, l50=43.25),units="NA")
lh=lhPar(lh)
eq=lhEql(lh)

gTime=c(round(gt(eq)))

fbar(eq)=refpts(eq)["msy","harvest"]%*%FLQuant(c(rep(.1,19),
                                              seq(.1,2,length.out=30),
                                              seq(2,1.0,length.out=gTime)[-1],
                                              rep(1.0,61)))[,1:105]
om=as(eq,"FLStock")
om=fwd(om,f=fbar(om)[,-1], sr=eq)


## ----turbot-ts-ref-------------------------------------------------------
plot(FLQuants(om, 
          "f" =   function(x) fbar(x)%/%refpts(eq)["msy","harvest"], 
          "ssb" = function(x) ssb(x)%/%refpts( eq)["msy","ssb"], 
          "catch"=function(x) landings(x)%/%refpts(eq)["msy","yield"],
          "rec" = function(x) rec(x)%/%refpts( eq)["msy","rec"])) + 
  geom_hline(aes(yintercept=1),col="red",linetype=2)+
  theme_bw() 


## ------------------------------------------------------------------------
library(popbio)
prior=popdyn(lh)

prior


## ------------------------------------------------------------------------
library(mpb)

## Bug, need to add to NAMESPACE
setMP=mpb:::setMP
p=mpb:::p

om=window(om,end=55)
mp=setMP(as(om,"biodyn"),
         r =median(prior["r"],na.rm=T),
         k =median(prior["v"],na.rm=T),
         b0=0.8,
         p =median(p(prior["bmsy"]/prior["v"]),na.rm=TRUE))


## ------------------------------------------------------------------------
mp=fit(mp,stock(mp)[,20:54])


## ------------------------------------------------------------------------
plot(mp)


## ------------------------------------------------------------------------
plot(as(list("MP"=mp,"OM"=as(om,"biodyn")),"biodyns"))


## ---- sra----------------------------------------------------------------
sra=mp
dplIndex =window(stock(mp,0.5)%/%params(sra)["k"],start=20,end=54)
dplIndex[,ac(c(22:52))]=NA

## change starting values and bounds for catchability
control(sra)["q1",2:4]=c(100,1000,10000)

sra=fit(sra,dplIndex)

## ------------------------------------------------------------------------
plot(as(list("SRA"=sra,"MP"=mp,"OM"=as(om,"biodyn")),"biodyns"))


## ---- turbot-------------------------------------------------------------
dplIndex[,"54"]=0.2
sra2=fit(sra,dplIndex)

## ------------------------------------------------------------------------
plot(as(list("SRA-2"=sra2,"SRA"=sra,"MP"=mp,"OM"=as(om,"biodyn")),"biodyns"))

