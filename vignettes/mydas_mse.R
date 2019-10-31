## ----knitr_init, echo=FALSE, results="hide"------------------------------
library(knitr)
## Global options
opts_chunk$set(cache     =!TRUE,
               cache.path='cache/mse/',
               echo      =TRUE,
               eval      =TRUE,
               prompt    =FALSE,
               comment   =NA,
               message   =FALSE,
               warning   =FALSE,
               tidy      =FALSE,
               fig.height=6,
               fig.width =8,
               fig.path  ='tex/mse-')

iFig=0


## ----install, eval=FALSE-------------------------------------------------
## install.packages("mpb", repos = "http://flr-project.org/R")


## ----init, echo=FALSE----------------------------------------------------
library(FLife)


## ----init-2, echo=FALSE--------------------------------------------------
library(FLCore)
library(FLBRP)
library(FLAssess)
library(FLXSA)

library(ggplotFL)

library(FLasher)
library(FLBRP)
library(FLife)
library(mpb)
library(plyr)

theme_set(theme_bw())


## ----turbot--------------------------------------------------------------
turbot=FLPar(c(linf= 59.1,  k=0.28, t0=-0.4, a=0.01111,b=3.15,a50=4.0, l50=43.25),units="NA")


## ------------------------------------------------------------------------
turbot=lhPar(turbot)

turbot


## ----eql-----------------------------------------------------------------
eq=lhEql(turbot)

range(eq)[c("minfbar","maxfbar")]=ceiling(mean(turbot["a1"]))


## ----vectors, echo=FALSE, fig.height=6-----------------------------------
sel<-function(x) 
  catch.sel(x)%/%fapex(catch.sel(x))

ggplot(FLQuants(eq,"m","catch.sel"=sel,"mat","catch.wt"))+
  geom_line(aes(age,data))+
  facet_wrap(~qname,scale="free")+
  scale_x_continuous(limits=c(0,20))+ 
  guides(colour=guide_legend(title="Species",title.position="top"))


## ----eqlcurves, echo=FALSE-----------------------------------------------
plot(eq)


## ----fbar----------------------------------------------------------------
fbar(eq)=refpts(eq)["msy","harvest"]%*%FLQuant(c(rep(.1,19),
                                              seq(.1,2,length.out = 30),
                                              seq(2,1.0,length.out = 10),
                                              rep(1.0,61)))[,1:105]
om=as(eq,"FLStock")

om=fwd(om,fbar=fbar(om)[,-1],sr=eq)


## ----stock-stochastic-u--------------------------------------------------
nits=10
set.seed(3321)
uDev =rlnorm(nits,setPlusGroup(stock.n(eq),10)*0,.2)


## ----stock-stochastic-rec------------------------------------------------
set.seed(1234)
srDev=rlnoise(nits,fbar(om)%=%0,.3,b=0.0)


## ----stock-stochastic-plot, echo=FALSE-----------------------------------
plot(srDev,iter=c(7,2,9))


## ----stock-stochastic-1--------------------------------------------------
om =propagate(om,nits)
oms=FLStocks("Projection"=fwd(om,fbar=fbar(om)[,-1],residuals=srDev,sr=eq))


## ----stock-stochastic-2, echo=FALSE--------------------------------------
plot(oms[["Projection"]],iter=1:3)+
  theme(legend.position="none")


## ----hcr,echo=TRUE-------------------------------------------------------
library(kobe)

hcr= data.frame(stock  =c(0.0 ,0.1 , 0.6,2.0), 
                harvest=c(0.01,0.01, 0.7,0.7))
kobePhase()+
  geom_line(aes(stock,harvest),data=hcr,col="orange",size=2)


## ----xsa-xtest-----------------------------------------------------------
library(FLXSA)
mp=window(setPlusGroup(oms[["Projection"]],10),end=80)

xsaControl=FLXSA.control(tol    =1e-09, maxit   =150, 
                         min.nse=0.3,   fse     =1.0, 
                         rage   =1,     qage    =6, 
                         shk.n  =TRUE,  shk.f   =TRUE, 
                         shk.yrs=1,     shk.ages=4, 
                         window =10,    tsrange =10, 
                         tspower= 0,
                         vpa    =FALSE)


## ----xsa-xtest-2---------------------------------------------------------
idx=FLIndex(index=stock.n(mp)%*%uDev[,dimnames(stock.n(mp))$year])
range(idx)[c("plusgroup","startf","endf")]=c(NA,0.1,.2)

xsa=FLXSA(mp,idx,
          control=xsaControl,diag.flag=FALSE)
range(xsa)[c("min","max","plusgroup")]=range(mp)[c("min","max","plusgroup")]
mp=mp+xsa 

sr=fmle(as.FLSR(mp,model="bevholt"),control=list(silent=TRUE))
rf=FLBRP(mp,sr) 


## ----xsa-xtest-plot------------------------------------------------------
plot(FLStocks("Stock\nAssessment"=mp,
              "Operating\nModel" =window(oms[["Projection"]],end=80)))


## ----xsa-mse-------------------------------------------------------------
library(mydas)
save(oms,eq,mp,xsaControl,rf,srDev,uDev,file="/home/laurence/tmp/oms.RData")  

## ------------------------------------------------------------------------
oms["Age"]=mseXSA(oms[["Projection"]],eq,  
                    mp,xsaControl,rf=rf,       
                    sr_deviances=srDev,u_deviances=uDev,   
                    start=75,end=103,maxF=1.0)      


## ----xsa-mse-plot, echo=FALSE--------------------------------------------
plot(oms[["Age"]],iter=1:3)+
  theme(legend.position="none")


## ---- echo=FALSE---------------------------------------------------------
library(plyr)
library(dplyr)
library(reshape)
library(ggplot2)

library(FLCore)
library(ggplotFL)
library(FLasher)
library(FLBRP)
library(FLife)
library(mpb)
library(mydas)


## ----biodyn-2------------------------------------------------------------
## MP
library(mpb)

mp=setMP(as(window(om,end=75),"biodyn"),
         r =   0.25,
         k =1000.0,
         b0=   0.9,
         p =  -0.6)


## ----biodyn-3------------------------------------------------------------
library(FLasher)
library(mpb)

nits=dims(om)$iter
set.seed(1234)
srDev =rlnoise(nits,FLQuant(0,dimnames=dimnames(iter(catch(om),1))),0.3,b=0.0)
uDev  =rlnoise(nits,FLQuant(0,dimnames=dimnames(iter(catch(om),1))),0.2,b=0.0)
selDev=rlnoise(nits,FLQuant(0,dimnames=dimnames(iter(    m(om),1))),0.2,b=0.0)

eq=FLCore::iter(eq,seq(nits))

oms["Biomass"]=mydas:::mseMPB2(om,eq,mp,start=75,end=100,ftar=0.5,sr_deviances=srDev,
                       u_deviances=uDev,sel_deviances=selDev)


## ----biodyn-mse-plot, echo=FALSE-----------------------------------------
plot(window(oms[["Biomass"]],end=100),iter=1:3)+
  theme(legend.position="none")


## ----emp-----------------------------------------------------------------
control=FLPar(k1=0.5,k2=0.5,gamma=1)
oms["Emprirical"]=mydas:::mseSBTD(om,eq,control=control,
                            sr_deviances=srDev,u_deviances=uDev,
                            start=75,end=100)


## ----emp-mse-plot, echo=FALSE--------------------------------------------
plot(window(oms[["Emprirical"]],end=100),iter=1:3)+
  theme(legend.position="none")


## ----eval=FALSE----------------------------------------------------------
## res=ldply(oms,omSmry,eq)

