## ---- ini, echo=FALSE, results='hide', message=FALSE, warnings=FALSE, cache=FALSE,eval=!TRUE----
## library(knitr)
## #source("R/ini.R")


## ----knitr_init, echo=FALSE, results="hide",eval=TRUE--------------------
library(knitr)
## Global options
opts_chunk$set(echo    =!TRUE,
               eval    =TRUE,
               cache   =TRUE,
               cache.path="cache/oem/",
               prompt  =FALSE,
               comment =NA,
               message =FALSE,
               tidy    =FALSE,
               warnings=FALSE,
               fig.height=4.5,
               fig.width =6,
               fig.path  ="tex/oem-")
iFig=0


## ------------------------------------------------------------------------
library(plyr)
library(reshape)
library(ggplot2)


## ---- eval=FALSE---------------------------------------------------------
## install.packages(c("FLCore","FLFishery","FLasher","FLBRP","mpb","FLife"),
##              repos="http://flr-project.org/R")


## ------------------------------------------------------------------------
library(FLCore)
library(FLasher)
library(FLBRP)
library(FLife)


## ---- eval=FALSE---------------------------------------------------------
## install.packages("devtools",dependencies=TRUE)


## ---- echo=TRUE, eval=FALSE----------------------------------------------
## library(devtools)
## 
## devtools::install_github("lauriekell/mydas-pkg")


## ---- echo=TRUE----------------------------------------------------------
library(mydas) 


## ---- om-lh--------------------------------------------------------------
lh=FLPar(c(linf= 59.1,  k=0.28, t0=-0.4, a=0.01111,b=3.15,a50=4.0,l50=43.25),units="NA")

lh=lhPar(lh)   


## ---- om-eql, echo=TRUE--------------------------------------------------
eq=lhEql(lh)


## ---- om-f, fig.height=3, fig.width=6, echo=TRUE-------------------------
gTime=c(round(gt(eq))) 

fbar(eq)=refpts(eq)["msy","harvest"]%*%FLQuant(c(rep(.1,19),
                                                 seq(.1,2,length.out=30),
                                                 seq(2,1.0,length.out=gTime)[-1],
                                                 rep(1.0,61)))[,1:105]
plot(fbar(eq))


## ------------------------------------------------------------------------
om=as(eq,"FLStock")

om=fwd(om,f=fbar(om)[,-1], sr=eq)

plot(om)


## ---- survey, fig.height=4-----------------------------------------------
plot(stock(om))


## ---- cpue, fig.height=6-------------------------------------------------
plot(FLQuants(om,Survey=stock,
                 CPUE=function(x) catch(x)/fbar(x)))


## ---- cpue2, fig.height=4, fig.width=6-----------------------------------
library(plyr)

dat=as.data.frame(FLQuants(om,Survey=stock,
                 CPUE=function(x) catch(x)/fbar(x)))
dat=ddply(dat,.(qname), transform, index=data/mean(data))
ggplot(dat)+geom_line(aes(year,index,col=qname))+
  theme(legend.position="bottom")


## ---- trend, eval=TRUE---------------------------------------------------
trend<-function(object,bias=0.02) 
  object%*%FLQuant(cumprod(1+rep(bias,dim(object)[2])),dimnames=dimnames(object))


## ------------------------------------------------------------------------
plot(FLQuants(stock(om)%=%1,trend(stock(om)%=%1)))


## ------------------------------------------------------------------------
plot(FLQuants(catch(om)%/%fbar(om),trend(catch(om)%/%fbar(om))))


## ---- hyperstability, eval=TRUE------------------------------------------
hyperstability<-function(object,omega=1,ref=apply(object,c(1,3:6),mean)) 
  ref%*%((object%/%ref)^omega)


## ---- hyper, fig.height=4------------------------------------------------
om2=window(om,start=25,end=50)

plot(catch(om2)%/%fbar(om2))+
  geom_line(aes(year,data),
            data=as.data.frame(hyperstability(catch(om2)%/%fbar(om2),omega=0.5)),col="red")


## ----plot,echo=TRUE, eval=TRUE-------------------------------------------
plot(FLQuants(om,"Stock"   =stock,
                 "Adults"   =ssb,
                 "Recruits" =rec,
                 "Juveniles"=function(x) stock(x)-ssb(x),
                 "CPUE"     =function(x) catch(x)%/%fbar(x)))


## ----example2, fig.height=4, echo=TRUE, eval=TRUE------------------------
cv=rlnorm(100,log(stock(om)),0.3)
plot(cv, iter=1)


## ----example3, fig.height=4, echo=TRUE, eval=TRUE------------------------
cv=rlnoise(100,log(stock(om)),0.3,0.8)
plot(cv, iter=1)


## ----example9,echo=TRUE,eval=TRUE----------------------------------------
set.seed(1234)

u =FLQuants("Unbiased"       =rlnorm(100,log(stock(om)),.3),
            "AR"             =rlnoise(100,log(stock(om)),.3,b=.7),
            "Hyperstability" =rlnorm(100,log(hyperstability(stock(om),0.5)),.3),
            "Trend"          =rlnorm(100,log(trend(stock(om),0.01)),.3),
            "Hetroscedascity"=rlnorm(100,log(stock(om)),.3)*trend(stock(om)%=%1,0.01),
            "Juvenile"       =rlnorm(100,log(stock(om)-ssb(om)),.3),
            "Mature"         =rlnorm(100,log(ssb(om)),.3),
            "Numbers"        =rlnorm(100,log(apply(stock.n(om),2:6,sum)),.3))

u=FLQuants(llply(u,function(x) x/mean(x)))
u=ldply(u,as.data.frame)

u.=ddply(u,.(year,.id), with, quantile(data,na.rm=T))
ggplot()+
  geom_line(aes(year,data,col=factor(iter)),
            data=subset(u,iter%in%c(2,11)))+
  geom_ribbon(aes(year,ymin=`25%`,ymax=`75%`),data=u.,col="grey",alpha=.5)+
  facet_wrap(~.id,ncol=2)+
  theme_bw()+theme(legend.position="none")


## ------------------------------------------------------------------------
mnLn=apply(wt2len(stock.wt(om),lh)%*%stock.n(om),c(2,6),sum)%/%apply(stock.n(om),c(2.,6),sum)


## ------------------------------------------------------------------------
plot(FLQuants(mnLn,fbar(om)))


## ---- turbot-alk,echo=TRUE,eval=TRUE-------------------------------------
ak=alk(lh)  


## ---- lfd----------------------------------------------------------------
lfd=lenSample(catch.n(om)[,20:65],ak,nsample=500)


## ---- turbot-oemplot,echo=TRUE,eval=TRUE---------------------------------
ggplot(melt(lfd[,seq(1,45,10)]))+
  geom_histogram(aes(length,weight=value),binwidth=1)+
  facet_grid(year~iter,scale="free")+
  xlab("Length (cm)")+ylab("Frequency")+
  coord_cartesian(xlim=c(0,mean(lh["linf"])))


## ---- devtools, echo=TRUE, eval=FALSE------------------------------------
## 	library(devtools)
## 	install_github('flr/FLPKG')

