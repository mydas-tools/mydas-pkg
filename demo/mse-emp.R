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


theme_set(theme_bw())

dirMy ="/home/laurence/Desktop/sea++/mydas/tasks/task4"
dirDat=file.path(dirMy,"data")

## OM
load(file.path(dirDat,"turbot.RData"))

om=window(om,start=20,end=90)

## MP
#source('~/Desktop/flr/mpb/R/mseEMP.R')

nits=dims(om)$iter
set.seed(1234)
sr_deviates=FLife:::rlnoise(nits,FLQuant(0,dimnames=list(year=1:100)),0.3,b=0.0)
u_deviates =FLife:::rlnoise(nits,FLQuant(0,dimnames=list(year=1:100)),0.2,b=0.0)
eq=iter(eq,seq(nits))

msesbt1 =mseEMP(om,eq,control=c(k1=1.5,k2=2.0,gamma=1),start=60,end=100)  
