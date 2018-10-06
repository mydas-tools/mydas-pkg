
#http://www.ices.dk/sites/pub/Publication%20Reports/Expert%20Group%20Report/acom/2017/WKMSYCAT34/01.%20WKMSYCAT34%20REPORT%202017.pdf

hcr2o3<-function(yrs,
         catch,
         index,
         fproxy=1,
         bproxy=1,
         rule=c("2o3","slope","loess")){
  
  dmns=dimnames(catch)
  dmns$year=yrs
  
  if (("numeric"%in%is(fproxy))) fproxy =FLQuant(fproxy)
  if (("numeric"%in%is(bproxy))) bproxy =FLQuant(bproxy)

  catch=apply(catch,3:6,mean)
  index=index%/%apply(index,6,mean)
  
  if (rule[1]=="slope"){
    trend=as.FLQuant(ddply(transform(as.data.frame(index),iter=as.numeric(iter)),.(iter),with,
                data.frame(data=coefficients(lm(data~year,data=data.frame(year=year,val=log(data))))[2])))
    r    =exp(trend)}
  else if (rule[1]=="2o3"){
    r = apply(index[,dim(index)[2]-0:1],6,mean)%/%apply(index[,dim(index)[2]-2:4],6,mean)}
  else if (rule[1]=="loess"){
    r = as.FLQuant(ddply(transform(as.data.frame(index),iter=as.numeric(iter)),.(iter),with,{
                    hat=predict(loess(stock~year,data=subset(res,iter==1)))
                    data.frame(data=hat[length(hat)]/hat[length(hat)-1])}))}
    
  tac=catch%*%r%*%fproxy%*%bproxy
  tac=FLQuant(rep(c(tac),each=length(yrs)),dimnames=dmns)
  
  tac}

if (FALSE){
 load(file.path(dirDat,"turbot.RData"))
 res=omSmry(om,eq,lh)
 ggplot(res,aes(fbar,sln))+
   geom_point()+
   geom_smooth(lm="loess")+
   scale_x_log10()+
   scale_y_log10()
 
source('~/Desktop/flr/FLife/R/omOut.R')
 
object=omSmry(om,eq,lh)
object=subset(object,year>=57)

## Quantities
res=transmute(object,
                 iter=iter,
                 year=year,
                 fmsy=fbar/msy_harvest,
                 bmsy=ssb/msy_ssb,
                 msy =catch/msy_yield,
                 pri =rec_hat/virgin_rec)

## Variation
### CV
### AV
### Runs

## Discount Rates

#aav
ddply(res,.(iter), with, aav(msy))
ddply(res,.(iter), with, aav(fmsy))

#green quadrant
ddply(res,.(iter), with, (bmsy>1)*(fmsy<1))


}