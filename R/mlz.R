mlzFn<-function(year,len,lc,linf,k,ss=500){
  
  dat=new("MLZ_data", 
          Year         =unique(year), 
          MeanLength   =len,
          ss           =rep(ss,length(year)),
          length.units ="cm")
  
  dat@Lc    =lc
  dat@vbLinf=linf
  dat@vbK   =k
  
  hat=ML(dat,ncp=1,figure=FALSE)
  
  res=hat@estimates
  res}

mlz<-function(object,params,ss=500){
  nits=max(dim(object)[6],dims(params)$iter) 
  
  rtn=NULL  
  for (i in seq(nits)){
    
    dat=as.data.frame(iter(object,i))
    
    res=mlzFn(
         year=dat$year,
         len =dat$data,
         lc  =c(iter(prior["lc"],  i)),
         linf=c(iter(prior["linf"],i)),
         k   =c(iter(prior["k"],   i)),
         ss=ss)
    
    rtn=rbind(rtn,cbind(iter=i,res))}
  
  dimnames(rtn)[[2]][2:3]=c("hat","se")
  
  z=data.frame(subset(rtn,substr(dimnames(rtn)[[1]],1,1)=="Z"))
  z=ddply(z,.(iter), transform, block=seq(length(hat)))
  if (!("iter"%in%names(z))) z=data.frame(iter=1,z)
  z=as(daply(z,.(block), with, as(data.frame(hat=hat,se=se,iter=iter),"FLPar")),"FLPar")
  
  y=subset(rtn,substr(dimnames(rtn)[[1]],1,1)!="Z")
  y=data.frame(y,name=substr(dimnames(y)[[1]],1,1))
  y=ddply(y,.(iter,name), transform, block=seq(length(hat)))
  if (!("iter"%in%names(y))) y=data.frame(iter=1,y)
  y=as(daply(y,.(block,name), with, as(data.frame(hat=hat,se=se,iter=iter),"FLPar")),"FLPar")
  
  FLPars("z"=z,"year"=y)}  
