simCovar<-function(x,cv,cor,nits=100){
  
  mn=aaply(x,1,mean,na.rm=TRUE)
  y =exp(rmvnorm(nits,       log(mn[dimnames(cor)[[1]]]),
                 cor2cov(cor,log(mn[dimnames(cor)[[1]]])*cv*cv)))
  
  res=FLPar(t(array(unlist(c(y)),c(nits,dim(cor)[1]),
                    dimnames=list(iter=seq(100),params=dimnames(cor)[[1]]))))
  
  mn=propagate(FLPar(mn),nits)
  mn[dimnames(cor)[[1]]]=res
  mn["l50"]=res["linf"]*res["l50linf"]
  lhPar(mn)}

