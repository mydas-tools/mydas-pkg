sra<-function(object,depletion,nyr=1){
  index=FLQuant(NA,dimnames=dimnames(catch(object)))
  
  index[,seq(nyr)]                        =params(object)["b0"]
  index[,dim(catch(object))[2]+1-seq(nyr)]=depletion
  
  params(object)=params(object)[1:4]
  setParams(object) =index
  setControl(object)=params(object)
  
  res=fit(object,index)
  
  #dpl =stock(bd,0.5)%/%params(bd)["k"]
  #dpl[,ac(c(22:52))]=NA
  
  #control(bd)["q1",2:4]=c(100,1000,10000)
  
  #sra=fit(bd,dpl)
  
  res}