lbsprFn<-function(len,params,species="",units="cm"){
  
  pars        =new("LB_pars")
  pars@Linf   =c(params["linf"]) 
  pars@L50    =vonB(c(params["a50"]),params) 
  pars@L95    =pars@L50+vonB(c(params["ato95"]),params)
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

lbspr<-function(object,params){
  
  nits=max(dim(object)[6],dim(params)[2])
  
  if (!(dim(object)[6]%in%c(1,nits)|(dim(params)[2]%in%c(1,nits))))
    stop("iters should be equal to n or 1")
  
  res=mdply(data.frame(iter=seq(nits)), function(iter)
    lbsprFn(iter(object,iter)[drop=T],iter(params,iter)))
  res=data.frame(year=dimnames(object)$year,res)
  
  rtn=FLPar(cast(melt(res,id=c("year","iter")),variable~year~iter),units="NA")
  
  rtn}
