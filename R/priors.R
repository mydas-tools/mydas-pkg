priors<-function(object,eq=lhEql(lhPar(object))){
  
  par=lhPar(object)

  ## return
  rtn=par[c("linf","k","t0","a50","ato95","a","b","s","v")]
  
  ## SRR
  sr=params(eq)
  dimnames(sr)$params
  dimnames(sr)$params=c("alpha","beta")
  rtn=rbind(rtn,sr)
  
  spr=FLPar(spr0=array(c(spr0(eq)),c(1,dim(rtn)[2])))
  rtn=rbind(rtn,spr)
  
  #msy  
  rfs=FLPar(refpts(eq)["msy",c("harvest","yield","ssb"),drop=T])
  names(rfs)[1]="params"
  dimnames(rfs)[[1]]=c("fmsy","msy","bmsy")
  rtn=rbind(rtn,rfs)
  
  #lopt
  growth=vonB
  lop=lopt(par)
  rtn=rbind(rtn,lop)
  
  #lc
  lc=vonB(as.FLQuant(c(par["a50"]-par["ato95"]),dimnames=list(iter=dimnames(par)$iter)),par)
  lc=FLPar(lc=array(c(lc),c(1,length(c(lc)))))
  
  rtn=rbind(rtn,lc)
  
  r=maply(seq(dims(eq)$iter), function(x) 
    log(lambda(leslie(iter(eq,x),fbar=c(refpts(eq)["crash","harvest",x]))[drop=TRUE])))
  r=FLPar(r=array(c(r),c(1,length(c(r)))))
  rtn=rbind(rtn,r)
  
  #"fmsy/m"
  eq@fbar=fbar(eq)[,1]
  eq@fbar[]=refpts(eq)["msy","harvest"]
  fm=apply(fbar(eq)%/%m(eq),6,mean)
  fm=FLPar(fm=array(c(fm),c(1,length(c(fm)))))
  rtn=rbind(rtn,fm)
  
  #"m/k",priors<-function(object,eq=lhEql(lhPar(object))){
  
  par=lhPar(object)
  
  ## return
  rtn=par[c("linf","k","t0","a50","ato95","a","b","s","v")]
  
  ## SRR
  sr=params(eq)
  dimnames(sr)$params
  dimnames(sr)$params=c("alpha","beta")
  rtn=rbind(rtn,sr)
  
  spr=FLPar(spr0=array(c(spr0(eq)),c(1,dim(rtn)[2])))
  rtn=rbind(rtn,spr)
  
  #msy  
  rfs=FLPar(refpts(eq)["msy",c("harvest","yield","ssb"),drop=T])
  names(rfs)[1]="params"
  dimnames(rfs)[[1]]=c("fmsy","msy","bmsy")
  rtn=rbind(rtn,rfs)
  
  #lopt
  growth=vonB
  lop=lopt(par)
  rtn=rbind(rtn,lop)
  
  #lc
  lc=vonB(as.FLQuant(c(par["a50"]-par["ato95"]),dimnames=list(iter=dimnames(par)$iter)),par)
  lc=FLPar(lc=array(c(lc),c(1,length(c(lc)))))
  
  rtn=rbind(rtn,lc)
  
  r=maply(seq(dims(eq)$iter), function(x) 
    log(lambda(leslie(iter(eq,x),fbar=c(refpts(eq)["crash","harvest",x]))[drop=TRUE])))
  r=FLPar(r=array(c(r),c(1,length(c(r)))))
  rtn=rbind(rtn,r)
  
  rc=maply(seq(dims(eq)$iter), function(x) 
    log(lambda(leslie(iter(eq,x),fbar=c(refpts(eq)["msy","harvest",x]))[drop=TRUE])))
  rc=FLPar(rc=array(c(rc),c(1,length(c(rc)))))
  rtn=rbind(rtn,rc)
  
  #"fmsy/m"
  eq@fbar=fbar(eq)[,1]
  eq@fbar[]=refpts(eq)["msy","harvest"]
  fm=apply(fbar(eq)%/%m(eq),6,mean)
  fm=FLPar(fm=array(c(fm),c(1,length(c(fm)))))
  rtn=rbind(rtn,fm)
  
  #"m/k",
  mk=apply(m(eq)%/%par["k"],6,mean)
  mk=FLPar(mk=array(c(mk),c(1,length(c(mk)))))
  rtn=rbind(rtn,mk)
  
  rtn}