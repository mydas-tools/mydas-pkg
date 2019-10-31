utils::globalVariables(c("lhEql","lhPar","lopt","maply","lambda"))

#' popdyn
#' 
#' @title popdyn 
#' 
#' @description Calculates popdyn, i.e. r and reference points, based on life history parameters
#' @author Laurence Kell, Sea++
#' 
#' @name popdyn 
#' 
#' @param object blah,blah,blah,...
#' @param eq blah,blah,blah,...
#' 
#' @export popdyn
#' @docType methods
#' 
#' @rdname popdyn
#' 
#' @examples
#' \dontrun{
#' par=FLPar(linf=100)
#' popdyn(par)
#' }
  
  
popdyn<-function(object,eq=lhEql(lhPar(object))){

  par=lhPar(object)

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
  
  rtn=rbind(rtn,FLPar(fm=fm))

  #generation time
  eq@fbar[]=0
  n=stock.n(eq)[,1]
  gt=apply((stock.wt(eq)%*%mat(eq)%*%n)%*%ages(stock.wt(eq)),6,sum)%/%
       apply(stock.wt(eq)%*%mat(eq)%*%n,6,sum)
  gt=FLPar(gt=array(c(gt),c(1,length(c(gt)))))
  rtn=rbind(rtn,gt)
  
  if (!("l50"%in%dimnames(rtn)[[1]])){
     l50=FLPar("l50"=vonB(age=c(par["a50"]),par))
     rtn=rbind(rtn,l50)}

  lns=lenFn(eq,par)
  if (dim(lns)[1]==1) lns=as(lns,"FLPar") else lns=as(lns[,-1],"FLPar")

  rtn=rbind(rtn,lns["sln"])

  dimnames(rtn)$params[dimnames(rtn)$params=="sln"]="slv"
  warning("add age etc")
  
 
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
    
  #LF=M
  lfm=FLPar("lfm"=c(0.75*rtn["lc"]+0.25*rtn["linf"]))
  rtn=rbind(rtn,lfm)
    
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
  
  avirgin<-FLPar(apply(stock.n(eq)%*%ages(stock.n(eq)),6,sum)%/%
                   apply(stock.n(eq),2:6,sum))
  wvirgin<-FLPar(apply(stock.n(eq)%*%stock.wt(eq),6,sum)%/%
                   apply(stock.n(eq),2:6,sum))
  
  eq@fbar[]=refpts(eq)["msy","harvest"]
  amsy<-FLPar(apply(stock.n(eq)%*%ages(stock.n(eq)),6,sum)%/%
                apply(stock.n(eq),2:6,sum))
  wmsy<-FLPar(apply(stock.n(eq)%*%stock.wt(eq),6,sum)%/%
                apply(stock.n(eq),6,sum))
  
  rtn=rbind(rtn,avirgin)
  rtn=rbind(rtn,wvirgin)
  rtn=rbind(rtn,amsy)
  rtn=rbind(rtn,wmsy)
  
  nms=c("fmsy","msy","bmsy","v","spr0","l50","lc","lopt","gt","slmsy","clmsy","r","rc","fm","mk")
  
  rtn=rbind(rtn[dimnames(rtn)$params%in%nms][nms[nms%in%dimnames(rtn)$params]],
            rtn[!dimnames(rtn)$params%in%nms])
  
  return(rtn)}
