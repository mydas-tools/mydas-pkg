#' @title refs
#' 
#' @description 
#' 
#' Returns a variety of reference points based on spawner per recruit (S/R) and stock recruitment
#' relationships.
#' 
#' $F_{MSY}$ corresponds to the level of exploitation that provides the maximum, derived from a 
#' spawner and yield curves combined with a stock recruitment relationship, while $F_{Crash}$ is 
#' the level of F that will drive the stock to extinction. 
#' 
#' Reference points can also be derived from the spawner and yield curves alone. For example
#' $F_{0.1}$ is a proxy for $F_{MSY}$ and is the fishing mortality that corresponds to a point 
#' on the yield per recruit curve where the slope is 10\% of that at the origin; $F_{Max}$ is 
#' the maximum of the yield per recruit and $SPR_O$ is the spawner per recruit at virgin biomass, 
#' and $SPR_30$ corressponds to the point on the curve where SPR is 30\% of $SPR_0$. In these 
#' cases the biomass, ssb and yield values are derived by multiplying the per recruit values by 
#' the average recruitment. 
#'
#' @param object \code{FLStock} 
#' @param ... other arguments
#' @aliases refs refs-method refs,FLStock-method
#' @return \code{data.frame} 
#' 
#' @export
#' @rdname refs
#' 
#' @examples
#' \dontrun{
#' 
#'  library(FLCore)
#'  library(FLBRP)
#'  library(ggplotFL)
#'  library(FLife)
#'  library(popbio)
#'  library(plyr)
#'  
#'  data(ple4)
#'  
#'  refs(ple4)
#' } 

setGeneric('refs', function(object,...) standardGeneric('refs'))

old<-function(object,s=NULL){
  
  if (!requireNamespace("popbio", quietly = TRUE)) {
    stop("Package \"popbio\" needed for this function to work. Please install it.",
         call. = FALSE)}
  
  warn=options()$warn
  options(warn=-1)
  
  ## msy
  eql       =FLBRP(object)
  model(eql)=bevholt()$model
  
  #steepness fixed
  if (!is.null(s)){
    spr0=mean(spr0(FLBRP(object)))
    srr =as.FLSR(object,model="bevholtSV")
    
    upper(srr)[1:2]=1e12
    srr=fmle(srr,
             fixed=list(s=s,spr0=spr0),
             control=list(silent=TRUE),
             method="Brent")
    params(eql)=ab(params(srr),"bevholt")[c("a","b")]
  }else{ 
    srr=as.FLSR(object,model="bevholt")
    lower(srr)[1:2]=0.00001
    srr=fmle(srr,control=list(trace=FALSE),method="L-BFGS-B")
    params(eql)=params(srr)
    }
  
  eql=brp(eql)

  ## per recruit regime shift
  eql1        =eql
  model(eql1) =geomean()$model
  params(eql1)=FLPar(a=NA)
  params(eql1)=propagate(params(eql1),dim(refpts(eql))[3])
  params(eql1)[]=unlist(c(ddply(rod(rec.obs(eql1)),.(iter), function(x)
    mean(subset(x,regime==max(as.numeric(regime)))$data))["V1"]))
  refpts(eql1)=refpts(eql1)[c("f0.1","fmax","spr.30","virgin"),]
  dimnames(refpts(eql1))[[1]][4]="spr.100"

  eql1=brp(eql1)

  ## per recruit stationarity
  eql2        =eql
  model(eql2) =geomean()$model
  params(eql2)=FLPar(a=NA)
  params(eql2)=propagate(params(eql2),dim(refpts(eql))[3])
  params(eql2)[]=apply(rec.obs(eql2),6,mean)
  refpts(eql2)=refpts(eql2)[c("f0.1","fmax","spr.30","virgin"),]
  dimnames(refpts(eql2))[[1]][4]="spr.100"

  eql2=brp(eql2)

  res =model.frame(refpts(eql )[c("msy","crash","virgin"), c("harvest","yield","rec","ssb","biomass")])
  res1=model.frame(refpts(eql1)[c("f0.1","fmax","spr.30","spr.100"),c("harvest","yield","rec","ssb","biomass")])
  res2=model.frame(refpts(eql2)[c("f0.1","fmax","spr.30","spr.100"),c("harvest","yield","rec","ssb","biomass")])
  res=cbind(res[,4:5],res[1:3],res1[,1:4],res2[,1:4])
  
  names(res)[rev(length(names(res))-0:3)]=paste(names(res)[rev(length(names(res))-0:3)],"_",sep="")
  res=transform(res,iter=as.numeric(iter))
  names(res)[names(res)=="quant"]="quantity"
  
  if (dims(refpts(eql))$iter>1){
    catch(object)=propagate(catch(object),dims(refpts(eql))$iter)
    dimnames(catch(object))$iter=dimnames(rec(object))$iter
  }
  
  current=as.data.frame(mcf(FLQuants(object,"harvest"=FLCore::fbar,
                                            "yield"  =catch,
                                            "rec"    =rec, 
                                            "ssb"    =ssb,
                                            "biomass"=computeStock)),drop=TRUE)[,-1]
  
  current=subset(current,year==max(year)&!is.na(data))[,-1]
  names(current)[names(current)=="data"]="current"
  names(current)[names(current)=="qname"]="quantity"
  if (!("iter"%in%names(current))) 
    current=cbind(iter=1,current)
  current=transform(current,iter=as.numeric(as.factor(iter)))
  
  res=merge(res,current,by=c("iter","quantity"))
  
  #Warning message:
  #  In .nextMethod(x = x, MARGIN = MARGIN, STATS = STATS, FUN = FUN) :
  #  length(STATS) or dim(STATS) do not match dim(x)[MARGIN]
  #r =log(aaply(leslie(eql, f=c(refpts(eql)["crash","harvest"])),3,
  #             function(x) lambda(x[drop=T])))
  r =log(maply(seq(dims(eql)$iter), 
               function(i) {res=leslie(iter(eql,i), f=c(refpts(iter(eql,i))["crash","harvest"]))
                            lambda(res[drop=TRUE])}))
  
  names(r)=dimnames(refpts(eql))$iter

  rc=log(maply(seq(dims(eql)$iter), 
               function(i) {res=leslie(iter(eql,i), f=c(refpts(iter(eql,i))["msy","harvest"]))
               lambda(res[drop=TRUE])}))
  names(rc)=dimnames(refpts(eql))$iter
  
  r =data.frame(r=r, rc=rc, iter=seq(dim(stock.n(object))[6]))
  
  res=cbind(iter=subset(res,quantity=="ssb")[ ,  1 ],
            b   =subset(res,quantity=="biomass")[,-(1:2)],
            s   =subset(res,quantity=="ssb"    )[,-(1:2)],
            r   =subset(res,quantity=="rec"    )[,-(1:2)],
            f   =subset(res,quantity=="harvest")[,-(1:2)],
            y   =subset(res,quantity=="yield"  )[,-(1:2)])
  
  res=transform(merge(res,r,by="iter"),
                rt=log(b.msy/b.current)/r)
  
  options(warn=warn)
  
  res=res[, !(names(res)%in%c("b.crash","s.crash","r.crash","y.crash",
                              "b.year", "s.year", "r.year", "y.year","f.year",
                              "f.virgin","f.spr.100","f.spr.100_","y.virgin","y.spr.100"))]
  
  res[,"r.current"]=c(params(eql1))
  mf2FLPar(res)}

### r conditional on maturity so comparable across stocks
new<-function(object,s=NULL,recent=3){

  if (!requireNamespace("popbio", quietly = TRUE)) {
    stop("Package \"popbio\" needed for this function to work. Please install it.",
         call. = FALSE)}
  
  warn=options()$warn
  options(warn=-1)
  
  #SRR
  if (!is.null(s)){
    spr0=mean(spr0(FLBRP(object)))
    srr =as.FLSR(object,model="bevholtSV")
    
    upper(srr)[1:2]=1e12
    srr=fmle(srr,
             fixed=list(s=s,spr0=spr0),
             control=list(silent=TRUE),
             method="Brent")
    eql=FLBRP(object,sr=list(model=bevholt()$model,params=ab(params(srr),"bevholt")[c("a","b")]))
  }else{ 
    srr=as.FLSR(object,model="bevholt")
    lower(srr)[1:2]=0.00001
    srr=fmle(srr,control=list(trace=FALSE),method="L-BFGS-B")
    eql=FLBRP(object,sr=srr)
    }
  
  dimnames(refpts(eql))$refpt[7]="spr.100"

  ## recruits
  shift=unlist(c(ddply(rod(rec.obs(eql)),.(iter), function(x)
    mean(subset(x,regime==max(as.numeric(regime)))$data))["V1"]))
  shift=FLPar(rshift=array(shift,c(1,length(shift))))
  mn=apply(rec.obs(eql),6,mean)
  mn=FLPar(rmean=array(mn,c(1,length(mn))))
  rtn=rbind(mn,shift)

  ## latest
  res=FLQuants(object[,dim(object)[2]],
               "f"      =FLCore::fbar,
               "yield"  =catch,
               "rec"    =rec, 
               "ssb"    =ssb,
               "biomass"=computeStock)
  late=FLPar(laply(res,FLPar))
  names(late)=c("params","iter")
  dimnames(late)$params=c("f","yield","rec","ssb","stock")
  rtn=rbind(rtn,late)
  dimnames(rtn)$params[3:7]=paste(dimnames(rtn)$params[3:7],"latest",sep=".")
  
  ## recent
  res=FLQuants(object[,dim(object)[2]-seq(recent)+1],
               "f"      =FLCore::fbar,
               "yield"  =catch,
               "rec"    =rec, 
               "ssb"    =ssb,
               "biomass"=computeStock)
  rcnt=laply(res,function(x) FLPar(array(aaply(x,6,mean),c(1,dim(x)[6]))))
  rcnt=FLPar(rcnt)
  names(rcnt)=c("params","iter")
  dimnames(rcnt)$params=c("f","yield","rec","ssb","stock")
  rtn=rbind(rtn,rcnt)
  dimnames(rtn)$params[8:12]=paste(dimnames(rtn)$params[8:12],"recent",sep=".")
  
  ## growth rate
  r =log(maply(seq(dims(eql)$iter), 
               function(i) {res=leslie(iter(eql,i), f=c(refpts(eql)["crash","harvest",i]))
               lambda(res[drop=TRUE])}))
  r =FLPar("r"=array(r,c(1,length(r))))
  rtn=rbind(rtn,r)
  
  rc=log(maply(seq(dims(eql)$iter), 
               function(i) {res=leslie(iter(eql,i), f=c(refpts(eql)["msy","harvest",i]))
               lambda(res[drop=TRUE])}))
  rc=FLPar("rc"=array(rc,c(1,length(rc))))
  rtn=rbind(rtn,rc)
  
  spr=FLPar(spr0=array(c(spr0(eql)),c(1,dim(rtn)[2])))
  rtn=rbind(rtn,spr)
  
  rfs=computeRefpts(eql)[,c("harvest","yield","rec","ssb","biomass")]
  names(rfs)[2]="params"
  dimnames(rfs)$params[1]="f"
  
  nms=c("yield","rec","ssb","biomass")
  rfs["f0.1",   nms]=sweep(rfs["f0.1",   nms],c(1,3),rfs["f0.1",   "rec"],"/")
  rfs["fmax",   nms]=sweep(rfs["fmax",   nms],c(1,3),rfs["fmax",   "rec"],"/")
  rfs["spr.30", nms]=sweep(rfs["spr.30", nms],c(1,3),rfs["spr.30", "rec"],"/")
  rfs["spr.100",nms]=sweep(rfs["spr.100",nms],c(1,3),rfs["spr.100","rec"],"/")
  
  rtn=rbind(rtn,FLPar(rfs["msy",                            drop=TRUE]))
  dimnames(rtn)[[1]][16:20]=paste(dimnames(rtn)[[1]][16:20],"msy",sep=".")
  
  rtn=rbind(rtn,FLPar(crash=array(rfs["crash","f",drop=TRUE],c(1,dim(rtn)[2]))))

  rtn=rbind(rtn,FLPar(rfs["virgin",c("ssb","biomass","rec"),drop=TRUE]))
  dimnames(rtn)[[1]][22:24]=paste(dimnames(rtn)[[1]][22:24],"virgin",sep=".")
  
  rtn=rbind(rtn,FLPar(rfs["f0.1",                           drop=TRUE]))
  dimnames(rtn)[[1]][25:29]=paste(dimnames(rtn)[[1]][25:29],"f0.1",sep=".")
  
  rtn=rbind(rtn,FLPar(rfs["fmax",                           drop=TRUE]))
  dimnames(rtn)[[1]][30:34]=paste(dimnames(rtn)[[1]][30:34],"fmax",sep=".")
  
  rtn=rbind(rtn,FLPar(rfs["spr.30",                         drop=TRUE]))
  dimnames(rtn)[[1]][35:39]=paste(dimnames(rtn)[[1]][35:39],"spr.30",sep=".")
  
  rtn=rbind(rtn,FLPar(rfs["spr.100",                        drop=TRUE]))
  dimnames(rtn)[[1]][40:44]=paste(dimnames(rtn)[[1]][40:44],"spr.100",sep=".")
  
  ssb.kobe    =rtn["ssb.latest"]/rtn["ssb.msy"]
  f.kobe      =rtn["f.latest"]/rtn["f.msy"]
  biomass.kobe=rtn["stock.latest"]/rtn["biomass.msy"]
  ssb.pri     =rtn["ssb.latest"]/rtn["ssb.virgin"]

  dimnames(ssb.kobe)$params="ssb.kobe"
  dimnames(f.kobe)$params="f.kobe"
  dimnames(biomass.kobe)$params="biomass.kobe"
  dimnames(ssb.pri)$params="ssb.pri"
  
  rt          =log(rtn["ssb.msy"]/rtn["ssb.latest"])/rtn["rc"]
  dimnames(rt)$params="rt"
  
  rtn=rbind(rtn,rt)
  rtn=rbind(rtn,ssb.kobe)
  rtn=rbind(rtn,f.kobe)
  rtn=rbind(rtn,biomass.kobe)
  rtn=rbind(rtn,ssb.pri)

  sr=params(eql)
  dimnames(sr)$params
  dimnames(sr)$params=c("alpha","beta")

  rtn=rbind(rtn,sr)
  
  x=rtn[c("alpha","beta","spr0")]
  dimnames(x)$params=c("a","b","spr0")
  rtn=rbind(rtn,sv(x,"bevholt")[1:2])

  shape=rtn["ssb.msy"]/rtn["ssb.virgin"]
  dimnames(shape)$params="shape"
  rtn=rbind(rtn,shape)
  
  prod =rtn["yield.msy"]/rtn["ssb.msy"]
  dimnames(prod)$params="prod"
  rtn=rbind(rtn,prod)
  
  rtn=rtn[!dimnames(rtn)$params%in%c("rec.f0.1","rec.fmax","rec.spr.30","rec.spr.100")]

  options(warn=warn)

  rtn}

setMethod("refs", signature(object="FLStock"), function(object,s=NULL) new(object,s))
          