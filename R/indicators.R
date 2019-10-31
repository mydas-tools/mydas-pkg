#' @export indicators
#' 
indicators<-function(params,
                     m   ="gislason",
                     f    =1,
                     srDev=FLQuant(rep(1,121)),
                     fbar  =srDev%=%1){
  
  ##set up equilibrium object
  if ("numeric"%in%is(f))
    f=FLPar(f=array(f,c(1,length(f))))
   
  ## need to add interactions for f and par
  if (dim(params)[2]>1&dim(f)[2]>1){
    npar=dim(params)[2]
    params=as(mdply(seq(dim(f)[2]), with, 
                    cbind(model.frame(params)))[,-c(1,dim(params)[1]+2)],"FLPar")
    f  =rep(c(f),each=npar)
    f  =FLPar(array(f,c(1,length(f))))
  }
  
  eql=lhEql(params,m=m)
  
  ## convert to FLStock with constant F
  eq=eql
  
  fbar(eq)=fbar
  mlt=FLPar(f=array(c(f)*c(eq@refpts["msy","harvest"]),
                    c(1,length(c(f)*c(eq@refpts["msy","harvest"])))))
  fbar(eq)=fbar(eq)%*%mlt

  stk=fwd(propagate(as(eq,"FLStock"),dim(srDev)[6]),fbar=fbar(eq)[,-1],
          sr=eq,residuals=srDev)
  
  ## Other stuff
  #srr=model.frame(FLQuants(eql,"ssb"=ssb,"rec"=rec),drop=TRUE)
  #srp=model.frame(FLQuants("rp"=setPlusGroup(stock.n(eq)[,1]*stock.wt(eq)[,1]*eq@mat[,1],12)),drop=T)
  #ts =model.frame(FLQuants(stk,"ssb"=ssb,"biomass"=stock,"rec"=rec,"catch"=catch,
  #                             "dev"=function(x) propagate(srDev,dim(x)[6])),drop=T)
  
  ## Summary stats
  ind=mydas:::omSmry(stk,eql,params)
  
  refs=model.frame(popdyn(params,eq=lhEql(params,m=m)))
  
  key=cbind(model.frame(params),f=c(f))
  
  list(ind=ind,refs=refs,key=key,
       ctn=catch.n(stk),
       cln=exp(log(catch.wt(stk)%/%params["a"])%/%params["b"]))
}




#' $L_{max5\%}$ Mean length of largest 5\%
#' $L_{95\%}$ $95^{th}$ 95thpercentile
#' $P_{mega}$ Proportion of individuals above $L_{opt} + 10\%$
#' $L_{25\%}$ $25^{th}$ percentile of length distribution
#' $L_{c}$ Length at $50\%$ of modal abundance
#' $L_{mean}$ Mean length of individuals $> L_c$
#' $L_{max_{y}}$ Length class with maximum biomass in catch
#' $L_{mean}$ Meanlength of individuals $> L$
    