utils::globalVariables(c("prior","ddply",".", "hat","aply","block","se","ML","mlzFn","daply"))

#' mlz
#' 
#' @title mlz 
#' 
#' @description A wrapper function for the MLZ package https://cran.r-project.org/web/packages/MLZ/index.html
#' to allow i be simulation tested.
#' 
#' @author Laurence Kell, Sea++
#'  
#' @param year blah,blah,blah,...
#' @param len blah,blah,blah,...
#' @param lc blah,blah,blah,...
#' @param linf v blah,blah,blah,...
#' @param k blah,blah,blah,...
#' @param ncp blah,blah,blah,...
#' @param ss blah,blah,blah,...
#' 
#' @export
#' @docType methods
#' 
#' @rdname mlz
#' 
#' @examples
#' \dontrun{
#' res=mydas:::mlz(mnLen[,ac(40:60)],prior)
#' }

mlzFn<-function(year,len,lc,linf,k,ncp=2,ss=500){
  
  dat=new("MLZ_data", 
          Year         =unique(year), 
          MeanLength   =len,
          ss           =rep(ss,length(year)),
          length.units ="cm")
  
  dat@Lc    =lc
  dat@vbLinf=linf
  dat@vbK   =k
  
  hat=ML(dat,ncp=ncp,figure=FALSE)
  
  res=hat@estimates
  res}

mlz<-function(object,params,ncp=2,ss=500){
  
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
         ncp=ncp,
         ss=ss)
    
    rtn=rbind(rtn,cbind(iter=i,res))}
  
  return(rtn)
  
  dimnames(rtn)[[2]][2:3]=c("hat","se")
  
  z=data.frame(subset(rtn,substr(dimnames(rtn)[[1]],1,1)=="Z"))
  z=ddply(z,.(iter), transform, block=seq(length(hat)))
  z=as(daply(z,.(block), with, as(data.frame(hat=hat,se=se,iter=iter),"FLPar")),"FLPar")
  
  y=subset(rtn,substr(dimnames(rtn)[[1]],1,1)!="Z")
  y=data.frame(y,name=substr(dimnames(y)[[1]],1,1))
  y=ddply(y,.(iter,name), transform, block=seq(length(hat)))
  y=as(daply(y,.(block,name), with, as(data.frame(hat=hat,se=se,iter=iter),"FLPar")),"FLPar")
  
  res=list("z"=z,"year"=y)
  
  return(res)
  
  names(res)=c("z","year")
  
  res}  

mlz.z<-function(object,data){
  year=as.numeric(dimnames(data)$year)
  brk =round(object[grep("year",dimnames(object)[[1]]),2])
  year=c(min(year),brk,max(year))
  
  rtn =mdply(data.frame(block=seq(length(year)-1)),function(block)
    data.frame(year=year[block]:year[block+1]))
  
  z   =object[substring(dimnames(object)[[1]],1,1)=="Z",2:3]
  dimnames(z)[[2]]=c("z","se")
  
  z   =data.frame(block=seq(dim(z)[1]),z)
  merge(rtn,z)[,c(2,1,3:4)]}

