#' hcrSBTD
#' 
#' @title hcrSBTD 
#' 
#' @description An emprirical Harvest Control Rule that sets a quota based on a trend in an index
#' 
#' @author Laurence Kell, Sea++
#'  
#' @anme hcrSBTD  
#' 
#' @param yrs  blah,blah,blah,...
#' @param control  blah,blah,blah,...FLPar(c(k1=0.25,k2=0.25,gamma=1))
#' @param index  blah,blah,blah,...
#' @param catch  blah,blah,blah,...
#' @param ... any additional arguments
#' 
#' @export hcrSBTD
#' @docType methods
#' 
#' @rdname hcrSBTD
#' 
#' @examples
#' \dontrun{
#' data(pl4)
#' }
#' @export hcrSBTD
#' 
hcrSBTD<-function(yrs,
                  control=FLPar(c(k1=0.25,k2=0.25,gamma=1)),
                  index,
                  catch,...){
  
  lambda=as.FLQuant(ddply(transform(as.data.frame(index%/%apply(index,6,mean)),iter=as.numeric(ac(iter))), 
                          .(iter), with, data.frame(data=coefficients(lm(data~year))[2])))

  flag =c(lambda<0)
  flag1=seq(dims(lambda)$iter)[ flag]
  flag2=seq(dims(lambda)$iter)[!flag]
  lambda =abs(lambda)
  control=propagate(FLPar(control),dims(lambda)$iter)

  gain=control["k1"]
  if (length(flag1)>0)
   gain[1,flag1]=-(control["k1",flag1])%*%(lambda[,,,,,flag1]%^%control["gamma",flag1])

  if (length(flag2)>0)
    gain[1,flag2]=  control["k2",flag2] %*% lambda[,,,,,flag2]
  #res   =1+ifelse(flag,-control["k1"],control["k2"])*exp(log(lambda)*ifelse(flag,control["gamma"],1))

  res   =(1+gain)%*%apply(catch,6,mean)
  
  dmns     =dimnames(catch)
  dmns$year=yrs
  dmns$iter=dimnames(index)$iter
  res      =FLQuant(rep(c(res),each=length(yrs)),dimnames=dmns)
  
  #chk=list(lambda=lambda,gain=gain,catch=catch,index=index,res=res) 
 
  return(res)}
