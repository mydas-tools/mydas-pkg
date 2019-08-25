## Observation Error Model
globalVariables("ctrl")
globalVariables("prrs")
globalVariables("cpue")
globalVariables("phaseQ")
globalVariables("bounds")
globalVariables("uCV")

#' setGeneric('oem',    function(object,...) standardGeneric('oem'))
#' 
#' #' @title oem
#' #'
#' #' @description Creates an \code{FLQuant} to represent an index of relative abundance
#' #' 
#' #' @param object \code{FLStock} for which an index is to be generated
#' #' @param sel \code{FLQuant} vector at age that shapes the catch or biomass.
#' #' additive, by default set to \code{TRUE},  
#' #' @param fish.depenpent \code{logical} that determines whether the index is proportional to the 
#' #' stock or catch/effort, by default set to \code{TRUE}
#' #' @param effort \code{character} c("h","f") that determines what proxy to use for effort,
#' #' i.e. harvest rate or instanteous fishing mortality 
#' #' @param mass \code{logical} default is TRUE, whether index is in mass or numbers
#' #' @param ... other arguments
#' #' 
#' #' @aliases 
#' #' oem,FLStock-method 
#' #' 
#' #' @export
#' #' @rdname oem
#' #' 
#' #' @return \code{FLQuant} with time series
#' #' 
#' #' @aliases oem-method oem,FLStock,ANY-method oem 
#' #' 
#' #' @export
#' #' 
#' #' @examples
#' #' \dontrun{
#' #'  data(ple4) 
#' #'  cpue=oem(ple4) 
#' #'  }
#' setMethod('oem',   signature(object='FLStock'),
#'            function(object,
#'                     sel       =FLQuant(FLQuant(1,dimnames=dimnames(harvest(object)))),
#'                     timing    =0.5,
#'                     fish.dependent=TRUE,
#'                     effort    =c("f","h"),
#'                     mass      =TRUE){
#' 
#'   timing=pmax(pmin(timing,1.0),0.0)
#'   stock=(stock(object)[,-dim(stock(object))[2]]*timing+stock(object)[,-1]*(1.0-timing))
#'   
#'   object=window(object,start=range(object)["minyear"],end=range(object)["maxyear"])
#'   
#'   if (fish.dependent) {
#'     if (effort[1]=="h")
#'       E=catch(object)%/%stock
#'     else  
#'       E=fbar(object)
#'     
#'     cpue=(catch.n(object)%*%sel)%/%E
#'   }else 
#'       cpue=apply((stock%*%sel),2:6,sum)
#'   
#'   if (mass)
#'     cpue=cpue%*%catch.wt(object)
#' 
#'   cpue})

# #setMethod('survey', signature(object='FLStock'),
# survey=function(object,timing=FLQuant(0,dimnames=dimnames(m(object))),wt=stock.wt(object),sd=0,log=FALSE,...){
#   
#   if(units(harvest(object)) == 'f'){
#     res      <- colSums(stock.n(object)*exp(-harvest(object)*timing - m(object)*timing)*wt, na.rm=FALSE)
#     dim(res) <- c(1, dim(res))
#     dmns     <-dimnames(stock(object))
#     dmns$iter<-dimnames(res)$iter
#     
#     res=FLQuant(res, dimnames=dmns)
#   } else if(units(harvest(object)) == 'hr') {
#     
#     res      = colSums(object@stock.n*(1-object@harvest*timing)*exp(-object@m*timing)*timing*wt)
#     dim(res) = c(1, dim(res))
#     
#     res=FLQuant(res, dimnames=dimnames(object@stock))
#   } else
#     stop('Correct units (f or hr) not specified in the harvest slot')
#   
#   if (sd>0){
#     if (!log) res=apply(res,2:6,function(x,sd) rnorm( 1,x,sd=sd),sd=sd) else 
#       res=apply(res,2:6,function(x,sd) rlnorm(1,x,sdlog=sd),sd=sd)
#   }
# 
#    res}
# 
# #setMethod('cpue', signature(object='FLStock'),
# u=function(object,partialf=FLQuant(1,dimnames=dimnames(m(object))),wt=catch.wt(object),sd=0,log=FALSE,...){
#   
#   if(units(harvest(object)) == 'f'){
#     res      <- colSums(stock.n(object)*harvest(object)*partialf/(harvest(object)+m(object))*exp(1-harvest(object)-m(object)), na.rm=FALSE)
#     dim(res) <- c(1, dim(res))
#     dmns     <-dimnames(stock(object))
#     dmns$iter<-dimnames(res)$iter
#     
#     res=FLQuant(res, dimnames=dmns)
#   } else if(units(harvest(object)) == 'hr') {
#     
#     res      = colSums(object@stock.n*(1-object@harvest*partialf)*exp(-object@m)*wt)
#     dim(res) = c(1, dim(res))
#     
#     res=FLQuant(res, dimnames=dimnames(object@stock))
#   } else
#     stop('Correct units (f or hr) not specified in the harvest slot')
#   
#   if (sd>0){
#     if (!log) res=apply(res,2:6,function(x,sd) rnorm( 1,x,sd=sd),sd=sd) else 
#       res=apply(res,2:6,function(x,sd) rlnorm(1,x,sdlog=sd),sd=sd)}
#   
#   res} 
# 
# cpueBiodym2Aspic=function(bd,type='B0',n=1,sd=0,log=TRUE){
#   
#   type=toupper(type)
#   effort=harvest(bd)
#   if (n>1)  effort=propagate(effort,n)
#   
#   dmns=dimnames(stock(bd))
#   dmns$iter=1              
#   dev=if (log) rlnorm(n,FLQuant(0,dimnames=dmns),sd) else rnorm( n,FLQuant(0,dimnames=dmns),sd)
#   
#   ctc=stock(bd)[,-dims(bd)$year]*effort
#   
#   switch(type,
#          
#          ## Error on Catch
#          CC={ res=cbind(name='CC sim',model.frame(mcf(FLQuants(effort=effort,catch=ctc,dev=dev)),drop=T, stringsAsFactors=FALSE))
#               res=transform(res,catch=catch*dev,index=catch*dev/effort)
#               res},
#          
#          ## Error on index
#          CE={ res=cbind(name='CE sim',model.frame(mcf(FLQuants(effort=effort,catch=ctc,dev=dev)),drop=T, stringsAsFactors=FALSE))
#               res=transform(res,index=catch*dev/effort)
#          },   
#          
#          ## Error on stock
#          B0=cbind(name='B0 sim',model.frame(mcf(FLQuants(stock=stock(bd)*dev)),drop=T, stringsAsFactors=FALSE)),
#          B1=cbind(name='I1 sim',model.frame(mcf(FLQuants(stock=(stock(bd)[,-dim(stock(bd))[2]]+stock(bd)[,-1]/2),dev=dev)),drop=T, stringsAsFactors=FALSE)),
#          B2={ res=cbind(name='I0 sim',model.frame(mcf(FLQuants(stock=stock(bd),dev=dev)),drop=T, stringsAsFactors=FALSE))
#               res},          
#          
#          ## Error on stock
#          I0=cbind(name='B0 sim',model.frame(mcf(FLQuants(index=.1*stock(bd)*dev)),drop=T, stringsAsFactors=FALSE)),
#          I1=cbind(name='I1 sim',model.frame(mcf(FLQuants(index=.1*(stock(bd)[,-dim(stock(bd))[2]]+stock(bd)[,-1]/2),dev=dev)),drop=T, stringsAsFactors=FALSE)),
#          I2={ res=cbind(name='I0 sim',model.frame(mcf(FLQuants(index=.1*stock(bd),dev=dev)),drop=T, stringsAsFactors=FALSE))
#               res}                
#   )}
# 
# oemFn<-function(flt,om,paa=paa,dev=uCV){
#   
#   res=apply(sweep(catch.n(om),1,paa[,flt],"*")*catch.wt(om),c(2,6),sum)
#   
#   if (!is.FLQuant(dev)){
#     
#     dev=rlnorm(dims(res)$iter,FLQuant(0,dimnames=dimnames(iter(res,1))),dev[flt])}
#   
#   res*dev[,dimnames(res)$year]}
# 
# 
# oemOld=function(om,cv,trendQ=FLQuant(1,dimnames=dimnames(stock(om))),
#               omega =1,refB=1,
#               fishDepend=FALSE){
#   
#   nits=max(dims(stock(om))$iter,dims(catch(om))$iter)
#   rnd=rlnorm(nits,FLQuant(0,dimnames=list(year=dims(om)$minyear:dims(om)$maxyear)),cv)
#   
#   if (fishDepend)
#     cpue=catch(om)/fbar(om)
#   else
#     cpue=computeStock(om)
#   
#   cpue=trendQ[,dimnames(cpue)$year]*cpue
#   cpue=cpue*(stock(om)%/%refB)^(omega)
#   cpue=rnd*cpue
#   
#   cpue}


# cpue
# /AvCatch(iAge,iYear,pTune->GetStartFishing(iFleet), pTune->GetEndFishing(iFleet));
# 
# AvCatch(iYear, double StartFishing, double EndFishing)
# {
#   //Equation (2) VPA User Guide
#   double Z;
#   
#   Z = (GetM(iAge, iYear)+F[iAge][iYear]);                   
#   
#   if (StartFishing == 0.0 && EndFishing == 0.0)
#     return 1.0;
#   else if (StartFishing == EndFishing)
#     EndFishing += 0.00001;
#   
#   return (exp(-StartFishing*Z) - exp(-EndFishing*Z))/((EndFishing-StartFishing)*Z); 
# }

# setMethod('oem<-', signature(object='biodyn',value="FLStock"), 
#           function(object,value, 
#                    sel           =FLQuant(FLQuant(1,dimnames=dimnames(harvest(value)))),
#                    timing        =0.5,
#                    fish.dependent=TRUE,
#                    effort        =c("f","h"),
#                    mass          =TRUE){
#           res=oem(value,sel,timing,fish.dependent,effort,mass)  
#           res=apply(res,2:6,sum)
#           object@indices=FLQuants(res)
#           object})



#' @title oem
#'
#' @description Creates an \code{FLQuant} to represent an index of relative abundance
#' 
#' @param object \code{FLStock} for which an index is to be generated
#' @param sel \code{FLQuant} vector at age that shapes the catch or biomass.
#' additive, by default set to \code{TRUE},  
#' @param fish.depenpent \code{logical} that determines whether the index is proportional to the 
#' stock or catch/effort, by default set to \code{TRUE}
#' @param effort \code{character} c("h","f") that determines what proxy to use for effort,
#' i.e. harvest rate or instanteous fishing mortality 
#' @param mass \code{logical} default is TRUE, whether index is in mass or numbers
#' @param ... other arguments
#' 
#' @aliases 
#' oem,FLStock-method 
#' 
#' @export
#' @rdname oem
#' 
#' @return \code{FLQuant} with time series
#' 
#' @aliases oem-method oem,FLStock,ANY-method oem 
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#'  data(ple4) 
#'  cpue <- oem(ple4) 
#'  }
setMethod('oem',   signature(object='FLStock'),
          function(object, sel=FLQuant(FLQuant(1, dimnames=dimnames(harvest(object)))),
                   timing = 0.5, fish.dependent = TRUE, effort = c("f","h"), mass = TRUE) {
            
            timing <- pmax(pmin(timing, 1.0), 0.0)
            stock <- (stock(object)[, -dim(stock(object))[2]] * timing +
                        stock(object)[, -1]*(1.0 - timing))
            
            if (fish.dependent) {
              if (effort[1] == "h")
                E <- catch(object) %/% stock
              else  
                E <- fbar(object)
              
              cpue <- (catch.n(object) %*% sel) %/% E
            } else { 
              cpue <- apply((stock%*%sel), 2:6, sum)
            }
            if (mass)
              cpue <- cpue %*% catch.wt(object)
            
            return(cpue)
          })

hyperstability<-function(object,omega=1,ref=apply(object,c(1,3:6),mean)) 
  ref%*%((object%/%ref)^omega)

bias<-function(object,bias=0.02) 
  FLQuant(cumprod(1+rep(bias,dim(object)[2])),dimnames=dimnames(object))




