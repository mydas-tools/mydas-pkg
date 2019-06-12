#' gt
#' 
#' @title gt 
#' 
#' @description 
#' @author Laurence Kell, Sea++
#'  
#' @name gt
#' @param object \code{FLBRP}
#' 
#' @aliases
#' 
#' @export gt
#' @docType method
#' 
#' @rdname gt
#' @seealso 
#' 
#' @examples
#' \dontrun{
#' data(pl4)
#' alk=gt(FLBRP(ple4))
#' }
fitSvr=function(dat, target, predictors, tuneLength=5){
  
  ctrl=caret::trainControl(
    method ="repeatedcv",
    repeats=5)
  
  print(paste(target, "~", predictors, sep=" "))
  
  svrFit = caret::train(
    as.formula(paste(target, "~", predictors, sep=" ")),
    data      =dat,
    method    ='svmRadial',
    tuneLength=tuneLength,
    trControl =ctrl,
    metric    ="MAE")
  
  return(svrFit)}

svrFn=function(ind,svr1,svr2){
  ind         =t(data.frame(ind))
  prediction1 =predict(svr1,ind)
  prediction2 =predict(svr2,ind)
  result      =c(prediction1, prediction2)
  
  return(result)}

svrFn2=function(ind,...)
  laply(list(...), function(x) predict(x,t(data.frame(ind))))
