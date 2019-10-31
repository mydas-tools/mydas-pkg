pmTimeseriesFn<-function(om,eql){
  
  res=model.frame(FLQuants(
    #stock status
    ssb=ssb( om)%/%refpts(eql)["msy",   "ssb"],
    f  =fbar(om)%/%refpts(eql)["msy",   "harvest"],
    #saftety
    rec=rec( om)%/%refpts(eql)["virgin","rec"],
    #yield
    yld= catch( om)%/%refpts(eql)["msy","yield"]),drop=TRUE)
  
  res}
