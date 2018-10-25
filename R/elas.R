
ProdFn<-function(x){
  res=lhPar(FLPar(x))
  eql=lhEql(res,m=function(wt,params){
    res=wt%^%params["m2"] 
    res=params["m1"]%*%res
    res})
  
  c(log(lambda(leslie(eql,fbar=c(refpts(eql)["crash","harvest"]))[drop=TRUE])),
    log(lambda(leslie(eql,fbar=c(refpts(eql)["msy",  "harvest"]))[drop=TRUE])))}

rFn=function(x){
  eql=lhEql(lhPar(FLPar(x)),m=function(wt,params){
    res=wt%^%params["m2"] 
    res=params["m1"]%*%res
    res})
  
  c("r"  =log(lambda(leslie(eql,fbar=c(refpts(eql)["crash","harvest"]))[drop=TRUE])),
    "r.c"=log(lambda(leslie(eql,fbar=c(refpts(eql)["msy","harvest"]))[drop=TRUE])))}

gtFn=function(x){   
  eql=lhEql(lhPar(FLPar(x)),m=function(wt,params){
    res=wt%^%params["m2"] 
    res=params["m1"]%*%res
    res})
  
  c("gt"    =genTime(eql),
    "age"   =aMsy(eql),
    "lopt"  =lopt(x,m=function(l,params){
      wt=params["a"]*(l^params["b"])
      res=exp(log(wt)*params["m2"]) 
      res=params["m1"]%*%res
      res
    }))}

refFn<-function(x){
  
  rf=refpts(lhEql(lhPar(FLPar(x)),m=function(wt,params){
    res=wt%^%params["m2"] 
    res=params["m1"]%*%res
    res}))
  
  res=c(hmsy  =c(rf["msy","yield"]/rf["msy","ssb"]),
        fcrash=c(rf["crash","harvest"]),
        fmsy  =c(rf["msy","harvest"]),
        bmsy  =c(rf["msy","ssb"]),
        msy   =c(rf["msy","yield"]),
        shape =c(rf["msy","ssb"]/rf["virgin","ssb"]))
  
  res}
