lfm<-function(params,m,lc){
   (params["linf"]%+%((2*m)%/%params["k"])%*%lc)%/%(1+2*(m%/%params["k"]))}


