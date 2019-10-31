
roc<-function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), 
             FPR=cumsum(!labels)/sum(!labels),
             labels,
             reference=sort(scores))}


