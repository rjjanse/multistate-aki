#Estimates of PDI and its components
pdiest<-function(data){
  
  y<-data$outcome
  ymin<-min(y)
  ymax<-max(y)
  noutcome<-ymax-ymin
  p<-prod(table(y))
  pdi<-c()
  
  for (i in 1:(noutcome+1)){
    
    predprob<-data[,(i+1)]  #READ predicted probabilities for level i
    t0<-table(predprob,y)   #CALCULATE frequencies of predicted probabilities for level i by outcome
    
    dim1<-dim(t0)[1]
    dim2<-dim(t0)[2]
    t<-cbind(t0[,i],t0[,-i]) #REORDER columns
    restrictt<- if (noutcome == 1){matrix(t[,2:(noutcome+1)],ncol=1)} else {t[,2:(noutcome+1)] } #REMOVE first column of t
    
    c<-apply(restrictt,2,cumsum) #CALCULATE cumulative frequencies of predicted probabilities for level i by outcome
    cnew<- if (noutcome == 1) {rbind(rep(0,noutcome),matrix(c[1:(dim(c)[1]-1),],ncol=))} else {rbind(rep(0,noutcome),c[1:(dim(c)[1]-1),])} #INTRODUCE a row of zeros at the begining of c
    
    mat<-c()                     #MATRIX of 0s and 1s of dimension 2^(noutcome) x noutcome
    for (j in 1:noutcome){
      mat0<-cbind(mat,0)
      mat1<-cbind(mat,1)
      mat<-rbind(mat0,mat1)}
    
    r<-0
    for (k in 1:dim(mat)[1]){
      dt<-t(apply(restrictt, 1, function(x) mat[k,]*x))
      dcnew<-t(apply(cnew, 1, function(x) (1-mat[k,])*x))
      dfinal<-if (noutcome == 1) {cbind(t[,1],t(dt+dcnew))} else {cbind(t[,1],dt+dcnew)} #TAKE all combinations of frequencies and cumulative frequencies
      r<-r+sum(apply(dfinal,1,prod))/(1+sum(mat[k,]))}                                   #MULTIPLYIES across rows
    
    r<-r/p     #PDI component for outcome i
    pdi<-rbind(pdi,r)
  }
  pdi<-rbind(mean(pdi),pdi)
  pdi
}
