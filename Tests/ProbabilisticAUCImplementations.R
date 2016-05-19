
binaryclass.Probabilisticauc = function(pred, truth) {
  ptm <- proc.time()
  # Parameters, names=FALSE
  K=dim(pred)[2]
  n=dim(pred)[1]
  levels=levels(truth)
  print(paste("levels",levels))
  
  
  # Design Matrix f
  f= model.matrix( ~ 0 + truth, truth)
  f.summed=apply(f,2,sum)
  
  # matrix AUC
  probAUC=matrix(rep(NA, K*K),K,K)
  
  for (j in c(1:K)) {
    for (k in c(1:K)) {
      if (k!=j) {
        probAUC[k,j]=((t(f[,j])%*%pred[,j])/f.summed[j]-(t(f[,k])%*%pred[,j])/f.summed[k]+1)/2
      }
    }
  }
    colnames(probAUC)<-levels
    rownames(probAUC)<-levels

  print(proc.time()-ptm)
  return(probAUC)
}



binaryclass.Probabilisticauc.permutations = function(pred, truth) {
  ptm <- proc.time()
  # Parameters, names=FALSE
  K=dim(pred)[2]
  n=dim(pred)[1]
  levels=levels(truth)
  print(paste("levels",levels))
  
  
  # Design Matrix f
  f= model.matrix( ~ 0 + truth, truth)
  f.summed=apply(f,2,sum)
  
  # matrix AUC
  probAUC=matrix(rep(NA, K*K),K,K)
  
  for (j in c(1:K)) {
    for (k in c(1:K)) {
      if (k!=j) {
        probAUC[k,j]=((t(f[,j])%*%pred[,j])/f.summed[j]-(t(f[,k])%*%pred[,j])/f.summed[k]+1)/2
      }
    }
  }
  colnames(probAUC)<-levels
  rownames(probAUC)<-levels
  
  print(proc.time()-ptm)
  return(probAUC)
}