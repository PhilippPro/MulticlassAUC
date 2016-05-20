
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



binaryclass.Probabilisticauc.conventions<-function (X, y, plotROC = FALSE, alg = c("Wilcoxon", "ROC")) {
y = as.factor(y)
X = as.matrix(X)
alg = match.arg(alg)
if (nrow(X) == 1) 
  X = t(X)
if (plotROC) 
  alg = "ROC"
nR = nrow(X)
nC = ncol(X)
nY = table(y)
uL = as.factor(rownames(nY))
nL = length(nY)
if (nL <= 1) 
  stop("colAUC: List of labels 'y' have to contain at least 2 class labels.")
if (!is.numeric(X)) 
  stop("colAUC: 'X' must be numeric")
if (nR != length(y)) 
  stop("colAUC: length(y) and nrow(X) must be the same")
L = matrix(rep(uL, each = nR), nR, nL)
per = rbind(combs(1:nL, 2),combs(nL:1, 2))
nP = nrow(per)
Auc = matrix(0.5, nP, nC)
rownames(Auc) = paste(uL[per[, 1]], " vs. ", uL[per[, 2]], 
                      sep = "")
colnames(Auc) = colnames(X)
if (plotROC) {
  plot(c(0, 1), c(0, 1), type = "n", xaxs = "i", yaxs = "i", 
       xlab = "probability of false alarm", sub = "(1-Specificity)", 
       ylab = "probability of detection\n(Sensitivity)")
  title("ROC Curves")
  abline(h = 0:10/10, v = 0:10/10, col = "lightgray")
  if (nC * nP < 20) {
    S = colnames(Auc)
    if (is.null(S)) 
      S = paste("col", 1:nC)
    if (nP > 1) 
      S = paste(rep(S, each = nP), "[", rownames(Auc), 
                "]")
    legend("bottomright", S, col = 1:(nC * nP), lty = 1, 
           lwd = 1, pch = 20, merge = TRUE, inset = 0.01, 
           bg = "white")
  }
  nClr = 1
}


for (j in 1:nC) {
  d = (matrix(rep(y, nL), nR, nL) == L)
  d2 = d*1
  for (i in 1:nL) d[, i] = cumsum(d[, i])
  nD = nrow(d)
  for (i in 1:nP) {
    c1 = per[i, 1]
    c2 = per[i, 2]
    Auc[i,j]=((t(d2[,c1])%*%X[,j])/d[nD,c1]-(t(d2[,c2])%*%X[,j])/d[nD,c2]+1)/2
  }
  
}

Auc = pmax(Auc, 1 - Auc)
return(Auc)
}