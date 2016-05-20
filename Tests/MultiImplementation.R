binaryclass.multiAUC<-function (X, y, plotROC = FALSE, TypeofAUC = c("auc", "sauc", "pauc")) 
{
  ptm <- proc.time()
  y = as.factor(y)
  X = as.matrix(X)
  TypeofAUC = match.arg(TypeofAUC)
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
  per = combs(1:nL, 2)
  if ((TypeofAUC=="sauc")||(TypeofAUC=="pauc"))
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
  
  # Case auc
  if (TypeofAUC=="auc") {
    
    for (j in 1:nC) {
      x = sort(X[, j], index = TRUE)
      nunq = which(diff(x$x) == 0)
      nTies = length(nunq)
      if (nTies < nR - 1) {
        idx = y[x$ix]
        d = (matrix(rep(idx, nL), nR, nL) == L)
        for (i in 1:nL) d[, i] = cumsum(d[, i])
        if (nTies) 
          d = d[-nunq, ]
        d = rbind(matrix(0, 1, nL), d)
        nD = nrow(d)
        for (i in 1:nP) {
          c1 = per[i, 1]
          c2 = per[i, 2]
          n = d[nD, c1] * d[nD, c2]
          if (n > 0) 
            Auc[i, j] = trapz(d[, c1], d[, c2])/n
          if (plotROC) {
            xx = if (n > 0) 
              d[, c1]/d[nD, c1]
            else c(0, 1)
            yy = if (n > 0) 
              d[, c2]/d[nD, c2]
            else c(0, 1)
            if (2 * Auc[i, j] < 1) {
              xx = 1 - xx
              yy = 1 - yy
            }
            lines(xx, yy, col = nClr, type = "o", pch = 20)
            nClr = nClr + 1
          }
        }
      }
    }
  }
  
  # Case sauc
  if (TypeofAUC=="sauc") {
    for (j in c(1:nC)) {                      # For all the classes probability predictions
      x = sort(X[, j], index = TRUE)              # sort with increasing probability for class j
      idx = y[x$ix]                               # get the resulting indexes
      d = (matrix(rep(idx, nL), nR, nL) == L)        # Transform the predicted vector in K columns which have TRUE/FALSE values
      d2=d*1
      d=apply(d,2,cumsum)                         # Do the cumulative sum of element for each class
      nD = nrow(d)                                                                    # Get the previously created matrix n*K with 1 and 0 instead of TRUE/FALSE
      
      for (i in 1:nP) {                           # For all the possible permutations
        c1 = per[i, 1]
        c2 = per[i, 2]
        number = d[nD, c1] * d[nD, c2]              # Number of elements in the cross product
        sum = 0
        
        positiveProbability.number = d2[,c1]*d[,c2]
        negativeProbability.number = d2[,c2]*(d[nD, c1]-d[,c1])
        combined.number = positiveProbability.number-negativeProbability.number
        res = sum(x$x*combined.number)
        Auc[i, j] = res/number 
      }
    }
  }
  
  # Case pauc
  if (TypeofAUC=="pauc") {
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
  }
  
  
  Auc = pmax(Auc, 1 - Auc)
  print(proc.time()-ptm)
  return(Auc)
}