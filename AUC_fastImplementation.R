function (X, y, plotROC = FALSE, alg = c("Wilcoxon", "ROC")) 
{
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
  per = combs(1:nL, 2)
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
  if (alg == "Wilcoxon") {
    idxL = vector(mode = "list", length = nL)
    for (i in 1:nL) idxL[[i]] = which(y == uL[i])
    for (j in 1:nC) {
      for (i in 1:nP) {
        c1 = per[i, 1]
        c2 = per[i, 2]
        n1 = as.numeric(nY[c1])
        n2 = as.numeric(nY[c2])
        if (n1 > 0 & n2 > 0) {
          r = rank(c(X[idxL[[c1]], j], X[idxL[[c2]], 
                                         j]))
          Auc[i, j] = (sum(r[1:n1]) - n1 * (n1 + 1)/2)/(n1 * 
                                                          n2)
        }
      }
    }
  }
  if (alg == "ROC") {
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
  Auc = pmax(Auc, 1 - Auc)
  return(Auc)
}