# File: R/measures_colAUC.R
# mainly copied from the caTools package
colAUC = function(probabilities, truth) {
  y = as.factor(truth)
  X = as.matrix(probabilities)
  if (nrow(X) == 1)
    X = t(X)
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
  rownames(Auc) = paste(uL[per[, 1]], " vs. ", uL[per[, 2]], sep = "")
  colnames(Auc) = colnames(X)
  # Wilcoxon AUC
  idxL = vector(mode = "list", length = nL)
  for (i in 1:nL) idxL[[i]] = which(y == uL[i])
  for (j in 1:nC) {
    for (i in 1:nP) {
      c1 = per[i, 1]
      c2 = per[i, 2]
      n1 = as.numeric(nY[c1])
      n2 = as.numeric(nY[c2])
      if (n1 > 0 & n2 > 0) {
        r = rank(c(X[idxL[[c1]], j], X[idxL[[c2]], j]))
        Auc[i, j] = (sum(r[1:n1]) - n1 * (n1 + 1)/2)/(n1 * n2)
      }
    }
  }
  Auc = pmax(Auc, 1 - Auc)
  return(Auc)
}

combs = function(v,k) {
  # combs(V,K) - finds all unordered combinations of K elements from vector V 
  #  V is a vector of length N
  #  K is a integer 
  # combs(V,K) creates a matrix with N!/((N-K)! K!) rows
  # and K columns containing all possible combinations of N elements taken K at a time.
  # example: combs(1:3,2) returns matrix with following rows (1 2), (1 3), (2 3)
  n = length(v)
  if      (n==k  ) P = matrix(v,1,n)
  else if (k==1  ) P = matrix(v,n,1)
  else if (k==n-1) P = matrix( rep(v, each=n-1), n, n-1)
  else if (k< n) {
    P = matrix(0,0,k)
    if (k < n & k > 1) {
      for (i in 1:(n-k+1)) {
        Q = combs(v[(i+1):n],k-1)
        j = nrow(Q)
        P = rbind(P, cbind(rep(v[i],j), Q))
      }
    }
  } else 
    stop("combs: number m has to be smaller or equal to length of vector v")
  return(P)
}


binaryclass.multiAUC<-function (X, y, plotROC = FALSE, TypeofAUC = c("auc", "sauc", )) 
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
  if (TypeofAUC=="sauc")
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
  
  
  Auc = pmax(Auc, 1 - Auc)
  print(proc.time()-ptm)
  return(Auc)
}


# File: R/measures.R (delete old multiclass.auc!)
#' @export multiclass.aunu
#' @rdname measures
#' @format none
multiclass.aunu = makeMeasure(id = "multiclass.aunu", minimize = FALSE, best = 0, worst = 1,
                               properties = c("classif", "classif.multi", "req.pred", "req.truth", "req.prob"),
                               name = "Average multiclass AUC",
                               note = "Following the definition in the Ferri et. al Paper: https://www.math.ucdavis.edu/~saito/data/roc/ferri-class-perf-metrics.pdf",                             
                               fun = function(task, model, pred, feats, extra.args) {
                                   measureAUNU(getPredictionProbabilities(pred), pred$data$truth)
                               }
)

#' @export measureAUNU
#' @rdname measures
#' @format none
measureAUNU = function(probabilities, truth) {
  if (nlevels(truth) == 2L) probabilities = cbind(probabilities,1-probabilities)
  mean(vnapply(1:nlevels(truth), function(i) colAUC(probabilities[, i], truth == levels(truth)[i])))
}

#' @export multiclass.aunp
#' @rdname measures
#' @format none
multiclass.aunp = makeMeasure(id = "multiclass.aunp", minimize = FALSE, best = 0, worst = 1,
                              properties = c("classif", "classif.multi", "req.pred", "req.truth", "req.prob"),
                              name = "Weighted average multiclass AUC",
                              note = "Following the definition in the Ferri et. al paper: https://www.math.ucdavis.edu/~saito/data/roc/ferri-class-perf-metrics.pdf",                             
                              fun = function(task, model, pred, feats, extra.args) {
                                measureAUNP(getPredictionProbabilities(pred), pred$data$truth)
                              }
)

#' @export measureAUNP
#' @rdname measures
#' @format none
measureAUNP = function(probabilities, truth) {
  if (nlevels(truth) == 2L) probabilities = cbind(probabilities,1-probabilities)
  sum(vnapply(1:nlevels(truth), function(i) mean(truth == levels(truth)[i]) * colAUC(probabilities[,i], truth == levels(truth)[i])))  
}

#' @export multiclass.au1u
#' @rdname measures
#' @format none
multiclass.au1u = makeMeasure(id = "multiclass.au1u", minimize = FALSE, best = 0, worst = 1,
                              properties = c("classif", "classif.multi", "req.pred", "req.truth", "req.prob"),
                              name = "Average multiclass AUC with each class against each other",
                              note = "Following the definition in the Ferri et. al paper: https://www.math.ucdavis.edu/~saito/data/roc/ferri-class-perf-metrics.pdf",                             
                              fun = function(task, model, pred, feats, extra.args) {
                                measureAU1U(getPredictionProbabilities(pred), pred$data$truth)
                              }
)

#' @export measureAU1U
#' @rdname measures
#' @format none
measureAU1U = function(probabilities, truth) {
  if (nlevels(truth) == 2L) probabilities = cbind(probabilities,1-probabilities)
  m = colAUC(probabilities, truth)
  c = c(combn(1:nlevels(truth), 2))
  mean(m[cbind(rep(1:nrow(m), each = 2), c)])
}

#' @export multiclass.au1p
#' @rdname measures
#' @format none
multiclass.au1p = makeMeasure(id = "multiclass.au1p", minimize = FALSE, best = 0, worst = 1,
                              properties = c("classif", "classif.multi", "req.pred", "req.truth", "req.prob"),
                              name = "Weighted average multiclass AUC with each class against each other",
                              note = "Following the definition in the Ferri et. al paper: https://www.math.ucdavis.edu/~saito/data/roc/ferri-class-perf-metrics.pdf",                             
                              fun = function(task, model, pred, feats, extra.args) {
                                measureAU1P(getPredictionProbabilities(pred), pred$data$truth)
                              }
)

#' @export measureAU1P
#' @rdname measures
#' @format none
measureAU1P = function(probabilities, truth) {
  if (nlevels(truth) == 2L) probabilities = cbind(probabilities,1-probabilities)
  m = colAUC(probabilities, truth)
  weights = table(truth) / length(truth)
  m = m * matrix(rep(weights, each = nrow(m)), ncol = length(weights))
  c = c(combn(1:nlevels(truth), 2))
  sum(m[cbind(rep(1:nrow(m), each = 2), c)]) / (nlevels(truth) - 1)
}

# File: tests/testthat/test_base_measures.R (delete old multiclass.auc)
# not sure about how the tests should look like!