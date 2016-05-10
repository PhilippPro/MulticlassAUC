rm(list = ls())


## The functions ------

binaryclass.auc = function(pred, truth,  names=FALSE) {
  ptm <- proc.time()
  # Parameters
  y<-truth
  X<-as.matrix(pred)
  ncol=ncol(X) #maj k
  n=nrow(X)
  levels.values=as.factor(levels(as.factor(truth)))
  levels.names=levels(as.factor(truth))
  K=length(levels.names)
  L = matrix(rep(levels.values, each = n), n, K)
  permutations = combs(1:K, 2)
  nP = nrow(permutations)
  Auc = matrix(0.5, nP, ncol) # Creates the AUC matrix
  
  
  for (j in c(1:ncol)) {
    x = sort(X[, j], index = TRUE) # sorting with increasing probability for class j
    
    idx = y[x$ix] # indexes by increasing order of probabilities
    d = (matrix(rep(idx, K), n, K) == L) # Transforming the predicted vector in K columns of TRUE/FALSE
    d1=d
    d=apply(d,2,cumsum) # Cumulative sum of number in class
    dd = rbind(matrix(0, 1, K), d) # Add one line for integral computation
    nD = nrow(d) # New number of lines for dj
    d2=d1*1
    
    for (i in 1:nP) {
      c1 = permutations[i, 1]
      c2 = permutations[i, 2]
      number = d[nD, c1] * d[nD, c2]
      sum=0
      for (p in c(1:n)) {
        if (d2[p,c1]==1) {
          #sum=sum+length(which(d2[c(1:p),c2]==1))
          sum=sum+d[p,c2]
          index=which(d2[c(1:p),c2]==1)
        }
      }
      Auc[i, j] = sum/number # En gros le nombre de non variations de c2 par rapport a c1 qui est constant
    }
  }
  
  Auc<-apply(Auc,c(1,2),function(x) return(max(x,1-x)))
  # Add the names
  if (names==TRUE) {
    rownames(Auc) = paste(levels.values[permutations[, 1]], " vs. ", levels.values[permutations[, 2]], sep = "")
    colnames(Auc) = colnames(X)
  }
  
  print(proc.time()-ptm)
  return(Auc)
}

## AUNU
multiclass.aunu= function(pred, truth,  names=FALSE) {
  mean(vnapply(1:nlevels(truth), function(i) binaryclass.auc(pred[,i], truth == levels(truth)[i])))
}

multiclass.aunu(pred, truth)

## AUNP
multiclass.aunp= function(pred, truth,  names=FALSE) {
  # Matrix f
  f= model.matrix( ~ 0 + truth, truth)
  n=dim(f)[1]
  p=dim(f)[2]
  f.proba=apply(f,2,sum)/n
  nClasses=length(levels(truth))
  
  # Computation
  nClasses*mean(vnapply(1:nlevels(truth), function(i) f.proba[i]*binaryclass.auc(pred[,i],  truth == levels(truth)[i])))
}

multiclass.aunp(pred,truth)

## AU1U
multiclass.au1u= function(pred, truth,  names=FALSE) {
  
  # Computes binary AUC
  AUC.matrix<-binaryclass.auc(pred, truth,  names=FALSE)
  
  # Retrieve the useful elements only
  K=nlevels(truth)
  permutations = combs(1:K, 2) # same permutation matrix as used in binaryclass.au1u
  nCol=ncol(AUC.matrix) # also number of classes
  nRow=nrow(AUC.matrix) # also number of vs 
  count= (nCol-1)*nCol # also number of element of matrix AUC
  sum=0 
  
  for (j in c(1:nCol)) {
    for (i in c(1:nRow)) {
      if (is.element(j,permutations[i,])) {
        sum=sum+AUC.matrix[i,j]
      }
    }
  }
  res=sum/count
  return(res)
}


## AU1P
multiclass.au1p= function(pred, truth,  names=FALSE) {
  
  # Computes binary AUC
  AUC.matrix<-binaryclass.auc(pred, truth,  names=FALSE)
  
  # Matrix f to get the class probabilities
  f= model.matrix( ~ 0 + truth, truth)
  n=dim(f)[1] # numberof observations
  K=dim(f)[2] # number of classes as number of columns in f
  f.proba=apply(f,2,sum)/n
  
  # Weights AUC matrix
  AUC.matrix.weightsproba<-t(t(AUC.matrix)*f.proba)
  
  # Retrieve the useful elements only
  permutations = combs(1:K, 2) # same permutation matrix as used in binaryclass.au1u
  nCol=ncol(AUC.matrix) # also number of classes
  nRow=nrow(AUC.matrix) # also number of vs 
  count= (nCol-1)*nCol # also number of element of matrix AUC
  sum=0 
  
  for (j in c(1:nCol)) {
    for (i in c(1:nRow)) {
      if (is.element(j,permutations[i,])) {
        sum=sum+AUC.matrix.weightsproba[i,j]
      }
    }
  }
  res=(sum/count)*K
  return(res)
}



## Scored AUC
binaryclass.Scoredauc = function(pred, truth,  names=FALSE) {
  ptm <- proc.time()
  # Parameters
  y<-truth
  X<-pred
  K=dim(pred)[2] #maj k
  n=dim(pred)[1]
  levels.values=as.factor(levels(truth))
  levels.names=levels(truth)
  L = matrix(rep(levels.values, each = n), n, K)
  permutations = combs(1:K, 2)
  nP = nrow(permutations)
  Auc = matrix(0.5, nP, K) # Creates the AUC matrixc
  count=0
  
  
  for (j in c(1:K)) {
    x = sort(X[, j], index = TRUE) # sorting with increasing probability for class j
    
    idx = y[x$ix] # indexes by increasing order of probabilities
    d = (matrix(rep(idx, K), n, K) == L) # Transforming the predicted vector in K columns of TRUE/FALSE
    d1=d
    d=apply(d,2,cumsum) # Cumulative sum of number in class
    dd = rbind(matrix(0, 1, K), d) # Add one line for integral computation
    nD = nrow(d) # New number of lines for dj
    d2=d1*1
    
    for (i in 1:nP) {
      c1 = permutations[i, 1]
      c2 = permutations[i, 2]
      number = d[nD, c1] * d[nD, c2]
      sum=0
      sumInt=0
      l<-NULL
      
      for (p in c(1:n)) {
        if (d2[p,c1]==1) {
          index=which(d2[c(1:p),c2]==1)
          if (length(index)>0) {
            for (q in c(1:length(index))) {
              indexq=index[q]
              #sum=sum+1
              sum=sum+pred[(x$ix[p]),c1]-pred[(x$ix[indexq]),c2]
              l<-c(l,pred[(x$ix[p]),c1]-pred[(x$ix[indexq]),c2])
              sumInt=sumInt+1
              count=count+1
            }
          }
        }
      }
      
      
      Auc[i, j] = sum/number # En gros le nombre de non variations de c2 par rapport a c1 qui est constant
    }
  }
  
  Auc<-abs(Auc)
  Auc<-apply(Auc,c(1,2),function(x) return(max(x,1-x)))
  # Add the names
  if (names==TRUE) {
    rownames(Auc) = paste(levels.values[permutations[, 1]], " vs. ", levels.values[permutations[, 2]], sep = "")
    colnames(Auc) = colnames(X)
  }
  
  print(proc.time()-ptm)
  print(paste("count",count))
  return(Auc)
}


## Scored SAUC
multiclass.Scoredauc= function(pred, truth,  names=FALSE) {
  
  # Computes the AUC scored matrix
  AUC.scored.matrix<-binaryclass.Scoredauc(pred, truth,  names=FALSE)
  
  # Retrieve the useful elements only
  K=nlevels(truth)
  permutations = combs(1:K, 2) # same permutation matrix as used in binaryclass.au1u
  nCol=ncol(AUC.matrix) # also number of classes
  nRow=nrow(AUC.matrix) # also number of vs 
  count= (nCol-1)*nCol # also number of element of matrix AUC
  sum=0 
  
  for (j in c(1:nCol)) {
    for (i in c(1:nRow)) {
      if (is.element(j,permutations[i,])) {
        sum=sum+AUC.scored.matrix[i,j]
      }
    }
  }
  res=sum/count
  return(res)
}


## Probabilistic AUC
binaryclass.Probabilisticauc = function(pred, truth,  names=FALSE) {
  ptm <- proc.time()
  
  # Parameters
  y<-truth
  X<-as.matrix(pred)
  ncol=ncol(X) 
  n=nrow(X)
  levels.values=as.factor(levels(as.factor(truth)))
  levels.names=levels(as.factor(truth))
  K=length(levels.names)
  L = matrix(rep(levels.values, each = n), n, K)
  permutations = combs(1:K, 2)
  nP = nrow(permutations)
  Auc.probabilistic = matrix(0.5, nP, ncol) 
  
  # Design Matrix f
  f= model.matrix( ~ 0 + truth, truth)
  f.summed=apply(f,2,sum)
  
  # Iterative
  for (j in c(1:ncol)) {
    for (i in 1:nP) {
      c1 = permutations[i, 1]
      c2 = permutations[i, 2]
      Auc.probabilistic[i,j]=((t(f[,c2])%*%pred[,j])/f.summed[c2]-(t(f[,c1])%*%pred[,j])/f.summed[c1]+1)/2
    }
  }
  
  Auc.probabilistic<-apply(Auc.probabilistic,c(1,2),function(x) return(max(x,1-x)))
  # Add the names
  if (names==TRUE) {
    rownames(Auc.probabilistic) = paste(levels.values[permutations[, 1]], " vs. ", levels.values[permutations[, 2]], sep = "")
    colnames(Auc.probabilistic) = colnames(X)
  }
  
  print(proc.time()-ptm)
  return(Auc.probabilistic)
}


## multiclass PAUC one v one
multiclass.Probabilisticauc = function(pred, truth,  names=FALSE) {
  
  # Computation of the probabilistic AUC
  AUC.probabilistic.matrix<-binaryclass.Probabilisticauc(pred, truth,  names=FALSE)
  
  # Retrieve the useful elements only
  K=nlevels(truth)
  permutations = combs(1:K, 2) # same permutation matrix as used in binaryclass.au1u
  nCol=ncol(AUC.matrix) # also number of classes
  nRow=nrow(AUC.matrix) # also number of vs 
  count= (nCol-1)*nCol # also number of element of matrix AUC
  sum=0 
  
  for (j in c(1:nCol)) {
    for (i in c(1:nRow)) {
      if (is.element(j,permutations[i,])) {
        sum=sum+AUC.probabilistic.matrix[i,j]
      }
    }
  }
  res=sum/count
  return(res)
}



###############################################
## Examples -----
###############################################

library(HandTill2001)
library(mlr)
library(caTools)

## Data Iris
learner <- makeLearner('classif.lda', predict.type="prob")
task <- makeClassifTask(data=iris, target="Species")
mod <- train(learner=learner, task=task)
pred.obj <- predict(mod, newdata=iris)
predicted <- as.matrix(pred.obj$data[,paste("prob.", levels(pred.obj$data$response), sep="")])
colnames(predicted)<-levels(pred.obj$data$response)
truth=iris$Species

## Computation of AUC measures------

# Classic AUC
binaryclass.auc(predicted, truth, names = TRUE)
binaryclass.auc(predicted, truth, names = TRUE)
multiclass.aunu(predicted, truth, names = TRUE)
multiclass.aunp(predicted, truth, names = TRUE)
multiclass.au1u(predicted, truth, names = TRUE)
multiclass.au1p(predicted, truth, names = TRUE)

#Scored AUC
binaryclass.Scoredauc(predicted, truth, names = TRUE)
binaryclass.ScoredaucOld(predicted, truth, names = TRUE)
multiclass.Scoredauc(predicted, truth, names = TRUE)

# Probabilistic AUC
binaryclass.Probabilisticauc(predicted, truth, names = TRUE)
multiclass.Probabilisticauc(predicted, truth, names = TRUE)


# tend to confirm the 2 first lines are in fact 1 because no problem between
i=1
test=cbind(predicted[,i],as.factor(iris$Species))
test[order(test[,1]),]
