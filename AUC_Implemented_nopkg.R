rm(list = ls())


## The functions ------

## Basic AUC -----
binaryclass.auc = function(pred, truth,  names=FALSE) {
  ptm <- proc.time()
  # Parameters
  K=dim(pred)[2]
  n=dim(pred)[1]
  levels=levels(truth)
  
  # Making the matrix of probabilities p(i,j) is pred
  
  # function I
  I=function(a,b){
    res=0
    if (a>b){res=1}
    if (a==b){res=0.5}
    return(res)
  }
  
  # function makeMatrixI
  makeMatrixI=function(j, k, pred){
    n=dim(pred)[1]
    matrix.I.temp=matrix(NA,nrow = n,ncol = n)
    for (i.i in c(1:n)) {
      for (i.t in c(1:n)) {
        matrix.I.temp[i.i,i.t]=I(pred[i.i,j],pred[i.t,k])
      }
    }
    return(matrix.I.temp)
  }
  
  
  # Design Matrix f
  f= model.matrix( ~ 0 + truth, truth)
  f.summed=apply(f,2,sum)
  
  # matrix AUC
  AUC=matrix(rep(0.5, K*K),K,K)
  
  for (j in c(1:K)) {
    matrixI.tempj=makeMatrixI(j,j,pred)
    for (k in c(1:K)) {
      if (k!=j) {
        value=t(f[,j])%*%matrixI.tempj%*%f[,k]/(f.summed[j]*f.summed[k])
        AUC[j,k]=max(value,1-value)
      }
    }
  }
  
  if (names==TRUE) {
    colnames(AUC)<-levels
    rownames(AUC)<-levels
  }
  print(proc.time()-ptm)
  return(AUC)
}


## AUNU
multiclass.aunu= function(pred, truth,  names=FALSE) {
  mean(vnapply(1:nlevels(truth), function(i) binaryclass.auc(cbind(pred[,i],1-pred[,i]), as.factor(truth == levels(truth)[i]))[1,2]))
}

## AUNP
multiclass.aunp= function(pred, truth,  names=FALSE) {
  # Matrix f
  f= model.matrix( ~ 0 + truth, truth)
  n=dim(f)[1]
  p=dim(f)[2]
  f.proba=apply(f,2,sum)/n
  
  # Computation
  dim(pred)[2]*mean(vnapply(1:nlevels(truth), function(i) f.proba[i]*binaryclass.auc(cbind(pred[,i],1-pred[,i]), as.factor(truth == levels(truth)[i]))[1,2]))
}



## AU1U
multiclass.au1u= function(pred, truth,  names=FALSE) {
  AUC.matrix<-binaryclass.auc(pred, truth,  names=FALSE)
  mean(AUC.matrix[row(AUC.matrix)!=col(AUC.matrix)])
}

## AU1P
multiclass.au1p= function(pred, truth,  names=FALSE) {
  # Matrix f
  f= model.matrix( ~ 0 + truth, truth)
  n=dim(f)[1]
  p=dim(f)[2]
  f.proba=apply(f,2,sum)/n
  
  # Computation
  AUC.matrix<-binaryclass.auc(pred, truth,  names=FALSE)
  AUC.matrix.weightsproba<-t(t(AUC.matrix)*f.proba)
  return(mean(AUC.matrix.weightsproba[row(AUC.matrix.weightsproba)!=col(AUC.matrix.weightsproba)])*p)
}


## Scored AUC

binaryclass.Scoredauc = function(pred, truth,  names=FALSE) {
  ptm <- proc.time()
  # Parameters
  K=dim(pred)[2]
  n=dim(pred)[1]
  levels=levels(truth)
  print(paste("levels",levels))
  
  # Making the matrix of probabilities p(i,j) is pred
  
  # function I
  I=function(a,b){
    res=0
    if (a>b){res=1}
    if (a==b){res=0.5}
    return(res)
  }
  
  # function makeMatrixI
  makeMatrixI.scored=function(j, k, pred){
    n=dim(pred)[1]
    matrix.I.temp=matrix(NA,nrow = n,ncol = n)
    for (i.i in c(1:n)) {
      for (i.t in c(1:n)) {
        matrix.I.temp[i.i,i.t]=I(pred[i.i,j],pred[i.t,k])*(pred[i.i,j]-pred[i.t,k])
      }
    }
    return(matrix.I.temp)
  }
  
  
  # Design Matrix f
  f= model.matrix( ~ 0 + truth, truth)
  f.summed=apply(f,2,sum)
  
  # matrix AUC
  AUC=matrix(rep(NA, K*K),K,K)
  
  for (j in c(1:K)) {
    matrixI.tempj=makeMatrixI.scored(j,j,pred)
    for (k in c(1:K)) {
      if (k!=j) {
        AUC[j,k]=t(f[,j])%*%matrixI.tempj%*%f[,k]/(f.summed[j]*f.summed[k])
      }
    }
  }
  
  if (names==TRUE) {
    colnames(AUC)<-levels
    rownames(AUC)<-levels
  }
  print(proc.time()-ptm)
  return(AUC)
}


## Scored SAUC

multiclass.Scoredauc= function(pred, truth,  names=FALSE) {
  AUC.scored.matrix<-binaryclass.Scoredauc(pred, truth,  names=FALSE)
  mean(AUC.scored.matrix[row(AUC.scored.matrix)!=col(AUC.scored.matrix)])
}


## Probabilistic AUC


binaryclass.Probabilisticauc = function(pred, truth, names=FALSE) {
  ptm <- proc.time()
  # Parameters
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
        probAUC[j,k]=((t(f[,j])%*%pred[,j])/f.summed[j]-(t(f[,k])%*%pred[,j])/f.summed[k]+1)/2
      }
    }
  }
  
  if (names==TRUE) {
    colnames(probAUC)<-levels
    rownames(probAUC)<-levels
  }
  print(proc.time()-ptm)
  return(probAUC)
}

## multiclass PAUC one v one


multiclass.Probabilisticauc = function(pred, truth,  names=FALSE) {
  pAUC.scored.matrix<-binaryclass.Probabilisticauc(pred, truth,  names=FALSE)
  mean(pAUC.scored.matrix[row(pAUC.scored.matrix)!=col(pAUC.scored.matrix)])
}


## Fast Implementations -------


## Fast AUC -----
binaryclass.auc.fast = function(pred, truth,  names=FALSE) {
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
  Auc = matrix(0.5, K, K) # Creates the AUC matrix
  Auc2 = matrix(NA, nP, K) # Creates the AUC matrix
  
  
  for (j in c(1:K)) {
    x = sort(X[, j], index = TRUE) # sorting with increasing probability for class j
    nunq = which(diff(x$x) == 0) # Which probability values are equal ?
    nTies = length(nunq) # How many are equl ?
    
    if (nTies < n - 1) { # Do all the probabilities are equal ?
      idx = y[x$ix] # indexes by increasing order of probabilities
      d = (matrix(rep(idx, K), n, K) == L) # Transforming the predicted vector in K columns of TRUE/FALSE
      d=apply(d,2,cumsum) # Cumulative sum of number in class
      if (nTies) # If at least one of the probability equals one another
        d = d[-nunq, ] # We remove the doublons
      d = rbind(matrix(0, 1, K), d) # Add one line for integral computation
      nD = nrow(d) # New number of lines for dj
      
      
      for (i in 1:K) {
        
        if (i!=j) {
          number = d[nD, i] * d[nD, j]
          if (number > 0) 
            Auc[i, j] = trapz(x=d[, i], y=d[, j])/number # Number of non-variations of j regarding constant values of i
          # --> gives the number of false negative when trying to predict j
        }
      }
      
      
      for (i in 1:nP) {
        c1 = permutations[i, 1]
        c2 = permutations[i, 2]
        number = d[nD, c1] * d[nD, c2]
        if ((number > 0)&&((j==c1)||(j==c2))) 
          Auc2[i, j] = trapz(x=d[, c1], y=d[, c2])/number # En gros le nombre de non variations de c2 par rapport a c1 qui est constant
      }
      
      
      
    }
  }
  
  
  # Add the names
  if (names==TRUE) {
    colnames(Auc)<-levels.names
    rownames(Auc)<-levels.names
  }
  Auc<-apply(Auc,c(1,2),function(x) return(max(x,1-x)))
  Auc2<-apply(Auc2,c(1,2),function(x) return(max(x,1-x)))
  rownames(Auc2) = paste(levels.values[permutations[, 1]], " vs. ", levels.values[permutations[, 2]], 
                        sep = "")
  colnames(Auc2) = colnames(X)
  
  print(proc.time()-ptm)
  return(Auc)
}


## AUNU fast
multiclass.aunu.fast= function(pred, truth,  names=FALSE) {
  mean(vnapply(1:nlevels(truth), function(i) binaryclass.auc.fast(cbind(pred[,i],1-pred[,i]), as.factor(truth == levels(truth)[i]))[1,2]))
}

## AUNP fast
multiclass.aunp.fast= function(pred, truth,  names=FALSE) {
  # Matrix f
  f= model.matrix( ~ 0 + truth, truth)
  n=dim(f)[1]
  p=dim(f)[2]
  f.proba=apply(f,2,sum)/n
  
  # Computation
  dim(pred)[2]*mean(vnapply(1:nlevels(truth), function(i) f.proba[i]*binaryclass.auc.fast(cbind(pred[,i],1-pred[,i]), as.factor(truth == levels(truth)[i]))[1,2]))
}



## AU1U fast
multiclass.au1u.fast= function(pred, truth,  names=FALSE) {
  AUC.matrix<-binaryclass.auc.fast(pred, truth,  names=FALSE)
  mean(AUC.matrix[row(AUC.matrix)!=col(AUC.matrix)])
}

## AU1P fast
multiclass.au1p.fast= function(pred, truth,  names=FALSE) {
  # Matrix f
  f= model.matrix( ~ 0 + truth, truth)
  n=dim(f)[1]
  p=dim(f)[2]
  f.proba=apply(f,2,sum)/n
  
  # Computation
  AUC.matrix<-binaryclass.auc.fast(pred, truth,  names=FALSE)
  AUC.matrix.weightsproba<-t(t(AUC.matrix)*f.proba)
  return(mean(AUC.matrix.weightsproba[row(AUC.matrix.weightsproba)!=col(AUC.matrix.weightsproba)])*p)
}


## Fast Scored auc
binaryclass.Scoredauc.fast = function(pred, truth,  names=FALSE) {
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
  Auc = matrix(0.5, K, K) # Creates the AUC matrix
  Auc2 = matrix(NA, nP, K) # Creates the AUC matrix
  
  
  for (j in c(1:K)) {
    x = sort(X[, j], index = TRUE) # sorting with increasing probability for class j
    nunq = which(diff(x$x) == 0) # Which probability values are equal ?
    nTies = length(nunq) # How many are equl ?
    
    if (nTies < n - 1) { # Do all the probabilities are equal ?
      idx = y[x$ix] # indexes by increasing order of probabilities
      d = (matrix(rep(idx, K), n, K) == L) # Transforming the predicted vector in K columns of TRUE/FALSE
      d=apply(d,2,cumsum) # Cumulative sum of number in class
      if (nTies) # If at least one of the probability equals one another
        d = d[-nunq, ] # We remove the doublons
      d = rbind(matrix(0, 1, K), d) # Add one line for integral computation
      nD = nrow(d) # New number of lines for dj
      
      
      for (i in 1:K) {
        
        if (i!=j) {
          number = d[nD, i] * d[nD, j]
          if (number > 0) 
            Auc[i, j] = trapz(x=d[, i], y=d[, j])/number # Number of non-variations of j regarding constant values of i
          # --> gives the number of false negative when trying to predict j
        }
      }
      
      
      for (i in 1:nP) {
        c1 = permutations[i, 1]
        c2 = permutations[i, 2]
        number = d[nD, c1] * d[nD, c2]
        if ((number > 0)&&((j==c1)||(j==c2))) 
          Auc2[i, j] = trapz(x=d[, c1], y=d[, c2])/number # En gros le nombre de non variations de c2 par rapport a c1 qui est constant
      }
      
      
      
    }
  }
  
  
  # Add the names
  if (names==TRUE) {
    colnames(Auc)<-levels.names
    rownames(Auc)<-levels.names
  }
  Auc<-apply(Auc,c(1,2),function(x) return(max(x,1-x)))
  Auc2<-apply(Auc2,c(1,2),function(x) return(max(x,1-x)))
  rownames(Auc2) = paste(levels.values[permutations[, 1]], " vs. ", levels.values[permutations[, 2]], 
                         sep = "")
  colnames(Auc2) = colnames(X)
  
  print(proc.time()-ptm)
  return(Auc)
}




###############################################
## Examples -----
###############################################

library(HandTill2001)
library(mlr)

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
multiclass.aunu(predicted, truth, names = TRUE)
multiclass.aunp(predicted, truth, names = TRUE)
multiclass.au1u(predicted, truth, names = TRUE)
multiclass.au1p(predicted, truth, names = TRUE)

#Scored AUC
binaryclass.Scoredauc(predicted, truth, names = TRUE)
multiclass.Scoredauc(predicted, truth, names = TRUE)

# Probabilistic AUC
binaryclass.Probabilisticauc(predicted, truth, names = TRUE)
multiclass.Probabilisticauc(predicted, truth, names = TRUE)

## Fast computation of AUC measures ------
# (trick from ca_tool package of considering only the errors false negative)

# Classic AUC fast
binaryclass.auc.fast(predicted, truth, names = TRUE)
multiclass.aunu.fast(predicted, truth, names = TRUE)
multiclass.aunp.fast(predicted, truth, names = TRUE)
multiclass.au1u.fast(predicted, truth, names = TRUE)
multiclass.au1p.fast(predicted, truth, names = TRUE)

# Peut être plutot placer un NA dans les endroits ou ca n'a pas de sens ??? et faire le même tableau qu'eux


