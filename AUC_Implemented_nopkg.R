rm(list = ls())


## The functions ------

## Basic AUC -----
binaryclass.auc = function(pred, truth,  names=FALSE) {
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
  AUC=matrix(rep(NA, K*K),K,K)
  
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


## AUNU ------
multiclass.aunu= function(pred, truth,  names=FALSE) {
  mean(vnapply(1:nlevels(truth), function(i) binaryclass.auc(cbind(pred[,i],1-pred[,i]), as.factor(truth == levels(truth)[i]))[1,2]))
}

## AUNP ------
multiclass.aunp= function(pred, truth,  names=FALSE) {
  # Matrix f
  f= model.matrix( ~ 0 + truth, truth)
  n=dim(f)[1]
  p=dim(f)[2]
  f.proba=apply(f,2,sum)/n
  
  # Computation
  dim(pred)[2]*mean(vnapply(1:nlevels(truth), function(i) f.proba[i]*binaryclass.auc(cbind(pred[,i],1-pred[,i]), as.factor(truth == levels(truth)[i]))[1,2]))
}



## AU1U ------
multiclass.au1u= function(pred, truth,  names=FALSE) {
  AUC.matrix<-binaryclass.auc(pred, truth,  names=FALSE)
  mean(AUC.matrix[row(AUC.matrix)!=col(AUC.matrix)])
}

## AU1P ------
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


## Scored AUC ------

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


## Scored SAUC ------

multiclass.Scoredauc= function(pred, truth,  names=FALSE) {
  AUC.scored.matrix<-binaryclass.Scoredauc(pred, truth,  names=FALSE)
  mean(AUC.scored.matrix[row(AUC.scored.matrix)!=col(AUC.scored.matrix)])
}


## Probabilistic AUC ------


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

## multiclass PAUC one v one -----


multiclass.Probabilisticauc = function(pred, truth,  names=FALSE) {
  pAUC.scored.matrix<-binaryclass.Probabilisticauc(pred, truth,  names=FALSE)
  mean(pAUC.scored.matrix[row(pAUC.scored.matrix)!=col(pAUC.scored.matrix)])
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

## Computation of AUC measures

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
