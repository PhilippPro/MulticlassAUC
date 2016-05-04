rm(list = ls())

library(mlr)
library(HandTill2001)
library(caTools)

# Learner
learner <- makeLearner('classif.lda', predict.type="prob")

## Data classification example binary-------
n=1e2
Beta1=c(2, 5, -2)
X1<-rnorm(n)
X2<-rnorm(n)
X=cbind(1,X1,X2)
product=X%*%Beta1
Y<-as.factor(rbinom(n,1,prob = plogis(product)))
df1<-data.frame(X1,X2,Y)

task.2d<-makeClassifTask(data=df1, target = "Y")
mod.2d<- train(learner=learner, task=task.2d)
pred.obj.2d <- predict(mod.2d, newdata=df1)
library(HandTill2001)
predicted.2d <- as.matrix(pred.obj.2d$data[,paste("prob.", levels(pred.obj.2d$data$response), sep="")])
colnames(predicted.2d)<-levels(pred.obj.2d$data$response)
truth.2d=df1$Y


## Data example multiclass-----
task <- makeClassifTask(data=iris, target="Species")
mod <- train(learner=learner, task=task)
pred.obj <- predict(mod, newdata=iris)
library(HandTill2001)
predicted <- as.matrix(pred.obj$data[,paste("prob.", levels(pred.obj$data$response), sep="")])
colnames(predicted)<-levels(pred.obj$data$response)
#auc(multcap(response=pred.obj$data$response, predicted=predicted))
truth=iris$Species




## Implementations AUC -----

# imp 1 ----

multiclass.auc1 = function(pred, truth){
  ptm<-proc.time()
  colAUC(pred,truth, plotROC = TRUE)
  print(proc.time()-ptm)
}


#imp 2 pROC -----
multiclass.auc2 = function(pred, truth){
  if(nlevels(truth) == 2) {pred = cbind(pred,1-pred)}
  predP = pred
  # choose the probablity of the truth
  predV = vnapply(seq_row(pred), function(i) {
    pred[i, truth[i]]
  })
  auc = pROC::multiclass.roc(response = truth, predictor = predV)$auc
  as.numeric(auc)
}


multiclass.auc2 = function(pred, truth){
  auc = pROC::multiclass.roc(response = truth, predictor = pred)$auc
  as.numeric(auc)
}



#imp3 implemented by hand # -----
multiclass.auc3 = function(pred, truth,  names=FALSE) {
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



## Scored AUC

#imp.scored3 implemented by hand # -----
multiclass.auc.scored.3 = function(pred, truth,  names=FALSE) {
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
    matrixI.tempj=makeMatrixI(j,j,pred)
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

multiclass.sauc.scored.3 = function(pred, truth,  names=FALSE) {
  AUC.scored.matrix<-multiclass.auc.scored.3(pred, truth,  names=FALSE)
  mean(AUC.scored.matrix[row(AUC.scored.matrix)!=col(AUC.scored.matrix)])
}


## Probabilistic AUC ------


ProbabilisticAUC = function(pred, truth, names=FALSE) {
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
    matrixI.tempj=makeMatrixI(j,j,pred)
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


multiclassPAUC = function(pred, truth,  names=FALSE) {
  pAUC.scored.matrix<-ProbabilisticAUC(pred, truth,  names=FALSE)
  mean(pAUC.scored.matrix[row(pAUC.scored.matrix)!=col(pAUC.scored.matrix)])
}






## AU1U (Hand and Till) -----
m = colAUC(getPredictionProbabilities(pred), truth)
mean(m[col(m) != ncol(m) + 1 - row(m)])
# or
c = t(combn(1:nlevels(truth), 2, simplify = T))
mean(m[cbind(1:nrow(c), c[, 2])])

## Computation AUC binary ------
res1.2d=multiclass.auc1(predicted.2d, truth.2d)
res3.2d=multiclass.auc3(predicted.2d, truth.2d, names = TRUE)




## Computation AUC ----

#AUC1
res1=multiclass.auc1(predicted, truth)

#AUC2
res2=multiclass.auc2(predicted[,1], truth)

#AUC3
res3=multiclass.auc3(predicted, truth, names = TRUE)

res.scored3=multiclass.auc.scored.3(predicted, truth, names = TRUE)
res.scored3.sauc=multiclass.sauc.scored.3(predicted, truth, names = TRUE)


res.probAUC<-ProbabilisticAUC(pred, truth, names=FALSE)

res.pAUC.mean<-multiclassPAUC(pred, truth, names=FALSE)

## Explication Philipp
res1.1col=multiclass.auc1(predicted[,2], truth)


i=3
res3.1col=multiclass.auc3(cbind(data.matrix(predicted[,i]),(1-data.matrix(predicted[,i]))), truth)
res3.1col


# Hypothese : le auc va chercher pour telle colonne seulemetn deux classes, et va calculer l'AUC

# prenons la premiere colonne de predicted : setosa
predicted.setosa=predicted[,1]

# apres troisieme ligne par exemple c'est versicolor vs virginia, donc prenons ces obs la
truth.versiversusvirg=truth[51:150]
predicted.setosa.versiversusvirg=predicted.setosa[51:150]

res=multiclass.auc3(cbind(data.matrix(predicted.setosa.versiversusvirg),1-data.matrix(predicted.setosa.versiversusvirg)),as.factor(data.matrix(truth.versiversusvirg)))


# Prenons la deuxieme colonnes de predicted : versicolor
predicted.versicolor=predicted[,2]

# apres troisieme ligne par exemple c'est versicolor vs virginia, donc prenons ces obs la
predicted.versicolor.versiversusvirg=predicted.versicolor[51:150]

res=multiclass.auc3(cbind(data.matrix(predicted.versicolor.versiversusvirg),1-data.matrix(predicted.versicolor.versiversusvirg)),as.factor(data.matrix(truth.versiversusvirg)))


# Prenons la deuxieme colonnes de predicted : versicolor
predicted.virginica=predicted[,3]

# apres troisieme ligne par exemple c'est versicolor vs virginia, donc prenons ces obs la
predicted.versicolor.versiversusvirg=predicted.versicolor[51:150]

res=multiclass.auc3(cbind(data.matrix(predicted.versicolor.versiversusvirg),1-data.matrix(predicted.versicolor.versiversusvirg)),as.factor(data.matrix(truth.versiversusvirg)))
