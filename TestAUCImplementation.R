rm(list = ls())

library(mlr)
library(HandTill2001)

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
  colAUC(pred,truth)
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


#imp3 implemented by hand # -----
multiclass.auc3 = function(pred, truth) {
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
  makeMatrixI=function(j, pred){
    n=dim(pred)[1]
    matrix.I.temp=matrix(NA,nrow = n,ncol = n)
    for (i.i in c(1:n)) {
      for (i.t in c(1:n)) {
        matrix.I.temp[i.i,i.t]=I(pred[i.i,j],pred[i.t,j])
      }
    }
    return(matrix.I.temp)
  }
  
  makeMatrixI(1,pred)
  
  # Design Matrix f
  f= model.matrix( ~ 0 + truth, truth)
  f.summed=apply(f,2,sum)
  
  # matrix AUC
  AUC=matrix(rep(NA, K*K),K,K)
  
  for (j in c(1:K)) {
    for (k in c(1:K)) {
      AUC[j,k]=t(f[,j])%*%makeMatrixI(j,pred)%*%f[,k]/(f.summed[j]*f.summed[k])
      
    }
  }
  colnames(AUC)<-levels
  rownames(AUC)<-levels
  return(AUC)
}


## Computation AUC binary ------
res1.2d=multiclass.auc1(predicted.2d, truth.2d)
res3.2d=multiclass.auc3(predicted.2d, truth.2d)


## Computation AUC ----

#AUC1
res1=multiclass.auc1(predicted, truth)

#AUC2
res2=multiclass.auc2(predicted, truth)

#AUC3
res3=multiclass.auc3(predicted, truth)


