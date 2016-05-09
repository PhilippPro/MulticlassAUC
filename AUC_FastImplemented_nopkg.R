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
        }
      }
      Auc[i, j] = sum/number # En gros le nombre de non variations de c2 par rapport a c1 qui est constant
    }
  }
  
  Auc<-apply(Auc,c(1,2),function(x) return(max(x,1-x)))
  # Add the names
  if (names==TRUE) {
  rownames(Auc) = paste(levels.values[permutations[, 1]], " vs. ", levels.values[permutations[, 2]], 
                        sep = "")
  colnames(Auc) = colnames(X)
  }
  print(proc.time()-ptm)
  return(Auc)
}

i=3
#pred=pred[,i]
#truth=as.factor(truth == levels(truth)[i])

binaryclass.auc(pred, truth, names = TRUE)
colAUC(pred, truth)

i=3
pred=probabilities[,i]
truth=truth == levels(truth)[i]
binaryclass.auc(pred, truth)

## AUNU
multiclass.aunu= function(pred, truth,  names=FALSE) {
  mean(vnapply(1:nlevels(truth), function(i) binaryclass.auc(probabilities[,i], truth == levels(truth)[i])))
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
  y<-truth
  X<-pred
  K=dim(pred)[2] #maj k
  n=dim(pred)[1]
  levels.values=as.factor(levels(truth))
  levels.names=levels(truth)
  L = matrix(rep(levels.values, each = n), n, K)
  permutations = combs(1:K, 2)
  nP = nrow(permutations)
  Auc = matrix(0.5, K, K) # Creates the AUC matrixc
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
      
      for (p in c(1:n)) {
        if (d2[p,c1]==1) {
          index=which(d2[c(1:p),c2]==1)
          if (length(index)>1) {
            for (q in c(1:length(index)))
              indexq=index[q]
            sum=sum+pred[p,c1]-pred[q,c2]
            count=count+1
          }
        }
      }
      Auc[i, j] = sum/number # En gros le nombre de non variations de c2 par rapport a c1 qui est constant
    }
  }
  
  # Add the names
  if (names==TRUE) {
    colnames(Auc)<-levels.names
    rownames(Auc)<-levels.names
    
    Auc<-apply(Auc,c(1,2),function(x) return(max(x,1-x)))
    rownames(Auc) = paste(levels.values[permutations[, 1]], " vs. ", levels.values[permutations[, 2]], 
                          sep = "")
    colnames(Auc) = colnames(X)
    
  }
  
  print(proc.time()-ptm)
  print(paste("count",count))
  return(Auc)
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




# tend to confirm the 2 first lines are in fact 1 because no problem between
i=1
test=cbind(predicted[,i],as.factor(iris$Species))
test[order(test[,1]),]
