rm(list = ls())



## Libraries

## Get the functions
source("Tests/ScoredAUCImplementations.R")

## ========================
## Example I : iris dataset 
## ========================

#### Part 1

# Data Iris ----
learner <- makeLearner('classif.lda', predict.type="prob")
task <- makeClassifTask(data=iris, target="Species")
mod <- train(learner=learner, task=task)
pred.obj <- predict(mod, newdata=iris)
predicted <- as.matrix(pred.obj$data[,paste("prob.", levels(pred.obj$data$response), sep="")])
colnames(predicted)<-levels(pred.obj$data$response)
truth=iris$Species
pred = predicted


# Run the sAUC ----
binaryclass.Scoredauc(predicted, truth, names = TRUE) # fast implementation
binaryclass.Scoredauc.naive(predicted, truth, names = TRUE) # simple implementation, like definition
binaryclass.Scoredauc.naive.permutations(predicted, truth, names = TRUE) # Matrix products


#### Part 2 : verification of values -----

# sAUC.coeff is coeff computed with explicit formula, m comes from fast computation

m<-binaryclass.Scoredauc(predicted, truth, names = TRUE)
m.naive<-binaryclass.Scoredauc.naive.permutations(predicted, truth, names = TRUE) 
m.final<-binaryclass.Scoredauc.final(predicted, truth, names = TRUE)
m.naive-m

j=3

# c1=1 c=2
sAUC1 = sAUC.coeff(pred, truth, c1=1, c2=2)
m[1,3]
sAUC1

# c1=2 c2=1
sAUC2 = sAUC.coeff(pred, truth, c1=2, c2=1)
m[6,3]
sAUC2

# c1=1 c2=3
sAUC3 = sAUC.coeff(pred, truth, c1=1, c2=3)
m[2,3]
sAUC3

# c1=3 c2=2
sAUC4 = sAUC.coeff(pred, truth, c1=3, c2=1)
m[5,3]
sAUC4

# c1=2 c2=3
sAUC5 = sAUC.coeff(pred, truth,c1=2, c2=3)
m[3,3]
sAUC5

# c1=3 c2=2
sAUC6 = sAUC.coeff(pred, truth, c1=3, c2=2)
m[4,3]
sAUC6

# Difference
m[1,3]-sAUC1
m[6,3]-sAUC2
m[2,3]-sAUC3
m[5,3]-sAUC4
m[3,3]-sAUC5
m[4,3]-sAUC6

m.naive[1,3]-sAUC1
m.naive[6,3]-sAUC2
m.naive[2,3]-sAUC3
m.naive[5,3]-sAUC4
m.naive[3,3]-sAUC5
m.naive[4,3]-sAUC6

# Still some 1e-15 unexplained difference...





#### Part 3 : verification of formula


## test for 1-2----
# Here we have j=3, thus probabilitie values according to virginica prediction
# We compare the Setosa vs versicolor

number=50*50

# Part 2.a : Virginica vs Setosa, c1=3, c2=1, sum of probabilities difference ----
predVirginicapart1=pred[1:50,3]
predvirginicapart2=pred[51:100,3]

sum=0
for (i in c(1:50)) {
  for (j in c(1:50)) {
    sum=sum+predvirginicapart2[j]-predVirginicapart1[i]
  }
}
sum.diffproba.12=sum/number
sum.diffproba.12

## test for 1-3 ----
# Here we have j=3, thus probabilitie values according to virginica prediction
# We compare the Setosa vsvirginica


# Part 2.b : Virginica vs Setosa, c1=3, c2=1, sum of probabilities difference ----
predVirginicapart1=pred[1:50,3]
predvirginicapart3=pred[101:150,3]

sum=0
for (i in c(1:50)) {
  for (j in c(1:50)) {
    sum=sum+predvirginicapart3[j]-predVirginicapart1[i]
  }
}
sum.diffproba.13=sum/number
sum.diffproba.13

## test for 2-3 ----
# Here we have j=3, thus probabilitie values according to virginica prediction
# We compare the versicolor vs virginica


# Part 2.c : Virginica vs Setosa, c1=3, c2=1 sum of probabilities difference ----
predVirginicapart2=pred[51:100,3]
predvirginicapart3=pred[101:150,3]

sum=0
for (i in c(1:50)) {
  for (j in c(1:50)) {
    sum=sum+predvirginicapart3[j]-predVirginicapart2[i]
  }
}
sum.diffproba.23=sum/number
sum.diffproba.23



## Verification ----

m[6,3]-(m[1,3]+sum.diffproba.12)

m[5,3]-(m[2,3]+sum.diffproba.13)

m[4,3]-(m[3,3]+sum.diffproba.23)



## ==========================================
## Example II : 2 class example
## ==========================================

predicted=c(0.1,0.12,0.14,0.26,0.45,0.65,0.74,0.8,0.85)
truth=as.factor(c(0,0,0,1,0,1,0,1,1))

# Run the sAUC ----
binaryclass.Scoredauc(predicted, truth, names = TRUE) # fast implementation
binaryclass.Scoredauc.naive.permutations(predicted, truth, names = TRUE) # Matrix products
binaryclass.Scoredauc.final(predicted, truth, names = TRUE)

# Computing the values of the the matrix with definition ----
sAUC.coeff(predicted, truth, c1=1, c2=2, j=1)
sAUC.coeff(predicted, truth, c1=2, c2=1, j=1)

# Computing the sum of crossed difference
pred0 = pred[which(truth==0)]
pred1 = pred[which(truth==1)]
number=length(pred0)*length(pred1)

sum=0
for (i in c(1:length(pred0))) {
  for (j in c(1:length(pred1))) {
    sum=sum+pred0[i]-pred1[j]
  }
}
sum.diff=sum/number
sum.diff

sAUC.coeff(predicted, truth, c1=2, c2=1, j=1)
sAUC.coeff(predicted, truth, c1=1, c2=2, j=1)
sum.diff

sAUC.coeff(predicted, truth, c1=2, c2=1, j=1)-sAUC.coeff(predicted, truth, c1=1, c2=2, j=1)+sum.diff

## ========================================================
## Example III : case where some probabilities are the same
## ========================================================
# But thus with the definition it is simplier that for basic AUC ??? let's see...


predicted=c(0.1,0.2,0.2,0.2,0.2,0.2,0.6,0.7,0.7,0.7,0.7,0.7)
truth=as.factor(c(0,0,0,1,1,0,1,1,1,0,1,1))

# Run the sAUC ----
binaryclass.Scoredauc(predicted, truth, names = TRUE) # fast implementation
binaryclass.Scoredauc.naive.permutations(predicted, truth, names = TRUE) # Matrix products
binaryclass.Scoredauc.final(predicted, truth, names = TRUE)

# Computing the values of the the matrix with definition ----
sAUC.coeff(predicted, truth, c1=1, c2=2, j=1)
sAUC.coeff(predicted, truth, c1=2, c2=1, j=1)

