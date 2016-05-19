rm(list = ls())



## Libraries

## Get the functions
source("Tests/ProbabilisticAUCImplementations.R")

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


# Run the pAUC ----
binaryclass.Probabilisticauc(pred, truth)
