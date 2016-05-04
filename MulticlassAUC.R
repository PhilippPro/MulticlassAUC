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

library(mlr)

lrn = makeLearner("classif.lda", predict.type = "prob")
model = train(lrn, iris.task)
pred = predict(model, iris.task)

library(caTools)
m = colAUC(getPredictionProbabilities(pred), pred$data$truth)
# First column: j setosa
colAUC(getPredictionProbabilities(pred)[,1], pred$data$truth)

# Comparison of colAUC with pROC
library(pROC)
# Syntax (response, predictor):
result_pROC = auc(droplevels(pred$data$truth[pred$data$truth %in% c("versicolor", "virginica")]), getPredictionProbabilities(pred)[,2][pred$data$truth %in% c("versicolor", "virginica")])
result_pROC == m[3,2]
result_pROC = auc(droplevels(pred$data$truth[pred$data$truth %in% c("setosa", "virginica")]), getPredictionProbabilities(pred)[,3][pred$data$truth %in% c("setosa", "virginica")])
result_pROC == m[2,3]
# Seems to be reliable!

multiclass.auc3(getPredictionProbabilities(pred), pred$data$truth)


truth = pred$data$truth
probabilities = getPredictionProbabilities(pred)
# AUNU
mean(vnapply(1:nlevels(truth), function(i) colAUC(getPredictionProbabilities(pred)[,i], truth == levels(truth)[i])))

# AUNP
sum(vnapply(1:nlevels(truth), function(i) mean(truth == levels(truth)[i]) * colAUC(getPredictionProbabilities(pred)[,i], truth == levels(truth)[i])))

# AU1U (Hand and Till)
m = colAUC(getPredictionProbabilities(pred), truth)
mean(m[col(m) != ncol(m) + 1 - row(m)])
# or
c = t(combn(1:nlevels(truth), 2, simplify = T))
mean(m[cbind(1:nrow(c), c[, 2])])

# AU1P
m = colAUC(getPredictionProbabilities(pred), truth)
weights = table(truth[-1]) / length(truth[-1])
m = m * matrix(rep(weights, each = length(weights)), ncol = length(weights))
sum(m[col(m) != ncol(m) + 1 - row(m)]) / (ncol(m) - 1)

# ensure that it works, also when no observation is present!!