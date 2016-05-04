library(mlr)

lrn = makeLearner("classif.lda", predict.type = "prob")
model = train(lrn, iris.task)
pred = predict(model, iris.task)

truth = pred$data$truth
probabilities = getPredictionProbabilities(pred)

library(caTools)
m = colAUC(probabilities, pred$data$truth)
# First column: j setosa
colAUC(probabilities[,1], pred$data$truth)

# Comparison of colAUC with pROC
library(pROC)
# Syntax (response, predictor):
result_pROC = auc(droplevels(pred$data$truth[pred$data$truth %in% c("versicolor", "virginica")]), probabilities[,2][pred$data$truth %in% c("versicolor", "virginica")])
result_pROC == m[3,2]
result_pROC = auc(droplevels(pred$data$truth[pred$data$truth %in% c("setosa", "virginica")]), probabilities[,3][pred$data$truth %in% c("setosa", "virginica")])
result_pROC == m[2,3]
# Seems to be reliable!

# multiclass.auc3(getPredictionProbabilities(pred), pred$data$truth)

# AUNU
mean(vnapply(1:nlevels(truth), function(i) colAUC(probabilities[,i], truth == levels(truth)[i])))

# AUNP
sum(vnapply(1:nlevels(truth), function(i) mean(truth == levels(truth)[i]) * colAUC(probabilities[,i], truth == levels(truth)[i])))

# AU1U (Hand and Till)
m = colAUC(probabilities, truth)
c = c(combn(1:nlevels(truth), 2))
mean(m[cbind(rep(1:nrow(m), each = 2), c)])

# AU1P
m = colAUC(probabilities, truth)
weights = table(truth) / length(truth)
m = m * matrix(rep(weights, each = nrow(m)), ncol = length(weights))
sum(m[cbind(rep(1:nrow(m), each = 2), c)]) / (nlevels(truth) - 1)

# ensure that it works, also when no observation is present!!
