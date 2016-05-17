library(mlr)

lrn = makeLearner("classif.lda", predict.type = "prob")
model = train(lrn, sonar.task)
pred = predict(model, sonar.task)

performance(pred, measures = list(multiclass.aunu, multiclass.aunp, multiclass.au1u, multiclass.au1p))

# Multiclass Example
library(OpenML)
data = getOMLTask(task.id = 3548, verbosity=0)$input$data.set

multiclass.task = makeClassifTask(id = "abc", data = data$data, target = "Type")
model = train(lrn, multiclass.task)
pred = predict(model, multiclass.task)
performance(pred, measures = list(multiclass.aunu, multiclass.aunp, multiclass.au1u, multiclass.au1p, logloss, multiclass.brier))

probabilities = getPredictionProbabilities(pred)

library(caTools)
m = colAUC(probabilities, pred$data$truth)
# First column: j 
colAUC(probabilities[,1], pred$data$truth)

truth = pred$data$truth

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
c = c(combn(1:nlevels(truth), 2))
weights = table(truth) / length(truth)
m = m * matrix(rep(weights, each = nrow(m)), ncol = length(weights))
sum(m[cbind(rep(1:nrow(m), each = 2), c)]) / (nlevels(truth) - 1)

# ensure that it works, also when no observation is present and in higher class case!