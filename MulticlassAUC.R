multiclass.auc2 = function(pred, truth){
  if(nlevels(truth) == 2) {pred = cbind(pred,1-pred)}
  predP = pred
  # choose the probablity of the truth
  predV = vnapply(seq_row(pred), function(i) {
    pred[i, resp[i]]
  })
  auc = pROC::multiclass.roc(response = truth, predictor = predV)$auc
  as.numeric(auc)
}

library(mlr)

lrn = makeLearner("classif.randomForest", predict.type = "prob")
model = train(lrn, iris.task)
pred = predict(model, iris.task)

library(caTools)
colAUC(getPredictionProbabilities(pred), pred$data$truth)
