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
m = colAUC(getPredictionProbabilities(pred)[,1], pred$data$truth)
# First column: j setosa, 
colAUC(getPredictionProbabilities(pred)[,1], pred$data$truth)

# Is this Hand and Till?
mean(m[col(m)!=4-row(m)])
colMeans(m)

m = colAUC(getPredictionProbabilities(pred), pred$data$truth)
c = t(combn(1:nlevels(pred$data$truth), 2, simplify = T))
m[cbind(1:nrow(c),c[,2])]
mean(m[cbind(1:nrow(c),c[,2])])

# AUNU
for(j in 1:nlevels(pred$data$truth)){
  
}

