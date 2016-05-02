
#' @export multiclass.auc
#' @rdname measures
#' @format none
multiclass.auc = makeMeasure(id = "multiclass.auc", minimize = FALSE, best = 1, worst = 0,
                             properties = c("classif", "classif.multi", "req.pred", "req.truth", "req.prob"),
                             name = "Multiclass area under the curve",
                             note = "Calls `pROC::multiclass.roc`.",
                             fun = function(task, model, pred, feats, extra.args) {
                               requirePackages("pROC", why = "multiclass.auc", default.method = "load")
                               resp = pred$data$response
                               predP = getPredictionProbabilities(pred)
                               # choose the probablity of the choosen response
                               predV = vnapply(seq_row(predP), function(i) {
                                 predP[i, resp[i]]
                               })
                               auc = pROC::multiclass.roc(response = resp, predictor = predV)$auc
                               as.numeric(auc)
                             }
)

###############################################################################
### classif binary ###
###############################################################################
#' @export auc
#' @rdname measures
#' @format none
auc = makeMeasure(id = "auc", minimize = FALSE, best = 1, worst = 0,
                  properties = c("classif", "req.pred", "req.truth", "req.prob"),
                  name = "Area under the curve",
                  fun = function(task, model, pred, feats, extra.args) {
                    # ROCR does not work with NAs
                    if (anyMissing(pred$data$response) || length(unique(pred$data$truth)) == 1L)
                      return(NA_real_)
                    measureAUC(getPredictionProbabilities(pred), pred$data$truth, pred$task.desc$negative, pred$task.desc$positive)
                  }
)

#' @export measureAUC
#' @rdname measures
#' @format none
measureAUC = function(probabilites, truth, negative, positive) {
  rpreds = asROCRPredictionIntern(probabilites, truth, negative, positive)
  ROCR::performance(rpreds, "auc")@y.values[[1L]]
}
