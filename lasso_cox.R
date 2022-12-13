library(glmnet)
library(Hmisc)
library(survcomp)
library(survivalROC)
library(survival)
lasso_cox <- function(x.train,fold){

time = as.double(x.train$time)
status = as.double(x.train$status)

set.seed(4321)
fit = cv.glmnet(as.matrix(x.train[,3:ncol(x.train)]),as.matrix(x.train[,1:2]) , 
                family = "cox",alpha=1,nfolds = 5)

pred_cox = predict(fit,type='link',
                   newx=as.matrix(x.test[,3:ncol(x.test)]),
                   s = fit$lambda.min)
pred_cox <- as.numeric(pred_cox)
cindex_cox = 1-rcorr.cens(pred_cox,Surv(t.test, s.test))[[1]]
print(sprintf("cox_%d",fold))
print(cindex_cox)


result <- list()
result$best_param <- list( s = fit$lambda.min)
result$model = fit
result$cindex = 1-rcorr.cens(pred_cox,Surv(t.test, s.test))[[1]]
result$pred = pred_cox 

### auc km ###

cutoff = 12*1
if ( min(x.test$time) < cutoff )
{
  y <- survivalROC(Stime = t.test, status = s.test, marker = pred_cox,
                   predict.time = cutoff, method = "KM")
  result$km_fp_1 = y$FP
  result$km_tp_1 = y$TP
  result$km_auc_1 = y$AUC
}else
{
  result$km_fp_1 = NA
  result$km_tp_1 = NA
  result$km_auc_1 = NA
}

cutoff=12*3
y <- survivalROC(Stime = t.test, status = s.test, marker = pred_cox, 
                 predict.time = cutoff,method = "KM")
result$km_fp_3 = y$FP
result$km_tp_3 = y$TP
result$km_auc_3 = y$AUC

cutoff=12*5
y <- survivalROC(Stime = t.test, status = s.test, marker = pred_cox, 
                 predict.time = cutoff,method = "KM")
result$km_fp_5 = y$FP
result$km_tp_5 = y$TP
result$km_auc_5 = y$AUC

dd_ext <- data.frame("time"=x.test$time, "event"=x.test$status, "score"= pred_cox)
Brier_score <- sbrier.score2proba(data.tr=dd_ext, data.ts=dd_ext, method="cox")
result$brier_Score <- Brier_score


save("result", file = sprintf("%d_lasso_cox_result.RData", fold))

}




