library(plsRcox)
library(survival)
library(survivalROC)
library(Hmisc)
library(survcomp)

plsR_cox <- function(x.train,fold){
  set.seed(4321)
  cv.plsRcox.res=cv.plsRcox(list(x=x.train[,3:ncol(x.train)],
                                 time=x.train$time,
                                 status=x.train$status),
                            nfold = 5,
                            nt=10,
                            verbose = FALSE)
  
  
  fit <- plsRcox(x.train[,3:ncol(x.train)],
                 time=x.train$time,
                 event=x.train$status,
                 nt=as.numeric(cv.plsRcox.res[5]))
  
  pred_pls <- predict(fit,type="lp",newdata=x.test[,-c(1,2)])
  pred_pls <- as.numeric(pred_pls)
  cindex_pls = 1-rcorr.cens(pred_pls,Surv(t.test, s.test))[[1]]
  print(sprintf("plscox_%d",fold))
  print(cindex_pls)
  
  result <- list()
  result$best_param =list(nt=as.numeric(cv.plsRcox.res[5]))
  result$model = fit
  result$cindex = cindex_pls
  result$pred = pred_pls 
  
  
  ### auc km ###
  
  cutoff = 12*1
  if ( min(x.test$time) < cutoff )
  {
    y <- survivalROC(Stime = t.test, status = s.test, marker = pred_pls,
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
  y <- survivalROC(Stime = t.test, status = s.test, marker = pred_pls, 
                   predict.time = cutoff,method = "KM")
  result$km_fp_3 = y$FP
  result$km_tp_3 = y$TP
  result$km_auc_3 = y$AUC
  
  cutoff=12*5
  y <- survivalROC(Stime = t.test, status = s.test, marker = pred_pls, 
                   predict.time = cutoff,method = "KM")
  result$km_fp_5 = y$FP
  result$km_tp_5 = y$TP
  result$km_auc_5 = y$AUC
  

  dd_ext <- data.frame("time"=x.test$time, "event"=x.test$status, "score"= pred_pls)
  Brier_score <- survcomp::sbrier.score2proba(data.tr=dd_ext, data.ts=dd_ext, method="cox")
  result$brier_Score <- Brier_score
  
  save("result", file = sprintf("%d_plsRcox_result.RData", fold))
  
  
}

