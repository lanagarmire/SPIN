
rm(list = ls())
suppressPackageStartupMessages({
  library(e1071)
  library(glmnet)
  library(pROC)
  library(readxl)
})

# ---------- data preparation ----------
setwd("D:\\Research_project\\For Test\\Elastic_net_Lasso\\Data")
Y_ls = readRDS("hdWGCNA_co_expression_list_uni.RDS")
meta = read_excel("meta_uni.xlsx")
cov_df = meta[,2]
y = meta$Cognitive_Status


compute_metrics <- function(y_true, prob) {
  y_true <- as.integer(y_true)
  prob   <- as.numeric(prob)
  label  <- as.integer(prob > 0.5)
  
  TP <- sum(y_true == 1 & label == 1)
  TN <- sum(y_true == 0 & label == 0)
  FP <- sum(y_true == 0 & label == 1)
  FN <- sum(y_true == 1 & label == 0)
  safe_div <- function(a, b) if (b == 0) NA_real_ else a / b
  
  acc  <- safe_div(TP + TN, TP + TN + FP + FN)
  prec <- safe_div(TP, TP + FP)            
  rec  <- safe_div(TP, TP + FN)             
  spec <- safe_div(TN, TN + FP)        
  f1   <- if (is.na(prec) || is.na(rec) || (prec + rec) == 0) NA_real_ else 2 * prec * rec / (prec + rec)
  
  roc_obj <- tryCatch(pROC::roc(y_true, prob, quiet = TRUE), error = function(e) NULL)
  auc_val <- if (is.null(roc_obj)) NA_real_ else as.numeric(pROC::auc(roc_obj))
  
  list(
    confusion = c(TP = TP, TN = TN, FP = FP, FN = FN),
    accuracy = acc, precision = prec, recall = rec, specificity = spec, F1 = f1, auc = auc_val,
    label = label, prob = prob
  )
}

prep_edges <- function(Y_ls) {
  stopifnot(length(Y_ls) > 0)
  N <- nrow(Y_ls[[1]])
  lt <- lower.tri(matrix(0, N, N), diag = FALSE)
  p  <- sum(lt)
  edge_mat <- t(vapply(Y_ls, function(A) A[lt], numeric(p)))
  colnames(edge_mat) <- paste0("e", seq_len(p))
  edge_mat
}

prep_cov <- function(cov_df) {
  model.matrix(~ . - 1, data = cov_df)
}

svm_rank_edges <- function(Xedge, y, cost = 1) {
  fit <- svm(x = Xedge, y = as.factor(y),
             kernel = "linear", type = "C-classification",
             cost = cost, scale = TRUE)
  w <- drop(t(fit$coefs) %*% fit$SV)  
  order(abs(w), decreasing = TRUE)
}

# ---------- SVM  ----------
tune_svm_topk <- function(Xedge, Xcov, y,
                          topk_grid = c(10, 15, 20, 25, 30),
                          cost_grid = c(1),
                          repeats = 5, val_frac = 0.2, seed = 1) {
  set.seed(seed)
  n <- length(y)
  best <- list(score = -Inf, top_k = topk_grid[1], cost = cost_grid[1])
  for (Cval in cost_grid) {
    for (k in topk_grid) {
      aucs <- numeric(repeats)
      for (r in seq_len(repeats)) {
        idx <- sample.int(n)
        nval <- max(1, round(n * val_frac))
        val_idx <- idx[seq_len(nval)]
        tr_idx  <- idx[-seq_len(nval)]
        
        ord  <- svm_rank_edges(Xedge[tr_idx, , drop = FALSE], y[tr_idx], cost = Cval)
        take <- ord[seq_len(min(k, length(ord)))]
        
        df_tr  <- data.frame(y = y[tr_idx],
                             cbind(Xcov[tr_idx, , drop = FALSE], Xedge[tr_idx, take, drop = FALSE]))
        df_val <- data.frame(y = y[val_idx],
                             cbind(Xcov[val_idx, , drop = FALSE], Xedge[val_idx, take, drop = FALSE]))
        fit <- glm(y ~ ., data = df_tr, family = "binomial")
        prob <- as.vector(predict(fit, newdata = df_val, type = "response"))
        ro <- tryCatch(pROC::roc(df_val$y, prob, quiet = TRUE), error = function(e) NULL)
        aucs[r] <- if (is.null(ro)) NA_real_ else as.numeric(pROC::auc(ro))
      }
      mauc <- mean(aucs, na.rm = TRUE)
      if (!is.nan(mauc) && mauc > best$score) best <- list(score = mauc, top_k = k, cost = Cval)
    }
  }
  best[c("top_k", "cost")]
}

# ---------- Lasso ----------
fit_lasso_cv <- function(Xcov, Xedge, y, nfolds = 5, use_1se = FALSE) {
  X  <- cbind(Xcov, Xedge)
  pf <- c(rep(0, ncol(Xcov)), rep(1, ncol(Xedge)))
  cvfit <- cv.glmnet(X, y, family = "binomial",
                     alpha = 1, penalty.factor = pf,
                     nfolds = nfolds, standardize = TRUE)
  lam <- if (use_1se) cvfit$lambda.1se else cvfit$lambda.min
  glmnet(X, y, family = "binomial",
         alpha = 1, lambda = lam,
         penalty.factor = pf, standardize = TRUE)
}

# ---------- Elastic Net ----------
fit_enet_cv <- function(Xcov, Xedge, y,
                        alpha_grid = c(0.2, 0.5, 0.8),
                        nfolds = 5, use_1se = FALSE) {
  X  <- cbind(Xcov, Xedge)
  pf <- c(rep(0, ncol(Xcov)), rep(1, ncol(Xedge)))  ## 临床项不惩罚
  best <- list(score = Inf, alpha = alpha_grid[1], lambda = NA)
  for (a in alpha_grid) {
    cvfit <- cv.glmnet(X, y, family = "binomial",
                       alpha = a, penalty.factor = pf,
                       nfolds = nfolds, standardize = TRUE)
    cur_score <- min(cvfit$cvm)
    cur_lam   <- if (use_1se) cvfit$lambda.1se else cvfit$lambda.min
    if (cur_score < best$score) best <- list(score = cur_score, alpha = a, lambda = cur_lam)
  }
  fit <- glmnet(X, y, family = "binomial",
                alpha = best$alpha, lambda = best$lambda,
                penalty.factor = pf, standardize = TRUE)
  list(fit = fit, alpha = best$alpha, lambda = best$lambda)
}

# ---------- loocv ----------
nested_loocv_run <- function(Y_ls, cov_df, y,
                             topk_grid = c(10, 15, 20, 25, 30),
                             cost_grid = c(1),                
                             enet_alpha_grid = c(0.2, 0.5, 0.8),
                             use_1se = FALSE,
                             inner_repeats = 5, inner_val_frac = 0.2,
                             outer_seed = 2025) {
  set.seed(outer_seed)
  M <- length(Y_ls); stopifnot(M == length(y), nrow(cov_df) == M)
  
  Xedge_full <- prep_edges(Y_ls)
  Xcov_full  <- prep_cov(cov_df)
  
  prob_svm  <- numeric(M)
  prob_las  <- numeric(M)
  prob_enet <- numeric(M)
  
  for (i in seq_len(M)) {
    tr <- setdiff(seq_len(M), i); te <- i
    
    Xcov_tr <- Xcov_full[tr, , drop = FALSE]
    Xcov_te <- Xcov_full[te, , drop = FALSE]
    Xedg_tr <- Xedge_full[tr, , drop = FALSE]
    Xedg_te <- Xedge_full[te, , drop = FALSE]
    y_tr <- y[tr]; y_te <- y[te]
    
    # ---- A) SVM + Logistic  ----
    ksel <- tune_svm_topk(Xedg_tr, Xcov_tr, y_tr,
                          topk_grid = topk_grid, cost_grid = cost_grid,
                          repeats = inner_repeats, val_frac = inner_val_frac, seed = outer_seed + i)
    ord  <- svm_rank_edges(Xedg_tr, y_tr, cost = ksel$cost)
    take <- ord[seq_len(min(ksel$top_k, length(ord)))]
    df_tr <- data.frame(y = y_tr, cbind(Xcov_tr, Xedg_tr[, take, drop = FALSE]))
    df_te <- data.frame(        cbind(Xcov_te, Xedg_te[, take, drop = FALSE]))
    fit_logit <- glm(y ~ ., data = df_tr, family = "binomial")
    prob_svm[i] <- as.vector(predict(fit_logit, newdata = df_te, type = "response"))
    
    # ---- B) Lasso ----
    fit_las <- fit_lasso_cv(Xcov_tr, Xedg_tr, y_tr, nfolds = 5, use_1se = use_1se)
    prob_las[i] <- as.vector(predict(fit_las, newx = cbind(Xcov_te, Xedg_te), type = "response"))
    
    # ---- C) Elastic Net  ----
    en <- fit_enet_cv(Xcov_tr, Xedg_tr, y_tr,
                      alpha_grid = enet_alpha_grid, nfolds = 5, use_1se = use_1se)
    prob_enet[i] <- as.vector(predict(en$fit, newx = cbind(Xcov_te, Xedg_te), type = "response"))
  }
  
  res_svm  <- compute_metrics(y, prob_svm)
  res_las  <- compute_metrics(y, prob_las)
  res_enet <- compute_metrics(y, prob_enet)
  
  list(svm = res_svm, lasso = res_las, enet = res_enet)
}

print_summary <- function(res, name) {
  cat("\n====", name, "====\n")
  cat("AUC        :", sprintf("%.3f", res$auc), "\n")
  cat("Accuracy   :", sprintf("%.3f", res$accuracy), "\n")
  cat("Precision  :", sprintf("%.3f", res$precision), "\n")
  cat("Recall     :", sprintf("%.3f", res$recall), "\n")
  cat("Specificity:", sprintf("%.3f", res$specificity), "\n")
  cat("F1         :", sprintf("%.3f", res$F1), "\n")
  cat("Confusion  : TP=", res$confusion["TP"], " TN=", res$confusion["TN"],
      " FP=", res$confusion["FP"], " FN=", res$confusion["FN"], "\n")
}

 res <- nested_loocv_run(
   Y_ls, cov_df, y,
   topk_grid = seq(5,20),
   cost_grid = c(1),                 
   enet_alpha_grid = seq(0.1,0.9,0.05),
   use_1se = FALSE,                  
   inner_repeats = 5, inner_val_frac = 0.2,
   outer_seed = 2025
 )
setwd("D:\\Research_project\\For Test\\Elastic_net_Lasso\\For_benchmark")
saveRDS(res,file = "benchmark_res_uni.RDS")


