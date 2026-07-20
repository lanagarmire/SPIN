#!/usr/bin/env Rscript

.libPaths(c("~/R_libs_amd", .libPaths()))
Sys.setenv(PKG_CXXFLAGS="-O3 -march=native -mtune=native")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(coda)
  library(Rcpp)
  library(RcppArmadillo)
  library(readxl)
})

# Load BSNMani Utils & Diagnostics 
base_src_dir <- "/home/koe/BSNMani_application-main/Github_BSNMani/SEA-AD/BSNMani/SEA-AD/BSNMani-dev"
source("/home/koe/BSNMani_application-main/Github_BSNMani/IMC/Preprocessed/Utils.R")
source(file.path(base_src_dir, "g1_diagnostics_helper.R"))

# NEW: Compile the C++ functions (including polar_expansion)
sourceCpp(file.path(base_src_dir, "hybrid_M0_MALA_LR_FAST_v2.cpp"))
sourceCpp(file.path(base_src_dir, "hybrid_M0_MALA_LR_g2.cpp"))

# ==============================================================================
# 1. OUT-OF-SAMPLE PROJECTION FUNCTION
# ==============================================================================
project_lambda_test <- function(Y_test, U_train, s2_train, t2_lambda_train) {
  d_UtYU <- diag(t(U_train) %*% Y_test %*% U_train)
  shrinkage_factor <- 1 / ((1 / t2_lambda_train) + (1 / s2_train))
  lambda_test <- shrinkage_factor * (d_UtYU / s2_train)
  return(as.numeric(lambda_test))
}

# ==============================================================================
# 2. CONFIGURATION & LOAD DATA
# ==============================================================================
## check point                     ####
# cell_type = "WGCNA"
# cell_type = "Smoothie"
cell_type = "Oli"
# cell_type = "hdWGCNA"
# cell_type = "SpaceX"

q_candidates <- 2:5
r_val <- 2 # Intercept + Atherosclerosis
n_chains <- 2

cat("========================================================\n")
cat(" Compiling Unbiased LOOCV Results & Computing Metrics\n")
cat("========================================================\n\n")

# Load BSNMani Utils & Diagnostics 
base_src_dir <- "/home/koe/BSNMani_application-main/Github_BSNMani/SEA-AD/BSNMani/SEA-AD/BSNMani-dev"
source("/home/koe/BSNMani_application-main/Github_BSNMani/IMC/Preprocessed/Utils.R")
source(file.path(base_src_dir, "g1_diagnostics_helper.R"))
## check point                     ####
# Load Clinical & Network Data
setwd("/home/koe/BSNMani_application-main/Github_BSNMani/SEA-AD")
# clinical_df = readRDS("SEA-AD_Data/filtered_baseline_meta.rds")
clinical_df = read_excel("SEA-AD_Data/meta.xlsx")

# FIX 2: Safe conversion from factor to numeric
clinical_df$Cognitive_Status <- as.numeric(as.character(clinical_df$Cognitive_Status))

cov_df <- as.matrix(clinical_df %>% select(Atherosclerosis))
## check point                     ####
# wgcna_transformed_ls = readRDS("SEA-AD_Data/Filtered_baseline_co_expression_list.rds")
# wgcna_transformed_ls = readRDS("SEA-AD_Data/Filtered_Smoothie_co_expression_list.rds")
wgcna_transformed_ls = readRDS("SEA-AD_Data/Smoothie_co_expression_list_Oli.RDS")
# wgcna_transformed_ls = readRDS("SEA-AD_Data/hdWGCNA_co_expression_list.RDS")
# wgcna_transformed_ls = readRDS("SEA-AD_Data/Shared_co_expression_list.RDS")

n <- nrow(clinical_df)
n_roi <- nrow(wgcna_transformed_ls[[1]])

FC_arr <- array(NA, dim=c(n_roi, n_roi, n))
for(i in 1:n){
  tmp_mat <- as.matrix(wgcna_transformed_ls[[i]])
  diag(tmp_mat) <- 0 # FIX: Zero out the diagonal to match your training script!
  FC_arr[,,i] <- tmp_mat
}

# Tracking Variables
test_predicted_probs <- numeric(n)
test_predicted_labels <- numeric(n)
test_true_labels <- numeric(n)
selected_qs <- numeric(n)
dic_log <- matrix(0, nrow = n, ncol = length(q_candidates))
colnames(dic_log) <- paste0("q=", q_candidates)

# ==============================================================================
# 3. COMPILE RESULTS ACROSS FOLDS
# ==============================================================================
for (i in 1:n) {
  cat(sprintf("Processing Fold %02d / %02d...\n", i, n))
  
  train_idx <- setdiff(1:n, i)
  FC_arr_train <- FC_arr[,,train_idx]
  clinical_df_train <- clinical_df[train_idx, ]
  cov_df_train <- cov_df[train_idx, , drop=FALSE]
  
  FC_test <- FC_arr[,,i]
  clinical_test <- clinical_df[i, ]
  
  best_q <- NA
  best_dic <- Inf
  best_model_params <- list()
  
  fold_str <- sprintf("Fold_%02d", i)
  fold_dir <- fs::path("/home/koe/BSNMani_application-main/Github_BSNMani/result", cell_type, fold_str)
  
  # A. Evaluate all candidate models for this fold
  for (q in q_candidates) {
    hybrid_res <- list()
    
    # --- UPDATED: Loop over both chains ---
    for (c_idx in 1:n_chains) {
      mh_file <- fs::path(fold_dir, paste("q", q, sep="_"), "MH", paste0("MH_res_GQN_KPN_chain_", c_idx, "_q_", q, ".RDS"))
      if(!file.exists(mh_file)) {
        stop(sprintf("ERROR: Missing MH file for Fold %d, q = %d, chain = %d.", i, q, c_idx))
      }
      
      # Safely Un-nest the MCMC object to prevent the dim(X) apply error ---
      tmp_res <- readRDS(mh_file)
      # Safely Un-nest the MCMC object
      if ("mcmc" %in% names(tmp_res)) {
        hybrid_res[[c_idx]] <- list(mcmc = tmp_res$mcmc)
      } else {
        hybrid_res[[c_idx]] <- list(mcmc = tmp_res)
      }
    }
    
    # 2. Extract Posterior Means using ALL chains
    dat_ls <- list(N = n_roi, q = q, n_chains = n_chains, n_samps = length(hybrid_res[[1]]$mcmc$s2), n_burnin = 30000)
    
    var_df_full <- data.frame(
      var = c("Lambda_flat","s2","t2_lambda","X","U","Y_C_llk","Y_llk","C_llk","d","t2_alpha","t2_beta","t2"),
      dim = c(3,1,1,3,3,1,1,1,2,1,1,1)
    )
    
    MALA_samps_ls <- gather_samples(mcmc_res_ls = hybrid_res, dat = dat_ls, var_df = var_df_full, sign_flip=FALSE)
    
    MALA_subset_ls <- subset_samps(samps_ls = MALA_samps_ls, n_samps = dat_ls$n_samps, n_burnin = dat_ls$n_burnin)
    MALA_pos_mean <- posterior_mean(MALA_subset_ls)
    # Map the likelihood variable so DIC() can find it! ---
    MALA_pos_mean$dat_llk <- MALA_pos_mean$Y_C_llk
    
    # 3. Compute DIC
    Y_flat_train <- matrix(NA, n_roi*n_roi, (n-1))
    for(k in 1:(n-1)) Y_flat_train[,k] <- c(FC_arr_train[,,k])
    
    dic_val <- DIC(mean_ls=MALA_pos_mean, Y_flat=t(Y_flat_train), C=as.matrix(clinical_df_train[,"Cognitive_Status"]), Z=as.matrix(cbind(1,cov_df_train)), M=(n-1), N=n_roi, q=q, r=r_val)
    dic_log[i, paste0("q=", q)] <- dic_val
    
    # 4. Parsimonious Selection
    if (!is.na(dic_val) && (dic_val < best_dic)) {
      best_dic <- dic_val
      best_q <- q
      best_model_params <- list(
        U_train = matrix(MALA_pos_mean$U, n_roi, q),
        s2_train = MALA_pos_mean$s2,
        t2_lambda_train = MALA_pos_mean$t2_lambda,
        Lambda_train = matrix(MALA_pos_mean$Lambda_flat, (n-1), q)
      )
    }
  }
  
  selected_qs[i] <- best_q
  cat(sprintf("Selected q = %d (DIC = %.2f)\n", best_q, best_dic))
  
  # B. Final Prediction for Held-Out Sample
  lambda_test_vec <- project_lambda_test(
    Y_test = FC_test, 
    U_train = best_model_params$U_train, 
    s2_train = best_model_params$s2_train, 
    t2_lambda_train = best_model_params$t2_lambda_train
  )
  
  # Refit GLM with selected q
  train_df_glm <- data.frame(
    Cognitive_Status = clinical_df_train$Cognitive_Status,
    Atherosclerosis = clinical_df_train$Atherosclerosis
  )
  for(k in 1:best_q) train_df_glm[[paste0("lambda_", k)]] <- best_model_params$Lambda_train[, k]
  
  test_df_glm <- data.frame(Atherosclerosis = clinical_test$Atherosclerosis)
  for(k in 1:best_q) test_df_glm[[paste0("lambda_", k)]] <- lambda_test_vec[k]
  
  # Min-Max Scaling to prevent numerical instability in glm()
  train_mins <- numeric(best_q)
  train_maxs <- numeric(best_q)
  for(k in 1:best_q) {
    col_name <- paste0("lambda_", k)
    train_mins[k] <- min(train_df_glm[[col_name]])
    train_maxs[k] <- max(train_df_glm[[col_name]])
    
    # Scale train
    train_df_glm[[col_name]] <- (train_df_glm[[col_name]] - train_mins[k]) / (train_maxs[k] - train_mins[k])
    # Scale test strictly using training boundaries
    test_df_glm[[col_name]] <- (test_df_glm[[col_name]] - train_mins[k]) / (train_maxs[k] - train_mins[k])
  }
  
  form <- as.formula(paste0("Cognitive_Status ~ Atherosclerosis + ", paste0("lambda_", 1:best_q, collapse = " + ")))
  suppressWarnings(final_fit <- glm(form, data = train_df_glm, family = "binomial"))
  
  # Predict
  suppressWarnings(test_prob <- predict(final_fit, newdata = test_df_glm, type = "response"))
  
  test_predicted_probs[i] <- test_prob
  test_predicted_labels[i] <- ifelse(test_prob > 0.5, 1, 0)
  test_true_labels[i] <- clinical_test$Cognitive_Status
  
  gc()
}

# ==============================================================================
# 4. SAVE INTERMEDIATE RESULTS FOR PART 2
# ==============================================================================
report_dir <- file.path("/home/koe/BSNMani_application-main/Github_BSNMani/result", cell_type, "Unbiased_LOOCV_Report")
if(!dir.exists(report_dir)) dir.create(report_dir, recursive = TRUE)

save(test_predicted_probs, test_predicted_labels, test_true_labels, selected_qs, dic_log, n, cell_type, 
     file = file.path(report_dir, "LOOCV_Intermediate_Predictions.RData"))

cat(sprintf("\n[Part 1 Complete] Raw predictions successfully saved to:\n   %s\n", file.path(report_dir, "LOOCV_Intermediate_Predictions.RData")))
