#!/usr/bin/env Rscript

.libPaths(c("~/R_libs_amd", .libPaths()))
Sys.setenv(PKG_CXXFLAGS="-O3 -march=native -mtune=native")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(Rcpp)
  library(RcppArmadillo)
  library(SparseGrid)
  library(abind)
  library(matrixStats)
  library(data.table) 
  library(pracma)
  library(coda)
  library(readxl)
  library(fs)
})

# Parse command line arguments from the Slurm script
args <- commandArgs(trailingOnly = TRUE)
fold_idx <- as.numeric(args[1])
q_val <- as.numeric(args[2])
# chain_idx <- 1
cat(sprintf("Running Fold %d with q = %d\n", fold_idx, q_val))

# Configuration
# cell_type = "WGCNA"
cell_type = "Smoothie"
# cell_type = "Oli"
# cell_type = "hdWGCNA"
# cell_type = "SpaceX"
r_val <- 2

# Load sources
source("/home/koe/BSNMani_application-main/Github_BSNMani/IMC/Preprocessed/Utils.R")
base_src_dir <- "/home/koe/BSNMani_application-main/Github_BSNMani/SEA-AD/BSNMani/SEA-AD/BSNMani-dev"
sourceCpp(file.path(base_src_dir, "hybrid_M0_MALA_LR_FAST_v2.cpp"), rebuild = FALSE)
sourceCpp(file.path(base_src_dir, "hybrid_M0_MALA_LR_A_lambda_new.cpp"), rebuild = FALSE)
sourceCpp(file.path(base_src_dir, "hybrid_M0_MALA_LR_g2.cpp"), rebuild = FALSE)
source(file.path(base_src_dir, "g1_BFGS_init.R"))
source(file.path(base_src_dir, "g1_diagnostics_helper.R"))
source(file.path(base_src_dir, "two_stage_train_pipeline.R"))

# Load Full Data
## check point                     ####
setwd("/home/koe/BSNMani_application-main/Github_BSNMani/SEA-AD")
clinical_df = readRDS("SEA-AD_Data/filtered_baseline_meta.rds")
# clinical_df = read_excel("SEA-AD_Data/meta.xlsx")

# wgcna_transformed_ls = readRDS("SEA-AD_Data/Filtered_baseline_co_expression_list.rds")
wgcna_transformed_ls = readRDS("SEA-AD_Data/Filtered_Smoothie_co_expression_list.rds")
# wgcna_transformed_ls = readRDS("SEA-AD_Data/Smoothie_co_expression_list_Oli.RDS")
# wgcna_transformed_ls = readRDS("SEA-AD_Data/hdWGCNA_co_expression_list.RDS")
# wgcna_transformed_ls = readRDS("SEA-AD_Data/Shared_co_expression_list.RDS")

n <- nrow(clinical_df)
n_roi <- nrow(wgcna_transformed_ls[[1]])

# Added diagonal zeroing loop to fix network topology mapping
for(i in 1:length(wgcna_transformed_ls)) {
  tmp_mat <- as.matrix(wgcna_transformed_ls[[i]])
  diag(tmp_mat) <- 0
  wgcna_transformed_ls[[i]] <- tmp_mat
}

FC_arr <- array(NA, dim=c(n_roi, n_roi, n))
for(i in 1:n){
  FC_arr[,,i] <- as.matrix(wgcna_transformed_ls[[i]])
}

# Ensure Atherosclerosis is strictly numeric
clinical_df$Atherosclerosis <- as.numeric(as.character(clinical_df$Atherosclerosis))
cov_df <- as.matrix(clinical_df %>% select(Atherosclerosis))

# Split Data: Hold out 'fold_idx'
train_idx <- setdiff(1:n, fold_idx)
FC_arr_train <- FC_arr[,,train_idx]
clinical_df_train <- clinical_df[train_idx, ]
cov_df_train <- cov_df[train_idx, , drop=FALSE]

# Define save paths isolated by fold
fold_dir <- fs::path("/home/koe/BSNMani_application-main/Github_BSNMani/result", cell_type, paste0("Fold_", sprintf("%02d", fold_idx)))

for (chain_idx in 1:2) {
  suppressWarnings({
    # Run the training pipeline strictly on n-1 samples
    two_stage_single_chain_train(
      n_roi = n_roi, q_val = q_val, q_val_0 = q_val, n_pol = length(train_idx), 
      FC_dat = FC_arr_train, 
      idx_ls = list(module_1 = 1:n_roi/2, module_2=(1+n_roi/2):n_roi), 
      mask = TRUE, 
      clinical_df = as.matrix(clinical_df_train[,"Cognitive_Status"]), 
      cov_df = cov_df_train, r_val = r_val,    
      noise = 0.001, t2_lambda_0 = 1000, 
      nu0 = 2, s20 = 1, eta0=2, t20=1, SGI_g1 = "GQN", k1 = 25, stepsize=1e-3, 
      acpt_step = 50, target_acpt = 0.25, tune = TRUE, fix_ls = c(""), 
      rho0=2, psi20=1, gamma0=2, kappa20=1, omega0=2, phi20=1,
      k2 = 10, SGI_g2 = "KPN", g2_weighted = TRUE, 
      seed = 12345 + fold_idx + chain_idx, burn_in = 100000, mcmc_sample = 200000, chain_idx = chain_idx, 
      save_path_g1 = fs::path(fold_dir, paste("q", q_val, sep="_")),
      save_path_g2 = fs::path(fold_dir, paste("q", q_val, sep="_")),
      save_path_MH = fs::path(fold_dir, paste("q", q_val, sep="_")),
      fname_suffix_g1 = paste("q", q_val, sep="_"), 
      fname_suffix_g2 = paste("q", q_val, sep="_"),  
      fname_suffix_MH = paste("q", q_val, sep="_"),  
      run_g1=TRUE, run_g2=TRUE, run_MH=TRUE,
      load_g1=FALSE, load_g2=FALSE, load_MH=FALSE,
      save_g1=TRUE, save_g2=TRUE, save_MH=TRUE,
      dir_tags = c("result","q"),
      basename_tags = c("q","chain")
    )
  })
}

cat(sprintf("Successfully completed Fold %d with q = %d\n", fold_idx, q_val))