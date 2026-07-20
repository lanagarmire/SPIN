#' Run BSNMani on BC data
#' Modified from /nfs/dcmb-lgarmire/haowenwu/Breast_Cancer_IMC_BSNMani.R

source("/nfs/dcmb-lgarmire/yangiwen/workspace/common/Utils.R")

# Get running options
spec <- matrix(
  c(
    "seed", "s", 1, "integer",
    "qval", "q", 1, "integer",
    "output", "o", 1, "character"
  ),
  byrow = T,
  ncol = 4
)

opt <- getopt::getopt(spec)
seed <- opt[["seed"]]
q_val <- opt[["qval"]]
output_dir <- opt[["output"]]
fold <- 5

set.seed(seed)
mkdir(output_dir)

#########################################
#### modeling fitting for SEA-AD
#### GL environment
#### module load Bioinformatics
#### module load Rgiotto
#########################################

#########################################
##### set parameters
#########################################
library(dplyr)
library(tidyr)
library(abind)
library(Rcpp)
library(RcppArmadillo)
library(SparseGrid)
library(matrixStats)
library(data.table)
library(pracma)
library(coda)
library(survival)

## pull the scripts from github
sourceCpp("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SEA-AD/BSNMani-dev/hybrid_M0_MALA_LR_FAST_v2.cpp")
sourceCpp("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SEA-AD/BSNMani-dev/hybrid_M0_MALA_LR_A_lambda_new.cpp")
sourceCpp("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SEA-AD/BSNMani-dev/hybrid_M0_MALA_LR_g2.cpp")
source("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SEA-AD/BSNMani-dev/g1_BFGS_init.R")
source("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SEA-AD/BSNMani-dev/g1_diagnostics_helper.R")
source("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SEA-AD/BSNMani-dev/two_stage_train_pipeline.R")

fix_ls <- ""
SGI_g1 <- "GQN"
SGI_g2 <- "KPN"
dat_setting <- "all"
r_val <- 4

#########################################
### load data
#########################################
folds <- readRDS("/nfs/dcmb-lgarmire/haowenwu/5CV_folds.RDS")
clinical_df <- readRDS("/nfs/dcmb-lgarmire/haowenwu/clinical_df_3cov.RDS")
# clinical_df <- clinical_df[clinical_df$clinical_type == "HR-HER2-",]
# clinical_df$clinical_type_int <- as.integer(factor(clinical_df$clinical_type, levels = c("HR+HER2-", "HR+HER2+", "HR-HER2+", "HR-HER2-"))) - 1
clinical_df$clinical_type <- factor(clinical_df$clinical_type, levels = c("HR+HER2-", "HR+HER2+", "HR-HER2+", "HR-HER2-"))
clinical_df$grade <- factor(clinical_df$grade, levels = c(1, 2, 3))
clinical_df_train <- clinical_df[clinical_df$core %in% folds[[fold]],]
clinical_df_test <- clinical_df[!clinical_df$core %in% folds[[fold]],]

wgcna_transformed_ls <- readRDS("/nfs/dcmb-lgarmire/haowenwu/WGCNA_transformed_BC_core.RDS")
wgcna_transformed_ls_train <- wgcna_transformed_ls[clinical_df_train$core]
wgcna_transformed_ls_test <- wgcna_transformed_ls[clinical_df_test$core]

n_pol_train <- length(wgcna_transformed_ls_train)
n_pol_test <- length(wgcna_transformed_ls_test)
n_roi <- nrow(wgcna_transformed_ls_train[[1]])

## network input: n_roi x n_roi x n_pol array
FC_arr_train <- array(NA, dim = c(n_roi, n_roi, n_pol_train))
for (i in 1:n_pol_train) {
  FC_arr_train[, , i] <- wgcna_transformed_ls_train[[i]]
}

FC_arr_test <- array(NA, dim = c(n_roi, n_roi, n_pol_test))
for (i in 1:n_pol_test) {
  FC_arr_test[, , i] <- wgcna_transformed_ls_test[[i]]
}

cov_df_train <- as.matrix(clinical_df_train[c("age", "grade", "clinical_type")])
cov_df_test <- as.matrix(clinical_df_test[c("age", "grade", "clinical_type")])

### load g1 results
n_burnin <- 2e+3
n_samps <- 5e+3
fix_ls_train <- ""
fix_ls_test <- c("U", "X", "s2")

#########################################
### run BSNMani stage 1
#########################################
n_chains <- 2
save_path_train <- file.path(output_dir, "train")
mkdir(save_path_train)
for (chain_idx in 1:n_chains) {
  two_stage_single_chain_train(
    n_roi = n_roi, q_val = q_val, q_val_0 = q_val, n_pol = n_pol_train,
    FC_dat = FC_arr_train, idx_ls = list(module_1 = 1:n_roi / 2, module_2 = (1 + n_roi / 2):n_roi), mask = TRUE,
    clinical_df = as.matrix(clinical_df_train[, "Patientstatus"]), cov_df = cov_df_train, r_val = r_val,    ## data
    noise = 0.001, t2_lambda_0 = 1000, ## g1-initialize
    nu0 = 2, s20 = 1, eta0 = 2, t20 = 1, SGI_g1 = "GQN", k1 = 25, stepsize = 1e-3, acpt_step = 50, target_acpt = 0.4, tune = TRUE,
    fix_ls = fix_ls_train, ## g1
    rho0 = 2, psi20 = 1, gamma0 = 2, kappa20 = 1, omega0 = 2, phi20 = 1, ## g2
    k2 = 10, SGI_g2 = "KPN", g2_weighted = TRUE, ## MH
    seed = seed + chain_idx, burn_in = n_burnin, mcmc_sample = n_samps, chain_idx = chain_idx,
    save_path_g1 = save_path_train,
    save_path_g2 = save_path_train,
    save_path_MH = save_path_train,
    fname_suffix_g1 = "train",
    fname_suffix_g2 = "train",
    fname_suffix_MH = "train",
    run_g1 = TRUE, run_g2 = FALSE, run_MH = FALSE,
    load_g1 = FALSE, load_g2 = FALSE, load_MH = FALSE,
    save_g1 = TRUE, save_g2 = FALSE, save_MH = FALSE,
    dir_tags = "train",
    basename_tags = "train"
  )
}

############################
### diagnostics on stage 1
############################
g1_train_res <- list()
fname_suffix_g1_train <- "train"
for (c in 1:n_chains) {
  tmp <- readRDS(file = fs::path(save_path_train, "g1", paste("g1_res", SGI_g1, "chain", c, fname_suffix_g1_train, sep = "_"), ext = "RDS"))
  g1_train_res[[c]] <- tmp
}
str(g1_train_res)

dat_ls <- list(N = n_roi, q = q_val, n_chains = n_chains, n_samps = length(g1_train_res[[1]]$mcmc$s2), n_burnin = 2e+3)
str(dat_ls)

###### gather samples (ALL: 1e+5)
print("formatting samples")
MALA_samps_ls_train <- gather_samples(mcmc_res_ls = g1_train_res, dat = dat_ls,
                                      var_df = data.frame(var = c("Lambda_flat", "s2", "t2_lambda", "X", "U"),
                                                          dim = c(3, 1, 1, 3, 3)), sign_flip = FALSE)
heidel_ls <- lapply(MALA_samps_ls_train, FUN = function(x) { heidel.diag(mcmc.list(x)) })

########## posterior mean
print("computing posterior mean")
MALA_subset_ls_train <- subset_samps(samps_ls = MALA_samps_ls_train[c("Lambda_flat", "s2", "t2_lambda", "X", "U")], n_samps = dat_ls$n_samps, n_burnin = dat_ls$n_burnin)
MALA_pos_mean_train <- posterior_mean(MALA_subset_ls_train[c("Lambda_flat", "s2", "t2_lambda", "X", "U")])

### save diagnostics results
if (!dir.exists(fs::path(save_path_train, "g1", "diagnostics"))) {
  dir.create(fs::path(save_path_train, "g1", "diagnostics"), recursive = TRUE)
}
saveRDS(MALA_pos_mean_train, file = fs::path(save_path_train, "g1", "diagnostics", paste("g1_res", SGI_g1, fname_suffix_g1_train, sep = "_"), ext = "RDS"))

############################
### run frequentist survival analysis (change this stratified survival analysis)
############################
Lambda_mean_train <- matrix(MALA_pos_mean_train$Lambda_flat, n_pol_train, q_val)

lambda_vars <- paste0("lambda_", 1:q_val)
colnames(Lambda_mean_train) <- lambda_vars
clinical_df_train <- cbind(clinical_df_train, Lambda_mean_train)

# vars <- c("age", "grade", "clinical_type_int", lambda_vars)
clinical_df_train[, lambda_vars] <- lapply(clinical_df_train[, lambda_vars], minmax)
formula <- as.formula(paste("Surv(OSmonth, Patientstatus) ~ grade + age + clinical_type +", paste(lambda_vars, collapse = " + ")))
mod <- survival::coxph(formula, data = clinical_df_train)

if (!dir.exists(fs::path(save_path_train, "g2_survival"))) {
  dir.create(fs::path(save_path_train, "g2_survival"), recursive = TRUE)
}
saveRDS(mod, file = fs::path(save_path_train, "g2_survival", paste("g2_survival_res", fname_suffix_g1_train, sep = "_"), ext = "RDS"))

###########################
### make prediction
###########################
X_mean <- matrix(MALA_pos_mean_train$X, n_roi, q_val)
U_mean <- polar_expansion(X_mean)
s2_mean <- c(MALA_pos_mean_train$s2)
t2_lambda_mean <- c(MALA_pos_mean_train$t2_lambda)

init_ls <- list()
sgrid_g1 <- createSparseGrid(type = SGI_g1, dimension = 1, k = 25)
region_idx <- sort(sapply(list(module_1 = 1:n_roi / 2, module_2 = (1 + n_roi / 2):n_roi), FUN = function(x) { sample(x, 1) }))
n_region_pairs <- length(region_idx) * (length(region_idx) - 1) / 2
save_path_test <- file.path(output_dir, "test")
mkdir(save_path_test)
if (!dir.exists(fs::path(save_path_test, "g1"))) {
  dir.create(fs::path(save_path_test, "g1"), recursive = TRUE)
}
for (chain_idx in 1:n_chains) {
  set.seed(seed + chain_idx)

  init_ls$X <- X_mean
  init_ls$Lambda <- matrix(rnorm(n_pol_test * q_val, 0, 0.001), n_pol_test, q_val) + Lambda_cls(U_mean, FC_arr_test, n_pol_test, q_val)
  init_ls$s2 <- s2_mean
  init_ls$t2_lambda <- t2_lambda_mean

  start <- Sys.time()
  g1_res_test <- hybrid_MALA_g1(X_0 = init_ls$X, Lambda_0_flat = init_ls$Lambda, s2_0 = init_ls$s2, t2_lambda_0 = init_ls$t2_lambda,
                                nu0 = 2, s20 = 1, eta0 = 2, t20 = 1,
                                Y = FC_arr_test, M = n_pol_test, N = n_roi, q = q_val, SGI_nodes = sgrid_g1$nodes, SGI_wts = sgrid_g1$weights, weighted = (SGI_g1 == "GQN"), region_select = region_idx,
                                mcmc_sample = n_samps, stepsize = 1e-3, acpt_step = 50, Gibbs_step = 1, MALA_step = 1, target_acpt = 0.4, tune = TRUE, fixed = fix_ls_test)
  end <- Sys.time()
  runtime <- difftime(end, start, units = "secs")
  g1_res_test$runtime <- runtime
  g1_res_test$init <- init_ls

  ## save
  fname_suffix_g1_test <- "test"
  saveRDS(g1_res_test, file = fs::path(save_path_test, "g1", paste("g1_res", SGI_g1, "chain", chain_idx, fname_suffix_g1_test, sep = "_"), ext = "RDS"))
}

g1_test_res <- list()
fname_suffix_g1_test <- "test"
for (c in 1:n_chains) {
  tmp <- readRDS(file = fs::path(save_path_test, "g1", paste("g1_res", SGI_g1, "chain", c, fname_suffix_g1_test, sep = "_"), ext = "RDS"))
  g1_test_res[[c]] <- tmp
}
str(g1_test_res)

###### gather samples (ALL: 1e+5)
print("formatting samples")
MALA_samps_ls_test <- gather_samples(mcmc_res_ls = g1_test_res, dat = dat_ls,
                                     var_df = data.frame(var = c("Lambda_flat", "t2_lambda"),
                                                         dim = c(3, 1)), sign_flip = FALSE)
########## posterior mean
print("computing posterior mean")
MALA_subset_ls_test <- subset_samps(samps_ls = MALA_samps_ls_test[c("Lambda_flat", "t2_lambda")], n_samps = dat_ls$n_samps, n_burnin = dat_ls$n_burnin)
MALA_pos_mean_test <- posterior_mean(MALA_subset_ls_test[c("Lambda_flat", "t2_lambda")])

##### make predictions
Lambda_mean_test <- matrix(MALA_pos_mean_test$Lambda_flat, n_pol_test, q_val)
colnames(Lambda_mean_test) <- lambda_vars
clinical_df_test <- cbind(clinical_df_test, Lambda_mean_test)

clinical_df_test[, lambda_vars] <- lapply(clinical_df_test[, lambda_vars], minmax)
mod_pred <- predict(mod, clinical_df_test, type = "lp")
pred_res <- list(mod_pred = mod_pred, Lambda_mean = Lambda_mean_test)

if (!dir.exists(fs::path(save_path_test, "g2_survival"))) {
  dir.create(fs::path(save_path_test, "g2_survival"), recursive = TRUE)
}
saveRDS(pred_res, fs::path(save_path_test, "g2_survival", paste("g2_survival", fname_suffix_g1_test, sep = "_"), ext = "RDS"))

#save test information
saveRDS(mod, file = fs::path(save_path_train, "g2_survival", paste("g2_survival_res", fname_suffix_g1_train, sep = "_"), ext = "RDS"))

lp_train <- predict(mod, type = "lp")

saveRDS(lp_train, fs::path(save_path_train, "g2_survival", paste0("lp_train_", fname_suffix_g1_train, ".RDS")))
saveRDS(clinical_df_train, fs::path(save_path_train, "g2_survival", paste0("clinical_df_train_", fname_suffix_g1_train, ".RDS")))

saveRDS(mod_pred, fs::path(save_path_test, "g2_survival", paste0("lp_test_", fname_suffix_g1_test, ".RDS")))
saveRDS(clinical_df_test, fs::path(save_path_test, "g2_survival", paste0("clinical_df_test_", fname_suffix_g1_test, ".RDS")))

