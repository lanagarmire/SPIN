#########################################
#### modeling fitting for SEA-AD
#### GL environment
#### module load Bioinformatics
#### module load Rgiotto
#########################################

#########################################
##### set parameters
#########################################
rm(list=ls())
library(dplyr)
library(tidyr)
library(readxl)
library(abind)
library(Rcpp)
library(RcppArmadillo)
library(SparseGrid)
library(matrixStats)
library(data.table) 
library(pracma)
library(coda)
library(readxl)

setwd("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SEA-AD/BSNMani-dev") ##"/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes"
## pull the scripts from github
sourceCpp("hybrid_M0_MALA_LR_FAST_v2.cpp")
sourceCpp("hybrid_M0_MALA_LR_A_lambda_new.cpp")
sourceCpp("hybrid_M0_MALA_LR_g2.cpp")
source("g1_BFGS_init.R")
source("g1_diagnostics_helper.R")
source("two_stage_train_pipeline.R")

########################################################################
#For array calculation
########################################################################
args <- commandArgs(trailingOnly = TRUE)
q_val <- as.numeric(args[1])
chain_idx <- as.numeric(args[2])
cat("Running for q_val =", q_val, "chain_idx =", chain_idx, "\n")
#######################################################################

## check point 1/2                          ####
#q_val = 5   ## Change the number from 3 to 8 by Tong
##


fix_ls = c("")
SGI_g1 = "GQN"
SGI_g2 = "KPN"
dat_setting = "all"
r_val = 2   ## change by TONG, since we have Arteriolosclerosis, RIN, PMI

## check point                   ####
cell_type = "Oli"
#########################################
### load data
#########################################

## check point                     ####
clinical_df = read_excel("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SpaceX_BSNMani/data/meta.xlsx")
wgcna_transformed_ls = readRDS("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SpaceX_BSNMani/data/Smoothie_co_expression_list_Oli.RDS")
##

n_pol = length(wgcna_transformed_ls)
n_roi = nrow(wgcna_transformed_ls[[1]])

## network input: n_roi x n_roi x n_pol array
FC_arr = array(NA,dim=c(n_roi,n_roi,n_pol))
for(i in 1:n_pol){
  wgcna_transformed_ls[[i]] = as.matrix(wgcna_transformed_ls[[i]])
  FC_arr[,,i] = wgcna_transformed_ls[[i]]
}

## Add other covariates By TONG 16 MAY 2025
cov_df = as.matrix(clinical_df %>% select(Atherosclerosis))

### load g1 results
n_burnin = 9e+4   ## change from 9e+4 to 50
n_samps = 1e+5    ## change from 1e+5 to 1000 for faster computaion
fix_ls = c("")

#########################################
### run training two-stage model
#########################################
## n_roi: number of ROIs
## q_val: number of latent dimensions (explore multiple values for this parameter)
## q_val_0: the real number of latent dimensions if available
## n_pol: number of patients
## FC_dat: n_roi x n_roi x n_pol array (network input)
## idx_ls: named list of ROI indices whose subnetwork edges to track 
## mask: whether to mask the diagonal entries of the input networks (default: TRUE)
## clinical_df: clinical outcome dataframe
## cov_df: clinical covariate dataframe
## r_val: number of clinical covariates
## noise: sampling noise parameter for initialization (default: 0.001)
## t2_lambda_0: initialiation parameter (default: 1000)
## nu0, s20: hyperparameters for the s2 prior (default: (2,1))
## eta0, t20: hyperparameters for the t2_lambda prior (default: (2,1))
## SGI_g1: sparse grid integration scheme for stage-one sampling (default: GQN)
## k1: sparse grid integration parameter for stage-one sampling (default: 25)
## stepsize: step size parameter for MALA sampling in stage 1 (default: 1e-3)
## acpt_step: number of steps before adjusting acceptance rate in MALA sampling in stage 1 (default: 50)
## target_acpt: target acceptance rate for MALA sampling in stage 1 (the higher the number of ROIs the lower the target acceptance rate. for lower dimension networks with around 100 ROIs: around 0.4; for higher dimension networks with around 400 ROIs: around 0.15)
## tune: whether to tune the stepsize (default: TRUE)
## fix_ls: which parameters to set as fixed; keep empty unless debugging (default: c(""))
## rho0, psi20: hyperparameters for the t2 prior (default: (2,1))
## gamma0, kappa20: hyperparameters for the t2_beta prior (default: (2,1))
## omega0, phi20: hyperparameters for the t2_alpha prior (default: (2,1))
## k2: sparse grid integration parameter for stage-two sampling (default: 10)
## SGI_g2: sparse grid integration scheme for stage-two sampling (default: KPN)
## g2_weighted: sparse grid integration parameter for stage-two sampling (default: TRUE)
## seed: random seed for MCMC sampling
## burn_in: number of burn-in samples
## mcmc_sample: number of MCMC samples
## chain_idx: MCMC chain index
## save_path_g1: path for saving stage-one sampling results
## save_path_g2: path for saving stage-two sampling results
## save_path_MH: path for saving Metropolis Hastings sampling results
## fname_suffix_g1: stage-one sampling results filename suffix
## fname_suffix_g2: stage-two sampling results filename suffix
## fname_suffix_MH: Metropolis Hastings sampling results filename suffix
## run_g1: whether to run stage-one sampling (default: TRUE)
## run_g2: whether to run stage-two sampling (default: TRUE)
## run_MH: whether to run Metropolis Hastings sampling (default: TRUE)
## load_g1: whether to load existing stage-one sampling results (default: FALSE)
## load_g2: whether to load existing stage-two sampling results (default: FALSE)
## load_MH: whether to load Metropolis Hastings sampling results (default: FALSE)
## save_g1: whether to save stage-one sampling results (default: TRUE)
## save_g2: whether to save stage-two sampling results (default: TRUE)
## save_MH: whether to save Metropolis Hastings sampling results (default: TRUE)
## dir_tags: keywords for sanity checking the results directory path
## basename_tags: keywords for sanity checking the results file names

## change "clinical_df = as.matrix(clinical_df[,"MMSE"])" TO "Cognitive_Status"
two_stage_single_chain_train(n_roi = n_roi, q_val = q_val, q_val_0 = q_val, n_pol = n_pol, FC_dat = FC_arr, idx_ls = list(module_1 = 1:n_roi/2,module_2=(1+n_roi/2):n_roi), mask = TRUE, clinical_df = as.matrix(clinical_df[,"Cognitive_Status"]), cov_df = cov_df, r_val = r_val,    ## data
                             noise = 0.001, t2_lambda_0 = 1000, ## g1-initialize
                             nu0 = 2, s20 = 1, eta0=2, t20=1, SGI_g1 = "GQN", k1 = 25, stepsize=1e-3, acpt_step = 50, target_acpt = 0.25, tune = TRUE, fix_ls = fix_ls, ## g1
                             rho0=2, psi20=1, gamma0=2, kappa20=1, omega0=2, phi20=1,## g2
                             k2 = 10, SGI_g2 = "KPN", g2_weighted = TRUE, ## MH
                             seed = 12345+chain_idx, burn_in = n_burnin, mcmc_sample = n_samps, chain_idx = chain_idx, 
                             save_path_g1 = fs::path("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SpaceX_BSNMani/result",cell_type,paste("q",q_val,sep="_")),
                             save_path_g2 = fs::path("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SpaceX_BSNMani/result",cell_type,paste("q",q_val,sep="_")),
                             save_path_MH = fs::path("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SpaceX_BSNMani/result",cell_type,paste("q",q_val,sep="_")),
                             fname_suffix_g1 = paste("q",q_val,sep="_"), 
                             fname_suffix_g2 = paste("q",q_val,sep="_"),  
                             fname_suffix_MH = paste("q",q_val,sep="_"),  
                             run_g1=TRUE, run_g2=TRUE, run_MH=TRUE,
                             load_g1=FALSE,load_g2=FALSE,load_MH=FALSE,
                             save_g1=TRUE, save_g2=TRUE, save_MH=TRUE,
                             dir_tags = c("SpaceX_BSNMani","q"),
                             basename_tags = c("q","chain"))
