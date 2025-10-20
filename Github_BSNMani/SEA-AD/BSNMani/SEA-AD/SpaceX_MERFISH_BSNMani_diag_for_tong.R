#########################################
##### diagnostics: model fitting
##### GL environment
##### module load Bioinformatics
##### module load Rgiotto
#########################################

######## set environment
rm(list=ls())
library(data.table)
library(tidyr)
library(readxl)
library(Rcpp)
library(RcppArmadillo)
library(abind)
library(coda)
library(gplots)
library(RColorBrewer)
library(dplyr)
library(brainconn)
library(gridExtra)
library(hrbrthemes)
library(readxl)

setwd("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SEA-AD/BSNMani-dev") ##"/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes"

#### load functions and set parameters
source("g1_diagnostics_helper.R")
sourceCpp("hybrid_M0_MALA_LR_FAST_v2.cpp")
sourceCpp("hybrid_M0_MALA_LR_g2.cpp")

fix_ls = c("")
SGI_g1 = "GQN"
SGI_g2 = "KPN"
dat_setting = "all"

########################################################################
#For array calculation
########################################################################
args <- commandArgs(trailingOnly = TRUE)
q_val <- as.numeric(args[1])  
cat("Running with q_val =", q_val, "\n")
#######################################################################

## check point 1/5                                      ####
cell_type = "Oli"
##

cat("q_val is:",q_val)
cat("celltype is:",cell_type)
r_val = 2
n_chains=2

#########################################
### posterior summary and diagnostics
#########################################
## check point 2/5                                      ####
clinical_df = read_excel("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SpaceX_BSNMani/data/meta.xlsx")
wgcna_transformed_ls = readRDS("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SpaceX_BSNMani/data/Smoothie_co_expression_list_Oli.RDS")
##

n_pol = length(wgcna_transformed_ls)
n_roi = nrow(wgcna_transformed_ls[[1]])

FC_arr = array(NA,dim=c(n_roi,n_roi,n_pol))
FC_flat_mask_full = matrix(NA,n_roi*n_roi,n_pol)
for(i in 1:n_pol){
  wgcna_transformed_ls[[i]] = as.matrix(wgcna_transformed_ls[[i]])
  FC_arr[,,i] = wgcna_transformed_ls[[i]]
  FC_flat_mask_full[,i] = c(wgcna_transformed_ls[[i]])
}

## modified covariates By TONG 09 MAY 2025
cov_df = as.matrix(clinical_df %>% select(Atherosclerosis))

######## load real data results (whole)
save_path = fs::path("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SpaceX_BSNMani/result",cell_type,paste("q",q_val,sep="_"))

hybrid_res = list()
fname_suffix_MH = paste("q",q_val,sep="_")
for(c in 1:n_chains){
  tmp = readRDS(file = fs::path(save_path,"MH",paste("MH_res",SGI_g1,SGI_g2,"chain",c,fname_suffix_MH,sep="_"),ext="RDS"))
  hybrid_res[[c]] = list(mcmc = tmp)
}

dat_ls = list(N = n_roi, q = q_val, n_chains = n_chains, n_samps = length(hybrid_res[[1]]$mcmc$s2), n_burnin = 30000)   ## change from 30000 to 8000 by tong

###### gather samples (ALL: 1e+5)
print("formatting samples")
MALA_samps_ls = gather_samples(mcmc_res_ls = hybrid_res, dat = dat_ls,
                               var_df = data.frame(var = c("l_lambda","Lambda_flat","s2","t2_lambda","X","U","Y_C_llk","Y_llk","C_llk","d","t2_alpha","t2_beta","t2"),
                                                   dim = c(1,3,1,1,3,3,1,1,1,2,1,1,1)), sign_flip=FALSE)
heidel_ls = lapply(MALA_samps_ls,FUN=function(x){heidel.diag(mcmc.list(x))})

########## convergence: Rhat, Geweke, neff
print("obtaining convergence diagnostics")
MALA_subset_ls = subset_samps(samps_ls = MALA_samps_ls[c("Lambda_flat","s2","t2_lambda","X","U","Y_C_llk","Y_llk","C_llk","d","t2_alpha","t2_beta","t2")], n_samps = dat_ls$n_samps, n_burnin = dat_ls$n_burnin)
saveRDS(MALA_subset_ls,file="/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SpaceX_BSNMani/MALA_subset_ls_debug.RDS")
rm(MALA_samps_ls)

MALA_Rhat = Rhat(samps_subset = MALA_subset_ls[c("Lambda_flat","s2","t2_lambda","X","U","Y_C_llk","Y_llk","C_llk","d","t2_alpha","t2_beta","t2")])

########## model selection DIC
print("computing posterior mean and DIC")
MALA_pos_mean = posterior_mean(MALA_subset_ls[c("Lambda_flat","s2","t2_lambda","X","U","Y_C_llk","Y_llk","C_llk","d","t2_alpha","t2_beta","t2")])
MALA_pos_mean$dat_llk = MALA_pos_mean$Y_C_llk
DIC = DIC(mean_ls=MALA_pos_mean, Y_flat=t(FC_flat_mask_full), C=as.matrix(clinical_df[,"Cognitive_Status"]), Z=as.matrix(cbind(1,cov_df)), M=n_pol, N=n_roi, q=q_val,r=r_val)

########## training R2
print("training R2")

## change merFISH_clinical.RDS (the clinical file that we open) to clinical_df_final

## check point 3/5                            ####
clinical_df = read_excel("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SpaceX_BSNMani/data/meta.xlsx")
##

n_pol = nrow(clinical_df)

## Add covariates By TONG 09 MAY 2025
cov_df = as.matrix(clinical_df %>% select(Atherosclerosis))
Lam_mean = matrix(MALA_pos_mean$Lambda_flat,n_pol,q_val)
C_pred = cbind(Lam_mean,1,cov_df) %*% MALA_pos_mean$d
train_R2 = cor(clinical_df$Cognitive_Status,C_pred)^2
train_RMSE = sqrt(mean((clinical_df$Cognitive_Status-C_pred)^2)/mean(clinical_df$Cognitive_Status^2))
train_rMSE = mean((clinical_df$Cognitive_Status-C_pred)^2)/var(clinical_df$Cognitive_Status)

########## 90% CI for coefficients
d_CI = apply(Reduce(rbind,MALA_subset_ls$d),2,FUN=function(x){quantile(x,c(0.05,0.95))})
d_CI_80 = apply(Reduce(rbind,MALA_subset_ls$d),2,FUN=function(x){quantile(x,c(0.1,0.9))})

##########  generate PPD
print("generate PPD")
mcmc_d = MALA_subset_ls$d
mcmc_t2 = MALA_subset_ls$t2
mcmc_lambda = MALA_subset_ls$Lambda_flat
C_PPD_ls = d_PPD_ls = list()
n_samps = dat_ls$n_samps - dat_ls$n_burnin #### should match the dimensions for mcmc lambda samples [02/27/2025]
set.seed(456)

########## save diagnostics results
dir.create(fs::path(save_path,"MH","diagnostics"), recursive=TRUE, showWarnings = FALSE)

diag = list(posterior_mean = MALA_pos_mean)
saveRDS(diag, file = fs::path(save_path,"MH","diagnostics",paste("MH_diag_res",SGI_g1,SGI_g2,"q",q_val,sep="_"), ext="RDS"))
#########################################
### visualization
#########################################
dir.create(fs::path(save_path,"MH","viz"))
col.heat = colorRampPalette(c("blue", "white", "red"))(256)

## subnetworks (prior to clustering genes into functional modules)
X_est = matrix(diag$posterior_mean$X, n_roi, q_val)
U_est = polar_expansion(X_est)
U_est_file = paste0("U_est_",q_val)
saveRDS(U_est, file = fs::path(save_path,"MH","diagnostics",U_est_file,ext="RDS"))

for(i in 1:q_val){
   print(i)
    U_sub = U_est[,i,drop=FALSE]%*%t(U_est[,i,drop=FALSE])
    diag(U_sub)=NA

    png(fs::path(save_path,"MH","viz",paste("MH_diag",SGI_g1,SGI_g2,"U_subnet_q",q_val,"i",i,sep="_"),ext="png"),width=12,height=12,units="in",res=300)
    heatmap.2(U_sub, col=col.heat, scale="none",dendrogram='none',main=paste("q =",i),lhei=c(1,5),breaks = NULL, na.color="grey", lwid=c(1,5), trace='none', cexRow=1/sqrt(n_roi)*2, cexCol=1/sqrt(n_roi)*2)
    dev.off()
}

rm(list = ls())
print("#######grouping genes###################")

## check point 4/5                            #### genes' name vector may need to be changed when changing the cell type
args <- commandArgs(trailingOnly = TRUE)
q_val <- as.numeric(args[1])
cell_type = "Oli"
##
source("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SpaceX_BSNMani/building/MerFISH_methods_of_grouping_genes_TL/MerFISH_methods_of_grouping_genes_TL.R")


rm(list = ls())
print("#######clinical model construction & LOOCV###################")

## check point 5/5                           ####
args <- commandArgs(trailingOnly = TRUE)
q_val <- as.numeric(args[1]) 
cell_type = "Oli"
print("no lambda filtered\n")
##
                                               ####
## need to change code inside this R script when changing the cell type   ####
source("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SpaceX_BSNMani/building/clincal_model_construction/clinical_model_construction.R")



