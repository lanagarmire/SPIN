#########################################
#### BSNMani model fitting for SEA-AD all-slide data
#### Donor-level LOOCV version
####
#### Usage:
####   Rscript BSNMani_allslides_fit_donorLOOCV.R <q_val> <chain_idx> <fold_id> [cell_type]
#########################################

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(abind)
  library(Rcpp)
  library(RcppArmadillo)
  library(data.table)
  library(coda)
  library(fs)
})

## =========================
## User-configurable paths
## =========================
BSNMANI_CODE_DIR <- "/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SEA-AD/BSNMani-dev"
DATA_DIR <- "/nfs/turbo/umms-lgarmire/liutong/BSNMani/SpaceX_BSNMani/data"
OUT_ROOT <- "/nfs/turbo/umms-lgarmire/liutong/BSNMani/SpaceX_BSNMani/BSNMani_output_result"
CLINICAL_FILE <- "/home/koe/BSNMani_application-main/Tong_results/BSNMani_for_transfer/BSNMani_for_transfer/BSNMani/SEA-AD_Updated_Data/meta_uni.xlsx"
SAMPLE_META_FILE <- "/home/koe/BSNMani_application-main/Tong_results/BSNMani_for_transfer/BSNMani_for_transfer/BSNMani/SEA-AD_Updated_Data/patient_sample_ID_all_slides.csv"

setwd(BSNMANI_CODE_DIR)

sourceCpp("hybrid_M0_MALA_LR_FAST_v2.cpp")
sourceCpp("hybrid_M0_MALA_LR_A_lambda_new.cpp")
sourceCpp("hybrid_M0_MALA_LR_g2.cpp")
source("g1_BFGS_init.R")
source("g1_diagnostics_helper.R")
source("two_stage_train_pipeline.R")

## =========================
## Arguments
## =========================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript BSNMani_allslides_fit_donorLOOCV.R <q_val> <chain_idx> <fold_id> [cell_type]")
}

q_val <- as.numeric(args[1])
chain_idx <- as.numeric(args[2])
fold_id <- as.numeric(args[3])
cell_type <- ifelse(length(args) >= 4, args[4], "SpaceX_uni")

cat("Running donor-level LOOCV fit | q =", q_val,
    "| chain =", chain_idx,
    "| fold =", fold_id,
    "| cell_type =", cell_type, "\n")

fix_ls <- c("")
SGI_g1 <- "GQN"
SGI_g2 <- "KPN"
r_val <- 2
n_burnin <- 9e4
n_samps <- 1e5

## =========================
## Helper functions
## =========================
extract_donor_id <- function(x) {
  m <- regexpr("H[0-9]+\\.[0-9]+\\.[0-9]+", x)
  ifelse(m > 0, regmatches(x, m), NA_character_)
}

safe_name <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)

standardize_clinical_cols <- function(clinical_df) {
  cn <- colnames(clinical_df)
  donor_candidates <- c("Donor ID", "donor_id", "Donor_ID", "patient", "Patient", "patient_ID")
  donor_col <- donor_candidates[donor_candidates %in% cn][1]
  if (is.na(donor_col)) stop("Cannot find donor ID column in clinical_df.")

  if (!"Cognitive_Status" %in% cn) {
    if ("Dementia" %in% cn) clinical_df$Cognitive_Status <- clinical_df$Dementia
    else if ("Dementia_v2" %in% cn) clinical_df$Cognitive_Status <- clinical_df$Dementia_v2
    else stop("Cannot find Cognitive_Status, Dementia, or Dementia_v2 in clinical_df.")
  }

  if (!"Atherosclerosis" %in% colnames(clinical_df)) stop("Cannot find Atherosclerosis column in clinical_df.")

  clinical_df$donor_id_std <- as.character(clinical_df[[donor_col]])
  clinical_df
}

load_networks <- function(cell_type) {
  file_map <- c(
    Baseline_uni = "Baseline_co_expression_list_uni_all_slides.RDS",
    smoothie_uni = "Smoothie_co_expression_list_uni_all_slides.RDS",
    SpaceX_uni   = "SpaceX_co_expression_list_uni_all_slides.RDS",
    hdWGCNA_uni  = "hdWGCNA_co_expression_list_uni_all_slides.RDS",
    Oli_uni      = "hdWGCNA_co_express_list_uni_Oli_all_slides.RDS"
  )
  if (!cell_type %in% names(file_map)) stop("Unknown cell_type: ", cell_type)
  f <- file.path(DATA_DIR, file_map[[cell_type]])
  if (!file.exists(f)) stop("Network file not found: ", f)
  readRDS(f)
}

build_sample_meta <- function(wgcna_ls) {
  sample_names <- names(wgcna_ls)
  if (is.null(sample_names)) {
    sample_names <- paste0("sample_", seq_along(wgcna_ls))
    names(wgcna_ls) <- sample_names
  }

  if (file.exists(SAMPLE_META_FILE)) {
    meta <- fread(SAMPLE_META_FILE)
    cn <- colnames(meta)
    sample_col <- c("patient_sample_ID", "sample_id", "sample_name", "sample")[c("patient_sample_ID", "sample_id", "sample_name", "sample") %in% cn][1]
    donor_col <- c("patient", "donor_id", "Donor ID", "Donor_ID")[c("patient", "donor_id", "Donor ID", "Donor_ID") %in% cn][1]

    if (!is.na(sample_col) && !is.na(donor_col)) {
      meta2 <- data.frame(
        sample_id = as.character(meta[[sample_col]]),
        donor_id = as.character(meta[[donor_col]]),
        stringsAsFactors = FALSE
      )
      meta2 <- meta2[meta2$sample_id %in% sample_names, , drop = FALSE]
      if (nrow(meta2) == length(sample_names)) {
        meta2 <- meta2[match(sample_names, meta2$sample_id), ]
        return(meta2)
      }
    }
  }

  data.frame(
    sample_id = sample_names,
    donor_id = extract_donor_id(sample_names),
    stringsAsFactors = FALSE
  )
}

align_clinical_to_slides <- function(sample_meta, clinical_df) {
  clinical_df <- standardize_clinical_cols(clinical_df)
  merged <- sample_meta %>% left_join(clinical_df, by = c("donor_id" = "donor_id_std"))
  if (any(is.na(merged$Cognitive_Status))) stop("Missing Cognitive_Status after donor-level merge.")
  if (any(is.na(merged$Atherosclerosis))) stop("Missing Atherosclerosis after donor-level merge.")
  merged
}

make_FC_array <- function(wgcna_ls) {
  n_pol <- length(wgcna_ls)
  n_roi <- nrow(wgcna_ls[[1]])
  FC_arr <- array(NA_real_, dim = c(n_roi, n_roi, n_pol))
  for (i in seq_len(n_pol)) {
    mat <- as.matrix(wgcna_ls[[i]])
    if (nrow(mat) != n_roi || ncol(mat) != n_roi) stop("Network dimensions differ at index ", i)
    FC_arr[, , i] <- mat
  }
  FC_arr
}

## =========================
## Load and align all-slide data
## =========================
clinical_df_donor <- read_excel(CLINICAL_FILE)
wgcna_transformed_ls <- load_networks(cell_type)
sample_meta <- build_sample_meta(wgcna_transformed_ls)

if (any(is.na(sample_meta$donor_id))) stop("Some donor IDs could not be parsed/found.")

clinical_slide_df <- align_clinical_to_slides(sample_meta, clinical_df_donor)
unique_donors <- sort(unique(clinical_slide_df$donor_id))

if (fold_id < 1 || fold_id > length(unique_donors)) stop("fold_id out of range.")

heldout_donor <- unique_donors[fold_id]
train_idx <- which(clinical_slide_df$donor_id != heldout_donor)
test_idx <- which(clinical_slide_df$donor_id == heldout_donor)

cat("Total slides:", length(wgcna_transformed_ls), "\n")
cat("Total donors:", length(unique_donors), "\n")
cat("Held-out donor:", heldout_donor, "\n")
cat("Training slides:", length(train_idx), "| Test slides:", length(test_idx), "\n")

fold_dir <- fs::path(
  OUT_ROOT,
  cell_type,
  "donor_LOOCV",
  paste0("fold_", sprintf("%02d", fold_id), "_", safe_name(heldout_donor)),
  paste0("q_", q_val)
)
dir.create(fold_dir, recursive = TRUE, showWarnings = FALSE)

fold_info <- list(
  q_val = q_val,
  chain_idx = chain_idx,
  fold_id = fold_id,
  heldout_donor = heldout_donor,
  train_idx = train_idx,
  test_idx = test_idx,
  sample_meta = sample_meta,
  clinical_slide_df = clinical_slide_df,
  cell_type = cell_type
)
saveRDS(fold_info, file = fs::path(fold_dir, "fold_info.rds"))

## Use training slides only for model fitting.
FC_arr_train <- make_FC_array(wgcna_transformed_ls[train_idx])
n_pol_train <- length(train_idx)
n_roi <- dim(FC_arr_train)[1]

cov_df_train <- as.matrix(clinical_slide_df[train_idx, "Atherosclerosis", drop = FALSE])
C_train <- as.matrix(clinical_slide_df[train_idx, "Cognitive_Status", drop = FALSE])

cat("Starting BSNMani model fitting...\n")

two_stage_single_chain_train(
  n_roi = n_roi,
  q_val = q_val,
  q_val_0 = q_val,
  n_pol = n_pol_train,
  FC_dat = FC_arr_train,
  idx_ls = list(module_1 = 1:(n_roi/2), module_2 = (1 + n_roi/2):n_roi),
  mask = TRUE,
  clinical_df = C_train,
  cov_df = cov_df_train,
  r_val = r_val,
  noise = 0.001,
  t2_lambda_0 = 1000,
  nu0 = 2, s20 = 1,
  eta0 = 2, t20 = 1,
  SGI_g1 = SGI_g1,
  k1 = 25,
  stepsize = 1e-3,
  acpt_step = 50,
  target_acpt = 0.25,
  tune = TRUE,
  fix_ls = fix_ls,
  rho0 = 2, psi20 = 1,
  gamma0 = 2, kappa20 = 1,
  omega0 = 2, phi20 = 1,
  k2 = 10,
  SGI_g2 = SGI_g2,
  g2_weighted = TRUE,
  seed = 12345 + chain_idx + 1000 * fold_id,
  burn_in = n_burnin,
  mcmc_sample = n_samps,
  chain_idx = chain_idx,
  save_path_g1 = fold_dir,
  save_path_g2 = fold_dir,
  save_path_MH = fold_dir,
  fname_suffix_g1 = paste("q", q_val, "fold", fold_id, sep = "_"),
  fname_suffix_g2 = paste("q", q_val, "fold", fold_id, sep = "_"),
  fname_suffix_MH = paste("q", q_val, "fold", fold_id, sep = "_"),
  run_g1 = TRUE, run_g2 = TRUE, run_MH = TRUE,
  load_g1 = FALSE, load_g2 = FALSE, load_MH = FALSE,
  save_g1 = TRUE, save_g2 = TRUE, save_MH = TRUE,
  dir_tags = c("donor_LOOCV", "q"),
  basename_tags = c("q", "chain")
)

cat("Finished donor-level LOOCV fit.\n")
