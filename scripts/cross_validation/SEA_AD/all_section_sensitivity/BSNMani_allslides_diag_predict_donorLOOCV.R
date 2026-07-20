#########################################
#### BSNMani diagnostics and donor-level held-out prediction
#### SEA-AD all-slide donor-level LOOCV version
####
#### Usage:
####   Rscript BSNMani_allslides_diag_predict_donorLOOCV.R <q_val> <fold_id> [cell_type] [n_chains]
#########################################

rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readxl)
  library(Rcpp)
  library(RcppArmadillo)
  library(coda)
  library(fs)
})

## =========================
## User-configurable paths
## =========================
BSNMANI_CODE_DIR <- "/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SEA-AD/BSNMani-dev"
DATA_DIR <- "/nfs/turbo/umms-lgarmire/liutong/BSNMani/SpaceX_BSNMani/data"
OUT_ROOT <- "/nfs/turbo/umms-lgarmire/liutong/BSNMani/SpaceX_BSNMani/BSNMani_output_result"
CLINICAL_FILE <- "/nfs/turbo/umms-lgarmire/liutong/BSNMani/SEA-AD Updated Data/meta_uni.xlsx"
SAMPLE_META_FILE <- "/nfs/turbo/umms-lgarmire/liutong/BSNMani/SEA-AD Updated Data/patient_sample_ID_all_slides.csv"

setwd(BSNMANI_CODE_DIR)

source("g1_diagnostics_helper.R")
sourceCpp("hybrid_M0_MALA_LR_FAST_v2.cpp")
sourceCpp("hybrid_M0_MALA_LR_g2.cpp")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript BSNMani_allslides_diag_predict_donorLOOCV.R <q_val> <fold_id> [cell_type] [n_chains]")
}

q_val <- as.numeric(args[1])
fold_id <- as.numeric(args[2])
cell_type <- ifelse(length(args) >= 3, args[3], "SpaceX_uni")
n_chains <- ifelse(length(args) >= 4, as.numeric(args[4]), 2)

SGI_g1 <- "GQN"
SGI_g2 <- "KPN"
r_val <- 2
n_burnin <- 30000

cat("Diagnostics/prediction | q =", q_val,
    "| fold =", fold_id,
    "| cell_type =", cell_type,
    "| n_chains =", n_chains, "\n")

extract_donor_id <- function(x) {
  m <- regexpr("H[0-9]+\\.[0-9]+\\.[0-9]+", x)
  ifelse(m > 0, regmatches(x, m), NA_character_)
}
safe_name <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)

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
  if (!"Atherosclerosis" %in% colnames(clinical_df)) stop("Cannot find Atherosclerosis column.")
  clinical_df$donor_id_std <- as.character(clinical_df[[donor_col]])
  clinical_df
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

make_FC_flat <- function(wgcna_ls) {
  n_pol <- length(wgcna_ls)
  n_roi <- nrow(wgcna_ls[[1]])
  FC_flat <- matrix(NA_real_, n_roi * n_roi, n_pol)
  for (i in seq_len(n_pol)) FC_flat[, i] <- c(as.matrix(wgcna_ls[[i]]))
  list(n_roi = n_roi, FC_flat = FC_flat)
}

project_lambda <- function(Y, U) {
  ## For Y ≈ U diag(lambda) U^T, lambda_k = u_k^T Y u_k.
  q <- ncol(U)
  out <- numeric(q)
  for (k in seq_len(q)) {
    u <- U[, k, drop = FALSE]
    out[k] <- as.numeric(t(u) %*% Y %*% u)
  }
  out
}

## =========================
## Reconstruct donor-level LOOCV split
## =========================
clinical_df_donor <- read_excel(CLINICAL_FILE)
wgcna_transformed_ls <- load_networks(cell_type)
sample_meta <- build_sample_meta(wgcna_transformed_ls)
if (any(is.na(sample_meta$donor_id))) stop("Some donor IDs could not be parsed/found.")

clinical_slide_df <- align_clinical_to_slides(sample_meta, clinical_df_donor)
unique_donors <- sort(unique(clinical_slide_df$donor_id))
heldout_donor <- unique_donors[fold_id]

train_idx <- which(clinical_slide_df$donor_id != heldout_donor)
test_idx <- which(clinical_slide_df$donor_id == heldout_donor)

fold_dir <- fs::path(
  OUT_ROOT,
  cell_type,
  "donor_LOOCV",
  paste0("fold_", sprintf("%02d", fold_id), "_", safe_name(heldout_donor)),
  paste0("q_", q_val)
)
diag_dir <- fs::path(fold_dir, "diagnostics")
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

cat("Held-out donor:", heldout_donor, "\n")
cat("Training slides:", length(train_idx), "| Test slides:", length(test_idx), "\n")

flat_train <- make_FC_flat(wgcna_transformed_ls[train_idx])
n_roi <- flat_train$n_roi
FC_flat_train <- flat_train$FC_flat
n_pol_train <- length(train_idx)

C_train <- as.matrix(clinical_slide_df[train_idx, "Cognitive_Status", drop = FALSE])
Z_train <- as.matrix(cbind(1, clinical_slide_df[train_idx, "Atherosclerosis", drop = FALSE]))

## =========================
## Load chains
## =========================
hybrid_res <- list()
fname_suffix_MH <- paste("q", q_val, "fold", fold_id, sep = "_")

for (c in seq_len(n_chains)) {
  f <- fs::path(
    fold_dir,
    "MH",
    paste("MH_res", SGI_g1, SGI_g2, "chain", c, fname_suffix_MH, sep = "_"),
    ext = "RDS"
  )
  if (!file.exists(f)) stop("Missing chain file: ", f)
  hybrid_res[[c]] <- list(mcmc = readRDS(f))
}

dat_ls <- list(
  N = n_roi,
  q = q_val,
  n_chains = n_chains,
  n_samps = length(hybrid_res[[1]]$mcmc$s2),
  n_burnin = n_burnin
)

cat("Formatting posterior samples...\n")

MALA_samps_ls <- gather_samples(
  mcmc_res_ls = hybrid_res,
  dat = dat_ls,
  var_df = data.frame(
    var = c("l_lambda", "Lambda_flat", "s2", "t2_lambda", "X", "U",
            "Y_C_llk", "Y_llk", "C_llk", "d", "t2_alpha", "t2_beta", "t2"),
    dim = c(1, 3, 1, 1, 3, 3, 1, 1, 1, 2, 1, 1, 1)
  ),
  sign_flip = FALSE
)

MALA_subset_ls <- subset_samps(
  samps_ls = MALA_samps_ls[c("Lambda_flat", "s2", "t2_lambda", "X", "U",
                             "Y_C_llk", "Y_llk", "C_llk", "d",
                             "t2_alpha", "t2_beta", "t2")],
  n_samps = dat_ls$n_samps,
  n_burnin = dat_ls$n_burnin
)
rm(MALA_samps_ls)

MALA_Rhat <- Rhat(MALA_subset_ls[c("Lambda_flat", "s2", "t2_lambda", "X", "U",
                                   "Y_C_llk", "Y_llk", "C_llk", "d",
                                   "t2_alpha", "t2_beta", "t2")])

MALA_pos_mean <- posterior_mean(MALA_subset_ls[c("Lambda_flat", "s2", "t2_lambda", "X", "U",
                                                 "Y_C_llk", "Y_llk", "C_llk", "d",
                                                 "t2_alpha", "t2_beta", "t2")])
MALA_pos_mean$dat_llk <- MALA_pos_mean$Y_C_llk

DIC_val <- DIC(
  mean_ls = MALA_pos_mean,
  Y_flat = t(FC_flat_train),
  C = C_train,
  Z = Z_train,
  M = n_pol_train,
  N = n_roi,
  q = q_val,
  r = r_val
)

## Posterior mean U and coefficients
X_est <- matrix(MALA_pos_mean$X, n_roi, q_val)
U_est <- polar_expansion(X_est)
d_mean <- as.matrix(MALA_pos_mean$d)

## Held-out section-level predictions
test_lambda <- matrix(NA_real_, nrow = length(test_idx), ncol = q_val)
for (j in seq_along(test_idx)) {
  Y_test <- as.matrix(wgcna_transformed_ls[[test_idx[j]]])
  test_lambda[j, ] <- project_lambda(Y_test, U_est)
}

test_design <- cbind(test_lambda, 1, clinical_slide_df[test_idx, "Atherosclerosis"])
test_score <- as.numeric(test_design %*% d_mean)
test_prob <- plogis(test_score)

section_pred_df <- data.frame(
  fold_id = fold_id,
  q_val = q_val,
  heldout_donor = heldout_donor,
  sample_id = clinical_slide_df$sample_id[test_idx],
  donor_id = clinical_slide_df$donor_id[test_idx],
  y_true = clinical_slide_df$Cognitive_Status[test_idx],
  atherosclerosis = clinical_slide_df$Atherosclerosis[test_idx],
  score = test_score,
  prob = test_prob,
  stringsAsFactors = FALSE
)

for (k in seq_len(q_val)) {
  section_pred_df[[paste0("lambda_", k)]] <- test_lambda[, k]
}

donor_pred_df <- section_pred_df %>%
  group_by(fold_id, q_val, donor_id) %>%
  summarise(
    y_true = unique(y_true)[1],
    n_sections = n(),
    score_mean = mean(score, na.rm = TRUE),
    prob_mean = mean(prob, na.rm = TRUE),
    .groups = "drop"
  )

diag <- list(
  posterior_mean = MALA_pos_mean,
  Rhat = MALA_Rhat,
  DIC = DIC_val,
  U_est = U_est,
  train_idx = train_idx,
  test_idx = test_idx,
  heldout_donor = heldout_donor,
  section_prediction = section_pred_df,
  donor_prediction = donor_pred_df
)

saveRDS(diag, file = fs::path(diag_dir, paste0("fold_", fold_id, "_q_", q_val, "_diag_pred.rds")))
write.csv(section_pred_df, file = fs::path(diag_dir, paste0("fold_", fold_id, "_q_", q_val, "_section_predictions.csv")), row.names = FALSE)
write.csv(donor_pred_df, file = fs::path(diag_dir, paste0("fold_", fold_id, "_q_", q_val, "_donor_prediction.csv")), row.names = FALSE)

cat("Saved diagnostics and donor-level held-out predictions.\n")
