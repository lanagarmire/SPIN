#########################################
#### Collect donor-level LOOCV predictions across folds
####
#### Usage:
####   Rscript collect_donorLOOCV_predictions.R <q_val> [cell_type]
#########################################

rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript collect_donorLOOCV_predictions.R <q_val> [cell_type]")

q_val <- as.numeric(args[1])
cell_type <- ifelse(length(args) >= 2, args[2], "SpaceX_uni")

OUT_ROOT <- "/nfs/turbo/umms-lgarmire/liutong/BSNMani/SpaceX_BSNMani/BSNMani_output_result"
root <- file.path(OUT_ROOT, cell_type, "donor_LOOCV")

files <- list.files(
  root,
  pattern = paste0("fold_.*_q_", q_val, "_donor_prediction.csv$"),
  recursive = TRUE,
  full.names = TRUE
)

if (length(files) == 0) stop("No donor prediction files found.")

pred <- bind_rows(lapply(files, fread))

out_dir <- file.path(root, paste0("summary_q_", q_val))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(pred, file = file.path(out_dir, paste0("all_donor_predictions_q_", q_val, ".csv")), row.names = FALSE)

y <- as.numeric(pred$y_true)
score <- as.numeric(pred$prob_mean)
yhat <- as.numeric(score >= 0.5)

tp <- sum(yhat == 1 & y == 1)
tn <- sum(yhat == 0 & y == 0)
fp <- sum(yhat == 1 & y == 0)
fn <- sum(yhat == 0 & y == 1)

accuracy <- mean(yhat == y)
precision <- ifelse((tp + fp) == 0, NA, tp / (tp + fp))
recall <- ifelse((tp + fn) == 0, NA, tp / (tp + fn))
specificity <- ifelse((tn + fp) == 0, NA, tn / (tn + fp))
f1_pos <- ifelse(is.na(precision + recall) || (precision + recall) == 0, NA, 2 * precision * recall / (precision + recall))

precision0 <- ifelse((tn + fn) == 0, NA, tn / (tn + fn))
recall0 <- specificity
f1_neg <- ifelse(is.na(precision0 + recall0) || (precision0 + recall0) == 0, NA, 2 * precision0 * recall0 / (precision0 + recall0))
macro_f1 <- mean(c(f1_pos, f1_neg), na.rm = TRUE)

denom <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
mcc <- ifelse(denom == 0, NA, (tp * tn - fp * fn) / denom)

## AUROC without requiring extra packages
calc_auc <- function(labels, scores) {
  o <- order(scores, decreasing = TRUE)
  labels <- labels[o]
  pos <- sum(labels == 1)
  neg <- sum(labels == 0)
  if (pos == 0 || neg == 0) return(NA_real_)
  ranks <- rank(scores)
  (sum(ranks[labels == 1]) - pos * (pos + 1) / 2) / (pos * neg)
}
auroc <- calc_auc(y, score)

metrics <- data.frame(
  q_val = q_val,
  n_donors = length(unique(pred$donor_id)),
  accuracy = accuracy,
  AUROC = auroc,
  MCC = mcc,
  macro_F1 = macro_f1,
  precision = precision,
  recall = recall,
  specificity = specificity
)

write.csv(metrics, file = file.path(out_dir, paste0("donor_LOOCV_metrics_q_", q_val, ".csv")), row.names = FALSE)
print(metrics)
