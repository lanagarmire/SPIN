library(readxl)
library(car)
library(PRROC)

###############################
library(readxl)
library(car)
library(fs)


print(paste0("q value in clinical model construction is ", q_val))
print(paste0("cell type is ", cell_type))

################################
print(paste0("q value in clinical model construction is ", q_val))
print(paste0("Current threshold selection method is ", threshold_method))
clinical_df = read_excel("/nfs/turbo/umms-lgarmire/liutong/BSNMani/SEA-AD Updated Data/meta_uni.xlsx")

##

lambda_file_direction = fs::path("/nfs/turbo/umms-lgarmire/liutong/BSNMani/SpaceX_BSNMani/BSNMani_output_result",cell_type,paste("q",q_val,sep = "_"),"MH/diagnostics")
setwd(lambda_file_direction)
lambda_file = readRDS(paste0("MH_diag_res_GQN_KPN_q_",q_val,".RDS"))

lambda = lambda_file$posterior_mean$Lambda_flat


lambda_list = list()

for (i in 1:q_val) {
  lambda_list[[i]] = lambda[seq(i, length(lambda), by = q_val)]
}

lambda_df = as.data.frame(lambda_list)
colnames(lambda_df) = paste0("lambda_", 1:q_val)

clinical_df_final = cbind(clinical_df, lambda_df)


formula_all = paste0("Cognitive_Status ~ Atherosclerosis+",paste0("lambda_", seq(1, q_val), collapse = "+"))
formula_lambda_only = paste0("Cognitive_Status ~ ",paste0("lambda_", seq(1, q_val), collapse = "+"))
#if(cell_type == "IT" && q_val == 4) {formula = Cognitive_Status ~ lambda_1+lambda_2+lambda_4+Atherosclerosis}
## need to be changed according to heatmaps
##  if (q_val == 3)
 ## {formula_all = "Cognitive_Status ~ Atherosclerosis+lambda_1+lambda_3"}   ## lambda2 p value = 0.6
  #if (q_val == 4)
  #{formula_all = "Cognitive_Status ~ Atherosclerosis+lambda_1+lambda_2+lambda_3"}  ## lambda2 p value = 0.7
#  if (q_val == 5)
#  {formula_all = "Cognitive_Status ~ Atherosclerosis+lambda_2+lambda_3+lambda_4"}  ## lambda1 p value = 0.6
#  if (q_val == 2)
#  {formula_all = "Cognitive_Status ~ Atherosclerosis+lambda_1+lambda_2"}

model_lambda_all = glm(as.formula(formula_all), data = clinical_df_final, family = "binomial")
model_lambda_only = glm(as.formula(formula_lambda_only), data = clinical_df_final, family = "binomial")
print(summary(model_lambda_all))
print("\n")
print(summary(model_lambda_only))

print("check multicollinearity\n")
print(vif(model_lambda_only))

cat("################Begin LOOCV################\n")

##################
## using LOOCV algorithm for logistic regression construction
#################
library(pscl)
library(pROC)
library(caret)

n = nrow(clinical_df_final)
errors_all = numeric(n)
r2 = numeric(n)
aic = numeric(n)

predicted_probs = numeric(n)
predicted_labels = numeric(n)
true_labels = numeric(n)
test_patient_ids = character(n)

for (i in 1:n) {
  train_data = clinical_df_final[-i, ]
  test_data = clinical_df_final[i, , drop = FALSE]
  
  formula_all = paste0("Cognitive_Status ~ Atherosclerosis + ", 
                       paste0("lambda_", seq(1, q_val), collapse = " + "))

  model_all = glm(as.formula(formula_all), data = train_data, family = "binomial")
  
  pred_prob_all = predict(model_all, newdata = test_data, type = "response")
  #pred_class_all = ifelse(pred_prob_all > predictive_threshold, 1, 0)
  pred_class_all = NA
  
  predicted_probs[i] = pred_prob_all
  #redicted_labels[i] = pred_class_all
  true_labels[i] = test_data$Cognitive_Status
  test_patient_ids[i] = as.character(test_data[[1]])

  #errors_all[i] = ifelse(pred_class_all != test_data$Cognitive_Status, 1, 0)
  r2[i] = pR2(model_all)[["McFadden"]]
  aic[i] = AIC(model_all)
}

#mean_error_all = mean(errors_all)
mean_r2 = mean(r2)
mean_aic = mean(aic)

true_labels = as.numeric(true_labels)
predicted_probs = as.numeric(predicted_probs)

if (threshold_method == "fixed_0.5") {
  
  predictive_threshold = fixed_threshold
  
} else if (threshold_method == "prevalence") {
  
  predictive_threshold = prevalence_threshold
  
} else if (threshold_method == "roc") {
  
  roc_for_threshold = roc(
    response = true_labels,
    predictor = predicted_probs,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
  
  best_threshold_info = coords(
    roc_for_threshold,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )
  
  predictive_threshold = as.numeric(best_threshold_info["threshold"])
  
} else {
  
  stop("Invalid threshold_method. Please use 'fixed_0.5', 'prevalence', or 'roc'.")
}

cat(sprintf("Selected threshold method = %s\n", threshold_method))
cat(sprintf("Selected predictive threshold = %.6f\n", predictive_threshold))

predicted_labels = ifelse(predicted_probs >= predictive_threshold, 1, 0)
errors_all = ifelse(predicted_labels != true_labels, 1, 0)

mean_error_all = mean(errors_all)

loocv_prediction_df = data.frame(
  patient_id = test_patient_ids,
  true_label = as.numeric(true_labels),
  predicted_prob = as.numeric(predicted_probs),
  predicted_label = as.numeric(predicted_labels),
  error = as.numeric(errors_all),
  threshold_method = threshold_method,
  predictive_threshold = predictive_threshold,
  stringsAsFactors = FALSE
)

write.csv(
  loocv_prediction_df,
  #paste0("LOOCV_prediction_results_q_", q_val, "_", cell_type, ".csv"),
  paste0("LOOCV_prediction_results_q_", q_val, "_", cell_type, "_", threshold_method, ".csv"),
  row.names = FALSE
)
# conf_matrix = confusionMatrix(
#   factor(predicted_labels),
#   factor(true_labels),
#   positive = "1"
# )

# precision = conf_matrix$byClass["Precision"]
# recall = conf_matrix$byClass["Recall"]
# f1 = conf_matrix$byClass["F1"]
# specificity = conf_matrix$byClass["Specificity"]
# accuracy = conf_matrix$overall["Accuracy"]

# roc_obj = roc(true_labels, predicted_probs)
# auc_val = auc(roc_obj)
library(PRROC)

## make sure labels are 0/1 numeric
true_labels = as.numeric(true_labels)
predicted_labels = as.numeric(predicted_labels)
predicted_probs = as.numeric(predicted_probs)

## confusion matrix with fixed level order
conf_matrix = confusionMatrix(
  factor(predicted_labels, levels = c(0, 1)),
  factor(true_labels, levels = c(0, 1)),
  positive = "1"
)

## basic metrics from caret
precision = unname(conf_matrix$byClass["Precision"])
recall = unname(conf_matrix$byClass["Recall"])
f1 = unname(conf_matrix$byClass["F1"])
specificity = unname(conf_matrix$byClass["Specificity"])
accuracy = unname(conf_matrix$overall["Accuracy"])

## extract confusion matrix counts
cm_table = conf_matrix$table
TN = cm_table["0", "0"]
FP = cm_table["1", "0"]
FN = cm_table["0", "1"]
TP = cm_table["1", "1"]

## MCC
mcc_denominator = sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
mcc = ifelse(mcc_denominator == 0, NA,
             (TP * TN - FP * FN) / mcc_denominator)

## Macro F1
## F1 for class 1
precision_pos = ifelse((TP + FP) == 0, NA, TP / (TP + FP))
recall_pos    = ifelse((TP + FN) == 0, NA, TP / (TP + FN))
f1_pos = ifelse(is.na(precision_pos) || is.na(recall_pos) || (precision_pos + recall_pos) == 0,
                NA,
                2 * precision_pos * recall_pos / (precision_pos + recall_pos))

## F1 for class 0
precision_neg = ifelse((TN + FN) == 0, NA, TN / (TN + FN))
recall_neg    = ifelse((TN + FP) == 0, NA, TN / (TN + FP))
f1_neg = ifelse(is.na(precision_neg) || is.na(recall_neg) || (precision_neg + recall_neg) == 0,
                NA,
                2 * precision_neg * recall_neg / (precision_neg + recall_neg))

macro_f1 = mean(c(f1_pos, f1_neg), na.rm = TRUE)

## AUROC
#roc_obj = roc(response = true_labels, predictor = predicted_probs, quiet = TRUE)
roc_obj = roc(
  response = true_labels,
  predictor = predicted_probs,
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)
auroc_val = as.numeric(auc(roc_obj))

## AUPRC
## scores.class0 should be scores for positive class
pr_obj = pr.curve(
  scores.class0 = predicted_probs[true_labels == 1],
  scores.class1 = predicted_probs[true_labels == 0],
  curve = TRUE
)
auprc_val = pr_obj$auc.integral

cat("========== LOOCV Logistic Regression Evaluation ==========\n")
cat(sprintf("Mean Classification Error = %.4f\n", mean_error_all))
cat(sprintf("Mean McFadden R²          = %.4f\n", mean_r2))
cat(sprintf("Mean AIC                  = %.4f\n\n", mean_aic))

cat(sprintf("Accuracy     = %.4f\n", accuracy))
cat(sprintf("AUROC        = %.4f\n", auroc_val))
cat(sprintf("AUPRC        = %.4f\n", auprc_val))
cat(sprintf("MCC          = %.4f\n", mcc))
cat(sprintf("F1 Score     = %.4f\n", f1))
cat(sprintf("Macro F1     = %.4f\n", macro_f1))
cat(sprintf("Precision    = %.4f\n", precision))
cat(sprintf("Recall       = %.4f\n", recall))
cat(sprintf("Specificity  = %.4f\n", specificity))

############################################
## Bootstrap based on 27 LOOCV test results
############################################

set.seed(123)
B = 5000. # previously 3000, now 5000

calculate_metrics_from_predictions = function(true_labels, predicted_probs, threshold) {
  
  true_labels = as.numeric(true_labels)
  predicted_probs = as.numeric(predicted_probs)
  predicted_labels = ifelse(predicted_probs >= threshold, 1, 0)
  
  conf_matrix = confusionMatrix(
    factor(predicted_labels, levels = c(0, 1)),
    factor(true_labels, levels = c(0, 1)),
    positive = "1"
  )
  
  accuracy = unname(conf_matrix$overall["Accuracy"])
  precision = unname(conf_matrix$byClass["Precision"])
  recall = unname(conf_matrix$byClass["Recall"])
  f1 = unname(conf_matrix$byClass["F1"])
  specificity = unname(conf_matrix$byClass["Specificity"])
  
  cm_table = conf_matrix$table
  TN = cm_table["0", "0"]
  FP = cm_table["1", "0"]
  FN = cm_table["0", "1"]
  TP = cm_table["1", "1"]
  
  mcc_denominator = sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc = ifelse(mcc_denominator == 0, NA,
               (TP * TN - FP * FN) / mcc_denominator)
  
  precision_pos = ifelse((TP + FP) == 0, NA, TP / (TP + FP))
  recall_pos = ifelse((TP + FN) == 0, NA, TP / (TP + FN))
  f1_pos = ifelse(is.na(precision_pos) || is.na(recall_pos) || 
                    (precision_pos + recall_pos) == 0,
                  NA,
                  2 * precision_pos * recall_pos / (precision_pos + recall_pos))
  
  precision_neg = ifelse((TN + FN) == 0, NA, TN / (TN + FN))
  recall_neg = ifelse((TN + FP) == 0, NA, TN / (TN + FP))
  f1_neg = ifelse(is.na(precision_neg) || is.na(recall_neg) || 
                    (precision_neg + recall_neg) == 0,
                  NA,
                  2 * precision_neg * recall_neg / (precision_neg + recall_neg))
  
  macro_f1 = mean(c(f1_pos, f1_neg), na.rm = TRUE)
  
  auroc_val = NA
  if (length(unique(true_labels)) == 2) {
    roc_obj = roc(response = true_labels, predictor = predicted_probs, quiet = TRUE)
    auroc_val = as.numeric(auc(roc_obj))
  }
  
  auprc_val = NA
  if (sum(true_labels == 1) > 0 && sum(true_labels == 0) > 0) {
    pr_obj = pr.curve(
      scores.class0 = predicted_probs[true_labels == 1],
      scores.class1 = predicted_probs[true_labels == 0],
      curve = FALSE
    )
    auprc_val = pr_obj$auc.integral
  }
  
  mean_classification_error = mean(predicted_labels != true_labels)
  
  return(c(
    classification_error = mean_classification_error,
    accuracy = accuracy,
    auroc = auroc_val,
    auprc = auprc_val,
    mcc = mcc,
    f1 = f1,
    macro_f1 = macro_f1,
    precision = precision,
    recall = recall,
    specificity = specificity
  ))
}

bootstrap_metrics = matrix(
  NA,
  nrow = B,
  ncol = 10
)

colnames(bootstrap_metrics) = c(
  "classification_error",
  "accuracy",
  "auroc",
  "auprc",
  "mcc",
  "f1",
  "macro_f1",
  "precision",
  "recall",
  "specificity"
)

for (b in 1:B) {
  
  boot_index = sample(seq_len(nrow(loocv_prediction_df)),
                      size = nrow(loocv_prediction_df),
                      replace = TRUE)
  
  boot_data = loocv_prediction_df[boot_index, ]
  
  bootstrap_metrics[b, ] = calculate_metrics_from_predictions(
    true_labels = boot_data$true_label,
    predicted_probs = boot_data$predicted_prob,
    threshold = predictive_threshold
  )
}

bootstrap_metrics_df = as.data.frame(bootstrap_metrics)

bootstrap_summary_df = data.frame(
  metric = colnames(bootstrap_metrics_df),
  mean = sapply(bootstrap_metrics_df, mean, na.rm = TRUE),
  variance = sapply(bootstrap_metrics_df, var, na.rm = TRUE),
  sd = sapply(bootstrap_metrics_df, sd, na.rm = TRUE),
  min = sapply(bootstrap_metrics_df, min, na.rm = TRUE),
  max = sapply(bootstrap_metrics_df, max, na.rm = TRUE),
  q025 = sapply(bootstrap_metrics_df, quantile, probs = 0.025, na.rm = TRUE),
  q975 = sapply(bootstrap_metrics_df, quantile, probs = 0.975, na.rm = TRUE),
  row.names = NULL
)

write.csv(
  bootstrap_metrics_df,
  #paste0("Bootstrap_3000_all_metrics_q_", q_val, "_", cell_type, ".csv"),
  paste0("Bootstrap_",B,"_all_metrics_q_", q_val, "_", cell_type, "_", threshold_method, ".csv"),
  row.names = FALSE
)

write.csv(
  bootstrap_summary_df,
  #paste0("Bootstrap_3000_summary_q_", q_val, "_", cell_type, ".csv"),
  paste0("Bootstrap_",B,"_all_metrics_q_", q_val, "_", cell_type, "_", threshold_method, ".csv"),
  row.names = FALSE
)

cat("\n========== Bootstrap Summary ==========\n")
print(paste0("The number of boostrap is ",B))
print(bootstrap_summary_df)

