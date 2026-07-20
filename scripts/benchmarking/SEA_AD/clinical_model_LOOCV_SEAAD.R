library(readxl)
library(car)

print(paste0("q value in clinical model construction is ", q_val))

## check point 1                                     ####
clinical_df = read_excel("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SpaceX_BSNMani/data/meta.xlsx")
##

lambda_file_direction = fs::path("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SpaceX_BSNMani/result",cell_type,paste("q",q_val,sep = "_"),"MH/diagnostics")
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

for (i in 1:n) {
  train_data = clinical_df_final[-i, ]
  test_data = clinical_df_final[i, , drop = FALSE]
  
  formula_all = paste0("Cognitive_Status ~ Atherosclerosis + ", 
                       paste0("lambda_", seq(1, q_val), collapse = " + "))

  model_all = glm(as.formula(formula_all), data = train_data, family = "binomial")
  
  pred_prob_all = predict(model_all, newdata = test_data, type = "response")
  pred_class_all = ifelse(pred_prob_all > 0.5, 1, 0)
  
  predicted_probs[i] = pred_prob_all
  predicted_labels[i] = pred_class_all
  true_labels[i] = test_data$Cognitive_Status
  
  errors_all[i] = ifelse(pred_class_all != test_data$Cognitive_Status, 1, 0)
  r2[i] = pR2(model_all)[["McFadden"]]
  aic[i] = AIC(model_all)
}

mean_error_all = mean(errors_all)
mean_r2 = mean(r2)
mean_aic = mean(aic)

conf_matrix = confusionMatrix(
  factor(predicted_labels),
  factor(true_labels),
  positive = "1"
)

precision = conf_matrix$byClass["Precision"]
recall = conf_matrix$byClass["Recall"]
f1 = conf_matrix$byClass["F1"]
specificity = conf_matrix$byClass["Specificity"]
accuracy = conf_matrix$overall["Accuracy"]

roc_obj = roc(true_labels, predicted_probs)
auc_val = auc(roc_obj)


cat("========== LOOCV Logistic Regression Evaluation ==========\n")
cat(sprintf("Mean Classification Error = %.4f\n", mean_error_all))
cat(sprintf("Mean McFadden R²          = %.4f\n", mean_r2))
cat(sprintf("Mean AIC                  = %.4f\n\n", mean_aic))

cat(sprintf("Accuracy     = %.4f\n", accuracy))
cat(sprintf("Precision    = %.4f\n", precision))
cat(sprintf("Recall       = %.4f\n", recall))
cat(sprintf("F1 Score     = %.4f\n", f1))
cat(sprintf("Specificity  = %.4f\n", specificity))
cat(sprintf("AUC          = %.4f\n", auc_val))


