#' Run baseline model on clinical only

source("/nfs/dcmb-lgarmire/yangiwen/workspace/common/Utils.R")

# Create output dir
output_dir <- "/nfs/dcmb-lgarmire/yangiwen/workspace/bsnmani/output/bc/baseline"
mkdir(output_dir)
# Read clinical df inputs
clinical_df <- readRDS("/nfs/dcmb-lgarmire/haowenwu/clinical_df_3cov.RDS")
folds <- readRDS("/nfs/dcmb-lgarmire/haowenwu/5CV_folds.RDS")
fold <- 3
clinical_df$clinical_type <- factor(clinical_df$clinical_type, levels = c("HR+HER2-", "HR+HER2+", "HR-HER2+", "HR-HER2-"))
clinical_df$grade <- factor(clinical_df$grade, levels = c(1, 2, 3))
clinical_df_train <- clinical_df[clinical_df$core %in% folds[[fold]],]

# Train model on clinical variables only
formula <- as.formula("survival::Surv(OSmonth, Patientstatus) ~ age + grade + clinical_type")
model <- survival::coxph(formula, data = clinical_df_train)

# Run prediction for train
clinical_df_train$pred <- predict(model, clinical_df_train, type = "lp")

# Save result
saveRDS(model, file.path(output_dir, "model.rds"))
saveRDS(clinical_df_train, file.path(output_dir, "clinical_df_train.rds"))
