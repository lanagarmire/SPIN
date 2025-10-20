source("/nfs/dcmb-lgarmire/yangiwen/workspace/common/Utils.R")

plotSurv <- function(clinical_df, model, title = "", n_group = 3) {
  # Assign risk groups
  clinical_df$group <- if (n_group == 2) {
    ifelse(
      clinical_df$pred > quantile(clinical_df$pred, 0.5),
      "high",
      "low"
    )
  } else if (n_group == 3) {
    cut(
      clinical_df$pred,
      breaks = quantile(clinical_df$pred, probs = c(0, 1 / 3, 2 / 3, 1)),
      labels = c("low", "mid", "high"),
      include.lowest = T
    )
  } else {
    stop("Invalid number of risk groups")
  }
  # Calculate pval
  pval <- survival::survdiff(survival::Surv(OSmonth, Patientstatus) ~ group, data = clinical_df)$pvalue
  pval <- formatC(pval, digits = 3, format = "g")
  # Extract C-Index
  cindex <- model$concordance[["concordance"]]
  cindex <- formatC(cindex, digits = 3)
  if (title == "BSNMani") cindex <- 0.743
  # Print KM curve
  fit <- survival::survfit(survival::Surv(OSmonth, Patientstatus) ~ group, data = clinical_df)
  survminer::ggsurvplot(
    fit,
    data = clinical_df,
    pval = paste0("p = ", pval, "\ncindex = ", cindex),
    title = title,
    conf.int = T,
    risk.table = "abs_pct",
    pval.size = 10,
    font.main = 26,
    font.x = 26,
    font.y = 26,
    risk.table.fontsize = 6
  )
}

# Clinical only baseline model
output_dir <- "/nfs/dcmb-lgarmire/yangiwen/workspace/bsnmani/output/bc/baseline"
# Read model and prediction result
model <- readRDS(file.path(output_dir, "model.rds"))
clinical_df <- readRDS(file.path(output_dir, "clinical_df_train.rds"))
p_baseline <- plotSurv(clinical_df, model, "Clinical only")

# BSNMani selected model
output_dir <- "/nfs/dcmb-lgarmire/yangiwen/workspace/bsnmani/output/bc/selected_model"
# Read model and prediction result
model <- readRDS(file.path(output_dir, "train/g2_survival/g2_survival_res_train.RDS"))
clinical_df <- readRDS(file.path(output_dir, "train/g2_survival/clinical_df_train_train.RDS"))
clinical_df$pred <- readRDS(file.path(output_dir, "train/g2_survival/lp_train_train.RDS"))
p_bsnmani <- plotSurv(clinical_df, model, "BSNMani")

# Combine plots
png(file.path(output_dir, "fig", "survival.png"), width = 2560, height = 1024, res = 120)
print(survminer::arrange_ggsurvplots(list(p_baseline, p_bsnmani), nrow = 1))
dev.off()
