#' Print table of C-index and log rank p for WGCNA, hdWGCNA, Smoothie

q_vals <- 2:8
seeds <- c(0, 42, 64, 123, 894)
output_dirs <- list(
  "WGCNA" = "/nfs/dcmb-lgarmire/yangiwen/workspace/bsnmani/output/bc/output_23",
  "hdWGCNA" = "/nfs/dcmb-lgarmire/yangiwen/workspace/bsnmani/output/bc/output_21",
  "Smoothie" = "/nfs/dcmb-lgarmire/yangiwen/workspace/bsnmani/output/bc/output_12"
)

#' Read BSNMani prediction, calculate c-index and log rank p
loadResults <- function(output_dir, q_val, seed, method, n_group = 3) {
  output_dir <- file.path(output_dir, paste0("q", q_val), paste0("s", seed))
  # Read train/test set and prediction
  clinical_df_train <- readRDS(file.path(output_dir, "train", "g2_survival", "clinical_df_train_train.RDS"))
  clinical_df_test <- readRDS(file.path(output_dir, "test", "g2_survival", "clinical_df_test_test.RDS"))
  clinical_df_train$pred <- readRDS(file.path(output_dir, "train", "g2_survival", "lp_train_train.RDS"))
  clinical_df_test$pred <- readRDS(file.path(output_dir, "test", "g2_survival", "lp_test_test.RDS"))
  # Calculate C-Index
  surv_obj_test <- survival::Surv(clinical_df_test$OSmonth, clinical_df_test$Patientstatus)
  cindex_test <- survival::concordance(surv_obj_test ~ pred, data = clinical_df_test, reverse = T)$concordance
  # Assign risk groups
  clinical_df_train$group <- if (n_group == 2) {
    ifelse(
      clinical_df_train$pred > quantile(clinical_df_train$pred, 0.5),
      "high",
      "low"
    )
  } else if (n_group == 3) {
    cut(
      clinical_df_train$pred,
      breaks = quantile(clinical_df_train$pred, probs = c(0, 1 / 3, 2 / 3, 1)),
      labels = c("low", "mid", "high"),
      include.lowest = T
    )
  } else {
    stop("Invalid number of risk groups")
  }
  # Calculate pval
  pval <- survival::survdiff(survival::Surv(OSmonth, Patientstatus) ~ group, data = clinical_df_train)$pvalue
  data.frame(q = q_val, seed = seed, cindex = cindex_test, method = method, pval = pval)
}

# Read all results
results_list <- list()
for (method in names(output_dirs)) {
  for (q_val in q_vals) {
    for (seed in seeds) {
      res <- loadResults(output_dirs[[method]], q_val, seed, method)
      results_list[[length(results_list) + 1]] <- res
    }
  }
}
results <- as.data.frame(do.call(rbind, results_list))

# Aggregate cindex and pval separately
agg_cindex <- aggregate(cindex ~ method + q, data = results, FUN = mean)
agg_pval <- aggregate(pval ~ method + q, data = results, FUN = function(x) exp(mean(log(as.numeric(x)))))
agg_df <- merge(agg_cindex, agg_pval, by = c("method", "q"))
write.csv(agg_df, "/nfs/dcmb-lgarmire/yangiwen/workspace/bsnmani/output/bc/benchmark_spatial_coexpr.csv")
