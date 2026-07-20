#!/usr/bin/env Rscript

# Ensure correct libraries on the cluster
.libPaths(c("~/R_libs_amd", .libPaths()))

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(boot) # Added for bootstrapping
})

cat("========================================================\n")
cat(" Evaluating C-Index with Bootstrap CIs (Single Output)\n")
cat("========================================================\n\n")

# --- 1. Define Paths and Iterators ---
base_dir <- "/home/koe/BSNMani_application-main/Github_BSNMani/IMC/bc"

# >>> SET THIS BEFORE RUNNING <<<
# cell_type = "hdWGCNA"
# cell_type = "Smoothie"
cell_type = "WGCNA"

# output_folder <- "output_21" # Single target hdWGCNA output folder
# output_folder <- "output_22" # Single target Smoothie output folder
output_folder <- "output_23" # Single target WGCNA output folder

q_folders <- paste0("q", 2:8)
seed_folders <- c("s0", "s123", "s42", "s64", "s894")

# Number of bootstrap resamples 
B_resamples <- 2000 

# Define the bootstrap function for C-Index
cindex_boot_fn <- function(data, indices) {
  d <- data[indices, ]
  # Use tryCatch to prevent errors if a bootstrap sample randomly draws only one patient status
  fit <- tryCatch({
    survival::concordance(Surv(OSmonth, Patientstatus) ~ lp_test, data = d, reverse = TRUE)$concordance
  }, error = function(e) NA)
  return(fit)
}

summary_results_list <- list()

# --- 2. Master Loop: q -> Seed ---
for (q in q_folders) {
  for (seed in seed_folders) {
    
    cat(sprintf("\n--- Processing Data for %s | %s ---\n", q, seed))
    
    target_dir <- file.path(base_dir, output_folder, q, seed, "test", "g2_survival")
    
    clinical_file <- file.path(target_dir, "clinical_df_test_test.RDS")
    lp_file <- file.path(target_dir, "lp_test_test.RDS")
    
    # Safety Check
    if (!file.exists(clinical_file) || !file.exists(lp_file)) {
      message(sprintf("  [Warning] Missing files in %s/%s/%s. Skipping...", output_folder, q, seed))
      next
    }
    
    # Load Data
    clinical_df <- readRDS(clinical_file)
    lp_test <- readRDS(lp_file)
    
    # Ensure numeric formats
    clinical_df$OSmonth <- as.numeric(clinical_df$OSmonth)
    clinical_df$Patientstatus <- as.numeric(clinical_df$Patientstatus)
    clinical_df$lp_test <- as.numeric(lp_test)
    
    # --- 3. Bootstrap the Test Predictions ---
    if (nrow(clinical_df) > 0) {
      cat(sprintf("▶️ Bootstrapping predictions (N = %d patients, B = %d)...\n", nrow(clinical_df), B_resamples))
      
      # Calculate the original C-Index on the test set
      cindex_overall <- survival::concordance(Surv(OSmonth, Patientstatus) ~ lp_test, data = clinical_df, reverse = TRUE)$concordance
      
      # Run Bootstrap
      boot_out <- boot(data = clinical_df, statistic = cindex_boot_fn, R = B_resamples)
      
      # Calculate 95% Percentile Confidence Intervals
      ci <- tryCatch({ boot.ci(boot_out, type = "perc", conf = 0.95) }, error = function(e) NULL)
      
      ci_lower <- ifelse(!is.null(ci), ci$percent[4], NA)
      ci_upper <- ifelse(!is.null(ci), ci$percent[5], NA)
      
      # Store the summarized bootstrap results
      summary_results_list[[length(summary_results_list) + 1]] <- data.frame(
        Seed = seed,
        q_value = q,
        N_Patients = nrow(clinical_df),
        Test_C_Index = round(cindex_overall, 4),
        Boot_Mean_C_Index = round(mean(boot_out$t, na.rm = TRUE), 4),
        CI_95_Lower = round(ci_lower, 4),
        CI_95_Upper = round(ci_upper, 4)
      )
    }
  }
}

# --- 4. Aggregate and Print ---
if (length(summary_results_list) > 0) {
  
  summary_df <- bind_rows(summary_results_list) %>% arrange(Seed, q_value)
  
  cat("\n========================================================\n")
  cat(" FINAL SUMMARY: C-Index and 95% Bootstrap CIs\n")
  cat("========================================================\n")
  print(as.data.frame(summary_df))
  cat("========================================================\n\n")
  
  # --- 5. Save Results ---
  save_dir <- "/home/koe/BSNMani_application-main/Github_BSNMani/IMC/Results"
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  summary_file <- file.path(save_dir, paste0(cell_type, "_single_run_bootstrap_summary.csv"))
  
  write.csv(summary_df, file = summary_file, row.names = FALSE)
  
  cat(sprintf("[Success] Bootstrap summary table saved to:\n   %s\n", summary_file))
  
} else {
  cat("\n[Error] No valid results were found. Please verify the folder structure.\n")
}