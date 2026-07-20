#!/usr/bin/env Rscript

# Ensure correct libraries on the cluster
.libPaths(c("~/R_libs_amd", .libPaths()))

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(tidyr) # Added tidyr to pivot the table into wide format
})

cat("========================================================\n")
cat(" Evaluating C-Index (Wide Format: Seeds as Columns)\n")
cat("========================================================\n\n")

# --- 1. Define Paths and Iterators ---
base_dir <- "/home/koe/BSNMani_application-main/Github_BSNMani/IMC/bc"

# >>> SET THIS BEFORE RUNNING <<<
cell_type = "Smoothie"

# Define the Smoothie output folders you want to evaluate.
# (Based on your image, output_10 to output_14 are the 5-fold Smoothie folders. 
# You can also add output_22 here if you want to reproduce your specific screenshot).
fold_folders <- c(
  "Fold_1" = "output_10",
  "Fold_2" = "output_11",
  "Fold_3" = "output_12",
  "Fold_4" = "output_13",
  "Fold_5" = "output_14"
)

# If you only want to test output_22 like in your screenshot, uncomment this:
# fold_folders <- c("Fold_3_Norm" = "output_22")

q_folders <- paste0("q", 2:8)
seed_folders <- c("s0", "s123", "s42", "s64", "s894")

results_list <- list()

# --- 2. Master Loop: Fold -> q -> Seed ---
for (fold_name in names(fold_folders)) {
  output_folder <- fold_folders[[fold_name]]
  cat(sprintf("Extracting raw data for %s (%s)...\n", fold_name, output_folder))
  
  for (q in q_folders) {
    for (seed in seed_folders) {
      
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
      
      # Ensure numeric survival variables
      clinical_df$OSmonth <- as.numeric(clinical_df$OSmonth)
      clinical_df$Patientstatus <- as.numeric(clinical_df$Patientstatus)
      
      # Calculate C-Index
      cindex_obj <- survival::concordance(
        Surv(OSmonth, Patientstatus) ~ lp_test,
        data = clinical_df,
        reverse = TRUE
      )
      
      # Store the result in a long format
      results_list[[length(results_list) + 1]] <- data.frame(
        Output_Folder = output_folder,
        q_value = as.numeric(gsub("q", "", q)), # Convert "q2" to numeric 2 for clean sorting
        Seed = seed,
        Test_C_Index = cindex_obj$concordance
      )
    }
  }
}

# --- 3. Pivot to Wide Format and Calculate Averages ---
if (length(results_list) > 0) {
  all_results_df <- bind_rows(results_list)
  
  # Reshape data: q_value as rows, Seeds as columns, and calculate the average
  wide_summary_df <- all_results_df %>%
    pivot_wider(
      names_from = Seed, 
      values_from = Test_C_Index
    ) %>%
    # Ensure columns are in the specific order from your screenshot
    select(Output_Folder, q = q_value, any_of(c("s0", "s123", "s42", "s64", "s894"))) %>%
    rowwise() %>%
    # Calculate the row average across all seed columns
    mutate(Average = mean(c_across(starts_with("s")), na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(Output_Folder, q)
  
  # --- 4. Print Grids to Console and Save ---
  save_dir <- "/home/koe/BSNMani_application-main/Github_BSNMani/IMC/Results"
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Loop through unique output folders to print them beautifully (like your screenshot)
  unique_folders <- unique(wide_summary_df$Output_Folder)
  
  for (folder in unique_folders) {
    cat("\n====================================================================\n")
    cat(sprintf(" %s | %s Output\n", folder, cell_type))
    cat("====================================================================\n")
    
    # Filter for just this folder and remove the Output_Folder column for clean printing
    folder_df <- wide_summary_df %>% 
      filter(Output_Folder == folder) %>%
      select(-Output_Folder)
    
    # Round numbers to 4 decimal places for clean display
    print(folder_df %>% mutate(across(where(is.numeric), ~round(.x, 4))))
    
    # Save individual CSVs for each output folder
    folder_file <- file.path(save_dir, paste0(cell_type, "_", folder, "_seed_grid.csv"))
    write.csv(folder_df, file = folder_file, row.names = FALSE)
  }
  
  cat("\n====================================================================\n")
  cat(sprintf("[Success] Individual grid CSVs saved to:\n   %s\n", save_dir))
  
} else {
  cat("\n[Error] No valid results were found. Please verify the folder structure.\n")
}