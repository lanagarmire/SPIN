rm(list = ls())
library(DescTools)

## SpaceX_res_all_3celltype_update file
setwd("/home/koe/BSNMani_application-main/Tong_results/BSNMani_for_transfer/BSNMani_for_transfer/BSNMani/SpaceX_BSNMani/data")

## read 4 data
baseline = readRDS("Baseline_co_expression_list_uni.RDS")
hdwgcna = readRDS("hdWGCNA_co_expression_list_uni.RDS")
spacex = readRDS("SpaceX_co_expression_list_uni.RDS")
smoothie= readRDS("Smoothie_co_expression_list_uni.RDS")

 
baseline_filter = baseline
hdwgcna_filter = hdwgcna
spacex_filter = spacex
smoothie_filter = smoothie


safe_fisherz <- function(mat) {
  FisherZ(mat)
}
baseline_z  <- lapply(baseline_filter, safe_fisherz)
hdwgcna_z   <- lapply(hdwgcna_filter, safe_fisherz)
spacex_z    <- lapply(spacex_filter, safe_fisherz)
smoothie_z  <- lapply(smoothie_filter, safe_fisherz)

set_diag_zero <- function(lst) {
  lapply(lst, function(mat) {
    diag(mat) <- 0
    return(mat)
  })
}

baseline_z  <- set_diag_zero(baseline_z)
hdwgcna_z   <- set_diag_zero(hdwgcna_z)
spacex_z    <- set_diag_zero(spacex_z)
smoothie_z  <- set_diag_zero(smoothie_z)

keep_vars <- c("baseline_z", "hdwgcna_z", "spacex_z", "smoothie_z")
rm(list = setdiff(ls(), keep_vars))


################
## begin to draw heatmaps
################

library(pheatmap)
library(gridExtra)
library(grid)

setwd("/home/koe/BSNMani_application-main/Tong_results/BSNMani_for_transfer/BSNMani_for_transfer/BSNMani/4heatmap_comparison/Visualization")

for (i in 1:27) {
  print(i)
  patient_name <- names(baseline_z)[i]
  print(paste0("we are dealing with: ", patient_name))
  
  baseline_temp <- as.matrix(baseline_z[[i]])
  hdwgcna_temp <- as.matrix(hdwgcna_z[[i]])
  spacex_temp <- as.matrix(spacex_z[[i]])
  smoothie_temp <- as.matrix(smoothie_z[[i]])
  
  all_vals <- c(baseline_temp, hdwgcna_temp, spacex_temp, smoothie_temp)
  all_vals <- all_vals[is.finite(all_vals)]
  val_range <- max(abs(quantile(all_vals, probs = c(0.01, 0.99))))
  val_range <- ceiling(val_range * 10) / 10  
  
  breaks <- seq(-val_range, val_range, length.out = 100)
  colors <- colorRampPalette(c("blue", "white", "red"))(99)
  
  p1_baseline <- pheatmap(baseline_temp,
                          color = colors,
                          breaks = breaks,
                          main = "Baseline",
                          cluster_rows = TRUE,
                          cluster_cols = TRUE,
                          fontsize_row = 1.5,
                          fontsize_col = 1.5,
                          treeheight_row = 0,
                          treeheight_col = 0,
                          border_color = NA,
                          silent = TRUE)
  
  row_order <- p1_baseline$tree_row$order
  col_order <- p1_baseline$tree_col$order
  
  p2_hdwgcna <- pheatmap(hdwgcna_temp[row_order, col_order],
                         color = colors,
                         breaks = breaks,
                         main = "hdWGCNA",
                         cluster_rows = FALSE,
                         cluster_cols = FALSE,
                         fontsize_row = 1.5,
                         fontsize_col = 1.5,
                         border_color = NA,
                         silent = TRUE)
  
  p3_spacex <- pheatmap(spacex_temp[row_order, col_order],
                        color = colors,
                        breaks = breaks,
                        main = "SpaceX",
                        cluster_rows = FALSE,
                        cluster_cols = FALSE,
                        fontsize_row = 1.5,
                        fontsize_col = 1.5,
                        border_color = NA,
                        silent = TRUE)
  
  p4_smoothie <- pheatmap(smoothie_temp[row_order, col_order],
                          color = colors,
                          breaks = breaks,
                          main = "smoothie",
                          cluster_rows = FALSE,
                          cluster_cols = FALSE,
                          fontsize_row = 1.5,
                          fontsize_col = 1.5,
                          border_color = NA,
                          silent = TRUE)
  
  png(paste0("compare_4_methods_", patient_name, ".png"), width = 5000, height = 3000, res = 600)
  grid.arrange(p1_baseline$gtable,
               p2_hdwgcna$gtable,
               p3_spacex$gtable,
               p4_smoothie$gtable,
               ncol = 2,
               top = textGrob(paste0("Patient ", patient_name,"'s comparison heatmap"), gp = gpar(fontsize = 15, fontface = "bold")))
  dev.off()
}


