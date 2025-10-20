##############################
## This script is used for lambda exploration
##############################

rm(list = ls())
library(readxl)
setwd("D:\\Research_project\\For Test\\about_lambda_itself\\data\\q_5")
lambda_file = readRDS("MH_diag_res_GQN_KPN_q_5.RDS")
lambda = lambda_file$posterior_mean$Lambda_flat
lambda_list = list()

# 5 means we have 5 subnetwork(q_Val = 5)
for (i in 1:5) {
  lambda_list[[i]] = lambda[seq(i, length(lambda), by = 5)]
}

lambda_df = as.data.frame(lambda_list)
colnames(lambda_df) = paste0("lambda_", 1:5)

## normalization for lambdas
log_lambda_df = log(lambda_df)

## find patient id vector
patient_id = names(readRDS("D:\\Research_project\\For Test\\about_lambda_itself\\data\\q_5\\Smoothie_co_expression_list_Oli.RDS"))
meta = read_excel("D:\\Research_project\\For Test\\about_lambda_itself\\data\\q_5\\meta.xlsx")
data_final = cbind(meta[,c(1,2,4)],log_lambda_df)


#############################
## draw - with clustering of patient
############################
library(pheatmap)
library(dplyr)

ordered_vars = c("Atherosclerosis", "lambda_4","lambda_3", "lambda_5", "lambda_2","lambda_1")

heatmap_mat = as.data.frame(data_final[, ordered_vars])
heatmap_scaled = as.data.frame(sapply(heatmap_mat, function(col) (col - min(col)) / (max(col) - min(col))))

rownames(heatmap_scaled) = data_final$'Donor ID'


annotation_row = data.frame(
  Dementia = factor(
    ifelse(data_final$Cognitive_Status == 1, "Dementia", "No Dementia"),
    levels = c("No Dementia", "Dementia")
  )
)
rownames(annotation_row) = data_final$'Donor ID'


ann_colors = list(
  Dementia = c("No Dementia" = "#A6CEE3",  
               "Dementia" = "#FB9A99")  
)

attr(heatmap_scaled, "name") = "Value"

pheatmap(
  mat = heatmap_scaled,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  scale = "none",
  fontsize = 12,
  fontsize_row = 6,
  fontsize_col = 10,
  border_color = NA,
  color = colorRampPalette(c("#4575b4", "#e0f3f8", "#ffffbf", "#fee090", "#d73027"))(100)
)


#############################
## re-draw - WITHOUT clustering patients
############################
library(pheatmap)
library(dplyr)

data_final_without = data_final[order(data_final$Cognitive_Status),]
data_final_without = data_final_without[,c("Donor ID","Cognitive_Status")]

heatmap_scaled_group = heatmap_scaled[data_final_without$`Donor ID`,]

annotation_row <- data.frame(
  Dementia = factor(
    ifelse(data_final_without$Cognitive_Status == 1, "Dementia", "No Dementia"),
    levels = c("No Dementia", "Dementia")
  )
)
rownames(annotation_row) = data_final_without$'Donor ID'


ann_colors = list(
  Dementia = c("No Dementia" = "#A6CEE3",  
               "Dementia" = "#FB9A99") 
)

attr(heatmap_scaled, "name") = "Value"

pheatmap(
  mat = heatmap_scaled_group,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  scale = "none",
  fontsize = 12, 
  fontsize_row = 6,
  fontsize_col = 10,
  border_color = NA,
  color = colorRampPalette(c("#4575b4", "#e0f3f8", "#ffffbf", "#fee090", "#d73027"))(100)
)


