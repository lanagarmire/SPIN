
########################
## generate input for Smoothie (expression csv file and spatial csv file for each patient)
########################
rm(list = ls())
#setwd("D:\\Research_project\\For Test\\Cell_type_specific_BSNMani\\data")
#load("final_data_for_SPACEX.RData")

setwd("D:\\Research_project\\For Test\\H20.33.001\\coexpression-generation\\data_uni")
load("data_for_uni_27.RData")
count_list = express_matrix_temp[["H20.33.001"]]  # split out H20.33.001 for smoothie construction
loc_list = cell_loc_temp[["H20.33.001"]]

#setwd("D:\\Research_project\\For Test\\smoothie\\smoothie_input\\celltype_specific\\csv")
setwd("D:\\Research_project\\For Test\\smoothie\\smoothie_input\\csv")
#for (i in 1:length(count_list)) {
#  sample_id <- names(count_list)[i]
#  counts_df <- as.matrix(count_list[[i]])
#  rownames(counts_df) = seq(1,nrow(counts_df))
#  write.csv(counts_df, file = paste0("sample_", sample_id, "_counts.csv"), row.names = TRUE)
#  write.csv(loc_list[[i]][, c("sdimx", "sdimy")], file = paste0("sample_", sample_id, "_spatial.csv"), row.names = FALSE)
#}   # since we only use H20.33.001, then list -> dataframe

# only for H20.33.001
sample_id = "H20.33.001"
counts_df = as.matrix(count_list)
rownames(counts_df) = seq(1,nrow(counts_df))
write.csv(counts_df, file = paste0("sample_", sample_id, "_counts.csv"), row.names = TRUE)
write.csv(loc_list[, c("sdimx", "sdimy")], file = paste0("sample_", sample_id, "_spatial.csv"), row.names = FALSE)

count_temp = read.csv("sample_H20.33.001_counts.csv")
loc_temp = read.csv("sample_H20.33.001_spatial.csv")


names(loc_list)

setwd("D:\\Research_project\\For Test\\smoothie\\smoothie_input")
count_temp_1 = read.csv("sample_H21.33.028_counts.csv")
loc_temp_1 = read.csv("sample_H21.33.028_spatial.csv")

smooth_temp = read.csv("D:\\Research_project\\For Test\\smoothie\\smoothie_output\\celltype_specific\\csv\\smoothed_expression_matrix.csv")
#########################
## generate data for BSNMani input
#########################
rm(list = ls())
folder_path <- "D:\\Research_project\\For Test\\smoothie\\smoothie_output\\csv"
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)
patient_ids <- gsub("_pearsonR\\.csv$", "", basename(csv_files))
csv_list <- lapply(csv_files, read.csv)
names(csv_list) <- patient_ids
gene_name = readRDS("D:\\Research_project\\For Test\\smoothie\\common_genes.RDS")

for (i in seq_along(csv_list)) {
  rownames(csv_list[[i]]) <- gene_name
  colnames(csv_list[[i]]) <- gene_name
}
setwd("D:\\Research_project\\For Test\\smoothie")
saveRDS(csv_list,file = "Smoothie_co_expression_list.RDS")

rm(list = ls())
setwd("D:\\Research_project\\For Test\\smoothie")
temp = readRDS("Smoothie_co_expression_list_Astrocyte_update.RDS")

############################################
# generate BSNMani input for H20.33.001
setwd("D:\\Research_project\\For Test\\smoothie\\smoothie_output\\csv")
expr = read.csv("H20.33.001_pearsonR.csv")
gene_name = readRDS("D:\\Research_project\\For Test\\smoothie\\common_genes.RDS")
rownames(expr) = gene_name
colnames(expr) = gene_name

former_expr = readRDS("D:\\Research_project\\For Test\\smoothie\\Smoothie_co_expression_list.RDS")
meta_uni = readxl::read_excel("C:\\Users\\liutong\\Desktop\\SEA-AD Updated Data\\meta_uni.xlsx")

new_expr = former_expr[names(former_expr) != "H20.33.002"]
new_expr[["H20.33.001"]] = expr

new_expr_ordered <- new_expr[meta_uni$`Donor ID`]
setwd("D:\\Research_project\\For Test\\smoothie")
saveRDS(new_expr_ordered, file = "Smoothie_co_expression_list_uni.RDS")
