###################
## try to open raw csv data from Tay
###################
rm(list = ls())
library(bit64)
library(data.table)
library(data.table)
library(openxlsx)  
#options(scipen = 999)
setwd("D:\\For Test\\For Cell Annotation")
tay = fread("SEAAD_MTG_MERFISH.2024-12-11_metadata.csv")
tay$`Cell ID` = as.character(tay$`Cell ID`)

cell_type = tay[,c("Donor ID","Subclass","Specimen Barcode","Cell ID")]
cell.type_filter = unique(cell_type)
length(unique(cell.type_filter$`Cell ID`))   ## Has 67965 unique cells

dornor_ID = unlist(unique(cell.type_filter$`Donor ID`))
dornor_spe_info = list()
for (i in 1:27)
{
  print(i)
  print(dornor_ID[i])
  dornor_spe_info[[i]] = cell.type_filter[cell.type_filter$`Donor ID` == dornor_ID[i], ]
  minbarcode = min(dornor_spe_info[[i]]$`Specimen Barcode`)
  if (dornor_ID[i] == "H20.33.012") { minbarcode = 1218552479 } ## DONE
  dornor_spe_info[[i]] = dornor_spe_info[[i]][dornor_spe_info[[i]]$`Specimen Barcode` == minbarcode, ]
  dornor_spe_info[[i]] = dornor_spe_info[[i]][,c("Subclass","Cell ID")]
  names(dornor_spe_info)[i] = dornor_ID[i]
}

saveRDS(dornor_spe_info,file = "dornor_spe_info.RDS")



###################
## load MerFish data
##################
rm(list = ls())
dornor_spe_info = readRDS("dornor_spe_info.RDS")

##################
## begin to match
##################
library(dplyr)
load("for_annotation.RData")     ## include norm_scaled_exprs_merFISH and spatial_locs_merFISH, 
                                 ## and label each corresponding sub-dataframe with the respective patient ID.

names_ref = names(dornor_spe_info)
names_now = names(spatial_locs_merFISH)

## find common dornor name and filter
common_names = intersect(names_ref, names_now)

dornor_spe_info = dornor_spe_info[common_names]
spatial_locs_merFISH = spatial_locs_merFISH[common_names]
count_list = count_list[common_names]

## start matching
patient_ID = unlist(names(dornor_spe_info))

express_matrix = list()
cell_loc = list()

# process 26 patients' cell annotation
for (i in 1:26)
{
  print(i)
  spatial_temp = spatial_locs_merFISH[[patient_ID[i]]]
  reference_temp = dornor_spe_info[[patient_ID[i]]]
  reference_temp = reference_temp[!duplicated(reference_temp$`Cell ID`), ]
  
  
  spatial_temp$cell_ID = as.character(spatial_temp$cell_ID)
  reference_temp$`Cell ID` = as.character(reference_temp$`Cell ID`)
  
  spatial_temp = left_join(
    spatial_temp,
    reference_temp,
    by = c("cell_ID" = "Cell ID")  
  )
  spatial_temp = spatial_temp[!is.na(spatial_temp$Subclass), ]
  
  express_temp = count_list[[patient_ID[i]]]
  express_temp = t(express_temp)
  express_temp = as.data.frame(express_temp)
  
  
  spatial_temp$cell_ID = as.character(spatial_temp$cell_ID)
  rownames(express_temp) = as.character(rownames(express_temp))
  
  express_temp = express_temp[match(spatial_temp$cell_ID, rownames(express_temp)), ]
  
  rownames(express_temp) = seq_len(nrow(express_temp))
  rownames(spatial_temp) = seq_len(nrow(spatial_temp))
  spatial_temp = spatial_temp[,-3]
  
  ## save in list
  express_matrix[[i]] = express_temp
  cell_loc[[i]] = spatial_temp
  names(cell_loc)[i] = patient_ID[i]
  names(express_matrix)[i] = patient_ID[i]
}

save(cell_loc, express_matrix, file = "data_for_SPACEX.RData")



###################
## Try to summarize cell type
###################
rm(list = ls())
setwd("D:\\For Test\\For Cell Annotation")
load("data_for_SPACEX.RData")

class_counts_list = lapply(cell_loc, function(dt) {
  if (is.data.table(dt)) {
    dt[, .N, by = Subclass]
  } else if (is.list(dt) && is.data.table(dt[[1]]) && "Subclass" %in% names(dt[[1]])) {
    dt[[1]][, .N, by = Subclass]
  } else {
    NULL
  }
})


merged_list_temp = list()

for (sample_name in names(class_counts_list)) {
  dt <- class_counts_list[[sample_name]]
  setnames(dt, c("Subclass", sample_name))  
  merged_list_temp[[sample_name]] <- dt
}

merged_dt = Reduce(function(x, y) merge(x, y, by = "Subclass", all = TRUE), merged_list_temp)

merged_dt[is.na(merged_dt)] <- 0
merged_dt = t(merged_dt)
colnames(merged_dt) = merged_dt[1,]
merged_dt = merged_dt[-1,]
merged_dt_filter = merged_dt[,c("Chandelier","L5 ET","L5/6 NP","Pax6","Sncg","Sst Chodl")]

write.xlsx(merged_dt, file = "Subclass_counts_summary.xlsx", rowNames = TRUE)
write.xlsx(merged_dt_filter, file = "Subclass_counts_summary_filter.xlsx", rowNames = TRUE)



