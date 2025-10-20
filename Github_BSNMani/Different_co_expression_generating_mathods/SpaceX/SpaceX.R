
library(MSFA)
library(SpaceX)
library(dplyr)
library(doParallel)
setwd("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SpaceX/data")

#####################################
## For array computing
#####################################
args = commandArgs(trailingOnly = TRUE)
i = as.numeric(args[1])  
cat("Running analysis for i =", i, "\n")



data = readRDS("final_data_for_SPACEX.RDS")
patient_name = names(data)[i]
cat(paste0("we are dealing with_", patient_name))
data_temp = data[[i]]
    

BC_fit = SpaceX(data_temp$count, data_temp$loc[,1:2], data_temp$loc[,3],sPMM=FALSE,Post_process = TRUE,numCore = 16)

cat(paste0("Complete patient_",patient_name))

setwd("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SpaceX/data/SpaceX_res_all_3celltype_update")
saveRDS(BC_fit, file = paste0("SpaceX_res_", i, ".RDS"))

cat("\n")
str(BC_fit)
