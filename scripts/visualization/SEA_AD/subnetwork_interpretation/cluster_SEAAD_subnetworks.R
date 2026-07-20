#############################
library(cluster)   
library(gplots)    
library(factoextra) 
library(ggplot2)
library(mclust)
library(umap)
library(dplyr)
library(clue)
library(SNFtool) # For SNF clustering algorithm
############################
source("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SEA-AD/BSNMani-dev/g1_diagnostics_helper.R")
###########################

wd_now = getwd()     

Genes_name_vec = readRDS("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/SEA_AD/preprocessed/common_genes.RDS")


print(paste0("q value in grouping genes is ",q_val))

#U_est = readRDS("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/SEA_AD/preprocessed/U_est.RDS")
U_est = readRDS(fs::path("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SpaceX_BSNMani/result",cell_type,paste0("q_",q_val),"MH","diagnostics",
                paste0("U_est_",q_val),ext = "RDS"))

##################################
## Method updated version: Gaussian Mixture Model - GMM update
##################################

Clustering_GMM_update = function(X){

    # Clustering using GMM
    gmm_model = Mclust(X)
    summary(gmm_model)

    # Filter the grouping information with high confidence result
    threshold = 0.95
    probs = gmm_model$z
    high_confidence_genes = which(apply(probs, 1, max) > threshold)

    # Gaining filtered matrix
    X = X[high_confidence_genes,high_confidence_genes]

    # Using UMAP to lower the dimension
    umap_result = umap(X)
    clusters = gmm_model$classification[high_confidence_genes]  
    print(clusters)

    # preparation for heatmap
    final_Grouping_result = data.frame(Gene = names(clusters),  
                            Cluster = as.numeric(clusters),  
                            stringsAsFactors = FALSE)
    GMM_classified_genes = split(final_Grouping_result$Gene, final_Grouping_result$Cluster)
    assign("GMM_classified_genes", GMM_classified_genes, envir = .GlobalEnv)
    
}


##############################################
## Similarity Network Fusion (SNF) algorithm
##############################################
cat("Similarity Network Fusion (SNF) algorithm")
## Function preparation
# Save q_val subnetwork dataframe in advance: U_sub_1/U_sub_2/U_sub_3/...
generate_subnetwork = function(X, i)
{
    U_sub = U_est[,i,drop=FALSE]%*%t(U_est[,i,drop=FALSE])
    colnames(U_sub) = Genes_name_vec
    rownames(U_sub) = Genes_name_vec
    assign(paste0("U_sub_", i), U_sub, envir = .GlobalEnv)

}
## data preparation
for (i in 1:q_val)
{
   generate_subnetwork(X, i)
}

setwd("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SpaceX_BSNMani/building/MerFISH_methods_of_grouping_genes_TL/Visualization/GMM_clustering")
dir.create(file.path(cell_type, paste("q_val", q_val, sep = "_")), recursive = TRUE)
setwd(fs::path("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SpaceX_BSNMani/building/MerFISH_methods_of_grouping_genes_TL/Visualization/GMM_clustering",cell_type,paste("q_val",q_val,sep = "_")))

# SNF integration
name_vec = paste0("U_sub_", 1:q_val)   
U_list = mget(name_vec)         
fused_matrix = SNF(U_list, K = 20)

#fused_matrix = SNF(list(U_sub_1, U_sub_2), K = 20)
Clustering_GMM_update(fused_matrix)
GMM_classified_genes_update = GMM_classified_genes
# Heatmap
for(i in 1:q_val){
    U_sub = U_est[,i,drop=FALSE]%*%t(U_est[,i,drop=FALSE])
    colnames(U_sub) = Genes_name_vec
    rownames(U_sub) = Genes_name_vec


    png(filename = paste0("final_Heatmap_update_", i, ".png"),
                width=12,height=12,units="in",res=300)

    idx_ls_whole = GMM_classified_genes_update
    col.side = c("#845EC2","#D65DB1","#FF6F91","#FF9671","#FFC75F","#F9F871",  
                "#2C73D2","#008F7A","#D5CABD","#AF5C00","#00C9A7","#CA4362","#4E2B00","#FF8066")
    col.side = col.side[1:length(idx_ls_whole)]
    col.heat = colorRampPalette(c("blue", "white", "red"))(256)
    subnetwork_heatmap_v2(idx_ref=idx_ls_whole, U_sub, col.side=col.side, col.heat=col.heat, title=paste("q =",i),
                             na.color="grey", breaks_vec = NULL, lhei_ratio=c(3,5), lwid_ratio=c(3,5))
}


