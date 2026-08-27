################
# Trying other methods to group genes
# Tong Liu 2/6/2026
# Try two methods of grouping genes: Hierarchy Clustering and Gaussian Mixture Model
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
U_est = readRDS(fs::path("/nfs/turbo/umms-lgarmire/liutong/BSNMani/SpaceX_BSNMani/BSNMani_output_result",cell_type,paste0("q_",q_val),"MH","diagnostics",
                paste0("U_est_",q_val),ext = "RDS"))

####################
## Method 1 Hierarchy clustering (abandoned)
####################
clustering_hclust = function(X, i)
{
    dist_matrix = as.dist(1 - abs(X))  # Generating distance matrix 
    result_hcluster = hclust(dist_matrix, method = "ward.D2") # Generating hc

    # Plot
    png(fs::path(wd_now,"Visualization","Hierarchy_clustering",paste("HC Dendrogram","i",i,sep="_"),ext="png"),width=12,height=12,units="in",res=300)
    plot(result_hcluster, hang = -0.5, cex = 0.4, main = "Hierarchical Clustering Dendrogram", sub = "", xlab = "", ylab = "Height")
    dev.off()  

    #Using Silhouette Coefficiency to find best group number
    sil_scores = numeric()
    for (k in 2:10) {  
        clusters <- cutree(result_hcluster, k)
        sil_scores[k] <- mean(silhouette(clusters, dist_matrix)[, 3])  
        }
    best_k = which.max(sil_scores)  
    cat("best_k is", best_k)

    png(fs::path(wd_now,"Visualization","Hierarchy_clustering",paste("Best_K","i",i,sep="_"),ext="png"),width=12,height=12,units="in",res=300)
    plot(2:10, sil_scores[2:10], type="b", pch=19, main="Silhouette Score for K Selection")
    abline(v=best_k, col="red", lty=2)

    # Grouping
    clusters = cutree(result_hcluster, k = best_k)
    final_Grouping_result = data.frame(Gene=rownames(X), Cluster=clusters)
    HC_classified_genes = split(final_Grouping_result$Gene, final_Grouping_result$Cluster)
    names(HC_classified_genes) = paste0("Group ", names(HC_classified_genes))
    assign(paste0("HC_classified_genes_", i), HC_classified_genes, envir = .GlobalEnv)


    png(fs::path(wd_now,"Visualization","Hierarchy_clustering",paste("HC_Dendrogram_optimized","i",i,sep="_"),ext="png"),width=12,height=12,units="in",res=300)
    fviz_dend(result_hcluster, 
          k = best_k, 
          rect = FALSE, 
          show_labels = TRUE, 
          cex = 0.4, 
          label_cols = "black",
          lwd = 0.8
          ) + 
        ggtitle("Hierarchical Clustering of Genes") +  
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))   
    dev.off()  

    png(fs::path(wd_now,"Visualization","Hierarchy_clustering",paste("Heatmap","i",i,sep="_"),ext="png"),
                width=12,height=12,units="in",res=300)
    idx_ls_whole = HC_classified_genes
    col.side = c("#845EC2","#D65DB1","#FF6F91","#FF9671","#FFC75F","#F9F871",  
                "#2C73D2","#008F7A","#D5CABD","#AF5C00","#00C9A7","#CA4362","#4E2B00","#FF8066")
    col.side = col.side[1:length(idx_ls_whole)]
    col.heat = colorRampPalette(c("blue", "white", "red"))(256)
    subnetwork_heatmap_v2(idx_ref=idx_ls_whole, X, col.side=col.side, col.heat=col.heat, title=paste("q =",i),
                             na.color="grey", breaks_vec = NULL, lhei_ratio=c(3,5), lwid_ratio=c(3,5))
}

##################################
## Method 2: Gaussian Mixture Model - GMM
##################################

Clustering_GMM = function(X, i){

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

    # Save grouping information
    assign(paste0("GMM_classified_genes_", i), clusters, envir = .GlobalEnv)

    # Visualization
    umap_data = data.frame(
        X = umap_result$layout[,1], 
        Y = umap_result$layout[,2], 
        Cluster = as.factor(clusters)
    )

    png(fs::path(wd_now,"Visualization","GMM_clustering",paste("UMAP","i",i,sep="_"),ext="png"),
        width=12,height=12,units="in",res=300)
    print(
        ggplot(umap_data, aes(x = X, y = Y, color = Cluster)) +
        geom_point(alpha = 0.8, size = 3) +
        theme_minimal() +
        labs(title = "UMAP Visualization of GMM Clustering") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        )
    dev.off()

    # preparation for heatmap
    final_Grouping_result = data.frame(Gene = names(clusters),  
                            Cluster = as.numeric(clusters),  
                            stringsAsFactors = FALSE)
    GMM_classified_genes = split(final_Grouping_result$Gene, final_Grouping_result$Cluster)
    names(GMM_classified_genes) = paste0("Group ", names(GMM_classified_genes))
    assign(paste0("final_Grouping_result_", i), final_Grouping_result, envir = .GlobalEnv)
    assign(paste0("GMM_classified_genes_", i), GMM_classified_genes, envir = .GlobalEnv)
}

##################################
## Method 2 updated version: Gaussian Mixture Model - GMM update
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
# Save 3 subnetwork dataframe in advance: U_sub_1/U_sub_2/U_sub_3
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

setwd("/nfs/turbo/umms-lgarmire/liutong/BSNMani/SpaceX_BSNMani/building/MerFISH_methods_of_grouping_genes_TL/Visualization/GMM_clustering")
dir.create(file.path(cell_type, paste("q_val", q_val, sep = "_")), recursive = TRUE)
setwd(fs::path("/nfs/turbo/umms-lgarmire/liutong/BSNMani/SpaceX_BSNMani/building/MerFISH_methods_of_grouping_genes_TL/Visualization/GMM_clustering",cell_type,paste("q_val",q_val,sep = "_")))

# SNF integration

name_vec = paste0("U_sub_", 1:q_val)   
U_list = mget(name_vec)         
fused_matrix = SNF(U_list, K = 20)

#fused_matrix = SNF(list(U_sub_1, U_sub_2), K = 20)
Clustering_GMM_update(fused_matrix)
GMM_classified_genes_update = GMM_classified_genes
saveRDS(GMM_classified_genes_update, file = "GMM_classified_genes_update.RDS")
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
##############################################



