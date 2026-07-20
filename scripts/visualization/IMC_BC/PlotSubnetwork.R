#' Plot heatmaps of q subnetworks

source("/nfs/dcmb-lgarmire/yangiwen/workspace/common/Utils.R")
Rcpp::sourceCpp("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SEA-AD/BSNMani-dev/hybrid_M0_MALA_LR_FAST_v2.cpp")

subnetwork_heatmap_v2 <- function(idx_ref, U_sub, col.side, col.heat, ignore.diag = TRUE, title, breaks_vec, na.color, legend_position = "left", lhei_ratio = c(2, 5), lwid_ratio = c(2, 5)) {
  modules <- names(idx_ref)
  modules_sz <- sapply(idx_ref, FUN = function(x) { length(x) })
  idx_vec <- unname(unlist(idx_ref))
  u_arrange <- U_sub[idx_vec, idx_vec]
  borders <- cumsum(modules_sz)
  borders <- borders[-length(borders)]
  borders <- borders

  if (ignore.diag) {
    diag(u_arrange) <- NA
  }

  gplots::heatmap.2(
    u_arrange, col = col.heat, scale = "none", dendrogram = 'none', main = title, breaks = breaks_vec, na.color = na.color,
    lwid = lwid_ratio, lhei = lhei_ratio, Rowv = FALSE, Colv = FALSE, trace = 'none', colsep = borders, rowsep = borders, sepcolor = "black",
    RowSideColors = unname(unlist(mapply(FUN = function(x, y) { rep(x, y) }, col.side, modules_sz))),
    ColSideColors = unname(unlist(mapply(FUN = function(x, y) { rep(x, y) }, col.side, modules_sz))),
    cexRow = 2, cexCol = 2, margins = c(10, 10)
  )
  legend(legend_position, legend = modules, col = col.side, lty = 1, lwd = 5, cex = 1, border = NA)
}

q_vals <- 2:8
seeds <- c(0, 42, 64, 123, 894)
n_roi <- 44

# for (q in q_vals) {
#   for (seed in seeds) {
#     output_dir <- file.path("/nfs/dcmb-lgarmire/yangiwen/workspace/bsnmani/output/bc/output_10", paste0("q", q), paste0("s", seed))
q <- 2
output_dir <- file.path("/nfs/dcmb-lgarmire/yangiwen/workspace/bsnmani/output/bc/selected_model")
fig_output_dir <- file.path(output_dir, "fig")
mkdir(fig_output_dir)

pos_mean <- readRDS(file.path(output_dir, "train", "g1", "diagnostics", "g1_res_GQN_train.RDS"))

x_mean <- matrix(pos_mean$X, n_roi, q)
u_mean <- polar_expansion(x_mean)

palette <- colorRampPalette(c("blue", "white", "red"))(256)

marker_data <- read.csv("/nfs/dcmb-lgarmire/yangiwen/workspace/bsnmani/data/bc/Filtered_and_Renamed_Markers.csv")
metal_to_marker <- setNames(marker_data$Gene, marker_data$Original.Metal)
selected_metals <- marker_data$Original.Metal

genes_name_vec_all <- c("Ar80", "Dy162", "Dy163", "Dy164", "Er166", "Er167", "Er168", "Eu151", "Gd155", "Gd156", "Gd158", "Gd160", "Hg202", "In113", "In115", "Ir191", "Ir193", "La139", "Lu175", "Nd142", "Nd143", "Nd144", "Nd145", "Nd146", "Nd148", "Nd150", "Pb204", "Pb206", "Pr141", "Ru100", "Ru101", "Ru102", "Ru104", "Ru96", "Ru98", "Ru99", "Sm147", "Sm149", "Sm152", "Tm169", "Xe134", "Yb172", "Yb174", "Yb176")

rownames(u_mean) <- genes_name_vec_all
u_mean <- u_mean[selected_metals, , drop = F]

u_list <- list()
for (i in 1:q) {
  ui <- u_mean[, i, drop = F]
  u_sub <- ui %*% t(ui)
  rownames(u_sub) <- selected_metals
  colnames(u_sub) <- selected_metals
  u_list[[i]] <- u_sub
}
names(u_list) <- paste0("u_sub_", 1:q)

fused_matrix <- SNFtool::SNF(u_list, K = 20)
rownames(fused_matrix) <- metal_to_marker[rownames(fused_matrix)]
colnames(fused_matrix) <- metal_to_marker[colnames(fused_matrix)]

library(mclust)
gmm_model <- mclust::Mclust(fused_matrix, G = 3:14)
clusters <- gmm_model$classification
final_grouping_result <- data.frame(
  gene = rownames(fused_matrix),
  cluster = clusters
)
gmm_classified_genes <- split(final_grouping_result$gene, final_grouping_result$cluster)
names(gmm_classified_genes) <- paste0("group_", names(gmm_classified_genes))

color <- c("#845EC2", "#D65DB1", "#FF6F91", "#FF9671", "#FFC75F", "#F9F871",
           "#2C73D2", "#008F7A", "#D5CABD", "#AF5C00", "#00C9A7", "#CA4362",
           "#4E2B00", "#FF8066")[seq_along(gmm_classified_genes)]

for (i in 1:q) {
  u_sub_metal <- u_list[[i]]

  rownames(u_sub_metal) <- metal_to_marker[rownames(u_sub_metal)]
  colnames(u_sub_metal) <- metal_to_marker[colnames(u_sub_metal)]

  ordered_genes <- unlist(gmm_classified_genes)
  u_sub_metal <- u_sub_metal[ordered_genes, ordered_genes]

  png(file.path(fig_output_dir, paste0("subnetwork_", i, ".png")),
      width = 1920, height = 1920, res = 120)

  subnetwork_heatmap_v2(
    gmm_classified_genes,
    u_sub_metal,
    color,
    palette,
    paste("SNF-GMM Clustered U", i),
    na.color = "grey",
    breaks_vec = NULL,
    lhei_ratio = c(2, 8),
    lwid_ratio = c(2, 8),
    ignore.diag = TRUE
  )
  dev.off()
}
#   }
# }
