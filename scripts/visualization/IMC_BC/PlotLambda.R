source("/nfs/dcmb-lgarmire/yangiwen/workspace/common/Utils.R")

output_dir <- file.path("/nfs/dcmb-lgarmire/yangiwen/workspace/bsnmani/output/bc/selected_model")
q <- 2
clinical_df <- read.csv(file.path(output_dir, "clinical_df_coxph.csv"))
clinical_df$clusters <- clinical_df$group
categorical_cols <- c("grade", "clinical_type", "Patientstatus", "clusters")
for (col in categorical_cols) {
  clinical_df[[col]] <- factor(clinical_df[[col]], levels = sort(unique(clinical_df[[col]])))
}
# clinical_types <- unique(clinical_df$clinical_type)

# for (clinical_type in clinical_types) {
#   df <- clinical_df[clinical_df$clinical_type == clinical_type,]
df <- clinical_df
ann <- df[, c("clusters", "OSmonth", "Patientstatus"), drop = F]
ann$Patientstatus <- ifelse(ann$Patientstatus == "0", "Alive", "Dead")
levels(ann$Patientstatus) <- c("Alive", "Dead")

mtx <- as.matrix(sapply(df[, c("age", "grade", paste0("lambda_", 1:q))], as.numeric))
rownames(mtx) <- rownames(df)
mtx <- apply(mtx, 2, minmax)

# Order rows by cluster and subcluster (hierarchical clustering within each cluster)
order_idx <- unlist(lapply(levels(ann$Patientstatus), function(status) {
  idx <- which(ann$Patientstatus == status)
  sub_hclust <- hclust(dist(mtx[idx,]))
  idx[sub_hclust$order]
}))
order_idx <- rev(order(ann$OSmonth))

ann <- ann[order_idx, , drop = FALSE]
mtx <- mtx[order_idx, , drop = FALSE]

# n_clusters <- length(levels(ann$clusters))
ann_colors <- list(
  Patientstatus = c("Alive" = "lightblue", "Dead" = "lightcoral")
  # clinical_type = setNames(RColorBrewer::brewer.pal(length(levels(ann$clinical_type)), "Set2"), levels(ann$clinical_type)),
  # grade = setNames(RColorBrewer::brewer.pal(length(levels(ann$grade)), "Set3"), levels(ann$grade)),
  # clusters = setNames(RColorBrewer::brewer.pal(n_clusters, "Set1")[1:n_clusters], levels(ann$clusters))
)

# png(file.path(output_dir, "fig", paste0("lambda_cluster_", clinical_type, ".png")), width = 1024, height = 1024, res = 120)
png(file.path(output_dir, "fig", "lambda_cluster.png"), width = 1024, height = 1024, res = 120)
print(pheatmap::pheatmap(
  mtx,
  cluster_cols = F,
  cluster_rows = T,
  annotation_row = ann[, c("OSmonth", "Patientstatus")],
  annotation_colors = ann_colors,
  show_rownames = F,
  fontsize = 18
  # main = paste0("Heatmap for ", clinical_type, " q = ", q)
))
dev.off()
# }
