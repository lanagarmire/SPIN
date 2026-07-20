#' Plot heatmaps of smoothie co-expression matrix for each patients

source("/nfs/dcmb-lgarmire/yangiwen/workspace/common/Utils.R")

# Read raw data for BC
marker_data <- read.csv("/nfs/dcmb-lgarmire/yangiwen/workspace/bsnmani/data/bc/Filtered_and_Renamed_Markers.csv")
spatial <- read.csv("/nfs/dcmb-lgarmire/yangiwen/workspace/bsnmani/data/bc/Basel_SC_locations.csv")
clusters <- read.csv("/nfs/dcmb-lgarmire/yangiwen/workspace/bsnmani/data/bc/Cluster_labels/Basel_metaclusters.csv")
cluster_annotation <- read.csv("/nfs/dcmb-lgarmire/yangiwen/workspace/bsnmani/data/bc/Cluster_labels/Metacluster_annotations.csv")
clinical_df_coxph <- read.csv("/nfs/dcmb-lgarmire/yangiwen/workspace/bsnmani/output/bc/output_05/q5/s123/clinical_df_coxph.csv")

# Assign annotation to the cells
cluster_annotation <- data.frame(do.call(rbind, strsplit(cluster_annotation$Metacluster..Cell.type.Class, ";")))
clusters$annotation <- cluster_annotation[match(clusters$cluster, cluster_annotation$X1), "X3"]
spatial$annotation <- clusters[match(spatial$id, clusters$id), "annotation"]

# Create output dir
output_dir <- "/nfs/dcmb-lgarmire/yangiwen/workspace/bsnmani/output/bc/coexpr"
mkdir(output_dir)

# Read processed coexpression matrix from Haowen
coexpr_smoothie_all <- readRDS("/home/yangiwen/n/workspace/bsnmani/r/R/extern/haowen/BSNMANi/BC/WGCNA_transformed_smoothed_padded44.RDS")
# coexpr_wgcna_all <- readRDS("/home/yangiwen/n/workspace/bsnmani/r/R/extern/haowen/BSNMANi/BC/WGCNA_transformed_BC_core.RDS")
# coexpr_hdwgcna_all <- readRDS("/home/yangiwen/n/workspace/bsnmani/r/R/extern/haowen/BSNMANi/BC/hdwgcna_44matched_list.RDS")
# coexpr_spacex_all <- readRDS("/home/yangiwen/n/workspace/bsnmani/r/R/extern/haowen/BSNMANi/BC/spacex_44matched_list.RDS")

applyGeneName <- function(mtx) {
  mtx <- mtx[marker_data$Original.Metal, marker_data$Original.Metal]
  rownames(mtx) <- marker_data$Gene
  colnames(mtx) <- marker_data$Gene
  mtx
}

coexpr_smoothie_all <- lapply(coexpr_smoothie_all, applyGeneName)
# coexpr_wgcna_all <- lapply(coexpr_wgcna_all, applyGeneName)
# coexpr_hdwgcna_all <- lapply(coexpr_hdwgcna_all, applyGeneName)

plotCoexpression <- function(
  core,
  plot_pearson = T,
  plot_smoothie = F,
  plot_spatial = T
) {
  coexpr_list <- list()
  if (plot_pearson) {
    coexpr_pearson <- read.csv(file.path("/nfs/dcmb-lgarmire/yangiwen/workspace/bsnmani/r/R/extern/haowen/BSNMANi/BC/correlation_matrix", paste0("corr_", core, ".csv")), row.names = 1)
    coexpr_pearson <- applyGeneName(as.matrix(coexpr_pearson))
    coexpr_list$pearson <- coexpr_pearson
  }
  if (plot_smoothie) {
    coexpr_list$smoothie <- coexpr_smoothie_all[[core]]
  }
  n_plots <- length(coexpr_list)
  plots <- if (n_plots > 0) {
    order <- hclust(dist(coexpr_list[[1]]))$order
    lapply(names(coexpr_list), function(name) {
      coexpr <- coexpr_list[[name]]
      pheatmap::pheatmap(coexpr[order, order], cluster_rows = F, cluster_cols = F, main = name)$gtable
    })
  } else list()
  spatial_plot <- if (plot_spatial) {
    color_list <- c(
      "Immune" = "royalblue",
      "Stroma" = "seagreen",
      "Tumor" = "lightcoral",
      "Vessel" = "purple"
    )
    spatial_df <- spatial[spatial$core == core,]
    patient_status <- clinical_df_coxph[clinical_df_coxph$core == core, "Patientstatus"]
    predicted_group <- clinical_df_coxph[clinical_df_coxph$core == core, "group"]
    osmonth <- clinical_df_coxph[clinical_df_coxph$core == core, "OSmonth"]
    age <- clinical_df_coxph[clinical_df_coxph$core == core, "age"]
    grade <- clinical_df_coxph[clinical_df_coxph$core == core, "grade"]
    clinical_type <- clinical_df_coxph[clinical_df_coxph$core == core, "clinical_type"]
    title <- paste0(
      "Patient ",
      "<span style='color:",
      if (patient_status == 0) "blue" else "red",
      "'>",
      if (patient_status == 0) "alive" else "dead",
      "</span>, predicted <span style='color:",
      if (predicted_group == "low") "blue" else "red",
      "'>",
      predicted_group,
      " risk</span>, ",
      osmonth,
      " OS months ",
      age,
      " grade ",
      grade,
      " ",
      clinical_type
    )
    ggplot2::ggplot(spatial_df, ggplot2::aes(Location_Center_X, Location_Center_Y, color = annotation)) +
      ggplot2::geom_point(size = 3) +
      ggplot2::scale_color_manual(values = color_list) +
      ggplot2::theme_void() +
      ggplot2::ggtitle(title) +
      ggplot2::theme(plot.title = ggtext::element_markdown(hjust = 0.5, face = "bold"))
  } else NULL
  if (!is.null(spatial_plot)) {
    n_plots <- n_plots + 1
    plots <- c(list(spatial_plot), plots)
  }
  if (n_plots > 0) {
    patchwork::wrap_plots(plots, ncol = n_plots)
  } else NULL
}

for (core in names(coexpr_smoothie_all)) {
  png(file.path(output_dir, paste0(core, ".png")), width = 3072, height = 1080, res = 120)
  print(plotCoexpression(core, plot_smoothie = T))
  dev.off()
}
