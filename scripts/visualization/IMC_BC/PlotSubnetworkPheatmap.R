#' Plot heatmaps of q subnetworks

source("/nfs/dcmb-lgarmire/yangiwen/workspace/common/Utils.R")
# Source functions to calculate U
Rcpp::sourceCpp("/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SEA-AD/BSNMani-dev/hybrid_M0_MALA_LR_FAST_v2.cpp")

output_dir <- file.path("/nfs/dcmb-lgarmire/yangiwen/workspace/bsnmani/output/bc/selected_model")
fig_output_dir <- file.path(output_dir, "fig")
mkdir(fig_output_dir)

q <- 2
n_roi <- 44
# All 44 metal tags
metals_all <- c("Ar80", "Dy162", "Dy163", "Dy164", "Er166", "Er167", "Er168", "Eu151", "Gd155", "Gd156", "Gd158", "Gd160", "Hg202", "In113", "In115", "Ir191", "Ir193", "La139", "Lu175", "Nd142", "Nd143", "Nd144", "Nd145", "Nd146", "Nd148", "Nd150", "Pb204", "Pb206", "Pr141", "Ru100", "Ru101", "Ru102", "Ru104", "Ru96", "Ru98", "Ru99", "Sm147", "Sm149", "Sm152", "Tm169", "Xe134", "Yb172", "Yb174", "Yb176")
# Ordered 29 metal tags
selected_metals <- c("Lu175", "Er167", "Eu151", "Nd143", "Dy163", "Nd144", "Gd158", "Gd156", "In113", "Er166", "La139", "Sm152", "Dy164", "Er168", "Yb172", "Sm147", "Gd155", "Pr141", "Tm169", "Yb174", "Nd145", "Nd150", "Yb176", "Dy162", "Gd160", "Nd148", "Nd146", "Sm149", "Nd142")

# Calculate U from training output
pos_mean <- readRDS(file.path(output_dir, "train", "g1", "diagnostics", "g1_res_GQN_train.RDS"))
x_mean <- matrix(pos_mean$X, n_roi, q)
u_mean <- polar_expansion(x_mean)
# Assign row names of U to be metal tags
rownames(u_mean) <- metals_all
u_mean <- u_mean[selected_metals, , drop = F]
# Map metal tag to target genes
marker_data <- read.csv("/nfs/dcmb-lgarmire/yangiwen/workspace/bsnmani/data/bc/Filtered_and_Renamed_Markers.csv")
rownames(u_mean) <- marker_data[match(rownames(u_mean), marker_data$Original.Metal), "Gene"]

# Plot one heatmap for each U (subnetwork)
plots <- lapply(1:q, function(i) {
  u <- u_mean[, i, drop = F]
  u <- u %*% t(u)
  diag(u) <- NA
  pheatmap::pheatmap(u, cluster_cols = F, cluster_rows = F, main = paste0("Subnetwork ", i), fontsize = 20)$gtable
})

# Save the combined plot
png(file.path(fig_output_dir, "subnetwork.png"), width = 1024 * length(plots), height = 1024, res = 120)
print(patchwork::wrap_plots(plots))
dev.off()
