# Select and plot important (hub) genes from BSNMani subnetworks
#
# Required objects before running this script:
#   u_mean : matrix of gene/marker loadings (genes x subnetworks)
#   output_dir : directory where results will be saved

library(dplyr)
library(ggplot2)
library(patchwork)

# Number of important genes to show for each subnetwork
top_n <- 10

# -----------------------------------------------------------------------------
# 1. Select top hub genes
# -----------------------------------------------------------------------------
# For subnetwork k, the adjacency matrix is u_k %*% t(u_k).
# The hub score is the sum of absolute off-diagonal connection weights.
# This can be calculated directly from u_k without constructing the full matrix:
#   HubScore_i = |u_i| * (sum_j |u_j| - |u_i|)

get_top_hub_genes <- function(u_mean, subnetwork, top_n = 10) {
  u <- u_mean[, subnetwork]
  abs_u <- abs(u)

  hub_score <- abs_u * (sum(abs_u) - abs_u)

  tibble(
    Subnetwork = subnetwork,
    Gene = rownames(u_mean),
    HubScore = hub_score
  ) %>%
    arrange(desc(HubScore)) %>%
    slice_head(n = top_n)
}

# Select top genes from every subnetwork
top_hub_genes <- bind_rows(
  lapply(seq_len(ncol(u_mean)), function(i) {
    get_top_hub_genes(u_mean, i, top_n)
  })
)

print(top_hub_genes)

# -----------------------------------------------------------------------------
# 2. Plot top hub genes
# -----------------------------------------------------------------------------
hub_plots <- lapply(seq_len(ncol(u_mean)), function(i) {
  plot_df <- top_hub_genes %>%
    filter(Subnetwork == i)

  ggplot(plot_df, aes(x = reorder(Gene, HubScore), y = HubScore, fill = HubScore)) +
    geom_col(color = "black", alpha = 0.8) +
    coord_flip() +
    scale_fill_gradient(low = "#377EB8", high = "#E41A1C") +
    labs(
      title = paste0("Top ", top_n, " Hub Markers (Subnetwork ", i, ")"),
      x = "Marker / Gene",
      y = "Total Connectivity (Hub Score)"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "none"
    )
})

combined_plot <- wrap_plots(hub_plots)
print(combined_plot)

# -----------------------------------------------------------------------------
# 3. Save results
# -----------------------------------------------------------------------------
fig_output_dir <- file.path(output_dir, "fig")
dir.create(fig_output_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  top_hub_genes,
  file.path(fig_output_dir, "Top_Hub_Markers.csv"),
  row.names = FALSE
)

ggsave(
  filename = file.path(fig_output_dir, "Top_Hub_Markers.png"),
  plot = combined_plot,
  width = 8 * ncol(u_mean),
  height = 7,
  dpi = 300
)
