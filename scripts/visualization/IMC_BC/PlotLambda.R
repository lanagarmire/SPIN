source("/nfs/dcmb-lgarmire/yangiwen/workspace/common/Utils.R")

output_dir <- file.path("/nfs/dcmb-lgarmire/yangiwen/workspace/bsnmani/output/bc/selected_model")
q <- 2
clinical_df <- read.csv(file.path(output_dir, "clinical_df_coxph.csv"))
# 1. Prepare Annotations (Columns/Top)
ann <- data.frame(
  "Patient status" = factor(
    ifelse(clinical_df$Patientstatus %in% c("1", 1, "Dead"), "Dead", "Alive"),
    levels = c("Alive", "Dead")),
  "OS month" = as.numeric(clinical_df$OSmonth),
  check.names = FALSE
)
rownames(ann) <- rownames(clinical_df)

# 2. Prepare Matrix and Transpose
mtx <- as.matrix(sapply(clinical_df[, c("age", "grade", paste0("lambda_", 1:q))], as.numeric))
rownames(mtx) <- rownames(clinical_df)
mtx <- apply(mtx, 2, minmax)
mtx <- t(mtx)   # features as rows, patients as columns

# 3. SORT PATIENTS BY OUTCOME (replaces the dendrogram / cluster_cols = TRUE)
order_idx <- switch(sort_mode,
  # os        = order(ann[["OS month"]]),                                  # short -> long
  os = order(ann[["OS month"]], decreasing = TRUE),   # long -> short
  status_os = order(match(ann[["Patient status"]], c("Dead", "Alive")),  # Dead | Alive
                    ann[["OS month"]])
)
mtx <- mtx[, order_idx, drop = FALSE]
ann <- ann[order_idx, , drop = FALSE]

gaps <- if (sort_mode == "status_os") {
  as.integer(table(ann[["Patient status"]])[["Dead"]])
} else NULL

# 4. Colors (unchanged from your version)
ann_colors <- list(
  `Patient status` = c("Alive" = "lightblue", "Dead" = "#F8766D"),
  `OS month` = c("white", "#2BA25F")
)
heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

# 5. Generate Plot — no clustering anywhere; column order IS the survival order
png(file.path(output_dir, "lambda_by_survival.png"), width = 1600, height = 800, res = 120)
pheatmap(
  mtx,
  cluster_cols = FALSE,   # CRITICAL CHANGE: was TRUE -> hclust discarded any order
  cluster_rows = FALSE,   # features stay in fixed order: age, grade, lambda_1..q
  annotation_col = ann,
  annotation_colors = ann_colors,
  color = heatmap_colors,
  gaps_col = gaps,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize = 18
)
dev.off()
# }
