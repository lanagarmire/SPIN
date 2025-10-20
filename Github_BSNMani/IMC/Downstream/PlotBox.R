#' Plot comparison of different q-values in BSNMani

# Plot using all patients

q_vals <- 2:8
seeds <- c(0, 42, 64, 123, 894)
output_dir <- "/nfs/dcmb-lgarmire/yangiwen/workspace/bsnmani/output/bc/output_21"

# Read BSNMani prediction for test sets and calculate c-index
results <- setNames(lapply(q_vals, function(q_val) {
  setNames(lapply(seeds, function(seed) {
    clinical_df_test <- readRDS(file.path(output_dir, paste0("q", q_val), paste0("s", seed), "test", "g2_survival", "clinical_df_test_test.RDS"))
    lp_test <- readRDS(file.path(output_dir, paste0("q", q_val), paste0("s", seed), "test", "g2_survival", "lp_test_test.RDS"))
    clinical_df_test$pred <- lp_test
    surv_obj_test <- survival::Surv(clinical_df_test$OSmonth, clinical_df_test$Patientstatus)
    cindex_test <- survival::concordance(surv_obj_test ~ pred, data = clinical_df_test, reverse = T)$concordance
    data.frame(q = q_val, seed = seed, cindex = cindex_test)
  }), seeds)
}), q_vals)
results <- do.call(rbind, lapply(names(results), function(q) {
  do.call(rbind, results[[q]])
}))

# Save box plot
png(file.path(output_dir, "q_comparison.png"), width = 1024, height = 1024, res = 120)
aggregate(cindex ~ q, data = results, FUN = max)
selected_q <- 2
pairs <- lapply(setdiff(2:8, selected_q), function(x) c(selected_q, x))
print(ggplot2::ggplot(results, ggplot2::aes(x = q, y = cindex, group = q)) +
  ggplot2::geom_boxplot(fill = "#69b3a2", color = "#1f3552") +
  # ggplot2::geom_boxplot(fill = "#dce9cb", color = "#1f3552") +
  ggsignif::geom_signif(
    comparisons = pairs,
    map_signif_level = T,
    test = "wilcox.test",
    test.args = list(exact = F),
    step_increase = 0.1
  ) +
  ggplot2::scale_x_discrete(limits = 2:8) +
  # ggplot2::ylim(0.65, 0.78) +
  ggplot2::labs(x = "q", y = "C-index") +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 26, face = "bold"),
    axis.text = ggplot2::element_text(size = 20),
    plot.title = ggplot2::element_text(size = 26, face = "bold")
  )
)
dev.off()
