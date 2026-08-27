##########################
## This R script is used for benchmark BSNMani - Lasso - Elastic Net
## And draw bar plots with error term
##########################

rm(list = ls())
library(readxl)
library(readxl)
library(car)
library(pROC)
library(dplyr)
setwd("D:\\Research_project\\For Test\\Elastic_net_Lasso\\For_benchmark")

benchmark_res = readRDS("benchmark_res_uni.RDS")
bsnmani = read.csv("LOOCV_prediction_results_q_3_hdWGCNA_uni_roc.csv")
clinical_df = read_excel("meta_uni.xlsx")
lasso_prob_vec = benchmark_res$lasso$prob
Elastic_prob_vec = benchmark_res$enet$prob
bsnmani_prob_vec = bsnmani$predicted_prob


############################################
## Now we have obtain prob of 27 patients with 3 methods
############################################
predicted_probs = bsnmani_prob_vec
prob_df = cbind(lasso_prob_vec,Elastic_prob_vec,predicted_probs,clinical_df$Cognitive_Status)
colnames(prob_df) = c("Lasso","Elastic_net","BSNMani","Cognitive_Status")

######################################
## Boostrap and error term calculation

metrics_once <- function(y, prob, thr = 0.50) {
  pred <- as.integer(prob >= thr)
  TP <- sum(pred == 1 & y == 1)
  FP <- sum(pred == 1 & y == 0)
  TN <- sum(pred == 0 & y == 0)
  FN <- sum(pred == 0 & y == 1)
  acc  <- (TP + TN) / length(y)
  prec <- if ((TP + FP) > 0) TP/(TP + FP) else NA
  rec  <- if ((TP + FN) > 0) TP/(TP + FN) else NA
  spec <- if ((TN + FP) > 0) TN/(TN + FP) else NA
  f1   <- if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0) 2*prec*rec/(prec+rec) else NA
  auc  <- as.numeric(auc(roc(y, prob, quiet=TRUE)))
  c(Accuracy=acc, Precision=prec, Recall=rec, F1=f1, Specificity=spec, AUC=auc)
}

metrics_once_uni <- function(y, prob, thr = 0.70) {
  pred <- as.integer(prob >= thr)
  TP <- sum(pred == 1 & y == 1)
  FP <- sum(pred == 1 & y == 0)
  TN <- sum(pred == 0 & y == 0)
  FN <- sum(pred == 0 & y == 1)
  acc  <- (TP + TN) / length(y)
  prec <- if ((TP + FP) > 0) TP/(TP + FP) else NA
  rec  <- if ((TP + FN) > 0) TP/(TP + FN) else NA
  spec <- if ((TN + FP) > 0) TN/(TN + FP) else NA
  f1   <- if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0) 2*prec*rec/(prec+rec) else NA
  auc  <- as.numeric(auc(roc(y, prob, quiet=TRUE)))
  c(Accuracy=acc, Precision=prec, Recall=rec, F1=f1, Specificity=spec, AUC=auc)
}

bootstrap_metrics <- function(y, prob, B = 5000, thr = 0.50, conf = 0.9, keep_samples = TRUE, seed = 123) {
  n <- length(y)
  set.seed(seed)
  boot_mat <- t(replicate(B, {
    idx <- sample.int(n, replace=TRUE)
    metrics_once(y[idx], prob[idx], thr)
  }))
  vals <- colMeans(boot_mat, na.rm=TRUE)
  lohi <- apply(boot_mat, 2, function(v) quantile(v, c((1-conf)/2, 1-(1-conf)/2), na.rm=TRUE))
  summary_df <- data.frame(
    metric = names(vals),
    value  = as.numeric(vals),
    lo     = as.numeric(lohi[1,]),
    hi     = as.numeric(lohi[2,]),
    row.names = NULL
  )
  original_vec <- metrics_once(y, prob, thr)
  
  out <- list(summary = summary_df, original = original_vec)
  
  if (keep_samples) {
    samples_wide <- as.data.frame(boot_mat)
    colnames(samples_wide) <- names(original_vec)
    samples_wide$replicate <- seq_len(B)
    samples_long <- reshape(
      samples_wide,
      direction = "long",
      varying = names(original_vec),
      v.names = "value",
      times = names(original_vec),
      timevar = "metric",
      idvar = "replicate"
    )
    rownames(samples_long) <- NULL
    samples_long$metric <- as.character(samples_long$metric)
    out$samples <- samples_long[, c("replicate","metric","value")]
  }
  return(out)
}

bootstrap_metrics_uni <- function(y, prob, B = 5000, thr = 0.70, conf = 0.9, keep_samples = TRUE, seed = 123) {
  n <- length(y)
  set.seed(seed)
  boot_mat <- t(replicate(B, {
    idx <- sample.int(n, replace=TRUE)
    metrics_once_uni(y[idx], prob[idx], thr)
  }))
  vals <- colMeans(boot_mat, na.rm=TRUE)
  lohi <- apply(boot_mat, 2, function(v) quantile(v, c((1-conf)/2, 1-(1-conf)/2), na.rm=TRUE))
  summary_df <- data.frame(
    metric = names(vals),
    value  = as.numeric(vals),
    lo     = as.numeric(lohi[1,]),
    hi     = as.numeric(lohi[2,]),
    row.names = NULL
  )
  original_vec <- metrics_once_uni(y, prob, thr)
  
  out <- list(summary = summary_df, original = original_vec)
  
  if (keep_samples) {
    samples_wide <- as.data.frame(boot_mat)
    colnames(samples_wide) <- names(original_vec)
    samples_wide$replicate <- seq_len(B)
    samples_long <- reshape(
      samples_wide,
      direction = "long",
      varying = names(original_vec),
      v.names = "value",
      times = names(original_vec),
      timevar = "metric",
      idvar = "replicate"
    )
    rownames(samples_long) <- NULL
    samples_long$metric <- as.character(samples_long$metric)
    out$samples <- samples_long[, c("replicate","metric","value")]
  }
  return(out)
}

y <- prob_df[,"Cognitive_Status"]

res_las  <- bootstrap_metrics(y, prob_df[,"Lasso"])
res_enet <- bootstrap_metrics(y, prob_df[,"Elastic_net"])
res_bsn  <- bootstrap_metrics_uni(y, prob_df[,"BSNMani"])

res_las_df  <- transform(res_las$summary,  method = "Lasso")
res_enet_df <- transform(res_enet$summary, method = "Elastic Net")
res_bsn_df  <- transform(res_bsn$summary,  method = "BSNMani")

plot_df <- rbind(res_las_df, res_enet_df, res_bsn_df)
rownames(plot_df) <- NULL

#############################
## pair test
#############################
res_bsn_sample = res_bsn$samples
res_las_sample = res_las$samples
res_enet_sample = res_enet$samples
metrics <- c("Accuracy","Precision","Recall","F1","Specificity","AUC")


pairwise_pvals <- function(samples_A, samples_B, name_A, name_B,
                           metrics, paired = TRUE, method = "wilcox") {
  joined <- samples_A %>%
    inner_join(samples_B, by = c("replicate","metric"), suffix = c(".A",".B")) %>%
    filter(metric %in% metrics)
  
  res <- lapply(metrics, function(m) {
    dat <- joined %>% filter(metric == m)
    x <- dat$value.A
    y <- dat$value.B
    # Wilcoxon/Mann-Whitney
    if (method == "wilcox") {
      tst <- wilcox.test(x, y, paired = paired, exact = FALSE)
    } else if (method == "ttest") {
      tst <- t.test(x, y, paired = paired)
    } else {
      stop("Unsupported method.")
    }
    data.frame(
      metric = m,
      comparison = paste(name_A, "vs", name_B),
      n_pairs = nrow(dat),
      p_value = unname(tst$p.value)
    )
  })
  do.call(rbind, res)
}

p_bsn_vs_enet <- pairwise_pvals(res_bsn_sample, res_enet_sample,
                                "BSNMani", "Elastic Net",
                                metrics, paired = TRUE, method = "ttest")

p_bsn_vs_lasso <- pairwise_pvals(res_bsn_sample, res_las_sample,
                                 "BSNMani", "Lasso",
                                 metrics, paired = TRUE, method = "ttest")
p_all <- bind_rows(p_bsn_vs_enet, p_bsn_vs_lasso) %>%
  group_by(comparison) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(match(metric, metrics), comparison)

p_all


p2star <- function(p) {
  ifelse(p < 0.001, "*",
         ifelse(p < 0.01,  "*",
                ifelse(p < 0.05,  "*", "ns")))
}
p_all$stars <- p2star(p_all$p_value)

library(dplyr)
y_pos <- plot_df %>%
  group_by(metric) %>%
  summarise(y = max(hi, na.rm = TRUE) + 0.04, .groups = "drop")

metric_levels <- unique(plot_df$metric)                 
x_index <- setNames(seq_along(metric_levels), metric_levels)
dodge_width <- 0.8
n_groups <- 3
offset <- dodge_width / n_groups                        # ≈ 0.2667

p_plot <- p_all %>% left_join(y_pos, by = "metric")
p_enet <- p_plot %>%
  filter(comparison == "BSNMani vs Elastic Net") %>%
  mutate(xmid = x_index[metric],
         xmin = xmid - offset,
         xmax = xmid + 0,
         y_position = y + 0.00,
         pair = "ENET") %>%    
  select(metric, xmin, xmax, y_position, stars, pair)

p_lasso <- p_plot %>%
  filter(comparison == "BSNMani vs Lasso") %>%
  mutate(xmid = x_index[metric],
         xmin = xmid - offset,
         xmax = xmid + offset,
         y_position = y + 0.05,
         pair = "LASSO") %>%   
  select(metric, xmin, xmax, y_position, stars, pair)

y_global <- max(plot_df$hi, na.rm = TRUE) + 0.04
gap <- 0.05  

p_enet <- p_enet %>% mutate(y_position = y_global)
p_lasso <- p_lasso %>% mutate(y_position = y_global + gap)

sign_df <- bind_rows(p_enet, p_lasso)
sign_df$stars <- as.character(sign_df$stars)



library(ggsignif)
ylim_up <- min(1.2, max(sign_df$y_position) + 0.02)

p <- ggplot(plot_df, aes(x = metric, y = value, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = lo, ymax = hi),
                position = position_dodge(width = 0.8), width = 0.25) +
  labs(x = NULL, y = "Score", fill = "Method",
       title = "LOOCV Benchmark") +
  theme_minimal(base_size = 14) +
  scale_y_continuous(limits = c(0, ylim_up), expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(clip = "off") +
  ggsignif::geom_signif(
    data = sign_df,
    mapping = aes(xmin = xmin, xmax = xmax, y_position = y_position,
                  annotations = stars, group = interaction(metric, pair)),
    manual = TRUE, inherit.aes = FALSE,
    tip_length = 0.01, textsize = 5, vjust = 0.3
  )
p

setwd("D:\\Research_project\\For Test\\Elastic_net_Lasso\\For_benchmark")
ggsave("LOOCV_benchmark_uni_newest.png", plot = p,
       width = 180, height = 120, units = "mm",
       dpi = 600, device = ragg::agg_png, bg = "white")
