#' Plot hazard ratio
source("/nfs/dcmb-lgarmire/yangiwen/workspace/common/Utils.R")

q_vals <- 2:8
seeds <- c(0, 42, 64, 123, 894)
output_dir <- "/nfs/dcmb-lgarmire/yangiwen/workspace/bsnmani/output/bc/output_22"

for (q in q_vals) {
  for (seed in seeds) {
    seed_output_dir <- file.path(output_dir, paste0("q", q), paste0("s", seed))
    g2_survival_dir <- file.path(seed_output_dir, "train", "g2_survival")
    fig_dir <- file.path(seed_output_dir, "fig")
    mkdir(fig_dir)
    model <- readRDS(file.path(g2_survival_dir, "g2_survival_res_train.RDS"))
    clinical_df_train <- readRDS(file.path(g2_survival_dir, "clinical_df_train_train.RDS"))
    png(file.path(fig_dir, "hazard.png"), width = 768, height = 768, res = 120)
    print(survminer::ggforest(model, data = clinical_df_train))
    dev.off()
  }
}
