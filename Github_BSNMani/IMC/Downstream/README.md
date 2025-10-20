# Pipeline of breast cancer IMC data analysis using BSNMani

1. `BSNMani.R`: Run BSNMani on all 253 patients
2. `run_bsnmani.sh`: Slurm job to run `BSNMani.R`, allow for different q values and seeds
3. `ClusterLambda.R`: Do cluster (kmeans or coxph high-low risk stratification) using lambdas (or with clinical confounders)
4. `PlotBox.R`: Plot boxplot for comparison of different q
5. `PlotLambda.R`: Plot heatmap for lambda matrix vs survival groups
6. `PlotSubnetwork.R`: Plot heatmap for subnetwork correlation matrix
7. `PlotSurv.R`: Plot KM curve for survival groups derived from `ClusterLambda.R`
8. `Preprocess.R`: [**deprecated**] Preprocess BC data
9. `Clinical.R`: [**deprecated**] Run CoxPH on test split using trained model from Haowen,  summarize C-Index and pval
10. `Gsea.R`: Use EnrichR to do GSEA