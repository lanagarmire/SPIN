#!/bin/bash
#SBATCH --job-name=SEA_AD_R_job        
#SBATCH --output=/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SpaceX_BSNMani/building/codes/Diagnosis/Diagnosis_log/Oli/output_%A_%a.log  
#SBATCH --error=/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SpaceX_BSNMani/building/codes/Diagnosis/Diagnosis_log/Oli/error_%A_%a.log   
#SBATCH --time=2:00:00                  
#SBATCH --mem=64G                       
#SBATCH --cpus-per-task=4              
#SBATCH --partition=standard           
#SBATCH --account=lgarmire0            
#SBATCH --array=2-8   # <--- 设置数组任务，从 q_val = 2 到 8

module load R/4.4.0 
module load Bioinformatics

cd /nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SpaceX_BSNMani/building/codes

# 将任务索引传给 R 脚本
Rscript SpaceX_MERFISH_BSNMani_diag_for_tong.R ${SLURM_ARRAY_TASK_ID}
