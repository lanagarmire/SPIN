#!/bin/bash
#SBATCH --job-name=SEA_AD_R_job        
#SBATCH --output=/nfs/turbo/umms-lgarmire/liutong/BSNMani/SpaceX_BSNMani/building/codes/Diagnosis/Diagnosis_log/Oli_uni/v4output_%A_%a.log  
#SBATCH --error=/nfs/turbo/umms-lgarmire/liutong/BSNMani/SpaceX_BSNMani/building/codes/Diagnosis/Diagnosis_log/Oli_uni/v4error_%A_%a.log   
#SBATCH --time=2:00:00                  
#SBATCH --mem=32G                       
#SBATCH --cpus-per-task=2              
#SBATCH --partition=standard            
#SBATCH --account=lgarmire1            
#SBATCH --array=2-5   

module load R/4.4.0 
module load Bioinformatics

cd /nfs/turbo/umms-lgarmire/liutong/BSNMani/SpaceX_BSNMani/building/codes


Rscript SpaceX_MERFISH_BSNMani_diag_for_tong.R ${SLURM_ARRAY_TASK_ID}
