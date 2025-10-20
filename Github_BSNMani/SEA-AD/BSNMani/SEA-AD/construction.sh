#!/bin/bash
#SBATCH --job-name=array_R_job        
#SBATCH --output=/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SpaceX_BSNMani/building/codes/construction/construction_log/Oli/output_%A_%a.out
#SBATCH --error=/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SpaceX_BSNMani/building/codes/construction/construction_log/Oli/error_%A_%a.err
#SBATCH --array=1-14                  
#SBATCH --time=4:00:00               
#SBATCH --mem=16G
#SBATCH --partition=standard 
#SBATCH --account=lgarmire1
#SBATCH --cpus-per-task=8             

module load R/4.4.0

# 明确指定每个 array task 对应的参数组合
q_val_list=(2 2 3 3 4 4 5 5 6 6 7 7 8 8)
chain_list=(1 2 1 2 1 2 1 2 1 2 1 2 1 2)

q_val=${q_val_list[$SLURM_ARRAY_TASK_ID - 1]}
chain_idx=${chain_list[$SLURM_ARRAY_TASK_ID - 1]}

echo "Running job with q_val=$q_val and chain_idx=$chain_idx"

cd /nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/codes/SpaceX_BSNMani/building/codes

Rscript SpaceX_MERFISH_BSNMani_for_tong.R $q_val $chain_idx
