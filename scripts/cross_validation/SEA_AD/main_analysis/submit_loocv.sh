#!/bin/bash
#SBATCH --job-name=BSNMani_LOOCV
#SBATCH --output=/home/koe/BSNMani_application-main/Github_BSNMani/Output/output_%A_%a.out
#SBATCH --error=/home/koe/BSNMani_application-main/Github_BSNMani/Error/error_%A_%a.err
#SBATCH --array=1-104                 # 26 folds * 4 q values
#SBATCH --time=150:00:00               # Adjust based on expected MCMC run time
#SBATCH --cpus-per-task=8             # MCMC is single-threaded
#SBATCH --mem=192G                      # Memory per job
#SBATCH --partition=amd-hdr100          # Update to your specific cluster partition

# Load necessary modules (adjust to your cluster's R module)
module load R/4.4.0-foss-2022b

# Calculate Fold (1 to 26) and Q (2 to 5) from the SLURM_ARRAY_TASK_ID
# Bash arithmetic uses 0-indexing for arrays
TASK_ID=$((SLURM_ARRAY_TASK_ID - 1))

# 4 candidate models per fold
FOLD=$(( (TASK_ID / 4) + 1 ))
Q_INDEX=$(( TASK_ID % 4 ))

# Array of candidate q values
Q_ARRAY=(2 3 4 5)
Q=${Q_ARRAY[$Q_INDEX]}

echo "Starting Array Task: $SLURM_ARRAY_TASK_ID"
echo "Assigned Fold: $FOLD | Assigned q: $Q"

# Execute the R script, passing the fold and q as arguments
Rscript /home/koe/BSNMani_application-main/Github_BSNMani/SEA-AD/BSNMani/SEA-AD/run_single_fold.R $FOLD $Q