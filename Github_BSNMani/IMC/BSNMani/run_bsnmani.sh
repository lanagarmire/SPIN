#!/usr/bin/env bash

#SBATCH --job-name=bsnmani_bc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --mem=4g
#SBATCH --mail-user=yangiwen@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=/nfs/dcmb-lgarmire/yangiwen/workspace/bsnmani/logs/%x-%j-%A-%a.log
#SBATCH --account=lgarmire99
#SBATCH --partition=standard,largemem
#SBATCH --array=0-34
#SBATCH --mail-type=ARRAY_TASKS

dir='/nfs/dcmb-lgarmire/yangiwen/workspace/bsnmani/r'
cd $dir || exit

qs=({2..8})
seeds=(0 42 64 123 894)
q=${qs[$((SLURM_ARRAY_TASK_ID / 5))]}
seed=${seeds[$((SLURM_ARRAY_TASK_ID % 5))]}
output_dir="/nfs/dcmb-lgarmire/yangiwen/workspace/bsnmani/output/bc/output_23/q${q}/s${seed}"

/sw/pkgs/arc/stacks/gcc/10.3.0/R/4.3.1/bin/Rscript R/bc/Bsnmani.R \
  -s "$seed" \
  -q "$q" \
  -o "$output_dir"
