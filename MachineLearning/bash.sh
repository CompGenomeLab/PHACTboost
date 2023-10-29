#!/bin/bash
#SBATCH --job-name=PHACTboost
#SBATCH --account=investor
#SBATCH --nodes=1
#SBATCH --array=1
#SBATCH --mem=16G
#SBATCH --qos=long_investor
#SBATCH --partition=long_investor
#SBATCH --time=7-00:00:00
#SBATCH --output=./output/%A_%a-slurm.out	
#SBATCH --mail-type=ALL

module load r-4.1.0-gcc-9.2.0-wq6f22q 
echo "======================================================================================"
Rscript lgb_Final.R ${SLURM_ARRAY_TASK_ID} "CountNodes_3" "PHACTboost_Result" 
