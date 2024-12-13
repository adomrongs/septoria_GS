#!/bin/bash
#SBATCH --job-name=predictions
#SBATCH --output=predictions_%j.log
#SBATCH --error=predictions_%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --mem=30G
#SBATCH --array=1-3

module load R
srun Rscript code/R/11_predictions.R
