#!/bin/bash
#SBATCH --job-name=predictions
#SBATCH --output=predictions_%j.log
#SBATCH --error=predictions_%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=05:00:00
#SBATCH --mem=60G

module load R
srun Rscript code/R/10_predictions.R
