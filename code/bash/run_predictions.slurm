#!/bin/bash
#SBATCH --job-name=predictions
#SBATCH --output=predictions_%j.log
#SBATCH --error=predictions_%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=12:00:00
#SBATCH --mem=30G
#SBATCH --array=1-3

module load R/4.2.0-foss-2021b
srun Rscript code/R/11_predictions.R
