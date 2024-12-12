#!/bin/bash
#SBATCH --job-name=AIC
#SBATCH --output=AIC_%j.log
#SBATCH --error=AIC_%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=05:00:00
#SBATCH --mem=40G
#SBATCH --array=1-3

module load R
srun Rscript code/R/10_modelAIC.R