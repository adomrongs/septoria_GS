#!/bin/bash
#SBATCH --job-name=combine
#SBATCH --output=combine_%j.log
#SBATCH --error=combine_%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --time=05:00:00
#SBATCH --mem=60G

module load R
srun Rscript code/R/9_prepTest.R
