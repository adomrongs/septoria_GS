#!/bin/bash
#SBATCH --job-name=cv
#SBATCH --output=cv_%j.log
#SBATCH --error=cv_%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=20:00:00
#SBATCH --mem=50G
#SBATCH --array=1-10

module load R
srun Rscript code/R/7_runCV_all.R
