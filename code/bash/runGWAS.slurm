#!/bin/bash
#SBATCH --job-name=GWAS
#SBATCH --output=GWAS_%j.log
#SBATCH --error=GWAS_%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=10:00:00
#SBATCH --mem=40G
#SBATCH --array=0-5

module load R
srun Rscript code/R/2_runGWAS_sep.R

