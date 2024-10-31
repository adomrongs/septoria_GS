#!/bin/bash
#SBATCH --job-name=GWAS  
#SBATCH --output=GWAS_%j.log
#SBATCH --error=GWAS_%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=15:00:00
#SBATCH --mem=60G

module load R
srun Rscript code/R/03_runGWAS_septoria.R 
  