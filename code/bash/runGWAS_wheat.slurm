#!/bin/bash
#SBATCH --job-name=GWAS_wheat
#SBATCH --output=GWAS_wheat_%j.log
#SBATCH --error=GWAS_wheat_%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=20:00:00
#SBATCH --mem=40G
#SBATCH --array=0-4

module load R
srun Rscript code/R/4_runGWAS_wheat.R