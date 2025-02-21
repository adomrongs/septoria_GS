#!/bin/bash
#SBATCH --job-name=GWAS_wheat
#SBATCH --output=GWAS_wheat_%j.log
#SBATCH --error=GWAS_wheat_%j.err
#SBATCH --partition=fast
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --mem=40G
#SBATCH --array=1-4

module load R/4.2.0-foss-2021b
srun Rscript code/R/4_runGWAS_mix.R