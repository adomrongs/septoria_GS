#!/bin/bash
#SBATCH --job-name=cv
#SBATCH --output=cv_%j.log
#SBATCH --error=cv_%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=20:00:00
#SBATCH --mem=30G
#SBATCH --array=1-30

module load R/4.2.0-foss-2021b
srun Rscript code/R/14_septoriaCV.R
