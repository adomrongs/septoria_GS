library(tidyverse)
library(here)
source("http://zzlab.net/GAPIT/gapit_functions.txt")

load("data/1_septoria_data.Rdata")
dim(septoria_pheno); dim(septoria_geno); dim(septoria_map); dim(septoria_kinship)
str(septoria_pheno); str(septoria_geno[1:5, 1:5]); str(septoria_map); str(septoria_kinship[1:5, 1:5])

# Retrieve iteration ID from SLURM
iter <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
dir <- paste0("outputs/GWAS_sep/PC_", iter)
if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
setwd(dir)

print(paste0("Working on PC_", iter))

  myGAPIT <- GAPIT(
    Y  = septoria_pheno,
    GD = septoria_geno,
    GM = septoria_map,
    KI = septoria_kinship,
    CV = NULL,
    PCA.total = iter,
    model = c('Blink', 'FarmCPU')
  )

# Restore the working directory
setwd(here())
