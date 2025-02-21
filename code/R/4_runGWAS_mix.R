library(tidyverse)
library(here())
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("code/R/function_septoria_GS.R")

load("data/modified_data/3_wheat_GWAS.Rdata")
i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
out_dir <- "outputs/GWAS_wheat"

# adjust output directory
out_dir_final <- file.path(out_dir, paste0("mix", i))
dir.create(out_dir_final, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir_final)
blues <- mix_blues[[i]]
genotype <- genotype_wheat[genotype_wheat[[1]] %in% blues[[1]], ]
k <- mix_k[[i]]

print(paste("Working on mix", i))

myGAPIT <- GAPIT(
  Y = blues,
  GD = genotype,
  GM = map_wheat,
  KI = k,
  CV = NULL,
  PCA.total = 3,
  model = c("Blink", 'FarmCPU', 'MLM')
)
# Reset working directory
setwd(here())
# Confirm progress
print(paste("Completed mix", i))