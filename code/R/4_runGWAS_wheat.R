library(tidyverse)
library(here())
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("code/R/function_septoria_GS.R")

load("data/modified_data/3_wheat_GWAS.Rdata")
i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
out_dir <- "outputs/GWAS_wheat"

for(j in seq_along(blues_wheat_list)){
  # select blues
  blues <- blues_wheat_list[[j]]
  name <- names(blues_wheat_list)[[j]]
  
  # adjust output directory
  out_dir_final <- file.path(out_dir, name, paste0("PC", i))
  dir.create(out_dir_final, recursive = TRUE, showWarnings = FALSE)
  setwd(out_dir_final)
  
  # match genotype and kinship
  genotype <- genotype_wheat[genotype_wheat[,1] %in% blues[,1], ]
  kinship <- k_wheat[
    k_wheat[, 1] %in% genotype[, 1], 
    c(TRUE, colnames(k_wheat)[-1] %in% genotype[, 1]) # Always keep the first column, filter remaining columns
  ]
  
  print(paste("Working on model", i, "and PC", j))

    myGAPIT <- GAPIT(
      Y = blues,
      GD = genotype,
      GM = map_wheat,
      KI = kinship,
      CV = NULL,
      PCA.total = i,
      model = c('Blink', 'FarmCPU', 'MLM')
    )
  # Reset working directory
  setwd(here())
  # Confirm progress
  print(paste("Completed model", name, "and PC", j))
}
