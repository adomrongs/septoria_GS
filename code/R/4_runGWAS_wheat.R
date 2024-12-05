library(tidyverse)
library(here())
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("code/R/function_septoria_GS.R")

load("data/modified_data/3_wheat_GWAS.Rdata")
i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
out_dir <- "outputs/GWAS_wheat"

# adjust output directory
out_dir_final <- file.path(out_dir, paste0("PC", i))
dir.create(out_dir_final, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir_final)
  
  
print(paste("Working on PC", i))

myGAPIT <- GAPIT(
      Y = blues_wheat,
      GD = genotype_wheat,
      GM = map_wheat,
      KI = k_wheat,
      CV = NULL,
      PCA.total = i,
      model = c("Blink", "MLM")
)
  # Reset working directory
setwd(here())
  # Confirm progress
print(paste("Completed model PC", i))
