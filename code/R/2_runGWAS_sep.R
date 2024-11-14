library(tidyverse)
library(here)
source("http://zzlab.net/GAPIT/gapit_functions.txt")

load("data/1_septoria_data.Rdata")

# Retrieve iteration ID from SLURM
iter <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
dir <- paste0("outputs/GWAS_sep/PC_", iter)
if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
setwd(dir)

print(paste0("Working on PC_", iter))

# Run GAPIT with error handling
tryCatch({
  myGAPIT <- GAPIT(
    Y  = septoria_pheno,
    GD = septoria_geno,
    GM = septoria_map,
    CV = NULL,
    PCA.total = iter,
    model = c('Blink', 'FarmCPU')
  )
}, error = function(e) {
  message(paste("Error in iteration", iter, ":", e$message))
})

# Restore the working directory
setwd(here())