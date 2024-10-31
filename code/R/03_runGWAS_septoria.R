library(tidyverse)
library(here)
source("http://zzlab.net/GAPIT/gapit_functions.txt")

load("data/septoria_data.Rdata")

for(i in 1:5){
  new_dir <- paste0("outputs/GWAS_PC", i)
  dir.create(new_dir)
  setwd(new_dir)
  tmp <- capture.output({
    scores <- GAPIT(Y = BLUPs_df,
                    GD = geno_imputed,
                    GM = map,
                    KI = K_sep,
                    CV = NULL,
                    PCA.total = i,
                    model = "Blink",
                    file.output = T)
  })
  rm(tmp)
  gc()
  setwd(here())
}
