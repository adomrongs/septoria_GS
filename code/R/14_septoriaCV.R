library(tidyverse)
library(sommer)
library(here)
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("code/R/function_septoria_GS.R")

i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
load("data/modified_data/6_CV_septoria.Rdata")

traits <- c("PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion")
test <- sample(rownames(kinship), ceiling(0.2 * nrow(kinship)))

dir <- "data/modified_data/cv_septoria/"
dir.create(dir, recursive = T, showWarnings = FALSE)

normal_list <- list()
weighted_list <- list()
for(j in seq_along(traits)){
  trait <- traits[[j]]
  
  message("Processing trait: ", trait, ", iteration: ", i)
  
  normal_list[[paste0("iter_", i)]][[trait]] <- cv_septoria(genotype, phenotype, kinship,
                        map_septoria, test, trait, blues_all, wModel = FALSE)
  
  weighted_list[[paste0("iter_", i)]][[trait]] <- cv_septoria(genotype, phenotype, kinship,
                          map_septoria, test, trait, blues_all, wModel = TRUE)
  
  message("Analysis Ready")
}

results_lits <- list(normal_list, weighted_list)
save(results_lits, file = paste0(dir, "iter_", i, ".Rdata"))
