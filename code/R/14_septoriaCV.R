library(tidyverse)
library(sommer)
library(here)
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("code/R/function_septoria_GS.R")

load("data/modified_data/6_CV_septoria.Rdata")

i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
traits <- c("PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion")
test <- sample(rownames(kinship), ceiling(0.2 * nrow(kinship)))

for(j in seq_along(traits)){
  trait <- traits[[j]]
  dir <- paste0("data/modified_data/cv_septoria/", trait)
  dir.create(dir, recursive = T)
  
  message("Processing trait: ", trait, ", iteration: ", i)
  
  normal <- cv_septoria(genotype, phenotype, kinship, map, test, trait, blues_all, wModel = FALSE)
  weighted <- cv_septoria(genotype, phenotype, kinship, map, test, trait, blues_all, wModel = TRUE)
  
  save(normal, weighted, file = paste0(dir, "/iter_", i, ".Rdata"))
}
