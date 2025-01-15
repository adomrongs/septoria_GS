library(tidyverse)
library(sommer)
library(here)
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("code/R/function_septoria_GS.R")

i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
load('data/modified_data/7_CV_sep_119.Rdata')

traits <- c("PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion")
test <- sample(rownames(k_all_clean), ceiling(0.2 * nrow(k_all_clean)))

dir <- "data/modified_data/cv_septoria_119/"
dir.create(dir, recursive = T, showWarnings = FALSE)

normal_list <- list()
weighted_list <- list()
for(j in seq_along(traits)){
  trait <- traits[[j]]
  
  message("Processing trait: ", trait, ", iteration: ", i)
  
  normal_list[[paste0("iter_", i)]][[trait]] <- cv_septoria2(genotype = genotype_all_clean,
                                                             blues_all = blues_combined_clean,
                                                             kinship = k_all_clean,
                                                             map = map_septoria,
                                                             test = test,
                                                             trait = trait,
                                                             wModel = FALSE)
  
  weighted_list[[paste0("iter_", i)]][[trait]] <- cv_septoria2(genotype = genotype_all_clean,
                                                               blues_all = blues_combined_clean,
                                                               kinship = k_all_clean,
                                                               map = map_septoria,
                                                               test = test,
                                                               trait = trait
                                                               , wModel = TRUE)
  
  message("Analysis Ready")
}

results_lits <- list(normal_list, weighted_list, test)
save(results_lits, file = paste0(dir, "iter_", i, ".Rdata"))