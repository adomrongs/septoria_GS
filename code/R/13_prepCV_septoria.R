library(tidyverse)
library(sommer)
library(here)
source("code/R/function_septoria_GS.R")

load("data/modified_data/5_predictions.Rdata")
load("data/modified_data/1_septoria.Rdata")

genotype <- genotype_septoria
phenotype <- cleaned_septoria_phenotype_1 |> 
  filter(Leaf == 2) |> 
  droplevels()
kinship <- k_septoria[,-1]
colnames(kinship) <- rownames(kinship) <- k_septoria[,1]
map <- map_septoria
formula <- "~ -1 + Isolate + Line + Trial + Year + BRep"
blues_all <- extract_blues_df_adapted(cleaned_septoria_phenotype_1,
                                      trait = c("PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion"),
                                      formula = formula, 
                                      colname = "Isolate")

save(genotype, map_septoria,phenotype, kinship, blues_all, file = "data/modified_data/6_CV_septoria.Rdata")
