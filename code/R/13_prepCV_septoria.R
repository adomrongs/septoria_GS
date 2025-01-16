library(tidyverse)
library(sommer)
library(here)
source("code/R/function_septoria_GS.R")

load("data/modified_data/5_predictions.Rdata")
load("data/modified_data/1_septoria.Rdata")

phenotype <- cleaned_septoria_phenotype_1 |> 
  filter(Leaf == 2) |> 
  droplevels()
genotype <- genotype_septoria[genotype_septoria[,1] %in% unique(phenotype$Isolate), ]
k_inter <- k_septoria
rownames(k_inter) <- k_inter[,1]
k_inter <- k_inter[,-1]
colnames(k_inter) <- rownames(k_inter) 
kinship <- k_inter[rownames(k_inter) %in% unique(phenotype$Isolate), 
                      colnames(k_inter) %in% unique(phenotype$Isolate)]
map <- map_septoria

dim(kinship); dim(phenotype); dim(genotype); dim(map_septoria)

formula <- "~ -1 + Isolate + Line + Trial + Year + BRep"
blues_all <- extract_blues_df_adapted(phenotype,
                                      trait = c("PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion"),
                                      formula = formula, 
                                      colname = "Isolate")

save(genotype, map_septoria,phenotype, kinship, blues_all, file = "data/modified_data/6_CV_septoria.Rdata")
