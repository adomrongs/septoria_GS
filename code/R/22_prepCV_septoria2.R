library(tidyverse)
source('code/R/function_septoria_GS.R')

load("data/modified_data/1_septoria.Rdata")

info <- read_csv('data/modified_data/info_strain_complete.csv')
gh1_pheno <- read_csv('data/modified_data/gh1_clean_pheno.csv') |> 
  left_join(info) |> 
  dplyr::select(Isolate:pycnidiaPerCm2Lesion, Region) |> 
  filter(Leaf == 2)

cultivars_phenotype <- gh1_pheno |> 
  split(gh1_pheno$Line)

formula <- "~ -1 + Isolate + Trial + Year + BRep + Region"
cultivars_blues <- map(cultivars_phenotype, \(x) extract_blues_df_adapted(phenotype = x,
                                                                          traits = c('PLACL', 'pycnidiaPerCm2Leaf', 'pycnidiaPerCm2Lesion'),
                                                                          formula = formula,
                                                                          colname = 'Isolate'))
k_inter <- k_septoria
rownames(k_inter) <- k_inter[,1]
k_inter <- k_inter[,-1]
colnames(k_inter) <- rownames(k_inter) 
kinship <- k_inter[rownames(k_inter) %in% unique(gh1_pheno$Isolate), 
                   colnames(k_inter) %in% unique(gh1_pheno$Isolate)]


save(cultivars_phenotype, cultivars_blues, kinship, file = 'data/modified_data/9_cvseptoria2.Rdata')

