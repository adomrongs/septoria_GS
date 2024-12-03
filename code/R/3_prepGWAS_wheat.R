library(tidyverse)
library(lme4)
source("code/R/function_septoria_GS.R")

load("data/modified_data/2_wheat.Rdata")
phenotype <- read_csv("data/modified_data/clean_phenotype.csv")
summary(phenotype)

traits <- c("PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion")
cleaned_phenotype <- remove_outliers(df = phenotype, cols = traits)
write_csv(cleaned_phenotype, file = "data/modified_data/clean_phenotype_no_outliers.csv")
dim(phenotype); dim(cleaned_phenotype)

# First we will calculate the blups including the Leaf in the model
formula <- "~ -1 + Plant + Strain + Rep + Leaf "
blues_wheat <- extract_blues_df(cleaned_phenotype, traits, formula)

# Additionally we will also calcualte the blups when using just leaf 2
cleaned_phenotype_leaf2 <- cleaned_phenotype %>% 
  filter(Leaf == 2)
formula2 <- "~ -1 + Plant + Strain + Rep"
blues_leaf2 <- extract_blues_df(cleaned_phenotype_leaf2, traits, formula2)

# prepare the data to run the GWAS
blues_wheat_list <- list(blues_wheat, blues_leaf2)
colnames(k_wheat) <- c("GenoID", k_wheat[,1])
map_wheat$Position <- as.numeric(map_wheat$Position)

names(blues_wheat_list) <- c("All_leaves", "Leaf2")
save(genotype_wheat, map_wheat, k_wheat, blues_wheat_list,
     file = "data/modified_data/3_wheat_GWAS.Rdata")

