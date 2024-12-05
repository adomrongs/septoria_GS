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

# prepare the data to run the GWAS
colnames(k_wheat) <- c("GenoID", k_wheat[,1])
map_wheat$Position <- as.numeric(map_wheat$Position)
k_wheat <- k_wheat[k_wheat$GenoID %in% blues_wheat$GenoID, ]
k_wheat <- k_wheat[, colnames(k_wheat) %in% c("GenoID", blues_wheat$GenoID)]
genotype_wheat <- genotype_wheat[genotype_wheat[, "GenoID"] %in% blues_wheat$GenoID,]

dim(genotype_wheat); dim(map_wheat); dim(k_wheat); dim(blues_wheat)
str(blues_wheat)
save(genotype_wheat, map_wheat, k_wheat, blues_wheat,
     file = "data/modified_data/3_wheat_GWAS.Rdata")

