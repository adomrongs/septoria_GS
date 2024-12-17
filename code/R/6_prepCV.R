library(tidyverse)
library(rrBLUP)
source("code/R/function_septoria_GS.R")

# what do we need for the cross validation that we dont have yet
# - GRM mixes
# - I mix
# - I wheat

# lets prepare the GRM of the mixxes
load("data/modified_data/1_septoria.Rdata")
mix1 <- c("22_EcijaSec83Ica_L2", "22_EcijaSecCris_L1", "22_EcijaRegTej_L1")
mix2 <- c("22_CorKiko_L1", "22_CorCale_L1", "22_Cor3927_L1")
mix3 <- c("22_ConilAmi_L1", "22_Conil3806_L1", "22_Jerez3927_L1")
mix4 <- "22_CorVal_L1"

# Combine the mixes into a list
mixes <- list(mix1, mix2, mix3, mix4)

# Define the function to extract means
extractM_mix <- function(genotype, mix){
  tmp_M_mix <- genotype[genotype[, 1] %in% mix, ]
  means <- colMeans(tmp_M_mix[,-1])
  return(means)
}

# Apply the function to each mix using purrr::map
mixes_means <- purrr::map(mixes, ~ extractM_mix(genotype_septoria, .x))
M_mix <- bind_rows(mixes_means)

# calculate GRM of mixes
k_mixes <- A.mat(M_mix)
# calculate the identity of the mixes 
I_mixes <- diag(nrow(k_mixes))
rownames(k_mixes) <- colnames(k_mixes) <- rownames(I_mixes) <- colnames(I_mixes)<- c("mix1", "mix2", "mix3", "mix4")

# adjuist phenotype to stract blues for Plant:Strain
phenotype <- read_csv("data/modified_data/clean_phenotype_no_outliers.csv")
phenotype <- phenotype %>%
  mutate(PS = factor(paste(Plant,  Strain, sep = ":")))
formula_blues <- " ~ -1 + PS + Set + Rep + Leaf"
blues_PS <- extract_blues_df_adapted(phenotype,
                                     c("PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion"),
                                     formula_blues,
                                     "PS")
adjusted_phenotype <- blues_PS %>% 
  separate(col = PS, into = c("Plant", "Strain"), sep = ":") %>% 
  filter(!Strain %in% c("check1", "check2", "MOCK"))

# adjust k_wheat for GS
load("data/modified_data/2_wheat.Rdata")
intersect_lines <- intersect(genotype_wheat[,1], unique(adjusted_phenotype$Plant))
genotype_wheat <- genotype_wheat[genotype_wheat[,1] %in% intersect_lines, ]

rownames(k_wheat) <- k_wheat[,1]
k_wheat <- k_wheat[,-1]
colnames(k_wheat) <- gsub("\\.", "_", colnames(k_wheat))
k_wheat <- k_wheat %>% # adjust kinship to blups lines
  filter(rownames(k_wheat) %in% genotype_wheat$GenoID) %>% 
  dplyr::select(which(colnames(k_wheat) %in% genotype_wheat$GenoID))

# calculate the identity of the wheat lines
I_wheat <- diag(nrow(k_wheat))
rownames(I_wheat) <- colnames(I_wheat) <- rownames(k_wheat)

# veryfy dimensions
dim(k_mixes); dim(I_mixes); dim(k_wheat); dim(I_wheat)

#we have everything we need for the CV
save(k_mixes, I_mixes, k_wheat, I_wheat, adjusted_phenotype, map_wheat, genotype_wheat,
     file = "data/modified_data/4_CV_data.Rdata")
