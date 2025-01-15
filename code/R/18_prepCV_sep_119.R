library(tidyverse)
source('code/R/function_septoria_GS.R')


# I need blues of all individuals, genotype, kinship, map and phenotype

# Lets start by combining the two phenotypes, lets see how can we do this
load('data/modified_data/1_septoria.Rdata')
load("data/modified_data/5_predictions.Rdata")

phenotype_train <- cleaned_septoria_phenotype_1 |> 
  filter(Leaf == 2) |> 
  mutate(across(Isolate:BRep, as.factor)) |> 
  droplevels()
phenotype_test <- clean_test_pheno |> 
  dplyr::rename(Trial = TRep) |> 
  mutate(Leaf = 2,
         Year = 23) |> 
  dplyr::select(Isolate, Line, Trial, Leaf, Year, BRep, PLACL, pycnidiaPerCm2Leaf, pycnidiaPerCm2Lesion) |> 
  mutate(across(Isolate:BRep, as.factor),
         Trial = as.factor(as.numeric(Trial)),
         Isolate = dplyr::recode(Isolate,
                                 "CorVal_22" = "22_CorVal_L1",
                                 "ConilOrt_22" = "22_ConilOrt_L1"))

phenotypes_combined <- bind_rows(phenotype_train, phenotype_test) |> 
  relocate(Isolate, Line, Year, Trial, BRep, Leaf, PLACL, pycnidiaPerCm2Leaf, pycnidiaPerCm2Lesion)

# now we should extract blues for all the lines and exclude mock
formula <- "~ -1 + Isolate + Line + Trial + Year + BRep"
blues_combined <- extract_blues_df_adapted(phenotypes_combined,
                         trait = c("PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion"),
                         formula = formula, 
                         colname = "Isolate")
blues_combined_clean <- blues_combined |> filter(Isolate != "Mock")

names_df <- data.frame(
  Complete_Name = genotype_all[,1], 
  Name = unlist(map(str_split(genotype_all[,1], "_"), \(x) x[[2]])),
  Year = as.numeric(str_extract(genotype_all[,1], "([0-9]{2})_", group = T))
) 

diff <- setdiff(blues_combined_clean$Isolate, genotype_all[,1])
blues_combined_clean[blues_combined_clean$Isolate %in% diff, "Isolate"] <- names_df[names_df$Name %in% diff & names_df$Year == 23, "Complete_Name"]
genotype_all_clean <- genotype_all[genotype_all[,1] %in% blues_combined_clean$Isolate, ]

k_all_clean <- k_all[rownames(k_all) %in% blues_combined_clean$Isolate, 
                     rownames(k_all) %in% blues_combined_clean$Isolate]

# save the information for the analysis
save(blues_combined_clean,
     genotype_all_clean,
     map_septoria,
     k_all_clean,
     file = 'data/modified_data/7_CV_sep_119.Rdata')
