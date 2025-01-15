library(tidyverse)
source('code/R/function_septoria_GS.R')


# I need blues of all individuals, genotype, kinship, map and phenotype

# Lets start by combining the two phenotypes, lets see how can we do this
load('data/modified_data/1_septoria.Rdata')
load("data/modified_data/5_predictions.Rdata")

phenotype_train <- cleaned_septoria_phenotype_1 |> 
  mutate(across(Isolate:BRep, as.factor))
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
formula <- "~ -1 + Isolate + Line + Trial + Leaf + Year + BRep"
blues_combined <- extract_blues_df_adapted(phenotypes_combined,
                         trait = c("PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion"),
                         formula = formula, 
                         colname = "Isolate")
blues_combined_clean <- blues_combined |> filter(Isolate != "Mock")


# save the information for the analysis
save(blues_combined_clean,
     genotype_all,
     map_septoria,
     k_all,
     file = 'data/modified_data/7_CV_sep_119.Rdata')
