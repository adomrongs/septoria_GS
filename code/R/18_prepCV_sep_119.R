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
  filter(Isolate != "Mock") |>  # Replace "Mock" with NA
  droplevels() |> 
  relocate(Isolate, Line, Year, Trial, BRep, Leaf, PLACL, pycnidiaPerCm2Leaf, pycnidiaPerCm2Lesion)

# now we should extract blues for all the lines and exclude mock
names_df <- data.frame(
  Complete_Name = genotype_all[,1],
  Name = unlist(map(str_split(genotype_all[,1], "_"), \(x) x[[2]])),
  Year = as.numeric(str_extract(genotype_all[,1], "([0-9]{2})_", group = T))
)
diff <- setdiff(unique(phenotypes_combined$Isolate), genotype_all[,1])

phenotypes_combined <- phenotypes_combined %>%
  left_join(
    names_df %>% filter(Year == 23) %>% dplyr::select(Name, Complete_Name), 
    by = c("Isolate" = "Name")
  ) %>%
  mutate(Isolate = ifelse(Isolate %in% diff, Complete_Name, Isolate)) %>% 
  dplyr::select(-Complete_Name)  

genotype_all_clean <- genotype_all[genotype_all[,1] %in% unique(phenotypes_combined$Isolate), ]
k_all_clean <- k_all[rownames(k_all) %in% unique(phenotypes_combined$Isolate), 
                      colnames(k_all) %in% unique(phenotypes_combined$Isolate)]

formula <- "~ -1 + Isolate + Line + Trial + Year + BRep"
blues_combined <- extract_blues_df_adapted(phenotypes_combined,
                                           trait = c("PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion"),
                                           formula = formula,
                                           colname = "Isolate")
# blues_combined_clean <- blues_combined |> filter(Isolate != "Mock")

dim(phenotypes_combined); dim(k_all_clean); dim(genotype_all_clean); dim(map_septoria); dim(blues_combined)
length(intersect(rownames(k_all_clean), unique(phenotypes_combined$Isolate)))

# save the information for the analysis
save(phenotypes_combined,
     genotype_all_clean,
     map_septoria,
     k_all_clean,
     blues_combined, 
     file = 'data/modified_data/7_CV_sep_119.Rdata')
