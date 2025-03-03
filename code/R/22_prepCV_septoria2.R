library(tidyverse)
source('code/R/function_septoria_GS.R')

load("data/modified_data/1_septoria.Rdata")
load("data/modified_data/5_predictions.Rdata")

phenotype_test <- clean_test_pheno |> 
  dplyr::rename(Trial = TRep) |> 
  mutate(Leaf = 2,
         Year = 23) |> 
  dplyr::select(Isolate, Line, Trial, Leaf, Year, BRep, PLACL, pycnidiaPerCm2Leaf, pycnidiaPerCm2Lesion) |> 
  mutate(across(Isolate:BRep, as.factor),
         Trial = as.factor(as.numeric(Trial)),
         Isolate = dplyr::recode(Isolate,
                                 "CorVal_22" = "22_CorVal_L1",
                                 "ConilOrt_22" = "22_ConilOrt_L1"), 
         Region = 'Córdoba', 
         Line = case_when(
           Line == 'Don_Ricardo' ~ 'Don Ricardo',
           TRUE ~ Line  # Mantiene los demás valores sin cambios
         ))

info <- read_csv('data/modified_data/info_strain_complete.csv')
gh1_pheno <- read_csv('data/modified_data/gh1_clean_pheno.csv') |> 
  left_join(info) |> 
  dplyr::select(Isolate:pycnidiaPerCm2Lesion, Region) |> 
  filter(Leaf == 2) |> 
  mutate(across(Isolate:BRep, as.factor))

phenotypes_combined <- bind_rows(gh1_pheno, phenotype_test) |> 
  filter(Isolate != "Mock") |>  # Replace "Mock" with NA
  droplevels() |> 
  relocate(Isolate, Line, Year, Trial, BRep, Leaf, PLACL, pycnidiaPerCm2Leaf, pycnidiaPerCm2Lesion)
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


cultivars_phenotype <- phenotypes_combined |> 
  split(phenotypes_combined$Line)

formula <- "~ -1 + Isolate + Trial + Year + BRep + Region"
cultivars_blues <- map(cultivars_phenotype, \(x) extract_blues_df_adapted(phenotype = x,
                                                                          traits = c('PLACL', 'pycnidiaPerCm2Leaf', 'pycnidiaPerCm2Lesion'),
                                                                          formula = formula,
                                                                          colname = 'Isolate'))

kinship <- k_all[rownames(k_all) %in% unique(phenotypes_combined$Isolate), 
                     colnames(k_all) %in% unique(phenotypes_combined$Isolate)]


save(cultivars_phenotype, cultivars_blues, kinship, file = 'data/modified_data/9_cvseptoria2.Rdata')

