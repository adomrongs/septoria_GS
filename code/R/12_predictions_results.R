library(tidyverse)
library(lme4)
source("code/R/function_septoria_GS.R")

load("data/modified_data/5_predictions.Rdata")

# retrieve information about the strains name
info_strains <- read_csv("data/raw_data/INFO_STRAINS.csv")
test <- info_strains %>% 
  filter(Partition == "Test") %>% 
  pull(Isolate)
test_name <- unlist(map(strsplit(test, "_"), ~.x[2]))

isolate_conversion <- info_strains %>%
  filter(Partition == "Test") %>%
  dplyr::select(Isolate_complete = Isolate) %>% 
  mutate(Isolate_incomplete = map_chr(strsplit(Isolate_complete, "_"), ~ .x[2]))

# extract the blues for the test lines 
formula <- "~ -1 + Isolate + Line + TRep + BRep"
blues_test <- extract_blues_df_adapted(clean_test_pheno,
                                   traits = c("PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion"),
                                   formula = formula,
                                   colname = "Isolate")
blues_test <- blues_test %>% 
  filter(Isolate %in% test_name) %>% 
  left_join(isolate_conversion, by = c("Isolate" = "Isolate_incomplete")) %>% 
  dplyr::select(Isolate = Isolate_complete, PLACL, pycnidiaPerCm2Leaf, pycnidiaPerCm2Lesion) %>% 
  arrange(Isolate)
  

# extract blups from the different models test as well as pohenotypes
dirs <- list.dirs("data/modified_data/predictions", full.names = T, recursive = F)
blups_list <- list()
for(dir in dirs){
  subdirs <- list.dirs(dir, full.names = T, recursive = F)
  for(subdir in subdirs){
    files <- list.files(subdir, full.names = T, recursive = F)
    for(file in files){
        load(file)
        blups_list[[basename(dir)]][[basename(subdir)]][[gsub(".Rdata", "", basename(file))]] <- blups_df
    }
  }
}

# Crear los data frames para cada trait
# Apply the function to each trait
PLACL_df <- process_and_filter("PLACL", test)
pycnidiaLeaf_df <- process_and_filter("pycnidiaPerCm2Leaf", test)
pycnidiaLesion_df <- process_and_filter("pycnidiaPerCm2Lesion", test)


cor_PLACL <- correlation_dfs(blues_test, "PLACL", PLACL_df)
cor_pycnidiaPerCm2Leaf <- correlation_dfs(blues_test, "pycnidiaPerCm2Leaf", pycnidiaLeaf_df)
cor_pycnidiaPerCm2Lesion <- correlation_dfs(blues_test, "pycnidiaPerCm2Lesion", pycnidiaLesion_df)

# based on the AIC comparison, wqe will select the model 2 (including BRep) and the phenotype1 (inclusing all leaves)
cor_traits <- list(cor_PLACL, cor_pycnidiaPerCm2Leaf, cor_pycnidiaPerCm2Lesion)
pheno1_model2 <- unlist(map(cor_traits, \(x) x[2,2]))

prediction_df <- data.frame(Traits =  c("PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion"),
                            Correlation = pheno1_model2)

write_csv(prediction_df, file = 'outputs/plots/prediction.csv')


