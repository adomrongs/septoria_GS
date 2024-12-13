library(tidyverse)
library(lme4)
source("code/R/function_septoria_GS.R")

load("data/modified_data/5_predictions.Rdata")

info_strains <- read_csv("data/raw_data/INFO_STRAINS.csv")
test <- info_strains %>% 
  filter(Partition == "Test") %>% 
  pull(Isolate)
test_name <- unlist(map(strsplit(test, "_"), ~.x[2]))

# extract the blues for the test lines 
formula <- "~ -1 + strain + dpi + cultivar + replicate"
blues_test <- extract_blues_df_adapted(clean_test_pheno,
                                   traits = c("PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion"),
                                   formula = formula,
                                   colname = "strain")
blues_test <- blues_test %>% 
  filter(strain %in% test_name) %>% 
  rename(Isolate = strain)

files <- list.files("data/modified_data/predictions/", full.names = T)
blups_list <- list()
for(i in seq_along(files)){
  rdata <- files[[i]]
  load(rdata)
  blups_list[[i]] <- blups_df
}

blups_df <- do.call(cbind, blups_list)
blups_df <- blups_df[, c(1, 2, 5, 8)]
colnames(blups_df) <- c("Isolate", "PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion")
rownames(blups_df) <- NULL

blups_test <- blups_df %>% 
  filter(Isolate %in% test)

blups_test <- blups_test %>% arrange(Isolate)
blues_test <- blues_test %>% arrange(Isolate)

correlations <- corCalculation(blups_test, blues_test)

