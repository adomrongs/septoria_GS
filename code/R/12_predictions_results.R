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

# extract the blues for the test lines 
formula <- "~ -1 + Isolate + Line + TRep + BRep"
blues_test <- extract_blues_df_adapted(clean_test_pheno,
                                   traits = c("PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion"),
                                   formula = formula,
                                   colname = "Isolate")
blues_test <- blues_test %>% 
  filter(Isolate %in% test_name) 

# load the calculated BLUPs calculated like this " trait ~ Line + Year + Trial + BRep + (1|Isolate)"
files <- list.files("data/modified_data/predictions/", full.names = T)
blups_list <- list()
for(i in seq_along(files)){
  rdata <- files[[i]]
  load(rdata)
  blups_list[[i]] <- blups_df
}
# bind the list to create a df containing the result for the 3 traits
blups_df <- do.call(cbind, blups_list)
blups_df <- blups_df[, c(1, 2, 5, 8)]
colnames(blups_df) <- c("Isolate", "PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion")
rownames(blups_df) <- NULL
#subset the blups for the test
blups_test <- blups_df %>% 
  filter(Isolate %in% test)
# sort them to fairly correlate each isolate
blups_test <- blups_test %>% arrange(Isolate)
blues_test <- blues_test %>% arrange(Isolate)

# calculate correlation
# this formula just compute the correlation between each column except for the first one for 2 dfs
correlations <- corCalculation(blups_test, blues_test)

# THIS IS THE SAME BUT FOR THE "OLD BLUPS" " trait ~ Line + Year + Trial + (1|Isolate)"
# predictions without BRep
files <- list.files("data/modified_data/predictions2/", full.names = T)
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

blups_test2 <- blups_test %>% arrange(Isolate)
blues_test <- blues_test %>% arrange(Isolate)

correlations2 <- corCalculation(blups_test2, blues_test)

# THIS IS THE SAME BUT FOR THE "OLD new BLUPS" " trait ~ Line + Year + Trial + (1|Isolate)"
# predictions without BRep
files <- list.files("data/modified_data/predictions3/", full.names = T)
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

blups_test2 <- blups_test %>% arrange(Isolate)
blues_test <- blues_test %>% arrange(Isolate)

correlations3 <- corCalculation(blups_test2, blues_test)

# THIS IS THE SAME BUT FOR THE "OLD new BLUPS" " trait ~ Line + Year + Trial + (1|Isolate)"
# predictions without BRep
files <- list.files("data/modified_data/predictions4/", full.names = T)
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

blups_test4 <- blups_test %>% arrange(Isolate)
blues_test <- blues_test %>% arrange(Isolate)

correlations4 <- corCalculation(blups_test4, blues_test)

