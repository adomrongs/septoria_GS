library(tidyverse)
library(sommer)
source("code/R/function_septoria_GS.R")

load("data/modified_data/5_predictions.Rdata")

# extract the blues for the test lines 
formula <- "~ -1 + Isolate + Line + N"
blues <- extract_blues_df_isolates(clean_test_pheno,
                                   traits = colnames(clean_test_pheno)[2:4],
                                   formula = formula)


# First stage model
model <- mmer(PLACL ~ Line + Year + REP + Leaf,
              random = ~ vsr(Isolate, Gu = k_all),
              rcov = ~ units,
              data = cleaned_septoria_phenotype,
              verbose = TRUE)




