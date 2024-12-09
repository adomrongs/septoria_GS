library(tidyverse)
library(sommer)
source("code/R/function_septoria_GS.R")

load("data/modified_data/5_predictions.Rdata")
test <- unique(clean_test_pheno$Isolate)

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

blups <- model$U$`u:Isolate`
blups_df <- data.frame(do.call(cbind, blups)) %>% 
  rownames_to_column("Isolate") %>%
  dplyr::select(Isolate, everything())

save(blues, blups_df, file = "data/modified_data/6_pred_results.Rdata")


