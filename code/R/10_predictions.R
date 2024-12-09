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
traits <- c("PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion")
blups_list <- list()

# Loop over traits to fit the model and extract BLUPs
for (i in seq_along(traits)) {
  trait <- traits[i]
  formula <- as.formula(paste0(trait, " ~ Line + Year + REP + Leaf"))
  
  model <- mmer(
    formula,
    random = ~ vsr(Isolate, Gu = k_all),
    rcov = ~ units,
    data = cleaned_septoria_phenotype,
    verbose = TRUE
  )
  
  # Store BLUPs in the list, naming by trait
  blups_list[[trait]] <- model$U$`u:Isolate`[[trait]]
}

# Combine BLUPs into a single data frame
blups_df <- data.frame(do.call(cbind, blups_list)) %>% 
  rownames_to_column("Isolate") %>% 
  dplyr::select(Isolate, everything())

save(blues, blups_df, file = "data/modified_data/6_pred_results.Rdata")


