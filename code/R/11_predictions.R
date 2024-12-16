library(tidyverse)
library(sommer)
source("code/R/function_septoria_GS.R")

# Load necessary data
load("data/modified_data/5_predictions.Rdata")

# Retrieve task ID
i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Ensure task ID is valid
traits <- c("PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion")

# Select trait for this run
trait <- traits[i]

phenotypes <- list(cleaned_septoria_phenotype_1, cleaned_septoria_phenotype_2)
# Define formula
formula_1 <- as.formula(paste0(trait, " ~ Line + Year + Trial + Leaf"))
formula_2 <- as.formula(paste0(trait, " ~ Line + Year + Trial"))
formulas <- list(formula_1, formula_2)

for(j in seq_along(phenotypes)){
  phenotype <- phenotypes[[j]]
  formula <- formulas[[j]]
  
  sub_f1 <- formula
  sub_f2 <- update(sub_f1, . ~ . + BRep)
  sub_formulas <- list(sub_f1, sub_f2)
  
  for(z in seq_along(sub_formulas)){
    model_formula <- sub_formulas[[z]]
    model <- mmer(
      model_formula,
      random = ~ vsr(Isolate, Gu = k_all),
      rcov = ~ units,
      data = phenotype,
      verbose = TRUE
    )
    
    new_output <- paste0("data/modified_data/predictions/phenotype_", j,"/model_", z, "/", trait, ".Rdata")
    new_dir <- paste0("data/modified_data/predictions/phenotype_", j,"/model_", z)
    dir.create(new_dir, recursive = T)
    
    blups_df <- data.frame(Isolate = names(model$U$`u:Isolate`[[trait]]), 
                           BLUP = model$U$`u:Isolate`[[trait]]) %>% 
      mutate(Trait = trait)
    save(blups_df, file = new_output)
  }
}



