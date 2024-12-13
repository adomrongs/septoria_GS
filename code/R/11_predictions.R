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

# Define formula
formula <- as.formula(paste0(trait, " ~ Line + Year + Trial + Leaf + BRep"))

# Fit model
model <- mmer(
  formula,
  random = ~ vsr(Isolate, Gu = k_all),
  rcov = ~ units,
  data = cleaned_septoria_phenotype,
  verbose = TRUE
)

# Extract BLUPs and create data frame
blups_list <- list()
blups_list[[trait]] <- model$U$`u:Isolate`[[trait]]

blups_df <- data.frame(Isolate = names(blups_list[[trait]]), 
                       BLUP = blups_list[[trait]]) %>% 
  mutate(Trait = trait)

# Save results
dir.create("data/modified_data/predictions4/", recursive = TRUE, showWarnings = FALSE)
save(blups_df, file = paste0("data/modified_data/predictions4/pred_", trait, ".Rdata"))



