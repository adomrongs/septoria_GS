library(tidyverse)
library(sommer)
source("code/R/function_septoria_GS.R")

load("data/modified_data/5_predictions.Rdata")

traits <- c("PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion")
i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
trait <- traits[[i]]
# model1: REP 
formula1 <- as.formula(paste0(trait, " ~ -1 + Line + Year + Trial"))
model1 <- mmer(formula1, 
               random = ~ vsr(Isolate, Gu = k_all),
               rcov = ~ units,
               data = cleaned_septoria_phenotype,
               verbose = TRUE)

# model1: REP + BRep
formula2 <- as.formula(paste0(trait, " ~ -1 + Line + Year + Trial + BRep"))
model2 <- mmer(formula2, 
               random = ~ vsr(Isolate, Gu = k_all),
               rcov = ~ units,
               data = cleaned_septoria_phenotype,
               verbose = TRUE)

AIC_df <- data.frame(Model = c("Model1", "Model2"),
                     AIC = c(model1$AIC, model2$AIC),
                     Trait = trait)
results_dir <- "data/modified_data/AIC"
dir.create(results_dir, recursive = T) 
save(AIC_df, file = paste0(results_dir, "/",trait, ".Rdata"))
