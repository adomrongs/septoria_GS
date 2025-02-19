library(tidyverse)
library(BGLR)
library(here)
library(lme4)
library(rrBLUP)
source("code/R/function_septoria_GS.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

load("data/modified_data/4_CV_data.Rdata")

Kw_list <- list(G = k_wheat, Diagonal = I_wheat)
Ks_list <- list(G = k_mixes, Diagonal = I_mixes)

# Get iteration number from SLURM
iter <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Initialize lists to store results
allResults <- list()
wtest_lines <- list()  # Initialize wtest_lines before using it

# Select test lines for wheat
wheat_test <- sample(rownames(k_wheat), ceiling(0.2 * nrow(k_wheat)))
wtest_lines[[iter]] <- wheat_test
formula <- "~ -1 + Plant + Strain"

# Loop over wheat and mix kernels
for (wheatKey in names(Kw_list)) {
  for (mixKey in names(Ks_list)) {
    message(paste("Iteration", iter, "- Wheat:", wheatKey, "- Mix:", mixKey))
    
    # Create a list to store results for this iteration
    scenarioResults <- list()
    
    message("Running strategy 1")
    # Run Scenario 1
    ST1Models <- runS1(trait = "PLACL",
                       Kw = Kw_list[[wheatKey]],
                       Kmix = Ks_list[[mixKey]],
                       pheno = adjusted_phenotype,
                       genoW = genotype_wheat,
                       map = map_wheat,
                       wtest = wheat_test,
                       formula = formula, 
                       name_kw = wheatKey,
                       name_km = mixKey)
    message("Strategy 1 completed")
    message("Running evaluation 1")
    scenarioResults$Scenario1 <- eval_S1(strategy = ST1Models,
                                         phenotype = adjusted_phenotype,
                                         trait = "PLACL")
    message("Evaluation 1 completed")
    message("Running strategy 1 weighted")
    # Run Scenario 1 with weighting (wModel = TRUE)
    weightedModel <- runS1(trait = "PLACL",
                           Kw = Kw_list[[wheatKey]],
                           Kmix = Ks_list[[mixKey]],
                           pheno = adjusted_phenotype,
                           genoW = genotype_wheat,
                           map = map_wheat,
                           wtest = wheat_test,
                           wModel = TRUE,
                           formula = formula, 
                           name_kw = wheatKey,
                           name_km = mixKey)
    message("Strategy 1 weighted completed")
    message("Running evaluation 1 weighted")
    scenarioResults$Scenario1w <- eval_S1(strategy = weightedModel,
                                          phenotype = adjusted_phenotype,
                                          trait = "PLACL")
    message("Evaluation 1 weighted completed")
    
    for (mix in unique(adjusted_phenotype$Strain)) {
      ST2strategy <- runS2(trait = "PLACL",
                           Kw = Kw_list[[wheatKey]],
                           Kmix = Ks_list[[mixKey]],
                           phenotype = adjusted_phenotype,
                           genoW = genotype_wheat,
                           map = map_wheat,
                           sMix = mix, 
                           formula = formula)
      scenarioResults[[paste("Scenario2", mix)]] <- eval_S2(strategy = ST2strategy,
                                                            phenotype = adjusted_phenotype,
                                                            trait = "PLACL")
      
      weightedModel2 <- runS2(trait = "PLACL",
                              Kw = Kw_list[[wheatKey]],
                              Kmix = Ks_list[[mixKey]],
                              phenotype = adjusted_phenotype,
                              genoW = genotype_wheat,
                              map = map_wheat,
                              sMix = mix,
                              wModel = TRUE,
                              formula = formula)
      scenarioResults[[paste("Scenario2w", mix)]] <- eval_S2(strategy = weightedModel2,
                                                             phenotype = adjusted_phenotype,
                                                             trait = "PLACL")
    }
    
    for (mix in unique(adjusted_phenotype$Strain)) {
      ST3Models <- run_S3(trait = "PLACL",
                          Kw = Kw_list[[wheatKey]],
                          Kmix = Ks_list[[mixKey]],
                          phenotype = adjusted_phenotype, 
                          sMix = mix,
                          genoW = genotype_wheat,
                          map = map_wheat,
                          wtest = wheat_test, 
                          formula = formula)
      scenarioResults[[paste("Scenario3", mix)]] <- eval_S3(strategy = ST3Models,
                                                            phenotype = adjusted_phenotype,
                                                            trait = "PLACL")
      
      weightedModel3 <- run_S3(trait = "PLACL",
                               Kw = Kw_list[[wheatKey]],
                               Kmix = Ks_list[[mixKey]],
                               phenotype = adjusted_phenotype, 
                               sMix = mix,
                               genoW = genotype_wheat,
                               map = map_wheat,
                               wtest = wheat_test,
                               wModel = TRUE,
                               formula = formula)
      scenarioResults[[paste("Scenario3w", mix)]] <- eval_S3(strategy = weightedModel3,
                                                             phenotype = adjusted_phenotype,
                                                             trait = "PLACL")
    }
    
    # Store results in a structured format
    allResults[[paste("Iteration", iter, "Wheat", wheatKey, "Mix", mixKey)]] <- scenarioResults
  }
}

# Save allResults and wtest_lines for this iteration in a separate file
dir.create("data/modified_data/cv", recursive = T)
save(allResults, wtest_lines, file = paste0("data/modified_data/cv/iter_", iter, ".Rdata"))
