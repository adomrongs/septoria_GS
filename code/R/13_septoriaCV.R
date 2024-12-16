library(tidyverse)
library(sommer)
library(here())
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("code/R/function_septoria_GS.R")

load("data/modified_data/5_predictions.Rdata")
load("data/modified_data/1_septoria.Rdata")

genotype <- genotype_septoria
phenotype <- cleaned_septoria_phenotype
kinship <- k_septoria[,-1]
colnames(kinship) <- rownames(kinship) <- k_septoria[,1]
map <- map_septoria
trait <- "PLACL"
wModel <- F
test <- sample(rownames(kinship), ceiling(0.2 * nrow(kinship)))
formula <- "~ Isolate + Line + Trial + Year + BRep"

blues_all <- extract_blues_df_adapted(cleaned_septoria_phenotype,
                                   trait = c("PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion"),
                                   formula = formula, 
                                   colname = "Isolate")

cv_septoria <- function(genotype, phenotype, kinship, map, test, trait, blues_all,  wModel = FALSE){
  #===============================================
  # Create data for train and test
  # ==============================================
  train <- setdiff(rownames(kinship), test)
  #===============================================
  # Run GWAS on train set
  # ==============================================
  bonferroni <- 0.05/nrow(map)
  ptrain <- phenotype %>% filter(Isolate %in% train, )
  blues <- extract_blues_adapted(ptrain, trait, formula, "Isolate")
  gtrain <- genotype %>% filter(genotype[,1] %in% blues$Isolate)
  ktrain <- kinship %>%
    filter(rownames(.) %in% blues$Isolate) %>%
    dplyr::select(all_of(blues$Isolate)) %>%
    rownames_to_column("Isolate") %>%
    dplyr::select(Isolate, everything())
  
  dim(blues); dim(gtrain); dim(map); dim(ktrain)
  #===============================================
  # Run predictions with/withouth markers
  # ==============================================
  
  if (!wModel) {
    formula_blups <- as.formula(paste(trait, "~ Line + Year + Trial + Leaf + BRep"))
  } 
  if (wModel) {
    tmp <- capture.output({
      scores <- GAPIT(Y = blues,
                      GD = gtrain,
                      GM = map,
                      KI = ktrain,
                      CV = NULL,
                      PCA.total = 3,
                      model = "Blink",
                      file.output = F)
    })
    rm(tmp)
    gc()
    setwd(here())
    
    results <- scores[["GWAS"]] %>% 
      arrange(P.value)
    hits_bonferroni <- results %>% filter(P.value <= bonferroni)
    
    if (nrow(hits_bonferroni) == 0) {
      sSNPs <- results$SNP[1]
    }
    if (nrow(hits_bonferroni) >= 3) {
      sSNPs <- results$SNP[1:3]
    }
    if (nrow(hits_bonferroni) == 2) {
      sSNPs <- results$SNP[2]
    }
    if (nrow(hits_bonferroni) == 1) {
      sSNPs <- results$SNP[2]
    }
    sSNPs_data <- gtrain[, sSNPs, drop = FALSE] %>%
      mutate(Isolate = gtrain[,1])
    
    ptrain <- ptrain %>% 
      left_join(sSNPs_data)
    
    formula_blups <- as.formula(
      paste0(trait, "~ Line + Year + Trial + Leaf + BRe", paste(sSNPs, collapse = " + "))
    )
  }
  
  model <- mmer(formula_blups,
                random = ~ vsr(Isolate, Gu = as.matrix(kinship)),
                rcov = ~ units,
                data = ptrain,
                verbose = TRUE)
  H2 <- h2_sommer(model, n = 12)
  
  blups_test <- data.frame(Isolate = names(model$U[[1]][[1]])) %>%
    mutate(!!trait := model$U[[1]][[1]]) %>% 
    filter(Isolate %in% test) %>% 
    arrange(Isolate)
  
  blues_test <- blues_all %>% 
    filter(Isolate %in% test) %>% 
    arrange(Isolate)
  
  ability <- cor(blups_test[,2], blues_test[,2])
  accuracy <- ability/sqrt(H2)
  
  results <- list(ability = ability, accuracy = accuracy)
  return(results)

}