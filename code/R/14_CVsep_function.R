library(tidyverse)
library(sommer)
library(here)
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("code/R/function_septoria_GS.R")

load("data/modified_data/6_CV_septoria.Rdata")

genotype <- genotype
phenotype <- phenotype 
kinship <- kinship
map <- map_septoria
test <- sample(rownames(kinship), ceiling(0.2 * nrow(kinship)))
trait <- "PLACL"
wModel <- T


cv_septoria <- function(genotype, phenotype, kinship, map, test, trait, blues_all,  wModel = FALSE){
  #===============================================
  # Create data for train and test
  # ==============================================
  train <- setdiff(rownames(kinship), test)
  #===============================================
  # Run GWAS on train set
  # ==============================================
  bonferroni <- 0.05/nrow(map)
  ptrain <- phenotype %>% filter(Isolate %in% train)
  blues <- blues_all %>% 
    filter(Isolate %in% train) %>% 
    dplyr::select(Isolate, trait)
  gtrain <- genotype %>% filter(genotype[,1] %in% blues$Isolate)
  ktrain <- kinship %>%
    filter(rownames(.) %in% blues$Isolate) %>%
    dplyr::select(all_of(blues$Isolate)) %>%
    rownames_to_column("Isolate") %>%
    dplyr::select(Isolate, everything())
  
  message("Data Ready")
  dim(blues); dim(gtrain); dim(map); dim(ktrain)
  #===============================================
  # Run predictions with/without markers
  # ==============================================
  
  if (!wModel) {
    formula_blups <- as.formula(paste(trait, "~ Line + Year + Trial + BRep"))
    results <- data.frame()
    message("Results created")
  } 
  if (wModel) {
    message("Running GWAS")
    tmp <- capture.output({
      scores <- GAPIT(Y = blues,
                      GD = gtrain,
                      GM = map,
                      KI = NULL,
                      CV = NULL,
                      PCA.total = 3,
                      model = "Blink",
                      file.output = F)
    })
    rm(tmp)
    gc()
    setwd(here())
    message("GWAS completed")
    
    results <- scores[["GWAS"]] %>% 
      arrange(P.value)
    message("Results created")
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
      paste0(trait, " ~ Line + Year + Trial + BRep + ", 
             paste0("`", sSNPs, "`", collapse = " + "))
    )
  }
  
  model <- mmer(formula_blups,
                random = ~ vsr(Isolate, Gu = as.matrix(kinship)),
                rcov = ~ units,
                data = ptrain,
                verbose = TRUE)
  message("model correctly created")
  
  H2 <- h2_sommer(model, n = 12)
  message("hertability calculated")
  
  blups_test <- data.frame(Isolate = names(model$U[[1]][[1]])) %>%
    mutate(!!trait := model$U[[1]][[1]]) %>% 
    filter(Isolate %in% test) %>% 
    arrange(Isolate)
  message("blups extracted")
  
  blues_test <- blues_all %>% 
    filter(Isolate %in% test) %>% 
    arrange(Isolate)
  message("blues extracted")
  
  ability <- cor(blups_test[,trait], blues_test[,trait])
  accuracy <- ability/sqrt(H2)
  message("results extracted")
  
  results <- list(ability = ability, accuracy = accuracy, results = results)
  return(results)
  
}
