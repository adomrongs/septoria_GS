library(tidyverse)
library(BGLR)
library(here)
library(lme4)
library(rrBLUP)
source("code/R/function_septoria_GS.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

load("data/modified_data/4_CV_data.Rdata")

trait <- "PLACL"
Kw <- k_wheat
Kmix <- k_mixes
pheno <- phenotype
genoW <- genotype_wheat
map <- map_wheat
wtest <- sample(rownames(Kw), ceiling(0.2 * nrow(Kw)))
wModel <- TRUE
sMix <- unique(phenotype$Strain)[1]
formula <- "~ -1 + Plant + Strain + Rep + Leaf "

runS1 <- function(trait, Kw, Kmix, pheno, genoW, map, wtest, formula, wModel = FALSE){
  #===============================================
  # Create data for train and test
  # ==============================================
  wlines <- rownames(Kw)
  wtrain <- setdiff(wlines, wtest)
  
  ptrain <- pheno
  ptrain[ptrain$Plant %in% wtest, trait] <- NA
  #===============================================
  # Run GWAS on train set
  # ==============================================
  bonferroni <- 0.05/nrow(map)
  # prepare train data
  gtrain <- data.frame(genoW[genoW[,1] %in% wtrain, ]) # subset genotype
  ptrain_GWAS <- pheno[pheno$Plant %in% gtrain[,1], ] # subset pheno
  Ktrain <- data.frame(A.mat(gtrain[,-1])) 
  colnames(Ktrain) <- gtrain[,1]
  Ktrain <- Ktrain %>% #subset K
    mutate(GenoID = gtrain[,1]) %>% 
    dplyr::select(GenoID, everything())
  
  BLUEs <- extract_blues_df(ptrain_GWAS, "PLACL", formula) # obtain BLUEs
  gtrain_GWAS <- gtrain %>% # adjust genotype to BLUEs lines
    filter(GenoID %in% BLUEs$GenoID)
  Ktrain_GWAS <- Ktrain %>% # adjust kinship to BLUEs lines
    filter(GenoID %in% BLUEs$GenoID) %>% 
    dplyr::select(c(GenoID, which(colnames(Ktrain) %in% BLUEs$GenoID)))
  
  dim(BLUEs); dim(gtrain_GWAS); dim(map); dim(Ktrain_GWAS)
  # run GWAS
  tmp <- capture.output({
    scores <- GAPIT(Y = BLUEs,
                    GD = gtrain_GWAS,
                    GM = map,
                    KI = Ktrain_GWAS,
                    CV = NULL,
                    PCA.total = 2,
                    model = "Blink",
                    file.output = F)
  })
  rm(tmp)
  gc()
  setwd(here())
  
  results <- scores[["GWAS"]] %>% 
    arrange(P.value)
  
  #===============================================
  # Run predictions with/withouth markers
  # ==============================================
  
  Zwtrain <- model.matrix(~0 + Plant, data = ptrain)
  Zmixtrain <- model.matrix(~0 + Strain, data = ptrain)
  if (!wModel) {
    Xwtrain <- model.matrix(~ Rep + Leaf, data = ptrain)
  } 
  if (wModel) {
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
    sSNPs_data <- genoW[, sSNPs, drop = FALSE] %>%
      data.frame() %>% 
      mutate(ID = genoW[,1])
    
    Xwtrain <- model.matrix(~ Rep + Leaf, data = ptrain) %>%
      data.frame() %>% 
      mutate(ID = ptrain$Plant) %>% 
      left_join(sSNPs_data, by = "ID") %>% 
      dplyr::select(-ID) %>% 
      as.matrix()
  }
  
  K12_wheat <- Zwtrain %*% as.matrix(Kw) %*% t(Zwtrain)
  K12_mix <- Zmixtrain %*% Kmix %*% t(Zmixtrain)
  K12_combined <- K12_mix * K12_wheat
  
  model1 <- BGLR(y = as.numeric(ptrain[[trait]]), ETA = list(
    list(X = Xwtrain, model = "FIXED"),
    list(K = K12_wheat, model = "RKHS"),
    list(K = K12_mix, model = "RKHS")
  ), nIter = 10000, burnIn = 2000, verbose = FALSE)
  
  model1_I <- BGLR(y = as.numeric(ptrain[[trait]]), ETA = list(
    list(X = Xwtrain, model = "FIXED"),
    list(K = K12_wheat, model = "RKHS"),
    list(K = K12_mix, model = "RKHS"),
    list(K = K12_combined, model = "RKHS")
  ), nIter = 10000, burnIn = 2000, verbose = FALSE)
  
  return(list(
    S1 = model1, # results model without interaction term
    S1_I = model1_I, # results with interaction term
    wtest = wtest, # partition test 
    wtrain = wtrain,
    hits = results
  ))
}

eval_S1 <- function(strategy, phenotype, trait) {
  predictions <- predict(strategy$S1)
  predictions_I <- predict(strategy$S1_I)
  
  lines_test <- phenotype$Plant %in% strategy$wtest
  
  pheno_test <- phenotype[lines_test,]
  
  cor <- cor(predictions[lines_test],
             pheno_test[[trait]],
             use = "complete.obs")
  cor_I <- cor(predictions_I[lines_test],
               pheno_test[[trait]],
               use = "complete.obs")
  
  withinStrainCorrelations <- function(predictions, phenotypeData, testLocations) {
    data <- data.frame(predictions = predictions[testLocations],
                       traits = phenotypeData[[trait]][testLocations],
                       Strain = phenotypeData$Strain[testLocations])
    
    listData <- split(data, data$Strain)
    lapply(listData, function(df) cor(df$predictions, df$traits, use = "complete.obs"))
  }
  
  cor_withinStrain <- withinStrainCorrelations(predictions, phenotype, lines_test)
  cor_withinStrain_I <- withinStrainCorrelations(predictions_I, phenotype, lines_test)
  
  cor_results <- list(
    cor = cor,
    cor_I = cor_I,
    cor_withinStrain = cor_withinStrain,
    cor_withinStrain_I = cor_withinStrain_I,
    hits <- strategy$hits
  )
  
  return(cor_results)
}

runS2 <- function(trait, Kw, Kmix, phenotype, genoW, map, sMix, formula, wModel = FALSE) {
  #===============================================
  # Run GWAS 
  # ==============================================
  bonferroni <- 0.05/nrow(map)
  Kw <- data.frame(Kw) %>% 
    rownames_to_column("GenoID") %>% 
    dplyr::select(GenoID, everything())
    
  BLUEs <- extract_blues_df(phenotype, "PLACL", formula) 
  
  dim(BLUEs); dim(genoW); dim(map); dim(Kw)
  tmp <- capture.output({
    scores <- GAPIT(Y = BLUEs,
                    GD = genoW,
                    GM = map,
                    KI = Kw,
                    CV = NULL,
                    PCA.total = 2,
                    model = "Blink",
                    file.output = F)
  })
  rm(tmp)
  gc()
  setwd(here())
  
  results <- scores[["GWAS"]] %>% 
    arrange(P.value)
  
  #===============================================
  # Run predictions with/withouth markers
  # ==============================================
  ptrain <- phenotype
  ptrain[ptrain$Strain %in% sMix, trait] <- NA
  
  if (!wModel) {
    Xwtrain <- model.matrix(~ Rep + Leaf, data = ptrain)
  } 
  if (wModel) {
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
    sSNPs_data <- genoW[, sSNPs, drop = FALSE] %>%
      data.frame() %>% 
      mutate(ID = genoW[,1])
    
    Xwtrain <- model.matrix(~ Rep + Leaf, data = ptrain) %>%
      data.frame() %>% 
      mutate(ID = ptrain$Plant) %>% 
      left_join(sSNPs_data, by = "ID") %>% 
      dplyr::select(-ID) %>% 
      as.matrix()
  }
  
  Zwtrain <- model.matrix(~0 + Plant, data = ptrain)
  Zmixtrain <- model.matrix(~0 + Strain, data = ptrain)
  
  K12_mix <- Zmixtrain %*% as.matrix(Kmix) %*% t(Zmixtrain)
  K12_wheat <- Zwtrain %*% as.matrix(Kw[,-1]) %*% t(Zwtrain)
  K12_combined <- K12_mix * K12_wheat
  
  model2 <- BGLR(y = as.numeric(ptrain[[trait]]), ETA = list(
    list(X = Xwtrain, model = "FIXED"),
    list(K = K12_wheat, model = "RKHS"),
    list(K = K12_mix, model = "RKHS")
  ), nIter = 10000, burnIn = 2000, verbose = FALSE)
  
  model2_I <- BGLR(y = as.numeric(ptrain[[trait]]), ETA = list(
    list(X = Xwtrain, model = "FIXED"),
    list(K = K12_mix, model = "RKHS"),
    list(K = K12_wheat, model = "RKHS"),
    list(K = K12_combined, model = "RKHS")
  ), nIter = 10000, burnIn = 2000, verbose = FALSE)
  
  return(list(model2 = model2,
              model2_I = model2_I,
              sMix = sMix))
}

eval_S2 <- function(strategy, phenotype, trait) {
  predictions <- predict(strategy$model2)
  predictions_I <- predict(strategy$model2_I)
  
  mix_test <- phenotype$Strain %in% strategy$sMix
  
  ptest <- phenotype[mix_test,]
  ptestTrait <- ptest[[trait]]
  
  cor <- cor(predictions[mix_test], ptestTrait)
  cor_I <- cor(predictions_I[mix_test], ptestTrait)
  
  correlationResults <- list(
    cor = cor, 
    corInteraction = cor_I 
  )
  
  return(list(CorrelationResults = correlationResults))
}

run_S3 <- function(trait, Kw, Kmix, phenotype, genoW, map, sMix, formula, wtest, wModel = FALSE) {
  #===============================================
  # Create data for train and test
  # ==============================================
  wlines <- rownames(Kw)
  wtrain <- setdiff(wlines, wtest)
  
  #===============================================
  # Run GWAS on train set
  # ==============================================
  bonferroni <- 0.05/nrow(map)
  gtrain <- data.frame(genoW[genoW[,1] %in% wtrain, ])
  
  ptrain_GWAS <- phenotype[phenotype$Plant %in% gtrain[,1], ]
  Ktrain <- data.frame(A.mat(gtrain[,-1]))
  colnames(Ktrain) <- rownames(Ktrain) <- gtrain[,1]
  Ktrain <- data.frame(Ktrain) %>% 
    rownames_to_column("GenoID") %>% 
    dplyr::select(GenoID, everything())
  rownames(Ktrain) <- NULL
    
  BLUEs <- extract_blues_df(ptrain_GWAS, "PLACL", formula) 
  gtrain_GWAS <- gtrain %>% # adjust genotype to BLUEs lines
    filter(GenoID %in% BLUEs$GenoID)
  Ktrain_GWAS <- Ktrain %>% # adjust kinship to BLUEs lines
    filter(GenoID %in% BLUEs$GenoID) %>% 
    dplyr::select(c(GenoID, which(colnames(Ktrain) %in% BLUEs$GenoID)))
  
  dim(BLUEs); dim(gtrain_GWAS); dim(map); dim(Ktrain_GWAS)
  tmp <- capture.output({
    scores <- GAPIT(Y = BLUEs,
                    GD = gtrain_GWAS,
                    GM = map,
                    KI = Ktrain_GWAS,
                    CV = NULL,
                    PCA.total = 2,
                    model = "Blink",
                    file.output = F)
  })
  rm(tmp)
  gc()
  setwd(here())
  
  results <- scores[["GWAS"]] %>% 
    arrange(P.value)
  #===============================================
  # Run predictions with/withouth markers
  # ==============================================
  
  ptrain <- phenotype
  ptrain[ptrain$Strain %in% sMix, trait] <- NA
  ptrain[ptrain$Plant %in% wtest, trait] <- NA
  
  if (!wModel) {
    Xwtrain <- model.matrix(~ Rep + Leaf, data = ptrain)
  } 
  if (wModel) {
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
    sSNPs_data <- genoW[, sSNPs, drop = FALSE] %>%
      data.frame() %>% 
      mutate(ID = genoW[,1])
    
    Xwtrain <- model.matrix(~ Rep + Leaf, data = ptrain) %>%
      data.frame() %>% 
      mutate(ID = ptrain$Plant) %>% 
      left_join(sSNPs_data, by = "ID") %>% 
      dplyr::select(-ID) %>% 
      as.matrix()
  }
  
  Zwtrain <- model.matrix(~0 + Plant, data = ptrain)
  Zmixtrain <- model.matrix(~0 + Strain, data = ptrain)
  
  K12_mix <- Zmixtrain %*% as.matrix(Kmix) %*% t(Zmixtrain)
  K12_wheat <- Zwtrain %*% as.matrix(Kw) %*% t(Zwtrain)
  K12_combined <- K12_mix * K12_wheat
  
  model3 <- BGLR(y = as.numeric(ptrain[[trait]]), ETA = list(
    list(X = Xwtrain, model = "FIXED"),
    list(K = K12_wheat, model = "RKHS"),
    list(K = K12_mix, model = "RKHS")
  ), nIter = 10000, burnIn = 2000, verbose = FALSE)
  
  model3_I <- BGLR(y = as.numeric(ptrain[[trait]]), ETA = list(
    list(X = Xwtrain, model = "FIXED"),
    list(K = K12_mix, model = "RKHS"),
    list(K = K12_wheat, model = "RKHS"),
    list(K = K12_combined, model = "RKHS")
  ), nIter = 10000, burnIn = 2000, verbose = FALSE)
  
  return(list(model3 = model3,
              model3_I = model3_I,
              sMix = sMix,
              wtest = wtest))
}

eval_S3 <- function(strategy, phenotype, trait) {
  predictions <- predict(strategy$model3)
  predictions_I <- predict(strategy$model3_I)
  
  mix_test <- phenotype$Strain %in% strategy$sMix
  wheat_test <- phenotype$Plant %in% strategy$wtest
  
  all_test <- mix_test & wheat_test
  
  ptest <- phenotype[all_test,]
  ptestTrait <- ptest[[trait]]
  
  cor <- cor(predictions[all_test], ptestTrait)
  cor_I <- cor(predictions_I[all_test], ptestTrait)
  
  correlationResults <- list(
    cor = cor, 
    corInteraction = cor_I
  )
  
  return(list(CorrelationResults = correlationResults))
}