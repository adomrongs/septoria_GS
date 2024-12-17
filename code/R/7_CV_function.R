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
pheno <- phenotype <- adjusted_phenotype
genoW <- genotype_wheat
map <- map_wheat
wtest <- sample(rownames(Kw), ceiling(0.2 * nrow(Kw)))
wModel <- TRUE
sMix <- unique(adjusted_phenotype$Strain)[1]

runS1 <- function(trait, Kw, Kmix, pheno, genoW, map, wtest, formula, wModel = FALSE){
  #===============================================
  # Create data for train and test
  # ==============================================
  wlines <- rownames(Kw)
  wtrain <- setdiff(wlines, wtest)
  
  ptrain <- pheno
  ptrain[ptrain$Plant %in% wtest, trait] <- NA
  
  #===============================================
  # Run predictions with/withouth markers
  # ==============================================
  
  Zwtrain <- model.matrix(~0 + Plant, data = ptrain)
  Zmixtrain <- model.matrix(~0 + Strain, data = ptrain)
  if (!wModel) {
    Xwtrain <- model.matrix(~ 1, data = ptrain)
  } 
  if (wModel) {
    bonferroni <- 0.05/nrow(map)
    # prepare train data
    gwas_geno <- data.frame(genoW[genoW[,1] %in% wtrain, ]) # subset genotype
    gwas_pheno <- pheno[pheno$Plant %in% gwas_geno[,1], c("Plant", trait)] # subset pheno
    gwas_k <- Kw %>% 
      filter(rownames(.) %in% gwas_geno[,1]) %>% 
      dplyr::select(all_of(gwas_geno[,1])) %>% 
      rownames_to_column("GenoID") %>% 
      dplyr::select(GenoID, everything())
    
    dim(gwas_geno); dim(gwas_pheno); dim(map); dim(gwas_k)
    # run GWAS
    tmp <- capture.output({
      scores <- GAPIT(Y = gwas_pheno,
                      GD = gwas_geno,
                      GM = map,
                      KI = gwas_k,
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
    
    Xwtrain <- model.matrix(~ 1, data = ptrain) %>%
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
  
  H2 <- computeH2(model1)

  model1_I <- BGLR(y = as.numeric(ptrain[[trait]]), ETA = list(
    list(X = Xwtrain, model = "FIXED"),
    list(K = K12_wheat, model = "RKHS"),
    list(K = K12_mix, model = "RKHS"),
    list(K = K12_combined, model = "RKHS")
  ), nIter = 10000, burnIn = 2000, verbose = FALSE)
  
  H2_I <- computeH2(model1_I, interaction = T)
  
  return(list(
    S1 = model1, # results model without interaction term
    S1_I = model1_I, # results with interaction term
    wtest = wtest, # partition test 
    wtrain = wtrain,
    hits = results,
    H_list = list(H2, H2_I)
  ))
}

eval_S1 <- function(strategy, phenotype, trait) {
  predictions <- predict(strategy$S1)
  predictions_I <- predict(strategy$S1_I)
  
  total_results <- data.frame(Plant = phenotype$Plant, 
                      Strain = phenotype$Strain,
                      Observed = phenotype[[trait]], 
                      Predicted = predictions,
                      Predicted_I = predictions_I)
  test_results <- total_results %>% 
    filter(Plant %in% strategy$wtest)
  # overall correlation
  lines_test <- phenotype$Plant %in% strategy$wtest
  pheno_test <- phenotype[lines_test,]
  cor <- cor(predictions[lines_test],
             pheno_test[[trait]],
             use = "complete.obs")
  cor_I <- cor(predictions_I[lines_test],
               pheno_test[[trait]],
               use = "complete.obs")
  
  # within strain correlation for both model and model_I
  cor_withinStrain <- test_results %>%
    split(.$Strain) %>%  # Divide el dataframe por cada valor único de 'Strain'
    map(~ cor(.x$Observed, .x$Predicted, use = "complete.obs"))
  cor_withinStrain_I <- test_results %>%
    split(.$Strain) %>%  # Divide el dataframe por cada valor único de 'Strain'
    map(~ cor(.x$Observed, .x$Predicted_I, use = "complete.obs"))
  # prediction accuracy calculation as the prediction ability divided by srt of h2
  accuracy_withinStrain <- map(cor_withinStrain, ~ .x/sqrt(strategy$H_list[[1]]))
  accuracy_withinStrain_I <- map(cor_withinStrain, ~ .x/sqrt(strategy$H_list[[2]]))
  # save resulst
  cor_results <- list(
    cor = cor,
    cor_I = cor_I,
    cor_withinStrain = cor_withinStrain,
    cor_withinStrain_I = cor_withinStrain_I,
    hits <- strategy$hits,
    accuracy = accuracy_withinStrain,
    accuracy_I = accuracy_withinStrain_I
  )
  
  return(cor_results)
}

runS2 <- function(trait, Kw, Kmix, phenotype, genoW, map, sMix, formula, wModel = FALSE) {
#===============================================
  # Run predictions with/withouth markers
  # ==============================================
  ptrain <- phenotype
  ptrain[ptrain$Strain %in% sMix, trait] <- NA
  
  if (!wModel) {
    Xwtrain <- model.matrix(~ 1, data = ptrain)
  } 
  if (wModel) {
    bonferroni <- 0.05/nrow(map)
    gwas_geno <- genoW
    gwas_pheno <- phenotype %>% 
      dplyr::select(Plant, trait)
    gwas_k <- data.frame(Kw) %>% 
      rownames_to_column("GenoID") %>% 
      dplyr::select(GenoID, everything())
    
    dim(gwas_pheno); dim(gwas_geno); dim(map); dim(gwas_k)
    tmp <- capture.output({
      scores <- GAPIT(Y = gwas_pheno,
                      GD = gwas_geno,
                      GM = map,
                      KI = gwas_k,
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
    
    Xwtrain <- model.matrix(~ 1, data = ptrain) %>%
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
  
  model2 <- BGLR(y = as.numeric(ptrain[[trait]]), ETA = list(
    list(X = Xwtrain, model = "FIXED"),
    list(K = K12_wheat, model = "RKHS"),
    list(K = K12_mix, model = "RKHS")
  ), nIter = 10000, burnIn = 2000, verbose = FALSE)
  
  H2 <- computeH2(model2)
  
  model2_I <- BGLR(y = as.numeric(ptrain[[trait]]), ETA = list(
    list(X = Xwtrain, model = "FIXED"),
    list(K = K12_mix, model = "RKHS"),
    list(K = K12_wheat, model = "RKHS"),
    list(K = K12_combined, model = "RKHS")
  ), nIter = 10000, burnIn = 2000, verbose = FALSE)
  
  H2_I <- computeH2(model2_I, interaction = T)
  
  return(list(model2 = model2,
              model2_I = model2_I,
              sMix = sMix,
              H_list = list(H2, H2_I)))
}

eval_S2 <- function(strategy, phenotype, trait) {
  predictions <- predict(strategy$model2)
  predictions_I <- predict(strategy$model2_I)
  
  total_results <- data.frame(Plant = phenotype$Plant, 
                              Strain = phenotype$Strain, 
                              Observed = phenotype[[trait]], 
                              Predicted = predictions, 
                              Predicted_I = predictions_I)
  test_results <- total_results %>% 
    filter(Strain %in% strategy$sMix)
  
  cor <- cor(test_results$Observed, test_results$Predicted)
  cor_I <- cor(test_results$Observed, test_results$Predicted_I)
  
  accuracy <- cor/sqrt(strategy$H_list[[1]])
  accuracy_I <- cor_I/sqrt(strategy$H_list[[2]])
  
  correlationResults <- list(
    cor = cor, 
    corInteraction = cor_I,
    accuracy = accuracy,
    accuracy_I = accuracy_I
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
  # Run predictions with/withouth markers
  # ==============================================
  ptrain <- phenotype
  ptrain[ptrain$Strain %in% sMix, trait] <- NA
  ptrain[ptrain$Plant %in% wtest, trait] <- NA
  
  if (!wModel) {
    Xwtrain <- model.matrix(~ 1, data = ptrain)
  } 
  if (wModel) {
    bonferroni <- 0.05/nrow(map)
    gwas_pheno <- phenotype %>% 
      filter(Plant %in% wtrain, 
             Strain %in% sMix) %>% 
      dplyr::select(Plant, trait)
    gwas_geno <- data.frame(genoW[genoW[,1] %in% gwas_pheno$Plant, ])
    gwas_k <- Kw %>% 
      filter(rownames(.) %in% gwas_geno[,1]) %>% 
      dplyr::select(all_of(gwas_geno[,1])) %>% 
      rownames_to_column("GenoID") %>% 
      dplyr::select(GenoID, everything())
    
    dim(gwas_geno); dim(gwas_pheno); dim(gwas_k); dim(map)
    tmp <- capture.output({
      scores <- GAPIT(Y = gwas_pheno,
                      GD = gwas_geno,
                      GM = map,
                      KI = gwas_k,
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
    
    Xwtrain <- model.matrix(~ 1, data = ptrain) %>%
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
  
  H2 <- computeH2(model3)
  
  model3_I <- BGLR(y = as.numeric(ptrain[[trait]]), ETA = list(
    list(X = Xwtrain, model = "FIXED"),
    list(K = K12_mix, model = "RKHS"),
    list(K = K12_wheat, model = "RKHS"),
    list(K = K12_combined, model = "RKHS")
  ), nIter = 10000, burnIn = 2000, verbose = FALSE)
  
  H2_I <- computeH2(model3_I, interaction = T)
  
  return(list(model3 = model3,
              model3_I = model3_I,
              sMix = sMix,
              wtest = wtest,
              H_list = list(H2, H2_I)))
}

eval_S3 <- function(strategy, phenotype, trait) {
  predictions <- predict(strategy$model3)
  predictions_I <- predict(strategy$model3_I)
  
  total_results <- data.frame(Plant = phenotype$Plant, 
                              Strain = phenotype$Strain, 
                              Observed = phenotype[[trait]], 
                              Predicted = predictions, 
                              Predicted_I = predictions_I)
  test_results <- total_results %>% 
    filter(Plant %in% strategy$wtest, 
          Strain %in% strategy$sMix)

  cor <- cor(test_results$Observed, test_results$Predicted)
  cor_I <- cor(test_results$Observed, test_results$Predicted)
  
  accuracy <- cor/sqrt(strategy$H_list[[1]])
  accuracy_I <- cor_I/sqrt(strategy$H_list[[2]])

  
  correlationResults <- list(
    cor = cor, 
    corInteraction = cor_I,
    accuracy = accuracy,
    accuracy_I = accuracy_I
  )
  
  return(list(CorrelationResults = correlationResults))
}