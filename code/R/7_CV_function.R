library(tidyverse)
library(BGLR)
library(here)
library(lme4)
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


runS1 <- function(trait, Kw, Kmix, pheno, genoW, map, wtest, wModel = NULL){
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
  
  formula <- "~ Strain + Rep + Leaf + (1|Plant)"
  BLUPs <- extract_blups_df(ptrain_GWAS, "PLACL", formula) # obtain blups
  gtrain_GWAS <- gtrain %>% # adjust genotype to blups lines
    filter(GenoID %in% BLUPs$GenoID)
  Ktrain_GWAS <- Ktrain %>% # adjust kinship to blups lines
    filter(GenoID %in% BLUPs$GenoID) %>% 
    dplyr::select(c(GenoID, which(colnames(Ktrain) %in% BLUPs$GenoID)))
  
  dim(BLUPs); dim(gtrain_GWAS); dim(map); dim(Ktrain_GWAS)
  # run GWAS
  tmp <- capture.output({
    scores <- GAPIT(Y = BLUPs,
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
    Xwtrain <- model.matrix(~1, data = ptrain) 
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
    
    Xwtrain <- model.matrix(~1, data = ptrain) %>%
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