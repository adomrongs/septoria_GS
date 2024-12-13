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
wModel <- TRUE
test <- sample(rownames(kinship), ceiling(0.2 * nrow(kinship)))
formula <- "~ Isolate + Line + REP + Year"

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
  #===============================================
  # Run predictions with/withouth markers
  # ==============================================
  
  if (!wModel) {
    formula_blups <- as.formula(trait, "~ Line + REP + Year")
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
    sSNPs_data <- gtrain[, sSNPs, drop = FALSE] %>%
      mutate(Isolate = gtrain[,1])
    
    ptrain <- ptrain %>% 
      left_join(sSNPs_data)
    
    formula_blups <- as.formula(
      paste0(trait, "~ Line + REP + Year + ", paste(sSNPs, collapse = " + "))
    )
  }
  
  model <- mmer(formula_blups,
                random = ~ vsr(Isolate, Gu = as.matrix(kinship)),
                rcov = ~ units,
                data = ptrain,
                verbose = TRUE)
  H2 <- 
  
  blups <- data.frame(Isolate = names(model$U$`u:Isolate`$PLACL)) %>%
    mutate(!!trait := model$U$`u:Isolate`$PLACL)
  rownames(blups) <- NULL
  blups_test <- blups %>% 
    filter(Isolate %in% test) %>% 
    arrange(Isolate)
  
  blues_test <- blues_all %>% 
    filter(Isolate %in% test) %>% 
    arrange(Isolate)
  
  ability <- cor(blups_test[,2], blues_test[,2])
  accuracy <- ability/sqrt(H2)
  
  
  
  
}