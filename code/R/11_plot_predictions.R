library(tidyverse)
library(lme4)
source("code/R/function_septoria_GS.R")

load("data/modified_data/5_predictions.Rdata")
test <- unique(clean_test_pheno$Isolate)

# extract the blues for the test lines 
formula <- "~ -1 + strain + dpi + cultivar + replicate"
blues_test <- extract_blues_df_isolates(clean_test_pheno,
                                   traits = colnames(clean_test_pheno)[2:4],
                                   formula = formula)

files <- list.files("data/modified_data/predictions/", full.names = T)
blups_list <- list()
for(i in seq_along(files)){
  rdata <- files[[i]]
  load(rdata)
  blups_list[[i]] <- blups_df
}

blups_df <- do.call(cbind, blups_list)
blups_df <- blups_df[, c(1, 2, 5, 8)]
colnames(blups_df) <- c("Isolate", "PLACL", "pycnidiaPerCm2Leaf", "pycnidiaPerCm2Lesion")
rownames(blups_df) <- NULL

blups_test <- blups_df %>% 
  filter(Isolate %in% blues_test$Isolate)

blups_test <- blups_test %>% arrange(Isolate)
blues_test <- blues_test %>% arrange(Isolate)

corCalculation <- function(df1, df2) {
  # Asegurarse de que ambos data frames tienen el mismo número de columnas
  if (ncol(df1) != ncol(df2)) {
    stop("Los data frames deben tener el mismo número de columnas.")
  }
  # Calcular la correlación para cada par de columnas (excluyendo la primera columna)
  correlations <- map2_dbl(df1[, -1], df2[, -1], ~ cor(.x, .y))
  # Crear un data frame con los resultados
  results <- data.frame(
    Column = colnames(df1)[-1],  # Nombres de las columnas correlacionadas
    Correlation = correlations
  )
  return(results)
}

correlations <- corCalculation(blups_test, blues_test)

