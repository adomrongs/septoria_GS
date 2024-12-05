library(tidyverse)
library(hrbrthemes)
library(ggplot2)
source("code/R/function_septoria_GS.R")

S1_df <- data.frame()
S2_df <- data.frame()
S3_df <- data.frame()
numIter <- 10
for (i in 1:numIter) {
  load(file = paste0("data/modified_data/cv/all/iter_", i, ".Rdata"))  # Cargar los datos
  
  # Aplicar las variaciones de la función con diferentes combinaciones de 'weighted' e 'interaction'
  S1 <- functS1(list = allResults, weighted = F, interaction = F)
  S1w <- functS1(list = allResults, weighted = T, interaction = F)
  S1_I <- functS1(list = allResults, weighted = F, interaction = T)
  S1w_I <- functS1(list = allResults, weighted = T, interaction = T)
  
  # Combinar los resultados de las cuatro variaciones
  S1_all <- rbind(S1, S1w, S1_I, S1w_I)
  S1_df <- rbind(S1_df, S1_all)
  
  S2 <- functS2(list = allResults, weighted = F, interaction = F)
  S2w <- functS2(list = allResults, weighted = T, interaction = F)
  S2_I <- functS2(list = allResults, weighted = F, interaction = T)
  S2w_I <- functS2(list = allResults, weighted = T, interaction = T)
  
  S2_all <- rbind(S2, S2w, S2_I, S2w_I)
  S2_df <- rbind(S2_df, S2_all)
  
  S3 <- functS3(list = allResults, weighted = F, interaction = F)
  S3w <- functS3(list = allResults, weighted = T, interaction = F)
  S3_I <- functS3(list = allResults, weighted = F, interaction = T)
  S3w_I <- functS3(list = allResults, weighted = T, interaction = T)
  
  S3_all <- rbind(S3, S3w, S3_I, S3w_I)
  S3_df <- rbind(S3_df, S3_all)
}

# Agrupar por 'mix', 'Model' y 'Strategy' y calcular la media y el error estándar
S1_df <- S1_df %>% 
  group_by(mix, Model, Strategy) %>% 
  summarise(
    Mean = mean(value),
    SE = sd(value, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'  # Para evitar advertencias sobre agrupamiento
  )

S2_df <- S2_df %>%
  group_by(mix, Model, strategy) %>%                          # Agrupar por la columna 'mix'
  summarize(
    Mean = mean(value, na.rm = TRUE),        # Calcular la media por grupo
    SE = sd(value, na.rm = TRUE) / sqrt(n()),# Calcular el error estándar por grupo
    .groups = 'drop'                         # Eliminar la agrupación después del resumen
  )

S3_df <- S3_df %>%
  group_by(mix, Model, strategy) %>%                          # Agrupar por la columna 'mix'
  summarize(
    Mean = mean(value, na.rm = TRUE),        # Calcular la media por grupo
    SE = sd(value, na.rm = TRUE) / sqrt(n()),# Calcular el error estándar por grupo
    .groups = 'drop'                         # Eliminar la agrupación después del resumen
  ) 

colnames(S1_df) <- colnames(S2_df) <- colnames(S3_df)

S1_plot <- plotCV(S1_df)
S2_plot <- plotCV(S2_df)
S3_plot <- plotCV(S3_df)

plots <- list(S1_plot, S2_plot, S3_plot)
names <- list("Strategy1", "Strategy2",  "Strategy3")
for(i in seq_along(plots)){
  plot <- plots[[i]]
  name <- names[[i]]
  png(paste0("outputs/plots/", name, ".png"), width = 4000, height = 2000, res = 400)
  plot(plot)
  dev.off()
}







